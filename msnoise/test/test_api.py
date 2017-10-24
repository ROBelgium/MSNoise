from unittest import makeSuite, TestCase, TestSuite, TextTestRunner
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from datetime import date

from msnoise.api import build_movstack_datelist
from msnoise.msnoise_table_def import Base, Config


class TestBuildMovstackDatelist(TestCase):
    """
    Test class to perform tests of msnoise.api.build_movstack_datelist method.
    """
    proper_startdate_str = '2017-01-01'
    proper_enddate_str = '2017-01-10'

    proper_startdate = date(2017, 1, 1)
    proper_enddate = date(2017, 1, 10)
    proper_list_of_dates = [date(2017, 1, 1), date(2017, 1, 2),
                            date(2017, 1, 3), date(2017, 1, 4),
                            date(2017, 1, 5), date(2017, 1, 6),
                            date(2017, 1, 7), date(2017, 1, 8),
                            date(2017, 1, 9), date(2017, 1, 10)]

    improper_dates = [('2017.01.01', '2017.01.10'),
                      ('2017.01.01', '2017.02.30'),
                      ('01.01.2017', '12.31.2017'),
                      ('2017 01 01', '2017 01 10'),
                      ('2017;01;01', '2017;01;10')]

    def setUp(self):
        # Create a in memory database only once for test suite
        engine = create_engine('sqlite:///')
        Base.metadata.create_all(engine)

        make_session = sessionmaker(bind=engine)
        self.session = make_session()

        self.session.add(Config(name='startdate',
                                value=self.proper_startdate_str))
        self.session.add(Config(name='enddate',
                                value=self.proper_enddate_str))
        self.session.commit()

    def test_custom_dates_proper_strings(self):
        start, end, dates = build_movstack_datelist(self.session,
                                                    self.proper_startdate_str,
                                                    self.proper_enddate_str)
        self.assertEqual(start, self.proper_startdate)
        self.assertEqual(end, self.proper_enddate)
        self.assertEqual(dates, self.proper_list_of_dates)

    def test_custom_dates_improper_strings(self):
        for pair in self.improper_dates:
            with self.assertRaises(ValueError):
                build_movstack_datelist(self.session, pair[0], pair[1])

    def test_dates_from_db(self):
        start, end, dates = build_movstack_datelist(self.session)

        self.assertEqual(start, self.proper_startdate)
        self.assertEqual(end, self.proper_enddate)
        self.assertEqual(dates, self.proper_list_of_dates)

    def tearDown(self):
        self.session.close()


class TestGetConfig(TestCase):
    def test_get_config(self):
        self.fail()

class TestGetResults(TestCase):
    def test_get_results_fail_non_esisting_path(self):
        self.fail()


def suite():
    testsuite = TestSuite()
    testsuite.addTest(makeSuite(TestBuildMovstackDatelist))
    testsuite.addTest(makeSuite(TestGetConfig))

    return testsuite


if __name__ == '__main__':
    runner = TextTestRunner()
    runner.run(suite())
