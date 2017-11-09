from unittest import makeSuite, TestCase, TestSuite, TextTestRunner
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from datetime import date

from msnoise.api import build_movstack_datelist, get_config, update_config
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


class TestUpdateConfigAndGetConfig(TestCase):
    """
    Test for msnoise.api.testconfig method

    Those tests depend on both UpdateConfig and GetConfig methods.
    It should be separated in the future.
    """

    startdate_name = "startdate"
    analysis_duration_name = "analysis_duration"
    cc_sampling_rate_name = "cc_sampling_rate"

    startdate_str = "2017-01-01"
    analysis_duration_int = str(3600)
    cc_sampling_rate_float = str(20.0)

    startdate_str_final = "2013-12-21"
    analysis_duration_int_final = str(45000)
    cc_sampling_rate_float_final = str(48.3)

    def setUp(self):
        # Create a in memory database only once for test suite
        engine = create_engine('sqlite:///')
        Base.metadata.create_all(engine)

        make_session = sessionmaker(bind=engine)
        self.session = make_session()

        self.session.add(Config(name=self.startdate_name,
                                value=self.startdate_str))
        self.session.add(Config(name=self.analysis_duration_name,
                                value=self.analysis_duration_int))
        self.session.add(Config(name=self.cc_sampling_rate_name,
                                value=self.cc_sampling_rate_float))
        self.session.commit()

    def test_update__and_get_config(self):
        update_config(self.session, self.startdate_name,
                      self.startdate_str_final)
        self.assertEqual(get_config(self.session, self.startdate_name),
                         self.startdate_str_final)

        update_config(self.session, self.analysis_duration_name,
                      self.analysis_duration_int_final)
        self.assertEqual(get_config(self.session, self.analysis_duration_name),
                         self.analysis_duration_int_final)

        update_config(self.session, self.cc_sampling_rate_name,
                      self.cc_sampling_rate_float_final)
        self.assertEqual(get_config(self.session, self.cc_sampling_rate_name),
                         self.cc_sampling_rate_float_final)

    def tearDown(self):
        self.session.close()


def suite():
    testsuite = TestSuite()
    testsuite.addTest(makeSuite(TestBuildMovstackDatelist))
    testsuite.addTest(makeSuite(TestUpdateConfigAndGetConfig))

    return testsuite


if __name__ == '__main__':
    runner = TextTestRunner()
    runner.run(suite())
