from unittest import TestCase, TestSuite, TextTestRunner
from datetime import date

from msnoise.api import build_movstack_datelist


class TestBuild_movstack_datelist(TestCase):
    """
    Test class to perform tests of build_movstack_datelist method.
    """
    proper_startdate_str = '2017-01-01'
    proper_enddate_str = '2017-01-10'

    proper_startdate = date(2017, 1, 1)
    proper_enddate = date(2017, 1, 10)

    improper_dates = [('2017.01.01', '2017.01.10'),
                      ('2017.01.01', '2017.02.30'),
                      ('01.01.2017', '12.31.2017'),
                      ('2017 01 01', '2017 01 10'),
                      ('2017;01;01', '2017;01;10')]

    def setUp(self):
        pass

    def test_custom_dates_proper_strings(self):
        start, end, list = build_movstack_datelist('session',
                                                   self.proper_startdate_str,
                                                   self.proper_enddate_str)
        self.assertEqual(start, self.proper_startdate)
        self.assertEqual(end, self.proper_enddate)

    def test_custom_dates_improper_strings(self):
        for pair in self.improper_dates:
            with self.assertRaises(ValueError):
                build_movstack_datelist('session', pair[0], pair[1])

    def test_dates_from_db(self):
        self.fail()

    def tearDown(self):
        pass


class TestGet_config(TestCase):
    def test_get_config(self):
        self.fail()


def suite():
    suite = TestSuite()
    suite.addTest(
        TestBuild_movstack_datelist('test_custom_dates_proper_strings'))
    suite.addTest(
        TestBuild_movstack_datelist('test_custom_dates_improper_strings'))
    suite.addTest(TestBuild_movstack_datelist('test_dates_from_db'))
    suite.addTest(TestGet_config('test_get_config'))

    return suite

if __name__ == '__main__':
    runner = TextTestRunner()
    runner.run(suite())