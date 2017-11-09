from unittest import makeSuite, TestCase, TestSuite, TextTestRunner
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from datetime import date

from msnoise.api import build_movstack_datelist, get_config, update_config,\
    get_results
from msnoise.msnoise_table_def import Base, Config


class TestGetResults(TestCase):
    def setUp(self):
        # Create a in memory database only once for test suite
        engine = create_engine('sqlite:///')
        Base.metadata.create_all(engine)

        make_session = sessionmaker(bind=engine)
        self.session = make_session()

        self.session.add(Config(name="maxlag",
                                value=120.0))
        self.session.add(Config(name='cc_sampling_rate',
                                value=12.0))


    def test_get_results_fail_non_esisting_path(self):
        with self.assertRaises(NotADirectoryError):
            get_results(self.session, "AABBGGRR", "BBCCDD", 1, "ZZ",
                        [0,1,2,3,4])


def suite():
    testsuite = TestSuite()
    testsuite.addTest(makeSuite(TestGetResults))

    return testsuite


if __name__ == '__main__':
    runner = TextTestRunner()
    runner.run(suite())
