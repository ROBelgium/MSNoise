import unittest
from msnoise.msnoise_table_def import Config


class TestConfig(unittest.TestCase):

    fields_to_skip = ["data_folder", "output_folder", "response_path"]

    before_stripping = [
        "  xxx  xxx   ",
        "  xxx  xxx",
        "xxx  xxx   ",
        "  xxxxxx   ",
        "xxx  xxx",
    ]

    after_stripping =  [
        "xxx  xxx",
        "xxx  xxx",
        "xxx  xxx",
        "xxxxxx",
        "xxx  xxx",
    ]

    after_replacing = "xxxxxx"


    def test_strip_value_path_field(self):

        for name in self.fields_to_skip:
            for before, after in zip(self.before_stripping,
                                     self.after_stripping):
                self.assertEqual(after, Config(name, before))

    def test_strip_value_non_path_field(self):
        name = "definitely_a_name_not_in_config"

        for before in self.before_stripping:
            self.assertEqual(self.after_replacing, Config(name, before))

def suite():
    testsuite = unittest.TestSuite()
    testsuite.addTest(unittest.makeSuite(TestConfig))
    return testsuite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
