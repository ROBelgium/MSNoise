import datetime
import glob
import logging
import os
import re
import shutil
import traceback
import unittest

logger = logging.getLogger('matplotlib')
# set WARNING for Matplotlib
logger.setLevel(logging.CRITICAL)
import tempfile
from click.testing import CliRunner
from obspy import read
from .. import FatalError, s01scan_archive
from ..scripts import msnoise as msnoise_script


class MSNoiseTests(unittest.TestCase):
    prefix = ""

    def setUp(self):
        # Copy test/data directory to ./data
        path = os.path.abspath(os.path.dirname(__file__))
        data_folder = os.path.join(path, 'data')
        if not os.path.isdir("data"):
            shutil.copytree(data_folder, "data/")
        self.data_folder = "data"
        # Create a click runner
        self.runner = CliRunner()

    def test_001_S01installer(self):
        import os
        from ..s000installer import main
        if "PREFIX" in os.environ:
            self.prefix=os.environ["PREFIX"]
        try:
            ret = main(tech=1, prefix=self.prefix,
                       filename='testmsnoise.sqlite')
            self.failUnlessEqual(ret, 0)
        except:
            traceback.print_exc()
            self.fail()

    def test_003_set_and_config(self):
        from ..api import connect, get_config, update_config
        db = connect()
        totests = []
        totests.append(['data_folder', self.data_folder])
        totests.append(['data_structure', 'PDF'])
        totests.append(['network', 'YA'])
        totests.append(['components_to_compute', 'ZZ'])
        totests.append(['mov_stack', '1,2,5'])

        for test in totests:
            update_config(db, test[0], test[1])
            d = get_config(db, test[0])
            self.failUnlessEqual(test[1], d)

        db.close()

    def test_004_set_and_get_filters(self):
        from ..api import connect, update_filter, get_filters, Filter
        db = connect()
        filters = []
        f = Filter()
        f.low = 0.1
        f.mwcs_low = 0.12
        f.high = 1.0
        f.mwcs_high = 0.98
        f.mwcs_wlen = 10
        f.mwcs_step = 5
        f.used = True
        filters.append(f)
        for f in filters:
            update_filter(db, f.ref, f.low, f.mwcs_low, f.high, f.mwcs_high,
                          f.mwcs_wlen, f.mwcs_step, f.used)

    def test_005_populate_station_table(self):
        from ..s002populate_station_table import main
        try:
            ret = main()
            self.failUnlessEqual(ret, True)
        except:
            self.fail()

    def test_007_update_stations(self):
        from ..api import connect, get_stations, update_station
        import pandas as pd
        db = connect()
        path = os.path.split(os.path.abspath(os.path.dirname(__file__)))[0]
        stations = pd.read_csv(os.path.join(path, 'test', 'extra',
                                            'stations.csv'),
                               header=None, index_col=0,
                               names=['X', 'Y', 'altitude'])
        for station in get_stations(db):
            fullname = "%s.%s" % (station.net, station.sta,)
            try:
                s = stations.loc[fullname]
                update_station(db, station.net, station.sta, s['X'],
                               s['Y'], s['altitude'])
            except:
                traceback.print_exc()
                self.fail()
                continue
            del s

    def test_008_scan_archive(self):
        from ..s01scan_archive import main
        try:
            main(init=True, threads=1)
        except:
            traceback.print_exc()
            self.fail()

    def test_008b_add_loc_chan_to_stations(self):
        result = self.runner.invoke(msnoise_script.da_stations_update_loc_chan)

    def test_010_new_jobs(self):
        from ..s02new_jobs import main

        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_012b_hack_noresample(self):
        from ..api import connect, update_config
        db = connect()
        update_config(db, 'resampling_method', 'Decimate')
        db.close()

    def test_013_s03compute_cc(self):
        from ..s03compute_no_rotation import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_014_check_results(self):
        # TO DO CHECK CCF RESULTS ARRAYS!
        pass

    def test_023_stack(self):
        from ..api import connect, update_config, reset_jobs
        from ..s04stack import main
        db = connect()
        update_config(db, 'ref_begin', '2009-01-01')
        update_config(db, 'ref_end', '2011-01-01')
        update_config(db, 'startdate', '2009-01-01')
        update_config(db, 'enddate', '2011-01-01')

        interval = 1.
        main('ref', interval)
        reset_jobs(db, "STACK", alljobs=True)
        main('mov', interval)
        reset_jobs(db, "STACK", alljobs=True)
        main('step', interval)
        db.close()

    def test_024_mwcs(self):
        from ..s05compute_mwcs import main
        main()

    def test_025_dtt(self):
        from ..s06compute_dtt import main
        main()


def main(prefix=""):
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")

    import os
    import sys
    test_dir = tempfile.mkdtemp()
    os.chdir(test_dir)
    print("Tests will be executed in %s" % test_dir)
    # c = len(os.listdir(os.getcwd()))
    # if c > 0:
    #     print("Directory is not empty, can't run tests here!")
    #     sys.exit()
    os.environ["PREFIX"] = prefix
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MSNoiseTests)
    runner = unittest.TextTestRunner(verbosity=4)
    result = runner.run(suite)
    print("Tests executed in %s" % test_dir)
    if not result.wasSuccessful():
        sys.exit(1)


if __name__ == '__main__':
    main()
