import unittest
import traceback
import os
import datetime


# Here's our "unit tests".
class MSNoiseTests(unittest.TestCase):
    def test_001_S01installer(self):
        from ..s000installer import main
        try:
            ret = main(tech=1)
            msg = "Installation Done! - Go to Configuration Step!"
            self.failUnlessEqual(ret, msg)
        except:
            traceback.print_exc()
            self.fail()

    def test_002_ConnectToDB(self):
        from ..api import connect
        try:
            db = connect()
            db.close()
        except:
            self.fail("Can't connect to MSNoise DB")

    def test_003_set_and_config(self):
        from ..api import connect, get_config, update_config
        db = connect()
        totests = []
        path = os.path.abspath(os.path.dirname(__file__))
        totests.append(['data_folder', os.path.join(path, 'data')])
        totests.append(['data_structure', 'PDF'])
        totests.append(['network', 'YA'])
        totests.append(['ZR', 'N'])
        totests.append(['ZT', 'N'])
        totests.append(['TZ', 'N'])
        totests.append(['TR', 'N'])
        totests.append(['TT', 'N'])
        totests.append(['RZ', 'N'])
        totests.append(['RR', 'N'])
        totests.append(['RT', 'N'])

        for test in totests:
            update_config(db, test[0], test[1])
            d = get_config(db, test[0])
            self.failUnlessEqual(test[1], d)

        db.close()

    def test_004_set_and_get_filters(self):
        from ..msnoise_table_def import Filter
        from ..api import connect, update_filter, get_filters
        db = connect()
        filters = []
        f = Filter()
        f.low = 0.01
        f.mwcs_low = 0.12
        f.high = 1.0
        f.mwcs_high = 0.98
        f.rms_threshold = 0
        f.mwcs_wlen = 10
        f.mwcs_step = 5
        f.used = True
        filters.append(f)
        f = Filter()
        f.low = 0.1
        f.mwcs_low = 0.12
        f.high = 1.0
        f.mwcs_high = 0.98
        f.rms_threshold = 0
        f.mwcs_wlen = 10
        f.mwcs_step = 5
        f.used = True
        filters.append(f)

        for f in filters:
            update_filter(db, f.ref, f.low, f.mwcs_low, f.high, f.mwcs_high,
                          f.rms_threshold, f.mwcs_wlen, f.mwcs_step, f.used)

        dbfilters = get_filters(db)
        for i, filter in enumerate(dbfilters):
            for param in ['low', 'mwcs_low', 'high', 'mwcs_high',
                          'rms_threshold', 'mwcs_wlen', 'mwcs_step', 'used']:
                self.failUnlessEqual(eval("filter.%s" % param),
                                     eval("filters[i].%s" % param))

    def test_005_populate_station_table(self):
        from ..s002populate_station_table import main
        try:
            ret = main()
            self.failUnlessEqual(ret, True)
        except:
            self.fail()

    def test_006_get_stations(self):
        from ..api import connect, get_stations
        db = connect()
        stations = get_stations(db).all()
        self.failUnlessEqual(len(stations), 3)
        db.close()

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

    def test_009_control_data_availability(self):
        from ..api import connect, get_new_files, get_data_availability,\
            count_data_availability_flags, get_stations

        db = connect()
        files = get_new_files(db)
        self.failUnlessEqual(len(files), 3)

        flags = count_data_availability_flags(db)
        self.failUnlessEqual(len(flags), 1)

        for station in get_stations(db):
            da = get_data_availability(db, net=station.net, sta=station.sta,
                                       comp='HHZ')
            self.failUnlessEqual(len(da), 1)

    def test_010_new_jobs(self):
        from ..s02new_jobs import main

        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_011_control_jobs(self):
        from ..api import connect, is_next_job, get_next_job
        from ..msnoise_table_def import Job
        db = connect()

        self.failUnlessEqual(is_next_job(db), True)
        jobs = get_next_job(db)
        self.failUnlessEqual(isinstance(jobs[0], Job), True)

    def test_012_reset_jobs(self):
        from ..api import connect, reset_jobs
        db = connect()
        reset_jobs(db, 'CC')
        db.close()

    def test_012b_hack_noresample(self):
        from ..api import connect, update_config
        db = connect()
        update_config(db,'resampling_method','Decimate')
        db.close()

    def test_013_s03compute_cc(self):
        from ..s03compute_cc import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_014_check_done_jobs(self):
        from ..api import connect, get_job_types
        db = connect()
        jobs = get_job_types(db, 'CC')
        self.failUnlessEqual(jobs[0][0], 3)
        self.failUnlessEqual(jobs[0][1], 'D')

    def test_015_check_cc_files(self, format="MSEED"):
        from ..api import connect, get_filters, get_station_pairs, \
            get_components_to_compute
        db = connect()
        for filter in get_filters(db):
            for components in get_components_to_compute(db):
                for (sta1, sta2) in get_station_pairs(db):
                    pair = "%s_%s_%s_%s" % (sta1.net, sta1.sta,
                                            sta2.net, sta2.sta)
                    tmp = os.path.join("STACKS",
                                       "%02i" % filter.ref,
                                       "001_DAYS",
                                       components,
                                       pair,
                                       "2010-09-01.%s" % format)
                    print("checking", tmp)
                    if not os.path.isfile(tmp):
                        self.fail()

    def test_016_update_config(self):
        from ..api import connect, update_config
        db = connect()
        update_config(db, "export_format", "SAC")

    def test_017_reset_cc_jobs(self):
        from ..api import connect, reset_jobs
        db = connect()
        reset_jobs(db, 'CC', alljobs=True)

    def test_018_recompute_cc(self):
        self.test_013_s03compute_cc()

    def test_019_check_SACS(self):
        self.test_015_check_cc_files(format='SAC')

    def test_020_update_config(self):
        from ..api import connect, update_config
        import shutil
        shutil.rmtree("STACKS")
        db = connect()
        update_config(db, "export_format", "BOTH")

    def test_021_reprocess_BOTH(self):
        self.test_017_reset_cc_jobs()
        self.test_013_s03compute_cc()
        self.test_015_check_cc_files(format='MSEED')
        self.test_015_check_cc_files(format='SAC')

    def test_022_check_content(self):
        from obspy.core import read
        from numpy.testing import assert_allclose
        from ..api import connect, get_filters, get_station_pairs, \
            get_components_to_compute
        db = connect()
        for filter in get_filters(db):
            for components in get_components_to_compute(db):
                for (sta1, sta2) in get_station_pairs(db):
                    pair = "%s_%s_%s_%s" % (sta1.net, sta1.sta,
                                            sta2.net, sta2.sta)
                    tmp1 = os.path.join("STACKS",
                                        "%02i" % filter.ref,
                                        "001_DAYS",
                                        components,
                                        pair,
                                        "2010-09-01.MSEED")
                    tmp2 = os.path.join("STACKS",
                                        "%02i" % filter.ref,
                                        "001_DAYS",
                                        components,
                                        pair,
                                        "2010-09-01.SAC")
                    tmp1 = read(tmp1)
                    tmp2 = read(tmp2)
                    assert_allclose(tmp1[0].data, tmp2[0].data)

    def test_023_stack(self):
        from ..api import connect, update_config
        from ..s04stack import main
        db = connect()
        update_config(db, 'ref_begin', '2009-01-01')
        update_config(db, 'ref_end', '2011-01-01')
        update_config(db, 'startdate', '2009-01-01')
        update_config(db, 'enddate', '2011-01-01')
        db.close()
        interval = 1.
        main('ref', interval)
        main('mov', interval)
        main('step', interval)

    def test_024_mwcs(self):
        from ..s05compute_mwcs import main
        main()

    def test_025_dtt(self):
        from ..s06compute_dtt import main
        main()

    def test_026_build_ref_datelist(self):
        from ..api import connect, build_ref_datelist
        db = connect()
        start, end, datelist = build_ref_datelist(db)
        self.failUnlessEqual(start, datetime.date(2009, 1, 1))
        self.failUnlessEqual(end, datetime.date(2011, 1, 1))
        self.failUnlessEqual(len(datelist), 731)

    def test_027_build_movstack_datelist(self):
        from ..api import connect, build_movstack_datelist
        db = connect()
        start, end, datelist = build_movstack_datelist(db)
        self.failUnlessEqual(start, datetime.date(2009, 1, 1))
        self.failUnlessEqual(end, datetime.date(2011, 1, 1))
        self.failUnlessEqual(len(datelist), 731)

    def test_028_S01installer(self):
        if not "TRAVIS" in os.environ:
            print("Seems to be running on local machine, skipping MySQL test")
            return
        import shutil
        shutil.move('db.ini','db.bak')
        from ..s000installer import main
        try:
            ret = main(tech=2, username="root", password="",
                       hostname="localhost", database="msnoise")
            msg = "Installation Done! - Go to Configuration Step!"
            self.failUnlessEqual(ret, msg)
        except:
            traceback.print_exc()
            self.fail()


def main():
    import os
    import sys
    c = len(os.listdir(os.getcwd()))
    if c > 0:
        print("Directory is not empty, can't run tests here!")
        sys.exit()

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MSNoiseTests)
    runner = unittest.TextTestRunner(verbosity=2)
    result = runner.run(suite)
    if not result.wasSuccessful():
        sys.exit(1)

if __name__ == '__main__':
    main()
