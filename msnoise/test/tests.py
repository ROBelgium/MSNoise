import datetime
import glob
import logging
import os
import re
import shutil
import traceback
import unittest

import pooch

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

        test_files = {
            "data/2010/UV05/HHZ.D/YA.UV05.00.HHZ.D.2010.244": "17034091285d485f7c2d4797f435228c408d6940db943be63f1769ec09854f4f",
            "data/2010/UV06/HHZ.D/YA.UV06.00.HHZ.D.2010.244": "51bfd1e735696e83ee6dba136c9e740c59120fac9f74b386eac75062eb9ca382",
            "data/2010/UV10/HHZ.D/YA.UV10.00.HHZ.D.2010.244": "530cc7f4a57fe69a8a5cedeb18e64773055c146e4ae4676012f6618dd0c92e82",
            "extra/DATA.RESIF_Jun_10,14_21_05_20264.RESIF": "95a6d007132fc41b6107d258aeee1170614d234cdd3eb4a6d5652e4661a6adcd",
            # "extra/RESP.GD.PUT..SHZ.L4C.5500.8905.1.171": "bubu",
            "extra/stations.csv": "057152c2823c5457bce879146d78984af422973ab313e1d1cd8baaa7f7a1d6b3",
            "extra/test_inventory.xml": "48bf3261a9c23f1782e452583306b73643a4d7f205df4338e5525df3c06eccb9",
        }

        BRIAN = pooch.create(
            # This is still the default
            path=pooch.os_cache("msnoise"),
            base_url="https://github.com/ROBelgium/MSNoise/raw/{version}/msnoise/test",
            version="1.6",
            version_dev="main",
            registry=test_files,
            # The name of an environment variable that can overwrite the path
            env="MSNOISE_DATA_DIR",
        )
        for fn in test_files:
            path = BRIAN.fetch(fn)
        # Copy test/data directory to ./data
        # path = os.path.abspath(os.path.dirname(__file__))
        # data_folder = os.path.join(path, 'data')

        data_folder = os.path.join(os.environ.get("MSNOISE_DATA_DIR",pooch.os_cache("msnoise")), "1.6", "data")
        self.response_folder = os.path.join(os.environ.get("MSNOISE_DATA_DIR",pooch.os_cache("msnoise")), "1.6", "extra")
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
            if os.environ["TECH"] == "1":
                ret = main(tech=1, prefix=self.prefix,
                           hostname="localhost",
                           username="root",
                           password="SECRET",
                           database=os.environ["hash"])
            elif os.environ["TECH"] == "2":
                ret = main(tech=2, prefix=self.prefix,
                           hostname=os.environ["MARIADB_HOSTNAME"],
                           username=os.environ["MARIADB_USERNAME"],
                           password=os.environ["MARIADB_PASSWORD"],
                           database=os.environ["hash"])

            self.assertEqual(ret, 0)
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
        totests.append(['data_folder', self.data_folder])
        totests.append(['data_structure', 'PDF'])
        totests.append(['network', 'YA'])
        totests.append(['components_to_compute', 'ZZ'])

        for test in totests:
            update_config(db, test[0], test[1])
            d = get_config(db, test[0])
            self.assertEqual(test[1], d)

        db.close()

    def test_004_set_and_get_filters(self):
        from ..api import connect, update_filter, get_filters, Filter
        db = connect()
        filters = []
        f = Filter()
        f.low = 0.01
        f.mwcs_low = 0.12
        f.high = 1.0
        f.mwcs_high = 0.98
        f.mwcs_wlen = 10
        f.mwcs_step = 5
        f.used = True
        filters.append(f)
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

        dbfilters = get_filters(db)
        for i, filter in enumerate(dbfilters):
            for param in ['low', 'mwcs_low', 'high', 'mwcs_high',
                          'mwcs_wlen', 'mwcs_step', 'used']:
                self.assertEqual(eval("filter.%s" % param),
                                     eval("filters[i].%s" % param))

    def test_005_populate_station_table(self):
        from ..s002populate_station_table import main
        try:
            ret = main()
            self.assertEqual(ret, True)
        except:
            self.fail()

    def test_006_get_stations(self):
        from ..api import connect, get_stations
        db = connect()
        stations = get_stations(db).all()
        self.assertEqual(len(stations), 3)
        db.close()

    def test_007_update_stations(self):
        from ..api import connect, get_stations, update_station
        import pandas as pd
        db = connect()
        stations = pd.read_csv(os.path.join(self.response_folder,
                                            'stations.csv'),
                               header=None, index_col=0,
                               names=['X', 'Y', 'altitude'])
        for station in get_stations(db):
            fullname = "%s.%s" % (station.net, station.sta,)
            try:
                s = stations.loc[fullname]
                update_station(db, net=station.net, sta=station.sta,
                               X=s['X'], Y=s['Y'], altitude=s['altitude'])
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
        result = self.runner.invoke(msnoise_script.db_da_stations_update_loc_chan)

    def test_009_control_data_availability(self):
        from ..api import connect, get_new_files, get_data_availability,\
            count_data_availability_flags, get_stations

        db = connect()
        files = get_new_files(db)
        self.assertEqual(len(files), 3)

        flags = count_data_availability_flags(db)
        self.assertEqual(len(flags), 1)

        for station in get_stations(db):
            for loc in station.locs():
                for chan in station.chans():
                    da = get_data_availability(db, net=station.net,
                                               sta=station.sta,
                                               loc=loc,
                                               chan=chan)
                    self.assertEqual(len(da), 1)

    def test_010_new_jobs(self):
        from ..s02new_jobs import main

        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_011_control_jobs(self):
        from ..api import connect, is_next_job, get_next_job, Job
        db = connect()

        self.assertEqual(is_next_job(db), True)
        jobs = get_next_job(db)
        self.assertEqual(isinstance(jobs[0], Job), True)

    def test_012_reset_jobs(self):
        from ..api import connect, reset_jobs
        db = connect()
        reset_jobs(db, 'CC', alljobs=True)
        db.close()

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

    def test_013_X_check_results(self):
        # TO DO CHECK CCF RESULTS ARRAYS!
        pass

    def test_014_check_done_jobs(self):
        from ..api import connect, get_job_types
        db = connect()
        jobs = get_job_types(db, 'CC')
        self.assertEqual(jobs[0][0], 3)
        self.assertEqual(jobs[0][1], 'D')

    def test_015_check_cc_files(self, format="MSEED"):
        from ..api import connect, get_filters, get_station_pairs, \
            get_components_to_compute
        db = connect()
        for filter in get_filters(db):
            for components in get_components_to_compute(db):
                for (sta1, sta2) in get_station_pairs(db):
                    for loc1 in sta1.locs():
                        for loc2 in sta2.locs():
                            pair = "%s.%s.%s_%s.%s.%s" % (sta1.net, sta1.sta,
                                                          loc1, sta2.net,
                                                          sta2.sta, loc2)
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
        db.close()

    def test_017_reset_cc_jobs(self):
        from ..api import connect, reset_jobs
        db = connect()
        reset_jobs(db, 'CC', alljobs=True)
        db.close()

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
        db.close()

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
                    for loc1 in sta1.locs():
                        for loc2 in sta2.locs():
                            sta1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                            sta2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                            pair = "%s_%s" % (sta1, sta2)
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
        db.close()

    def test_023_stack(self):
        from ..api import connect, update_config, reset_jobs
        from ..s04_stack2 import main
        db = connect()
        update_config(db, 'ref_begin', '2009-01-01')
        update_config(db, 'ref_end', '2011-01-01')
        update_config(db, 'startdate', '2009-01-01')
        update_config(db, 'enddate', '2011-01-01')
        update_config(db, 'mov_stack', "(('1d','1d'),('2d','1d'),('5d','1d'))")

        interval = 1.
        main('ref', interval)
        reset_jobs(db, "STACK", alljobs=True)
        main('mov', interval)
        reset_jobs(db, "STACK", alljobs=True)
        main('step', interval)
        db.close()

    def test_024_mwcs(self):
        from ..s05compute_mwcs2 import main
        main()

    def test_025_dtt(self):
        from ..s06compute_dtt2 import main
        main()

    def test_026_dvv(self):
        from ..s07_compute_dvv import main
        main()

    def test_027_build_ref_datelist(self):
        from ..api import connect, build_ref_datelist
        db = connect()
        start, end, datelist = build_ref_datelist(db)
        self.assertEqual(start, datetime.date(2009, 1, 1))
        self.assertEqual(end, datetime.date(2011, 1, 1))
        self.assertEqual(len(datelist), 731)
        db.close()

    def test_028_build_movstack_datelist(self):
        from ..api import connect, build_movstack_datelist
        db = connect()
        start, end, datelist = build_movstack_datelist(db)
        self.assertEqual(start, datetime.date(2009, 1, 1))
        self.assertEqual(end, datetime.date(2011, 1, 1))
        self.assertEqual(len(datelist), 731)
        db.close()

    def test_029_stretching(self):
        from ..api import connect, update_config, reset_jobs
        db = connect()
        update_config(db, "export_format", "MSEED")
        reset_jobs(db, "MWCS", alljobs=True)
        db.close()

        from ..stretch import main
        main()

    def test_030_create_fake_new_files(self):
        for f in sorted(glob.glob(os.path.join(self.data_folder, "2010", "*",
                                               "HHZ.D", "*"))):
            st = read(f)
            for tr in st:
                tr.stats.starttime += datetime.timedelta(days=1)
            out = f.replace("244", "245")
            st.write(out, format="MSEED")

        from ..s01scan_archive import main
        main(init=False, threads=1)

        from ..s02new_jobs import main
        main()

        from ..api import connect, get_job_types
        db = connect()
        jobs = get_job_types(db, 'CC')
        print(jobs)
        self.assertEqual(jobs[0][0], 3)
        self.assertEqual(jobs[0][1], 'D')
        self.assertEqual(jobs[1][0], 3)
        self.assertEqual(jobs[1][1], 'T')

    def test_031_instrument_response(self):
        from ..api import connect, update_config
        db = connect()
        update_config(db, 'response_path', self.response_folder)
        update_config(db, 'response_format', "dataless")
        update_config(db, 'remove_response', "Y")
        db.close()
        self.test_013_s03compute_cc()


    # def test_031_compute_cc_rot(self):
    #     import shutil
    #     shutil.rmtree("STACKS")
    #     from ..api import connect, reset_jobs
    #     db = connect()
    #     reset_jobs(db, "CC", alljobs=True)
    #     db.close()
    #     from ..s03compute_no_rotation import main
    #     main()

    def test_032_wct(self):
        from ..api import connect, read_db_inifile
        from sqlalchemy import text
        db = connect()
        dbini = read_db_inifile()
        prefix = (dbini.prefix + '_') if dbini.prefix != '' else ''
        db.execute(text("INSERT INTO {prefix}jobs (pair, day, jobtype, flag) "
                       "SELECT pair, day, '{right_type}', 'T' FROM {prefix}jobs "
                       "WHERE jobtype='{left_type}' AND flag='D';"
                       .format(prefix=prefix, right_type="WCT", left_type="STACK")))
        db.commit()

        from ..s08compute_wct import main
        main()
        db.close()
    
    # PLOTS

    def test_100_plot_cctfime(self):
        from ..api import connect, get_station_pairs, get_filters
        from ..plots.ccftime import main
        db = connect()
        for sta1, sta2 in get_station_pairs(db):
            for loc1 in sta1.locs():
                for loc2 in sta2.locs():
                    sta1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                    sta2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                    for filter in get_filters(db):
                        main(sta1, sta2, filter.ref, "ZZ",  1, show=False,
                             outfile="?.png")
                        fn = 'ccftime %s-%s-f%i-m%s_%s.png' % \
                             ("%s-%s" % (sta1, sta2), "ZZ", filter.ref, "1d", "1d")
                        self.assertTrue(os.path.isfile(fn), msg="%s doesn't exist" % fn)

    def test_100_plot_interferogram(self):
        from ..api import connect, get_station_pairs, get_filters
        from ..plots.interferogram import main
        db = connect()
        for sta1, sta2 in get_station_pairs(db):
            for loc1 in sta1.locs():
                for loc2 in sta2.locs():
                    sta1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                    sta2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                    for filter in get_filters(db):
                        main(sta1, sta2, filter.ref, "ZZ",  1, show=False,
                             outfile="?.png")
                        fn = 'interferogram %s-%s-f%i-m%s_%s.png' % \
                             ("%s-%s" % (sta1,
                                         sta2),
                              "ZZ", filter.ref, "1d", "1d")
                        self.assertTrue(os.path.isfile(fn), msg="%s doesn't exist" % fn)

    def test_101_plot_spectime(self):
        from ..api import connect, get_station_pairs, get_filters
        from ..plots.spectime import main
        db = connect()
        for sta1, sta2 in get_station_pairs(db):
            for loc1 in sta1.locs():
                for loc2 in sta2.locs():
                    sta1 = "%s.%s.%s" % (sta1.net, sta1.sta, loc1)
                    sta2 = "%s.%s.%s" % (sta2.net, sta2.sta, loc2)
                    for filter in get_filters(db):
                        main(sta1, sta2, filter.ref, "ZZ", 1, show=False,
                             outfile="?.png")
                        fn = 'spectime %s-%s-f%i-m%s_%s.png' % \
                             ("%s-%s" % (sta1,
                                         sta2),
                              "ZZ", filter.ref, "1d", "1d")
                        self.assertTrue(os.path.isfile(fn), msg="%s doesn't exist" % fn)

    def test_102_plot_distance(self):
        from ..plots.distance import main
        main(filterid=1, components="ZZ", show=False, outfile="?.png")
        fn = "distance ZZ-f1.png"
        self.assertTrue(os.path.isfile(fn),
                        msg="%s doesn't exist" % fn)

        main(filterid=1, components="ZZ", show=False, outfile="?_refilter.png",
             refilter="0.2:0.9")
        fn = "distance ZZ-f1_refilter.png"
        self.assertTrue(os.path.isfile(fn),
                        msg="%s doesn't exist" % fn)

    def test_103_plot_dvv(self):
        from ..plots.dvv import main
        main(filterid=1, components="ZZ", show=False, outfile="?.png")
        fn = "dvv ['ZZ']-f1-MM.png"
        self.assertTrue(os.path.isfile(fn),
                        msg="%s doesn't exist" % fn)


    # def test_103_plot_mwcs(self):
    #     from ..plots.mwcs import main
    #     main("YA.UV05", "YA.UV06", filterid=1, components="ZZ",
    #          mov_stack=5, show=False, outfile="?.png")
    #     fn = "mwcs YA_UV05_YA_UV06-ZZ-f1-m5.png"
    #     self.assertTrue(os.path.isfile(fn),
    #                     msg="%s doesn't exist" % fn)

    def test_104_plot_data_availability(self):
        import glob
        from ..plots.data_availability import main
        main(chan="HHZ", show=False, outfile="?.png")
        fn = glob.glob("data availability on*.png")
        self.assertEqual(len(fn), 1)

    def test_105_db_dump(self):
        """
        Tests the dump of the database and the creation of csv files
        """
        os.system("msnoise db dump")
        self.assertTrue(os.path.isfile("config.csv"))
        self.assertTrue(os.path.isfile("stations.csv"))
        self.assertTrue(os.path.isfile("filters.csv"))
        self.assertTrue(os.path.isfile("jobs.csv"))
        self.assertTrue(os.path.isfile("data_availability.csv"))

    def test_106_plot_wct(self):
        from ..plots.wct_dvv import main
        main(filterid=1, components="ZZ", show=False, outfile="?.png")
        fn = "wct ZZ-f1-dvv.png"
        self.assertTrue(os.path.isfile(fn),
                        msg="%s doesn't exist" % fn)
                        
    def test_107_plot_aurelogram(self):
        from ..plots.aurelogram import main
        main(filterid=1, components="ZZ", show=False, outfile="aurelogram.png")
        fn = "aurelogram.png"
        self.assertTrue(os.path.isfile(fn),
                        msg="%s doesn't exist" % fn)
  
    ### A few click CLI interface tests

    def test_201_config_get_unknown_param(self):
        result = self.runner.invoke(msnoise_script.config_get,
                                    ['inexistant_param'])
        self.assertEqual(result.exit_code, 0)
        self.assertIn('unknown parameter', result.output)

    def test_202_config_set_unknown_param(self):
        result = self.runner.invoke(msnoise_script.config_set,
                                    ['inexistant_param=value'])
        self.assertEqual(result.exit_code, 0)
        self.assertIn('unknown parameter', result.output)

    def test_203_config_set_param(self):
        result = self.runner.invoke(msnoise_script.config_set,
                                    ['channels=XXX'])
        self.assertEqual(result.exit_code, 0)
        result = self.runner.invoke(msnoise_script.config_get,
                                    ['channels'])
        self.assertEqual(result.exit_code, 0)
        self.assertIn('XXX', result.output)
        result = self.runner.invoke(msnoise_script.config_set,
                                    ['channels=*'])



    def test_301_compute_psd(self):
        from ..ppsd_compute import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_302_psd2hdf(self):
        from ..psd_to_hdf import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_303_hdf2rms(self):
        from ..psd_compute_rms import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_304_export_rms(self):
        from ..psd_export_rms import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()

    def test_400_run_manually(self):
        os.system("msnoise reset STACK --all")
        os.system("msnoise cc stack -m")
        os.system("msnoise cc dvv compute_mwcs")
        os.system("msnoise cc dvv compute_dtt")
        os.system("msnoise cc dvv compute_dvv")

    def test_99210_crondays_positive_float(self):
        """
        The crondays parameter can be a positive float that represents a
        number of days.
        """
        parsed_crondays = s01scan_archive.parse_crondays('2.5')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=2.5))

    def test_99211_crondays_negative_float(self):
        """
        A negative crondays parameter representing a number of days should be
        accepted for backward compatibility with MSNoise < 1.6.
        """
        parsed_crondays = s01scan_archive.parse_crondays('-3')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=3))

    def test_99212_crondays_weeks(self):
        """
        The crondays parameter can designate a number of weeks (7 days) using a
        string in the format 'Xw'.
        """
        parsed_crondays = s01scan_archive.parse_crondays('2w')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=7*2))

    def test_99213_crondays_days(self):
        """
        The crondays parameter can designate a number of days using a string in
        the format 'Xd'.
        """
        parsed_crondays = s01scan_archive.parse_crondays('5d')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=5))

    def test_99214_crondays_hours(self):
        """
        The crondays parameter can designate a number of hours using a string
        in the format 'Xh'.
        """
        parsed_crondays = s01scan_archive.parse_crondays('12h')
        self.assertEqual(parsed_crondays, datetime.timedelta(seconds=12*3600))

    def test_99215_crondays_weeks_days_hours(self):
        """
        The crondays parameter can designate a number of weeks, days, and hours
        in the same string.
        """
        parsed_crondays = s01scan_archive.parse_crondays('2w 3d 12h')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=2*7+3, seconds=12*3600))

    def test_99216_crondays_weeks_hours(self):
        """
        The crondays parameter does not have to designate the three units
        weeks, days and hours in the string.
        """
        parsed_crondays = s01scan_archive.parse_crondays('1w 6h')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=1*7, seconds=6*3600))

    def test_99217_crondays_weeks_days_hours_order_matters(self):
        """
        If the crondays parameter designates any weeks, days or hours in the
        string, they must be in the right order.
        """
        with self.assertRaises(FatalError):
            s01scan_archive.parse_crondays('16h 3d')

    def test_99218_crondays_weeks_days_hours_alone(self):
        """
        The crondays parameter can designate any weeks, days and hours, but the
        format must be [Xw][Xd][Xh] only.
        """
        with self.assertRaises(FatalError):
            s01scan_archive.parse_crondays('about 16h')

    def test_99219_crondays_weeks_days_hours_optional_blank(self):
        """
        If the crondays parameter designates any weeks, days and hours in
        the string, the separation blank is optional.
        """
        parsed_crondays = s01scan_archive.parse_crondays('3w4d12h')
        self.assertEqual(parsed_crondays, datetime.timedelta(days=3*7+4, seconds=12*3600))


def main(prefix="", tech=1):
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")

    import os
    import sys
    test_dir = tempfile.mkdtemp(prefix="msnoise_")
    os.chdir(test_dir)

    # c = len(os.listdir(os.getcwd()))
    # if c > 0:
    #     print("Directory is not empty, can't run tests here!")
    #     sys.exit()
    os.environ["PREFIX"] = prefix
    os.environ["hash"] = "h" + test_dir[-10:]
    os.environ["TECH"] = str(tech)
    print("Tests will be executed in %s" % test_dir)
    if tech == 2:
        print("The database localhost/%s was NOT deleted after this test!" % os.environ["HASH"])

    suite = unittest.defaultTestLoader.loadTestsFromTestCase(MSNoiseTests)
    runner = unittest.TextTestRunner(verbosity=4)
    result = runner.run(suite)

    print("Tests executed in %s" % test_dir)
    if tech == 2:
        print("The database localhost/%s was NOT deleted after this test!" % os.environ["HASH"])
    if not result.wasSuccessful():
        sys.exit(1)


if __name__ == '__main__':
    main()
