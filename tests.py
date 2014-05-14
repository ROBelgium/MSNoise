import unittest
import traceback

# Here's our "unit".
def IsOdd(n):
    return n % 2 == 1

# Here's our "unit tests".
class InitTests(unittest.TestCase):

    
    def test_001_S01installer(self):
        from s000installer import main
        try:
            ret = main(tech=1)
            msg = "Installation Done! - Go to Configuration Step!"
            self.failUnlessEqual(ret, msg)
        except:
            traceback.print_exc()
            self.fail()
    
    def test_002_ConnectToDB(self):
        from database_tools import connect
        try:
            db = connect()
            db.close()
        except:
            self.fail("Can't connect to MSNoise DB")
    
    def test_003_set_and_config(self):
        from database_tools import connect, get_config, update_config
        db = connect()
        totests = []
        totests.append(['data_folder', 'data'])
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
            update_config(db,test[0], test[1])
            d = get_config(db, test[0])
            self.failUnlessEqual(test[1], d)

        db.close()

    def test_004_set_and_get_filters(self):
        from msnoise_table_def import Filter
        from database_tools import connect, update_filter, get_filters
        db = connect()
        filters = []
        filters.append(Filter(0.01, 0.12, 0.98, 1.0, 0, 10, 5, 1))
        filters.append(Filter(0.1, 0.12, 0.98, 1.0, 0, 10, 5, 1))
        
        for f in filters:
            update_filter(db, f.ref, f.low, f.high, f.mwcs_low, f.mwcs_high,
                  f.rms_threshold, f.mwcs_wlen, f.mwcs_step, f.used)
        
        dbfilters = get_filters(db)
        
        for i, filter in enumerate(dbfilters):
            self.failUnlessEqual(filter.low, filters[i].low)
    
    def test_005_populate_station_table(self):
        from s002populate_station_table import main
        try:
            ret = main()
            self.failUnlessEqual(ret, True)
        except:
            self.fail()
    
    def test_006_get_stations(self):
        from database_tools import connect, get_stations
        db = connect()
        stations = get_stations(db).all()
        self.failUnlessEqual(len(stations), 3)
        db.close()
        
    def test_007_update_stations(self):
        from database_tools import connect, get_stations, update_station
        import pandas as pd
        db = connect()
        stations = pd.read_csv('extra/stations.csv',header=None, index_col = 0, names =['X','Y','altitude'])
        for station in get_stations(db):
            fullname = "%s.%s" % (station.net, station.sta,)
            try:
                s = stations.loc[fullname]
                update_station(db, station.net, station.sta, s['X'], s['Y'], s['altitude'])
            except:
                traceback.print_exc()
                self.fail()
                continue
            del s
        
    def test_008_scan_archive(self):
        from s01scan_archive import main
        try:
            main(init=True,threads=1)
        except:
            traceback.print_exc()
            self.fail()
    
    def test_009_control_data_availability(self):
        from database_tools import connect, get_new_files, get_data_availability, count_data_availability_flags, get_stations
        
        db = connect()
        files = get_new_files(db)
        self.failUnlessEqual(len(files), 3)
        
        flags = count_data_availability_flags(db)
        self.failUnlessEqual(len(flags), 1)
        
        for station in get_stations(db):
            da = get_data_availability(db, net=station.net, sta=station.sta, comp='HHZ')
            self.failUnlessEqual(len(da), 1)
        
    
    def test_010_new_jobs(self):
        from s02new_jobs import new_jobs
        
        try:
            new_jobs()
        except:
            traceback.print_exc()
            self.fail()
    
    def test_011_control_jobs(self):
        from database_tools import connect, is_next_job, get_next_job
        from msnoise_table_def import Job
        db = connect()
        
        self.failUnlessEqual(is_next_job(db), True)
        jobs = get_next_job(db)
        self.failUnlessEqual(isinstance(jobs[0], Job), True)
    
    def test_012_s03scan_archive(self):
        from s03compute_cc import main
        try:
            main()
        except:
            traceback.print_exc()
            self.fail()
        
    
    
def main():
    unittest.main()

if __name__ == '__main__':
    main()