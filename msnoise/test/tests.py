import datetime
import glob
import logging
import os
import shutil
import traceback
from click.testing import CliRunner
from obspy import read
from numpy.testing import assert_allclose
from sqlalchemy import text
import tempfile
import pandas as pd
import pooch
import pytest
from .. import s01scan_archive, FatalError
from ..scripts import msnoise as msnoise_script
from ..api import (connect, get_config, update_config, get_job_types,
                   get_new_files, get_filters, get_station_pairs,
                   get_components_to_compute, update_filter, Filter,
                   get_stations, update_station, get_data_availability,
                   count_data_availability_flags, is_next_job, get_next_job,
                   Job, reset_jobs, build_ref_datelist, build_movstack_datelist,
                   read_db_inifile)
from ..s01scan_archive import parse_crondays
from ..s02new_jobs import main as new_jobs_main
from ..s03compute_no_rotation import main as compute_cc_main
from ..s04_stack2 import main as stack_main
from ..s05compute_mwcs2 import main as compute_mwcs_main
from ..s06compute_dtt2 import main as compute_dtt_main
from ..s07_compute_dvv import main as compute_dvv_main
from ..psd_compute_rms import main as compute_rms_main
from ..psd_export_rms import main as export_rms_main
from ..ppsd_compute import main as ppsd_compute_main
from ..psd_to_hdf import main as psd_to_hdf_main

from ..plots.ccftime import main as ccftime_main
from ..plots.interferogram import main as interferogram_main
from ..plots.spectime import main as spectime_main
from ..plots.distance import main as distance_main
from ..plots.dvv import main as dvv_main
from ..plots.data_availability import main as data_availability_main
from ..plots.wct_dvv import main as wct_dvv_main

global logger
logger = logging.getLogger('matplotlib')
# set WARNING for Matplotlib
logger.setLevel(logging.CRITICAL)

@pytest.fixture(scope="session", autouse=True)
def setup_environment():
    # Switch Matplotlib backend
    import matplotlib.pyplot as plt
    plt.switch_backend("agg")

    # Create a temporary directory for tests
    test_dir = tempfile.mkdtemp(prefix="msnoise_")
    os.chdir(test_dir)

    # Set environment variables
    os.environ["PREFIX"] = ""
    os.environ["hash"] = "h" + test_dir[-10:]
    os.environ["TECH"] = "1"

    # Define test files and use pooch to handle fetching
    test_files = {
        "classic/data/2010/UV05/HHZ.D/YA.UV05.00.HHZ.D.2010.244": "17034091285d485f7c2d4797f435228c408d6940db943be63f1769ec09854f4f",
        "classic/data/2010/UV06/HHZ.D/YA.UV06.00.HHZ.D.2010.244": "51bfd1e735696e83ee6dba136c9e740c59120fac9f74b386eac75062eb9ca382",
        "classic/data/2010/UV10/HHZ.D/YA.UV10.00.HHZ.D.2010.244": "530cc7f4a57fe69a8a5cedeb18e64773055c146e4ae4676012f6618dd0c92e82",
        "classic/extra/DATA.RESIF_Jun_10,14_21_05_20264.RESIF": "95a6d007132fc41b6107d258aeee1170614d234cdd3eb4a6d5652e4661a6adcd",
        "classic/extra/stations.csv": "057152c2823c5457bce879146d78984af422973ab313e1d1cd8baaa7f7a1d6b3",
        "classic/extra/test_inventory.xml": "48bf3261a9c23f1782e452583306b73643a4d7f205df4338e5525df3c06eccb9",
    }

    BRIAN = pooch.create(
        path=pooch.os_cache("msnoise-testdata"),
        base_url="https://github.com/ROBelgium/msnoise-testdata/raw/{version}/",
        version="1.1",
        version_dev="main",
        registry=test_files,
        env="MSNOISE_DATA_DIR",
    )

    # Fetch all files
    for fn in test_files:
        BRIAN.fetch(fn)

    # Copy data folder
    data_folder = os.path.join(os.environ.get("MSNOISE_DATA_DIR", pooch.os_cache("msnoise-testdata")), "1.1", "classic", "data")
    response_path = os.path.join(os.environ.get("MSNOISE_DATA_DIR", pooch.os_cache("msnoise-testdata")), "1.1", "classic", "extra")

    if not os.path.isdir("data"):
        shutil.copytree(data_folder, "data/")

    runner = CliRunner()

    yield {
        "data_folder": "data",
        "response_path": response_path,
        "runner": runner
    }
@pytest.mark.order(1)
def test_001_S01installer(setup_environment):
    from ..s000installer import main
    prefix = os.environ.get("PREFIX", "")

    try:
        tech = int(os.environ.get("TECH", "1"))
        if tech == 1:
            ret = main(tech=1, prefix=prefix, hostname="localhost", username="root", password="SECRET",
                       database=os.environ["hash"])
        elif tech == 2:
            ret = main(tech=2, prefix=prefix, hostname=os.environ["MARIADB_HOSTNAME"],
                       username=os.environ["MARIADB_USERNAME"], password=os.environ["MARIADB_PASSWORD"],
                       database=os.environ["hash"])
        assert ret == 0
    except:
        traceback.print_exc()
        pytest.fail()


@pytest.mark.order(2)
def test_002_ConnectToDB():
    try:
        db = connect()
        db.close()
    except:
        pytest.fail("Can't connect to MSNoise DB")


@pytest.mark.order(3)
def test_003_set_and_config(setup_environment):
    db = connect()
    data_folder = setup_environment['data_folder']
    response_path = setup_environment['response_path']
    totests = [
        ['data_folder', data_folder],
        ['data_structure', 'PDF'],
        ['network', 'YA'],
        ['components_to_compute', 'ZZ'],
        ['response_path', response_path]
    ]
    for key, value in totests:
        update_config(db, key, value)
        config_value = get_config(db, key)
        assert config_value == value, f"Configuration parameter {key} did not set correctly."
    db.close()

@pytest.mark.order(4)
def test_004_set_and_get_filters():
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
        update_filter(db, f.ref, f.low, f.mwcs_low, f.high, f.mwcs_high, f.mwcs_wlen, f.mwcs_step, f.used)
    dbfilters = get_filters(db)
    for i, filter in enumerate(dbfilters):
        for param in ['low', 'mwcs_low', 'high', 'mwcs_high', 'mwcs_wlen', 'mwcs_step', 'used']:
            assert eval(f"filter.{param}") == eval(f"filters[i].{param}")

@pytest.mark.order(5)
def test_005_populate_station_table():
    from ..s002populate_station_table import main
    try:
        ret = main()
        assert ret is True
    except:
        pytest.fail()

@pytest.mark.order(6)
def test_006_get_stations():
    db = connect()
    stations = get_stations(db).all()
    assert len(stations) == 3
    db.close()

@pytest.mark.order(7)
def test_007_update_stations(setup_environment):
    db = connect()
    stations_file = os.path.join(setup_environment['response_path'], 'stations.csv')
    stations_df = pd.read_csv(stations_file, header=None, index_col=0, names=['X', 'Y', 'altitude'])
    for station in get_stations(db):
        fullname = f"{station.net}.{station.sta}"
        try:
            s = stations_df.loc[fullname]
            update_station(db, net=station.net, sta=station.sta, X=s['X'], Y=s['Y'], altitude=s['altitude'])
        except:
            traceback.print_exc()
            pytest.fail()

    runner = setup_environment['runner']
    result = runner.invoke(msnoise_script.config_sync)
    assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}"

    result = runner.invoke(msnoise_script.info)
    assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}"

@pytest.mark.order(8)
def test_008_scan_archive(setup_environment):
    from ..s01scan_archive import main
    try:
        main(init=True, threads=1)
    except:
        traceback.print_exc()
        pytest.fail()


    runner = setup_environment['runner']
    result = runner.invoke(msnoise_script.db_da_stations_update_loc_chan)
    assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}"
    # Optionally, add more assertions to check specific outputs in result.output


@pytest.mark.order(9)
def test_009_control_data_availability():
    db = connect()
    files = get_new_files(db)
    assert len(files) == 3
    flags = count_data_availability_flags(db)
    assert len(flags) == 1
    for station in get_stations(db):
        for loc in station.locs():
            for chan in station.chans():
                da = get_data_availability(db, net=station.net, sta=station.sta, loc=loc, chan=chan)
                assert len(da) == 1

@pytest.mark.order(10)
def test_010_new_jobs():
    try:
        new_jobs_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(11)
def test_011_control_jobs():
    db = connect()
    assert is_next_job(db) is True
    jobs = get_next_job(db)
    assert isinstance(jobs[0], Job)

@pytest.mark.order(12)
def test_012_reset_jobs():
    db = connect()
    reset_jobs(db, 'CC', alljobs=True)
    db.close()

@pytest.mark.order(13)
def test_013_s03compute_cc():
    try:
        compute_cc_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(14)
def test_014_check_done_jobs():
    db = connect()
    jobs = get_job_types(db, 'CC')
    assert jobs[0][0] == 3
    assert jobs[0][1] == 'D'

@pytest.mark.order(15)
def test_015_check_cc_files(format="MSEED"):
    db = connect()
    for filter in get_filters(db):
        for components in get_components_to_compute(db):
            for (sta1, sta2) in get_station_pairs(db):
                for loc1 in sta1.locs():
                    for loc2 in sta2.locs():
                        pair = f"{sta1.net}.{sta1.sta}.{loc1}_{sta2.net}.{sta2.sta}.{loc2}"
                        tmp = os.path.join("STACKS", f"{filter.ref:02}", "001_DAYS", components, pair, f"2010-09-01.{format}")
                        print("checking", tmp)
                        assert os.path.isfile(tmp), f"{tmp} does not exist"

@pytest.mark.order(16)
def test_016_update_config():
    db = connect()
    update_config(db, "export_format", "SAC")
    db.close()

@pytest.mark.order(17)
def test_017_reset_cc_jobs():
    db = connect()
    reset_jobs(db, 'CC', alljobs=True)
    db.close()

@pytest.mark.order(18)
def test_018_recompute_cc():
    test_013_s03compute_cc()

@pytest.mark.order(19)
def test_019_check_SACS():
    test_015_check_cc_files(format='SAC')


@pytest.mark.order(20)
def test_020_update_config():
    db = connect()
    shutil.rmtree("STACKS")
    update_config(db, "export_format", "BOTH")
    db.close()

@pytest.mark.order(21)
def test_021_reprocess_BOTH():
    test_017_reset_cc_jobs()
    test_013_s03compute_cc()
    test_015_check_cc_files(format='MSEED')
    test_015_check_cc_files(format='SAC')

@pytest.mark.order(22)
def test_022_check_content():
    db = connect()
    for filter in get_filters(db):
        for components in get_components_to_compute(db):
            for (sta1, sta2) in get_station_pairs(db):
                for loc1 in sta1.locs():
                    for loc2 in sta2.locs():
                        sta1_id = f"{sta1.net}.{sta1.sta}.{loc1}"
                        sta2_id = f"{sta2.net}.{sta2.sta}.{loc2}"
                        pair = f"{sta1_id}_{sta2_id}"
                        tmp1 = os.path.join("STACKS", f"{filter.ref:02}", "001_DAYS", components, pair, "2010-09-01.MSEED")
                        tmp2 = os.path.join("STACKS", f"{filter.ref:02}", "001_DAYS", components, pair, "2010-09-01.SAC")
                        data1 = read(tmp1)[0].data
                        data2 = read(tmp2)[0].data
                        assert_allclose(data1, data2)
    db.close()

@pytest.mark.order(23)
def test_023_stack():
    db = connect()
    update_config(db, 'ref_begin', '2009-01-01')
    update_config(db, 'ref_end', '2011-01-01')
    update_config(db, 'startdate', '2009-01-01')
    update_config(db, 'enddate', '2011-01-01')
    update_config(db, 'mov_stack', "(('1d','1d'),('2d','1d'),('5d','1d'))")
    interval = 1.
    stack_main('ref', interval)
    reset_jobs(db, "STACK", alljobs=True)
    stack_main('mov', interval)
    reset_jobs(db, "STACK", alljobs=True)
    stack_main('step', interval)
    
    update_config(db, 'wienerfilt', 'Y')
    reset_jobs(db, "STACK", alljobs=True)
    stack_main('mov', interval)
    stack_main('ref', interval)

    db.close()

@pytest.mark.order(24)
def test_024_mwcs():
    compute_mwcs_main()

@pytest.mark.order(25)
def test_025_dtt():
    compute_dtt_main()

@pytest.mark.order(26)
def test_026_dvv():
    compute_dvv_main()

@pytest.mark.order(27)
def test_027_build_ref_datelist():
    from ..api import build_ref_datelist
    db = connect()
    start, end, datelist = build_ref_datelist(db)
    assert start == datetime.date(2009, 1, 1)
    assert end == datetime.date(2011, 1, 1)
    assert len(datelist) == 731
    db.close()

@pytest.mark.order(28)
def test_028_build_movstack_datelist():
    from ..api import build_movstack_datelist
    db = connect()
    start, end, datelist = build_movstack_datelist(db)
    assert start == datetime.date(2009, 1, 1)
    assert end == datetime.date(2011, 1, 1)
    assert len(datelist) == 731
    db.close()

@pytest.mark.order(29)
def test_029_stretching():
    from ..stretch2 import main as stretch_main
    db = connect()
    update_config(db, "export_format", "MSEED")
    # reset_jobs(db, "MWCS", alljobs=True)
    db.close()
    stretch_main()

@pytest.mark.order(30)
def test_030_create_fake_new_files(setup_environment):
    data_folder = setup_environment['data_folder']
    for f in sorted(glob.glob(os.path.join(data_folder, "2010", "*", "HHZ.D", "*"))):
        st = read(f)
        for tr in st:
            tr.stats.starttime += datetime.timedelta(days=1)
        out = f.replace("244", "245")
        st.write(out, format="MSEED")

    from ..s01scan_archive import main
    try:
        main(init=False, threads=1)
    except:
        traceback.print_exc()
        pytest.fail()

    new_jobs_main()
    db = connect()
    jobs = get_job_types(db, 'CC')
    assert jobs[0][0] == 3
    assert jobs[0][1] == 'D'
    assert jobs[1][0] == 3
    assert jobs[1][1] == 'T'

@pytest.mark.order(31)
def test_031_instrument_response(setup_environment):
    db = connect()
    response_path = setup_environment['response_path']
    update_config(db, 'response_path', response_path)
    update_config(db, 'response_format', "dataless")
    update_config(db, 'remove_response', "Y")
    db.close()
    test_013_s03compute_cc()

@pytest.mark.order(32)
def test_032_wct():
    from ..s08compute_wct import main as compute_wct_main
    compute_wct_main()
  
@pytest.mark.order(33)
def test_033_validate_stack_data():
    from ..api import validate_stack_data
    import xarray as xr
    import numpy as np
    import pandas as pd
    
    # Test empty dataset
    ds = xr.Dataset()
    is_valid, message = validate_stack_data(ds, "reference")
    assert not is_valid
    assert "No data found for reference stack" in message
    
    # Test dataset without CCF
    ds = xr.Dataset({"wrong_var": 1})
    is_valid, message = validate_stack_data(ds, "reference") 
    assert not is_valid
    assert "Missing CCF data in reference stack" in message

    # Test empty CCF data
    times = pd.date_range('2020-01-01', periods=0)
    taxis = np.linspace(-50, 50, 100)
    data = np.random.random((0, len(taxis)))
    da = xr.DataArray(data, coords=[times, taxis], dims=['times', 'taxis'])
    ds = da.to_dataset(name='CCF')
    is_valid, message = validate_stack_data(ds, "reference")
    assert not is_valid 
    assert "Empty dataset in reference stack" in message

    # Test all NaN values
    times = pd.date_range('2020-01-01', periods=10)
    data = np.full((len(times), len(taxis)), np.nan)
    da = xr.DataArray(data, coords=[times, taxis], dims=['times', 'taxis'])
    ds = da.to_dataset(name='CCF')
    is_valid, message = validate_stack_data(ds, "reference")
    assert not is_valid
    assert "Reference stack contains only NaN values" in message

    # Test partial NaN values
    data = np.random.random((len(times), len(taxis)))
    data[0:5, :] = np.nan
    da = xr.DataArray(data, coords=[times, taxis], dims=['times', 'taxis'])
    ds = da.to_dataset(name='CCF')
    is_valid, message = validate_stack_data(ds, "reference")
    assert is_valid
    # We can still check if the warning message contains the percentage
    assert "50.0% NaN values" in message

    # Test valid data
    data = np.random.random((len(times), len(taxis)))
    da = xr.DataArray(data, coords=[times, taxis], dims=['times', 'taxis'])
    ds = da.to_dataset(name='CCF')
    is_valid, message = validate_stack_data(ds, "reference")
    assert is_valid
    assert message == "OK"
    
@pytest.mark.order(34)
def test_034_stack_validation_handling():
    from ..api import validate_stack_data
    import xarray as xr
    import numpy as np
    import pandas as pd
    
    # Create minimal test data
    times = pd.date_range('2020-01-01', periods=10)
    taxis = np.linspace(-50, 50, 100)
    
    # Test with actual code's variables and logic
    pairs = [('STA1', 'STA2'), ('STA3', 'STA3')]
    filters = [type('Filter', (), {'ref': '1'})]
    components = 'ZZ'
    
    for sta1, sta2 in pairs:
        for f in filters:
            filterid = int(f.ref)
            
            # Create test datasets that will trigger our paths
            if sta1 == 'STA1':
                # Test error path with empty dataset
                c = xr.Dataset()
                is_valid, message = validate_stack_data(c, "reference")
                if not is_valid:
                    logger.error(f"Invalid reference data for {sta1}:{sta2}-{components}-{filterid}: {message}")
                    continue
            else:
                # Test warning path with partial NaN data
                data = np.random.random((len(times), len(taxis)))
                data[0:5, :] = np.nan
                da = xr.DataArray(data, coords=[times, taxis], dims=['times', 'taxis'])
                c = da.to_dataset(name='CCF')
                
                is_valid, message = validate_stack_data(c, "reference")
                if "Warning" in message:
                    logger.warning(f"{sta1}:{sta2}-{components}-{filterid}: {message}")
                  
@pytest.mark.order(100)
def test_100_plot_interferogram():
    db = connect()
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                sta_id1 = f"{sta1.net}.{sta1.sta}.{loc1}"
                sta_id2 = f"{sta2.net}.{sta2.sta}.{loc2}"
                for filter in get_filters(db):
                    interferogram_main(sta_id1, sta_id2, filter.ref, "ZZ", 1, show=False, outfile="?.png")
                    fn = f'interferogram {sta_id1}-{sta_id2}-ZZ-f{filter.ref}-m1d_1d.png'
                    assert os.path.isfile(fn), f"{fn} doesn't exist"

@pytest.mark.order(101)
def test_101_plot_spectime():
    db = connect()
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                sta_id1 = f"{sta1.net}.{sta1.sta}.{loc1}"
                sta_id2 = f"{sta2.net}.{sta2.sta}.{loc2}"
                for filter in get_filters(db):
                    spectime_main(sta_id1, sta_id2, filter.ref, "ZZ", 1, show=False, outfile="?.png")
                    fn = f'spectime {sta_id1}-{sta_id2}-ZZ-f{filter.ref}-m1d_1d.png'
                    assert os.path.isfile(fn), f"{fn} doesn't exist"

@pytest.mark.order(102)
def test_102_plot_distance():
    distance_main(filterid=1, components="ZZ", show=False, outfile="?.png")
    fn = "distance ZZ-f1.png"
    assert os.path.isfile(fn), f"{fn} doesn't exist"

    distance_main(filterid=1, components="ZZ", show=False, outfile="?_refilter.png", refilter="0.2:0.9")
    fn = "distance ZZ-f1_refilter.png"
    assert os.path.isfile(fn), f"{fn} doesn't exist"

@pytest.mark.order(103)
def test_103_plot_dvv():
    dvv_main(filterid=1, components="ZZ", show=False, outfile="?.png")
    fn = "dvv ['ZZ']-f1-MM.png"
    assert os.path.isfile(fn), f"{fn} doesn't exist"

@pytest.mark.order(104)
def test_104_plot_data_availability():
    data_availability_main(chan="HHZ", show=False, outfile="?.png")
    fn = glob.glob("data availability on*.png")
    assert len(fn) == 1, "Data availability plot doesn't exist"

@pytest.mark.order(105)
def test_105_db_dump():
    """ Tests the dump of the database and the creation of csv files """
    os.system("msnoise db dump")
    assert os.path.isfile("config.csv")
    assert os.path.isfile("stations.csv")
    assert os.path.isfile("filters.csv")
    assert os.path.isfile("jobs.csv")
    assert os.path.isfile("data_availability.csv")

    # os.system("msnoise db import config --force")

@pytest.mark.order(106)
def test_106_plot_wct():
    wct_dvv_main(filterid=1, components="ZZ", show=False, outfile="?.png")
    fn = "wct ZZ-f1-dvv.png"
    assert os.path.isfile(fn), f"{fn} doesn't exist"


@pytest.mark.order(201)
def test_201_config_get_unknown_param(setup_environment):
    runner = setup_environment['runner']
    result = runner.invoke(msnoise_script.config_get, ['inexistant_param'])
    assert result.exit_code == 0
    assert 'unknown parameter' in result.output


@pytest.mark.order(202)
def test_202_config_set_unknown_param(setup_environment):
    runner = setup_environment['runner']
    result = runner.invoke(msnoise_script.config_set, ['inexistant_param=value'])
    assert result.exit_code == 0
    assert 'unknown parameter' in result.output


@pytest.mark.order(203)
def test_203_config_set_param(setup_environment):
    runner = setup_environment['runner']

    result = runner.invoke(msnoise_script.config_set, ['channels=XXX'])
    assert result.exit_code == 0

    result = runner.invoke(msnoise_script.config_get, ['channels'])
    assert result.exit_code == 0
    assert 'XXX' in result.output

    result = runner.invoke(msnoise_script.config_set, ['channels=*'])
    assert result.exit_code == 0

@pytest.mark.order(301)
def test_301_compute_psd():
    try:
        ppsd_compute_main()
    except:
        traceback.print_exc()
        pytest.fail("PSD Compute failed")

@pytest.mark.order(302)
def test_302_psd2hdf():
    try:
        psd_to_hdf_main()
    except:
        traceback.print_exc()
        pytest.fail("PSD to HDF conversion failed")

@pytest.mark.order(303)
def test_303_hdf2rms():
    try:
        compute_rms_main()
    except:
        traceback.print_exc()
        pytest.fail("HDF to RMS computation failed")

@pytest.mark.order(304)
def test_304_export_rms():
    try:
        export_rms_main()
    except:
        traceback.print_exc()
        pytest.fail("RMS export failed")

@pytest.mark.order(400)
def test_400_run_manually():
    os.system("msnoise reset STACK --all")
    os.system("msnoise cc stack -m")
    os.system("msnoise cc dvv compute_mwcs")
    os.system("msnoise cc dvv compute_dtt")
    os.system("msnoise cc dvv compute_dvv")
    os.system("msnoise cc dvv compute_stretching")
    os.system("msnoise cc dvv compute_wct")

def test_99210_crondays_positive_float():
    parsed_crondays = parse_crondays('2.5')
    assert parsed_crondays == datetime.timedelta(days=2.5)

def test_99211_crondays_negative_float():
    parsed_crondays = parse_crondays('-3')
    assert parsed_crondays == datetime.timedelta(days=3)

def test_99212_crondays_weeks():
    parsed_crondays = parse_crondays('2w')
    assert parsed_crondays == datetime.timedelta(days=7*2)

def test_99213_crondays_days():
    parsed_crondays = parse_crondays('5d')
    assert parsed_crondays == datetime.timedelta(days=5)

def test_99214_crondays_hours():
    parsed_crondays = parse_crondays('12h')
    assert parsed_crondays == datetime.timedelta(seconds=12*3600)

def test_99215_crondays_weeks_days_hours():
    parsed_crondays = parse_crondays('2w 3d 12h')
    assert parsed_crondays == datetime.timedelta(days=2*7+3, seconds=12*3600)

def test_99216_crondays_weeks_hours():
    parsed_crondays = parse_crondays('1w 6h')
    assert parsed_crondays == datetime.timedelta(days=1*7, seconds=6*3600)

def test_99217_crondays_weeks_days_hours_order_matters():
    with pytest.raises(FatalError):
        parse_crondays('16h 3d')

def test_99218_crondays_weeks_days_hours_alone():
    with pytest.raises(FatalError):
        parse_crondays('about 16h')

def test_99219_crondays_weeks_days_hours_optional_blank():
    parsed_crondays = parse_crondays('3w4d12h')
    assert parsed_crondays == datetime.timedelta(days=3*7+4, seconds=12*3600)

@pytest.mark.order(100000)
def test_100000_msnoise_admin():
    from ..msnoise_admin import get_app
    """
    GIVEN a Flask application configured for testing
    WHEN the '/' page is requested (GET)
    THEN check that the response is valid
    """
    # Set the Testing configuration prior to creating the Flask application
    os.environ['CONFIG_TYPE'] = 'config.TestingConfig'
    flask_app = get_app()

    # Create a test client using the Flask application configured for testing
    with flask_app.test_client() as test_client:
        response = test_client.get('admin/')
        assert response.status_code == 200, f"Error following route {route}"
        assert b"MSNoise Dashboard" in response.data
        for route in ["admin/config/","admin/stations/", "admin/filters/",
                      "admin/data_availability/", "admin/jobs/",
                      "admin/bugreport/"]:
            response = test_client.get(route)
            assert response.status_code == 200, f"Error following route {route}"

        route = "admin/stations/new/?url=/admin/stations/"
        response = test_client.get(route)
        assert response.status_code == 200, f"Error following route {route}"
        #
        for route in ["admin/filters/edit/?id=1",
                      "admin/stations/edit/?id=1&",
                      "admin/config/edit/?id=network",
                      "admin/data_availability/edit/?id=3",
                      "admin/jobs/edit/?id=9&url=/admin/jobs/&modal=True"]:
            response = test_client.get(route, follow_redirects=True)
            assert response.status_code == 200, f"Error following route {route}"

        for route in ["admin/pairs.json",
                      "admin/bugreport.json",
                      "admin/data_availability_flags.json",
                      ]:
            response = test_client.get(route, follow_redirects=True)
            assert response.status_code == 200, f"Error following route {route}"


@pytest.mark.order(200000)
def test_20000_invoke_script(setup_environment):
    runner = setup_environment['runner']

    for cmd in [
        msnoise_script.config_sync,
        msnoise_script.db_upgrade,
        msnoise_script.db_clean_duplicates,
        ]:

        result = runner.invoke(cmd)
        assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}"
