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


    test_files = {
        "3C/data/2014/PF/CSS/HHE.D/PF.CSS.00.HHE.D.2014.171": "1b69947c3481db5145864dcba33bf9f846ceff50ede2b0a66df00696a5d589f1",
        "3C/data/2014/PF/CSS/HHE.D/PF.CSS.00.HHE.D.2014.172": "4b1d1095337f9294180aad9b1ba3fdee40a472b34360ebd8b2f4b4e3b474754e",
        "3C/data/2014/PF/CSS/HHE.D/PF.CSS.00.HHE.D.2014.173": "cc8568adca6a6bdbc4d12fd8e6c5cfdaf10a0a7c98743262de4544114b30674e",
        "3C/data/2014/PF/CSS/HHN.D/PF.CSS.00.HHN.D.2014.171": "a0bf6049eb9608b0207308ed02ad555fad2071a1ca2a5415196d084af88d13b8",
        "3C/data/2014/PF/CSS/HHN.D/PF.CSS.00.HHN.D.2014.172": "d263f8eb19f4aa658eda893206c2affe8c665986e810fec20d7298d00c742661",
        "3C/data/2014/PF/CSS/HHN.D/PF.CSS.00.HHN.D.2014.173": "5bd1a4836933fa65d832fde70a05ab9f58d0cc763bd01e9cca4e3d059531a32e",
        "3C/data/2014/PF/CSS/HHZ.D/PF.CSS.00.HHZ.D.2014.171": "de8a962960517f1a3cff9142104c764a4811d9cd3a13754a2af3550fbc50d25c",
        "3C/data/2014/PF/CSS/HHZ.D/PF.CSS.00.HHZ.D.2014.172": "5a4b3da58b63264d5b8458b2e7730c6e15a4594903230e337b0a9064a24f7653",
        "3C/data/2014/PF/CSS/HHZ.D/PF.CSS.00.HHZ.D.2014.173": "872c7d54bedd34090b58812e39ff6404e537eb67f45cbc55c22dc1515fe49e53",
        "3C/data/2014/PF/FJS/HHE.D/PF.FJS.00.HHE.D.2014.171": "632273b8ae6f9f4a3905f3fc80e35a0d76b875aa04d47ce232ab9a2370fdf41b",
        "3C/data/2014/PF/FJS/HHE.D/PF.FJS.00.HHE.D.2014.172": "7177a0df3527282a9c0b701464799210ea9a58bbc9e9062d5da6297a4a25106b",
        "3C/data/2014/PF/FJS/HHE.D/PF.FJS.00.HHE.D.2014.173": "c5e1e5281fe6a9b6c3dead4d59ac0188b6b09341390cb5f1455cb74541b39dd5",
        "3C/data/2014/PF/FJS/HHN.D/PF.FJS.00.HHN.D.2014.171": "18f4de5676d9a2dd80dcf6a13b2d5e1c97d6d8901ca13fe9f623f1c2854e367b",
        "3C/data/2014/PF/FJS/HHN.D/PF.FJS.00.HHN.D.2014.172": "c1fcf39d2525d75cf53351e9d078a963a1680fd0d06eab01d7ada90a384b9c8e",
        "3C/data/2014/PF/FJS/HHN.D/PF.FJS.00.HHN.D.2014.173": "d185706cd11e55b3da1e9521204662da6a4f1b7f81cfed6f26e4ad626caad7e1",
        "3C/data/2014/PF/FJS/HHZ.D/PF.FJS.00.HHZ.D.2014.171": "47f4824f36601b2a569d75054390d539c58bf65ff39749bb9d27cc89991e6925",
        "3C/data/2014/PF/FJS/HHZ.D/PF.FJS.00.HHZ.D.2014.172": "1d91a2d1b176874d1f8cc6cf6087efe1bd743f93f16e9196cdb10c6d2befcdc0",
        "3C/data/2014/PF/FJS/HHZ.D/PF.FJS.00.HHZ.D.2014.173": "6bb3edb2be3b9a5907659345963464f3a22849be0b022475e57596b719f53a97",
        "3C/data/2014/PF/FOR/HHE.D/PF.FOR.00.HHE.D.2014.171": "ca32402b36c0ec04199e6fb317a64378a71a99b71c72c091ffc47eec0da6f7aa",
        "3C/data/2014/PF/FOR/HHE.D/PF.FOR.00.HHE.D.2014.172": "aed64ca298be9421887cbdd934078a71ced5806263f72f9f06eabf59aae6b18b",
        "3C/data/2014/PF/FOR/HHE.D/PF.FOR.00.HHE.D.2014.173": "f29e3e09f3b2b15215f976ef4dea35379b87a01a46dd335e46f29a40a29d2601",
        "3C/data/2014/PF/FOR/HHN.D/PF.FOR.00.HHN.D.2014.171": "553b67513dcc3a48a746adff68c59952d7af3b347c2bd6f91d62987c18ed05d6",
        "3C/data/2014/PF/FOR/HHN.D/PF.FOR.00.HHN.D.2014.172": "02bd07cf33bf8f778361de21f9354a1c3f8c51fb488b2fbefb34d53e0a33b312",
        "3C/data/2014/PF/FOR/HHN.D/PF.FOR.00.HHN.D.2014.173": "b4e9fe2719318f777414b7cc743da0e16e4abdb77c0d5b40ac0368516fbbb11a",
        "3C/data/2014/PF/FOR/HHZ.D/PF.FOR.00.HHZ.D.2014.171": "7d09f218e287b9fe5659e23ff10dfa50669ad10d4162ac804d56e528e61eb3b1",
        "3C/data/2014/PF/FOR/HHZ.D/PF.FOR.00.HHZ.D.2014.172": "d3707d62eb45bba1400b38c808aaf8a90477cb8be036eee77aadfb613afbadd9",
        "3C/data/2014/PF/FOR/HHZ.D/PF.FOR.00.HHZ.D.2014.173": "3e15cc7b49c0315fea44e52b848db1dffd65dca6a879cec4e0a2fe8ce8a90f1d",
        "3C/resp/PF.CSS.xml":"992dcdde5937c40f956d78b2264fb5b11b19dbe30ecf398434aba9af836db540",
        "3C/resp/PF.FJS.xml":"13d12cb5b46c9bce3b8aad08e89e6216c161d21f84ea213d7bec48e03cb39eca",
        "3C/resp/PF.FOR.xml":"257474d734adb13a8d64878119aad14c0e3121a4054a2f8e52dbf81983ec5284",
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
    data_folder = os.path.join(os.environ.get("MSNOISE_DATA_DIR", pooch.os_cache("msnoise-testdata")), "1.1", "3C", "data")
    response_path = os.path.join(os.environ.get("MSNOISE_DATA_DIR", pooch.os_cache("msnoise-testdata")), "1.1", "3C", "resp")

    # Create a temporary directory for tests
    test_dir = tempfile.mkdtemp(prefix="msnoise_")
    os.chdir(test_dir)

    # Set environment variables
    os.environ["PREFIX"] = ""
    os.environ["hash"] = "h" + test_dir[-10:]
    os.environ["TECH"] = "1"


    runner = CliRunner()

    yield {
        "data_folder": data_folder,
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

@pytest.mark.order(3)
def test_003_set_and_config(setup_environment):
    db = connect()
    data_folder = setup_environment['data_folder']
    response_path = setup_environment['response_path']
    totests = [
        ['data_folder', data_folder],
        ['data_structure', "SDS"],
        ['components_to_compute', 'ZZ,EN,NZ'],
        ['components_to_compute_single_station', 'ZZ,EN,EZ'],
        ['response_path', response_path],
        ['corr_duration', '7200'],
        ['overlap', '0']
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
    f.low = 1
    f.mwcs_low = 1
    f.high = 5
    f.mwcs_high = 5
    f.mwcs_wlen = 5
    f.mwcs_step = 2
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
    assert len(files) == 27
    flags = count_data_availability_flags(db)
    assert len(flags) == 1
    for station in get_stations(db):
        for loc in station.locs():
            for chan in station.chans():
                da = get_data_availability(db, net=station.net, sta=station.sta, loc=loc, chan=chan)
                assert len(da) == 3


@pytest.mark.order(10)
def test_010_new_jobs():
    try:
        new_jobs_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(13)
def test_013_s03compute_cc():
    try:
        compute_cc_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(23)
def test_023_stack():
    db = connect()
    update_config(db, 'mov_stack', "(('1h','1h'),('6h','1h'))")
    stack_main('ref')
    reset_jobs(db, "STACK", alljobs=True)
    stack_main('mov')
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

@pytest.mark.order(29)
def test_029_stretching():
    from ..stretch2 import main as stretch_main
    db = connect()
    reset_jobs(db, "MWCS", alljobs=True)
    db.close()
    stretch_main()

@pytest.mark.order(100)
def test_100_plot_cctfime():
    from ..plots.ccftime import main as ccftime_main
    db = connect()
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                for filter in get_filters(db):
                    ccftime_main("%s.%s.%s" % (sta1.net, sta1.sta, loc1), "%s.%s.%s" % (sta2.net, sta2.sta, loc2), filter.ref, "ZZ", 1, show=False, outfile="?.png")
                    fn = 'ccftime %s-%s-%s-f%i-m%s_%s.png' % ("%s.%s.%s" % (sta1.net, sta1.sta, loc1), "%s.%s.%s" % (sta2.net, sta2.sta, loc2),
                                                              "ZZ", filter.ref, "1h", "1h")
                    assert os.path.isfile(fn)

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
                    fn = f'interferogram {sta_id1}-{sta_id2}-ZZ-f{filter.ref}-m1h_1h.png'
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
                    fn = f'spectime {sta_id1}-{sta_id2}-ZZ-f{filter.ref}-m1h_1h.png'
                    assert os.path.isfile(fn), f"{fn} doesn't exist"

@pytest.mark.order(102)
def test_102_plot_distance():
    distance_main(filterid=1, components="ZZ", show=False, outfile="?.png")
    fn = "distance ZZ-f1.png"
    assert os.path.isfile(fn), f"{fn} doesn't exist"

    distance_main(filterid=1, components="ZZ", show=False, outfile="?_refilter.png", refilter="1.0:3.0")
    fn = "distance ZZ-f1_refilter.png"
    assert os.path.isfile(fn), f"{fn} doesn't exist"