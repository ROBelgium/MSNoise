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
from .. import s01_scan_archive, FatalError
from ..scripts import msnoise as msnoise_script
from ..api import *
                #(connect, get_config, update_config, get_job_types,
                #   get_new_files, get_filters, get_station_pairs,
                #   get_components_to_compute, update_filter, Filter,
                #   get_stations, update_station, get_data_availability,
                #   count_data_availability_flags, is_next_job, get_next_job,
                #   Job, reset_jobs, build_ref_datelist, build_movstack_datelist,
                #   read_db_inifile, DvvMwcs, update_dvv_mwcs, get_dvv_mwcs,
                #   DvvMwcsDtt, update_dvv_mwcs_dtt, get_dvv_mwcs_dtt,
                #   DvvStretching)
from ..s01_scan_archive import parse_crondays
from ..s02_new_jobs import main as new_jobs_main
from ..s03_compute_no_rotation import main as compute_cc_main
from ..s04_stack_mov import main as stack_mov
from ..s04_stack_ref import main as stack_ref  # deprecated
from ..s04_stack_refstack import main as stack_refstack_main
from ..s05_compute_mwcs import main as compute_mwcs_main
from ..s06_compute_mwcs_dtt import main as compute_dtt_main
from ..s07_compute_dvv import main as compute_dvv_main
from ..psd_compute_rms import main as compute_rms_main
from ..psd_export_rms import main as export_rms_main
from ..ppsd_compute import main as ppsd_compute_main
from ..psd_to_hdf import main as psd_to_hdf_main

from ..plots.ccftime import main as ccftime_main
from ..plots.interferogram import main as interferogram_main
from ..plots.spectime import main as spectime_main
from ..plots.distance import main as distance_main
from ..plots.mwcs_dtt import main as dvv_main
from ..plots.data_availability import main as data_availability_main
from ..plots.wct_dvv import main as wct_dvv_main
from ..s02_preprocessing import main as preprocess_main
from ..s08_compute_wct import main as compute_wct_main
from ..s09_compute_wct_dtt import main as wavelet_dtt_main

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
    print(f"Tests will be running in the {test_dir} folder")

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
    from ..s00_installer import main
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
def test_002b_create_workflow():
    db = connect()
    for category in ['preprocess', 'cc', 'filter', 'stack', 'refstack',
                     'mwcs', 'mwcs_dtt', 'stretching', 'wavelet', 'wavelet_dtt', 'qc']:
        set_number = create_config_set(db, category)
        assert set_number == 1, f"Expected set_number=1 for {category}, got {set_number}"
    created_steps, _, err = create_workflow_steps_from_config_sets(db)
    assert err is None, f"Error creating workflow steps: {err}"
    assert created_steps > 0
    created_links, _, err = create_workflow_links_from_steps(db)
    assert err is None, f"Error creating workflow links: {err}"
    assert created_links > 0
    db.close()


@pytest.mark.order(4)
def test_002c_verify_workflow():
    db = connect()
    steps = get_workflow_steps(db)
    step_names = {s.step_name for s in steps}
    for name in ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1',
                 'mwcs_1', 'mwcs_dtt_1', 'stretching_1', 'wavelet_1', 'wavelet_dtt_1']:
        assert name in step_names, f"Workflow step '{name}' not found"
    links = get_workflow_links(db)
    assert len(links) > 0, "No workflow links created"
    # Verify key refstack links: stack_1→refstack_1, refstack_1→mwcs_1
    step_map = {s.step_name: s.step_id for s in steps}
    link_pairs = {(lk.from_step_id, lk.to_step_id) for lk in links}
    assert (step_map['stack_1'], step_map['refstack_1']) in link_pairs, \
        "Missing link stack_1 → refstack_1"
    assert (step_map['refstack_1'], step_map['mwcs_1']) in link_pairs, \
        "Missing link refstack_1 → mwcs_1"
    db.close()


@pytest.mark.order(5)
def test_003_set_and_config(setup_environment):
    db = connect()
    data_folder = setup_environment['data_folder']
    response_path = setup_environment['response_path']
    totests = [
        ['data_folder', data_folder],
        ['data_structure', 'PDF'],
        ['network', 'YA'],
        ['response_path', response_path]
    ]
    for key, value in totests:
        update_config(db, name=key, value=value, category='global', set_number=1)
        config_value = get_config(db, name=key, category='global', set_number=1)
        assert config_value == value, f"Configuration parameter {key} did not set correctly."


    update_config(db, 'components_to_compute', 'ZZ', category='cc', set_number=1)

    db.close()


@pytest.mark.order(6)
def test_004_set_and_get_filters():
    db = connect()
    # filter_1 configset was created in test_002b; update its frequency parameters
    update_config(db, 'freqmin', '0.01', category='filter', set_number=1)
    update_config(db, 'freqmax', '1.0', category='filter', set_number=1)
    for param in ['CC', 'SC', 'AC']:
        update_config(db, param, 'Y', category='filter', set_number=1)
    # Create a second filter configset
    set_number = create_config_set(db, 'filter')
    assert set_number == 2
    update_config(db, 'freqmin', '0.1', category='filter', set_number=2)
    update_config(db, 'freqmax', '1.0', category='filter', set_number=2)
    for param in ['CC', 'SC', 'AC']:
        update_config(db, param, 'Y', category='filter', set_number=2)
    # Register the new WorkflowStep and links for filter_2
    create_workflow_steps_from_config_sets(db)
    create_workflow_links_from_steps(db)
    # Verify both filter configsets
    for set_num, expected_freqmin in [(1, '0.01'), (2, '0.1')]:
        details = {d['name']: d['value'] for d in get_config_set_details(db, 'filter', set_num)}
        assert details['freqmin'] == expected_freqmin, \
            f"freqmin for filter_{set_num} wrong: {details.get('freqmin')}"
        assert details['freqmax'] == '1.0'
        assert details['CC'] == 'Y'
    db.close()


@pytest.mark.order(7)
def test_005_populate_station_table():
    from ..s00_populate_station_table import main
    try:
        ret = main()
        assert ret is True
    except:
        pytest.fail()

@pytest.mark.order(8)
def test_006_get_stations():
    db = connect()
    stations = get_stations(db).all()
    assert len(stations) == 3
    db.close()

@pytest.mark.order(9)
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

@pytest.mark.order(10)
def test_008_scan_archive(setup_environment):
    from ..s01_scan_archive import main
    try:
        main(init=True, threads=1)
    except:
        traceback.print_exc()
        pytest.fail()


    runner = setup_environment['runner']
    result = runner.invoke(msnoise_script.db_da_stations_update_loc_chan)
    assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}"
    # Optionally, add more assertions to check specific outputs in result.output


@pytest.mark.order(11)
def test_009_control_data_availability():
    db = connect()
    files = get_new_files(db)
    assert len(files) == 3
    # count_data_availability_flags was removed; query distinct flags directly
    from sqlalchemy.sql.expression import func
    flags = db.query(func.count(DataAvailability.flag),
                     DataAvailability.flag).group_by(DataAvailability.flag).all()
    assert len(flags) == 1
    for station in get_stations(db):
        for loc in station.locs():
            for chan in station.chans():
                da = get_data_availability(db, net=station.net, sta=station.sta, loc=loc, chan=chan)
                assert len(da) == 1

@pytest.mark.order(12)
def test_010_new_jobs():
    try:
        new_jobs_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(13)
def test_010b_preprocess_and_propagate():
    try:
        preprocess_main()
        new_jobs_main(after='preprocess')
    except:
        traceback.print_exc()
        pytest.fail()


@pytest.mark.order(14)
def test_011_control_jobs():
    db = connect()
    assert is_next_job_for_step(db, step_category='cc') is True
    db.close()

@pytest.mark.order(15)
def test_012_reset_jobs():
    db = connect()
    reset_jobs(db, 'cc_1', alljobs=True)
    db.close()

@pytest.mark.order(16)
def test_013_s03compute_cc():
    try:
        compute_cc_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(17)
def test_014_check_done_jobs():
    db = connect()
    jobs = get_job_types(db, 'cc_1')
    counts = {flag: count for count, flag in jobs}
    assert counts.get('D', 0) == 3, f"Expected 3 done CC jobs, got: {counts}"
    db.close()

@pytest.mark.order(18)
def test_015_check_cc_files():
    db = connect()
    output_folder = get_config(db, 'output_folder') or 'OUTPUT'
    cc_params = get_config_set_details(db, 'cc', 1, format='AttribDict')
    components_to_compute = cc_params.components_to_compute.split(',')
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    for filter_step in filter_steps:
        for components in components_to_compute:
            for (sta1, sta2) in get_station_pairs(db):
                for loc1 in sta1.locs():
                    for loc2 in sta2.locs():
                        sta1_id = f"{sta1.net}.{sta1.sta}.{loc1}"
                        sta2_id = f"{sta2.net}.{sta2.sta}.{loc2}"
                        tmp = os.path.join(output_folder, "preprocess_1", "cc_1",
                                           filter_step.step_name, "_output", "daily",
                                           components, sta1_id, sta2_id,
                                           f"2010-09-01.MSEED")
                        print("checking", tmp)
                        assert os.path.isfile(tmp), f"{tmp} does not exist"
    db.close()


@pytest.mark.order(19)
def test_017_reset_cc_jobs():
    db = connect()
    reset_jobs(db, 'cc_1', alljobs=True)
    db.close()

@pytest.mark.order(20)
def test_018_recompute_cc():
    test_013_s03compute_cc()

@pytest.mark.order(21)
def test_023_stack():
    db = connect()
    # Configure MOV stack (no more ref_begin/ref_end here — moved to refstack)
    update_config(db, 'mov_stack', "(('6h','6h'),('1D','1D'))", category='stack', set_number=1)
    # startdate/enddate belong to global configset
    update_config(db, 'startdate', '2009-01-01', category='global', set_number=1)
    update_config(db, 'enddate', '2011-01-01', category='global', set_number=1)

    # Configure the refstack_1 configset (ref_begin/ref_end now live here)
    update_config(db, 'ref_begin', '2009-01-01', category='refstack', set_number=1)
    update_config(db, 'ref_end', '2011-01-01', category='refstack', set_number=1)
    db.close()

    # Propagate CC→stack (MOV day jobs)
    new_jobs_main(after='cc')
    db = connect()

    # Run MOV stack
    stack_mov('mov')

    # Propagate stack→refstack (REF jobs) then run refstack
    new_jobs_main(after='stack')
    stack_refstack_main()

    # Verify refstack_1 REF job is done
    jobs_ref = get_job_types(db, 'refstack_1')
    counts_ref = {flag: count for count, flag in jobs_ref}
    assert counts_ref.get('D', 0) >= 1, \
        f"Expected at least 1 done refstack_1 job, got: {counts_ref}"

    # Re-run with wiener filter enabled
    update_config(db, 'wienerfilt', 'Y', category='stack', set_number=1)
    reset_jobs(db, "stack_1", alljobs=True)
    stack_mov('mov')
    reset_jobs(db, "refstack_1", alljobs=True)
    stack_refstack_main()

    db.close()


@pytest.mark.order(22)
def test_024_mwcs_param_update():
    db = connect()
    details = get_config_set_details(db, 'mwcs', 1)
    assert len(details) > 0, "MWCS configset 1 not found"
    update_config(db, 'mwcs_wlen', '10', category='mwcs', set_number=1)
    update_config(db, 'mwcs_step', '5', category='mwcs', set_number=1)
    details = {d['name']: d['value'] for d in get_config_set_details(db, 'mwcs', 1)}
    assert details['mwcs_wlen'] == '10', f"mwcs_wlen not updated, got: {details.get('mwcs_wlen')}"
    assert details['mwcs_step'] == '5', f"mwcs_step not updated, got: {details.get('mwcs_step')}"
    db.close()

@pytest.mark.order(23)
def test_025_mwcs():
    # after='refstack' creates mwcs/stretching/wavelet day jobs
    # (stack jobs trigger refstack jobs; refstack completion triggers dvv jobs)
    new_jobs_main(after='refstack')
    db = connect()
    assert is_next_job_for_step(db, step_category='mwcs'), \
        "No mwcs jobs created after refstack propagation"
    db.close()
    compute_mwcs_main()

@pytest.mark.order(24)
def test_026_mwcs_dtt_param_update():
    db = connect()
    details = get_config_set_details(db, 'mwcs_dtt', 1)
    assert len(details) > 0, "MWCS DTT configset 1 not found"
    update_config(db, 'dtt_minlag', '5.0', category='mwcs_dtt', set_number=1)
    update_config(db, 'dtt_width', '30.0', category='mwcs_dtt', set_number=1)
    details = {d['name']: d['value'] for d in get_config_set_details(db, 'mwcs_dtt', 1)}
    assert details['dtt_minlag'] == '5.0', f"dtt_minlag not updated, got: {details.get('dtt_minlag')}"
    assert details['dtt_width'] == '30.0', f"dtt_width not updated, got: {details.get('dtt_width')}"
    db.close()

@pytest.mark.order(27)
def test_027_dtt():
    new_jobs_main(after='mwcs')
    compute_dtt_main()

# @pytest.mark.order(28)
# def test_028_dvv():
#     compute_dvv_main()

@pytest.mark.order(31)
def test_031_stretching_param_update():
    db = connect()
    details = get_config_set_details(db, 'stretching', 1)
    assert len(details) > 0, "Stretching configset 1 not found"
    update_config(db, 'stretching_max', '0.05', category='stretching', set_number=1)
    update_config(db, 'stretching_nsteps', '500', category='stretching', set_number=1)
    details = {d['name']: d['value'] for d in get_config_set_details(db, 'stretching', 1)}
    assert details['stretching_max'] == '0.05'
    assert details['stretching_nsteps'] == '500'
    db.close()

@pytest.mark.order(32)
def test_032_stretching():
    # stretching jobs were created by new_jobs_main(after='refstack') in test_025
    from ..s10_stretching import main as stretch_main
    stretch_main()

@pytest.mark.order(33)
def test_033_create_fake_new_files(setup_environment):
    data_folder = setup_environment['data_folder']
    for f in sorted(glob.glob(os.path.join(data_folder, "2010", "*", "HHZ.D", "*"))):
        st = read(f)
        for tr in st:
            tr.stats.starttime += datetime.timedelta(days=1)
        out = f.replace("244", "245")
        st.write(out, format="MSEED")

    from ..s01_scan_archive import main
    try:
        main(init=False, threads=1)
    except:
        traceback.print_exc()
        pytest.fail()

    new_jobs_main()
    preprocess_main()
    new_jobs_main(after='preprocess')
    db = connect()
    jobs = get_job_types(db, 'cc_1')
    counts = {flag: count for count, flag in jobs}
    assert counts.get('D', 0) == 3, f"Expected 3 done cc_1 jobs, got: {counts}"
    assert counts.get('T', 0) == 3, f"Expected 3 todo cc_1 jobs, got: {counts}"
    db.close()

@pytest.mark.order(34)
def test_034_instrument_response(setup_environment):
    db = connect()
    response_path = setup_environment['response_path']
    update_config(db, 'response_path', response_path)
    update_config(db, 'remove_response', "Y")
    db.close()
    test_013_s03compute_cc()


@pytest.mark.order(35)
def test_035_wct_param_update():
    """Test updating and retrieving wavelet configset parameters."""
    db = connect()
    details = get_config_set_details(db, 'wavelet', 1)
    assert len(details) > 0, "Wavelet configset 1 not found"
    update_config(db, 'wct_freqmin', '0.1', category='wavelet', set_number=1)
    update_config(db, 'wct_freqmax', '2.0', category='wavelet', set_number=1)
    details = {d['name']: d['value'] for d in get_config_set_details(db, 'wavelet', 1)}
    assert details['wct_freqmin'] == '0.1'
    assert details['wct_freqmax'] == '2.0'
    db.close()

@pytest.mark.order(36)
def test_036_wct_dtt_param_update():
    """Test updating and retrieving wavelet_dtt configset parameters."""
    db = connect()
    details = get_config_set_details(db, 'wavelet_dtt', 1)
    assert len(details) > 0, "Wavelet DTT configset 1 not found"
    update_config(db, 'wct_minlag', '5.0', category='wavelet_dtt', set_number=1)
    update_config(db, 'wct_mincoh', '0.5', category='wavelet_dtt', set_number=1)
    details = {d['name']: d['value'] for d in get_config_set_details(db, 'wavelet_dtt', 1)}
    assert details['wct_minlag'] == '5.0'
    assert details['wct_mincoh'] == '0.5'
    db.close()
  
@pytest.mark.order(37)
def test_037_validate_stack_data():
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
    
@pytest.mark.order(38)
def test_038_stack_validation_handling():
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

@pytest.mark.order(39)
def test_039_wct_pipeline():
    try:
        # WCT jobs were created by new_jobs_main(after='refstack') in test_025
        # (stretching already has its own test_032)
        compute_wct_main()
        new_jobs_main(after='wavelet')
        wavelet_dtt_main()
    except:
        traceback.print_exc()
        pytest.fail()

@pytest.mark.order(100)
def test_100_plot_interferogram():
    db = connect()
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                sta_id1 = f"{sta1.net}.{sta1.sta}.{loc1}"
                sta_id2 = f"{sta2.net}.{sta2.sta}.{loc2}"
                for filter_step in filter_steps:
                    interferogram_main(sta_id1, sta_id2,
                                      filter_id=filter_step.set_number,
                                      components="ZZ", stack_item=1,
                                      show=False, outfile="?.png")
                    fn = f'interferogram {sta_id1}-{sta_id2}-ZZ-f{filter_step.set_number}-m6h_6h.png'
                    assert os.path.isfile(fn), f"{fn} doesn't exist"

@pytest.mark.order(101)
def test_101_plot_spectime():
    db = connect()
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                sta_id1 = f"{sta1.net}.{sta1.sta}.{loc1}"
                sta_id2 = f"{sta2.net}.{sta2.sta}.{loc2}"
                for filter_step in filter_steps:
                    spectime_main(sta_id1, sta_id2, 1, 1,
                                  filter_step.set_number, 1, 1,
                                  show=False, outfile="?.png")
                    fn = f'spectime {sta_id1}_{sta_id2}-ZZ-f{filter_step.set_number}-m6h_6h.png'
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
    dvv_main(components="ZZ", show=False, outfile="?.png")
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
    assert os.path.isfile("jobs.csv")
    assert os.path.isfile("data_availability.csv")

    # os.system("msnoise db import config --force")

#@pytest.mark.order(106)
#def test_106_plot_wct():
#    wct_dvv_main(filterid=1, wctid=1, dttid=1, components="ZZ", show=False, outfile="?.png")
#    fn = "wct ZZ-f1-dvv.png"
#    assert os.path.isfile(fn), f"{fn} doesn't exist"


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
    # MOV stack
    os.system("msnoise reset stack_1 --all")
    os.system("msnoise cc stack -m")
    # Refstack (REF jobs triggered by stack completion)
    os.system("msnoise new_jobs --after stack")
    os.system("msnoise reset refstack_1 --all")
    os.system("msnoise cc stack_refstack")
    # DVV jobs triggered by refstack completion
    os.system("msnoise new_jobs --after refstack")
    # MWCS
    os.system("msnoise reset mwcs_1 --all")
    os.system("msnoise cc dtt compute_mwcs")
    os.system("msnoise reset mwcs_dtt_1 --all")
    os.system("msnoise cc dtt compute_dtt")
    os.system("msnoise cc dtt plot dvv -s 0 -o ?.png")
    # Stretching
    os.system("msnoise reset stretching_1 --all")
    os.system("msnoise cc dtt compute_stretching")
    os.system("msnoise cc dtt plot dvvs -s 0 -o ?.png")
    # Wavelet
    os.system("msnoise reset wavelet_1 --all")
    os.system("msnoise cc dtt compute_wct")
    os.system("msnoise new_jobs --after wavelet")
    os.system("msnoise cc dtt compute_wct_dtt")
    os.system("msnoise cc dtt plot wct -s 0 -o ?.png")

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
        assert response.status_code == 200, "Error following route admin/"
        assert b"MSNoise Admin" in response.data
        for route in ["admin/config/","admin/station/",
                      "admin/dataavailability/", "admin/job/"]:
            response = test_client.get(route)
            assert response.status_code == 200, f"Error following route {route}"

        route = "admin/station/new/?url=/admin/station/"
        response = test_client.get(route)
        assert response.status_code == 200, f"Error following route {route}"
        #
        for route in ["admin/station/edit/?id=1&",
                      "admin/dataavailability/edit/?id=3"]:
            response = test_client.get(route, follow_redirects=True)
            assert response.status_code == 200, f"Error following route {route}"
        #
        # for route in ["admin/pairs.json",
        #               "admin/data_availability_flags.json",
        #               ]:
        #     response = test_client.get(route, follow_redirects=True)
        #     assert response.status_code == 200, f"Error following route {route}"


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
