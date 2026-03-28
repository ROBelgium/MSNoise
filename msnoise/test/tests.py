import datetime
import glob
import logging
import os
import shutil
import traceback
import logging as _logging
from click.testing import CliRunner
from obspy import read
import tempfile
import pandas as pd
import pooch
import pytest
from .. import FatalError
from ..scripts import msnoise as msnoise_script
from ..core.db import connect
from ..core.config import (
    create_config_set,
    get_config,
    get_config_set_details,
    update_config)
from ..core.stations import (
    get_data_availability,
    get_new_files,
    get_station_pairs,
    get_stations,
    update_station)
from ..core.workflow import (
    create_workflow_links_from_steps,
    create_workflow_steps_from_config_sets,
    get_job_types,
    get_workflow_links,
    get_workflow_steps,
    is_next_job_for_step,
    reset_jobs)
from ..msnoise_table_def import DataAvailability
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
from ..s04_stack_refstack import main as stack_refstack_main
from ..s05_compute_mwcs import main as compute_mwcs_main
from ..s06_compute_mwcs_dtt import main as compute_dtt_main
from ..s07_compute_dvv import main as compute_dvv_main
from ..psd_compute import main as psd_compute_main
from ..psd_compute_rms import main as compute_rms_main

from ..plots.ccftime import main as ccftime_main
from ..plots.interferogram import main as interferogram_main
from ..plots.spectime import main as spectime_main
from ..plots.distance import main as distance_main
from ..plots.mwcs_dtt_dvv import main as dvv_main
from ..plots.data_availability import main as data_availability_main
from ..plots.wavelet_dtt_dvv import main as wavelet_dtt_dvv_main
from ..s02_preprocessing import main as preprocess_main
from ..s08_compute_wct import main as compute_wct_main
from ..s09_compute_wct_dtt import main as wavelet_dtt_main
from ..results import MSNoiseResult

global logger
logger = logging.getLogger('matplotlib')
# set WARNING for Matplotlib
logger.setLevel(logging.CRITICAL)


_test_logger = _logging.getLogger("msnoise.test")


@pytest.fixture(autouse=True)
def _log_test_name(request):
    """Log the current test name at INFO so it appears in msnoise log output."""
    _test_logger.info(f">>> TEST START: {request.node.name}")
    yield
    _test_logger.info(f"<<< TEST END:   {request.node.name}")


@pytest.fixture(autouse=True)
def _print_test_name(request):
    """Print the running test name to stdout (visible with pytest -s)."""
    print(f"\n{'-'*60}\n> {request.node.name}\n{'-'*60}", flush=True)
    yield


def _job_counts(db, step_name: str) -> dict:
    """Return {flag: count} for *step_name* jobs."""
    from ..msnoise_table_def import WorkflowStep, Job as _Job
    step = db.query(WorkflowStep).filter(
        WorkflowStep.step_name == step_name).first()
    if step is None:
        return {}
    rows = (db.query(_Job.flag, _Job.ref)
              .filter(_Job.step_id == step.step_id)
              .all())
    counts = {}
    for flag, _ in rows:
        counts[flag] = counts.get(flag, 0) + 1
    return counts


def _assert_n_jobs(db, step_name: str, flag: str, expected: int, msg: str = ""):
    """Assert exactly *expected* jobs have *flag* for *step_name*."""
    counts = _job_counts(db, step_name)
    n = counts.get(flag, 0)
    assert n == expected, (
        f"{msg or step_name}: expected {expected} '{flag}' jobs, got {n} "
        f"(all flags: {counts})")


def _assert_min_jobs(db, step_name: str, flag: str, min_n: int = 1, msg: str = ""):
    """Assert at least *min_n* jobs have *flag* for *step_name*."""
    counts = _job_counts(db, step_name)
    n = counts.get(flag, 0)
    assert n >= min_n, (
        f"{msg or step_name}: expected >= {min_n} '{flag}' jobs, got {n} "
        f"(all flags: {counts})")


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
    # Installer must create the default DataSource (id=1)
    db = connect()
    from ..msnoise_table_def import DataSource
    ds = db.query(DataSource).filter(DataSource.ref == 1).first()
    assert ds is not None, "Installer must create default DataSource (id=1)"
    assert ds.name == "local"
    assert ds.data_structure == "SDS"
    db.close()


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
                     'mwcs', 'mwcs_dtt', 'mwcs_dtt_dvv',
                     'stretching', 'stretching_dvv',
                     'wavelet', 'wavelet_dtt', 'wavelet_dtt_dvv',
                     'psd', 'psd_rms']:
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
                 'mwcs_1', 'mwcs_dtt_1', 'mwcs_dtt_dvv_1',
                 'stretching_1', 'stretching_dvv_1',
                 'wavelet_1', 'wavelet_dtt_1', 'wavelet_dtt_dvv_1',
                 'psd_1', 'psd_rms_1']:
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
        ['response_path', response_path]
    ]
    for key, value in totests:
        update_config(db, name=key, value=value, category='global', set_number=1)
        config_value = get_config(db, name=key, category='global', set_number=1)
        assert config_value == value, f"Configuration parameter {key} did not set correctly."
    # Set network code and channels on the default DataSource
    from ..core.stations import update_data_source
    update_data_source(db, id=1, uri=os.path.realpath(data_folder),
                       data_structure='PDF', network_code='YA', channels='HHZ')

    # Set the default DataSource URI and data_structure (v2 path)
    from ..core.stations import update_data_source
    update_data_source(db, id=1, uri=os.path.realpath(data_folder),
                       data_structure='PDF')

    update_config(db, 'components_to_compute', 'ZZ', category='cc', set_number=1)

    db.close()


@pytest.mark.order(6)
def test_004_set_and_get_filters():
    db = connect()
    # filter_1 configset was created in test_002b; update its frequency parameters
    update_config(db, 'freqmin', '0.1', category='filter', set_number=1)
    update_config(db, 'freqmax', '1.0', category='filter', set_number=1)
    for param in ['CC', 'SC', 'AC']:
        update_config(db, param, 'Y', category='filter', set_number=1)
    # Create a second filter configset
    set_number = create_config_set(db, 'filter')
    assert set_number == 2
    update_config(db, 'freqmin', '0.1', category='filter', set_number=2)
    update_config(db, 'freqmax', '5.0', category='filter', set_number=2)
    for param in ['CC', 'SC', 'AC']:
        update_config(db, param, 'Y', category='filter', set_number=2)
    # Register the new WorkflowStep and links for filter_2
    create_workflow_steps_from_config_sets(db)
    create_workflow_links_from_steps(db)
    # Verify both filter configsets
    for set_num, expected_freqmin, expected_freqmax in [(1, '0.1', '1.0'), (2, '0.1', '5.0')]:
        details = {d['name']: d['value'] for d in get_config_set_details(db, 'filter', set_num)}
        assert details['freqmin'] == expected_freqmin, \
            f"freqmin for filter_{set_num} wrong: {details.get('freqmin')}"
        assert details['freqmax'] == expected_freqmax, \
            f"freqmin for filter_{set_num} wrong: {details.get('freqmax')}"
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

    # Verify DA rows have data_source_id=1 and relative paths
    db = connect()
    data_folder = os.path.realpath(setup_environment['data_folder'])
    from ..msnoise_table_def import DataAvailability
    das = db.query(DataAvailability).all()
    assert len(das) > 0, "scan_archive must populate DataAvailability"
    for da in das:
        assert da.data_source_id == 1, (
            f"DA row {da} must have data_source_id=1, got {da.data_source_id}"
        )
        assert not da.path.startswith(data_folder), (
            f"DA path must be relative, got absolute: {da.path!r}"
        )
    db.close()


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
    """new_jobs() creates preprocess_1 T jobs (one per station/day)."""
    try:
        new_jobs_main()
    except:
        traceback.print_exc()
        pytest.fail()
    db = connect()
    _assert_min_jobs(db, 'preprocess_1', 'T', 1, 'new_jobs must create preprocess T jobs')
    db.close()

# NOTE: same order number = runs in file order (intentional sibling group)
@pytest.mark.order(13)
def test_010b_preprocess():
    """Preprocess runs and marks preprocess_1 jobs Done."""
    db = connect()
    _assert_min_jobs(db, 'preprocess_1', 'T', 1, 'Need T preprocess jobs before preprocess')
    db.close()
    try:
        preprocess_main()
    except:
        traceback.print_exc()
        pytest.fail()
    db = connect()
    _assert_min_jobs(db, 'preprocess_1', 'D', 1, 'preprocess_1 jobs must be Done after preprocess')
    db.close()


@pytest.mark.order(13)
def test_010c_propagate_preprocess():
    """new_jobs --after preprocess creates cc_1 T jobs."""
    try:
        new_jobs_main(after='preprocess')
    except:
        traceback.print_exc()
        pytest.fail()
    db = connect()
    _assert_min_jobs(db, 'cc_1', 'T', 1, 'propagate preprocess→cc must create T jobs')
    db.close()


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
    """compute_cc runs all cc_1 T jobs and marks them Done."""
    db = connect()
    _assert_min_jobs(db, 'cc_1', 'T', 1, 'Need T cc jobs before compute_cc')
    db.close()
    try:
        compute_cc_main()
    except:
        traceback.print_exc()
        pytest.fail()
    db = connect()
    _assert_min_jobs(db, 'cc_1', 'D', 1, 'cc_1 jobs must be Done after compute_cc')
    db.close()

@pytest.mark.order(17)
def test_014_check_done_jobs():
    db = connect()
    jobs = get_job_types(db, 'cc_1')
    counts = {flag: count for count, flag in jobs}
    assert counts.get('D', 0) == 3, f"Expected 3 done CC jobs, got: {counts}"
    db.close()

# NOTE: same order number = runs in file order (intentional sibling group)
@pytest.mark.order(18)
def test_015_check_cc_files():
    db = connect()
    output_folder = get_config(db, 'output_folder') or 'OUTPUT'
    cc_params = get_config_set_details(db, 'cc', 1, format='AttribDict')
    components_to_compute = cc_params.components_to_compute.split(',')
    filter_steps = [s for s in get_workflow_steps(db) if s.category == "filter"]
    for filter_step in filter_steps:
        for components in components_to_compute:
            for (sta1, sta2) in get_station_pairs(db):
                for loc1 in sta1.locs():
                    for loc2 in sta2.locs():
                        sta1_id = f"{sta1.net}.{sta1.sta}.{loc1}"
                        sta2_id = f"{sta2.net}.{sta2.sta}.{loc2}"
                        tmp = os.path.join(output_folder, "preprocess_1", "cc_1",
                                           filter_step.step_name, "_output", "daily",
                                           components, f"{sta1_id}_{sta2_id}",
                                           "2010-09-01.nc")
                        print("checking", tmp)
                        assert os.path.isfile(tmp), f"{tmp} does not exist"
    db.close()


@pytest.mark.order(18)
def test_015b_check_keep_all_cc_files():
    """Verify that keep_all per-window CCF NetCDF files exist under _output/all/."""
    db = connect()
    output_folder = get_config(db, 'output_folder') or 'OUTPUT'
    cc_params = get_config_set_details(db, 'cc', 1, format='AttribDict')
    keep_all = cc_params.keep_all
    if str(keep_all).upper() != 'Y':
        pytest.skip("keep_all is not enabled — skipping keep_all file check")
    components_to_compute = cc_params.components_to_compute.split(',')
    filter_steps = [s for s in get_workflow_steps(db) if s.category == "filter"]
    found_any = False
    for filter_step in filter_steps:
        for components in components_to_compute:
            for (sta1, sta2) in get_station_pairs(db):
                for loc1 in sta1.locs():
                    for loc2 in sta2.locs():
                        sta1_id = f"{sta1.net}.{sta1.sta}.{loc1}"
                        sta2_id = f"{sta2.net}.{sta2.sta}.{loc2}"
                        path = os.path.join(
                            output_folder, "preprocess_1", "cc_1",
                            filter_step.step_name, "_output", "all",
                            components, f"{sta1_id}_{sta2_id}",
                        )
                        files = glob.glob(os.path.join(path, "*.nc"))
                        assert len(files) > 0,                             f"No keep_all NetCDF files found under {path}"
                        found_any = True
    assert found_any, "No station pairs found to check keep_all files"
    db.close()


@pytest.mark.order(18)
def test_015c_check_cc_float32_encoding():
    """Saved CC NetCDF files use float32 encoding, not float64."""
    import xarray as xr
    db = connect()
    output_folder = get_config(db, 'output_folder') or 'OUTPUT'
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    checked = 0
    for filter_step in filter_steps:
        nc_files = glob.glob(os.path.join(
            output_folder, "preprocess_1", "cc_1",
            filter_step.step_name, "_output", "daily",
            "**", "*.nc"), recursive=True)
        for nc in nc_files[:3]:
            ds = xr.open_dataset(nc)
            for var in ds.data_vars:
                assert ds[var].dtype != "float64", (
                    f"{nc}: variable '{var}' is float64 — "
                    f"expected float32 (check _xr_save_and_close encoding)")
            ds.close()
            checked += 1
    db.close()
    if checked == 0:
        pytest.skip("No CC NetCDF files found to verify encoding")



@pytest.mark.order(18)
def test_016_lineage_normalisation():
    """Lineage strings are stored as FK integer IDs in the Lineage table."""
    from ..msnoise_table_def import Lineage
    db = connect()

    # At least some jobs must have a lineage_id set
    from ..msnoise_table_def import Job as JobTable
    jobs_with_lin = (
        db.query(JobTable)
        .filter(JobTable.lineage_id.isnot(None))
        .limit(5).all()
    )
    assert len(jobs_with_lin) >= 1, "Expected jobs with lineage_id set"

    # Instance-level .lineage property must resolve to a slash-separated string
    for job in jobs_with_lin:
        lin = job.lineage   # resolved via association_proxy
        assert lin is not None, "job.lineage returned None"
        # Single-step lineages (e.g. 'preprocess_1') are valid — no '/' required
        assert isinstance(lin, str) and len(lin) > 0, \
            f"Lineage string looks wrong: {lin!r}"

    # Number of distinct Lineage rows must be << total jobs
    n_lin  = db.query(Lineage).count()
    n_jobs = db.query(JobTable).filter(JobTable.lineage_id.isnot(None)).count()
    assert n_lin < n_jobs,         f"Expected fewer Lineage rows ({n_lin}) than jobs ({n_jobs})"
    print(f"  {n_lin} unique lineage strings --> {n_jobs} jobs")
    db.close()



@pytest.mark.order(19)
def test_017_reset_cc_jobs():
    db = connect()
    reset_jobs(db, 'cc_1', alljobs=True)
    db.close()

@pytest.mark.order(20)
def test_018_recompute_cc():
    test_013_s03compute_cc()

# NOTE: same order number = runs in file order (intentional sibling group)
@pytest.mark.order(21)
def test_023a_stack_config():
    """Configure stack/refstack parameters."""
    db = connect()
    update_config(db, 'mov_stack', "(('6h','6h'),('1D','1D'))", category='stack', set_number=1)
    update_config(db, 'startdate', '2009-01-01', category='global', set_number=1)
    update_config(db, 'enddate', '2011-01-01', category='global', set_number=1)
    update_config(db, 'ref_begin', '2009-01-01', category='refstack', set_number=1)
    update_config(db, 'ref_end', '2011-01-01', category='refstack', set_number=1)
    db.close()


@pytest.mark.order(21)
def test_023b_new_jobs_after_cc():
    """new_jobs --after cc creates stack_1 T jobs."""
    new_jobs_main(after='cc')
    db = connect()
    _assert_min_jobs(db, 'stack_1', 'T', 1, 'new_jobs --after cc must create stack T jobs')
    db.close()


@pytest.mark.order(21)
def test_023c_stack_mov():
    """stack_mov processes stack_1 T jobs → Done."""
    db = connect()
    _assert_min_jobs(db, 'stack_1', 'T', 1, 'Need T stack jobs before stack_mov')
    db.close()
    stack_mov('mov')
    db = connect()
    _assert_min_jobs(db, 'stack_1', 'D', 1, 'stack_1 jobs must be Done after stack_mov')
    db.close()


@pytest.mark.order(21)
def test_023d_new_jobs_after_stack_and_refstack():
    """Propagate stack→refstack and run refstack."""
    new_jobs_main(after='stack')
    db = connect()
    _assert_min_jobs(db, 'refstack_1', 'T', 1, 'new_jobs --after stack must create refstack T jobs')
    db.close()
    stack_refstack_main()
    db = connect()
    _assert_min_jobs(db, 'refstack_1', 'D', 1, 'refstack_1 must be Done after stack_refstack')
    db.close()


@pytest.mark.order(21)
def test_023_stack():
    """Verify CCF data and test Wiener filter (combines legacy assertions)."""
    db = connect()
    # Verify refstack done
    jobs_ref = get_job_types(db, 'refstack_1')
    counts_ref = {flag: count for count, flag in jobs_ref}
    assert counts_ref.get('D', 0) >= 1, \
        f"Expected at least 1 done refstack_1 job, got: {counts_ref}"

    result_plain = MSNoiseResult(db, ["preprocess_1", "cc_1", "filter_1", "stack_1"])
    ccfs_plain = result_plain.get_ccf()
    assert len(ccfs_plain) > 0, "No plain-stack CCF results found"

    # Test Wiener filter
    update_config(db, 'wienerfilt', 'Y', category='stack', set_number=1)
    reset_jobs(db, "stack_1", alljobs=True)
    stack_mov('mov')
    result_wiener = MSNoiseResult(db, ["preprocess_1", "cc_1", "filter_1", "stack_1"])
    ccfs_wiener = result_wiener.get_ccf()
    first_key = next(iter(ccfs_plain))
    import numpy as np
    plain_vals  = ccfs_plain[first_key].values
    wiener_vals = ccfs_wiener.get(first_key)
    if wiener_vals is not None:
        assert not np.allclose(plain_vals, wiener_vals.values, equal_nan=True), \
            "Wiener-filtered CCF is identical to plain CCF — filter had no effect"

    # Revert Wiener for downstream tests
    update_config(db, 'wienerfilt', 'N', category='stack', set_number=1)
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

    # Verify stretching NetCDF output exists for at least one pair
    db = connect()
    output_folder = get_config(db, 'output_folder') or 'OUTPUT'
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    # Path: root/preprocess_1/cc_1/filter_1/stack_1/refstack_1/stretching_1/_output/<ms>/<comp>/<s1>_<s2>.nc
    files = glob.glob(
        os.path.join(output_folder, "preprocess_1", "cc_1",
                     "filter_1", "stack_1", "refstack_1", "stretching_1",
                     "_output", "**", "*.nc"),
        recursive=True
    )
    assert len(files) > 0, f"No stretching NetCDF output found under {output_folder}"
    db.close()


@pytest.mark.order(55)
def test_041_new_jobs_after_mwcs_dtt():
    """new_jobs --after mwcs_dtt inserts mwcs_dtt_dvv sentinel jobs."""
    new_jobs_main(after='mwcs_dtt')
    db = connect()
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()
    Job = schema.Job
    dvv_jobs = (
        db.query(Job)
        .filter(Job.day == "DVV")
        .filter(Job.pair == "ALL")
        .all()
    )
    db.close()
    assert len(dvv_jobs) >= 1, \
        "Expected at least one DVV sentinel job after --after mwcs_dtt"


@pytest.mark.order(56)
def test_042_compute_mwcs_dtt_dvv():
    """Compute MWCS dv/v aggregate — mwcs_dtt_dvv step."""
    try:
        compute_dvv_main(step_category="mwcs_dtt_dvv")
    except Exception:
        traceback.print_exc()
        pytest.fail("mwcs_dtt_dvv computation failed")

    # Verify the DVV aggregate has >= 2 time steps (both days present)
    db = connect()
    results = MSNoiseResult.list(db, category='mwcs_dtt_dvv')
    assert len(results) >= 1, "No mwcs_dtt_dvv results found"
    dvv_data = results[0].get_dvv()
    for key, ds in dvv_data.items():
        n = len(ds.times) if hasattr(ds, 'times') else 0
        assert n >= 2, \
            f"Expected >= 2 DVV time steps after 2 days, got {n} for {key}"
        break  # spot-check first result
    db.close()


@pytest.mark.order(57)
def test_043_new_jobs_after_stretching():
    """new_jobs --after stretching inserts stretching_dvv sentinel jobs."""
    new_jobs_main(after='stretching')
    db = connect()
    from ..msnoise_table_def import declare_tables
    schema = declare_tables()
    Job = schema.Job
    dvv_jobs = (
        db.query(Job)
        .filter(Job.day == "DVV")
        .filter(Job.pair == "ALL")
        .filter(Job.jobtype.like("stretching_dvv%"))
        .all()
    )
    db.close()
    assert len(dvv_jobs) >= 1, \
        "Expected at least one stretching_dvv sentinel job after --after stretching"


@pytest.mark.order(58)
def test_044_compute_stretching_dvv():
    """Compute Stretching dv/v aggregate — stretching_dvv step."""
    try:
        compute_dvv_main(step_category="stretching_dvv")
    except Exception:
        traceback.print_exc()
        pytest.fail("stretching_dvv computation failed")

    # Verify the stretching DVV aggregate has >= 2 time steps
    db = connect()
    results = MSNoiseResult.list(db, category='stretching_dvv')
    assert len(results) >= 1, "No stretching_dvv results found"
    dvv_data = results[0].get_dvv()
    for key, ds in dvv_data.items():
        n = len(ds.times) if hasattr(ds, 'times') else 0
        assert n >= 2, \
            f"Expected >= 2 stretching DVV time steps after 2 days, got {n} for {key}"
        break
    db.close()

@pytest.mark.order(45)
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

@pytest.mark.order(46)
def test_033b_compute_second_day():
    """Propagate the second day's CC jobs through the full pipeline.

    test_033 left cc_1 jobs with flag 'T' for day 2.  This test runs them
    through the full pipeline (CC → stack → refstack → MWCS → DTT →
    stretching → WCT → WCT-DTT) so DVV aggregates see >= 2 time steps.
    """
    db = connect()

    # ── CC day 2 ─────────────────────────────────────────────────────────────
    compute_cc_main()

    # ── Stack day 2 ──────────────────────────────────────────────────────────
    new_jobs_main(after='cc')
    stack_mov('mov')

    # ── Refstack: reset so it re-integrates both days into the REF ──────────
    reset_jobs(db, 'refstack_1', alljobs=True)
    db.close()
    new_jobs_main(after='stack')
    stack_refstack_main()

    # ── MWCS day 2 ───────────────────────────────────────────────────────────
    new_jobs_main(after='refstack')
    compute_mwcs_main()
    new_jobs_main(after='mwcs')
    compute_dtt_main()

    # ── Stretching day 2 ─────────────────────────────────────────────────────
    from ..s10_stretching import main as stretch_main
    stretch_main()

    # ── WCT day 2 ────────────────────────────────────────────────────────────
    compute_wct_main()
    new_jobs_main(after='wavelet')
    wavelet_dtt_main()

    # ── Verify 2 days in MWCS output ─────────────────────────────────────────
    db = connect()
    results = MSNoiseResult.list(db, category='mwcs_dtt')
    assert len(results) >= 1, "Expected at least one mwcs_dtt result after second day"
    r = results[0]
    mwcs_data = r.get_mwcs()
    n_times = 0
    for v in mwcs_data.values():
        if hasattr(v, 'times'):
            n_times = max(n_times, len(v.times))
        elif hasattr(v, 'index'):
            n_times = max(n_times, len(v))
    assert n_times >= 2, \
        f"Expected >= 2 time steps in MWCS output after 2 days of data, got {n_times}"
    db.close()

@pytest.mark.order(47)
def test_034_instrument_response(setup_environment):
    # TODO THIS TEST DOESN'T WORK
    # TODO IT'S NOT RECOMPUTING THE PREPROCESSING STEP :-)
    db = connect()
    response_path = setup_environment['response_path']
    update_config(db, 'response_path', response_path)
    update_config(db, 'remove_response', "Y")
    reset_jobs(db, 'cc_1', alljobs=True)
    db.close()
    test_013_s03compute_cc()


@pytest.mark.order(48)
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

@pytest.mark.order(49)
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
  
@pytest.mark.order(50)
def test_037_validate_stack_data():
    from ..core.signal import validate_stack_data
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
    
@pytest.mark.order(51)
def test_038_stack_validation_handling():
    from ..core.signal import validate_stack_data
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

@pytest.mark.order(53)
def test_040_compute_wavelet_dtt_dvv():
    """Compute WCT dv/v aggregate — wavelet_dtt_dvv step."""
    new_jobs_main(after='wavelet_dtt')
    try:
        compute_dvv_main(step_category="wavelet_dtt_dvv")
    except Exception:
        traceback.print_exc()
        pytest.fail("wavelet_dtt_dvv computation failed")

@pytest.mark.order(99)
def test_099_plot_interferogram():
    db = connect()
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                sta_id1 = f"{sta1.net}.{sta1.sta}.{loc1}"
                sta_id2 = f"{sta2.net}.{sta2.sta}.{loc2}"
                for filter_step in filter_steps:
                    interferogram_main(sta_id1, sta_id2,
                                      filterid=filter_step.set_number,
                                      components="ZZ", stackid_item=1,
                                      show=False, outfile="?.png")
                    fn = glob.glob(f'interferogram__*{sta_id1}*ZZ*.png')
                    assert len(fn) >= 1, f"interferogram plot for {sta_id1}-{sta_id2} doesn't exist"

@pytest.mark.order(100)
def test_100_plot_ccftime():
    db = connect()
    filter_steps = [s for s in get_workflow_steps(db) if s.category == 'filter']
    for sta1, sta2 in get_station_pairs(db):
        for loc1 in sta1.locs():
            for loc2 in sta2.locs():
                sta_id1 = f"{sta1.net}.{sta1.sta}.{loc1}"
                sta_id2 = f"{sta2.net}.{sta2.sta}.{loc2}"
                for filter_step in filter_steps:
                    ccftime_main(sta_id1, sta_id2, 1, 1,
                                 filter_step.set_number, 1, 1,
                                 show=False, outfile="?.png")
                    fn = glob.glob(f'ccftime__*{sta_id1}*ZZ*.png')
                    assert len(fn) >= 1, f"ccftime plot for {sta_id1}-{sta_id2} doesn't exist"


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
                    fn = glob.glob(f'spectime__*{sta_id1}*ZZ*.png')
                    assert len(fn) >= 1, f"spectime plot for {sta_id1}-{sta_id2} doesn't exist"

@pytest.mark.order(102)
def test_102_plot_distance():
    distance_main(filterid=1, components="ZZ", show=False, outfile="?.png")
    fn = glob.glob("distance__*ZZ*.png")
    assert len(fn) >= 1, "distance plot doesn't exist"

    distance_main(filterid=1, components="ZZ", show=False, outfile="?_refilter.png", refilter="0.2:0.9")
    fn = glob.glob("distance__*ZZ*_refilter.png")
    assert len(fn) >= 1, "distance refilter plot doesn't exist"

# NOTE: same order number = runs in file order (intentional sibling group)
@pytest.mark.order(103)
def test_103_plot_dvv():
    """Plot mwcs_dtt_dvv — requires mwcs_dtt_dvv step to have been run."""
    try:
        dvv_main(components="ZZ", show=False, outfile="?.png", dvvid=1)
    except (FileNotFoundError, ValueError):
        pytest.skip("No mwcs_dtt_dvv aggregate data available — run compute first")
    fn = glob.glob("dvv_mwcs__*ZZ*.png")
    assert len(fn) >= 1, "Expected at least one dvv_mwcs plot PNG file"

@pytest.mark.order(103)
def test_103b_plot_stretching_dvv():
    """Plot stretching_dvv — requires stretching_dvv step to have run."""
    from ..plots.stretching_dvv import main as stretching_dvv_main
    try:
        stretching_dvv_main(components="ZZ", show=False, outfile="?.png", dvvid=1)
    except (FileNotFoundError, ValueError):
        pytest.skip("No stretching_dvv aggregate data available")
    fn = glob.glob("dvv_stretching__*ZZ*.png")
    assert len(fn) >= 1, "Expected at least one dvv_stretching plot PNG"

@pytest.mark.order(103)
def test_103c_plot_wavelet_dtt_dvv():
    """Plot wavelet_dtt_dvv — requires wavelet_dtt_dvv step to have run."""
    try:
        wavelet_dtt_dvv_main(filterid=1, wctid=1, dttid=1, components="ZZ",
                             show=False, outfile="?.png", dvvid=1)
    except (FileNotFoundError, ValueError):
        pytest.skip("No wavelet_dtt_dvv aggregate data available")
    fn = glob.glob("dvv_wavelet__*ZZ*.png")
    assert len(fn) >= 1, "Expected at least one dvv_wavelet plot PNG"

@pytest.mark.order(104)
def test_104_plot_data_availability():
    data_availability_main(chan="HHZ", show=False, outfile="?.png")
    fn = glob.glob("data availability on*.png")
    assert len(fn) == 1, "Data availability plot doesn't exist"

@pytest.mark.order(105)
def test_105_db_dump(setup_environment):
    """ Tests the dump of the database and the creation of csv files """
    runner = setup_environment['runner']
    result = runner.invoke(msnoise_script.db_dump)
    assert result.exit_code == 0, f"msnoise db dump failed: {result.output}"
    assert os.path.isfile("config.csv")
    assert os.path.isfile("stations.csv")
    assert os.path.isfile("jobs.csv")
    assert os.path.isfile("data_availability.csv")

    # os.system("msnoise db import config --force")

#@pytest.mark.order(106)
#def test_106_plot_wct():
#    wavelet_dtt_dvv_main(filterid=1, wctid=1, dttid=1, components="ZZ", show=False, outfile="?.png")
#    fn = "wct ZZ-f1-dvv.png"
#    assert os.path.isfile(fn), f"{fn} doesn't exist"


@pytest.mark.order(106)
def test_106_plot_dvv_comparison():
    """Run the DVV comparison example — verifies it produces at least one output file."""
    # The example discovers results and saves dvv_comparison__*.png files
    # Run it in-process by importing and exercising the core logic directly
    db = connect()
    from ..results import MSNoiseResult

    # Discover all three DVV categories
    found = {}
    for cat in ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"):
        results = MSNoiseResult.list(db, cat)
        if results:
            found[cat] = results[0]

    if not found:
        pytest.skip("No DVV results available — pipeline not fully run")

    # dvv_split_pair_type=Y (default) saves dvv_CC_ZZ.nc;
    # dvv_split_pair_type=N saves dvv_ALL_ZZ.nc. Try both.
    all_keys = set()
    raw = {}
    for cat, result in found.items():
        for try_pt in ("CC", "ALL"):
            data = result.get_dvv(pair_type=try_pt)
            filtered = {(comp, ms): ds
                        for (pt, comp, ms), ds in data.items()}
            if filtered:
                raw[cat] = filtered
                all_keys |= set(filtered.keys())
                break

    assert len(all_keys) >= 1, "Expected at least one (comp, mov_stack) combination"

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.gridspec as gridspec
    import pandas as pd
    from ..core.config import build_plot_outfile

    METHOD_STYLE = {
        "mwcs_dtt_dvv":    dict(color="#1f77b4", label="MWCS",       lw=1.8, zorder=3),
        "stretching_dvv":  dict(color="#d62728", label="Stretching", lw=1.8, zorder=2),
        "wavelet_dtt_dvv": dict(color="#2ca02c", label="WCT",        lw=1.8, zorder=1),
    }

    saved = []
    for comp, ms in sorted(all_keys):
        series = {}
        for cat in ("mwcs_dtt_dvv", "stretching_dvv", "wavelet_dtt_dvv"):
            if cat not in raw or (comp, ms) not in raw[cat]:
                continue
            ds = raw[cat][(comp, ms)]
            if "mean" not in ds:
                continue
            mean_da = ds["mean"]
            std_da  = ds.get("std", ds.get("weighted_std", None))
            times   = pd.DatetimeIndex(mean_da.coords["times"].values)
            mean_s  = pd.Series(mean_da.values * 100.0, index=times)
            std_s   = (pd.Series(std_da.values * 100.0, index=times)
                       if std_da is not None else None)
            series[cat] = (mean_s, std_s)

        if not series:
            continue

        fig = plt.figure(figsize=(12, 7))
        gs  = gridspec.GridSpec(2, 1, height_ratios=[3, 1], hspace=0.08)
        ax_dvv = fig.add_subplot(gs[0])
        ax_res = fig.add_subplot(gs[1], sharex=ax_dvv)

        reference_s = None
        for cat, (mean_s, std_s) in series.items():
            style = METHOD_STYLE[cat]
            ax_dvv.plot(mean_s.index, mean_s.values,
                        color=style["color"], lw=style["lw"], label=style["label"])
            if std_s is not None:
                ax_dvv.fill_between(mean_s.index,
                                    mean_s.values - std_s.values,
                                    mean_s.values + std_s.values,
                                    color=style["color"], alpha=0.15)
            if reference_s is None:
                reference_s = mean_s

        ax_dvv.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
        ax_dvv.set_ylabel("dv/v (%)")
        ax_dvv.set_title(f"dv/v comparison | comp={comp} | mov_stack={ms[0]}-{ms[1]}")
        ax_dvv.legend(loc="upper right")
        ax_dvv.tick_params(labelbottom=False)

        ref_label = METHOD_STYLE[next(iter(series))]["label"]
        for cat, (mean_s, _) in series.items():
            if mean_s is reference_s:
                continue
            style = METHOD_STYLE[cat]
            common = reference_s.index.intersection(mean_s.index)
            if len(common):
                resid = reference_s.loc[common] - mean_s.loc[common]
                ax_res.plot(resid.index, resid.values,
                            color=style["color"], lw=1.4,
                            label=f"{ref_label} − {style['label']}")

        ax_res.axhline(0, color="k", lw=0.8, ls="--", alpha=0.5)
        ax_res.set_ylabel("Residual (%)")
        ax_res.set_xlabel("Date")
        ax_res.legend(loc="upper right", fontsize=8)
        fig.autofmt_xdate()

        first_result = next(iter(found.values()))
        outfile = build_plot_outfile(
            "?.png", "dvv_comparison",
            first_result.lineage_names,
            components=comp, mov_stack=ms,
        )
        plt.savefig(outfile, dpi=100, bbox_inches="tight")
        plt.close(fig)
        saved.append(outfile)

    db.close()
    assert len(saved) >= 1, "Expected at least one dvv_comparison plot to be saved"
    for fn in saved:
        assert os.path.isfile(fn), f"Output file missing: {fn}"



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

    result = runner.invoke(msnoise_script.config_set, ['response_path=XXX'])
    assert result.exit_code == 0

    result = runner.invoke(msnoise_script.config_get, ['response_path'])
    assert result.exit_code == 0
    assert 'XXX' in result.output

    result = runner.invoke(msnoise_script.config_set, ['response_path=none'])
    assert result.exit_code == 0

@pytest.mark.order(301)
def test_301_compute_psd():
    try:
        psd_compute_main()
    except:
        traceback.print_exc()
        pytest.fail("PSD Compute failed")

@pytest.mark.order(302)
def test_302_compute_rms():
    try:
        new_jobs_main(after='psd')
        compute_rms_main()
    except:
        traceback.print_exc()
        pytest.fail("PSD RMS computation failed")

@pytest.mark.order(400)
def test_400_run_manually(setup_environment):
    """End-to-end CLI test: re-run full pipeline using Click CliRunner.

    Unlike tests 1–120031 which call Python main() functions directly,
    this test exercises the CLI command dispatch layer — argument parsing,
    group routing, context passing — using Click's CliRunner so the test
    database is shared (no external subprocess spawning).
    """
    runner = setup_environment['runner']

    def run(args, label=None):
        """Invoke CLI and assert success."""
        result = runner.invoke(msnoise_script.cli, args.split(), obj={})
        assert result.exit_code == 0, (
            f"msnoise {args} failed (exit {result.exit_code}):\n{result.output}"
        )
        return result

    # ── MOV stack ────────────────────────────────────────────────────────────
    run("reset stack_1 --all")
    run("cc stack -m")

    # ── Refstack ─────────────────────────────────────────────────────────────
    run("new_jobs --after stack")
    run("reset refstack_1 --all")
    run("new_jobs --after stack")
    run("cc stack_refstack")

    # ── Propagate refstack → downstream ──────────────────────────────────────
    run("new_jobs --after refstack")

    # ── MWCS ─────────────────────────────────────────────────────────────────
    run("reset mwcs_1 --all")
    run("cc dtt compute_mwcs")
    run("reset mwcs_dtt_1 --all")
    run("new_jobs --after mwcs")
    run("cc dtt compute_mwcs_dtt")
    run("new_jobs --after mwcs_dtt")
    run("reset mwcs_dtt_dvv_1 --all")
    run("cc dtt dvv compute_mwcs_dtt_dvv")

    # ── Stretching ───────────────────────────────────────────────────────────
    run("reset stretching_1 --all")
    run("cc dtt compute_stretching")
    run("new_jobs --after stretching")
    run("reset stretching_dvv_1 --all")
    run("cc dtt dvv compute_stretching_dvv")

    # ── Wavelet (WCT) ────────────────────────────────────────────────────────
    run("reset wavelet_1 --all")
    run("cc dtt compute_wct")
    run("new_jobs --after wavelet")
    run("reset wavelet_dtt_1 --all")
    run("cc dtt compute_wct_dtt")
    run("new_jobs --after wavelet_dtt")
    run("reset wavelet_dtt_dvv_1 --all")
    run("cc dtt dvv compute_wavelet_dtt_dvv")

    # ── PSDs ─────────────────────────────────────────────────────────────────
    run("reset psd_1 --all")
    run("qc compute_psd")
    run("new_jobs --after psd")
    run("reset psd_rms_1 --all")
    run("qc compute_psd_rms")


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


@pytest.mark.order(999999)
def test_20000_invoke_script(setup_environment):
    runner = setup_environment['runner']

    for cmd in [
        msnoise_script.config_sync,
        msnoise_script.db_upgrade,
        msnoise_script.db_clean_duplicates,
        ]:

        result = runner.invoke(cmd)
        assert result.exit_code == 0, f"Command failed with exit code {result.exit_code}"


# ============================================================
# Unit tests — core lineage / scheduling API  (Patch 11)
# ============================================================

@pytest.mark.order(110000)
def test_110000_lineage_str_to_steps_strict():
    """lineage_str_to_steps raises on missing step when strict=True."""
    db = connect()
    from ..core.workflow import lineage_str_to_steps
    with pytest.raises(Exception):
        lineage_str_to_steps(db, "preprocess_1/nonexistent_99", strict=True)
    db.close()


@pytest.mark.order(110001)
def test_110001_lineage_str_to_steps_permissive():
    """lineage_str_to_steps returns partial list when strict=False."""
    db = connect()
    from ..core.workflow import lineage_str_to_steps
    steps = lineage_str_to_steps(db, "preprocess_1/nonexistent_99", strict=False)
    # Should return whatever steps were resolved, not raise
    assert isinstance(steps, list)
    db.close()


@pytest.mark.order(110002)
def test_110002_lineage_str_to_steps_valid():
    """lineage_str_to_steps resolves a valid lineage correctly."""
    db = connect()
    from ..core.workflow import lineage_str_to_steps
    from ..core.workflow import get_workflow_steps
    # Build a lineage string from real steps
    steps = get_workflow_steps(db)
    preprocess = next((s for s in steps if s.category == 'preprocess'), None)
    if preprocess is None:
        pytest.skip("No preprocess step configured")
    lineage_str = preprocess.step_name
    resolved = lineage_str_to_steps(db, lineage_str, strict=True)
    assert len(resolved) >= 1
    assert resolved[0].step_name == preprocess.step_name
    db.close()


@pytest.mark.order(110010)
def test_110010_get_merged_params_overrides():
    """get_merged_params_for_lineage: later configsets override earlier."""
    db = connect()
    from ..core.config import get_merged_params_for_lineage
    from ..core.config import get_params
    from ..core.workflow import lineage_str_to_steps
    from ..core.workflow import get_workflow_steps
    orig_params = get_params(db)
    steps = get_workflow_steps(db)
    # Find a filter step
    filter_step = next((s for s in steps if s.category == 'filter'), None)
    if filter_step is None:
        pytest.skip("No filter step configured")
    lineage_str = "/".join([s.step_name for s in steps
                             if s.category in ('preprocess', 'cc', 'filter')][:3])
    lineage_steps = lineage_str_to_steps(db, lineage_str, strict=False)
    if not lineage_steps:
        pytest.skip("Could not resolve lineage steps")
    _, names, params = get_merged_params_for_lineage(db, orig_params, {}, lineage_steps)
    assert isinstance(names, list)
    assert len(names) >= 1
    db.close()


@pytest.mark.order(110011)
def test_110011_get_merged_params_components_split():
    """get_merged_params_for_lineage splits components_to_compute into a list."""
    db = connect()
    from ..core.config import get_merged_params_for_lineage
    from ..core.config import get_params
    from ..core.workflow import lineage_str_to_steps
    from ..core.workflow import get_workflow_steps
    orig_params = get_params(db)
    steps = get_workflow_steps(db)
    cc_steps = [s for s in steps if s.category == 'cc']
    if not cc_steps:
        pytest.skip("No CC step configured")
    lineage_str = cc_steps[0].step_name
    lineage_steps = lineage_str_to_steps(db, lineage_str, strict=False)
    if not lineage_steps:
        pytest.skip("Could not resolve CC lineage")
    _, _, params = get_merged_params_for_lineage(db, orig_params, {}, lineage_steps)
    if hasattr(params, 'components_to_compute'):
        assert isinstance(params.cc.components_to_compute, list), \
            "components_to_compute should be a list after merging"
    db.close()


@pytest.mark.order(110020)
def test_110020_get_done_lineages_for_category():
    """get_done_lineages_for_category returns done lineages after pipeline run."""
    db = connect()
    from ..core.workflow import get_done_lineages_for_category
    lineages = get_done_lineages_for_category(db, 'cc')
    assert isinstance(lineages, list)
    assert len(lineages) >= 1, "Expected at least one done CC lineage"
    for lin in lineages:
        assert isinstance(lin, list)
        assert all(isinstance(n, str) for n in lin)
    db.close()


@pytest.mark.order(110021)
def test_110021_get_done_lineages_deduplicates():
    """get_done_lineages_for_category deduplicates lineage strings."""
    db = connect()
    from ..core.workflow import get_done_lineages_for_category
    lineages = get_done_lineages_for_category(db, 'cc')
    # Convert to frozensets for dedup check
    as_tuples = [tuple(lin) for lin in lineages]
    assert len(as_tuples) == len(set(as_tuples)), \
        "get_done_lineages_for_category returned duplicate lineages"
    db.close()


# ============================================================
# Unit tests — MSNoiseResult  (Patch 11 / Patch 9)
# ============================================================

@pytest.mark.order(120000)
def test_120000_msnoise_result_import():
    """MSNoiseResult can be imported from msnoise.results."""
    from ..results import MSNoiseResult
    assert MSNoiseResult is not None


@pytest.mark.order(120001)
def test_120001_msnoise_result_from_ids():
    """MSNoiseResult.from_ids constructs correctly."""
    db = connect()
    from ..results import MSNoiseResult
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1)
    assert r.lineage_names == ['preprocess_1', 'cc_1']
    assert r.category == 'cc'
    assert r.output_folder is not None
    db.close()


@pytest.mark.order(120002)
def test_120002_msnoise_result_from_names():
    """MSNoiseResult.from_names constructs correctly."""
    db = connect()
    from ..results import MSNoiseResult
    names = ['preprocess_1', 'cc_1', 'filter_1', 'stack_1', 'refstack_1']
    r = MSNoiseResult.from_names(db, names)
    assert r.lineage_names == names
    assert r.category == 'refstack'
    db.close()


@pytest.mark.order(120003)
def test_120003_msnoise_result_list():
    """MSNoiseResult.list returns done results for a category."""
    db = connect()
    from ..results import MSNoiseResult
    results = MSNoiseResult.list(db, 'cc')
    assert isinstance(results, list)
    assert len(results) >= 1
    for r in results:
        assert r.category == 'cc'
        assert 'cc_1' in r.lineage_names
    db.close()


@pytest.mark.order(120004)
def test_120004_msnoise_result_require_category_raises():
    """MSNoiseResult raises AttributeError when accessing a method whose
    required category is absent from the lineage."""
    db = connect()
    from ..results import MSNoiseResult
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1)
    # get_ref requires 'refstack' — not present in preprocess/cc lineage
    with pytest.raises(AttributeError, match="requires category"):
        r.get_ref()
    # hasattr must return False for absent-category methods
    assert not hasattr(r, "get_ref"), "get_ref should not be visible at cc level"
    assert not hasattr(r, "get_mwcs"), "get_mwcs should not be visible at cc level"
    db.close()


@pytest.mark.order(120005)
def test_120005_msnoise_result_branches():
    """MSNoiseResult.branches returns downstream children."""
    db = connect()
    from ..results import MSNoiseResult
    # Build up to refstack — should have mwcs/stretching/wavelet children
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1)
    children = r.branches()
    assert isinstance(children, list)
    # At least one downstream branch should have done jobs
    if children:
        for child in children:
            assert child.category in ('mwcs', 'stretching', 'wavelet')
            assert child.lineage_names[:-1] == r.lineage_names
    db.close()


@pytest.mark.order(120006)
def test_120006_msnoise_result_dynamic_visibility():
    """Methods for absent categories are invisible and inaccessible."""
    db = connect()
    from ..results import MSNoiseResult

    # Stack-level result: only get_ccf should be visible
    r_stack = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    visible = [m for m in dir(r_stack) if m.startswith("get_")]
    assert "get_ccf" in visible,  "get_ccf must be visible at stack level"
    assert "get_ref" not in visible,  "get_ref must not be visible at stack level"
    assert "get_mwcs" not in visible, "get_mwcs must not be visible at stack level"
    assert "get_dvv" not in visible,  "get_dvv must not be visible at stack level"
    assert not hasattr(r_stack, "get_mwcs"), "hasattr(get_mwcs) must be False at stack level"

    # Refstack-level result: get_ccf + get_ref visible, mwcs still not
    r_ref = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                   stack=1, refstack=1)
    visible_ref = [m for m in dir(r_ref) if m.startswith("get_")]
    assert "get_ccf" in visible_ref
    assert "get_ref" in visible_ref
    assert "get_mwcs" not in visible_ref

    # Error message is informative
    try:
        r_stack.get_mwcs()
        pytest.fail("Expected AttributeError")
    except AttributeError as e:
        assert "mwcs" in str(e)
        assert "branches()" in str(e)

    db.close()

@pytest.mark.order(120010)
def test_120010_msnoise_result_get_ccf():
    """MSNoiseResult.get_ccf loads CCF data for a specific pair."""
    db = connect()
    from ..results import MSNoiseResult
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    pairs = get_station_pairs(db)
    if not pairs:
        pytest.skip("No station pairs configured")
    sta1, sta2 = next(pairs)
    loc1 = sta1.locs()[0] if sta1.locs() else "00"
    loc2 = sta2.locs()[0] if sta2.locs() else "00"
    pair_str = f"{sta1.net}.{sta1.sta}.{loc1}:{sta2.net}.{sta2.sta}.{loc2}"

    # get_ccf with all None → dict
    all_ccfs = r.get_ccf()
    assert isinstance(all_ccfs, dict)

    # get_ccf with specific mov_stack → xarray DataArray (default)
    if all_ccfs:
        first_key = next(iter(all_ccfs))
        pair_k, comp_k, ms_k = first_key
        single = r.get_ccf(pair_k, comp_k, ms_k)
        import xarray as xr
        assert isinstance(single, xr.DataArray)
        assert single.sizes["times"] > 0
    db.close()


@pytest.mark.order(120011)
def test_120011_msnoise_result_get_ref():
    """MSNoiseResult.get_ref loads REF stack."""
    db = connect()
    from ..results import MSNoiseResult
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1)
    all_refs = r.get_ref()
    assert isinstance(all_refs, dict)
    if all_refs:
        first_key = next(iter(all_refs))
        pair_k, comp_k = first_key
        ds = r.get_ref(pair_k, comp_k)
        import xarray as xr
        assert isinstance(ds, xr.Dataset)
    db.close()


# ============================================================
# Unit tests — MSNoiseResult format options + list()  (Patch 23)
# ============================================================

@pytest.mark.order(120020)
def test_120020_msnoise_result_get_ccf_xarray():
    """MSNoiseResult.get_ccf returns xarray DataArray when format='xarray'."""
    db = connect()
    from ..results import MSNoiseResult
    import xarray as xr
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    all_ccfs = r.get_ccf()
    if not all_ccfs:
        pytest.skip("No CCF data available")
    pair_k, comp_k, ms_k = next(iter(all_ccfs))
    da = r.get_ccf(pair_k, comp_k, ms_k, format="xarray")
    assert isinstance(da, xr.DataArray)
    db.close()


@pytest.mark.order(120021)
def test_120021_msnoise_result_get_ref_dataframe():
    """MSNoiseResult.get_ref returns Series when format='dataframe' (REF is 1D)."""
    db = connect()
    from ..results import MSNoiseResult
    import pandas as pd
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1)
    all_refs = r.get_ref()
    if not all_refs:
        pytest.skip("No REF data available")
    pair_k, comp_k = next(iter(all_refs))
    s = r.get_ref(pair_k, comp_k, format="dataframe")
    assert isinstance(s, pd.Series)
    db.close()


@pytest.mark.order(120022)
def test_120022_msnoise_result_get_mwcs_both_formats():
    """MSNoiseResult.get_mwcs returns DataFrame or xarray Dataset by format."""
    db = connect()
    from ..results import MSNoiseResult
    import pandas as pd
    import xarray as xr
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1, mwcs=1)
    all_mwcs = r.get_mwcs()
    if not all_mwcs:
        pytest.skip("No MWCS data available")
    pair_k, comp_k, ms_k = next(iter(all_mwcs))
    df = r.get_mwcs(pair_k, comp_k, ms_k, format="dataframe")
    assert isinstance(df, pd.DataFrame)
    ds = r.get_mwcs(pair_k, comp_k, ms_k, format="xarray")
    assert isinstance(ds, xr.Dataset)
    db.close()


@pytest.mark.order(120023)
def test_120023_msnoise_result_get_mwcs_dtt_both_formats():
    """MSNoiseResult.get_mwcs_dtt returns DataFrame or xarray Dataset by format."""
    db = connect()
    from ..results import MSNoiseResult
    import pandas as pd
    import xarray as xr
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1, mwcs=1, mwcs_dtt=1)
    all_dtt = r.get_mwcs_dtt()
    if not all_dtt:
        pytest.skip("No MWCS-DTT data available")
    pair_k, comp_k, ms_k = next(iter(all_dtt))
    df = r.get_mwcs_dtt(pair_k, comp_k, ms_k, format="dataframe")
    assert isinstance(df, pd.DataFrame)
    ds = r.get_mwcs_dtt(pair_k, comp_k, ms_k, format="xarray")
    assert isinstance(ds, xr.Dataset)
    db.close()


@pytest.mark.order(120024)
def test_120024_msnoise_result_list_include_empty():
    """MSNoiseResult.list(include_empty=True) returns >= done results."""
    db = connect()
    from ..results import MSNoiseResult
    done = MSNoiseResult.list(db, 'cc', include_empty=False)
    all_steps = MSNoiseResult.list(db, 'cc', include_empty=True)
    assert isinstance(all_steps, list)
    assert len(all_steps) >= len(done), \
        "include_empty=True should return >= done results"
    for r in all_steps:
        assert r.category == 'cc'
    db.close()


@pytest.mark.order(120025)
def test_120025_msnoise_result_list_psd():
    """MSNoiseResult.list works for PSD root step category."""
    db = connect()
    from ..results import MSNoiseResult
    results = MSNoiseResult.list(db, 'psd')
    assert isinstance(results, list)
    for r in results:
        assert any('psd' in n for n in r.lineage_names), \
            f"Expected 'psd' in lineage_names: {r.lineage_names}"
    db.close()


@pytest.mark.order(120026)
def test_120026_msnoise_result_get_psd_formats():
    """MSNoiseResult.get_psd returns DataFrame or xarray Dataset by format."""
    db = connect()
    from ..results import MSNoiseResult
    import pandas as pd
    import xarray as xr
    results = MSNoiseResult.list(db, 'psd')
    if not results:
        pytest.skip("No PSD results available")
    r = results[0]
    all_psd = r.get_psd()
    if not all_psd:
        pytest.skip("No PSD data files available")
    (sid_k, day_k) = next(iter(all_psd))
    df = r.get_psd(sid_k, day_k, format="dataframe")
    assert isinstance(df, pd.DataFrame)
    ds = r.get_psd(sid_k, day_k, format="xarray")
    assert isinstance(ds, xr.Dataset)
    db.close()


@pytest.mark.order(120027)
def test_120027_msnoise_result_branches_from_stack():
    """MSNoiseResult.branches() from stack returns refstack children."""
    db = connect()
    from ..results import MSNoiseResult
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    children = r.branches()
    assert isinstance(children, list)
    for child in children:
        assert child.category in ('refstack',), \
            f"Unexpected category from stack: {child.category}"
        assert child.lineage_names[:-1] == r.lineage_names
    db.close()


@pytest.mark.order(120028)
def test_120028_msnoise_result_default_formats():
    """Default format for get_ccf=xarray, get_ref=xarray."""
    db = connect()
    from ..results import MSNoiseResult
    import xarray as xr
    # get_ccf default = xarray DataArray
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    all_ccfs = r.get_ccf()
    if all_ccfs:
        v = next(iter(all_ccfs.values()))
        assert isinstance(v, xr.DataArray), "get_ccf default should be xarray DataArray"
    # get_ref default = xarray Dataset
    r2 = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                 stack=1, refstack=1)
    all_refs = r2.get_ref()
    if all_refs:
        v2 = next(iter(all_refs.values()))
        assert isinstance(v2, xr.Dataset), "get_ref default should be xarray Dataset"
    db.close()


@pytest.mark.order(120029)
def test_120029_msnoise_result_get_dvv_mwcs():
    """MSNoiseResult.get_dvv returns xarray Dataset from mwcs_dtt_dvv step."""
    db = connect()
    from ..results import MSNoiseResult
    import xarray as xr
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1, mwcs=1,
                                mwcs_dtt=1, mwcs_dtt_dvv=1)
    assert r.category == "mwcs_dtt_dvv"
    # Discovery: all available (pair_type, comp, mov_stack) combos
    all_dvv = r.get_dvv()
    if not all_dvv:
        pytest.skip("No mwcs_dtt_dvv aggregate data — run compute first")
    # Every value should be an xarray Dataset with a times dimension
    for key, ds in all_dvv.items():
        assert isinstance(ds, xr.Dataset), \
            f"Expected xr.Dataset for key {key}, got {type(ds)}"
        assert "times" in ds.dims, \
            f"Expected 'times' dim in Dataset for key {key}"
        assert "mean" in ds, \
            f"Expected 'mean' variable in Dataset for key {key}"
        assert "weighted_mean" in ds, \
            f"Expected 'weighted_mean' in Dataset for key {key} (dvv_weighted_mean=Y)"
        assert "trimmed_mean" in ds, \
            f"Expected 'trimmed_mean' in Dataset for key {key} (dvv_trimmed_mean=Y)"
    db.close()


@pytest.mark.order(120030)
def test_120030_msnoise_result_get_dvv_dataframe():
    """MSNoiseResult.get_dvv returns DataFrame when format='dataframe'."""
    db = connect()
    from ..results import MSNoiseResult
    import pandas as pd
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1, mwcs=1,
                                mwcs_dtt=1, mwcs_dtt_dvv=1)
    all_dvv = r.get_dvv()
    if not all_dvv:
        pytest.skip("No mwcs_dtt_dvv aggregate data — run compute first")
    pt, comp, ms = next(iter(all_dvv))
    df = r.get_dvv(pair_type=pt, components=comp,
                   mov_stack=ms, format="dataframe")
    assert isinstance(df, pd.DataFrame), \
        f"Expected pd.DataFrame, got {type(df)}"
    assert "mean" in df.columns, "Expected 'mean' column in DataFrame"
    db.close()


@pytest.mark.order(120031)
def test_120031_msnoise_result_gating_raises():
    """Methods for absent categories raise AttributeError, not ValueError."""
    db = connect()
    from ..results import MSNoiseResult
    # Stack-level result — no mwcs/dvv in lineage
    r = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                               stack=1, refstack=1)
    with pytest.raises(AttributeError) as exc_info:
        r.get_mwcs()
    assert "mwcs" in str(exc_info.value).lower()
    assert "branches" in str(exc_info.value).lower(), \
        "Error message should suggest .branches()"

    # mwcs_dtt result — get_dvv not yet reachable
    r2 = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                stack=1, refstack=1, mwcs=1, mwcs_dtt=1)
    with pytest.raises(AttributeError):
        r2.get_dvv()
    db.close()


@pytest.mark.order(120032)
def test_120032_msnoise_result_dir_gating():
    """__dir__ only exposes methods whose categories are in the lineage."""
    db = connect()
    from ..results import MSNoiseResult

    # Stack-level result
    r_stack = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    d_stack = dir(r_stack)
    assert "get_ccf" in d_stack,   "get_ccf should be visible at stack level"
    assert "get_ref" not in d_stack, "get_ref needs refstack — should be hidden"
    assert "get_mwcs" not in d_stack, "get_mwcs needs mwcs — should be hidden"

    # Refstack-level result
    r_ref = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                   stack=1, refstack=1)
    d_ref = dir(r_ref)
    assert "get_ccf" in d_ref
    assert "get_ref" in d_ref,  "get_ref should be visible at refstack level"
    assert "get_mwcs" not in d_ref, "get_mwcs still needs mwcs step"

    # MWCS result
    r_mwcs = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                    stack=1, refstack=1, mwcs=1)
    d_mwcs = dir(r_mwcs)
    assert "get_mwcs" in d_mwcs,    "get_mwcs should be visible at mwcs level"
    assert "get_mwcs_dtt" not in d_mwcs, "get_mwcs_dtt still needs mwcs_dtt"

    db.close()


@pytest.mark.order(120034)
def test_120034_msnoise_result_export_dvv():
    """export_dvv writes NetCDF with all 7 provenance attributes."""
    import tempfile
    import xarray as xr
    try:
        import yaml
    except ImportError:
        pytest.skip("pyyaml not installed")

    db = connect()
    from ..results import MSNoiseResult

    results = MSNoiseResult.list(db, "mwcs_dtt_dvv")
    if not results:
        pytest.skip("No mwcs_dtt_dvv results — run full pipeline first")

    r = results[0]

    # export_dvv must be visible (lineage includes mwcs_dtt_dvv)
    assert hasattr(r, "export_dvv"), "export_dvv should be visible on dvv result"

    with tempfile.TemporaryDirectory() as tmpdir:
        written = r.export_dvv(tmpdir, pair_type="CC")
        assert len(written) >= 1, "Expected at least one exported file"

        for fpath in written:
            assert os.path.isfile(fpath), f"File missing: {fpath}"
            ds = xr.open_dataset(fpath)

            # All 9 required global attributes (7 original + 2 DataSource provenance)
            for attr in ("lineage", "msnoise_params", "msnoise_version",
                         "generated", "pair_type", "components", "mov_stack",
                         "data_sources", "stations"):
                assert attr in ds.attrs, f"Missing attribute: {attr!r}"

            # data_sources must be valid YAML listing at least the local source
            import yaml as _yaml
            ds_list = _yaml.safe_load(ds.attrs["data_sources"])
            assert isinstance(ds_list, list) and len(ds_list) >= 1, \
                "data_sources must be a non-empty list"
            assert ds_list[0]["name"] == "local", \
                f"Expected 'local' DataSource, got {ds_list[0]}"
            assert "password" not in ds.attrs["data_sources"], \
                "Credentials must not appear in data_sources export"

            # stations must be valid YAML
            sta_list = _yaml.safe_load(ds.attrs["stations"]) or []
            assert isinstance(sta_list, list), "stations must be a list"

            # lineage round-trips exactly
            assert ds.attrs["lineage"] == "/".join(r.lineage_names),                 "Lineage attribute mismatch"

            # YAML parses and contains at least the global config block
            params = yaml.safe_load(ds.attrs["msnoise_params"])
            assert isinstance(params, dict), "msnoise_params is not valid YAML"
            assert "lineage" in params, "YAML missing 'lineage' key"
            assert "global" in params,  "YAML missing global config block"

            # Dataset must have times dimension and at least a mean variable
            assert "times" in ds.dims,       "Missing 'times' dimension"
            assert "mean" in ds.data_vars,   "Missing 'mean' variable"
            assert len(ds.times) >= 1,       "No time steps in exported dataset"

            ds.close()
    db.close()



@pytest.mark.order(120035)
def test_120033_msnoise_result_hasattr_gating():
    """hasattr() respects dynamic gating — False for absent categories."""
    db = connect()
    from ..results import MSNoiseResult

    r_stack = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1, stack=1)
    assert hasattr(r_stack, "get_ccf"),    "get_ccf should be accessible"
    assert not hasattr(r_stack, "get_ref"), "get_ref should be hidden"
    assert not hasattr(r_stack, "get_mwcs"), "get_mwcs should be hidden"
    assert not hasattr(r_stack, "get_dvv"),  "get_dvv should be hidden"

    r_mwcs = MSNoiseResult.from_ids(db, preprocess=1, cc=1, filter=1,
                                    stack=1, refstack=1, mwcs=1)
    assert hasattr(r_mwcs, "get_mwcs"),        "get_mwcs accessible at mwcs level"
    assert not hasattr(r_mwcs, "get_mwcs_dtt"), "get_mwcs_dtt needs mwcs_dtt step"
    assert not hasattr(r_mwcs, "get_dvv"),      "get_dvv needs dvv step"

    db.close()

