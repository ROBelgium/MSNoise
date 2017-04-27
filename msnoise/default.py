from collections import OrderedDict
default = OrderedDict()
default['data_folder'] = ["Data Folder",'']
default['output_folder'] = ["CC Output Folder",'CROSS_CORRELATIONS']
default['data_structure'] = ["Either a predefined acronym [SDS]/BUD/IDDS,<br> "
                             "or /-separated path (e.g. NET/STA/YEAR/NET.STA.YEAR.DAY.MSEED).", 'SDS']
default['network'] = ["Network to analyse [*]",'*']
default['channels'] = ["Channels need to match the value (ex: [\*], \*Z, BH\*, HHZ,...)",'*']

default['startdate'] = ["Start Date to process: [1970-01-01]='since beginning of the archive'","1970-01-01"]
default['enddate'] = ["End Date to process: [2100-01-01]='No end'","2018-01-01"]

default['analysis_duration'] = ["Duration of the Analysis (total in seconds : 3600, [86400])",'86400']
default['cc_sampling_rate'] = ["Sampling Rate for the CrossCorrelation [20.0]",'20.0']
default['resampling_method'] = ["Resampling method Resample/Decimate/[Lanczos]",'Lanczos']
default['preprocess_lowpass'] = ["Preprocessing Low-pass value in Hz [8.0]",'8.0']
default['preprocess_highpass'] = ["Preprocessing High-pass value in Hz [0.01]",'0.01']

default['remove_response'] = ["Remove instrument response Y/[N]",'N']
default['response_format'] = ["Remove instrument file format [dataless]/inventory/paz/resp",'dataless']
default['response_path'] = ["Instrument correction file(s) location (path relative to db.ini), defaults to './inventory', i.e."
                            " a subfolder in the current project folder.<br>All files in that folder will be parsed.",'inventory']
default['response_prefilt'] = ["Remove instrument correction **pre-filter** (0.005, 0.006, 30.0, 35.0)",'(0.005, 0.006, 30.0, 35.0)']

default['maxlag'] = ["Maximum lag (in seconds) [120.0]",'120.']
default['corr_duration'] = ["Data windows to correlate (in seconds) [1800.]",'1800.']
default['overlap'] = ["Amount of overlap between data windows [0:1[ [0.]",'0.0']
default['windsorizing'] = ["Windsorizing at N time RMS , 0 disables windsorizing, -1 enables 1-bit normalization [3]",'3']
default['whitening'] = ["Whiten Traces before cross-correlation: All (except for autocorr), None, or only if components are different [A]/N/C",'A']

default['stack_method'] = ["Stack Method: Linear Mean or Phase Weighted Stack: [linear]/pws ",'linear']
default['pws_timegate'] = ["If stack_method='pws', width of the smoothing in seconds : 10.0 ",'10.0']
default['pws_power'] = ["If stack_method='pws', Power of the Weighting: 2.0 ",'2.0']

default['crondays'] = ["Number of days to monitor with cron [-1], must be a negative number", '-1']

default['components_to_compute'] = ["List (comma separated) [ZZ]", 'ZZ']

default['autocorr'] = ["Compute Auto correlation [Y]/N",'N']
default['keep_all'] = ["Keep all cross-corr (length: corr_duration) [Y]/N",'N']
default['keep_days'] = ["Keep all daily cross-corr [Y]/N",'Y']

default['ref_begin'] = ["Beginning or REF stacks. Can be absolute (2012-01-01) or relative (-100) days",'1970-01-01']
default['ref_end'] = ["End or REF stacks. Same as ref_begin",'2018-01-01']

default['mov_stack'] = ["Number of days to stack for the Moving-window stacks ([5]= [day-4:day]), can be a comma-separated list 1,2,5,10","5"]

default['export_format'] = ["Export stacks in which format(s) ? SAC/MSEED/[BOTH]","MSEED"]
default['sac_format'] = ["Format for SAC stacks ? [doublets]/clarke","doublets"]

default['dtt_lag'] = ["How is the lag window defined [dynamic]/static","static"]
default['dtt_v'] = ["If dtt_lag=dynamic: what velocity to use to avoid ballistic waves [1.0]km/s","1.0"]
default['dtt_minlag'] = ["If dtt_lag=static: min lag time","5.0"]
default['dtt_width'] = ["Width of the time lag window [30]s","30.0"]
default['dtt_sides'] = ["Which sides to use [both]/left/right","both"]
default['dtt_mincoh'] = ["Minimum coherence on dt measurement, MWCS points with values lower than that will **not** be used in the WLS","0.65"]
default['dtt_maxerr'] = ["Maximum error on dt measurement, MWCS points with values larger than that will **not** be used in the WLS","0.1"]
default['dtt_maxdt'] = ["Maximum dt values, MWCS points with values larger than that will **not** be used in the WLS","0.1"]

default['plugins'] = ["Comma separated list of plugin names. Plugins names should be importable Python modules.",""]
