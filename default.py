from collections import OrderedDict
default = OrderedDict()
default['data_folder'] = ["Data Folder",'']
default['output_folder'] = ["CC Output Folder",'CROSS_CORRELATIONS']
default['data_structure'] = ["Data Structure [SDS]/BUD/IDDS",'SDS']
default['network'] = ["Network to analyse [*]",'*']
default['channels'] = ["Channels need to match the value (ex: [*], *Z, BH*, HHZ,...)",'*']

default['startdate'] = ["Start Date to process: [1970-01-01]='since beginning of the archive'","1970-01-01"]
default['enddate'] = ["End Date to process: [2100-01-01]='No end'","2100-01-01"]

default['analysis_duration'] = ["Duration of the Analysis (total in seconds : 3600, [86400])",'86400']
default['cc_sampling_rate'] = ["Sampling Rate for the CrossCorrelation [20.0]",'20.0']
default['resampling_method'] = ["Resampling method [Resample]/Decimate",'Resample']
default['decimation_factor'] = ["If Resampling mether=Decimate, decimation factor [5]",'5']
default['preprocess_lowpass'] = ["Preprocessing Low-pass value in Hz [8.0]",'8.0']
default['preprocess_highpass'] = ["Preprocessing High-pass value in Hz [0.01]",'0.01']


default['maxlag'] = ["Maximum lag (in seconds) [120.0]",'120.']
default['corr_duration'] = ["Data windows to correlate (in seconds) [1800.]",'1800.']
default['windsorizing'] = ["Windsorizing at N time RMS (in unit), 0 disables windsorizing [3]",'3']

default['crondays'] = ["Number of days to monitors with cron [-1]",'-1']


default['ZZ'] = ["Compute ZZ correlation [Y]/N",'Y']
default['ZR'] = ["Compute ZR correlation [Y]/N",'Y']
default['ZT'] = ["Compute ZT correlation [Y]/N",'Y']
default['RZ'] = ["Compute RZ correlation [Y]/N",'Y']
default['RR'] = ["Compute RR correlation [Y]/N",'Y']
default['RT'] = ["Compute RT correlation [Y]/N",'Y']
default['TZ'] = ["Compute TZ correlation [Y]/N",'Y']
default['TR'] = ["Compute TR correlation [Y]/N",'Y']
default['TT'] = ["Compute TT correlation [Y]/N",'Y']

default['autocorr'] = ["Compute Auto correlation [Y]/N",'N']
default['PAZ'] = ["Correct instrumental responce from paz [Y]/N",'N']
default['keep_all'] = ["Keep all 30 seconds cross-corr [Y]/N",'N']
default['keep_days'] = ["Keep all daily cross-corr [Y]/N",'Y']

default['ref_begin'] = ["Beginning or REF stacks. Can be absolute (2012-01-01) or relative (-100) days",'-100']
default['ref_end'] = ["End or REF stacks. Same as ref_begin",'0']

default['mov_stack'] = ["Number of days to stack for the Moving-window stacks ([5]= [day-4:day]), can be a comma-separated list 2,5,10","5"]

default['export_format'] = ["Export stacks in which format(s) ? SAC/MSEED/[BOTH]","MSEED"]
default['sac_format'] = ["Format for SAC stacks ? [doublets]/clarke","doublets"]

