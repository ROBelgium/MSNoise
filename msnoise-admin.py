from flask import Flask, redirect
from flask.ext.admin import Admin, BaseView, expose
from flask import render_template as template
import flask, os, sys, time, json, socket
import psutil
from subprocess import Popen, PIPE
from flask.ext.admin.contrib.sqla import ModelView
from flask import flash
import os
import sys
from wtforms.validators import ValidationError
from flask.ext.admin.actions import action
from flask.ext.admin.babel import gettext, ngettext, lazy_gettext
import webbrowser

from database_tools import *
from default import *
from msnoise_table_def import *

class GenericView(BaseView):
    name = "MSNoise"
    page = ""
    @expose('/')
    def index(self):
        return self.render('admin/%s.html'%self.page)

class FilterView(ModelView):
    view_title = "Filter Configuration"
    name = "filter"
    # Disable model creation
    def mwcs_low(form, field):
        if field.data <= form.data['low']:
            raise ValidationError("'mwcs_low' should be (at least slightly) greater than 'low'")
            
    def mwcs_high(form, field):
        if field.data <= form.data['mwcs_low']:
            raise ValidationError("'mwcs_high' should be greater than 'mwcs_low'")
    
    def high(form, field):
        if field.data <= form.data['mwcs_high']:
            raise ValidationError("'high' should be (at least slightly) greater than 'mwcs_high'")
    
    def mwcs_step(form, field):
        if field.data > form.data['mwcs_wlen']:
            raise ValidationError("'mwcs_step' should be smaller or equal to 'mwcs_wlen'")
    
    form_args = dict(
        mwcs_low=dict(validators=[mwcs_low]),
        mwcs_high=dict(validators=[mwcs_high]),
        high=dict(validators=[high]),
        mwcs_step=dict(validators=[mwcs_step]),
    )
    
    column_list = ('ref','low', 'mwcs_low', 'mwcs_high', 'high',
                   'rms_threshold', 'mwcs_wlen', 'mwcs_step', 'used')
    form_columns = ('low', 'mwcs_low', 'mwcs_high', 'high',
                   'rms_threshold', 'mwcs_wlen', 'mwcs_step', 'used')
    
    def __init__(self, session, **kwargs):
        # You can pass name and other parameters if you want to
        super(FilterView, self).__init__(Filter, session, **kwargs)
        
    @action('used',
            lazy_gettext('Toggle Used'),
            lazy_gettext('Are you sure you want to update selected models?'))
    def used(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            if s.used:
                s.used = False
            else:
                s.used = True
        self.session.commit()
        return  

class StationView(ModelView):
    view_title = "Station Configuration"
    column_filters = ('net', 'used')
    
    def __init__(self, session, **kwargs):
        super(StationView, self).__init__(Station, session, **kwargs)
    
    @action('used',
            lazy_gettext('Toggle Used'),
            lazy_gettext('Are you sure you want to update selected models?'))
    def used(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            if s.used:
                s.used = False
            else:
                s.used = True
        self.session.commit()
        return  


class DataAvailabilityView(ModelView):
    view_title = "Data Availability"
    can_create = False
    can_delete = False
    can_edit = True
    column_filters = ('net', 'sta', 'comp','data_duration','gaps_duration','samplerate','flag')
    def __init__(self, session, **kwargs):
        super(DataAvailabilityView, self).__init__(DataAvailability, session, **kwargs)
    
    @action('modified',
            lazy_gettext('Mark as (M)odified'),
            lazy_gettext('Are you sure you want to update selected models?'))
    def modified(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        count = 0
        for s in query.all():
            s.flag = 'M'
            count += 1
        self.session.commit()
        flash(ngettext('Model was successfully flagged (M)odified.',
               '%(count)s models were successfully flagged (M)odified.',
               count,
               count=count))
        return     
    

class JobView(ModelView):
    view_title = "Jobs"
    can_create = False
    can_delete = True
    can_edit = True
    column_filters = ('pair','type','flag')
    
    def __init__(self, session, **kwargs):
        super(JobView, self).__init__(Job, session, **kwargs)
    
    @action('todo',
        lazy_gettext('Mark as (T)odo'),
        lazy_gettext('Are you sure you want to update selected models?'))
    def todo(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            s.flag = 'T'
        self.session.commit()
        return
    
    @action('done',
        lazy_gettext('Mark as (D)one'),
        lazy_gettext('Are you sure you want to update selected models?'))
    def done(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            s.flag = 'D'
        self.session.commit()
        return
    
    @action('deletetype',
        lazy_gettext('Delete all Jobs of the same "Type"'),
        lazy_gettext('Are you sure you want to delete all those models?'))
    def deletetype(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            type_to_delete = s.type
        self.get_query().filter(Job.type == type_to_delete).delete()
        self.session.commit()
        return
    
    @action('massTodo',
        lazy_gettext('Mark all Jobs of the same Type as (T)odo'),
        lazy_gettext('Are you sure you want to update all those models?'))
    def massTodo(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            type_to_delete = s.type
        
        for s in self.get_query().filter(Job.type == type_to_delete).all():
            s.flag = 'T'
        self.session.commit()
        return
    
    
    
    

class ConfigView(ModelView):
    # Disable model creation
    view_title = "MSNoise General Configuration"
    def no_root_allowed(form, field):
        if field.data == 'root':
            raise ValidationError('"root" is not allowed')
 
    # inline_models = (Config,)
    form_args = dict(
        value=dict(validators=[no_root_allowed])
    )
    can_create = False
    can_delete = False
    page_size = 50
    # Override displayed fields
    column_list = ('name', 'value')

    def __init__(self, session, **kwargs):
        # You can pass name and other parameters if you want to
        super(ConfigView, self).__init__(Config, session, **kwargs)

class PSUtils(BaseView):
    name = "MSNoise"
    view_title = "Some Real Time Statistics"
    @expose('/')
    def index(self):
        cpucount = len(psutil.cpu_percent(0,True))
        return self.render('admin/psutils.html',cpucount=cpucount)

class ResultPlotter(BaseView):
    name = "MSNoise"
    view_title = "Result Plotter"
    @expose('/')
    def index(self):
        return self.render('admin/results.html')

class InterferogramPlotter(BaseView):
    name = "MSNoise"
    view_title = "Interferogram Plotter"
    @expose('/')
    def index(self):
        return self.render('admin/interferogram.html')


class DataAvailabilityPlot(BaseView):
    name = "MSNoise"
    view_title = "Data Availability"
    @expose('/')
    def index(self):
        return self.render('admin/data_availability.html')


class BugReport(BaseView):
    name = "MSNoise"
    view_title = "BugReport"
    @expose('/')
    def index(self):
        return self.render('admin/bugreport.html')

app = Flask(__name__, template_folder='templates')
app.secret_key = 'why would I tell you my secret key?'

@app.route('/admin/rawps')
def rawps():
    diskused = 0
    disktotal = 0
    for i in psutil.disk_partitions():
        try:
            x = psutil.disk_usage(i.mountpoint)
            diskused += x.used
            disktotal += x.total
        except OSError:
            pass
    data = {
        'uptime':	time.time() - psutil.BOOT_TIME,
        'fqdn':		socket.getfqdn(),
        'cpuusage':	psutil.cpu_percent(0,True),
        'cpu_count': len(psutil.cpu_percent(0,True)),
        'ramusage':	psutil.virtual_memory(),
        'diskio':	psutil.disk_io_counters(),
        'diskusage':	[diskused, disktotal],
        'netio':	psutil.net_io_counters(),
        'swapusage':	psutil.swap_memory()
    }
    o = json.dumps(data)
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/networks.json')
def networksJSON():
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    data = {}
    networks = get_networks(db)
    for network in networks:
        stations = get_stations(db, net=network)
        data[network] = [s.sta for s in stations]
    o = json.dumps(data)
    db.close()
    return flask.Response(o, mimetype='application/json')

@app.route('/admin/filters.json')
def filtersJSON():
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    data = {}
    filters = get_filters(db, all=False)
    for f in filters:
        data[f.ref] = "%.2f - %.2f"%(f.low, f.high)
    db.close()
    o = json.dumps(data)
    db.close()
    return flask.Response(o, mimetype='application/json')

@app.route('/admin/components.json')
def componentsJSON():
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    components = get_components_to_compute(db)
    data = {}
    for i,c in enumerate(components):
        data[i] = c
    db.close()
    o = json.dumps(data)
    db.close()
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/pairs.json')
def pairs():
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    stations = ["%s.%s" % (s.net, s.sta) for s in get_stations(db,all=False)]
    pairs = itertools.combinations(stations, 2)
    output = []
    for pair in pairs:
        output.append("%s - %s"%(pair[0],pair[1]))
    o = json.dumps(output)
    db.close()
    return flask.Response(o, mimetype='application/json')

@app.route('/admin/bugreport.json')
def bugreporter():
    proc = Popen(
            ["python", os.path.join(MSNoise_PATH,"bugreport.py"), "-m"], stdout=PIPE, stderr=PIPE)
    stdout, stderr = proc.communicate()
    output = {}
    for i, line in enumerate(stdout.split('\n')):
        output[i] = line
    
    o = json.dumps(output)
    db.close()
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/data_availability.json',methods=['GET','POST'])
def dataAvail():
    data = flask.request.get_json()
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    data = get_data_availability(db, net=data['net'], sta=data['sta'],comp='HHZ')
    o = {'dates':[o.starttime.strftime('%Y-%m-%d') for o in data]}
    db.close()
    o['result']='ok'
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')

@app.route('/admin/all_results.json',methods=['POST'])
def allresults():
    data = flask.request.get_json(force=True)
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    station1, station2 = flask.request.json['pair'].split("-")
    station1, station2 = data['pair'].replace('.','_').split(' - ')
    filterid = int(data['filter'])
    components = data['component']
    format = data['format']
    start, end, dates = build_ref_datelist(db)
    i, result = get_results(db,station1, station2, filterid, components, dates, format=format)
    
    data = {}
    if format == 'stack':
        if i != 0:
            maxlag = float(get_config(db, 'maxlag'))
            data['x'] = np.linspace(-maxlag, maxlag, get_maxlag_samples(db)).tolist()
            data['y'] = result.tolist()
            data["info"] = "Data OK"
        else:
            data["info"] = "No Data"
    else:
        if i != 0:
            x = []
            for y in range(len(dates)):
                r = result[y]
                r = np.nan_to_num(r)
                x.append( r.tolist() )
            data["image"] = x
            data["nx"] = len(r)
            data["ny"] = len(x)
            data["info"] = "Data OK"
        else:
            data["info"] = "No Data"
    o = json.dumps(data)
    return flask.Response(o, mimetype='application/json')

@app.route('/admin/new_jobs_TRIG.json')
def new_jobsTRIG():
    from s02new_jobs import new_jobs
    count = new_jobs()
    global db
    db.flush()
    db.commit()

    o = {}
    o['count']=count
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')
    

@app.route('/admin/jobs_list.json')
def joblists():
    jobtype = flask.request.args['type']
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    data = get_job_types(db,jobtype)
    db.close()
    o = {'T':0,'I':0,'D':0}
    for count, flag in data:
        o[flag] = count

    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')
 
@app.route('/admin/data_availability_flags.json')
def DA_flags():
    db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))
    data = count_data_availability_flags(db)
    db.close()
    o = {'N':0,'M':0,'A':0}
    for count, flag in data:
        o[flag] = count
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')
 
 


# Flask views
@app.route('/')
def index():
    return redirect("/admin/", code=302)


global db
db = connect(inifile=os.path.join(MSNoise_PATH,'db.ini'))

admin = Admin(app)
admin.name = "MSNoise"
admin.add_view(StationView(db,endpoint='stations', category='Configuration'))
admin.add_view(FilterView(db,endpoint='filters', category='Configuration'))
admin.add_view(ConfigView(db,endpoint='config', category='Configuration'))

admin.add_view(DataAvailabilityView(db,endpoint='data_availability',category='Analytics'))
admin.add_view(DataAvailabilityPlot(endpoint='data_availability_plot',category='Analytics'))
admin.add_view(JobView(db,endpoint='jobs',category='Analytics'))
admin.add_view(PSUtils(name='PS Utils', endpoint='psutils', category='Analytics'))
admin.add_view(BugReport(name='Bug Report', endpoint='bugreport', category='Analytics'))

admin.add_view(ResultPlotter(endpoint='results',category='Results'))
admin.add_view(InterferogramPlotter(endpoint='interferogram',category='Results'))
a = GenericView(endpoint='about',category='Help')
a.page = 'about'
a.name = 'About'
admin.add_view(a)

app.run(host='0.0.0.0', debug=True)
