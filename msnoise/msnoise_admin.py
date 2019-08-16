"""MSNoise Admin is a web interface that helps the user define the
configuration for all the processing steps, it also allows configuring the
stations and filters to be used in the processes. It gives a view on the
database tables.

To start the admin:

.. code-block:: sh

    $ msnoise admin

Which, by default, starts a web server listening on all interfaces on port
5000. This can be overriden by passing parameters to the command, e.g. for port
5099:

.. code-block:: sh

    $ msnoise admin -p 5099


Next step is to open a web browser and open the ip address of the machine,
by default on the current machine, it'll be http://localhost:5000/ or
http://127.0.0.1:5000/.

.. image:: ../.static/msnoise_admin_home.png
    :align: center


The top level menu shows four items:

Home
----

The index page shows

* The project location and its database
* Stats of the Data Availability, the CC, STACK, MWCS and DTT jobs.
* Quick action buttons for resetting or deleting jobs.

The name and the logo of the page can be overriden by setting an environment
variable with a name and the HTML tag of the logo image:

.. code:: sh
    
    set msnoise_brand="ROB|<img src='http://www.seismologie.be/img/oma/ROB-logo.svg' width=200 height=200>"
    
and then starting msnoise admin:

.. image:: ../.static/branding.png
    :align: center


Configuration
--------------

Station
~~~~~~~

Stations appear as a table and are editable.

Stations are defined as:

.. autoclass:: Station


Filter
~~~~~~

Filters appear as a table and are editable. The filter parameters are validated
before submission, so no errors should happen. Note: by default, the `used`
parameter is set to `False`, **don't forget to change it!**

Filters are defined as:

.. autoclass:: Filter

Config
~~~~~~

All configuration bits appear as a table and are editable. When editing one
configuration item, the Edit tab shows extra information about the
parameter, where it is used and its default value. Most of the configuration
bits are case-sensitive!

Example view:

.. image:: ../.static/msnoise_admin_config.png
    :align: center

The table below repeats this

.. include:: ../defaults.rst

Database
--------

Data Availability
~~~~~~~~~~~~~~~~~

Gives a view of the `data_availability` table. Allows to bulk edit/select rows.
Its main goal is to check that the `scan_archive` procedure has successfully
managed to list all files from one's archive.

Jobs
~~~~

Gives a view of the `jobs` table. Allows to bulk edit/select rows. Its main
goal is to check the `new_jobs` or any other workflow step (or Plugins)
successfully inserted/updated jobs.


Help
----

About
~~~~~

Shows some links and information about the package. Mostly the information
present on the github readme file.

Bug Report
~~~~~~~~~~

Web view of the `msnoise bugreport -m`, allows viewing if all required python
modules are properly installed and available for MSNoise.


"""

import json
from io import BytesIO

import flask
import jinja2
import markdown
from flask import Flask, redirect, request, render_template
from flask import Markup
from flask import flash
from flask_admin import Admin, BaseView, expose
from flask_admin.actions import action
from flask_admin.babel import ngettext, lazy_gettext
from flask_admin.contrib.sqla import ModelView
from flask_admin.model import typefmt
from wtforms.validators import ValidationError
from wtforms.fields import SelectField, StringField
from wtforms.utils import unset_value
from flask_wtf import Form

from .api import *
from .default import default

from .msnoise_table_def import Filter, Job, Station, Config, DataAvailability


class GenericView(BaseView):
    name = "MSNoise"
    page = "index"

    @expose('/')
    def index(self):
        return self.render('admin/%s.html'%self.page, msnoise_project="test")


class FilterView(ModelView):
    view_title = "Filter Configuration"
    name = "filter"

    def mwcs_low(form, field):
        if field.data < form.data['low']:
            raise ValidationError("'mwcs_low' should be greater or equal to"
                                  " 'low'")

    def mwcs_high(form, field):
        if field.data <= form.data['mwcs_low']:
            raise ValidationError("'mwcs_high' should be greater than"
                                  " 'mwcs_low'")

    def high(form, field):
        if field.data < form.data['mwcs_high']:
            raise ValidationError("'high' should be greater or equal than"
                                  " 'mwcs_high'")

    def mwcs_step(form, field):
        if field.data > form.data['mwcs_wlen']:
            raise ValidationError("'mwcs_step' should be smaller or equal to"
                                  " 'mwcs_wlen'")
    
    form_args = dict(
        mwcs_low=dict(validators=[mwcs_low]),
        mwcs_high=dict(validators=[mwcs_high]),
        high=dict(validators=[high]),
        mwcs_step=dict(validators=[mwcs_step]),
    )
    
    column_list = ('ref', 'low', 'mwcs_low', 'mwcs_high', 'high',
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


MY_DEFAULT_FORMATTERS = dict(typefmt.BASE_FORMATTERS)
MY_DEFAULT_FORMATTERS.update({
        type(None): typefmt.null_formatter,
    })


class StationView(ModelView):
    view_title = "Station Configuration"
    column_filters = ('net', 'used')
    column_type_formatters = MY_DEFAULT_FORMATTERS
    can_set_page_size = True

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
    can_delete = True
    can_edit = True
    page_size = 100
    can_set_page_size = True
    column_filters = ('net', 'sta', 'comp', 'data_duration', 'gaps_duration',
                      'samplerate', 'flag')

    def __init__(self, session, **kwargs):
        super(DataAvailabilityView, self).__init__(DataAvailability, session,
                                                   **kwargs)
    
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
              count, count=count))
        return
    

class JobView(ModelView):
    view_title = "Jobs"
    can_create = False
    can_delete = True
    can_edit = True
    column_filters = ('pair', 'jobtype', 'flag')
    page_size = 100
    edit_modal = True
    can_set_page_size = True

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
            type_to_delete = s.jobtype
        self.get_query().filter(Job.jobtype == type_to_delete).delete()
        self.session.commit()
        return
    
    @action('massTodo',
            lazy_gettext('Mark all Jobs of the same Type as (T)odo'),
            lazy_gettext('Are you sure you want to update all those models?'))
    def massTodo(self, ids):
        model_pk = getattr(self.model, self._primary_key)
        query = self.get_query().filter(model_pk.in_(ids))
        for s in query.all():
            type_to_delete = s.jobtype
        
        for s in self.get_query().filter(Job.jobtype == type_to_delete).all():
            s.flag = 'T'
        self.session.commit()
        return



class ConfigView(ModelView):
    # Disable model creation
    edit_template = 'admin/model/edit-config.html'

    view_title = "MSNoise General Configuration"
    can_create = False
    can_delete = False
    page_size = 100
    can_set_page_size = True
    column_sortable_list = ["name",]

    # Override displayed fields
    column_list = ('name', 'value')

    def __init__(self, session, **kwargs):
        super(ConfigView, self).__init__(Config, session, **kwargs)

    # def edit_form(self, obj=None):
    #     form = super(ModelView, self).edit_form(obj)
    #     nf = SelectField("Value", choices=[("E:\\", "YREAHHHH"), ("2", "KO")])
    #     nf = nf.bind(form, "value")
    #     nf.data = "E:\\"
    #     form._fields["value"] = nf
    #     return form

    @expose('/edit/', methods=['GET', 'POST'])
    def edit_view(self):
        id = request.args.get('id')
        helpstring = default[id][0]
        helpstring = Markup(markdown.markdown(helpstring))
        self._template_args['helpstring'] = helpstring
        self._template_args['helpstringdefault'] = default[id][1]
        return super(ConfigView, self).edit_view()

    def on_model_change(self, form, model, is_created):
        if form.data['value']:
            model.value = form.data['value'].strip()
        else:
            model.value = ''


def getitem(obj, item, default):
    if item not in obj:
        return default
    else:
        return obj[item]


def select_filter():
    db = connect()
    filters = []
    query = get_filters(db, all=False)
    for f in query:
        filters.append({'optid': f.ref,
                        'text': "%.2f - %.2f" % (f.low, f.high)})
    db.close()
    return filters


def select_pair():
    db = connect()
    stations = ["%s.%s" % (s.net, s.sta) for s in get_stations(db, all=False)]
    pairs = itertools.combinations(stations, 2)
    output = []
    i = 0
    for pair in pairs:
        output.append({'optid': i,
                       'text': "%s - %s" % (pair[0], pair[1])})
        i += 1
    db.close()
    return output


class ResultPlotter(BaseView):
    name = "MSNoise"
    view_title = "Result Plotter"

    @expose('/')
    def index(self):
        args = flask.request.args

        filters = select_filter()
        pairs = select_pair()

        # Get all the form arguments in the url with defaults
        filter = int(getitem(args, 'filter', 1))
        pair = int(getitem(args, 'pair', 0))
        component = getitem(args, 'component', 'ZZ')
        format = getitem(args, 'format', 'stack')

        db = connect()
        params = get_params(db)
        station1, station2 = pairs[pair]['text'].replace('.', '_').split(' - ')
        start, end, dates = build_movstack_datelist(db)
        print(station1, station2, filter, component, dates, format)
        i, result = get_results(db,station1, station2, filter, component, dates,
                                format=format, params=params)

        if format == 'stack':
            if i != 0:
                x = get_t_axis(db)
                y = result
        db.close()
        # print(result)
        # fig = figure(title=pairs[pair]['text'], plot_width=1000)
        # fig.line(x, y, line_width=2)
        #
        # plot_resources = RESOURCES.render(
        #     js_raw=CDN.js_raw,
        #     css_raw=CDN.css_raw,
        #     js_files=CDN.js_files,
        #     css_files=CDN.css_files,
        # )
        #
        # script, div = components(fig, INLINE)
        import matplotlib.pyplot as plt
        plt.figure(figsize=(12,5))
        plt.plot(x, y)

        from io import BytesIO
        figfile = BytesIO()
        plt.savefig(figfile, format='png')
        figfile.seek(0)  # rewind to beginning of file
        import base64
        # figdata_png = base64.b64encode(figfile.read())
        figdata_png = base64.b64encode(figfile.getvalue())
        result= figdata_png

        plot_div = ""

        return self.render(
            'admin/results.html',
            plot_script="", plot_div=plot_div, plot_resources="",
            filter_list=filters,
            pair_list=pairs, result=result.decode('utf8')
        )


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

app = Flask(__name__)
app.secret_key = 'why would I tell you my secret key?'


@app.route('/admin/networks.json')
def networksJSON():
    db = connect()
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
    db = connect()
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
    db = connect()
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
    db = connect()
    stations = ["%s.%s" % (s.net, s.sta) for s in get_stations(db, all=False)]
    pairs = itertools.combinations(stations, 2)
    output = []
    for pair in pairs:
        output.append("%s - %s" % (pair[0], pair[1]))
    o = json.dumps(output)
    db.close()
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/bugreport.json')
def bugreporter():
    from .bugreport import main
    output = main(modules=True, show=False)
    o = json.dumps(output)
    db.close()
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/data_availability.json',methods=['GET','POST'])
def dataAvail():
    data = flask.request.get_json()
    db = connect()
    data = get_data_availability(db, net=data['net'], sta=data['sta'],
                                 comp='HHZ')
    o = {'dates': [o.starttime.strftime('%Y-%m-%d') for o in data]}
    db.close()
    o['result'] = 'ok'
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/all_results.json',methods=['POST'])
def allresults():
    data = flask.request.get_json(force=True)
    db = connect()
    station1, station2 = data['pair'].replace('.','_').split(' - ')
    filterid = int(data['filter'])
    components = data['component']
    format = data['format']
    start, end, dates = build_ref_datelist(db)
    i, result = get_results(db,station1, station2, filterid, components, dates,
                            format=format)
    
    data = {}
    if format == 'stack':
        if i != 0:
            maxlag = float(get_config(db, 'maxlag'))
            data['x'] = np.linspace(-maxlag, maxlag, get_maxlag_samples(db)).\
                tolist()
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
                x.append(r.tolist())
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
    from .s02new_jobs import main
    count = main()
    global db
    db.flush()
    db.commit()
    o = {}
    o['count'] = count
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')
    

@app.route('/admin/jobs_list.json')
def joblists():
    jobtype = flask.request.args['type']
    db = connect()
    data = get_job_types(db, jobtype)
    db.close()
    o = {'T': 0, 'I': 0, 'D': 0}
    for count, flag in data:
        o[flag] = count
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/resetjobs.json')
def resetjobs():
    jobtype = flask.request.args['type']
    alljobs = flask.request.args['all']
    from msnoise.api import reset_jobs
    reset_jobs(db, jobtype, alljobs=alljobs)
    o = {}
    o["Done"] = "Jobs reset: Done."
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')


@app.route('/admin/data_availability.png')
def DA_PNG():
    from .plots.data_availability import main
    output = BytesIO()
    main(show=False, outfile=output)
    output.seek(0)
    return flask.Response(output.read(), mimetype='image/png')


@app.route('/admin/data_availability_flags.json')
def DA_flags():
    db = connect()
    data = count_data_availability_flags(db)
    db.close()
    o = {'N': 0, 'M': 0, 'A': 0}
    for count, flag in data:
        o[flag] = count
    o = json.dumps(o)
    return flask.Response(o, mimetype='application/json')


def shutdown_server():
    func = request.environ.get('werkzeug.server.shutdown')
    if func is None:
        raise RuntimeError('Not running with the Werkzeug Server')
    func()


@app.route('/shutdown', methods=['GET'])
def shutdown():
    shutdown_server()
    return 'Server shutting down...'


# Flask views
@app.route('/')
def index():
    db = connect()
    return redirect("/admin/", code=302)


def main(port=5000):
    global db
    db = connect()
    plugins = get_config(db, "plugins")
    db.close()

    admin = Admin(app, template_mode='bootstrap2')
    
    if "msnoise_brand" in os.environ:
        tmp = eval(os.environ["msnoise_brand"])
        name, logo = tmp.split("|")
        admin.logo = logo
        admin.name = name
    else:
        admin.name = "MSNoise"
    admin.project_folder = os.getcwd()
    dbini = read_db_inifile()
    if dbini.prefix != "":
        prefix = "/%s_*" % dbini.prefix
    else:
        prefix = ""
    if dbini.tech == 1:
        database = "SQLite: %s%s" % (dbini.hostname, prefix)
    else:
        database = "MySQL: %s@%s:%s%s" % (dbini.username, dbini.hostname,
                                          dbini.database, prefix)
    admin.project_database = database

    jobtypes = ["CC", "STACK", "MWCS", "DTT"]
    template_folders = []
    if plugins:
        for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.jobtypes'):
            module_name = ep.module_name.split(".")[0]
            if module_name in plugins:
                tmp = ep.load()()
                for t in tmp:
                    jobtypes.append(t["name"])
        for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.templates'):
            module_name = ep.module_name.split(".")[0]
            if module_name in plugins:
                tmp = ep.load()()
                for t in tmp:
                    template_folders.append(t)

    app.jinja_loader = jinja2.ChoiceLoader([
        app.jinja_loader,
        jinja2.FileSystemLoader(template_folders),
        ])


    admin.add_view(StationView(db,endpoint='stations', category='Configuration'))
    admin.add_view(FilterView(db,endpoint='filters', category='Configuration'))
    admin.add_view(ConfigView(db,endpoint='config', category='Configuration'))

    admin.add_view(DataAvailabilityView(db,endpoint='data_availability',
                                        category='Database'))

    admin.add_view(JobView(db,endpoint='jobs',category='Database'))


    admin.add_view(DataAvailabilityPlot(endpoint='data_availability_plot',
                                        category='Results'))
    admin.add_view(ResultPlotter(endpoint='results',category='Results'))
    admin.add_view(InterferogramPlotter(endpoint='interferogram',
                                        category='Results'))

    if plugins:
        plugins = plugins.split(',')
        for ep in pkg_resources.iter_entry_points(group='msnoise.plugins.admin_view'):
            module_name = ep.module_name.split(".")[0]
            if module_name in plugins:
                admin.add_view(ep.load()(db))

    a = GenericView(endpoint='about',category='Help', name='About')
    a.page = "about"
    admin.add_view(a)
    admin.add_view(BugReport(name='Bug Report', endpoint='bugreport',
                             category='Help'))

    print("MSNoise admin will run on all interfaces by default")
    print("access it via the machine's IP address or")
    print("via http://127.0.0.1:%i when running locally."%port)
    app.run(host='0.0.0.0', port=port, debug=False, reloader_interval=1,
            threaded=True)
