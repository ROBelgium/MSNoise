"""This is the most important tool of the whole package. Configurator is a GUI
that helps the user define the configuration for all the processing steps, plus
configure the stations and filters to be used in the processes.
The GUI presents a tree-like structure which contains three main items: the
global MSNoise config, the Networks/Stations config and the Filters config.

.. image:: .static/configurator1.png

To run this script:

.. code-block:: sh

    $ msnoise config


Default Global Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. include:: defaults.rst


Network-Station Parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~

.. seealso :: :class:`msnoise.msnoise_table_def.Station`

Filter Parameters
~~~~~~~~~~~~~~~~~

.. seealso :: :class:`msnoise.msnoise_table_def.Filter`

"""


from traits.api import HasTraits, Str, List, Instance, Float, CInt, CBool, Enum
from traitsui.api import Item, View, TreeEditor, TreeNode, TableEditor
from traitsui.table_column import ObjectColumn
from traitsui.extras.checkbox_column import CheckboxColumn

from .default import *
from .api import *
from .msnoise_table_def import *


class StationColumn (ObjectColumn):
    def get_text_color(self, obj):
        return ['red', 'black'][obj.used]


class TFilter(HasTraits):
    ref = CInt
    low = Float
    mwcs_low = Float
    high = Float
    mwcs_high = Float
    rms_threshold = Float
    mwcs_wlen = Float
    mwcs_step = Float
    used = CBool
    traits_view = View('low', 'high', 'rms_threshold', 'used',
                       buttons=['OK', 'Cancel'])


class TStation(HasTraits):
    """ Defines a Station"""
    ref = CInt
    net = Str
    sta = Str
    X = Float
    Y = Float
    altitude = Float
    coordinates = Enum('UTM', 'DEG')
    instrument = Str
    used = CBool

    view = View('net', 'sta', 'X', 'Y', 'altitude',
                'coordinates', 'instrument', 'used')

station_editor = TableEditor(
    sortable=True,
    configurable=False,
    auto_size=True,
    columns=[CheckboxColumn(name='used', label='Used ?'),
             StationColumn(name='net', editable=False, width=0.24,
                           horizontal_alignment='left'),
             StationColumn(name='sta'),
             StationColumn(name='X'),
             StationColumn(name='Y'),
             StationColumn(name='altitude'),
             StationColumn(name='coordinates'),
             StationColumn(name='instrument')])

filter_editor = TableEditor(
    sortable=True,
    configurable=True,
    auto_size=False,
    row_factory=TFilter,
    auto_add=False,
    show_toolbar=True,
    columns=[StationColumn(name='ref', editable=False),
             CheckboxColumn(name='used', label='Used ?'),
             StationColumn(name='low'),
             StationColumn(name='mwcs_low'),
             StationColumn(name='mwcs_high'),
             StationColumn(name='high'),
             StationColumn(name='rms_threshold'),
             StationColumn(name='mwcs_wlen'),
             StationColumn(name='mwcs_step'),
             ])

config_editor = TableEditor(
    sortable=True,
    configurable=True,
    auto_size=True,
    columns=[
        ObjectColumn(name='name', editable=False),
        ObjectColumn(name='value',),
        ObjectColumn(name='info', editable=False),
        ObjectColumn(name='default', editable=False),
    ])


class TNetwork(HasTraits):
    """ Defines a Network with stations."""

    name = Str('<unknown>')
    stations = List(TStation)


class ConfigItem(HasTraits):
    name = Str('<unknown>')
    value = Str('')
    info = Str('')
    default = Str('')


class Software(HasTraits):

    """ Defines an MSNoise configuration instance. """

    name = Str('<unknown>')
    networks = List(TNetwork)
    configs = List(ConfigItem)
    filters = List(TFilter)

no_view = View()

tree_editor = TreeEditor(
    nodes=[
        TreeNode(node_for=[Software],
                 auto_open=True,
                 children='',
                 label='name',
                 view=View(Item('configs',
                                show_label=False, editor=config_editor))),
        TreeNode(node_for=[Software],
                 auto_open=True,
                 children='networks',
                 label='=Networks',
                 view=no_view,
                 add=[TNetwork],),
        TreeNode(node_for=[TNetwork],
                 auto_open=False,
                 children='stations',
                 label='name',
                 view=View(Item('stations',
                                show_label=False,
                                editor=station_editor)
                           ),
                 add=[TStation]),
        TreeNode(node_for=[TStation],
                 auto_open=True,
                 label='sta',
                 view=View(
                     ['net', 'sta', 'X', 'Y', 'altitude', 'coordinates',
                      'used', 'instrument'])),

        TreeNode(node_for=[Software],
                 auto_open=True,
                 label='=Filters',
                 view=View(Item('filters',
                                show_label=False,
                                editor=filter_editor)),
                 add=[TFilter],),
    ]
)


class StationConfigurator(HasTraits):

    """ Defines a business partner."""

    name = Str('<unknown>')
    company = Instance(Software)

    view = View(
        Item(name='company',
             editor=tree_editor,
             show_label=False,
             ),
        title='Station Configuration GUI',
        buttons=['OK', 'Cancel'],
        resizable=True,
        style='custom',
        width=1024,
        height=768
    )


def main():
    db = connect()

    networks = []
    for network in get_networks(db, all=True):
        networks.append(TNetwork(name=network))
        for s in get_stations(db, net=network, all=True):
            station = TStation(ref=s.ref, net=s.net, sta=s.sta, X=s.X, Y=s.Y,
                               altitude=s.altitude, used=s.used,
                               instrument=s.instrument,
                               coordinates=s.coordinates)
            networks[-1].stations.append(station)

    filters = []
    for f in get_filters(db, all=True):
        filters.append(
            TFilter(
                ref=f.ref, low=f.low, high=f.high, mwcs_low=f.mwcs_low,
                mwcs_high=f.mwcs_high, rms_threshold=f.rms_threshold,
                mwcs_wlen=f.mwcs_wlen, mwcs_step=f.mwcs_step, used=f.used))

    configs = []
    for name in default.keys():
        value = get_config(db, name)
        configs.append(
            ConfigItem(name=name, value=value, info=default[name][0],
                       default=default[name][1]))

    db.close()
    demo = StationConfigurator(
        name='MSNoise',
        company=Software(
            name='MSNoise',
            networks=networks,
            configs=configs,
            filters=filters
        )
    )

    if demo.configure_traits():
        db = connect()

        print("Updating Config Table")
        for config in demo.company.configs:
            update_config(db, config.name, config.value)

        print("Updating Station Table")
        for network in demo.company.networks:
            for s in network.stations:
                update_station(
                    db, s.net, s.sta, s.X, s.Y, s.altitude, s.coordinates,
                    s.instrument, s.used)

        print("Updating Filter Table")
        for f in demo.company.filters:
            update_filter(db, f.ref, f.low, f.mwcs_low, f.high,
                          f.mwcs_high, f.rms_threshold, f.mwcs_wlen,
                          f.mwcs_step, f.used)

        db.close()
        print("Done !")

if __name__ == '__main__':
    main()

# EOF
