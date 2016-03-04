.. include:: ../configs.hrst

Extending MSNoise with Plugins
==============================

.. versionadded:: 1.4

Starting with :doc:`releasenotes/msnoise-1.4`, MSNoise supports Plugins, this
means the default workflow "from archive to dv/v" can be branched at any step!


Plugin Declaration
------------------

After installing a plugin, its package name must be declared in the ``plugins``
parameter in the configuration.