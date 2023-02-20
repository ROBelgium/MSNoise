Interaction Examples & Gallery
==============================

The following examples are meant to show you how to interact with MSNoise using
its API, thus avoiding having to dive into the folder structure.

Users should try examples while checking the :doc:`../api`. (application
programming interface) for understanding the calls to different functions.

In a nutshell, all examples start with the following Python code:

.. code:: python

    from msnoise.api import db
    db = connect()

This, if run in an MSNoise project folder (= a folder where you have already
run ``msnoise db init`` (as explained in :ref:`workflow`)), will provide a
``Session`` object, connected to the database. 