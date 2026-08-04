Optional Lineage API
====================

Bioscrape's lineage simulator is distributed from the source tree as an
optional extension. Build it from the repository with::

   python setup.py install lineage

The user-facing concepts are documented in 
:doc:`../user_guide/lineage_package`. The source for the extension lives
in ``lineage/lineage.pyx`` and is exposed as ``bioscrape.lineage``
when the optional extension is built.

.. automodule:: bioscrape.lineage
   :members:
   :show-inheritance: