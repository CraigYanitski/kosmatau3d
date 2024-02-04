{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline }}

.. currentmodule:: {{ module }}

.. autofunction:: {{objname }}
