{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline }}

.. current module:: {{ module }}

.. auto{{ objtype }}:: {{objname }}
