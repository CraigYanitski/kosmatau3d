{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline }}

.. currentmodule:: kosmatau3d

.. autofunction:: {{ shortname }}
