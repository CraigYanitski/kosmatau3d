{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline}}

:param {{ args[0] }}: {{ args[0] | escape }}{% if default_values[0] %} (default: {{ default_values[0] }}){% endif %}
:type {{ args[0] }}: {{ types[0] | escape }}
:returns: {{ rtype | escape }}{% if return_default %} (default: {{ return_default }}){% endif %}
:rtype: {{ return_type | escape }}
{% if options.docstring %}
{{ options.docstring | rst }}
{% endif %}