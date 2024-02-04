{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline }}

.. autofunction:: {{ fullname }}

   {% block parameters %}
   {% if args %}
   :param {{ args[0] }}: {{ args[0] | escape }}{% if default_values[0] %} (default: {{ default_values[0] }}){% endif %}
   :type {{ args[0] }}: {{ types[0] | escape }}
   {% endif %}
   {% if args[1:] %}
   {% for arg in args[1:] %}
   :param {{ arg }}: {{ arg | escape }}
   {% endfor %}
   {% endif %}
   {% endblock %}

   {% if returns %}
   :returns: {{ returns | escape }}
   {% endif %}
   {% if rtype %}
   :rtype: {{ rtype | escape }}
   {% endif %}

   {% if options.docstring %}
   {{ options.docstring | rst }}
   {% endif %}
