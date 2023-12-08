{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline}}

.. currentmodule:: {{ module }}

.. autoclass:: {{ objname }}
   :members:
   :show-inheritance:
   :inherited-members:
   :special-members: __call__, __add__, __mul__

   {% block methods %}
   {% if methods %}
   .. rubric:: {{ _('Methods') }}

   .. autosummary::
      :nosignatures:
   {% for item in methods %}
      {%- if not item.startswith('_') %}
      {%- set shortname = fullname.split('.')[-1] %}
      ~{{ shortname }}.{{ item }}
      {%- endif -%}
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: {{ _('Attributes') }}

   .. autosummary::
   {% for item in attributes %}
      {%- set shortname = fullname.split('.')[-1] %}
      ~{{ shortname }}.{{ item }}
   {%- endfor %}
   {% endif %}
   {% endblock %}
