{%- set shortname = '.'.join(fullname.split('.')[1:]) %}
.. _{{ shortname }}:

{{ shortname | escape | underline }}

.. automodule:: {{ fullname }}

   {% block attributes %}
   {% if attributes %}
   .. rubric:: Module attributes

   .. autosummary::
      :toctree:
   {% for item in attributes %}
      {%- set shortitem = '.'.join(item.split('.')[1:]) %}
      {{ item }} | shortitem
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block functions %}
   {% if functions %}
   .. rubric:: {{ _('Functions') }}

   .. autosummary::
      :toctree:
      :template: custom-function-template.rst
      :nosignatures:
   {% for item in functions %}
      {%- set shortitem = '.'.join(item.split('.')[1:]) %}
      {{ item }} | shortitem
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block classes %}
   {% if classes %}
   .. rubric:: {{ _('Classes') }}

   .. autosummary::
      :toctree:
      :template: custom-class-template.rst
      :nosignatures:
   {% for item in classes %}
      {%- set shortitem = '.'.join(item.split('.')[1:]) %}
      {{ item }} | shortitem
   {%- endfor %}
   {% endif %}
   {% endblock %}

   {% block exceptions %}
   {% if exceptions %}
   .. rubric:: {{ _('Exceptions') }}

   .. autosummary::
      :toctree:
   {% for item in exceptions %}
      {%- set shortitem = '.'.join(item.split('.')[1:]) %}
      {{ item }} | shortitem
   {%- endfor %}
   {% endif %}
   {% endblock %}

{% block modules %}
{% if modules %}
.. rubric:: {{ _('Modules') }}

.. autosummary::
   :toctree:
   :template: custom-module-template.rst
   :recursive:
{% for item in modules %}
   {%- set shortitem = '.'.join(item.split('.')[1:]) %}
   :ref:`{{ shortitem }} <{{ item }}>`
{%- endfor %}
{% endif %}
{% endblock %}
