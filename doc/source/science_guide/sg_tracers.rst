:tocdepth: 3

.. _tracers:

Tracers
=======

The basic conservation equations for ice area fraction :math:`a_{in}`,
ice volume :math:`v_{in}`, and snow volume :math:`v_{sn}` for each
thickness category :math:`n` are

.. math::
   {\frac{\partial}{\partial t}} (a_{in}) + \nabla \cdot (a_{in} {\bf u}) = 0,
   :label: transport-ai

.. math::
   \frac{\partial v_{in}}{\partial t} + \nabla \cdot (v_{in} {\bf u}) = 0,
   :label: transport-vi

.. math::
   \frac{\partial v_{sn}}{\partial t} + \nabla \cdot (v_{sn} {\bf u}) = 0.
   :label: transport-vs

The ice and snow volumes can be written equivalently in terms of
tracers, ice thickness :math:`h_{in}` and snow depth :math:`h_{sn}`:

.. math::
   \frac{\partial h_{in}a_{in}}{\partial t} + \nabla \cdot (h_{in}a_{in} {\bf u}) = 0,
   :label: transport-hi

.. math::
   \frac{\partial h_{sn}a_{in}}{\partial t} + \nabla \cdot (h_{sn}a_{in} {\bf u}) = 0.
   :label: transport-hs

Although we maintain ice and snow volume instead of the thicknesses as
state variables in CICE, the tracer form is used for volume transport
(section :ref:`horiz-trans`). There are many other tracers
available, whose values are contained in the `trcrn` array. Their
transport equations typically have one of the following three forms

.. math::
   \frac{\partial \left(a_{in} T_n\right)}{\partial t} + \nabla \cdot (a_{in} T_n {\bf u}) = 0,
   :label: transport-aT

.. math::
   \frac{\partial \left(v_{in} T_n\right)}{\partial t} + \nabla \cdot (v_{in} T_n {\bf u}) = 0,
   :label: transport-viT

.. math::
   \frac{\partial \left(v_{sn} T_n\right)}{\partial t} + \nabla \cdot (v_{sn} T_n {\bf u}) = 0.
   :label: transport-vsT

Equation :eq:`transport-aT` describes the transport of surface
temperature, whereas Equation :eq:`transport-viT` and Equation :eq:`transport-vsT`
describe the transport of ice and snow enthalpy, salt, and passive
tracers such as volume-weighted ice age and snow age. Each tracer field
is given an integer index, `trcr_depend`, which has the value 0, 1, or 2
depending on whether the appropriate conservation equation is
Equation :eq:`transport-aT`, Equation :eq:`transport-viT`, or Equation :eq:`transport-vsT`,
respectively. The total number of tracers is
:math:`N_{tr}\ge 1`. In the default configuration there are two
tracers: surface temperature and volume-weighted ice age. Tracers for
melt ponds, level ice area and level ice volume (used to diagnose ridged ice area and volume) 
are also available. Users may add any number of additional tracers that are transported 
conservatively provided that `trcr_depend` is defined appropriately. See Section :ref:`addtrcr` 
for guidance on adding tracers.

Please see the `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_ for additional information about tracers that depend on other tracers, age of the ice, aerosols, 
brine height, and the sea ice ecosystem.
