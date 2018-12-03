:tocdepth: 3

.. _modelcomps:

Fundamental Variables
=====================

The Arctic and Antarctic sea ice packs are mixtures of open water, thin
first-year ice, thicker multiyear ice, and thick pressure ridges. The
thermodynamic and dynamic properties of the ice pack depend on how much
ice lies in each thickness range. Thus the basic problem in sea ice
modeling is to describe the evolution of the ice thickness distribution
(ITD) in time and space.

The fundamental equation solved by CICE is :cite:`Thorndike75`:

.. math::
   \frac{\partial g}{\partial t} = -\nabla \cdot (g {\bf u}) 
    - \frac{\partial}{\partial h} (f g) + \psi,
   :label: transport-g

where :math:`{\bf u}` is the horizontal ice velocity,
:math:`\nabla = (\frac{\partial}{\partial x}, \frac{\partial}{\partial y})`,
:math:`f` is the rate of thermodynamic ice growth, :math:`\psi` is a
ridging redistribution function, and :math:`g` is the ice thickness
distribution function. We define :math:`g({\bf x},h,t)\,dh` as the
fractional area covered by ice in the thickness range :math:`(h,h+dh)`
at a given time and location.

Additional information about the ITD for CICE can be found in the
`Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_.

In addition to the fractional ice area, :math:`a_{in}`, we define the
following state variables for each category :math:`n`. In a change from
previous CICE versions, we no longer carry snow and ice energy as
separate variables; instead they and sea ice salinity are carried as
tracers on snow and ice volume.

-  :math:`v_{in}`, the ice volume, equal to the product of
   :math:`a_{in}` and the ice thickness :math:`h_{in}`.

-  :math:`v_{sn}`, the snow volume, equal to the product of
   :math:`a_{in}` and the snow thickness :math:`h_{sn}`.

-  :math:`e_{ink}`, the internal ice energy in layer :math:`k`, equal to
   the product of the ice layer volume, :math:`v_{in}/N_i`, and the ice
   layer enthalpy, :math:`q_{ink}`. Here :math:`N_i` is the total number
   of ice layers, with a default value :math:`N_i = 4`, and
   :math:`q_{ink}` is the negative of the energy needed to melt a unit
   volume of ice and raise its temperature to :math:`0\ ^{\circ}`\ C. 
   (NOTE: In the current code, :math:`e_i<0`
   and :math:`q_i<0` with :math:`e_i = v_iq_i`.)

-  :math:`e_{snk}`, the internal snow energy in layer :math:`k`, equal
   to the product of the snow layer volume, :math:`v_{sn}/N_s`, and the
   snow layer enthalpy, :math:`q_{snk}`, where :math:`N_s` is the number
   of snow layers. (Similarly, :math:`e_s<0` in the code.) CICE allows
   multiple snow layers, but the default value is :math:`N_s=1`.

-  :math:`S_i`, the bulk sea ice salt content in layer :math:`k`, equal
   to the product of the ice layer volume and the sea ice salinity
   tracer.

-  :math:`T_{sfn}`, the surface temperature.

Since the fractional area is unitless, the volume variables have units
of meters (i.e., m\ :math:`^3` of ice or snow per m\ :math:`^2` of grid
cell area), and the energy variables have units of J/m\ :math:`^2`.

The three terms on the right-hand side of Equation :eq:`transport-g` describe
three kinds of sea ice transport: (1) horizontal transport in
:math:`(x,y)` space; (2) transport in thickness space :math:`h` due to
thermodynamic growth and melting; and (3) transport in thickness space
:math:`h` due to ridging and other mechanical processes. We solve the
equation by operator splitting in three stages, with two of the three
terms on the right set to zero in each stage. We compute horizontal
transport using the incremental remapping scheme of
:cite:`Dukowicz00` as adapted for sea ice by
:cite:`Lipscomb04`; this scheme is discussed in
Section :ref:`horiz-trans`. Ice is transported in thickness space
using the remapping scheme of :cite:`Lipscomb01`. The mechanical
redistribution scheme, based on :cite:`Thorndike75`, :cite:`Rothrock75`,
:cite:`Hibler80`, :cite:`Flato95`, and :cite:`Lipscomb07` is outlined
in the `Icepack Documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_. 
To solve the horizontal transport and
ridging equations, we need the ice velocity :math:`{\bf u}`, and to
compute transport in thickness space, we must know the the ice growth
rate :math:`f` in each thickness category. We use the
elastic-viscous-plastic (EVP) ice dynamics scheme of
:cite:`Hunke97`, as modified by :cite:`Connolley04`,
:cite:`Hunke01`, :cite:`Hunke02` and
:cite:`Hunke03`, or a new elastic-anisotropic-plastic model
:cite:`Wilchinsky06,Weiss09,Tsamados13` to find the velocity, as
described in Section :ref:`dynam`. Finally, we use a thermodynamic
model to compute :math:`f`. The order in which
these computations are performed in the code itself was chosen so that
quantities sent to the coupler are consistent with each other and as
up-to-date as possible. The Delta-Eddington radiative scheme computes
albedo and shortwave components simultaneously, and in order to have the
most up-to-date values available for the coupler at the end of the
timestep, the order of radiation calculations is shifted. Albedo and
shortwave components are computed after the ice state has been modified
by both thermodynamics and dynamics, so that they are consistent with
the ice area and thickness at the end of the step when sent to the
coupler. However, they are computed using the downwelling shortwave from
the beginning of the timestep. Rather than recompute the albedo and
shortwave components at the beginning of the next timestep using new
values of the downwelling shortwave forcing, the shortwave components
computed at the end of the last timestep are scaled for the new forcing.
