:tocdepth: 3

.. _modelcomps:

Sea Ice Model Components
=========================

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
`Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/latest/science_guide/sg_itd.html#ice-thickness-distribution>`_.

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
   volume of ice and raise its temperature to :math:`0\ ^{\circ}`\ C; it is discussed in
   the `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/sg_thermo.html#thermodynamics>`_. (NOTE: In the current code, :math:`e_i<0`
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
using the remapping scheme of :cite:`Lipscomb01`, as
described in the `Icepack ITD documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/sg_itd.html#ice-thickness-distribution>`). The mechanical
redistribution scheme, based on :cite:`Thorndike75`, :cite:`Rothrock75`,
:cite:`Hibler80`, :cite:`Flato95`, and :cite:`Lipscomb07` is outlined
in Section :ref:`mech-red`. To solve the horizontal transport and
ridging equations, we need the ice velocity :math:`{\bf u}`, and to
compute transport in thickness space, we must know the the ice growth
rate :math:`f` in each thickness category. We use the
elastic-viscous-plastic (EVP) ice dynamics scheme of
:cite:`Hunke97`, as modified by :cite:`Connolley04`,
:cite:`Hunke01`, :cite:`Hunke02` and
:cite:`Hunke03`, or a new elastic-anisotropic-plastic model
:cite:`Wilchinsky06,Weiss09,Tsamados13` to find the velocity, as
described in Section :ref:`dynam`. Finally, we use a thermodynamic
model to compute :math:`f` (See the 
`Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/sg_thermo.html#thermodynamics>`_). The order in which
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

.. _tracers:

~~~~~~~
Tracers
~~~~~~~

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
melt ponds (See the Icepack documentation for
`Tracers <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/sg_tracers.html#tracers-that-depend-on-other-tracers>`_
and
'Melt Ponds <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/sg_thermo.html#melt-ponds>`_),
level ice
area and level ice volume (used to diagnose ridged ice area and volume,
Section :ref:`mech-red`) are also available. Users may add any number
of additional tracers that are transported conservatively provided that
`trcr_depend` is defined appropriately. See Section :ref:`addtrcr` for
guidance on adding tracers.


Additional information about tracers in CICE can be found in the Icepack documentation. 
Links the the Icepack online documentation are provided below and the Icepack
documentation can be accessed through the GitHub wiki

-  `Tracers that depend on other tracers <https://cice-consortium-icepack.readthedocs.io/en/latest/science_guide/sg_tracers.html#tracers-that-depend-on-other-tracers>`_

- `Age of the ice <https://cice-consortium-icepack.readthedocs.io/en/latest/science_guide/sg_tracers.html#ice-age>`_

- `Aerosols <https://cice-consortium-icepack.readthedocs.io/en/latest/science_guide/sg_bgc.html#aerosols>`_

- `Brine height <https://cice-consortium-icepack.readthedocs.io/en/latest/science_guide/sg_bgc.html#brine-height>`_

- `Sea ice ecosystem <https://cice-consortium-icepack.readthedocs.io/en/latest/science_guide/sg_bgc.html#sea-ice-ecosystem>`_


.. _horiz-trans:

~~~~~~~~~~~~~~~~~~~~
Horizontal transport
~~~~~~~~~~~~~~~~~~~~

We wish to solve the continuity or transport equation
(Equation :eq:`transport-ai`) for the fractional ice area in each
thickness category :math:`n`. Equation :eq:`transport-ai` describes
the conservation of ice area under horizontal transport. It is obtained
from Equation :eq:`transport-g` by discretizing :math:`g` and neglecting the
second and third terms on the right-hand side, which are treated
separately (Sections :ref:`mech-red` and in the
`Icepack ITD documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/sg_itd.html#ice-thickness-distribution>`).

There are similar conservation equations for ice volume
(Equation :eq:`transport-vi`), snow volume (Equation :eq:`transport-vs`), ice
energy and snow energy:

.. math::
   \frac{\partial e_{ink}}{\partial t} + \nabla \cdot (e_{ink} {\bf u}) = 0,
   :label: transport-ei

.. math::
   \frac{\partial e_{snk}}{\partial t} + \nabla \cdot (e_{snk} {\bf u}) = 0.
   :label: transport-es

By default, ice and snow are assumed to have constant densities, so that
volume conservation is equivalent to mass conservation. Variable-density
ice and snow layers can be transported conservatively by defining
tracers corresponding to ice and snow density, as explained in the
introductory comments in **ice\_transport\_remap.F90**. Prognostic
equations for ice and/or snow density may be included in future model
versions but have not yet been implemented.

Two transport schemes are available: upwind and the incremental
remapping scheme of :cite:`Dukowicz00` as modified for sea ice by
:cite:`Lipscomb04`. The remapping scheme has several desirable features:

-  It conserves the quantity being transported (area, volume, or
   energy).

-  It is non-oscillatory; that is, it does not create spurious ripples
   in the transported fields.

-  It preserves tracer monotonicity. That is, it does not create new
   extrema in the thickness and enthalpy fields; the values at
   time \ :math:`m+1` are bounded by the values at time \ :math:`m`.

-  It is second-order accurate in space and therefore is much less
   diffusive than first-order schemes (e.g., upwind). The accuracy may
   be reduced locally to first order to preserve monotonicity.

-  It is efficient for large numbers of categories or tracers. Much of
   the work is geometrical and is performed only once per grid cell
   instead of being repeated for each quantity being transported.

The time step is limited by the requirement that trajectories projected
backward from grid cell corners are confined to the four surrounding
cells; this is what is meant by incremental remapping as opposed to
general remapping. This requirement leads to a CFL-like condition,

.. math::
   {\max|{\bf u}|\Delta t\over\Delta x}
   \leq 1.

For highly divergent velocity fields the maximum time step must be
reduced by a factor of two to ensure that trajectories do not cross.
However, ice velocity fields in climate models usually have small
divergences per time step relative to the grid size.

The remapping algorithm can be summarized as follows:

#. Given mean values of the ice area and tracer fields in each grid
   cell, construct linear approximations of these fields. Limit the
   field gradients to preserve monotonicity.

#. Given ice velocities at grid cell corners, identify departure regions
   for the fluxes across each cell edge. Divide these departure regions
   into triangles and compute the coordinates of the triangle vertices.

#. Integrate the area and tracer fields over the departure triangles to
   obtain the area, volume, and energy transported across each cell
   edge.

#. Given these transports, update the state variables.

Since all scalar fields are transported by the same velocity field, step
(2) is done only once per time step. The other three steps are repeated
for each field in each thickness category. These steps are described
below.

.. _reconstruct:

*************************************
Reconstructing area and tracer fields
*************************************

First, using the known values of the state variables, the ice area and
tracer fields are reconstructed in each grid cell as linear functions of
:math:`x` and :math:`y`. For each field we compute the value at the cell
center (i.e., at the origin of a 2D Cartesian coordinate system defined
for that grid cell), along with gradients in the :math:`x` and :math:`y`
directions. The gradients are limited to preserve monotonicity. When
integrated over a grid cell, the reconstructed fields must have mean
values equal to the known state variables, denoted by :math:`\bar{a}`
for fractional area, :math:`\tilde{h}` for thickness, and
:math:`\hat{q}` for enthalpy. The mean values are not, in general, equal
to the values at the cell center. For example, the mean ice area must
equal the value at the centroid, which may not lie at the cell center.

Consider first the fractional ice area, the analog to fluid density
:math:`\rho` in :cite:`Dukowicz00`. For each thickness category
we construct a field :math:`a({\bf r})` whose mean is :math:`\bar{a}`,
where :math:`{\bf r} =
(x,y)` is the position vector relative to the cell center. That is, we
require

.. math::
   \int_A a \, dA = {\bar a} \, A,
   :label: mean-area

where :math:`A=\int_A dA` is the grid cell area.
Equation :eq:`mean-area` is satisfied if :math:`a({\bf r})` has the
form

.. math::
   a({\bf r}) = \bar{a} + \alpha_a \left<\nabla a\right> \cdot ({\bf r}-{\bf \bar{r}}),
   :label: recon-area

where :math:`\left<\nabla a\right>` is a centered estimate of the area
gradient within the cell, :math:`\alpha_a` is a limiting coefficient
that enforces monotonicity, and :math:`{\bf \bar{r}}` is the cell
centroid:

.. math:: 
   {\bf \bar{r}} = {1\over A} \int_A {\bf r} \, dA.

It follows from Equation :eq:`recon-area` that the ice area at the cell center
(:math:`\mathbf{r} = 0`) is

.. math:: 
   a_c = \bar{a} - a_x \overline{x} - a_y \overline{y},

where :math:`a_x = \alpha_a (\partial a / \partial x)` and
:math:`a_y = \alpha_a (\partial a / \partial y)` are the limited
gradients in the :math:`x` and :math:`y` directions, respectively, and
the components of :math:`{\bf \bar{r}}`,
:math:`\overline{x} = \int_A x \, dA / A` and
:math:`\overline{y} = \int_A y \, dA / A`, are evaluated using the
triangle integration formulas described in
Section :ref:`integ-flux`. These means, along with higher-order
means such as :math:`\overline{x^2}`, :math:`\overline{xy}`, and
:math:`\overline{y^2}`, are computed once and stored.

Next consider the ice and snow thickness and enthalpy fields. Thickness
is analogous to the tracer concentration :math:`T` in
:cite:`Dukowicz00`, but there is no analog in
:cite:`Dukowicz00` to the enthalpy. The reconstructed ice or snow
thickness :math:`h({\bf r})` and enthalpy :math:`q(\mathbf{r})` must
satisfy

.. math::
   \int_A a \, h \, dA       =  \bar{a} \, \tilde{h} \, A,
   :label: mean-thickness

.. math::
   \int_A a \, h \, q \, dA  =  \bar{a} \, \tilde{h} \, \hat{q} \, A,
   :label: mean-enthalpy

where :math:`\tilde{h}=h(\tilde{\bf r})` is the thickness at the center
of ice area, and :math:`\hat{q}=q(\hat{\bf r})` is the enthalpy at the
center of ice or snow volume. Equations :eq:`mean-thickness` and
:eq:`mean-enthalpy` are satisfied when :math:`h({\bf r})` and
:math:`q({\bf r})` are given by

.. math::
   h({\bf r}) = \tilde{h} + \alpha_h \left<\nabla h\right> \cdot
                                        ({\bf r}-{\bf \tilde{r}}),
   :label: recon-thickness

.. math::
   q({\bf r}) = \hat{q} + \alpha_q \left<\nabla q\right> \cdot
                                      ({\bf r}-{\bf \hat{r}}),
   :label: recon-enthalpy

where :math:`\alpha_h` and :math:`\alpha_q` are limiting coefficients.
The center of ice area, :math:`{\bf\tilde{r}}`, and the center of ice or
snow volume, :math:`{\bf \hat{r}}`, are given by

.. math:: 
   {\bf \tilde{r}} = {1\over\bar{a} \, A}\int_A a \, {\bf r} \, dA,

.. math::
   {\bf \hat{r}} =
           {1\over\bar{a} \, \tilde{h} \, A}\int_A a \, h \, {\bf r} \, dA.

Evaluating the integrals, we find that the components of
:math:`{\bf \tilde{r}}` are

.. math::
   \tilde{x} = \frac{a_c \overline{x}+a_x \overline{x^2}+a_y \overline{xy}}
                      {\bar{a}},

.. math::
   \tilde{y} = \frac{a_c \overline{y}+a_x \overline{xy} +a_y \overline{y^2}}
                      {\bar{a}},

and the components of :math:`{\bf \hat{r}}` are

.. math::
   \hat{x} = \frac { c_1 \overline{x}     + c_2 \overline{x^2}
                      + c_3 \overline{xy}    + c_4 \overline{x^3}
                      + c_5 \overline{x^2 y} + c_6 \overline{x y^2} }
                      {\bar{a} \, \tilde{h}},

.. math::
   \hat{y} = \frac { c_1 \overline{y}     + c_2 \overline{xy}
                      + c_3 \overline{y^2}   + c_4 \overline{x^2 y}
                      + c_5 \overline{x y^2} + c_6 \overline{y^3}   }
                       {\bar{a} \, \tilde{h}},

where

.. math::
   \begin{aligned}
    c_1 & \equiv & a_c h_c,            \\
    c_2 & \equiv & a_c h_x + a_x h_c,  \\
    c_3 & \equiv & a_c h_y + a_y h_c,  \\
    c_4 & \equiv & a_x h_x,            \\
    c_5 & \equiv & a_x h_y + a_y h_x,  \\
    c_6 & \equiv & a_y h_y.\end{aligned}

From Equation :eq:`recon-thickness` and Equation :eq:`recon-enthalpy`, 
the thickness and enthalpy at the cell center are given by

.. math:: 
   h_c = \tilde{h} - h_x \tilde{x} - h_y \tilde{y},

.. math:: 
   q_c = \hat{q}   - q_x \hat{x}   - q_y \hat{y},

where :math:`h_x`, :math:`h_y`, :math:`q_x` and :math:`q_y` are the
limited gradients of thickness and enthalpy. The surface temperature is
treated the same way as ice or snow thickness, but it has no associated
enthalpy. Tracers obeying conservation equations of the form Equation
:eq:`transport-viT` and Equation :eq:`transport-vsT` are treated in analogy
to ice and snow enthalpy, respectively.

We preserve monotonicity by van Leer limiting. If
:math:`\bar{\phi}(i,j)` denotes the mean value of some field in grid
cell :math:`(i,j)`, we first compute centered gradients of
:math:`\bar{\phi}` in the :math:`x` and :math:`y` directions, then check
whether these gradients give values of :math:`\phi` within cell
:math:`(i,j)` that lie outside the range of :math:`\bar{\phi}` in the
cell and its eight neighbors. Let :math:`\bar{\phi}_{\max}` and
:math:`\bar{\phi}_{\min}` be the maximum and minimum values of
:math:`\bar{\phi}` over the cell and its neighbors, and let
:math:`\phi_{\max}` and :math:`\phi_{\min}` be the maximum and minimum
values of the reconstructed :math:`\phi` within the cell. Since the
reconstruction is linear, :math:`\phi_{\max}` and :math:`\phi_{\min}`
are located at cell corners. If :math:`\phi_{\max} > \bar{\phi}_{\max}`
or :math:`\phi_{\min} < \bar{\phi}_{\min}`, we multiply the unlimited
gradient by :math:`\alpha = \min(\alpha_{\max}, \alpha_{\min})`, where

.. math::
   \alpha_{\max} =
     (\bar{\phi}_{\max} - \bar{\phi}) / (\phi_{\max} -\bar{\phi}),

.. math::
   \alpha_{\min} =
     (\bar{\phi}_{\min} - \bar{\phi}) / (\phi_{\min} -\bar{\phi}).

Otherwise the gradient need not be limited.

Earlier versions of CICE (through 3.14) computed gradients in physical
space. Starting in version 4.0, gradients are computed in a scaled space
in which each grid cell has sides of unit length. The origin is at the
cell center, and the four vertices are located at (0.5, 0.5),
(-0.5,0.5),(-0.5, -0.5) and (0.5, -0.5). In this coordinate system,
several of the above grid-cell-mean quantities vanish (because they are
odd functions of x and/or y), but they have been retained in the code
for generality.

.. _loc-dep-triangles:

****************************
Locating departure triangles
****************************

The method for locating departure triangles is discussed in detail by
:cite:`Dukowicz00`. The basic idea is illustrated in
:ref:`fig-deparr`, which shows a shaded quadrilateral departure region
whose contents are transported to the target or home grid cell, labeled
:math:`H`. The neighboring grid cells are labeled by compass directions:
:math:`NW`, :math:`N`, :math:`NE`, :math:`W`, and :math:`E`. The four
vectors point along the velocity field at the cell corners, and the
departure region is formed by joining the starting points of these
vectors. Instead of integrating over the entire departure region, it is
convenient to compute fluxes across cell edges. We identify departure
regions for the north and east edges of each cell, which are also the
south and west edges of neighboring cells. Consider the north edge of
the home cell, across which there are fluxes from the neighboring
:math:`NW` and :math:`N` cells. The contributing region from the
:math:`NW` cell is a triangle with vertices :math:`abc`, and that from
the :math:`N` cell is a quadrilateral that can be divided into two
triangles with vertices :math:`acd` and :math:`ade`. Focusing on
triangle :math:`abc`, we first determine the coordinates of vertices
:math:`b` and :math:`c` relative to the cell corner (vertex :math:`a`),
using Euclidean geometry to find vertex :math:`c`. Then we translate the
three vertices to a coordinate system centered in the :math:`NW` cell.
This translation is needed in order to integrate fields
(Section :ref:`integ-flux`) in the coordinate system where they
have been reconstructed (Section :ref:`reconstruct`). Repeating
this process for the north and east edges of each grid cell, we compute
the vertices of all the departure triangles associated with each cell
edge.

.. _fig-deparr:

.. figure:: ./figures/deparr.png
   :align: center
   :scale: 20%
 
   Departure Region

Figure :ref:`fig-deparr` shows that in incremental remapping, conserved quantities are
remapped from the shaded departure region, a quadrilateral formed by
connecting the backward trajectories from the four cell corners, to
the grid cell labeled :math:`H`. The region fluxed across the north
edge of cell :math:`H` consists of a triangle (:math:`abc`) in the
:math:`NW` cell and a quadrilateral (two triangles, :math:`acd` and
:math:`ade`) in the :math:`N` cell.


Figure :ref:`fig-triangles`, reproduced from :cite:`Dukowicz00`, shows
all possible triangles that can contribute fluxes across the north edge
of a grid cell. There are 20 triangles, which can be organized into five
groups of four mutually exclusive triangles as shown in
:ref:`tab-triangle`. In this table, :math:`(x_1, y_1)` and
:math:`(x_2,y_2)` are the Cartesian coordinates of the departure points
relative to the northwest and northeast cell corners, respectively. The
departure points are joined by a straight line that intersects the west
edge at :math:`(0,y_a)` relative to the northwest corner and intersects
the east edge at :math:`(0,y_b)` relative to the northeast corner. The
east cell triangles and selecting conditions are identical except for a
rotation through 90 degrees.

.. _fig-triangles:

.. figure:: ./figures/triangles.png
   :align: center
   :scale: 20%

   Triangles

Table :ref:`tab-triangle` show the evaluation of contributions from the 20
triangles across the north cell edge. The coordinates :math:`x_1`,
:math:`x_2`, :math:`y_1`, :math:`y_2`, :math:`y_a`, and :math:`y_b` are
defined in the text. We define :math:`\tilde{y}_1 =
y_1` if :math:`x_1>0`, else :math:`\tilde{y}_1 = y_a`. Similarly,
:math:`\tilde{y}_2
= y_2` if :math:`x_2<0`, else :math:`\tilde{y}_2 = y_b`.

.. _tab-triangle:

.. table:: Triangular Contributions

   +------------+------------+--------------------------------------------------------+----+
   | Triangle   | Triangle   | Selecting logical condition                            |    |
   | group      | label      |                                                        |    |
   +------------+------------+--------------------------------------------------------+----+
   | 1          | NW         | :math:`y_a>0` and :math:`y_1\geq0` and :math:`x_1<0`   |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | NW1        | :math:`y_a<0` and :math:`y_1\geq0` and :math:`x_1<0`   |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | W          | :math:`y_a<0` and :math:`y_1<0` and :math:`x_1<0`      |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | W2         | :math:`y_a>0` and :math:`y_1<0` and :math:`x_1<0`      |    |
   +------------+------------+--------------------------------------------------------+----+
   +------------+------------+--------------------------------------------------------+----+
   | 2          | NE         | :math:`y_b>0` and :math:`y_2\geq0` and :math:`x_2>0`   |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | NE1        | :math:`y_b<0` and :math:`y_2\geq0` and :math:`x_2>0`   |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | E          | :math:`y_b<0` and :math:`y_2<0` and :math:`x_2>0`      |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | E2         | :math:`y_b>0` and :math:`y_2<0` and :math:`x_2>0`      |    |
   +------------+------------+--------------------------------------------------------+----+
   +------------+------------+--------------------------------------------------------+----+
   | 3          | W1         | :math:`y_a<0` and :math:`y_1\geq0` and :math:`x_1<0`   |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | NW2        | :math:`y_a>0` and :math:`y_1<0` and :math:`x_1<0`      |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | E1         | :math:`y_b<0` and :math:`y_2\geq0` and :math:`x_2>0`   |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | NE2        | :math:`y_b>0` and :math:`y_2<0` and :math:`x_2>0`      |    |
   +------------+------------+--------------------------------------------------------+----+
   +------------+------------+--------------------------------------------------------+----+
   | 4          | H1a        | :math:`y_a y_b\geq 0` and :math:`y_a+y_b<0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | N1a        | :math:`y_a y_b\geq 0` and :math:`y_a+y_b>0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | H1b        | :math:`y_a y_b<0` and :math:`\tilde{y}_1<0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | N1b        | :math:`y_a y_b<0` and :math:`\tilde{y}_1>0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   +------------+------------+--------------------------------------------------------+----+
   | 5          | H2a        | :math:`y_a y_b\geq 0` and :math:`y_a+y_b<0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | N2a        | :math:`y_a y_b\geq 0` and :math:`y_a+y_b>0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | H2b        | :math:`y_a y_b<0` and :math:`\tilde{y}_2<0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   |            | N2b        | :math:`y_a y_b<0` and :math:`\tilde{y}_2>0`            |    |
   +------------+------------+--------------------------------------------------------+----+
   +------------+------------+--------------------------------------------------------+----+

This scheme was originally designed for rectangular grids. Grid cells in
CICE actually lie on the surface of a sphere and must be projected onto
a plane. The projection used in CICE maps each grid cell to a square
with sides of unit length. Departure triangles across a given cell edge
are computed in a coordinate system whose origin lies at the midpoint of
the edge and whose vertices are at (-0.5, 0) and (0.5, 0). Intersection
points are computed assuming Cartesian geometry with cell edges meeting
at right angles. Let CL and CR denote the left and right vertices, which
are joined by line CLR. Similarly, let DL and DR denote the departure
points, which are joined by line DLR. Also, let IL and IR denote the
intersection points (0, :math:`y_a`) and (0, :math:`y_b`) respectively,
and let IC = (:math:`x_c`, 0) denote the intersection of CLR and DLR. It
can be shown that :math:`y_a`, :math:`y_b`, and :math:`x_c` are given by

.. math::
   \begin{aligned}
    y_a &=& {x_{CL} (y_{DM}-y_{DL}) + x_{DM}y_{DL} - x_{DL}y_{DM}}\over{x_{DM} - x_{DL}}, \\
    y_b &=& {x_{CR} (y_{DR}-y_{DM}) - x_{DM}y_{DR} + x_{DR}y_{DM}}\over{x_{DR} - x_{DM}}, \\
    x_c &=& x_{DL} - y_{DL} \left({x_{DR} - x_{DL}} \over y_{DR} - y_{DL}\right)
    \end{aligned}

Each departure triangle is defined by three of the seven points (CL,
CR, DL, DR, IL, IR, IC).

Given a 2D velocity field **u**, the divergence
:math:`\nabla\cdot{\bf u}` in a given grid cell can be computed from the
local velocities and written in terms of fluxes across each cell edge:

.. math::
    \nabla\cdot{\bf u} = {1\over A}\left[\left({u_{NE}+u_{SE}}\over 2\right)L_E + \left({u_{NW}+u_{SW}}\over 2\right)L_W + \left({u_{NE}+u_{NW}}\over 2\right)L_N + \left({u_{SE}+u_{SW}}\over 2\right)L_S \right],
   :label: divergence

where :math:`L` is an edge length and the indices :math:`N, S, E, W`
denote compass directions. Equation :eq:`divergence` is equivalent to
the divergence computed in the EVP dynamics (Section :ref:`dynam`).
In general, the fluxes in this expression are not equal to those implied
by the above scheme for locating departure regions. For some
applications it may be desirable to prescribe the divergence by
prescribing the area of the departure region for each edge. This can be
done in CICE 4.0 by setting `l\_fixed\_area` = true in
**ice\_transport\_driver.F90** and passing the prescribed departure
areas (`edgearea\_e` and `edgearea\_n`) into the remapping routine. An extra
triangle is then constructed for each departure region to ensure that
the total area is equal to the prescribed value. This idea was suggested
and first implemented by Mats Bentsen of the Nansen Environmental and
Remote Sensing Center (Norway), who applied an earlier version of the
CICE remapping scheme to an ocean model. The implementation in CICE 4.0
is somewhat more general, allowing for departure regions lying on both
sides of a cell edge. The extra triangle is constrained to lie in one
but not both of the grid cells that share the edge. Since this option
has yet to be fully tested in CICE, the current default is
`l\_fixed\_area` = false.

We made one other change in the scheme of :cite:`Dukowicz00` for
locating triangles. In their paper, departure points are defined by
projecting cell corner velocities directly backward. That is,

.. math::
   \mathbf{x_D} = -\mathbf{u} \, \Delta t,
   :label: departure_points
  
where :math:`\mathbf{x}_D` is the location of the departure point
relative to the cell corner and :math:`\mathbf{u}` is the velocity at
the corner. This approximation is only first-order accurate. Accuracy
can be improved by estimating the velocity at the midpoint of the
trajectory.

.. _integ-flux:

******************
Integrating fields
******************

Next, we integrate the reconstructed fields over the departure triangles
to find the total area, volume, and energy transported across each cell
edge. Area transports are easy to compute since the area is linear in
:math:`x` and :math:`y`. Given a triangle with vertices
:math:`\mathbf{x_i} = (x_i,y_i)`, :math:`i\in\{1,2,3\}`, the triangle
area is

.. math::
   A_T = \frac{1}{2}\left|(x_2-x_1)(y_3-y_1) -
   (y_2-y_1)(x_3-x_1)\right|.

The integral :math:`F_a` of any linear function :math:`f(\mathbf{r})`
over a triangle is given by

.. math::
    F_a = A_T f(\mathbf{x_0}),
   :label: I1

where :math:`\mathbf{x}_0 = (x_0,y_0)` is the triangle midpoint,

.. math::
   \mathbf{x}_0={1\over 3}\sum_{i=1}^3\mathbf{x}_i.

To compute the area transport, we evaluate the area at the midpoint,

.. math:: 
   a(\mathbf{x}_0)  = a_c + a_x x_0 + a_y y_0,

and multiply by :math:`A_T`. By convention, northward and eastward
transport is positive, while southward and westward transport is
negative.

Equation :eq:`I1` cannot be used for volume transport, because the
reconstructed volumes are quadratic functions of position. (They are
products of two linear functions, area and thickness.) The integral of a
quadratic polynomial over a triangle requires function evaluations at
three points,

.. math::
    F_h = \frac{A_T}{3}\sum_{i=1}^3 f\left({\mathbf x}^\prime_i\right),
    :label: I2

where :math:`\mathbf{x}_i^\prime = (\mathbf{x}_0+\mathbf{x}_i)/2` are
points lying halfway between the midpoint and the three vertices.
:cite:`Dukowicz00` use this formula to compute transports of the
product :math:`\rho \, T`, which is analogous to ice volume. Equation
:eq:`I2` does not work for ice and snow energies, which are cubic
functions—products of area, thickness, and enthalpy. Integrals of a
cubic polynomial over a triangle can be evaluated using a four-point
formula :cite:`Stroud71`:

.. math::
    F_q = A_T \left[ -\frac{9}{16} f(\mathbf{x}_0) +
                 \frac{25}{48} \sum_{i=1}^3 f(\mathbf{x}_i^{\prime\prime})\right]
    :label: I3

where :math:`\mathbf{x_i}^{\prime\prime}=(3 \mathbf{x}_0 + 2
\mathbf{x}_i)/5`. To evaluate functions at specific points, we must
compute many products of the form :math:`a({\bf x}) \, h({\bf x})` and
:math:`a({\bf x}) \, h({\bf x}) \, q({\bf x})`, where each term in the
product is the sum of a cell-center value and two displacement terms. In
the code, the computation is sped up by storing some sums that are used
repeatedly.

.. _updating-state-var:

************************
Updating state variables
************************

Finally, we compute new values of the state variables in each ice
category and grid cell. The new fractional ice areas
:math:`a_{in}^\prime(i,j)` are given by

.. math::
   a_{in}^\prime(i,j) = a_{in}(i,j) +
                 \frac{F_{aE}(i-1,j) - F_{aE}(i,j)
                     + F_{aN}(i,j-1) - F_{aN}(i,j)}
                      {A(i,j)}
   :label: new-area

where :math:`F_{aE}(i,j)` and :math:`F_{aN}(i,j)` are the area
transports across the east and north edges, respectively, of cell
:math:`(i,j)`, and :math:`A(i,j)` is the grid cell area. All transports
added to one cell are subtracted from a neighboring cell; thus
Equation :eq:`new-area` conserves total ice area.

The new ice volumes and energies are computed analogously. New
thicknesses are given by the ratio of volume to area, and enthalpies by
the ratio of energy to volume. Tracer monotonicity is ensured because

.. math:: 
   h^\prime = {\int_A a \, h \, dA \over \int_A a \, dA},

.. math:: 
   q^\prime  = {\int_A a \, h \, q\,dA \over \int_A a \, h \ dA},

where :math:`h^\prime` and :math:`q^\prime` are the new-time thickness
and enthalpy, given by integrating the old-time ice area, volume, and
energy over a Lagrangian departure region with area :math:`A`. That is,
the new-time thickness and enthalpy are weighted averages over old-time
values, with non-negative weights :math:`a` and :math:`ah`. Thus the
new-time values must lie between the maximum and minimum of the old-time
values.

.. _mech-red:

*************************
Mechanical redistribution
*************************

The last term on the right-hand side of Equation :eq:`transport-g`
is :math:`\psi`, which describes the redistribution
of ice in thickness space due to ridging and other mechanical processes.
The mechanical redistribution scheme in CICE is based on
:cite:`Thorndike75`, :cite:`Rothrock75`,
:cite:`Hibler80`, :cite:`Flato95`, and
:cite:`Lipscomb07`. This scheme converts thinner ice to thicker
ice and is applied after horizontal transport. When the ice is
converging, enough ice ridges to ensure that the ice area does not
exceed the grid cell area.

First we specify the participation function: the thickness distribution
:math:`a_P(h) = b(h) \, g(h)` of the ice participating in ridging. (We
use “ridging” as shorthand for all forms of mechanical redistribution,
including rafting.) The weighting function :math:`b(h)` favors ridging
of thin ice and closing of open water in preference to ridging of
thicker ice. There are two options for the form of :math:`b(h)`. If
`krdg\_partic` = 0 in the namelist, we follow :cite:`Thorndike75`
and set

.. math::
   b(h) = \left\{\begin{array}{ll}  
          \frac{2}{G^*}(1-\frac{G(h)}{G^*}) & \mbox{if $G(h)<G^*$} \\
                    0                       & \mbox{otherwise}   
                 \end{array}  \right.
   :label: partic-old-contin

where :math:`G(h)` is the fractional area covered by ice thinner than
:math:`h`, and :math:`G^*` is an empirical constant. Integrating
:math:`a_P(h)` between category boundaries :math:`H_{n-1}` and
:math:`H_n`, we obtain the mean value of :math:`a_P` in category
:math:`n`:

.. math::
   a_{Pn} = \frac{2}{G^*} (G_n - G_{n-1})
            \left( 1 - \frac{G_{n-1}+G_n}{2 G^*} \right),
   :label: partic-old-discrete

where :math:`a_{Pn}` is the ratio of the ice area ridging (or open
water area closing) in category :math:`n` to the total area ridging and
closing, and :math:`G_n` is the total fractional ice area in categories
0 to :math:`n`. Equation :eq:`partic-old-discrete` applies to
categories with :math:`G_n < G^*`. If :math:`G_{n-1} < G^* < G_n`, then
Equation :eq:`partic-old-discrete` is valid with :math:`G^*` replacing
:math:`G_n`, and if :math:`G_{n-1} > G^*`, then :math:`a_{Pn} = 0`. If
the open water fraction :math:`a_0 > G^*`, no ice can ridge, because
“ridging” simply reduces the area of open water. As in
:cite:`Thorndike75` we set :math:`G^* = 0.15`.

If the spatial resolution is too fine for a given time step
:math:`\Delta t`, the weighting function Equation :eq:`partic-old-contin` can
promote numerical instability. For :math:`\Delta t = \mbox{1 hour}`,
resolutions finer than :math:`\Delta x \sim \mbox{10 km}` are typically
unstable. The instability results from feedback between the ridging
scheme and the dynamics via the ice strength. If the strength changes
significantly on time scales less than :math:`\Delta t`, the
viscous-plastic solution of the momentum equation is inaccurate and
sometimes oscillatory. As a result, the fields of ice area, thickness,
velocity, strength, divergence, and shear can become noisy and
unphysical.

A more stable weighting function was suggested by
:cite:`Lipscomb07`:

.. math::
   b(h) = \frac{\exp[-G(h)/a^*]}
               {a^*[1-\exp(-1/a^*)]}
   :label: partic-new-contin

When integrated between category boundaries, Equation :eq:`partic-new-contin`
implies

.. math::
   a_{Pn} = \frac {\exp(-G_{n-1}/a^*) - \exp(-G_{n}/a^*)}
                  {1 - \exp(-1/a^*)}
   :label: partic-new-discrete

This weighting function is used if `krdg\_partic` = 1 in the namelist.
From Equation :eq:`partic-new-contin`, the mean value of :math:`G` for ice
participating in ridging is :math:`a^*`, as compared to :math:`G^*/3`
for Equation :eq:`partic-old-contin`. For typical ice thickness distributions,
setting :math:`a^* = 0.05` with `krdg\_partic` = 1 gives participation
fractions similar to those given by :math:`G^* = 0.15` with `krdg\_partic`
= 0. See :cite:`Lipscomb07` for a detailed comparison of these
two participation functions.

Thin ice is converted to thick, ridged ice in a way that reduces the
total ice area while conserving ice volume and internal energy. There
are two namelist options for redistributing ice among thickness
categories. If `krdg\_redist` = 0, ridging ice of thickness :math:`h_n`
forms ridges whose area is distributed uniformly between
:math:`H_{\min} = 2 h_n` and :math:`H_{\max} = 2 \sqrt{H^* h_n}`, as in
:cite:`Hibler80`. The default value of :math:`H^*` is 25 m, as
in earlier versions of CICE. Observations suggest that
:math:`H^* = 50` m gives a better fit to first-year ridges
:cite:`Amundrud04`, although the lower value may be appropriate
for multiyear ridges :cite:`Flato95`. The ratio of the mean
ridge thickness to the thickness of ridging ice is
:math:`k_n = (H_{\min} + H_{\max}) / (2 h_n)`. If the area of category
:math:`n` is reduced by ridging at the rate :math:`r_n`, the area of
thicker categories grows simultaneously at the rate :math:`r_n/k_n`.
Thus the *net* rate of area loss due to ridging of ice in category
:math:`n` is :math:`r_n(1-1/k_n)`.

The ridged ice area and volume are apportioned among categories in the
thickness range :math:`(H_{\min}, H_{\max})`. The fraction of the new
ridge area in category :math:`m` is

.. math::
   f_m^{\mathrm{area}} = \frac{H_R - H_L} 
                              {H_{\max} - H_{\min}},
   :label: ridge-area-old

where :math:`H_L = \max(H_{m-1},H_{\min})` and
:math:`H_R= \min(H_m,H_{\max})`. The fraction of the ridge volume going
to category :math:`m` is

.. math::
   f_m^{\mathrm{vol}} = \frac{(H_R)^2 - (H_L)^2}
                             {(H_{\max})^2 - (H_{\min})^2}.
   :label: ridge-volume-old

This uniform redistribution function tends to produce too little ice in
the 3–5 m range and too much ice thicker than 10 m
:cite:`Amundrud04`. Observations show that the ITD of ridges is
better approximated by a negative exponential. Setting `krdg\_redist` = 1
gives ridges with an exponential ITD :cite:`Lipscomb07`:

.. math::
   g_R(h) \propto \exp[-(h - H_{\min})/\lambda]
   :label: redist-new

for :math:`h \ge H_{\min}`, with :math:`g_R(h) = 0` for
:math:`h < H_{\min}`. Here, :math:`\lambda` is an empirical *e*-folding
scale and :math:`H_{\min}=2h_n` (where :math:`h_n` is the thickness of
ridging ice). We assume that :math:`\lambda = \mu h_n^{1/2}`, where
:math:`\mu` (mu\_rdg) is a tunable parameter with units . Thus the mean
ridge thickness increases in proportion to :math:`h_n^{1/2}`, as in
:cite:`Hibler80`. The value :math:`\mu = 4.0`  gives
:math:`\lambda` in the range 1–4 m for most ridged ice. Ice strengths
with :math:`\mu = 4.0`  and `krdg\_redist` = 1 are roughly comparable to
the strengths with :math:`H^* = 50` m and `krdg\_redist` = 0.

From Equation :eq:`redist-new` it can be shown that the fractional area going
to category :math:`m` as a result of ridging is

.. math::
   f_m^{\mathrm{area}} = \exp[-(H_{m-1} - H_{\min}) / \lambda] 
                        - \exp[-(H_m - H_{\min}) / \lambda].
   :label: ridge-area-new

The fractional volume going to category :math:`m` is

.. math::
   f_m^{\mathrm{vol}} = \frac{(H_{m-1}+\lambda) \exp[-(H_{m-1}-H_{\min})/\lambda]
                              - (H_m + \lambda) \exp[-(H_m - H_{\min}) / \lambda]}
                                {H_{min} + \lambda}.
   :label: ridge-volume-new

Equations :eq:`ridge-area-new` and :eq:`ridge-volume-new` replace
Equations :eq:`ridge-area-old` and :eq:`ridge-volume-old` when `krdg\_redist`
= 1.

Internal ice energy is transferred between categories in proportion to
ice volume. Snow volume and internal energy are transferred in the same
way, except that a fraction of the snow may be deposited in the ocean
instead of added to the new ridge.

The net area removed by ridging and closing is a function of the strain
rates. Let :math:`R_{\mathrm{net}}` be the net rate of area loss for the
ice pack (i.e., the rate of open water area closing, plus the net rate
of ice area loss due to ridging). Following :cite:`Flato95`,
:math:`R_{\mathrm{net}}` is given by

.. math::
   R_{\mathrm{net}} = \frac{C_s}{2}
                    (\Delta - |D_D|) - \min(D_D,0),
   :label: Rnet

where :math:`C_s` is the fraction of shear dissipation energy that
contributes to ridge-building, :math:`D_D` is the divergence, and
:math:`\Delta` is a function of the divergence and shear. These strain
rates are computed by the dynamics scheme. The default value of
:math:`C_s` is 0.25.

Next, define :math:`R_{\mathrm{tot}} = \sum_{n=0}^N r_n`. This rate is
related to :math:`R_{\mathrm{net}}` by

.. math::
   R_{\mathrm{net}} =
      \left[ a_{P0} + \sum_{n=1}^N a_{Pn}\left(1-{1\over k_n}\right)\right]
       R_{\mathrm{tot}}.
   :label: Rtot-Rnet

Given :math:`R_{\mathrm{net}}` from Equation :eq:`Rnet`, we
use Equation :eq:`Rtot-Rnet` to compute :math:`R_{\mathrm{tot}}`. Then the area
ridged in category :math:`n` is given by :math:`a_{rn} = r_n \Delta t`,
where :math:`r_n = a_{Pn} R_{\mathrm{tot}}`. The area of new ridges is
:math:`a_{rn} / k_n`, and the volume of new ridges is :math:`a_{rn} h_n`
(since volume is conserved during ridging). We remove the ridging ice
from category :math:`n` and use Equations :eq:`ridge-area-old`
and :eq:`ridge-volume-old`: (or :eq:`ridge-area-new` and
:eq:`ridge-volume-new`) to redistribute the ice among thicker
categories.

Occasionally the ridging rate in thickness category :math:`n` may be
large enough to ridge the entire area :math:`a_n` during a time interval
less than :math:`\Delta t`. In this case :math:`R_{\mathrm{tot}}` is
reduced to the value that exactly ridges an area :math:`a_n` during
:math:`\Delta t`. After each ridging iteration, the total fractional ice
area :math:`a_i` is computed. If :math:`a_i > 1`, the ridging is
repeated with a value of :math:`R_{\mathrm{net}}` sufficient to yield
:math:`a_i = 1`.

Two tracers for tracking the ridged ice area and volume are available.
The actual tracers are for level (undeformed) ice area (`alvl`) and volume
(`vlvl`), which are easier to implement for a couple of reasons: (1) ice
ridged in a given thickness category is spread out among the rest of the
categories, making it more difficult (and expensive) to track than the
level ice remaining behind in the original category; (2) previously
ridged ice may ridge again, so that simply adding a volume of freshly
ridged ice to the volume of previously ridged ice in a grid cell may be
inappropriate. Although the code currently only tracks level ice
internally, both level ice and ridged ice are offered as history output.
They are simply related:

.. math::
   \begin{aligned}
   a_{lvl} + a_{rdg} &=& a_i, \\
   v_{lvl} + v_{rdg} &=& v_i.\end{aligned}

Level ice area fraction and volume increase with new ice formation and
decrease steadily via ridging processes. Without the formation of new
ice, level ice asymptotes to zero because we assume that both level ice
and ridged ice ridge, in proportion to their fractional areas in a grid
cell (in the spirit of the ridging calculation itself which does not
prefer level ice over previously ridged ice).

The ice strength :math:`P` may be computed in either of two ways. If the
namelist parameter kstrength = 0, we use the strength formula from
:cite:`Hibler79`:

.. math::
   P = P^* h \exp[-C(1-a_i)],
   :label: hib-strength

where :math:`P^* = 27,500 \, \mathrm {N/m}` and :math:`C = 20` are
empirical constants, and :math:`h` is the mean ice thickness.
Alternatively, setting kstrength = 1 gives an ice strength closely
related to the ridging scheme. Following
:cite:`Rothrock75`, the strength is assumed proportional
to the change in ice potential energy :math:`\Delta E_P` per unit area
of compressive deformation. Given uniform ridge ITDs (krdg\_redist = 0),
we have

.. math::
   P = C_f \, C_p \, \beta \sum_{n=1}^{N_C}
     \left[ -a_{Pn} \, h_n^2  + \frac{a_{Pn}}{k_n}
        \left( \frac{(H_n^{\max})^3 - (H_n^{\min})^3}
                    {3(H_n^{\max}-H_n^{\min})} \right) \right],
   :label: roth-strength0

where :math:`C_P = (g/2)(\rho_i/\rho_w)(\rho_w-\rho_i)`,
:math:`\beta =R_{\mathrm{tot}}/R_{\mathrm{net}} > 1`
from Equation :eq:`Rtot-Rnet`, and :math:`C_f` is an empirical parameter that
accounts for frictional energy dissipation. Following
:cite:`Flato95`, we set :math:`C_f = 17`. The first term in
the summation is the potential energy of ridging ice, and the second,
larger term is the potential energy of the resulting ridges. The factor
of :math:`\beta` is included because :math:`a_{Pn}` is normalized with
respect to the total area of ice ridging, not the net area removed.
Recall that more than one unit area of ice must be ridged to reduce the
net ice area by one unit. For exponential ridge ITDs (`krdg\_redist` = 1),
the ridge potential energy is modified:

.. math::
   P = C_f \, C_p \, \beta \sum_{n=1}^{N_C}
     \left[ -a_{Pn} \, h_n^2  + \frac{a_{Pn}}{k_n}
        \left( H_{\min}^2 + 2H_{\min}\lambda + 2 \lambda^2 \right) \right]
   :label: roth-strength1

The energy-based ice strength given by Equations :eq:`roth-strength0` or
:eq:`roth-strength1` is more physically realistic than the strength
given by Equation :eq:`hib-strength`. However, use of Equation :eq:`hib-strength` is
less likely to allow numerical instability at a given resolution and
time step. See :cite:`Lipscomb07` for more details.


.. _dynam:

~~~~~~~~
Dynamics
~~~~~~~~

There are now different rheologies available in the CICE code. The
elastic-viscous-plastic (EVP) model represents a modification of the
standard viscous-plastic (VP) model for sea ice dynamics
:cite:`Hibler79`. The elastic-anisotropic-plastic (EAP) model,
on the other hand, explicitly accounts for the observed sub-continuum
anisotropy of the sea ice cover :cite:`Wilchinsky06,Weiss09`. If
`kdyn` = 1 in the namelist then the EVP rheology is used (module
**ice\_dyn\_evp.F90**), while `kdyn` = 2 is associated with the EAP
rheology (**ice\_dyn\_eap.F90**). At times scales associated with the
wind forcing, the EVP model reduces to the VP model while the EAP model
reduces to the anisotropic rheology described in detail in
:cite:`Wilchinsky06,Tsamados13`. At shorter time scales the
adjustment process takes place in both models by a numerically more
efficient elastic wave mechanism. While retaining the essential physics,
this elastic wave modification leads to a fully explicit numerical
scheme which greatly improves the model’s computational efficiency.

The EVP sea ice dynamics model is thoroughly documented in
:cite:`Hunke97`, :cite:`Hunke01`,
:cite:`Hunke02` and :cite:`Hunke03` and the EAP
dynamics in :cite:`Tsamados13`. Simulation results and
performance of the EVP and EAP models have been compared with the VP
model and with each other in realistic simulations of the Arctic
respectively in :cite:`Hunke99` and
:cite:`Tsamados13`. Here we summarize the equations and
direct the reader to the above references for details. The numerical
implementation in this code release is that of :cite:`Hunke02`
and :cite:`Hunke03`, with revisions to the numerical solver as
in :cite:`Bouillon13`. The implementation of the EAP sea ice
dynamics into CICE is described in detail in
:cite:`Tsamados13`.

.. _momentum:

********
Momentum
********

The force balance per unit area in the ice pack is given by a
two-dimensional momentum equation :cite:`Hibler79`, obtained
by integrating the 3D equation through the thickness of the ice in the
vertical direction:

.. math::
   m{\partial {\bf u}\over\partial t} = \nabla\cdot{\bf \sigma}
   + \vec{\tau}_a+\vec{\tau}_w + \vec{\tau}_b - \hat{k}\times mf{\bf u} - mg\nabla H_\circ,
   :label: vpmom

where :math:`m` is the combined mass of ice and snow per unit area and
:math:`\vec{\tau}_a` and :math:`\vec{\tau}_w` are wind and ocean
stresses, respectively. The term :math:`\vec{\tau}_b` is a 
seabed stress (also referred to as basal stress) that represents the grounding of pressure
ridges in shallow water :cite:`Lemieux16`. The mechanical properties of the ice are represented by the
internal stress tensor :math:`\sigma_{ij}`. The other two terms on
the right hand side are stresses due to Coriolis effects and the sea
surface slope. The parameterization for the wind and ice–ocean stress
terms must contain the ice concentration as a multiplicative factor to
be consistent with the formal theory of free drift in low ice
concentration regions. A careful explanation of the issue and its
continuum solution is provided in :cite:`Hunke03` and
:cite:`Connolley04`.

The momentum equation is discretized in time as follows, for the classic
EVP approach. First, for clarity, the two components of Equation :eq:`vpmom` are

.. math::
   \begin{aligned}
   m{\partial u\over\partial t} &=& {\partial\sigma_{1j}\over\partial x_j} + \tau_{ax} + 
     a_i c_w \rho_w
     \left|{\bf U}_w - {\bf u}\right| \left[\left(U_w-u\right)\cos\theta - \left(V_w-v\right)\sin\theta\right]
     -C_bu +mfv - mg{\partial H_\circ\over\partial x}, \\
   m{\partial v\over\partial t} &=& {\partial\sigma_{2j}\over\partial x_j} + \tau_{ay} + 
     a_i c_w \rho_w
     \left|{\bf U}_w - {\bf u}\right| \left[\left(U_w-u\right)\sin\theta - \left(V_w-v\right)\cos\theta\right]
     -C_bv-mfu - mg{\partial H_\circ\over\partial y}. \end{aligned}

In the code,
:math:`{\tt vrel}=a_i c_w \rho_w\left|{\bf U}_w - {\bf u}^k\right|` and 
:math:`C_b=T_b \left( \sqrt{(u^k)^2+(v^k)^2}+u_0 \right)^{-1}`, 
where :math:`k` denotes the subcycling step. The following equations
illustrate the time discretization and define some of the other
variables used in the code.

.. math::
   \underbrace{\left({m\over\Delta t_e}+{\tt vrel} \cos\theta\ + C_b \right)}_{\tt cca} u^{k+1} 
   - \underbrace{\left(mf+{\tt vrel}\sin\theta\right)}_{\tt ccb}v^{k+1}
    =  \underbrace{{\partial\sigma_{1j}^{k+1}\over\partial x_j}}_{\tt strintx} 
    + \underbrace{\tau_{ax} - mg{\partial H_\circ\over\partial x} }_{\tt forcex}
     + {\tt vrel}\underbrace{\left(U_w\cos\theta-V_w\sin\theta\right)}_{\tt waterx}  + {m\over\Delta t_e}u^k,
   :label: umom

.. math::
    \underbrace{\left(mf+{\tt vrel}\sin\theta\right)}_{\tt ccb} u^{k+1} 
   + \underbrace{\left({m\over\Delta t_e}+{\tt vrel} \cos\theta + C_b \right)}_{\tt cca}v^{k+1}
    =  \underbrace{{\partial\sigma_{2j}^{k+1}\over\partial x_j}}_{\tt strinty} 
    + \underbrace{\tau_{ay} - mg{\partial H_\circ\over\partial y} }_{\tt forcey}
     + {\tt vrel}\underbrace{\left(U_w\sin\theta+V_w\cos\theta\right)}_{\tt watery}  + {m\over\Delta t_e}v^k,
   :label: vmom

and vrel\ :math:`\cdot`\ waterx(y) = taux(y).

We solve this system of equations analytically for :math:`u^{k+1}` and
:math:`v^{k+1}`. Define

.. math::
   \hat{u} = F_u + \tau_{ax} - mg{\partial H_\circ\over\partial x} + {\tt vrel} \left(U_w\cos\theta - V_w\sin\theta\right) + {m\over\Delta t_e}u^k 
   :label: cevpuhat

.. math::
   \hat{v} = F_v + \tau_{ay} - mg{\partial H_\circ\over\partial y} + {\tt vrel} \left(U_w\sin\theta + V_w\cos\theta\right) + {m\over\Delta t_e}v^k,
   :label: cevpvhat

where :math:`{\bf F} = \nabla\cdot\sigma^{k+1}`. Then

.. math::
   \begin{aligned}
   \left({m\over\Delta t_e} +{\tt vrel}\cos\theta\ + C_b \right)u^{k+1} - \left(mf + {\tt vrel}\sin\theta\right) v^{k+1} &=& \hat{u}  \\
   \left(mf + {\tt vrel}\sin\theta\right) u^{k+1} + \left({m\over\Delta t_e} +{\tt vrel}\cos\theta + C_b \right)v^{k+1} &=& \hat{v}.\end{aligned}

Solving simultaneously for :math:`u^{k+1}` and :math:`v^{k+1}`,

.. math::
   \begin{aligned}
   u^{k+1} = {a \hat{u} + b \hat{v} \over a^2 + b^2} \\
   v^{k+1} = {a \hat{v} - b \hat{u} \over a^2 + b^2}, \end{aligned}

where

.. math::
   a = {m\over\Delta t_e} + {\tt vrel}\cos\theta + C_b \\
   :label: cevpa

.. math::
   b = mf + {\tt vrel}\sin\theta.
   :label: cevpb

When the subcycling is finished for each (thermodynamic) time step, the
ice–ocean stress must be constructed from `taux(y)` and the terms
containing `vrel` on the left hand side of the equations.

The Hibler-Bryan form for the ice-ocean stress :cite:`Hibler87`
is included in **ice\_dyn\_shared.F90** but is currently commented out,
pending further testing.

.. _seabed-stress:

***************
Seabed stress
***************

The parameterization for the seabed stress is described in :cite:`Lemieux16`. The components of the basal seabed stress are 
:math:`\tau_{bx}=C_bu` and :math:`\tau_{by}=C_bv`, where :math:`C_b` is a coefficient expressed as

.. math::
   C_b= k_2 \max [0,(h_u - h_{cu})]  e^{-\alpha_b * (1 - a_u)} (\sqrt{u^2+v^2}+u_0)^{-1}, \\
   :label: Cb 

where :math:`k_2` determines the maximum seabed stress that can be sustained by the grounded parameterized ridge(s), :math:`u_0` 
is a small residual velocity and :math:`\alpha_b=20` is a parameter to ensure that the seabed stress quickly drops when 
the ice concentration is smaller than 1. In the code, :math:`k_2 \max [0,(h_u - h_{cu})]  e^{-\alpha_b * (1 - a_u)}` is defined as 
:math:`T_b`. The quantities :math:`h_u`, :math:`a_{u}` and :math:`h_{cu}` are calculated at 
the 'u' point based on local ice conditions (surrounding tracer points). They are respectively given by 

.. math::
   h_u=\max[v_i(i,j),v_i(i+1,j),v_i(i,j+1),v_i(i+1,j+1)], \\
   :label: hu 
   
.. math::
   a_u=\max[a_i(i,j),a_i(i+1,j),a_i(i,j+1),a_i(i+1,j+1)]. \\
   :label: au      
   
.. math::
   h_{cu}=a_u h_{wu} / k_1, \\
   :label: hcu

where the :math:`a_i` and :math:`v_i` are the total ice concentrations and ice volumes around the :math:`u` point :math:`i,j` and 
:math:`k_1` is a parameter that defines the critical ice thickness :math:`h_{cu}` at which the parameterized 
ridge(s) reaches the seafloor for a water depth :math:`h_{wu}=\min[h_w(i,j),h_w(i+1,j),h_w(i,j+1),h_w(i+1,j+1)]`. Given the formulation of :math:`C_b` in equation :eq:`Cb`, the seabed stress components are non-zero only 
when :math:`h_u > h_{cu}`. 

The maximum seabed stress depends on the weigth of the ridge 
above hydrostatic balance and the value of :math:`k_2`. It is, however, the parameter :math:`k_1` that has the most notable impact on the simulated extent of landfast ice. 
The value of :math:`k_1` can be changed at runtime using the namelist variable `k1`. The grounding scheme can be turned on or off using the namelist logical basalstress. 

Note that the user must provide a bathymetry field for using this grounding 
scheme. Grounding occurs up to water depth of ~25 m. It is suggested to have a bathymetry field with water depths larger than 5 m that represents well shallow water regions such as the Laptev Sea and the East Siberian Sea. 

   
.. _internal-stress:

***************
Internal stress
***************

For convenience we formulate the stress tensor :math:`\bf \sigma` in
terms of :math:`\sigma_1=\sigma_{11}+\sigma_{22}`,
:math:`\sigma_2=\sigma_{11}-\sigma_{22}`, and introduce the
divergence, :math:`D_D`, and the horizontal tension and shearing
strain rates, :math:`D_T` and :math:`D_S` respectively.

CICE now outputs the internal ice pressure which is an important field to support navigation in ice-infested water.
The internal ice pressure :math:`(sigP)` is the average of the normal stresses multiplied by :math:`-1` and 
is therefore simply equal to :math:`-\sigma_1/2`.

*Elastic-Viscous-Plastic*

In the EVP model the internal stress tensor is determined from a
regularized version of the VP constitutive law. Following the approach of :cite:`Konig10` (see also :cite:`Lemieux16`), the 
elliptical yield curve can be modified such that the ice has isotropic tensile strength. 
The tensile strength :math:`T_p` is expressed as a fraction of the ice strength :math:`P`, that is :math:`T_p=k_t P` 
where :math:`k_t` should be set to a value between 0 and 1. The constitutive law is therefore 

.. math::
   {1\over E}{\partial\sigma_1\over\partial t} + {\sigma_1\over 2\zeta} 
     + {P_R(1-k_t)\over 2\zeta} = D_D, \\
   :label: sig1 

.. math::
   {1\over E}{\partial\sigma_2\over\partial t} + {\sigma_2\over 2\eta} = D_T,
   :label: sig2

.. math::
   {1\over E}{\partial\sigma_{12}\over\partial t} + {\sigma_{12}\over
     2\eta} = {1\over 2}D_S,
   :label: sig12

where

.. math::
   D_D = \dot{\epsilon}_{11} + \dot{\epsilon}_{22}, 

.. math::
   D_T = \dot{\epsilon}_{11} - \dot{\epsilon}_{22}, 

.. math::
   D_S = 2\dot{\epsilon}_{12}, 

.. math::
   \dot{\epsilon}_{ij} = {1\over 2}\left({{\partial u_i}\over{\partial x_j}} + {{\partial u_j}\over{\partial x_i}}\right), 

.. math::
   \zeta = {P(1+k_t)\over 2\Delta}, 

.. math::
   \eta  = {P(1+k_t)\over {2\Delta e^2}}, 

.. math::
   \Delta = \left[D_D^2 + {1\over e^2}\left(D_T^2 + D_S^2\right)\right]^{1/2},

and :math:`P_R` is a “replacement pressure” (see :cite:`Geiger98`, for
example), which serves to prevent residual ice motion due to spatial
variations of :math:`P` when the rates of strain are exactly zero. The ice strength :math:`P` 
is a function of the ice thickness and concentration
as it is described in Section :ref:`mech-red`.

Viscosities are updated during the subcycling, so that the entire
dynamics component is subcycled within the time step, and the elastic
parameter :math:`E` is defined in terms of a damping timescale :math:`T`
for elastic waves, :math:`\Delta t_e < T < \Delta t`, as

.. math:: 
   E = {\zeta\over T},

where :math:`T=E_\circ\Delta t` and :math:`E_\circ` (eyc) is a tunable
parameter less than one. Including the modification proposed by :cite:`Bouillon13` for equations :eq:`sig2` and :eq:`sig12` in order to improve numerical convergence, the stress equations become

.. math::
   \begin{aligned}
   {\partial\sigma_1\over\partial t} + {\sigma_1\over 2T} 
     + {P_R(1-k_t)\over 2T} &=& {P(1+k_t)\over 2T\Delta} D_D, \\
   {\partial\sigma_2\over\partial t} + {\sigma_2\over 2T} &=& {P(1+k_t)\over
     2Te^2\Delta} D_T,\\
   {\partial\sigma_{12}\over\partial t} + {\sigma_{12}\over  2T} &=&
     {P(1+k_t)\over 4Te^2\Delta}D_S.\end{aligned}

Once discretized in time, these last three equations are written as

.. math::
   \begin{aligned}
   {(\sigma_1^{k+1}-\sigma_1^{k})\over\Delta t_e} + {\sigma_1^{k+1}\over 2T} 
     + {P_R^k(1-k_t)\over 2T} &=& {P(1+k_t)\over 2T\Delta^k} D_D^k, \\
   {(\sigma_2^{k+1}-\sigma_2^{k})\over\Delta t_e} + {\sigma_2^{k+1}\over 2T} &=& {P(1+k_t)\over
     2Te^2\Delta^k} D_T^k,\\
   {(\sigma_{12}^{k+1}-\sigma_{12}^{k})\over\Delta t_e} + {\sigma_{12}^{k+1}\over  2T} &=&
     {P(1+k_t)\over 4Te^2\Delta^k}D_S^k,\end{aligned}
   :label: sigdisc  
     

where :math:`k` denotes again the subcycling step. All coefficients on the left-hand side are constant except for
:math:`P_R`. This modification compensates for the decreased efficiency of including
the viscosity terms in the subcycling. (Note that the viscosities do not
appear explicitly.) Choices of the parameters used to define :math:`E`,
:math:`T` and :math:`\Delta t_e` are discussed in
Sections :ref:`revp` and :ref:`parameters`.

The bilinear discretization used for the stress terms
:math:`\partial\sigma_{ij}/\partial x_j` in the momentum equation is
now used, which enabled the discrete equations to be derived from the
continuous equations written in curvilinear coordinates. In this
manner, metric terms associated with the curvature of the grid are
incorporated into the discretization explicitly. Details pertaining to
the spatial discretization are found in :cite:`Hunke02`.

*Elastic-Anisotropic-Plastic*

In the EAP model the internal stress tensor is related to the
geometrical properties and orientation of underlying virtual diamond
shaped floes (see :ref:`fig-EAP`). In contrast to the isotropic EVP
rheology, the anisotropic plastic yield curve within the EAP rheology
depends on the relative orientation of the diamond shaped floes (unit
vector :math:`\mathbf r` in :ref:`fig-EAP`), with respect to the
principal direction of the deformation rate (not shown). Local
anisotropy of the sea ice cover is accounted for by an additional
prognostic variable, the structure tensor :math:`\mathbf{A}` defined
by

.. math:: 
   {\mathbf A}=\int_{\mathbb{S}}\vartheta(\mathbf r)\mathbf r\mathbf r d\mathbf r\label{structuretensor}.

where :math:`\mathbb{S}` is a unit-radius circle; **A** is a unit
trace, 2\ :math:`\times`\ 2 matrix. From now on we shall describe the
orientational distribution of floes using the structure tensor. For
simplicity we take the probability density function
:math:`\vartheta(\mathbf r )` to be Gaussian,
:math:`\vartheta(z)=\omega_{1}\exp(-\omega_{2}z^{2})`, where :math:`z`
is the ice floe inclination with respect to the axis :math:`x_{1}` of
preferential alignment of ice floes (see :ref:`fig-EAP`),
:math:`\vartheta(z)` is periodic with period :math:`\pi`, and the
positive coefficients :math:`\omega_{1}` and :math:`\omega_{2}` are
calculated to ensure normalization of :math:`\vartheta(z)`, i.e.
:math:`\int_{0}^{2\pi}\vartheta(z)dz=1`. The ratio of the principal
components of :math:`\mathbf{A}`, :math:`A_{1}/A_{2}`, are derived
from the phenomenological evolution equation for the structure tensor
:math:`\mathbf A`,

.. math:: 
   \frac{D\mathbf{A}}{D t}=\mathbf{F}_{iso}(\mathbf{A})+\mathbf{F}_{frac}(\mathbf{A},\boldsymbol\sigma),
   :label: evolutionA

where :math:`t` is the time, and :math:`D/Dt` is the co-rotational
time derivative accounting for advection and rigid body rotation
(:math:`D\mathbf A/Dt = d\mathbf A/dt -\mathbf W \cdot \mathbf A -\mathbf A \cdot \mathbf W^{T}`)
with :math:`\mathbf W` being the vorticity tensor.
:math:`\mathbf F_{iso}` is a function that accounts for a variety of
processes (thermal cracking, melting, freezing together of floes) that
contribute to a more isotropic nature to the ice cover.
:math:`\mathbf F_{frac}` is a function determining the ice floe
re-orientation due to fracture, and explicitly depends upon sea ice
stress (but not its magnitude). Following :cite:`Wilchinsky06`,
based on laboratory experiments by :cite:`Schulson01` we
consider four failure mechanisms for the Arctic sea ice cover. These
are determined by the ratio of the principal values of the sea ice
stress :math:`\sigma_{1}` and :math:`\sigma_{2}`: (i) under biaxial
tension, fractures form across the perpendicular principal axes and
therefore counteract any apparent redistribution of the floe
orientation; (ii) if only one of the principal stresses is
compressive, failure occurs through axial splitting along the
compression direction; (iii) under biaxial compression with a low
confinement ratio, (:math:`\sigma_{1}/\sigma_{2}<R`), sea ice fails
Coulombically through formation of slip lines delineating new ice
floes oriented along the largest compressive stress; and finally (iv)
under biaxial compression with a large confinement ratio,
(:math:`\sigma_{1}/\sigma_{2}\ge R`), the ice is expected to fail
along both principal directions so that the cumulative directional
effect balances to zero.

.. _fig-EAP:

.. figure:: ./figures/EAP.png
   :align: center
   :scale: 15%

   Diamond-shaped floes

Figure :ref:`fig-EAP` shows geometry of interlocking diamond-shaped floes (taken from
:cite:`Wilchinsky06`). :math:`\phi` is half of the acute angle
of the diamonds. :math:`L` is the edge length.
:math:`\boldsymbol n_{1}`, :math:`\boldsymbol n_{2}` and
:math:`\boldsymbol\tau_{1}`, :math:`\boldsymbol\tau_{2}` are
respectively the normal and tangential unit vectors along the diamond edges.
:math:`\mathbf v=L\boldsymbol\tau_{2}\cdot\dot{\boldsymbol\epsilon}`
is the relative velocity between the two floes connected by the
vector :math:`L \boldsymbol \tau_{2}`. :math:`\mathbf r` is the unit
vector along the main diagonal of the diamond. Note that the diamonds
illustrated here represent one possible realisation of all possible
orientations. The angle :math:`z` represents the rotation of the
diamonds’ main axis relative to their preferential orientation along
the axis :math:`x_1`.

The new anisotropic rheology requires solving the evolution
Equation :eq:`evolutionA` for the structure tensor in addition to the momentum
and stress equations. The evolution equation for :math:`\mathbf{A}` is
solved within the EVP subcycling loop, and consistently with the
momentum and stress evolution equations, we neglect the advection term
for the structure tensor. Equation :eq:`evolutionA` then reduces to the system
of two equations:

.. math::
   \begin{aligned}
   \frac{\partial A_{11}}{\partial t}&=&-k_{t}\left(A_{11}-\frac{1}{2}\right)+M_{11}  \mbox{,} \\ 
   \frac{\partial A_{12}}{\partial t}&=&-k_{t} A_{12}+M_{12}  \mbox{,}\end{aligned}

where the first terms on the right hand side correspond to the
isotropic contribution, :math:`F_{iso}`, and :math:`M_{11}` and
:math:`M_{12}` are the components of the term :math:`F_{frac}` in
Equation :eq:`evolutionA` that are given in :cite:`Wilchinsky06` and
:cite:`Tsamados13`. These evolution equations are
discretized semi-implicitly in time. The degree of anisotropy is
measured by the largest eigenvalue (:math:`A_{1}`) of this tensor
(:math:`A_{2}=1-A_{1}`). :math:`A_{1}=1` corresponds to perfectly
aligned floes and :math:`A_{1}=0.5` to a uniform distribution of floe
orientation. Note that while we have specified the aspect ratio of the
diamond floes, through prescribing :math:`\phi`, we make no assumption
about the size of the diamonds so that formally the theory is scale
invariant.

As described in greater detail in :cite:`Wilchinsky06`, the
internal ice stress for a single orientation of the ice floes can be
calculated explicitly and decomposed, for an average ice thickness
:math:`h`, into its ridging (r) and sliding (s) contributions

.. math::
   \boldsymbol \sigma^{b}(\mathbf r,h)=P_{r}(h) \boldsymbol \sigma_{r}^{b}(\mathbf r)+P_{s}(h) \boldsymbol \sigma_{s}^{b}(\mathbf r),
   :label: stress1

where :math:`P_{r}` and :math:`P_{s}` are the ridging and sliding
strengths and the ridging and sliding stresses are functions of the
angle :math:`\theta= \arctan(\dot\epsilon_{II}/\dot\epsilon_{I})`, the
angle :math:`y` between the major principal axis of the strain rate
tensor (not shown) and the structure tensor (:math:`x_1` axis in
:ref:`fig-EAP`, and the angle :math:`z` defined in :ref:`fig-EAP`. In
the stress expressions above the underlying floes are assumed parallel,
but in a continuum-scale sea ice region the floes can possess different
orientations in different places and we take the mean sea ice stress
over a collection of floes to be given by the average

.. math:: 
   \boldsymbol\sigma^{EAP}(h)=P_{r}(h)\int_{\mathbb{S}}\vartheta(\mathbf r)\left[\boldsymbol\sigma_{r}^{b}(\mathbf r)+ k \boldsymbol\sigma_{s}^{b}(\mathbf r)\right]d\mathbf r
   :label: stressaverage

where we have introduced the friction parameter :math:`k=P_{s}/P_{r}`
and where we identify the ridging ice strength :math:`P_{r}(h)` with the
strength :math:`P` described in section 1 and used within the EVP
framework.

As is the case for the EVP rheology, elasticity is included in the EAP
description not to describe any physical effect, but to make use of the
efficient, explicit numerical algorithm used to solve the full sea ice
momentum balance. We use the analogous EAP stress equations,

.. math::
   \frac{\partial \sigma_{1}}{\partial t}+\frac{\sigma_1}{2T} = \frac{\sigma^{EAP}_{1}}{2T}  \mbox{,}  
   :label: EAPsigma1

.. math::
   \frac{\partial \sigma_{2}}{\partial t}+\frac{\sigma_2}{2T} = \frac{\sigma^{EAP}_{2}}{2T} \mbox{,}  
   :label: EAPsigma2

.. math::
   \frac{\partial \sigma_{12}}{\partial t}+\frac{\sigma_{12}}{2T} = \frac{\sigma^{EAP}_{12}}{2T} \mbox{,}
   :label: EAPsigma12

where the anisotropic stress :math:`\boldsymbol\sigma^{EAP}` is defined
in a look-up table for the current values of strain rate and structure
tensor. The look-up table is constructed by computing the stress
(normalized by the strength) from Equations :eq:`EAPsigma1`–:eq:`EAPsigma12`
for discrete values of the largest eigenvalue of the structure tensor,
:math:`\frac{1}{2}\le A_{1}\le 1`, the angle :math:`0\le\theta\le2\pi`,
and the angle :math:`-\pi/2\le y\le\pi/2` between the major principal
axis of the strain rate tensor and the structure tensor
:cite:`Tsamados13`. The updated stress, after the elastic
relaxation, is then passed to the momentum equation and the sea ice
velocities are updated in the usual manner within the subcycling loop of
the EVP rheology. The structure tensor evolution equations are solved
implicitly at the same frequency, :math:`\Delta t_{e}`, as the ice
velocities and internal stresses. Finally, to be coherent with our new
rheology we compute the area loss rate due to ridging as
:math:`\vert\dot{\boldsymbol\epsilon}\vert\alpha_{r}(\theta)`, with
:math:`\alpha_r(\theta)` and :math:`\alpha_s(\theta)` given by
:cite:`Wilchinsky04`,

.. math::
   \begin{aligned}
   \alpha_{r}(\theta)=\frac{\sigma^{r}_{ij}\dot\epsilon_{ij}}{P_{r} \vert\dot{\boldsymbol\epsilon}\vert } , \qquad \alpha_{s}(\theta)=\frac{\sigma^{s}_{ij}\dot\epsilon_{ij}}{P_{s} \vert\dot{\boldsymbol\epsilon}\vert }.\label{alphas}\end{aligned}

Both ridging rate and sea ice strength are computed in the outer loop
of the dynamics.

.. _revp:

****************
Revised approach
****************

The revised EVP approach is based on a pseudo-time iterative scheme :cite:`Lemieux12`, :cite:`Bouillon13`, :cite:`Kimmritz15`. By construction, the revised EVP approach should lead to the VP solution 
(given the right numerical parameters and a sufficiently large number of iterations). To do so, the inertial term is formulated such that it matches the backward Euler approach of 
implicit solvers and there is an additional term for the pseudo-time iteration. Hence, with the revised approach, the discretized momentum equations :eq:`umom` and :eq:`vmom` become  

.. math::
    {\beta^*(u^{k+1}-u^k)\over\Delta t_e} + {m(u^{k+1}-u^n)\over\Delta t} + {\left({\tt vrel} \cos\theta + C_b \right)} u^{k+1} 
    - {\left(mf+{\tt vrel}\sin\theta\right)} v^{k+1}
    = {{\partial\sigma_{1j}^{k+1}\over\partial x_j}} 
    + {\tau_{ax} - mg{\partial H_\circ\over\partial x} }
    + {\tt vrel} {\left(U_w\cos\theta-V_w\sin\theta\right)},
    :label: umomr

.. math::
    {\beta^*(v^{k+1}-v^k)\over\Delta t_e} + {m(v^{k+1}-v^n)\over\Delta t} + {\left({\tt vrel} \cos\theta + C_b \right)}v^{k+1} 
    + {\left(mf+{\tt vrel}\sin\theta\right)} u^{k+1} 
    = {{\partial\sigma_{2j}^{k+1}\over\partial x_j}} 
    + {\tau_{ay} - mg{\partial H_\circ\over\partial y} }
    + {\tt vrel}{\left(U_w\sin\theta+V_w\cos\theta\right)},
    :label: vmomr

where :math:`\beta^*` is a numerical parameter and :math:`u^n, v^n` are the components of the previous time level solution. 
With :math:`\beta=\beta^* \Delta t \left(  m \Delta t_e \right)^{-1}` :cite:`Bouillon13`, these equations can be written as
 
.. math::
   \underbrace{\left((\beta+1){m\over\Delta t}+{\tt vrel} \cos\theta\ + C_b \right)}_{\tt cca} u^{k+1} 
   - \underbrace{\left(mf+{\tt vrel}\sin\theta\right)}_{\tt ccb}v^{k+1}
    =  \underbrace{{\partial\sigma_{1j}^{k+1}\over\partial x_j}}_{\tt strintx} 
    + \underbrace{\tau_{ax} - mg{\partial H_\circ\over\partial x} }_{\tt forcex}
     + {\tt vrel}\underbrace{\left(U_w\cos\theta-V_w\sin\theta\right)}_{\tt waterx}  + {m\over\Delta t}(\beta u^k + u^n),
   :label: umomr2

.. math::
    \underbrace{\left(mf+{\tt vrel}\sin\theta\right)}_{\tt ccb} u^{k+1} 
   + \underbrace{\left((\beta+1){m\over\Delta t}+{\tt vrel} \cos\theta + C_b \right)}_{\tt cca}v^{k+1}
    =  \underbrace{{\partial\sigma_{2j}^{k+1}\over\partial x_j}}_{\tt strinty} 
    + \underbrace{\tau_{ay} - mg{\partial H_\circ\over\partial y} }_{\tt forcey}
     + {\tt vrel}\underbrace{\left(U_w\sin\theta+V_w\cos\theta\right)}_{\tt watery}  + {m\over\Delta t}(\beta v^k + v^n),
   :label: vmomr2  

At this point, the solutions :math:`u^{k+1}` and :math:`v^{k+1}` are obtained in the same manner as for the standard EVP approach (see equations :eq:`cevpuhat` to :eq:`cevpb`).

Introducing another numerical parameter :math:`\alpha=2T \Delta t_e ^{-1}` :cite:`Bouillon13`, the stress equations in :eq:`sigdisc` become

.. math::
   \begin{aligned}
   {\alpha (\sigma_1^{k+1}-\sigma_1^{k})} + {\sigma_1^{k}} 
     + {P_R^k(1-k_t)} &=& {P(1+k_t)\over \Delta^k} D_D^k, \\
   {\alpha (\sigma_2^{k+1}-\sigma_2^{k})} + {\sigma_2^{k}} &=& {P(1+k_t)\over
     e^2\Delta^k} D_T^k,\\
   {\alpha (\sigma_{12}^{k+1}-\sigma_{12}^{k})} + {\sigma_{12}^{k}} &=&
     {P(1+k_t)\over 2e^2\Delta^k}D_S^k,\end{aligned}
   
where as opposed to the classic EVP, the second term in each equation is at iteration :math:`k` :cite:`Bouillon13`. Also, as opposed to the classic EVP, 
:math:`\Delta t_e` times the number of subcycles (or iterations) does not need to be equal to the advective time step :math:`\Delta t`. 
A last difference between the classic EVP and the revised approach is that the latter one initializes the stresses to 0 at the beginning of each time step, 
while the classic EVP approach uses the previous time level value. The revised EVP is activated by setting the namelist parameter `revised\_evp` = true. 
In the code :math:`\alpha = arlx` and :math:`\beta = brlx`. The values of :math:`arlx` and :math:`brlx` can be set in the namelist. 
It is recommended to use large values of these parameters and to set :math:`arlx=brlx` :cite:`Kimmritz15`.
