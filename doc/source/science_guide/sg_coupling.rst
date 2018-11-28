:tocdepth: 3

.. _coupl:

Coupling With Other Climate Model Components
============================================

The sea ice model exchanges information with the other model components
via a flux coupler. CICE has been coupled into numerous climate models
with a variety of coupling techniques. This document is oriented
primarily toward the CESM Flux Coupler :cite:`Kauffman02`
from NCAR, the first major climate model to incorporate CICE. The flux
coupler was originally intended to gather state variables from the
component models, compute fluxes at the model interfaces, and return
these fluxes to the component models for use in the next integration
period, maintaining conservation of momentum, heat, and fresh water.
However, several of these fluxes are now computed in the ice model
itself and provided to the flux coupler for distribution to the other
components, for two reasons. First, some of the fluxes depend strongly
on the state of the ice, and vice versa, implying that an implicit,
simultaneous determination of the ice state and the surface fluxes is
necessary for consistency and stability. Second, given the various ice
types in a single grid cell, it is more efficient for the ice model to
determine the net ice characteristics of the grid cell and provide the
resulting fluxes, rather than passing several values of the state
variables for each cell. These considerations are explained in more
detail below.

The fluxes and state variables passed between the sea ice model and the
CESM flux coupler are listed in the `Icepack documentation  <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_.
By convention,
directional fluxes are positive downward. In CESM, the sea ice model may
exchange coupling fluxes using a different grid than the computational
grid. This functionality is activated using the namelist variable
``gridcpl_file``. Another namelist variable ``highfreq``, allows the
high-frequency coupling procedure implemented in the Regional Arctic
System Model (RASM). In particular, the relative atmosphere-ice velocity
(:math:`\vec{U}_a-\vec{u}`) is used instead of the full atmospheric
velocity for computing turbulent fluxes in the atmospheric boundary
layer.

The ice fraction :math:`a_i` (aice) is the total fractional ice
coverage of a grid cell. That is, in each cell,

.. math::
   \begin{array}{cl}
                  a_{i}=0 & \mbox{if there is no ice} \\ 
                  a_{i}=1 & \mbox{if there is no open water} \\ 
                  0<a_{i}<1 & \mbox{if there is both ice and open water,}
   \end{array}

where :math:`a_{i}` is the sum of fractional ice areas for each category
of ice. The ice fraction is used by the flux coupler to merge fluxes
from the ice model with fluxes from the other components. For example,
the penetrating shortwave radiation flux, weighted by :math:`a_i`, is
combined with the net shortwave radiation flux through ice-free leads,
weighted by (:math:`1-a_i`), to obtain the net shortwave flux into the
ocean over the entire grid cell. The flux coupler requires the fluxes to
be divided by the total ice area so that the ice and land models are
treated identically (land also may occupy less than 100% of an
atmospheric grid cell). These fluxes are “per unit ice area" rather than
“per unit grid cell area."

For CICE run in stand-alone mode (i.e., uncoupled), the AOMIP shortwave
and longwave radiation formulas are available in **ice\_forcing.F90**.
In function *longwave\_rosati\_miyakoda*, downwelling longwave is
computed as

.. math:: 
   F_{lw\downarrow} = \epsilon\sigma T_s^4 - \epsilon\sigma T_a^4(0.39-0.05e_a^{1/2})(1-0.8f_{cld}) - 4\epsilon\sigma T_a^3(T_s-T_a)
   :label: lwflux

where the atmospheric vapor pressure (mb) is
:math:`e_a = 1000 Q_a/(0.622+0.378Q_a)`, :math:`\epsilon=0.97` is the
ocean emissivity, :math:`\sigma` is the Stephan-Boltzman constant,
:math:`f_{cld}` is the cloud cover fraction, and :math:`T_a` is the
surface air temperature (K). The first term on the right is upwelling
longwave due to the mean (merged) ice and ocean surface temperature,
:math:`T_s` (K), and the other terms on the right represent the net
longwave radiation patterned after :cite:`Rosati88`. 

The downwelling longwave formula of :cite:`Parkinson79` is also
available in function *longwave\_parkinson\_washington*:

.. math:: 
   F_{lw\downarrow} = \epsilon\sigma T_a^4 (1-0.261 \exp\left(-7.77\times 10^{-4}T_a^2\right)\left(1 + 0.275f_{cld}\right)
   :label: lwflux2

The value of :math:`F_{lw\uparrow}` is different for each ice thickness
category, while :math:`F_{lw\downarrow}` depends on the mean value of
surface temperature averaged over all of the thickness categories and
open water. The merged ice-ocean temperature in this formula creates a 
feedback between longwave radiation and sea surface temperature which is
unrealistic, resulting in erroneous model sensitivities to radiative changes, 
e.g. other emissivity values, when run in the stand-alone mode. Although our
stand-alone model test configurations are useful for model development 
purposes, we strongly recommend that scientific conclusions be drawn using 
the model only when coupled with other earth system components.

The AOMIP shortwave forcing formula (in subroutine *compute\_shortwave*)
incorporates the cloud fraction and humidity through the atmospheric
vapor pressure:

.. math:: 
   F_{sw\downarrow} = {\frac{1353 \cos^2 Z}{10^{-3}(\cos Z+2.7)e_a + 1.085\cos Z + 0.1}}\left(1-0.6 f_{cld}^3\right) > 0
   :label: swflux

where :math:`\cos Z` is the cosine of the solar zenith angle.

Many ice models compute the sea surface slope :math:`\nabla H_\circ`
from geostrophic ocean currents provided by an ocean model or other data
source. In our case, the sea surface height :math:`H_\circ` is a
prognostic variable in POP—the flux coupler can provide the surface
slope directly, rather than inferring it from the currents. (The option
of computing it from the currents is provided in subroutine
*dyn\_prep2*.) The sea ice model uses the surface layer currents
:math:`\vec{U}_w` to determine the stress between the ocean and the ice,
and subsequently the ice velocity :math:`\vec{u}`. This stress, relative
to the ice,

.. math::
   \begin{aligned}
   \vec{\tau}_w&=&c_w\rho_w\left|{\vec{U}_w-\vec{u}}\right|\left[\left(\vec{U}_w-\vec{u}\right)\cos\theta
   +\hat{k}\times\left(\vec{U}_w-\vec{u}\right)\sin\theta\right] \end{aligned}
   :label: tauw

is then passed to the flux coupler (relative to the ocean) for use by
the ocean model. Here, :math:`\theta` is the turning angle between
geostrophic and surface currents, :math:`c_w` is the ocean drag
coefficient, :math:`\rho_w` is the density of seawater, and
:math:`\hat{k}` is the vertical unit vector. The turning angle is
necessary if the top ocean model layers are not able to resolve the
Ekman spiral in the boundary layer. If the top layer is sufficiently
thin compared to the typical depth of the Ekman spiral, then
:math:`\theta=0` is a good approximation. Here we assume that the top
layer is thin enough.

Please see the `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_ for additional information about 
atmospheric and oceanic forcing and other data exchanged between the 
flux coupler and the sea ice model.
