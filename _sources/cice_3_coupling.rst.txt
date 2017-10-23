.. _coupl:

Coupling with other climate model components
============================================

The sea ice model exchanges information with the other model components
via a flux coupler. CICE has been coupled into numerous climate models
with a variety of coupling techniques. This document is oriented
primarily toward the CESM Flux Coupler :cite:`KL02`
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
CESM flux coupler are listed in :ref:`tab-flux-cpl`. By convention,
directional fluxes are positive downward. In CESM, the sea ice model may
exchange coupling fluxes using a different grid than the computational
grid. This functionality is activated using the namelist variable
``gridcpl_file``. Another namelist variable ``highfreq``, allows the
high-frequency coupling procedure implemented in the Regional Arctic
System Model (RASM). In particular, the relative atmosphere-ice velocity
(:math:`\vec{U}_a-\vec{u}`) is used instead of the full atmospheric
velocity for computing turbulent fluxes in the atmospheric boundary
layer.

:ref:`tab-flux-cpl`: *Data exchanged between the CESM flux coupler and the sea ice model*

.. _tab-flux-cpl:

.. table:: Table 1

   ===========================   ======================================   =======================================================================================
   Variable                       Description                              Interaction with flux coupler 
   ===========================   ======================================   =======================================================================================
   :math:`z_o`                    Atmosphere level height                  From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`\vec{U}_a`              Wind velocity                            From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`Q_a`                    Specific humidity                        From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`\rho_a`                 Air density                              From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`\Theta_a`               Air potential temperature                From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`T_a`                    Air temperature                          From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`F_{sw\downarrow}`       Incoming shortwave radiation             From *atmosphere model* via flux coupler **to** *sea ice model*
                                  (4 bands)

   :math:`F_{L\downarrow}`        Incoming longwave radiation              From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`F_{rain}`               Rainfall rate                            From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`F_{snow}`               Snowfall rate                            From *atmosphere model* via flux coupler **to** *sea ice model*

   :math:`F_{frzmlt}`             Freezing/melting potential               From *ocean model* via flux coupler **to** *sea ice model*

   :math:`T_w`                    Sea surface temperature                  From *ocean model* via flux coupler **to** *sea ice model*

   :math:`S`                      Sea surface salinity                     From *ocean model* via flux coupler **to** *sea ice model*

   :math:`\nabla H_o`             Sea surface slope                        From *ocean model* via flux coupler **to** *sea ice model*

   :math:`\vec{U}_w`              Surface ocean currents                   From *ocean model* via flux coupler **to** *sea ice model*

   :math:`\vec{\tau}_a`           Wind stress                              From *sea ice model* via flux coupler **to** *atmosphere model*

   :math:`F_s`                    Sensible heat flux                       From *sea ice model* via flux coupler **to** *atmosphere model*
 
   :math:`F_l`                    Latent heat flux                         From *sea ice model* via flux coupler **to** *atmosphere model*

   :math:`F_{L\uparrow}`          Outgoing longwave radiation              From *sea ice model* via flux coupler **to** *atmosphere model*

   :math:`F_{evap}`               Evaporated water                         From *sea ice model* via flux coupler **to** *atmosphere model*

   :math:`\alpha`                 Surface albedo (4 bands)                 From *sea ice model* via flux coupler **to** *atmosphere model*

   :math:`T_{sfc}`                Surface temperature                      From *sea ice model* via flux coupler **to** *atmosphere model*

   :math:`F_{sw\Downarrow}`       Penetrating shortwave radiation          From *sea ice model* via flux coupler **to** *ocean model*

   :math:`F_{water}`              Fresh water flux                         From *sea ice model* via flux coupler **to** *ocean model*

   :math:`F_{hocn}`               Net heat flux to ocean                   From *sea ice model* via flux coupler **to** *ocean model*

   :math:`F_{salt}`               Salt flux                                From *sea ice model* via flux coupler **to** *ocean model*

   :math:`\vec{\tau}_w`           Ice-ocean stress                         From *sea ice model* via flux coupler **to** *ocean model*

   :math:`F_{bio}`                Biogeochemical fluxes                    From *sea ice model* via flux coupler **to** *ocean model*

   :math:`a_{i}`                  Ice fraction                             From *sea ice model* via flux coupler **to** both *ocean and atmosphere models*

   :math:`T^{ref}_{a}`            2m reference temperature (diagnostic)    From *sea ice model* via flux coupler **to** both *ocean and atmosphere models*

   :math:`Q^{ref}_{a}`            2m reference humidity (diagnostic)       From *sea ice model* via flux coupler **to** both *ocean and atmosphere models*

   :math:`F_{swabs}`              Absorbed shortwave (diagnostic)          From *sea ice model* via flux coupler **to** both *ocean and atmosphere models*
   ===========================   ======================================   =======================================================================================

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

In some coupled climate models (for example, recent versions of the U.K.
Hadley Centre model) the surface air temperature and fluxes are computed
within the atmosphere model and are passed to CICE. In this case the
logical parameter ``calc_Tsfc`` in *ice_therm_vertical* is set to false.
The fields ``fsurfn`` (the net surface heat flux from the atmosphere), ``flatn``
(the surface latent heat flux), and ``fcondtopn`` (the conductive flux at
the top surface) for each ice thickness category are copied or derived
from the input coupler fluxes and are passed to the thermodynamic driver
subroutine, *thermo_vertical*. At the end of the time step, the surface
temperature and effective conductivity (i.e., thermal conductivity
divided by thickness) of the top ice/snow layer in each category are
returned to the atmosphere model via the coupler. Since the ice surface
temperature is treated explicitly, the effective conductivity may need
to be limited to ensure stability. As a result, accuracy may be
significantly reduced, especially for thin ice or snow layers. A more
stable and accurate procedure would be to compute the temperature
profiles for both the atmosphere and ice, together with the surface
fluxes, in a single implicit calculation. This was judged impractical,
however, given that the atmosphere and sea ice models generally exist on
different grids and/or processor sets.

.. _atmo:

Atmosphere
----------

The wind velocity, specific humidity, air density and potential
temperature at the given level height :math:`z_\circ` are used to
compute transfer coefficients used in formulas for the surface wind
stress and turbulent heat fluxes :math:`\vec\tau_a`, :math:`F_s`, and
:math:`F_l`, as described below. Wind stress is arguably the primary
forcing mechanism for the ice motion, although the ice–ocean stress,
Coriolis force, and slope of the ocean surface are also important
:cite:`SZRS97`. The sensible and latent heat fluxes,
:math:`F_s` and :math:`F_l`, along with shortwave and longwave
radiation, :math:`F_{sw\downarrow}`, :math:`F_{L\downarrow}`
and :math:`F_{L\uparrow}`, are included in the flux balance that
determines the ice or snow surface temperature when calc\_Tsfc = true.
As described in Section :ref:`thermo`, these fluxes depend nonlinearly
on the ice surface temperature :math:`T_{sfc}`. The balance
equation is iterated until convergence, and the resulting fluxes and
:math:`T_{sfc}` are then passed to the flux coupler.

The snowfall precipitation rate (provided as liquid water equivalent and
converted by the ice model to snow depth) also contributes to the heat
and water mass budgets of the ice layer. Melt ponds generally form on
the ice surface in the Arctic and refreeze later in the fall, reducing
the total amount of fresh water that reaches the ocean and altering the
heat budget of the ice; this version includes two new melt pond
parameterizations. Rain and all melted snow end up in the ocean.

Wind stress and transfer coefficients for the
turbulent heat fluxes are computed in subroutine
*atmo\_boundary\_layer* following :cite:`KL02`. For
clarity, the equations are reproduced here in the present notation.

The wind stress and turbulent heat flux calculation accounts for both
stable and unstable atmosphere–ice boundary layers. Define the
“stability”

.. math::
   \Upsilon = {\kappa g z_\circ\over u^{*2}}
   \left({\Theta^*\over\Theta_a\left(1+0.606Q_a\right)}  +
   {Q^*\over 1/0.606 + Q_a}\right),

where :math:`\kappa` is the von Karman constant, :math:`g` is
gravitational acceleration, and :math:`u^*`, :math:`\Theta^*` and
:math:`Q^*` are turbulent scales for velocity, temperature, and humidity,
respectively:

.. math::
   \begin{aligned}
   u^*&=&c_u \left|\vec{U}_a\right| \\
   \Theta^*&=& c_\theta\left(\Theta_a-T_{sfc}\right) \\
   Q^*&=&c_q\left(Q_a-Q_{sfc}\right).\end{aligned}
   :label: stars

The wind speed has a minimum value of 1 m/s. We have ignored ice motion
in :math:`u^*`, and :math:`T_{sfc}` and
:math:`Q_{sfc}` are the surface temperature and specific
humidity, respectively. The latter is calculated by assuming a saturated
surface, as described in Section :ref:`sfc-forcing`.

Neglecting form drag,the exchange coefficients :math:`c_u`,
:math:`c_\theta` and :math:`c_q` are initialized as

.. math:: 
   \kappa\over \ln(z_{ref}/z_{ice})

and updated during a short iteration, as they depend upon the turbulent
scales. The number of iterations is set by the namelist variable
`natmiter`. (For the case with form drag, see section :ref:`formdrag`.)
Here, :math:`z_{ref}` is a reference height of 10m and
:math:`z_{ice}` is the roughness length scale for the given
sea ice category. :math:`\Upsilon` is constrained to have magnitude less
than 10. Further, defining
:math:`\chi = \left(1-16\Upsilon\right)^{0.25}` and :math:`\chi \geq 1`,
the “integrated flux profiles” for momentum and stability in the
unstable (:math:`\Upsilon <0`) case are given by

.. math::
   \begin{aligned}
   \psi_m = &\mbox{}&2\ln\left[0.5(1+\chi)\right] +
            \ln\left[0.5(1+\chi^2)\right] -2\tan^{-1}\chi +
            {\pi\over 2}, \\
   \psi_s = &\mbox{}&2\ln\left[0.5(1+\chi^2)\right].\end{aligned}

In a departure from the parameterization used in
:cite:`KL02`, we use profiles for the stable case
following :cite:`JAM99`,

.. math::
   \psi_m = \psi_s = -\left[0.7\Upsilon + 0.75\left(\Upsilon-14.3\right)
            \exp\left(-0.35\Upsilon\right) + 10.7\right].

The coefficients are then updated as

.. math::
   \begin{aligned}
   c_u^\prime&=&{c_u\over 1+c_u\left(\lambda-\psi_m\right)/\kappa} \\
   c_\theta^\prime&=& {c_\theta\over 1+c_\theta\left(\lambda-\psi_s\right)/\kappa}\\
   c_q^\prime&=&c_\theta^\prime\end{aligned}

where :math:`\lambda = \ln\left(z_\circ/z_{ref}\right)`. The
first iteration ends with new turbulent scales from
equations :eq:`stars`. After five iterations the latent and sensible
heat flux coefficients are computed, along with the wind stress:

.. math::
   \begin{aligned}
   \nonumber
   C_l&=&\rho_a \left(L_{vap}+L_{ice}\right) u^* c_q \\
   C_s&=&\rho_a c_p u^* c_\theta^* + 1, \\
   \vec{\tau}_a&=&{\rho_a u^{*2}\vec{U}_a\over |\vec{U}_a|},\end{aligned}

where :math:`L_{vap}` and :math:`L_{ice}` are
latent heats of vaporization and fusion, :math:`\rho_a` is the density
of air and :math:`c_p` is its specific heat. Again following
:cite:`JAM99`, we have added a constant to the sensible
heat flux coefficient in order to allow some heat to pass between the
atmosphere and the ice surface in stable, calm conditions.

The atmospheric reference temperature :math:`T_a^{ref}` is computed from
:math:`T_a` and :math:`T_{sfc}` using the coefficients
:math:`c_u`, :math:`c_\theta` and :math:`c_q`. Although the sea ice
model does not use this quantity, it is convenient for the ice model to
perform this calculation. The atmospheric reference temperature is
returned to the flux coupler as a climate diagnostic. The same is true
for the reference humidity, :math:`Q_a^{ref}`.

Additional details about the latent and sensible heat fluxes and other
quantities referred to here can be found in
Section :ref:`sfc-forcing`.

For CICE run in stand-alone mode (i.e., uncoupled), the AOMIP shortwave
and longwave radiation formulas are available in **ice\_forcing.F90**.
In function *longwave\_rosati\_miyakoda*, downwelling longwave is
computed as

.. math:: 
   F_{lw\downarrow} = \epsilon\sigma T_s^4 - \epsilon\sigma T_a^4(0.39-0.05e_a^{1/2})(1-0.8f_{cld}) - 4\epsilon\sigma T_a^3(T_s-T_a)

where the atmospheric vapor pressure (mb) is
:math:`e_a = 1000 Q_a/(0.622+0.378Q_a)`, :math:`\epsilon=0.97` is the
ocean emissivity, :math:`\sigma` is the Stephan-Boltzman constant,
:math:`f_{cld}` is the cloud cover fraction, and :math:`T_a` is the
surface air temperature (K). The first term on the right is upwelling
longwave due to the mean (merged) ice and ocean surface temperature,
:math:`T_s` (K), and the other terms on the right represent the net
longwave radiation patterned after :cite:`RM88`. The
downwelling longwave formula of :cite:`PW79` is also
available in function *longwave\_parkinson\_washington*:

.. math:: 
   F_{lw\downarrow} = \epsilon\sigma T_a^4 (1-0.261 \exp\left(-7.77\times 10^{-4}T_a^2\right)\left(1 + 0.275f_{cld}\right)

The value of :math:`F_{lw\uparrow}` is different for each ice thickness
category, while :math:`F_{lw\downarrow}` depends on the mean value of
surface temperature averaged over all of the thickness categories and
open water.

The AOMIP shortwave forcing formula (in subroutine *compute\_shortwave*)
incorporates the cloud fraction and humidity through the atmospheric
vapor pressure:

.. math:: 
   F_{sw\downarrow} = {1353 \cos^2 Z \over {10^{-3}(\cos Z+2.7)e_a + 1.085\cos Z + 0.1}}\left(1-0.6 f_{cld}^3\right) > 0

where :math:`\cos Z` is the cosine of the solar zenith angle.

.. _ocean:

Ocean
-----

New sea ice forms when the ocean temperature drops below its freezing
temperature. In the Bitz and Lipscomb thermodynamics,
:cite:`BL99` :math:`T_f=-\mu S`, where :math:`S` is the
seawater salinity and :math:`\mu=0.054 \ ^\circ`/ppt is the ratio of the
freezing temperature of brine to its salinity (linear liquidus
approximation). For the mushy thermodynamics, :math:`T_f` is given by a
piecewise linear liquidus relation. The ocean model calculates the new
ice formation; if the freezing/melting potential
:math:`F_{frzmlt}` is positive, its value represents a certain
amount of frazil ice that has formed in one or more layers of the ocean
and floated to the surface. (The ocean model assumes that the amount of
new ice implied by the freezing potential actually forms.)

If :math:`F_{frzmlt}` is negative, it is used to heat already
existing ice from below. In particular, the sea surface temperature and
salinity are used to compute an oceanic heat flux :math:`F_w`
(:math:`\left|F_w\right| \leq \left|F_{frzmlt}\right|`) which
is applied at the bottom of the ice. The portion of the melting
potential actually used to melt ice is returned to the coupler in
:math:`F_{hocn}`. The ocean model adjusts its own heat budget
with this quantity, assuming that the rest of the flux remained in the
ocean.

In addition to runoff from rain and melted snow, the fresh water flux
:math:`F_{water}` includes ice melt water from the top surface
and water frozen (a negative flux) or melted at the bottom surface of
the ice. This flux is computed as the net change of fresh water in the
ice and snow volume over the coupling time step, excluding frazil ice
formation and newly accumulated snow. Setting the namelist option
update\_ocn\_f to true causes frazil ice to be included in the fresh
water and salt fluxes.

There is a flux of salt into the ocean under melting conditions, and a
(negative) flux when sea water is freezing. However, melting sea ice
ultimately freshens the top ocean layer, since the ocean is much more
saline than the ice. The ice model passes the net flux of salt
:math:`F_{salt}` to the flux coupler, based on the net change
in salt for ice in all categories. In the present configuration,
ice\_ref\_salinity is used for computing the salt flux, although the ice
salinity used in the thermodynamic calculation has differing values in
the ice layers.

A fraction of the incoming shortwave :math:`F_{sw\Downarrow}`
penetrates the snow and ice layers and passes into the ocean, as
described in Section :ref:`sfc-forcing`.

Many ice models compute the sea surface slope :math:`\nabla H_\circ`
from geostrophic ocean currents provided by an ocean model or other data
source. In our case, the sea surface height :math:`H_\circ` is a
prognostic variable in POP—the flux coupler can provide the surface
slope directly, rather than inferring it from the currents. (The option
of computing it from the currents is provided in subroutine
*evp\_prep*.) The sea ice model uses the surface layer currents
:math:`\vec{U}_w` to determine the stress between the ocean and the ice,
and subsequently the ice velocity :math:`\vec{u}`. This stress, relative
to the ice,

.. math::
   \begin{aligned}
   \vec{\tau}_w&=&c_w\rho_w\left|{\vec{U}_w-\vec{u}}\right|\left[\left(\vec{U}_w-\vec{u}\right)\cos\theta
   +\hat{k}\times\left(\vec{U}_w-\vec{u}\right)\sin\theta\right] \end{aligned}

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

For CICE run in stand-alone mode (i.e., uncoupled), a thermodynamic slab
ocean mixed-layer parameterization is available in **ice\_ocean.F90**.
The turbulent fluxes are computed above the water surface using the same
parameterizations as for sea ice, but with parameters appropriate for
the ocean. The surface flux balance takes into account the turbulent
fluxes, oceanic heat fluxes from below the mixed layer, and shortwave
and longwave radiation, including that passing through the sea ice into
the ocean. If the resulting sea surface temperature falls below the
salinity-dependent freezing point, then new ice (frazil) forms.
Otherwise, heat is made available for melting the ice.

.. _formdrag:

Variable exchange coefficients
------------------------------

In the default CICE setup, atmospheric and oceanic neutral drag
coefficients (:math:`c_u` and :math:`c_w`) are assumed constant in time
and space. These constants are chosen to reflect friction associated
with an effective sea ice surface roughness at the ice–atmosphere and
ice–ocean interfaces. Sea ice (in both Arctic and Antarctic) contains
pressure ridges as well as floe and melt pond edges that act as discrete
obstructions to the flow of air or water past the ice, and are a source
of form drag. Following :cite:`TFSFFKLB14` and based on
recent theoretical developments :cite:`LGHA12,LLCL11`, the
neutral drag coefficients can now be estimated from properties of the
ice cover such as ice concentration, vertical extent and area of the
ridges, freeboard and floe draft, and size of floes and melt ponds. The
new parameterization allows the drag coefficients to be coupled to the
sea ice state and therefore to evolve spatially and temporally. This
parameterization is contained in the subroutine *neutral\_drag\_coeffs*
and is accessed by setting `formdrag` = true in the namelist.

Following :cite:`TFSFFKLB14`, consider the general case of
fluid flow obstructed by N randomly oriented obstacles of height
:math:`H` and transverse length :math:`L_y`, distributed on a domain
surface area :math:`S_T`. Under the assumption of a logarithmic fluid
velocity profile, the general formulation of the form drag coefficient
can be expressed as

.. math:: 
   C_d=\frac{N c S_c^2 \gamma L_y  H}{2 S_T}\left[\frac{\ln(H/z_0)}{\ln(z_{ref}/z_0)}\right]^2,
   :label: formdrag

where :math:`z_0` is a roughness length parameter at the top or bottom
surface of the ice, :math:`\gamma` is a geometric factor, :math:`c` is
the resistance coefficient of a single obstacle, and :math:`S_c` is a
sheltering function that takes into account the shielding effect of the
obstacle,

.. math:: 
   S_{c}=\left(1-\exp(-s_l D/H)\right)^{1/2},
   :label: shelter

with :math:`D` the distance between two obstacles and :math:`s_l` an
attenuation parameter.

As in the original drag formulation in CICE (sections :ref:`atmo` and
:ref:`ocean`), :math:`c_u` and :math:`c_w` along with the transfer
coefficients for sensible heat, :math:`c_{\theta}`, and latent heat,
:math:`c_{q}`, are initialized to a situation corresponding to neutral
atmosphere–ice and ocean–ice boundary layers. The corresponding neutral
exchange coefficients are then replaced by coefficients that explicitly
account for form drag, expressed in terms of various contributions as

.. math::
   \tt{Cdn\_atm}  = \tt{Cdn\_atm\_rdg} + \tt{Cdn\_atm\_floe} + \tt{Cdn\_atm\_skin} + \tt{Cdn\_atm\_pond} ,
   :label: Cda

.. math::
   \tt{Cdn\_ocn}  =  \tt{Cdn\_ocn\_rdg} + \tt{Cdn\_ocn\_floe} + \tt{Cdn\_ocn\_skin}. 
   :label: Cdw

The contributions to form drag from ridges (and keels underneath the
ice), floe edges and melt pond edges can be expressed using the general
formulation of equation :eq:`formdrag` (see :cite:`TFSFFKLB14` for
details). Individual terms in equation :eq:`Cdw` are fully described in
:cite:`TFSFFKLB14`. Following :cite:`Arya75`
the skin drag coefficient is parametrized as

.. math:: 
   { \tt{Cdn\_(atm/ocn)\_skin}}=a_{i} \left(1-m_{(s/k)} \frac{H_{(s/k)}}{D_{(s/k)}}\right)c_{s(s/k)}, \mbox{       if  $\displaystyle\frac{H_{(s/k)}}{D_{(s/k)}}\ge\frac{1}{m_{(s/k)}}$,}
   :label: skindrag

where :math:`m_s` (:math:`m_k`) is a sheltering parameter that depends
on the average sail (keel) height, :math:`H_s` (:math:`H_k`), but is
often assumed constant, :math:`D_s` (:math:`D_k`) is the average
distance between sails (keels), and :math:`c_{ss}` (:math:`c_{sk}`) is
the unobstructed atmospheric (oceanic) skin drag that would be attained
in the absence of sails (keels) and with complete ice coverage,
:math:`a_{ice}=1`.

Calculation of equations :eq:`formdrag` – :eq:`skindrag` requires that small-scale geometrical
properties of the ice cover be related to average grid cell quantities
already computed in the sea ice model. These intermediate quantities are
briefly presented here and described in more detail in
:cite:`TFSFFKLB14`. The sail height is given by

.. math:: 
   H_{s} = \displaystyle 2\frac{v_{rdg}}{a_{rdg}}\left(\frac{\alpha\tan \alpha_{k} R_d+\beta \tan \alpha_{s} R_h}{\phi_r\tan \alpha_{k} R_d+\phi_k \tan \alpha_{s} R_h^2}\right),
   :label: Hs

and the distance between sails\ 

.. math:: 
   D_{s} = \displaystyle 2 H_s\frac{a_{i}}{a_{rdg}} \left(\frac{\alpha}{\tan \alpha_s}+\frac{\beta}{\tan \alpha_k}\frac{R_h}{R_d}\right),
   :label: Ds

where :math:`0<\alpha<1` and :math:`0<\beta<1` are weight functions,
:math:`\alpha_{s}` and :math:`\alpha_{k}` are the sail and keel slope,
:math:`\phi_s` and :math:`\phi_k` are constant porosities for the sails
and keels, and we assume constant ratios for the average keel depth and
sail height (:math:`H_k/H_s=R_h`) and for the average distances between
keels and between sails (:math:`D_k/D_s=R_d`). With the assumption of
hydrostatic equilibrium, the effective ice plus snow freeboard is
:math:`H_{f}=\bar{h_i}(1-\rho_i/\rho_w)+\bar{h_s}(1-\rho_s/\rho_w)`,
where :math:`\rho_i`, :math:`\rho_w` and :math:`\rho_s` are
respectively the densities of sea ice, water and snow, :math:`\bar{h_i}`
is the mean ice thickness and :math:`\bar{h_s}` is the mean snow
thickness (means taken over the ice covered regions). For the melt pond
edge elevation we assume that the melt pond surface is at the same level
as the ocean surface surrounding the floes
:cite:`FF07,FFT10,FSFH12` and use the simplification
:math:`H_p = H_f`. Finally to estimate the typical floe size
:math:`L_A`, distance between floes, :math:`D_F`, and melt pond size,
:math:`L_P` we use the parameterizations of :cite:`LGHA12`
to relate these quantities to the ice and pond concentrations. All of
these intermediate quantities are available as history output, along
with `Cdn\_atm`, `Cdn\_ocn` and the ratio `Cdn\_atm\_ratio\_n` between the
total atmospheric drag and the atmospheric neutral drag coefficient.

We assume that the total neutral drag coefficients are thickness
category independent, but through their dependance on the diagnostic
variables described above, they vary both spatially and temporally. The
total drag coefficients and heat transfer coefficients will also depend
on the type of stratification of the atmosphere and the ocean, and we
use the parameterization described in section :ref:`atmo` that accounts
for both stable and unstable atmosphere–ice boundary layers. In contrast
to the neutral drag coefficients the stability effect of the atmospheric
boundary layer is calculated separately for each ice thickness category.

The transfer coefficient for oceanic heat flux to the bottom of the ice
may be varied based on form drag considerations by setting the namelist
variable `fbot\_xfer\_type` to `Cdn\_ocn`; this is recommended when using
the form drag parameterization. Its default value of the transfer
coefficient is 0.006 (`fbot\_xfer\_type = ’constant’`).
