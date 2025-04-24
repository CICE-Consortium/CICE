:tocdepth: 3

.. _dynam:

Dynamics
========

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

For clarity, the two components of Equation :eq:`vpmom` are

.. math::
   \begin{aligned}
   m{\partial u\over\partial t} &= {\partial\sigma_{1j}\over\partial x_j} + \tau_{ax} +
     a_i c_w \rho_w
     \left|{\bf U}_w - {\bf u}\right| \left[\left(U_w-u\right)\cos\theta \mp \left(V_w-v\right)\sin\theta\right]
     -C_bu +mfv - mg{\partial H_\circ\over\partial x}, \\
   m{\partial v\over\partial t} &= {\partial\sigma_{2j}\over\partial x_j} + \tau_{ay} +
     a_i c_w \rho_w
     \left|{\bf U}_w - {\bf u}\right| \left[ \pm \left(U_w-u\right)\sin\theta + \left(V_w-v\right)\cos\theta\right]
     -C_bv-mfu - mg{\partial H_\circ\over\partial y}. \end{aligned}
   :label: momsys

On the B grid, the equations above are solved at the U point for the collocated u and v components (see figure :ref:`fig-Bgrid`). On the C grid, however, the two components are not collocated: the u component is at the E point while the v component is at the N point.

The B grid spatial discretization is based on a variational method described in :cite:`Hunke97` and :cite:`Hunke02`. A bilinear discretization is used for the stress terms
:math:`\partial\sigma_{ij}/\partial x_j`,
which enables the discrete equations to be derived from the
continuous equations written in curvilinear coordinates. In this
manner, metric terms associated with the curvature of the grid are
incorporated into the discretization explicitly. Details pertaining to
the spatial discretization are found in :cite:`Hunke02`

On the C grid, however, a finite difference approach is used for the spatial discretization. The C grid discretization is based on :cite:`Bouillon09`, :cite:`Bouillon13` and :cite:`Kimmritz16`.

There are different approaches in the CICE code for representing sea ice
rheology and for solving the sea ice momentum equation: the viscous-plastic (VP) rheology :cite:`Hibler79` with an implicit method,
the elastic-viscous-plastic (EVP) :cite:`Hunke97` model which represents a modification of the
VP model, the revised EVP (rEVP) approach :cite:`Lemieux12,Bouillon13` and the elastic-anisotropic-plastic (EAP) model which explicitly accounts for the sub-continuum
anisotropy of the sea ice cover :cite:`Wilchinsky06,Weiss09`. If
``kdyn`` = 1 in the namelist then the EVP model is used (module
**ice\_dyn\_evp.F90**), while ``kdyn`` = 2 is associated with the EAP
model (**ice\_dyn\_eap.F90**), and ``kdyn`` = 3 is associated with the
VP model (**ice\_dyn\_vp.F90**). The rEVP approach can be used by setting ``kdyn`` = 1 and  ``revised_evp`` = true in the namelist.

At times scales associated with the
wind forcing, the EVP model reduces to the VP model while the EAP model
reduces to the anisotropic rheology described in detail in
:cite:`Wilchinsky06,Tsamados13`. At shorter time scales the
adjustment process takes place in both models by a numerically more
efficient elastic wave mechanism. While retaining the essential physics,
this elastic wave modification leads to a fully explicit numerical
scheme which greatly improves the model’s computational efficiency. The rEVP is also a fully explicit scheme which by construction should lead to the VP solution. 

The EVP sea ice dynamics model is thoroughly documented in
:cite:`Hunke97`, :cite:`Hunke01`,
:cite:`Hunke02` and :cite:`Hunke03` and the EAP
dynamics in :cite:`Tsamados13`. Simulation results and
performance of the EVP and EAP models have been compared with the VP
model and with each other in realistic simulations of the Arctic
respectively in :cite:`Hunke99` and
:cite:`Tsamados13`.

The EVP numerical
implementation in this code release is that of :cite:`Hunke02`
and :cite:`Hunke03`, with revisions to the numerical solver as
in :cite:`Bouillon13`. Details about the rEVP solver can be found in  :cite:`Lemieux12`, :cite:`Bouillon13`, :cite:`Kimmritz15` and :cite:`Koldunov19`. The implementation of the EAP sea ice
dynamics into CICE is described in detail in
:cite:`Tsamados13`.

The VP solver implementation mostly follows :cite:`Lemieux08`, with
FGMRES :cite:`Saad93` as the linear solver and GMRES as the preconditioner.
Note that the VP solver has not yet been tested on the ``tx1`` grid.

The EVP, rEVP, EAP and VP approaches are all available with the B grid. However, at the moment, only the EVP and rEVP schemes are possible with the C grid.

The dynamics are solved for all gridcells with area concentration greater than ``dyn_area_min`` and mass
greater than ``dyn_mass_min``.  These parameters are respectively 0.001 and 0.01 by default but can be set in 
namelist.  Lower values can improve the solution but also lead to instabilities.

Here we summarize the equations and
direct the reader to the above references for details.

.. _momentumTS:

**********************
Momentum time stepping
**********************

.. _evp-momentum:

EVP time discretization and solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The momentum equation is discretized in time as follows, for the classic
EVP approach.
In the code,
:math:`{\tt vrel}=a_i c_w \rho_w\left|{\bf U}_w - {\bf u}^k\right|` and
:math:`C_b=T_b \left( \sqrt{(u^k)^2+(v^k)^2}+u_0 \right)^{-1}`,
where :math:`k` denotes the subcycling step. The following equations
illustrate the time discretization and define some of the other
variables used in the code.

.. math::
   \underbrace{\left({m\over\Delta t_e}+{\tt vrel} \cos\theta\ + C_b \right)}_{\tt cca} u^{k+1}
   - \underbrace{\left(mf \pm {\tt vrel}\sin\theta\right)}_{\tt ccb}v^{l}
    =  &\underbrace{{\partial\sigma_{1j}^{k+1}\over\partial x_j}}_{\tt strintx}
    + \underbrace{\tau_{ax} - mg{\partial H_\circ\over\partial x} }_{\tt forcex} \\
     &+ {\tt vrel}\underbrace{\left(U_w\cos\theta \mp V_w\sin\theta\right)}_{\tt waterx}  + {m\over\Delta t_e}u^k,
   :label: umom

.. math::
    \underbrace{\left(mf \pm {\tt vrel}\sin\theta\right)}_{\tt ccb} u^{l}
   + \underbrace{\left({m\over\Delta t_e}+{\tt vrel} \cos\theta + C_b \right)}_{\tt cca}v^{k+1}
    =  &\underbrace{{\partial\sigma_{2j}^{k+1}\over\partial x_j}}_{\tt strinty}
    + \underbrace{\tau_{ay} - mg{\partial H_\circ\over\partial y} }_{\tt forcey} \\
     &+ {\tt vrel}\underbrace{\left( \pm U_w\sin\theta+V_w\cos\theta\right)}_{\tt watery}  + {m\over\Delta t_e}v^k,
   :label: vmom

where :math:`{\tt vrel}\ \cdot\ {\tt waterx(y)}= {\tt taux(y)}` and the definitions of :math:`u^{l}` and :math:`v^{l}` vary depending on the grid.

As :math:`u` and :math:`v` are collocated on the B grid, :math:`u^{l}` and :math:`v^{l}` are respectively :math:`u^{k+1}` and :math:`v^{k+1}` such that this system of equations can be solved as follows. Define

.. math::
   \hat{u} = F_u + \tau_{ax} - mg{\partial H_\circ\over\partial x} + {\tt vrel} \left(U_w\cos\theta \mp V_w\sin\theta\right) + {m\over\Delta t_e}u^k
   :label: cevpuhat

.. math::
   \hat{v} = F_v + \tau_{ay} - mg{\partial H_\circ\over\partial y} + {\tt vrel} \left(\pm U_w\sin\theta + V_w\cos\theta\right) + {m\over\Delta t_e}v^k,
   :label: cevpvhat

where :math:`{\bf F} = \nabla\cdot\sigma^{k+1}`. Then

.. math::
   \begin{aligned}
   \left({m\over\Delta t_e} +{\tt vrel}\cos\theta\ + C_b \right)u^{k+1} - \left(mf \pm {\tt vrel}\sin\theta\right) v^{k+1} &= \hat{u}  \\
   \left(mf \pm {\tt vrel}\sin\theta\right) u^{k+1} + \left({m\over\Delta t_e} +{\tt vrel}\cos\theta + C_b \right)v^{k+1} &= \hat{v}.\end{aligned}

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
   b = mf \pm {\tt vrel}\sin\theta.
   :label: cevpb

Note that the time discretization and solution method for the EAP is exactly the same as for the B grid EVP. More details on the EAP model are given in Section :ref:`stress-eap`.

However, on the C grid, :math:`u` and :math:`v` are not collocated. When solving the :math:`u` momentum equation for :math:`u^{k+1}` (at the E point), :math:`v^{l}=v^{k}_{int}` where :math:`v^{k}_{int}` is :math:`v^{k}` from the surrounding N points interpolated to the E point. The same approach is used for the :math:`v` momentum equation. With this explicit treatment of the off-diagonal terms :cite:`Kimmritz16`, :math:`u^{k+1}` and :math:`v^{k+1}` are obtained by solving

.. math::
   \begin{aligned}
   u^{k+1} = {\hat{u} + b v^{k}_{int} \over a} \\
   v^{k+1} = {\hat{v} - b u^{k}_{int} \over a}. \end{aligned}

.. _revp-momentum:

Revised EVP time discretization and solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The revised EVP approach is based on a pseudo-time iterative scheme :cite:`Lemieux12`, :cite:`Bouillon13`, :cite:`Kimmritz15`. By construction, the revised EVP approach should lead to the VP solution
(given the right numerical parameters and a sufficiently large number of iterations). To do so, the inertial term is formulated such that it matches the backward Euler approach of
implicit solvers and there is an additional term for the pseudo-time iteration. Hence, with the revised approach, the discretized momentum equations :eq:`umom` and :eq:`vmom` become

.. math::
    {\beta^*(u^{k+1}-u^k)\over\Delta t_e} + {m(u^{k+1}-u^n)\over\Delta t} + {\left({\tt vrel} \cos\theta + C_b \right)} u^{k+1}
    - & {\left(mf \pm {\tt vrel}\sin\theta\right)} v^{l}
    =  {{\partial\sigma_{1j}^{k+1}\over\partial x_j}}
    + {\tau_{ax}} \\
      & - {mg{\partial H_\circ\over\partial x} }
    + {\tt vrel} {\left(U_w\cos\theta \mp V_w\sin\theta\right)},
    :label: umomr

.. math::
    {\beta^*(v^{k+1}-v^k)\over\Delta t_e} + {m(v^{k+1}-v^n)\over\Delta t} + {\left({\tt vrel} \cos\theta + C_b \right)}v^{k+1}
    + & {\left(mf \pm {\tt vrel}\sin\theta\right)} u^{l}
    =  {{\partial\sigma_{2j}^{k+1}\over\partial x_j}}
    + {\tau_{ay}} \\
     & - {mg{\partial H_\circ\over\partial y} }
    + {\tt vrel}{\left( \pm U_w\sin\theta+V_w\cos\theta\right)},
    :label: vmomr

where :math:`\beta^*` is a numerical parameter and :math:`u^n, v^n` are the components of the previous time level solution.
With :math:`\beta=\beta^* \Delta t \left(  m \Delta t_e \right)^{-1}` :cite:`Bouillon13`, these equations can be written as

.. math::
   \underbrace{\left((\beta+1){m\over\Delta t}+{\tt vrel} \cos\theta\ + C_b \right)}_{\tt cca} u^{k+1}
   - \underbrace{\left(mf \pm {\tt vrel} \sin\theta\right)}_{\tt ccb} & v^{l}
    = \underbrace{{\partial\sigma_{1j}^{k+1}\over\partial x_j}}_{\tt strintx}
    + \underbrace{\tau_{ax} - mg{\partial H_\circ\over\partial x} }_{\tt forcex} \\
    & + {\tt vrel}\underbrace{\left(U_w\cos\theta \mp V_w\sin\theta\right)}_{\tt waterx}  + {m\over\Delta t}(\beta u^k + u^n),
   :label: umomr2

.. math::
    \underbrace{\left(mf \pm {\tt vrel}\sin\theta\right)}_{\tt ccb} u^{l}
   + \underbrace{\left((\beta+1){m\over\Delta t}+{\tt vrel} \cos\theta + C_b \right)}_{\tt cca} & v^{k+1}
    = \underbrace{{\partial\sigma_{2j}^{k+1}\over\partial x_j}}_{\tt strinty}
    + \underbrace{\tau_{ay} - mg{\partial H_\circ\over\partial y} }_{\tt forcey} \\
    & + {\tt vrel}\underbrace{\left( \pm U_w\sin\theta+V_w\cos\theta\right)}_{\tt watery}  + {m\over\Delta t}(\beta v^k + v^n),
   :label: vmomr2

At this point, the solutions :math:`u^{k+1}` and :math:`v^{k+1}` for the B or the C grids are obtained in the same manner as for the standard EVP approach (see Section :ref:`evp-momentum` for details).

.. _vp-momentum:

Implicit (VP) time discretization and solution
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In the VP approach, equation :eq:`momsys` is discretized implicitly using a Backward Euler approach,
and stresses are not computed explicitly:

.. math::
  \begin{aligned}
  m\frac{(u^{n}-u^{n-1})}{\Delta t} &= \frac{\partial \sigma_{1j}^n}{\partial x_j}
  - \tau_{w,x}^n + \tau_{b,x}^n +  mfv^n
   + r_{x}^n,
  \\
  m\frac{(v^{n}-v^{n-1})}{\Delta t} &= \frac{\partial \sigma^{n} _{2j}}{\partial x_j}
  - \tau_{w,y}^n + \tau_{b,y}^n -  mfu^{n}
   + r_{y}^n
  \end{aligned}
  :label: u_sit

where :math:`r = (r_x,r_y)` contains all terms that do not depend on the velocities :math:`u^n, v^n` (namely the sea surface tilt and the wind stress).
As the water drag, seabed stress and rheology term depend on the velocity field, the only unknowns in equation :eq:`u_sit` are :math:`u^n` and :math:`v^n`.

Once discretized in space, equation :eq:`u_sit` leads to a system of :math:`N` nonlinear equations with :math:`N` unknowns that can be concisely written as

.. math::
  \mathbf{A}(\mathbf{u})\mathbf{u} = \mathbf{b}(\mathbf{u}),
  :label: nonlin_sys

where :math:`\mathbf{A}` is an :math:`N\times N` matrix and :math:`\mathbf{u}` and :math:`\mathbf{b}` are vectors of size :math:`N`.
Note that we have dropped the time level index :math:`n`.
The vector :math:`\mathbf{u}` is formed by stacking first the :math:`u` components, followed by the :math:`v` components of the discretized ice velocity.
The vector :math:`\mathbf{b}` is a function of the velocity vector :math:`\mathbf{u}` because of the water and seabed stress terms as well as parts of the rheology term that depend non-linearly on :math:`\mathbf{u}`.

The nonlinear system :eq:`nonlin_sys` is solved using a Picard iteration method.
Starting from a previous iterate :math:`\mathbf{u}_{k-1}`, the nonlinear system is linearized by substituting :math:`\mathbf{u}_{k-1}` in the expression of the matrix :math:`\mathbf{A}` and the vector :math:`\mathbf{b}`:

.. math::
  \mathbf{A}(\mathbf{u}_{k-1})\mathbf{u}_{k} =  \mathbf{b}(\mathbf{u}_{k-1})
  :label: picard

The resulting linear system is solved using the Flexible Generalized Minimum RESidual (FGMRES, :cite:`Saad93`) method and this process is repeated iteratively.

The maximum number of Picard iterations can be set using the namelist flag ``maxits_nonlin``.
The relative tolerance for the Picard solver can be set using the namelist flag ``reltol_nonlin``.
The Picard iterative process stops when :math:`\left\lVert \mathbf{u}_{k} \right\rVert_2 < {\tt reltol\_nonlin} \cdot \left\lVert\mathbf{u}_{0}\right\rVert_2` or when ``maxits_nonlin`` is reached.

Parameters for the FGMRES linear solver and the preconditioner can be controlled using additional namelist flags (see :ref:`dynamics_nml`).


.. _surfstress:

********************
Surface stress terms
********************

The formulation for the wind stress is described in `Icepack Documentation <https://cice-consortium-icepack.readthedocs.io/en/main/science_guide/index.html>`_. Below, some details about the ice-ocean stress and the seabed stress are given. 

Ice-Ocean stress
~~~~~~~~~~~~~~~~

At the end of each (thermodynamic) time step, the
ice–ocean stress must be constructed from :math:`{\tt taux(y)}` and the terms
containing :math:`{\tt vrel}` on the left hand side of the equations.
The water stress calculation has a hemispheric dependence on the sign of the
:math:`\pm {\tt vrel}\sin\theta` term.

The Hibler-Bryan form for the ice-ocean stress :cite:`Hibler87`
is included in **ice\_dyn\_shared.F90** but is currently commented out,
pending further testing.

.. _seabedstress:

Seabed stress
~~~~~~~~~~~~~

CICE includes two options for calculating the seabed stress,
i.e. the term in the momentum equation that represents the interaction
between grounded ice keels and the seabed. The seabed stress can be
activated by setting ``seabed_stress`` to true in the namelist. The seabed stress (or basal
stress) parameterization of :cite:`Lemieux16` is chosen if ``seabed_stress_method`` = ``LKD`` while the approach based on the probability of contact between the ice and the seabed is used if ``seabed_stress_method`` = ``probabilistic``.

For both parameterizations, the components of the seabed
stress are expressed as :math:`\tau_{bx}=C_bu` and
:math:`\tau_{by}=C_bv`, where :math:`C_b` is a seabed stress
coefficient.

The two parameterizations differ in their calculation of
the :math:`C_b` coefficients.

Note that the user must provide a bathymetry field for using these
grounding schemes. It is suggested to have a bathymetry field with water depths
larger than 5 m that represents well shallow water (less than 30 m) regions such as the Laptev Sea
and the East Siberian Sea.

**Seabed stress based on linear keel draft (LKD)**

This parameterization for the seabed stress is described in
:cite:`Lemieux16`. It assumes that the largest keel draft varies linearly with the mean thickness in a grid cell (i.e. sea ice volume). The :math:`C_b` coefficients are expressed as

.. math::
   C_b= k_2 \max [0,(h - h_{c})]  e^{-\alpha_b * (1 - a)} (\sqrt{u^2+v^2}+u_0)^{-1}, \\
   :label: Cb

where :math:`k_2` determines the maximum seabed stress that can be sustained by the grounded parameterized ridge(s), :math:`u_0`
is a small residual velocity and :math:`\alpha_b` is a parameter to ensure that the seabed stress quickly drops when
the ice concentration is smaller than 1. In the code, :math:`k_2 \max [0,(h - h_{c})]  e^{-\alpha_b * (1 - a)}` is defined as
:math:`T_b`. 

On the B grid, the quantities :math:`h`, :math:`a` and :math:`h_{c}` are calculated at
the U point and are referred to as :math:`h_u`, :math:`a_{u}` and :math:`h_{cu}`. They are respectively given by

.. math::
   h_u=\max[v_i(i,j),v_i(i+1,j),v_i(i,j+1),v_i(i+1,j+1)], \\
   :label: hu

.. math::
   a_u=\max[a_i(i,j),a_i(i+1,j),a_i(i,j+1),a_i(i+1,j+1)], \\
   :label: au

.. math::
   h_{cu}=a_u h_{wu} / k_1, \\
   :label: hcu

where the :math:`a_i` and :math:`v_i` are the total ice concentrations and ice volumes around the U point :math:`i,j` and
:math:`k_1` is a parameter that defines the critical ice thickness :math:`h_{cu}` at which the parameterized
ridge(s) reaches the seafloor for a water depth :math:`h_{wu}=\min[h_w(i,j),h_w(i+1,j),h_w(i,j+1),h_w(i+1,j+1)]`. Given the formulation of :math:`C_b` in equation :eq:`Cb`, the seabed stress components are non-zero only
when :math:`h_u > h_{cu}`.

As :math:`u` and :math:`v` are not collocated on the C grid, :math:`T_b` is calculated at E and N points. For example, at the E point, :math:`h_e`, :math:`a_{e}` and :math:`h_{ce}` are respectively

.. math::
   h_e=\max[v_i(i,j),v_i(i+1,j)], \\
   :label: he

.. math::
   a_e=\max[a_i(i,j),a_i(i+1,j)], \\
   :label: ae

.. math::
   h_{ce}=a_e h_{we} / k_1, \\
   :label: hce

where :math:`h_{we}=\min[h_w(i,j),h_w(i+1,j)]`. Similar calculations are done at the N points. 

To prevent unrealistic grounding, :math:`T_b` is set to zero when :math:`h_{wu}`
is larger than 30 m (same idea on the C grid depending on :math:`h_{we}` and :math:`h_{wn}`). This maximum value is chosen based on observations of large keels in the Arctic Ocean :cite:`Amundrud04`.

The maximum seabed stress depends on the weight of the ridge
above hydrostatic balance and the value of :math:`k_2`. It is, however, the parameter :math:`k_1` that has the most notable impact on the simulated extent of landfast ice.
The value of :math:`k_1` can be changed at runtime using the namelist variable ``k1``.

**Seabed stress based on probabilistic approach**

This more sophisticated grounding parameterization computes the seabed stress based
on the probability of contact between the ice thickness distribution
(ITD) and the seabed :cite:`Dupont22`. Multi-thickness category models such as CICE typically use a
few thickness categories (5-10). This crude representation of the ITD
does not resolve the tail of the ITD, which is crucial for grounding
events.

To represent the tail of the distribution, the simulated ITD is
converted to a positively skewed probability function :math:`f(x)`
with :math:`x` the sea ice thickness. The mean and variance are set
equal to the ones of the original ITD. A
log-normal distribution is used for :math:`f(x)`.

It is assumed that the bathymetry :math:`y` (at the 't' point) follows a normal
distribution :math:`b(y)`. The mean of :math:`b(y)` comes from the user's bathymetry field and the
standard deviation :math:`\sigma_b` is currently fixed to 2.5 m. Two
possible improvements would be to specify a distribution based on high
resolution bathymetry data and to take into account variations of the
water depth due to changes in the sea surface height.

Assuming hydrostatic balance and neglecting the impact of snow, the draft of floating ice of thickness
:math:`x` is :math:`D(x)=\rho_i x / \rho_w` where :math:`\rho_i` is the sea ice density. Hence, the probability of contact (:math:`P_c`) between the
ITD and the seabed is given by

.. math::
   P_c=\int_{0}^{\inf} \int_{0}^{D(x)} g(x)b(y) dy dx \label{prob_contact}.

:math:`T_b` is first calculated at the T point (referred to as :math:`T_{bt}`). :math:`T_{bt}` depends on the weight of the ridge in excess of hydrostatic balance. The parameterization first calculates

.. math::
   T_{bt}^*=\mu_s g \int_{0}^{\inf} \int_{0}^{D(x)} (\rho_i x - \rho_w
   y)g(x)b(y) dy dx, \\
   :label: Tbt

and then obtains :math:`T_{bt}` by multiplying :math:`T_{bt}^*` by :math:`e^{-\alpha_b * (1 - a_i)}` (similar to what is done for ``seabed_stress_method`` = ``LKD``).

To calculate :math:`T_{bt}^*` in equation :eq:`Tbt`, :math:`f(x)` and :math:`b(y)` are discretized using many small categories (100). :math:`f(x)` is discretized between 0 and 50 m while :math:`b(y)` is truncated at plus and minus three :math:`\sigma_b`. :math:`f(x)` is also modified by setting it to	zero after a certain percentile of the log-normal distribution. This percentile, which is currently set to 99.7%, notably affects the simulation of landfast ice and is used as a tuning parameter. Its impact is similar to the one of the parameter :math:`k_1` for the LKD method.

On the B grid, :math:`T_b` at the U point is calculated from the T point values around it according to

.. math::
   T_{bu}=\max[T_{bt}(i,j),T_{bt}(i+1,j),T_{bt}(i,j+1),T_{bt}(i+1,j+1)]. \\
   :label: Tb

Following again the LKD method, the seabed stress coefficients are finally expressed as

.. math::
   C_b= T_{bu} (\sqrt{u^2+v^2}+u_0)^{-1}. \\
   :label: Cb2

On the C grid, :math:`T_b` is needs to be calculated at the E and N points. :math:`T_{be}` and :math:`T_{bn}` are respectively given by

.. math::
   T_{be}=\max[T_{bt}(i,j),T_{bt}(i+1,j)], \\
   :label: Tbe

.. math::
   T_{bn}=\max[T_{bt}(i,j),T_{bt}(i,j+1)]. \\
   :label: Tbn

The :math:`C_{b}` are different at the E and N points and are respectively :math:`T_{be} (\sqrt{u^2+v^2_{int}}+u_0)^{-1}` and :math:`T_{bn} (\sqrt{u^2_{int} + v^2}+u_0)^{-1}` where :math:`v_{int}` (:math:`u_{int}`) is :math:`v` ( :math:`u`) interpolated to the E (N) point.

.. _internal-stress:

********
Rheology
********

For convenience we formulate the stress tensor :math:`\bf \sigma` in
terms of :math:`\sigma_1=\sigma_{11}+\sigma_{22}` (``stressp``),
:math:`\sigma_2=\sigma_{11}-\sigma_{22}` (``stressm``), and introduce the
divergence, :math:`D_D`, and the horizontal tension and shearing
strain rates, :math:`D_T` and :math:`D_S` respectively:

.. math::
   D_D = \dot{\epsilon}_{11} + \dot{\epsilon}_{22},

.. math::
   D_T = \dot{\epsilon}_{11} - \dot{\epsilon}_{22},

.. math::
   D_S = 2\dot{\epsilon}_{12},

where

.. math::
   \dot{\epsilon}_{ij} = {1\over 2}\left({{\partial u_i}\over{\partial x_j}} + {{\partial u_j}\over{\partial x_i}}\right)

Note that :math:`\sigma_1` and :math:`\sigma_2` are not to be confused with the normalized principal stresses,
:math:`\sigma_{n,1}` and :math:`\sigma_{n,2}` (``sig1`` and ``sig2``), which are defined as:

.. math::
   \sigma_{n,1}, \sigma_{n,2} = \frac{1}{P} \left( \frac{\sigma_1}{2} \pm \sqrt{\left(\frac{\sigma_2}{2}\right)^2 + \sigma_{12}^2} \right)

where :math:`P` is the ice strength.

In addition to the normalized principal stresses, CICE can output the internal ice pressure which is an important field to support navigation in ice-infested water.
The internal ice pressure (``sigP``) is the average of the normal stresses (:math:`\sigma_{11}`, :math:`\sigma_{22}`) multiplied by :math:`-1` and
is therefore simply equal to :math:`-\sigma_1/2`.

.. _stress-vp:

Viscous-Plastic
~~~~~~~~~~~~~~~

The VP constitutive law is given by

.. math::
   \sigma_{ij} = 2 \eta \dot{\epsilon}_{ij} + (\zeta - \eta) D_D - P_R\frac{\delta_{ij}}{2}
   :label: vp-const

where :math:`\eta` and :math:`\zeta` are the bulk and shear viscosities and
:math:`P_R` is a “replacement pressure” (see :cite:`Geiger98`, for example),
which serves to prevent residual ice motion due to spatial
variations of the ice strength :math:`P` when the strain rates are exactly zero.

An elliptical yield curve is used, with the viscosities given by

.. math::
   \zeta = {P(1+k_t)\over 2\Delta},
   :label: zeta

.. math::
   \eta  = e_g^{-2} \zeta,
   :label: eta

where

.. math::
   \Delta = \left[D_D^2 + {e_f^2\over e_g^4}\left(D_T^2 + D_S^2\right)\right]^{1/2}.
   :label: Delta

When the deformation :math:`\Delta` tends toward zero, the viscosities tend toward infinity. To avoid this issue, :math:`\Delta` needs to be limited and is replaced by :math:`\Delta^*` in equation :eq:`zeta`. Two methods for limiting :math:`\Delta` (or for capping the viscosities) are available in the code. If the namelist parameter ``capping_method`` is set to ``max``, :math:`\Delta^*=max(\Delta, \Delta_{min})` :cite:`Hibler79` while with ``capping_method`` set to ``sum``, the smoother formulation  :math:`\Delta^*=(\Delta + \Delta_{min})` of :cite:`Kreyscher00` is used. 

The ice strength :math:`P` is a function of the ice thickness distribution as
described in the `Icepack Documentation <https://cice-consortium-icepack.readthedocs.io/en/main/science_guide/index.html>`_.
 
Two other modifications to the standard VP rheology of :cite:`Hibler79` are available.
First, following the approach of :cite:`Konig10` (see also :cite:`Lemieux16`), the
elliptical yield curve can be modified such that the ice has isotropic tensile strength.
The tensile strength is expressed as a fraction of :math:`P`, that is :math:`k_t P`
where :math:`k_t` should be set to a value between 0 and 1 (this can
be changed at runtime with the namelist parameter ``Ktens``).

Second, while :math:`e_f` is the  ratio of the major and minor axes of the elliptical yield curve, the parameter
:math:`e_g` characterizes the plastic potential, i.e. another ellipse that decouples the flow rule from the
yield curve (:cite:`Ringeisen21`). :math:`e_f` and :math:`e_g` are respectively called ``e_yieldcurve`` and ``e_plasticpot`` in the code and
can be set in the namelist. The plastic potential can lead to more realistic fracture angles between linear kinematic features. :cite:`Ringeisen21` suggest to set :math:`e_f` to a value larger than 1 and to have :math:`e_g < e_f`.

By default, the namelist parameters are set to :math:`e_f=e_g=2` and :math:`k_t=0` which correspond to the standard VP rheology.

There are four options in the code for solving the sea ice momentum equation with a VP formulation: the standard EVP approach, a 1d EVP solver, the revised EVP approach and an implicit Picard solver. The choice of the capping method for the viscosities and the modifications to the yield curve and to the flow rule described above are available for these four different solution methods. Note that only the EVP and revised EVP methods are currently available if one chooses the C grid. 

.. _stress-evp:

Elastic-Viscous-Plastic
~~~~~~~~~~~~~~~~~~~~~~~

In the EVP model the internal stress tensor is determined from a
regularized version of the VP constitutive law :eq:`vp-const`.  The constitutive law is therefore

.. math::
   {1\over E}{\partial\sigma_1\over\partial t} + {\sigma_1\over 2\zeta}
     + {P_R\over 2\zeta} = D_D, \\
   :label: sig1

.. math::
   {1\over E}{\partial\sigma_2\over\partial t} + {\sigma_2\over 2\eta} = D_T,
   :label: sig2

.. math::
   {1\over E}{\partial\sigma_{12}\over\partial t} + {\sigma_{12}\over
     2\eta} = {1\over 2}D_S,
   :label: sig12


Viscosities are updated during the subcycling, so that the entire
dynamics component is subcycled within the time step, and the elastic
parameter :math:`E` is defined in terms of a damping timescale :math:`T`
for elastic waves, :math:`\Delta t_e < T < \Delta t`, as

.. math::
   E = {\zeta\over T},

where :math:`T=E_\circ\Delta t` and :math:`E_\circ` (elasticDamp) is a tunable
parameter less than one. Including the modification proposed by :cite:`Bouillon13` for equations :eq:`sig2` and :eq:`sig12` in order to improve numerical convergence, the stress equations become

.. math::
   \begin{aligned}
   {\partial\sigma_1\over\partial t} + {\sigma_1\over 2T}
     + {P_R\over 2T} &=& {\zeta \over T} D_D, \\
   {\partial\sigma_2\over\partial t} + {\sigma_2\over 2T} &=& {\eta \over
     T} D_T,\\
   {\partial\sigma_{12}\over\partial t} + {\sigma_{12}\over  2T} &=&
     {\eta \over 2T}D_S.\end{aligned}

Once discretized in time, these last three equations are written as

.. math::
   \begin{aligned}
   {(\sigma_1^{k+1}-\sigma_1^{k})\over\Delta t_e} + {\sigma_1^{k+1}\over 2T}
     + {P_R^k\over 2T} &=& {\zeta^k\over T} D_D^k, \\
   {(\sigma_2^{k+1}-\sigma_2^{k})\over\Delta t_e} + {\sigma_2^{k+1}\over 2T} &=& {\eta^k \over
     T} D_T^k,\\
   {(\sigma_{12}^{k+1}-\sigma_{12}^{k})\over\Delta t_e} + {\sigma_{12}^{k+1}\over  2T} &=&
     {\eta^k \over 2T}D_S^k,\end{aligned}
   :label: sigdisc


where :math:`k` denotes again the subcycling step. All coefficients on the left-hand side are constant except for
:math:`P_R`. This modification compensates for the decreased efficiency of including
the viscosity terms in the subcycling. Choices of the parameters used to define :math:`E`,
:math:`T` and :math:`\Delta t_e` are discussed in
Sections :ref:`revp` and :ref:`parameters`.

On the B grid, the stresses :math:`\sigma_{1}`, :math:`\sigma_{2}` and :math:`\sigma_{12}` are collocated at the U point. To calculate these stresses, the viscosities :math:`\zeta` and :math:`\eta` and the replacement pressure :math:`P_R` are also defined at the U point. 

However, on the C grid, :math:`\sigma_{1}` and :math:`\sigma_{2}` are collocated at the T point while :math:`\sigma_{12}` is defined at the U point. During a subcycling step, :math:`\zeta`, :math:`\eta` and :math:`P_R` are first calculated at the T point. To do so, :math:`\Delta` given by  equation :eq:`Delta` is calculated following the approach of :cite:`Bouillon13` (see also :cite:`Kimmritz16` for details). With this approach, :math:`D_S^2` at the T point is obtained by calculating :math:`D_S^2` at the U points and interpolating these values to the T point. As :math:`\sigma_{12}` is calculated at the U point, :math:`\eta` also needs to be computed as these locations. If ``visc_method`` in the namelist is set to ``avg_zeta`` (the default value), :math:`\eta` at the U point is obtained by interpolating T point values to this location. This corresponds to the approach used by :cite:`Bouillon13` and the one associated with the C1 configuration of :cite:`Kimmritz16`. On the other hand, if ``visc_method = avg_strength``, the strength :math:`P` calculated at T points is interpolated to the U point and :math:`\Delta` is calculated at the U point in order to obtain :math:`\eta` following equations :eq:`zeta` and :eq:`eta`. This latter approach is the one used in the C2 configuration of :cite:`Kimmritz16`.

.. _evp1d:

1d EVP solver
~~~~~~~~~~~~~

The standard EVP solver iterates hundreds of times, where each iteration includes a communication through MPI and a limited number of calculations. This limits how much the solver can be optimized as the speed is primarily determined by the communication. The 1d EVP solver avoids the communication by utilizing shared memory, which removes the requirement for calls to the MPI communicator. As a consequence of this the potential scalability of the code is improved. The performance is best on shared memory but the solver is also functional on MPI and hybrid MPI/OpenMP setups as it will run on the master processor alone.

The scalability of geophysical models is in general terms limited by the memory usage. In order to optimize this the 1d EVP solver solves the same equations that are outlined in the section :ref:`stress-evp` but it transforms all matrices to vectors (1d matrices) as this compiles better with the computer hardware. The vectorization and the contiguous placement of arrays in the memory makes it easier for the compiler to optimize the code and pass pointers instead of copying the vectors. The 1d solver is not supported for tripole grids and the code will abort if this combination is attempted.

.. _revp:

Revised EVP approach
~~~~~~~~~~~~~~~~~~~~

Introducing the numerical parameter :math:`\alpha=2T \Delta t_e ^{-1}` :cite:`Bouillon13`, the stress equations in :eq:`sigdisc` become

.. math::
   \begin{aligned}
   {\alpha (\sigma_1^{k+1}-\sigma_1^{k})} + {\sigma_1^{k}}
     + {P_R^k} &=& 2 \zeta^k D_D^k, \\
   {\alpha (\sigma_2^{k+1}-\sigma_2^{k})} + {\sigma_2^{k}} &=& 2 \eta^k D_T^k,\\
   {\alpha (\sigma_{12}^{k+1}-\sigma_{12}^{k})} + {\sigma_{12}^{k}} &=&
     \eta^k D_S^k,\end{aligned}

where as opposed to the classic EVP, the second term in each equation is at iteration :math:`k` :cite:`Bouillon13`. Also, contrary to the classic EVP,
:math:`\Delta t_e` times the number of subcycles (or iterations) does not need to be equal to the advective time step :math:`\Delta t`.
Finally, as with the classic EVP approach, the stresses are initialized using the previous time level values.
The revised EVP is activated by setting the namelist parameter ``revised_evp = true``.
In the code :math:`\alpha` is ``arlx`` and :math:`\beta` is ``brlx`` (introduced in Section :ref:`revp-momentum`). The values of ``arlx`` and ``brlx`` can be set in the namelist.
It is recommended to use large values of these parameters and to set :math:`\alpha=\beta` :cite:`Kimmritz15`.

.. _stress-eap:

Elastic-Anisotropic-Plastic
~~~~~~~~~~~~~~~~~~~~~~~~~~~

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
