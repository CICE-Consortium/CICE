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
available, whose values are contained in the ``trcrn`` array. Their
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
is given an integer index, ``trcr_depend``, which has the value 0, 1, or 2
depending on whether the appropriate conservation equation is
Equation :eq:`transport-aT`, Equation :eq:`transport-viT`, or Equation :eq:`transport-vsT`,
respectively. The total number of tracers is
:math:`N_{tr}\ge 1`.  Table :ref:`tracers` provides an overview of available tracers, 
including the namelist flags that turn them on and off, and their indices in the tracer 
arrays.  If any of the three explicit pond schemes is on, then ``tr_pond`` is true. 
Biogeochemistry tracers can be defined in the skeletal layer, dependent on the ice area
fraction, or through the full depth of snow and ice, in which case they utilize the 
bio grid and can depend on the brine fraction or the ice volume, if the brine fraction 
is not in use.

.. csv-table:: *Tracer flags and indices*
   :header: "flag", "num tracers", "dependency", "index (CICE grid)", "index (bio grid)"
   :widths: 12, 12, 18, 18, 18

   "default", "1", "aice", "nt_Tsfc=1", " "
   "default", "1", "vice", "nt_qice", " "
   "default", "1", "vsno", "nt_qsno", " "
   "default", "1", "vice", "nt_sice", " "
   "tr_iage", "1", "vice", "nt_iage", " "
   "tr_FY", "1", "aice", "nt_FY", " "
   "tr_lvl", "2", "aice", "nt_alvl", " "
   " ", " ", "vice", "nt_vlvl", " "
   "tr_pond_cesm", "2", "aice", "nt_apnd", " " 
   " ", " ", "apnd", "nt_vpnd", " "
   "tr_pond_lvl", "3", "aice", "nt_apnd", " " 
   " ", " ", "apnd", "nt_vpnd", " "
   " ", " ", "apnd", "nt_ipnd", " "
   "tr_pond_topo", "3", "aice", "nt_apnd", " " 
   " ", " ", "apnd", "nt_vpnd", " "
   " ", " ", "apnd", "nt_ipnd", " "
   "tr_aero", "n_aero", "vice, vsno", "nt_aero"," "
   "tr_brine", " ", "vice", "nt_fbri", " "
   "solve_zsal", "n_trzs", "fbri or (a,v)ice", "nt_bgc_S", " "
   "tr_bgc_N", "n_algae", "fbri or (a,v)ice", "nt_bgc_N", "nlt_bgc_N"
   "tr_bgc_Nit", " ", "fbri or (a,v)ice", "nt_bgc_Nit", "nlt_bgc_Nit"
   "tr_bgc_C", "n_doc", "fbri or (a,v)ice", "nt_bgc_DOC", "nlt_bgc_DOC"
   " ", "n_dic", "fbri or (a,v)ice", "nt_bgc_DIC", "nlt_bgc_DIC"
   "tr_bgc_chl", "n_algae", "fbri or (a,v)ice", "nt_bgc_chl", "nlt_bgc_chl"
   "tr_bgc_Am", " ", "fbri or (a,v)ice", "nt_bgc_Am", "nlt_bgc_Am"
   "tr_bgc_Sil", " ", "fbri or (a,v)ice", "nt_bgc_Sil", "nlt_bgc_Sil"
   "tr_bgc_DMS", " ", "fbri or (a,v)ice", "nt_bgc_DMSPp", "nlt_bgc_DMSPd"
   " ", " ", "fbri or (a,v)ice", "nt_bgc_DMSPd", "nlt_bgc_DMSPd"
   " ", " ", "fbri or (a,v)ice", "nt_bgc_DMS", "nlt_bgc_DMS"
   "tr_bgc_PON", " ", "fbri or (a,v)ice", "nt_bgc_PON", "nlt_bgc_PON"
   "tr_bgc_DON", " ", "fbri or (a,v)ice", "nt_bgc_DON", "nlt_bgc_DON"
   "tr_bgc_Fe", "n_fed", "fbri or (a,v)ice", "nt_bgc_Fed", "nlt_bgc_Fed"
   " ", "n_fep", "fbri or (a,v)ice", "nt_bgc_Fep", "nlt_bgc_Fep"
   "tr_bgc_hum", " ", "fbri or (a,v)ice", "nt_bgc_hum", "nlt_bgc_hum"
   "tr_zaero", "n_zaero", "fbri or (a,v)ice", "nt_zaero", "nlt_zaero"
   " ", "1", "fbri", "nt_zbgc_frac", " "



Users may add any number of additional tracers that are transported conservatively,
provided that the dependency ``trcr_depend`` is defined appropriately. 
See Section :ref:`addtrcr` for guidance on adding tracers.

Please see the `Icepack documentation <https://cice-consortium-icepack.readthedocs.io/en/master/science_guide/index.html>`_ for additional information about tracers that depend on other tracers, age of the ice, aerosols, 
brine height, and the sea ice ecosystem.
