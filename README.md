This repo is an extension of MPAS-Atmosphere v7.3 to include significant improvements to the potential vorticity (PV) diagnostics package. This revised package has greater functionality than the original package and contains several important fixes and modifications to close the PV budget, which are detailed below.


!=====================================================================================================================
! pv_diagnostics.F: MPAS code to compute Ertel's potential vorticity and tendency contributions to the Eulerian PV budget
!=====================================================================================================================
! Original diagnostics code was written by Nick Szapiro - 2016
! Significant changes to the code were made by Manda Chasteen (chasteen@ucar.edu) and May Wong (mwong@ucar.edu) - 2023
! Note: this revised PV diagnostics package is heavily reliant upon tendency calculations associated with the ITM
!       tendency package (also included in this repo)
!
! Added namelist options for ease of toggling on PV diagnostics calculations
! -- config_pv_diag          : flag for whether the 3D PV field and fields interpolated to dynamic tropopause are desired
! -- config_pv_tend          : flag for whether PV tendency diagnostics are desired (required for config_pv_microphys, config_pv_isobaric)
! -- config_pv_scalar        : flag for whether pv_scalar is initialized as PV and then transported as passive scalar
! -- config_pv_microphys     : flag for whether specific microphysics process PV tendencies are desired (Thompson only)
! -- config_pv_isobaric      : flag for whether isobaric interpolation of PV diagnostics variables is desired
!
! --------------------------------------------------------------------------------------------------------------------
! Subroutines contained in pv_diagnostics.F:
! ------------------------------------------
! pv_diagnostics_setup    : setup diagnostics package and performs initial check of PV config flags
! pv_diagnostics_reset    : calls store_previous_vars to save previous timestep variables and pv_diagnostics_dyn_init
!                         : to initialize and reset the computed theta and momentum tendencies as zeros
! pv_diagnostics_update   : calls atm_compute_pv_diagnostics and atm_compute_pvBudget_diagnostics to compute PV, all
!                           PV tendency variables, and interpolation of variables onto identified dynamic tropopause
! pv_diagnostics_init     : initializes PV scalar variable if desired and not a restart run. called in a
!                           new subroutine mpas_atm_diag_pv_init() that is then explicitly called in mpas_atm_core.F
! pv_diagnostics_compute  : calls MPAS_field_will_be_written for various quantities.
!                           MC_TODO NOTE: I'm not sure that this is necessary... variables are written if included in streams
!
!
! Order that subroutines are called in mpas_atm_core:
! -- At model intiialziation: pv_diagnostics_reset, pv_diagnostics_update, pv_diagnostics_init, pv_diagnostics_compute,
!                             pv_diagnostics_reset
!
! -- During time step integration: pv_diagnostics_update, pv_diagnostics_compute, pv_diagnostics_reset
! ********************************************************************************************************************
!
! Changes made from the original PV diagnostics code include:
  ! -----------------------------------------------------------
  ! * Different formulation for calculation of horizontal gradients on native MPAS grid. The updated method is based on
  !   Eq. 22 in Ringler et al. (2010) and is more robust than the previous method implemented by NS
  ! * Reconstruction of horizontal gradients on each cell's edges to the cell center following the same method as the
  !   horizontal wind reconstruction in mpas_vector_reconstruction.F
  ! * Changes to the calculation of the PV tendency terms to ensure that the correct time levels are used for the
  !   coefficients, as determined by discretizing the equation for PV. We employ consistent time levels for all relevant
  !   PV tendencies computed in MPAS:
  !   -- in diabatic PV tendencies, the 3D absolute vorticity vector from time level t
  !   -- in frictional PV tendencies, the 3D potential temperature gradient from time level t+dt
  !   -- density from t+dt is used in all relevant calculations 
  !   This important change requires storing fields from the beginning of the time step to be used in the PV tendency
  !   calculations because the model state and diagnostic fields are updated and assigned to time level 1 before the
  !   PV diagnostics are called at the end of the time step in mpas_atm_core.F. Thus, before this change was
  !   implemented, the updated variables from the end of the time step were incorrectly used alongside all these
  !    tendencies.
  ! * Update required to mpas_atm_core.F to ensure that diagnostic quantities theta and rho are updated at each time
  !   step. Previously, these were only calculated if alarm bell for writing an outfile was activated
  ! * Split frictional tendencies into components from explicit mixing, PBL+GWD schemes, and cumulus schemes, which
  !   are then summed to produce the full frictional tendency depv_dt_fric. This required the introduction of individual
  !   momentum tendency variables and renders the original tend_u_phys term obsolute, which has therefore been removed.
  !   These tendencies are derived from the coupled momentum tendencies rather than taking the uncoupled tendencies
  !   directly from physics. 
  ! * Corrections were made to the diffusion friction tendency terms, which had previously called tend_u_euler and
  !   tend_w_euler variables that comprised other momentum tendencies in addition to diffusion. These required calculating
  !   additional variables, u_tend_diff and w_tend_diff, in mpas_atm_time_integration.F that contain only the tendency
  !   contributions from diffusion.
  ! * The potential temperature tendency (dtheta_dt_mix) that is input into the diabatic diffusion tendency calculation
  !   was initially coupled to mass, which needed to be fixed. The tendency now is computed by decoupling the theta_m
  !   tendency associated with mixing from moisture (calculated in mpas_atm_time_integration), which is more accurate
  !   and enables closing the theta and PV budgets. 
  ! * All physics diabatic tendencies have been modified to use the derived theta tendencies by decoupling the associated
  !   theta_m tendencies from moisture, rather than the theta tendencies output directly from the physics schemes. Doing so
  !   is more accurate and enables closing the theta and PV budgets. 
  ! * Modified interpolation of PV tendencies to dynamic tropopause routine to interpolate to the DT identified at the
  !   beginning of the time step rather than at the end. This provides a better depiction of how processes may alter the
  !   height of the DT over the time step
  ! * Modified floodFill_tropo routine to better identify the dynamic tropopause in regions with low and/or negative
  !   PV values aloft.
  ! * Modified the DT interpolation routine (interp_pv) to mitigate prior issues of interpolating values to a falsely
  !   identified DT point where the bounding levels didn't change from (sign(f)*PV) < 2 PVU to (sign(f)*PV) >= 2 PVU.
  !   Interpolation weights assume this is true, leading to erroneous values of interpolated fields.
  !
  ! New additions include:
  ! ------------------------------------------
  ! * Inclusion of dynamics tendencies for all relevant variables, enabling the dynamics (advective) contributions to the PV
  !   budget to be accurately evaluated. The PV tendencies from dynamics do not include the effects of explicit diffusion,
  !   which are included as diabatic and frictional PV tendencies.
  ! * Incorporation of a PV passive scalar variable to advect initial PV field via the dynamics scalar transport routine
  !   throughout the model integration. Requires config_pv_scalar = .true.
  !   Note: using the PV scalar variable is a proxy for adiabatic PV transport and is not an adequate substitution
  !         for the dynamics tendencies (i.e., the PV budget will not close if scalar transport is used in lieu of the
  !         PV dynamics tendencies).
  ! * Accumulated PV tendencies were added to permit the evaluation of the net PV tendencies without outputting the model
  !   variables at each time step.
  ! * Added PV tendencies for specific microphysical processes in the Thompson scheme: net condensation/evaporation of cloud
  !   water, evaporation of rain water, net deposition/sublimation, melting, and freezing. Requires config_pv_microphys = .true.
  !   Note: these tendencies use the theta tendencies from the microphysics scheme directly, whereas depv_dt_mp is calculated
  !   using the derived theta tendency from the theta_m and qv tendencies. The differences in these approaches can be
  !   ascertained by comparing depv_dt_mp to depv_dt_mp_allproc
  ! * Incorporation of routine to interpolate PV diagnostics to isobaric levels (code also modified in isobaric_diagnostics.F)
  !   and then accumulate the interpolated tendencies to isobaric levels. Requires config_pv_isobaric = .true.
  !   Note: changes to this procedure requires making changes to isobaric_diagnostics.F and Registry_isobaric.xml

