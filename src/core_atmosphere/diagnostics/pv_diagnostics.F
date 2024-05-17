! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the
! LICENSE file
! distributed with this code, or at
! http://mpas-dev.github.com/license.html
!=================================================================================================================

module pv_diagnostics

    use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
    use mpas_kind_types, only : RKIND, StrKIND
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: tend
    type (MPAS_pool_type), pointer :: tend_physics
    type (MPAS_pool_type), pointer :: diag_physics
    type (MPAS_pool_type), pointer :: configs 

    type (MPAS_clock_type), pointer :: clock

    public :: pv_diagnostics_setup, &                          
              pv_diagnostics_compute, &                        
              pv_diagnostics_reset, &
              pv_diagnostics_update , &                         
              pv_diagnostics_init                              ! If config_pv_scalar, need to initialize as initial PV field 

    private

    ! MC_TODO: determine if these are needed; if so, make sure all relevant flags are included 
    logical :: need_ertel_pv, need_u_pv, need_v_pv, need_pres_pv, need_height_pv, need_theta_pv, &
               need_vort_pv, need_iLev_DT, &
               need_wcell, &                  
               need_tend_lw, need_tend_sw, need_tend_bl, need_tend_cu, need_tend_mix, need_tend_mp, &
               need_tend_dyn, &
               need_tend_diab, need_tend_fric, need_tend_diab_pv, need_tend_fric_pv, need_tend_dyn_pv, need_dtheta_mp, need_dtheta_cu, & 
               need_uzonal_prev, need_umerid_prev, need_wcell_prev, need_rho_prev, need_theta_prev, &                   
               need_qv_prev, & ! MW added for decoupling friction term 
               need_pvvert_prev, need_epv_prev, need_iLev_prev, need_epv_scalar, need_tends_mp_proc, &
               need_acc_tend_lw, need_acc_tend_sw, need_acc_tend_bl, need_acc_tend_cu, need_acc_tend_mix, &
               need_acc_tend_mp, need_acc_tend_diab, need_acc_tend_fric, need_acc_tend_diab_pv, need_acc_tend_fric_pv, &
               need_acc_tend_dyn_pv, need_acc_tend_dyn

    !=====================================================================================================================
    ! MPAS code to compute Ertel's potential vorticity and tendency contributions to the Eulerian PV budget
    !=====================================================================================================================
    ! Original diagnostics code was written by Nick Szapiro - 2016
    ! Significant changes to the code were made by Manda Chasteen (chasteen@ucar.edu) and May Wong (mwong@ucar.edu) - 2023
    ! Note: this revised PV diagnostics package is heavily reliant upon tendency calculations associated with the ITM
    !       tendency package. (MC_TODO: make sure this is stated accurately) 
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

    contains

    !*********************************************************************************************************************
    ! pv_diagnostics_setup: initialize the PV diagnostics when called in mpas_atm_diagnostics_manager.F
    !*********************************************************************************************************************

    subroutine pv_diagnostics_setup(configs_in, all_pools, simulation_clock)
        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type, MPAS_STREAM_OUTPUT, MPAS_STREAM_INPUT, &
                                       MPAS_STREAM_INPUT_OUTPUT, MPAS_LOG_WARN
        use mpas_pool_routines, only : mpas_pool_get_subpool, mpas_pool_get_config
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none
        
        type (MPAS_pool_type), pointer :: configs_in
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        logical, pointer :: config_tend, config_pv_diag, config_pv_tend, config_pv_scalar, &
                            config_pv_microphys, config_pv_isobaric, config_apply_lbcs
        character(len=StrKIND), pointer :: config_microp_scheme

        ! for zgrid_cell initialization
        integer, pointer :: nCells, nVertLevels
        real(kind=RKIND), dimension(:,:), pointer :: zgrid, zCell
        
        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'tend', tend)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)

        clock => simulation_clock
        configs => configs_in

        ! check configs.
        call mpas_pool_get_config(configs, 'config_apply_lbcs', config_apply_lbcs)
        call mpas_pool_get_config(configs,'config_microp_scheme',config_microp_scheme)
        call mpas_pool_get_config(configs, 'config_tend', config_tend)
        call mpas_pool_get_config(configs, 'config_pv_diag', config_pv_diag)  
        call mpas_pool_get_config(configs, 'config_pv_tend', config_pv_tend)
        call mpas_pool_get_config(configs, 'config_pv_scalar', config_pv_scalar)  
        call mpas_pool_get_config(configs, 'config_pv_microphys', config_pv_microphys)  
        call mpas_pool_get_config(configs, 'config_pv_isobaric', config_pv_isobaric)  

        ! if limited area simulation, disable all PV flags if activated
        if (config_apply_lbcs .and. (config_pv_diag .or. config_pv_tend .or. config_pv_scalar &
             .or. config_pv_microphys .or. config_pv_isobaric)) then

            call mpas_log_write('PV diagnostics are not supported for limited-area simulations. Disabling.', MPAS_LOG_WARN)
            config_pv_diag = .false.
            config_pv_tend = .false.
            config_pv_scalar = .false.
            config_pv_microphys = .false.
            config_pv_isobaric = .false.
        end if 

        ! if dependent config_pv flags are activated but parent config_pv_diag flag is not, deactivate them. 
        if ((.not. config_pv_diag) .and. (config_pv_tend .or. config_pv_scalar .or. config_pv_microphys .or. config_pv_isobaric)) then
            config_pv_tend = .false.
            config_pv_scalar = .false.
            config_pv_microphys = .false.
            config_pv_isobaric = .false.
            call mpas_log_write('config_pv_diag is not activated; deactivated all dependent PV configs.', MPAS_LOG_WARN)
        end if 

        ! if config_pv_tend is activated but config_tend is not, deactivate.
        if ((.not. config_tend) .and. (config_pv_tend .or. config_pv_microphys .or. config_pv_isobaric)) then
            config_pv_tend = .false.
            config_pv_microphys = .false.
            config_pv_isobaric = .false.
            call mpas_log_write('config_tend is not activated; deactivated all dependent PV configs.', MPAS_LOG_WARN)
        end if 
    
        ! if config_pv_microphys or config_pv_isobaric is activated but config_pv_tend is not, deactivate.
        if ((.not. config_pv_tend) .and. (config_pv_microphys .or. config_pv_isobaric)) then
            config_pv_microphys = .false.
            config_pv_isobaric = .false.
            call mpas_log_write('config_pv_tend is not activated; deactivated all dependent PV configs.', MPAS_LOG_WARN)
        end if

        if ((config_pv_microphys) .and. (config_microp_scheme /= 'mp_thompson')) then
            call mpas_log_write('config_pv_microphys is not compatible with   = '''//trim(config_microp_scheme)//''' -- disabling', MPAS_LOG_WARN)
            config_pv_microphys = .false.
        end if 

        ! If doing PV, need to initialize zCell array for vertical derivative calculations
        if (config_pv_diag) then
           call mpas_pool_get_dimension(mesh, 'nCells', nCells)
           call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

           call mpas_pool_get_array(mesh, 'zgrid', zgrid)
           call mpas_pool_get_array(diag, 'zgrid_cell', zCell)

           call interp_wLev_thetaLev(zgrid, nCells, nVertLevels, zCell)
           call mpas_log_write("Initialized zgrid_cell array in setup.")

        end if

  
    end subroutine pv_diagnostics_setup

    !*********************************************************************************************************************
    ! pv_diagnostics_init:   A subroutine was created in mpas_atm_diagnostics_manager.F that calls this subroutine and 
    !                        is then called explicitly in mpas_atm_core.F during the first time step so that initial PV
    !                        field is populated for scalar advection, if desired. Config flag set here to disable 
    !                        the re-initialization of pv_scalar if restart run.
    !*********************************************************************************************************************

    subroutine pv_diagnostics_init()
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config

        implicit none
        
        logical, pointer :: config_pv_scalar, config_do_restart
        ! TEST
        integer, pointer :: index_pv_scalar
        real(kind=RKIND), dimension(:,:,:), pointer :: pv_scalars

        call mpas_pool_get_config(configs, 'config_pv_scalar', config_pv_scalar)  
        call mpas_pool_get_config(configs, 'config_do_restart', config_do_restart) 

        if ((.not. config_do_restart) .and. config_pv_scalar) then 
          ! Initialize scalar PV field
           call mpas_log_write('--- config_do_restart = $l', logicArgs=(/config_do_restart/))
           call mpas_log_write("Calling pv_diagnostics_scalar_init().")
           call pv_diagnostics_scalar_init()
        else if (config_do_restart .and. config_pv_scalar) then 
            call mpas_log_write('config_do_restart = $l', logicArgs=(/config_do_restart/))
            call mpas_log_write('--- skipping pv_scalar initialization.')
            call mpas_pool_get_array(state, 'pv_scalars', pv_scalars)     ! MC may not be necessary 
            return
        end if 
       
    end subroutine pv_diagnostics_init

    !*********************************************************************************************************************
    ! pv_diagnostics_scalar_init: initialize PV scalar field to be advected by model's transport scheme, if desired
    !*********************************************************************************************************************

    subroutine pv_diagnostics_scalar_init()
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config

        ! Need to do halo exchange prior to calculating horizontal gradients in calc_epv
        use mpas_derived_types, only : field2DReal
        use mpas_pool_routines, only : mpas_pool_get_field
        use mpas_dmpar, only : mpas_dmpar_exch_halo_field

        implicit none 
        
        logical, pointer :: config_pv_diag, config_pv_scalar
        integer, pointer :: index_pv_scalar

        real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
        real(kind=RKIND), dimension(:,:,:), pointer :: pv_scalars
        type (field2DReal), pointer :: theta_f, uReconstructZonal_f, uReconstructMeridional_f, w_f, wCell_f, &  ! MC: needed for halo exchange
                                       rho_f, pv_vertex_f

        call mpas_pool_get_config(configs, 'config_pv_diag', config_pv_diag)
        call mpas_pool_get_config(configs, 'config_pv_scalar', config_pv_scalar)

        call mpas_pool_get_dimension(state, 'index_pv_scalar', index_pv_scalar)                     ! Note: PV scalars defined in Registry.xml
        call mpas_pool_get_array(state, 'pv_scalars', pv_scalars)                                   ! Note: PV scalars defined in Registry.xml
        call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
        
        call mpas_pool_get_field(diag, 'theta', theta_f)         
        call mpas_pool_get_field(diag, 'uReconstructZonal', uReconstructZonal_f)                     ! MC added
        call mpas_pool_get_field(diag, 'uReconstructMeridional', uReconstructMeridional_f)           ! MC added
        call mpas_pool_get_field(state, 'w', w_f, 1)                                                 ! MC added
        call mpas_pool_get_field(diag, 'wCell', wCell_f)                                             ! MC added 
        call mpas_pool_get_field(diag, 'rho', rho_f)                                                 ! MC added (MC_TODO: not sure if rho needed)
        call mpas_pool_get_field(diag, 'pv_vertex', pv_vertex_f)                                     ! MC added
                                         
        ! Need halo cells for everything w/ horizontal derivative 
        call mpas_dmpar_exch_halo_field(theta_f)
        call mpas_dmpar_exch_halo_field(uReconstructZonal_f)                               
        call mpas_dmpar_exch_halo_field(uReconstructMeridional_f)                         
        call mpas_dmpar_exch_halo_field(w_f)
        call mpas_dmpar_exch_halo_field(wCell_f)                                    
        call mpas_dmpar_exch_halo_field(rho_f)                                             ! MC_TODO: not sure if needed
        call mpas_dmpar_exch_halo_field(pv_vertex_f)                                       

        ! Calculate PV and initialize PV scalar variable as initial PV field 
        if (config_pv_diag) then
           call calc_epv(mesh, 1, state, diag) 
           
           if (config_pv_scalar) then
              pv_scalars(index_pv_scalar,:,:) = ertel_pv(:,:)  
           end if 
        end if 

    end subroutine pv_diagnostics_scalar_init

    !*********************************************************************************************************************
    ! pv_diagnostics_dyn_init: initialize/reset the dynamics and physics tendency fields that are used to compute the PV 
    !                          tendencies. This is called every time step in pv_diagnostics_reset()
    !
    ! MC_TODO: CLEAN UP 
    !*********************************************************************************************************************

    subroutine pv_diagnostics_dyn_init()
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
        use mpas_pool_routines, only : mpas_pool_get_config

        implicit none 
        
        logical, pointer :: config_pv_tend
    
        integer :: k
        integer, pointer :: nCells, nVertLevels
        real(kind=RKIND), dimension(:,:), pointer :: rho
        ! dynamics tendency variables
        real(kind=RKIND), dimension(:,:), pointer :: dthetam_dt_dyn, du_dt_dyn, dw_dt_dyn  
        ! physics tendency variables
        real(kind=RKIND), dimension(:,:), pointer :: tend_u_cu, tend_u_pbl, thmblten, thmcuten, &
                                                     thmswten, thmlwten, thmmpten

        ! NOTE: may not all be needed. clean up

        call mpas_pool_get_config(configs, 'config_pv_tend', config_pv_tend)
        
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_array(diag, 'rho', rho)
        
        call mpas_pool_get_array(diag, 'dthetam_dt_dyn', dthetam_dt_dyn)
        call mpas_pool_get_array(diag, 'du_dt_dyn', du_dt_dyn)
        call mpas_pool_get_array(diag, 'dw_dt_dyn', dw_dt_dyn)

        call mpas_pool_get_array(diag, 'tend_u_pbl', tend_u_pbl)
        call mpas_pool_get_array(diag, 'tend_u_cu', tend_u_cu)
        call mpas_pool_get_array(diag, 'thmblten', thmblten)
        call mpas_pool_get_array(diag, 'thmcuten', thmcuten)
        call mpas_pool_get_array(diag, 'thmswten', thmswten)
        call mpas_pool_get_array(diag, 'thmlwten', thmlwten)
        call mpas_pool_get_array(diag, 'thmmpten', thmmpten)

        if (config_pv_tend) then
           ! Dynamics tendencies
           dthetam_dt_dyn = 0.
           du_dt_dyn      = 0.
           dw_dt_dyn      = 0.

           ! Physics u tendencies 
           tend_u_pbl     = 0.
           tend_u_cu      = 0.

           ! Physics diabatic tendencies --- these may not be needed MC_TODO
           thmblten       = 0.
           thmcuten       = 0.
           thmlwten       = 0.
           thmswten       = 0.
           thmmpten       = 0.
        end if 

        ! Avoid FP errors caused by a potential division by zero below by 
        ! initializing the "garbage cell" of rho to a non-zero value
        do k=1,nVertLevels
           rho(k,nCells+1) = 1.0
        end do

    end subroutine pv_diagnostics_dyn_init

    !*********************************************************************************************************************
    ! pv_diagnostics_update: compute the PV diagnostics over each time step when called in mpas_atm_diagnostics_manager.F
    !*********************************************************************************************************************

    subroutine pv_diagnostics_update()
        use mpas_pool_routines, only: mpas_pool_get_config

        implicit none
        
        logical, pointer :: config_pv_diag, config_pv_tend

        call mpas_pool_get_config(configs, 'config_pv_diag', config_pv_diag)
        call mpas_pool_get_config(configs, 'config_pv_tend', config_pv_tend)

        if (config_pv_diag) then
           ! Calculate PV diagnostics at each time step
           call mpas_log_write("Computing Ertel's PV.")
           call atm_compute_pv_diagnostics(configs, state, 1, diag, mesh)
 
           if (config_pv_tend) then
              call mpas_log_write('Computing PV tendency diagnostics.')
              call atm_compute_pvBudget_diagnostics(configs, state, 1, diag, mesh, tend, tend_physics, diag_physics)
           end if
      end if 

    end subroutine pv_diagnostics_update

    !*********************************************************************************************************************
    ! pv_diagnostics_compute: compute subroutine is only called prior to writing fields to outfile.
    !
    ! MC -- CLEAN UP. NOT CLEAR THAT THIS DOES ANYTHING...
    !*********************************************************************************************************************

    subroutine pv_diagnostics_compute()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_pool_routines, only: mpas_pool_get_config

        implicit none
        
        logical, pointer :: config_pv_diag, config_pv_scalar

        call mpas_pool_get_config(configs, 'config_pv_diag', config_pv_diag)
        call mpas_pool_get_config(configs, 'config_pv_scalar', config_pv_scalar)

        ! MC 02/08/2024: I'm not sure these even do anything as long as variables are in output streams
        ! if we could remove them without repercussions, would clean up code a lot.
        ! MC_TODO:

        if (config_pv_diag) then
           ! Ertel's PV and state vars interpolated to dynamic tropopause
           need_ertel_pv = MPAS_field_will_be_written('ertel_pv')
           need_u_pv = MPAS_field_will_be_written('u_pv')
           need_v_pv = MPAS_field_will_be_written('v_pv')
           need_theta_pv = MPAS_field_will_be_written('theta_pv')
           need_pres_pv = MPAS_field_will_be_written('pres_pv')
           need_height_pv = MPAS_field_will_be_written('height_pv')
           need_vort_pv = MPAS_field_will_be_written('vort_pv')
           need_iLev_DT = MPAS_field_will_be_written('iLev_DT')
           need_wCell = MPAS_field_will_be_written('wCell')

           ! PV tendency terms 
           need_tend_lw = MPAS_field_will_be_written('depv_dt_lw')
           need_tend_sw = MPAS_field_will_be_written('depv_dt_sw')
           need_tend_bl = MPAS_field_will_be_written('depv_dt_bl')
           need_tend_cu = MPAS_field_will_be_written('depv_dt_cu')
           need_tend_mix = MPAS_field_will_be_written('depv_dt_mix')
           need_dtheta_mp = MPAS_field_will_be_written('dtheta_dt_mp')
           need_dtheta_cu = MPAS_field_will_be_written('dtheta_dt_cu')
           need_tend_mp = MPAS_field_will_be_written('depv_dt_mp')
           need_tend_dyn = MPAS_field_will_be_written('depv_dt_dyn')
           need_tend_diab = MPAS_field_will_be_written('depv_dt_diab')
           need_tend_fric = MPAS_field_will_be_written('depv_dt_fric')
           need_tend_diab_pv = MPAS_field_will_be_written('depv_dt_diab_pv')
           need_tend_fric_pv = MPAS_field_will_be_written('depv_dt_fric_pv')
           need_tend_dyn_pv = MPAS_field_will_be_written('depv_dt_dyn_pv')

           ! Accumulated PV tendency terms
           need_acc_tend_lw = MPAS_field_will_be_written('acc_depv_dt_lw')
           need_acc_tend_sw = MPAS_field_will_be_written('acc_depv_dt_sw')
           need_acc_tend_bl = MPAS_field_will_be_written('acc_depv_dt_bl')
           need_acc_tend_cu = MPAS_field_will_be_written('acc_depv_dt_cu')
           need_acc_tend_mp = MPAS_field_will_be_written('acc_depv_dt_mp')
           need_acc_tend_mix = MPAS_field_will_be_written('acc_depv_dt_mix')
           need_acc_tend_diab = MPAS_field_will_be_written('acc_depv_dt_diab')
           need_acc_tend_fric = MPAS_field_will_be_written('acc_depv_dt_fric')
           need_acc_tend_dyn = MPAS_field_will_be_written('acc_depv_dt_dyn')
           need_acc_tend_diab_pv = MPAS_field_will_be_written('acc_depv_dt_diab_pv')
           need_acc_tend_fric_pv = MPAS_field_will_be_written('acc_depv_dt_fric_pv')
           need_acc_tend_dyn_pv = MPAS_field_will_be_written('acc_depv_dt_dyn_pv')

           ! Specific microphysical process tendencies
           need_tends_mp_proc = MPAS_field_will_be_written('depv_dt_mp_evap_cw')
           need_tends_mp_proc = MPAS_field_will_be_written('depv_dt_mp_evap_rw')
           need_tends_mp_proc = MPAS_field_will_be_written('depv_dt_mp_depo_ice')
           need_tends_mp_proc = MPAS_field_will_be_written('depv_dt_mp_melt_ice')
           need_tends_mp_proc = MPAS_field_will_be_written('depv_dt_mp_frez_ice')
           need_tends_mp_proc = MPAS_field_will_be_written('depv_dt_mp_allproc')

           need_tends_mp_proc = MPAS_field_will_be_written('acc_depv_dt_mp_evap_cw')
           need_tends_mp_proc = MPAS_field_will_be_written('acc_depv_dt_mp_evap_rw')
           need_tends_mp_proc = MPAS_field_will_be_written('acc_depv_dt_mp_depo_ice')
           need_tends_mp_proc = MPAS_field_will_be_written('acc_depv_dt_mp_melt_ice')
           need_tends_mp_proc = MPAS_field_will_be_written('acc_depv_dt_mp_frez_ice')
           need_tends_mp_proc = MPAS_field_will_be_written('acc_depv_dt_mp_allproc')
           
           ! Save fields from beginning of time step to outfile before updating them
           call mpas_log_write('Saving beginning of time step variables to outfile.')
           need_uzonal_prev = MPAS_field_will_be_written('uReconstructZonal_prev')
           need_umerid_prev = MPAS_field_will_be_written('uReconstructMeridional_prev')
           need_wcell_prev = MPAS_field_will_be_written('wCell_prev')
           need_theta_prev = MPAS_field_will_be_written('theta_prev')
           need_qv_prev = MPAS_field_will_be_written('qv_prev')
           need_rho_prev = MPAS_field_will_be_written('rho_prev')
           need_pvvert_prev = MPAS_field_will_be_written('pv_vertex_prev')
           need_epv_prev = MPAS_field_will_be_written('ertel_pv_prev')
           need_iLev_prev = MPAS_field_will_be_written('iLev_DT_prev')

           if (config_pv_scalar) then
              ! Write pv_scalar 
              need_epv_scalar = MPAS_field_will_be_written('pv_scalars')
           end if
           
        end if
        
    end subroutine pv_diagnostics_compute

    !*********************************************************************************************************************
    ! pv_diagnostics_reset: Update beginning of time step fields after they have been written to outfile for use in next 
    !                       time step PV tendency calculations when called in mpas_atm_diagnostics_manager.F
    !*********************************************************************************************************************

    subroutine pv_diagnostics_reset()
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_pool_routines, only: mpas_pool_get_config

        implicit none
        
        logical, pointer :: config_pv_tend
        
        call mpas_pool_get_config(configs, 'config_pv_tend', config_pv_tend)

        if (config_pv_tend) then
            ! Called immediately after diagnostics have been written
            ! Update previous variables for next timestep calculation
            call mpas_log_write('Updating previous fields.')
            call store_previous_vars(mesh, 2, state, diag)
    
            ! Repopulate the fields used to calculate dynamics and physics tendencies with zeros 
            call pv_diagnostics_dyn_init()
        end if 
        
    end subroutine pv_diagnostics_reset
     
    !*********************************************************************************************************************
    ! NS: Below are two subroutines (floodFill_strato and floodFill_tropo), designed to determine the first model level 
    ! above the dynamic tropopause, iLev_DT, which is designated as the 2-PVU isosurface. Only one of these subroutines 
    ! is used (toggled with "call floodFill_strato(mesh, diag, pvuVal, stratoPV)" and "call floodFill_tropo(mesh,diag,pvuVal)"
    ! in the atm_compute_pv_diagnostics subroutine below. The routines *should* produce equivalent estimates for iLev_DT. 
    !
    ! MC - 2023: floodFill_tropo was updated to better ID the dynamic tropopause. This routine should be used in lieu of 
    !            floodFill_strato, and the two should no longer be expected to produce equivalent iLev_DT estimates.
    !********************************************************************************************************************* 
   
    subroutine floodFill_strato(mesh, diag, pvuVal, stratoPV)
       !To find model level of dynamic tropopause: 
       !Simply searching down from TOA within each column to find first
       !instance of 2-PVU surface (i.e., where the PV drops below values
       !characteristic of the stratosphere) 
       !is buggy due to stratospheric wave breaking, which may induce
       !regions of low PV (i.e., PV < 2 PVU) within the stratosphere and
       !thus yield artifically
       !high estimations of the tropopause height. This seems to be more
       !problematic as the mesh gets finer and the vertical vorticity
       !field exhibits greater variability
       !or jumps. 
       !Note that these low-PV anomalies in the stratosphere may persist
       !for long times w/ slow mixing downstream of mountains or deep
       !convection.
       !A few quicker fixes (e.g., make sure PV < 2 PVU for a number of
       !layers; search down from 10 PVU instead of TOA) are hacky and not
       !robust.
       
       !To (hopefully) alleviate the problems resulting from wave
       !breaking, we can flood fill from a known
       !stratosphere region (e.g., where the model top > 2 PVU) and
       !filter down and around any problematic regions.
       !The problem w/ using only the flood fill is that strong surface
       !PV anomalies can connect to the 2-PVU surface aloft,
       !and the resulting "flood-filled 2 PVU" region can have sizeable
       !areas that are located just at/near the surface, while there is
       !clearly a 
       !tropopause above (i.e., as evident in a vertical cross-section).
       !To address the large near-surface blobs of PV > 2 PVU, will take
       !the flood fill mask and try to move upward from near the surface
       !to 10 PVU within a vertical column. 
       !If this can be done, then the low-level PV anomaly extends to the
       !stratosphere. Else, remove the stratospheric designation to
       !disconnect the "surface blob".
       
       !The "output" is iLev_DT, which is the vertical index for the
       !model level just above the dynamic tropopause (i.e., where PV >=
       !pvuVal, which is set below in atm_compute_pv_diagnostics to 2
       !PVU). 
       !If iLev_DT > nVertLevels, then pvuVal is found only above the
       !column (i.e., entire column is in troposphere). If iLev_DT < 1,
       !PV >= pvuVal extends vertically through the entire column 
       !(i.e., the entire column is within the stratosphere).
       !Communication between blocks during the flood fill may be needed
       !to treat some edge cases appropriately.
 
       use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
       use mpas_dmpar, only : mpas_dmpar_max_int,mpas_dmpar_exch_halo_field
       use mpas_derived_types, only : dm_info, field2DInteger
      
       implicit none
       
       type (mpas_pool_type), intent(in) :: mesh
       type (mpas_pool_type), intent(inout) :: diag
       real(kind=RKIND), intent(in) :: pvuVal, stratoPV
       
       integer :: iCell, k, nChanged, iNbr, iCellNbr, levInd, haloChanged, global_haloChanged !INCORPORATE LEVEL INDEX FOR REMOVING SFC BLOB
       integer, pointer :: nCells, nVertLevels
       integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
       integer, dimension(:,:), pointer :: cellsOnCell, inStrato ! inStrato wasn't here in original procedure
       
       type (field2DInteger), pointer :: inStrato_f ! line added to match troposphere procedure workflow
       
       real(kind=RKIND) :: sgnHemi, sgn
       real(kind=RKIND),dimension(:),pointer:: latCell
       real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
       
       type (dm_info), pointer :: dminfo
       
       integer, dimension(:,:), allocatable :: candInStrato ! whether point is potentially inStrato
       
       call mpas_pool_get_dimension(mesh, 'nCells', nCells)
       call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
       call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(mesh, 'latCell', latCell)
 
       call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
       call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
       call mpas_pool_get_array(diag, 'inStrato', inStrato) ! was allocated below in original script
       
       allocate(candInStrato(nVertLevels, nCells+1))
       !allocate(inStrato(nVertLevels, nCells+1))
       candInStrato(:,:) = 0
       inStrato(:,:) = 0
       
       !store whether each grid point has |PV| >= pvuVal to avoid
       !repeating logic. we'll use candInStrato as a isVisited marker for
       !potential stratosphere grid points further below.
       do iCell=1,nCells
          sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
          if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
          do k=1,nVertLevels
             sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
             if (sgn .GE. 0) candInStrato(k,iCell) = 1
          end do
       end do
       
       !loop over cells and top 5 model levels to seed flood fill with
       !model top that's located above DT (i.e., where |PV| >= pvuVal).
       !can have model top with PV below 2 PVU (e.g., in tropics)
       nChanged = 0
       do iCell=1,nCells
          do k=nVertLevels-5,nVertLevels
             if (candInStrato(k,iCell) .GT. 0) then
                inStrato(k,iCell) = 1
                !candInStrato(k,iCell) = 0
                nChanged = nChanged+1
             end if
          end do
       end do
       
       !flood fill from the determined seeds. since I don't know enough
       !fortran,
       !we'll just brute force a continuing loop rather than queue.
       !here is where the changes to account for domain communication are
       !needed.
       
       call mpas_pool_get_field(diag, 'inStrato', inStrato_f)
       dminfo => inStrato_f % block % domain % dminfo
       global_haloChanged = 1
       
       do while(global_haloChanged .GT. 0) !any cell in a halo has changed, to propagate to other domains
         global_haloChanged = 0 !aggregate the number of changed cells w/in the loop below
         do while(nChanged .GT. 0)
           nChanged = 0
           do iCell=1,nCells !should we look for neighbors of halo cells?
             do k=nVertLevels,1,-1 ! loop over vertical levels from top down
                !update if candidate and neighbor in strato
                if ((candInStrato(k,iCell) .GT. 0) .AND. (inStrato(k,iCell).LT.1) ) then ! modified to match trop routine
                   !nbr above
                   if (k .LT. nVertLevels) then
                     if (inStrato(k+1,iCell) .GT. 0) then
                       inStrato(k,iCell) = 1
                       !candInStrato(k,iCell) = 0 ! commented out to be
                       !consistent with trop routine
                       nChanged = nChanged+1
                       cycle
                     end if
                   end if
                 
                   !side nbrs
                   do iNbr = 1, nEdgesOnCell(iCell)
                     iCellNbr = cellsOnCell(iNbr,iCell)
                     if (inStrato(k,iCellNbr) .GT. 0) then
                       inStrato(k,iCell) = 1
                       !candInStrato(k,iCell) = 0 ! commented out to be
                       !consistent with trop routine
                       nChanged = nChanged+1
                       !exit ! was cycle, but tropspheric loop has exit
                       !here. why?
                       cycle
                     end if
                   end do
                 
                   !nbr below
                   if (k .GT. 1) then
                     if (inStrato(k-1,iCell) .GT. 0) then
                       inStrato(k,iCell) = 1
                       !candInStrato(k,iCell) = 0 ! commented out to be
                       !consistent with trop routine
                       nChanged = nChanged+1
                       cycle
                     end if
                   end if
                     
                end if !candInStrato
             end do !levels
           end do !cells
           global_haloChanged = global_haloChanged+nChanged
         end do !while w/in domain
         
         !communicate to other domains for edge case where a chunk of a
         !block hasn't gotten to fill
         nChanged = global_haloChanged
         call mpas_dmpar_max_int(dminfo, nChanged, global_haloChanged)
         if (global_haloChanged .GT. 0) then !communicate inStrato everywhere
           call mpas_dmpar_exch_halo_field(inStrato_f)
         end if
         nChanged = global_haloChanged !so each block will iterate again if anything changed
       end do !while haloChanged
       deallocate(candInStrato)
       
       
       !Detach high surface PV blobs w/o vertical connection to
       !"stratosphere"
       do iCell=1,nCells
         if (inStrato(1,iCell) .GT. 0) then
           !see how high up we can walk in the column
           do k=2,nVertLevels
             if (inStrato(k,iCell) .LT. 1) then
               exit
             end if !k is highest connected level to sfc
             sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
             if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
             sgn = ertel_pv(k,iCell)*sgnHemi-stratoPV
             if (sgn .LT. 0) then !not actually connected to "stratosphere"
               inStrato(1:k,iCell) = 0
             end if
           end do !k
         end if !inStrato at sfc
       end do !iCell
       
       !Fill iLev_DT with the lowest level above the tropopause (If DT
       !above column, iLev>nVertLevels. If DT below column, iLev=0.
       nChanged = 0
       do iCell=1,nCells
         do k=1,nVertLevels
           if (inStrato(k,iCell) .GT. 0) then
             nChanged = 1
             exit
           end if
         end do !k
         if (nChanged .GT. 0) then !found lowest level
           if (k .EQ. 1) then 
             sgnHemi = sign(1.0_RKIND, latCell(iCell))
             sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
             if (sgn .GT. 0) then !whole column above DT
               iLev_DT(iCell) = 0
             end if
           else
             iLev_DT(iCell) = k
           end if
         else !whole column below DT
           iLev_DT(iCell) = nVertLevels+2
         end if
       end do !iCell
      
    end subroutine floodFill_strato
   
    !*********************************************************************************************************************
    ! Alternative approach: 
    !*********************************************************************************************************************

    subroutine floodFill_tropo(mesh, diag, pvuVal)
       !To find model level of dynamic tropopause: 
       !Simply searching down from TOA within each column to find first
       !instance of 2-PVU surface (i.e., where the PV drops below values
       !characteristic of the stratosphere) 
       !is buggy due to stratospheric wave breaking, which may induce
       !regions of low PV (i.e., PV < 2 PVU) within the stratosphere and
       !thus yield artifically
       !high estimations of the tropopause height. This seems to be more
       !problematic as the mesh gets finer and the vertical vorticity
       !field exhibits greater variability
       !or jumps. 
       !Note that these low-PV anomalies in the stratosphere may persist
       !for long times w/ slow mixing downstream of mountains or deep
       !convection.
       !A few quicker fixes (e.g., make sure PV < 2 PVU for a number of
       !layers; search down from 10 PVU instead of TOA) are hacky and not
       !robust.
       
       !Two flood fill subroutine options are to:
       ! (1) floodFill_strato: flood fill the stratosphere (PV >= 2 PVU)
       ! from stratosphere seeds placed near model top. Strong surface PV
       ! anomalies can connect to 2-PVU region aloft, 
       !     and the resulting "flood-filled 2 PVU" can have sizeable
       !     areas that are located just at/near the surface, while there
       !     is clearly a tropopause above
       !     (i.e., as evident in a vertical cross-section). To address
       !     the large near-surface blobs of PV > 2 PVU, will take the
       !     flood fill mask and try to move upward from
       !     near the surface to 10 PVU within a vertical column. If this
       !     can be done, then the low-level PV anomaly extends to the
       !     stratosphere. Else, remove the stratospheric
       !     designation to disconnect the "surface blob".
       ! (2) floodFill_tropo: flood fill the troposphere (PV < 2 PVU)
       ! from troposphere seeds placed near the surface.
       !
       ! Comparing the two procedures... Somewhat paradoxically, the
       ! bottom of the stratosphere is located lower than the top of the
       ! troposphere.  
      
       !The "output" is iLev_DT, which is the vertical index for the
       !model level just above the dynamic tropopause (i.e., where PV >=
       !pvuVal, which is set below in atm_compute_pv_diagnostics to 2
       !PVU). 
       !If iLev_DT > nVertLevels, then pvuVal is found only above the
       !column (i.e., entire column is in troposphere). If iLev_DT < 1,
       !PV >= pvuVal extends vertically through the entire column 
       !(i.e., the entire column is within the stratosphere).
       !Communication between blocks during the flood fill may be needed
       !to treat some edge cases appropriately.
       
       !Originally, it was assumed that each (MPI) domain would have > 0
       !cells with "right" DT found by flood filling.
       !However, for "small" domains (especially over the poles -- for
       !example, in the Arctic say during winter, when the entire surface
       !can be capped by high PV), 
       !this becomes problematic. So, we need to communicate between
       !domains during the flood fill procedure or else we will find the
       !DT located at/near the surface.
       !The extreme limiting case is if we had every cell as its own
       !domain; then, it's clear that there has to be communication.

       ! MC_TODO: clean up this code
 
       use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
       use mpas_dmpar, only : mpas_dmpar_max_int,mpas_dmpar_exch_halo_field
       use mpas_derived_types, only : dm_info, field2DInteger, field1DInteger
      
       implicit none
       
       type (mpas_pool_type), intent(in) :: mesh
       type (mpas_pool_type), intent(inout) :: diag
       real(kind=RKIND), intent(in) :: pvuVal

       integer :: nbr_count, intCounts, levUse, tropCounts 
       integer :: iCell, k, kk, nChanged, iNbr, iCellNbr, levInd, haloChanged, global_haloChanged
       integer, pointer :: nCells, nVertLevels
       integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
       integer, dimension(:,:), pointer :: cellsOnCell, inTropo, candInTropo, candInStrato
 
       type (field2DInteger), pointer :: inTropo_f
       type (field1DInteger), pointer :: iLev_DT_f
 
       real(kind=RKIND) :: sgnHemi, sgn_pv
       real(kind=RKIND),dimension(:),pointer:: latCell
       real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
       
       type (dm_info), pointer :: dminfo
 
       real(kind=RKIND), dimension(:,:), allocatable :: sgn, sgn2
       integer, dimension(:,:), allocatable :: oppSignPV, interfaceLev
 
       call mpas_pool_get_dimension(mesh, 'nCells', nCells)
       call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
       call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(mesh, 'latCell', latCell)
 
       call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
       call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
       call mpas_pool_get_array(diag, 'inTropo', inTropo)
       call mpas_pool_get_array(diag, 'candInTropo', candInTropo)
       call mpas_pool_get_array(diag, 'candInStrato', candInStrato)       

       allocate(oppSignPV(nVertLevels, nCells+1))
       allocate(sgn(nVertLevels, nCells+1))
       allocate(sgn2(nVertLevels, nCells+1))
       allocate(interfaceLev(nVertLevels, nCells+1))

       candInTropo(:,:) = 0
       candInStrato(:,:) = 0
       inTropo(:,:) = 0
       oppSignPV(:,:) = 0
       interfaceLev(:,:) = 0

       sgn(:,:) = 0.0
       sgn2(:,:) = 0.0

       ! Begin by looping over all cells and vertical levels and flagging cells as troposphere or stratosphere candidates
       do iCell=1,nCells
          sgnHemi = sign(1.0_RKIND, latCell(iCell))                            !at the equator, sign(0)=0
          if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
          do k=1,nVertLevels
             sgn(k,iCell) = ertel_pv(k,iCell)*sgnHemi-pvuVal                   ! quantity will be positive for |PV| > pvuVal
 
             !MC - need to account for pockets of inertial/symmetric instability that develop at upper levels in lee of mountains-- just excluding them as troposphere candidates
             !works well, but this will lead to low estimates of DT height if there's negative PV immediately below where PV drops below pvuVal (probably worse if fewer vertical levels)
             sgn_pv = sign(1.0_RKIND, ertel_pv(k,iCell))     
             if (sgn_pv .EQ. 0.0) sgn_pv = sgnHemi                             ! if sign(PV) = 0, set to sign of latitude
             sgn2(k,iCell) = sgnHemi*sgn_pv                                    ! if sgn2 > 0, PV and latitude are of same sign

             ! Add flags for PV interfaces (i.e., where PV switches from < 2 PVU*sgnHemi to >= 2 PVU*sgnHemi 
             if (k .LT. nVertLevels) then
                if (((sgnHemi .GT. 0) .AND. ((ertel_pv(k+1,iCell).GE.(pvuVal*sgnHemi)) .AND. (ertel_pv(k,iCell).LT.(pvuVal*sgnHemi)))) &
                   .OR. ((sgnHemi .LT. 0) .AND. ((ertel_pv(k+1,iCell).LE.(pvuVal*sgnHemi)) .AND. (ertel_pv(k,iCell).GT.(pvuVal*sgnHemi)))) ) then
                   interfaceLev(k+1,iCell) = 1                                    ! set level above as interface level. these are DT level candidates
                end if
             end if

             ! Assign as either tropo or strato candidates or cells with opposite-sign PV
             if ((sgn2(k,iCell) .GT. 0.0) .AND. (sgn(k,iCell) .LT. 0.0)) then          ! latitude and pv are same sign; |pv| < 2 PVU 
                 candInTropo(k,iCell) = 1

             else if ((sgn2(k,iCell) .GT. 0.0) .AND. (sgn(k,iCell) .GE. 0.0)) then     ! latitude and pv are same sign; but |pv| >= 2 PVU
                 candInStrato(k,iCell) = 1

             else if (sgn2(k,iCell) .LT. 0.0) then                                     ! latitude and pv are opposite sign 
                 oppSignPV(k,iCell) = 1
             end if

             !call mpas_log_write("sgnHemi, ertel_pv, sgn, sgn_pv, sgn2, candInTropo, candInStrato, and oppSignPV are $r, $r, $r, $r, $r, $i, $i, $i", &
             !      realArgs=(/ sgnHemi,ertel_pv(k,iCell),sgn(k,iCell), sgn_pv,sgn2(k,iCell)/), intArgs=(/candInTropo(k,iCell), candInStrato(k,iCell), oppSignPV(k,iCell)/))
         
          end do
       end do
       
       !seed flood fill with near surface that's below DT (can have
       !surface above 2pvu from pv anoms).
       !Note that this would be wrong if low PV "stratospheric" blobs are
       !right above the surface
       nChanged = 0
       levInd = min(nVertLevels, 3)
       do iCell=1,nCells
          lev_loop: do k=1,levInd                                                      !Assign points to troposphere in lowest 3 levels if they're troposphere candidates
             if (candInTropo(k,iCell) .GT. 0) then
                inTropo(k,iCell) = 1
                nChanged = nChanged+1

                !do assignment for cells above any low-level inTropo cells until stratospheric PV values are reached. This is partially redundant, but allows assigning 
                !negative PV cells as trop candidates if above regions of low PV
                if (k .EQ. levInd) then 
                   vert_loop: do kk=levInd+1,nVertLevels
                      if ((candInTropo(kk,iCell) .GT. 0) .AND. ((oppSignPV(kk,iCell) .LT. 1) .AND. (sgn(kk,iCell) .LT. 0.0))) then    ! cell identfied as trop candidate - pv same sign as latitude; pv < 2 PVU 
                         inTropo(kk,iCell) = 1                  
                         nChanged = nChanged+1                  
                         cycle vert_loop

                      else if (oppSignPV(kk,iCell) .GT. 0) then                                                                       ! pv is opposite sign as latitude -- assign as troposphere candidate
                         candInTropo(kk,iCell) = 1
                         cycle vert_loop

                      else if (candInStrato(kk,iCell) .GT. 0) then                                                                    ! if stratosphere candidate is reached while moving upward, exit loop
                         exit vert_loop
                      end if 
                   end do vert_loop 
                else 
                   cycle lev_loop
                end if   

             !MC - add a condition that allows inertially/symmetrically unstable points in the lowest 3 levels to be classified as in the troposphere
             else if (oppSignPV(k,iCell) .GT. 0) then                                  ! if opposite sign PV in lowest 3 levels
                inTropo(k,iCell) = 1
                candInTropo(k,iCell) = 1
                nChanged = nChanged+1

                if (k .EQ. levInd) then                                                ! if opposite sign PV extends up to 3rd level, look at points above
                   vertical_loop: do kk=levInd+1,nVertLevels
                      if (oppSignPV(kk,iCell) .GT. 0) then                             ! pv is opposite sign as latitude -- assign as troposphere candidate
                         candInTropo(kk,iCell) = 1
                         cycle vertical_loop

                      ! Assign regions with same sign PV as latitude as being in troposphere if above negative low-level PV and stratosphere candidate has not yet been reached
                      else if ((candInTropo(kk,iCell) .GT. 0) .AND. ((oppSignPV(kk,iCell) .LT. 1) .AND. (sgn(kk,iCell) .LT. 0.0))) then   ! cell identfied as trop candidate - pv same sign as latitude; pv < 2 PVU    
                         inTropo(kk,iCell) = 1                                                                                           
                         nChanged = nChanged+1                                                                                           
                         cycle vertical_loop

                      else if (candInStrato(kk,iCell) .GT. 0) then                     ! if candidate in stratosphere is reached moving upward, exit loop
                         exit vertical_loop                     
                      end if
                   end do vertical_loop
                end if
             end if
          end do lev_loop
       end do
 
       !flood fill from the given seeds. since I don't know enough
       !fortran,
       !we'll just brute force a continuing loop rather than queue.
       call mpas_pool_get_field(diag, 'inTropo', inTropo_f)
       dminfo => inTropo_f % block % domain % dminfo
       global_haloChanged = 1

       do while(global_haloChanged .GT. 0)                 !any cell in a halo has changed, to propagate to other domains
         global_haloChanged = 0                            !aggregate the number of changed cells w/in the loop below
         do while(nChanged .GT. 0)
           nChanged = 0

           do iCell=1,nCells                               !should we look for neighbors of halo cells?
             vert_loop_halo: do k=1,nVertLevels
                ! for points that are troposphere candidates and not yet assigned to troposphere:
                if ((candInTropo(k,iCell) .GT. 0) .AND. (inTropo(k,iCell) .LT. 1) ) then

                   ! evaluate whether cell below was identified as in the troposphere or not in stratosphere and not opposite sign PV
                   ! neighbor below:
                   if (k .GT. 1) then
                     if ((inTropo(k-1,iCell) .GT. 0) .OR. ((candInStrato(k-1,iCell) .LT. 1) .AND. oppSignPV(k-1,iCell) .LT. 1))  then
                       inTropo(k,iCell) = 1
                       nChanged = nChanged+1
                       cycle vert_loop_halo                 ! if so, move onto the next vertical level
                     end if
                   end if

                   !side neighbors -- original procedure
                   !do iNbr = 1, nEdgesOnCell(iCell)
                     !iCellNbr = cellsOnCell(iNbr,iCell)
                     !if (inTropo(k,iCellNbr) .GT. 0) then
                      !inTropo(k,iCell) = 1
                      !nChanged = nChanged+1
                      !!exit nbr_loop                ! if just exiting, will still do nbr above loop, inflating nChanged count.
                      !cycle vert_loop_halo
                     !end if
                   !end do
 
                   !side neighbors -- made this more stringent than original code, but it may be able to be relaxed due to added constraints in DT identification routine below.
                   nbr_loop: do iNbr = 1, nEdgesOnCell(iCell)
                     iCellNbr = cellsOnCell(iNbr,iCell)
                     if (inTropo(k,iCellNbr) .GT. 0) then
                       !add constraints that for neighbor to lead to trop classification, must be either bounded above or below by points meeting trop classification
                       !if ((k .LT. nVertLevels) .AND. (inTropo(k+1,iCell) .GT. 0)) then         ! try requiring that cell above is already assigned to trop... does this work with iteration?
                       if (k .LT. nVertLevels) then
                          if (inTropo(k+1,iCell) .GT. 0) then     
                             inTropo(k,iCell) = 1
                             nChanged = nChanged+1
                             !exit nbr_loop                                                        ! if just exiting, will still do nbr above loop, inflating nChanged count.
                             cycle vert_loop_halo
                          end if

                       !else if ((k .GT. 1) .AND. (candInStrato(k-1,iCell) .LT. 1)) then         ! if cell below was not assigned to stratosphere (this likely would've been established by neighbor below loop, but just in case...)
                       else if (k .GT. 1) then
                          if (candInStrato(k-1,iCell) .LT. 1) then
                             inTropo(k,iCell) = 1 
                             nChanged = nChanged+1
                             cycle vert_loop_halo
                          end if   
                       end if
                     end if
                   end do nbr_loop

                   !neighbor above
                   if (k .LT. nVertLevels) then
                     if (inTropo(k+1,iCell) .GT. 0) then
                       inTropo(k,iCell) = 1
                       nChanged = nChanged+1
                       cycle vert_loop_halo
                     end if
                   end if
 
                end if !candInTropo
             end do vert_loop_halo 
           end do !cells
           global_haloChanged = global_haloChanged+nChanged
         end do !while w/in domain

         !communicate to other domains for edge case where a chunk of a
         !block hasn't gotten to fill
         nChanged = global_haloChanged
         call mpas_dmpar_max_int(dminfo, nChanged, global_haloChanged)
         if (global_haloChanged .GT. 0) then !communicate inTropo everywhere
           call mpas_dmpar_exch_halo_field(inTropo_f)
         end if
         nChanged = global_haloChanged !so each block will iterate again if anything changed
       end do !while haloChanged
      
       !Moving downward, fill iLev_DT with the lowest level above the tropopause (If DT 
       !above column, iLev>nVertLevels. If DT below column, iLev=0.
       !Note: [previous routine] would find highest tropopause level in scenarios with a double tropopause
       do iCell=1,nCells
         !Keep a tally 
         nChanged = 0
         intCounts = 0     
         tropCounts = 0
         intCounts = COUNT(interfaceLev(:,iCell)==1)                   ! Number of interface levels in vertical column
         tropCounts = COUNT(inTropo(:,iCell)==1)                       ! Number of assigned troposphere cells

         ! First, deal with columns that are entirely in troposphere (i.e., in tropics) because many fail these routines...
         if (tropCounts .EQ. nVertLevels) then 
            nChanged = 1
            levUse = nVertLevels+1

         else
            ! Loop over vertical levels beginning at the top 
            lev_id: do k=nVertLevels,1,-1

              ! If PV interface exists in column, then look for those interfaces in vertical loop
              if (intCounts .GT. 0) then
                 !if DT candidate level identified
                 if (interfaceLev(k,iCell) .GT. 0) then

                    !if (k .GT. 3) then
                       !if ((inTropo(k-1,iCell) .GT. 0) .AND. (candInStrato(k-2,iCell) .LT. 1) .AND. (candInStrato(k-3,iCell) .LT. 1)) then ! one below is assigned as inTropo and two below that are not in stratosphere
                       !   nChanged = 1
                       !   levUse = k
                       !   exit lev_id

                    !Evaluate PV values and vertical continuity around interfaceLev:
                    if (k .GT. 5) then
                       ! if next level below was assigned inTropo and following 4 levels don't have PV characteristic of stratosphere
                       if ((inTropo(k-1,iCell) .GT. 0) .AND. (candInStrato(k-2,iCell) .LT. 1) .AND. (candInStrato(k-3,iCell) .LT. 1) .AND. & 
                          (candInStrato(k-4,iCell) .LT. 1) .AND. (candInStrato(k-5,iCell) .LT. 1))  then
                          nChanged = 1
                          levUse = k
                          exit lev_id

                       ! if next level below was IDed as troposphere candidate, one of following 2 levels is also troposphere candidate, and none of following 4 levels have PV characteristic of stratosphere
                       else if ((candInTropo(k-1,iCell).GT.0) .AND. (candInStrato(k-2,iCell) .LT. 1) .AND. (candInStrato(k-3,iCell) .LT. 1) .AND. &
                               (candInStrato(k-4,iCell) .LT. 1) .AND. (candInStrato(k-5,iCell) .LT. 1) .AND. & 
                               ((candInTropo(k-2,iCell).GT.0) .OR. (candInTropo(k-3,iCell).GT.0))) then 

                          nChanged = 1
                          levUse = k
                          exit lev_id

                       !if 5 consecutive levels beneath interface don't have characteristics of stratosphere and at least one of three levels beneath interface has characteristics of troposphere 
                       else if ((candInStrato(k-1,iCell) .LT. 1) .AND. (candInStrato(k-2,iCell) .LT. 1) .AND. (candInStrato(k-3,iCell) .LT. 1) .AND. &
                          (candInStrato(k-4,iCell) .LT. 1) .AND. (candInStrato(k-5,iCell) .LT. 1) .AND. ((candInTropo(k-1,iCell) .GT. 0) .OR. &
                          (candInTropo(k-2,iCell) .GT. 0) .OR. (candInTropo(k-3,iCell) .GT. 0))) then 
                          nChanged = 1
                          levUse = k
                          exit lev_id

                       end if
                    end if  ! k > 5

                    !otherwise, if no stratospheric characterstics for at least 8 consecutive levels
                    if (k .GT. 8) then 
                       if ((candInStrato(k-1,iCell) .LT. 1) .AND. (candInStrato(k-2,iCell) .LT. 1) .AND. (candInStrato(k-3,iCell) .LT. 1) .AND. &
                           (candInStrato(k-4,iCell) .LT. 1) .AND. (candInStrato(k-5,iCell) .LT. 1)  .AND. (candInStrato(k-6,iCell) .LT. 1) .AND. &
                           (candInStrato(k-7,iCell) .LT. 1) .AND. (candInStrato(k-8,iCell) .LT. 1) ) then
                           nChanged = 1
                           levUse = k
                           exit lev_id
                        end if
                    end if

                    !account for interfaces near the surface (e.g., in hurricanes)
                    if ((k .LE. 5) .AND. (k .GT. 1)) then
                       ! if next level below was assigned inTropo, candInTropo, or oppSignPV
                       if ((inTropo(k-1,iCell) .GT. 0) .OR. (candInTropo(k-1,iCell) .GT. 0) .OR. (oppSignPV(k-1,iCell) .GT. 0))  then
                          nChanged = 1
                          levUse = k
                          exit lev_id
                       end if
                    end if
 
                 end if ! interfaceLev > 0
    
              ! no interfaces identified in column 
              !else if (intCounts .EQ. 0) then

                 !Search for inTropo flags (this should be null -- addressed before do loop)
              !   if (inTropo(k,iCell) .GT. 0) then
              !      nChanged = 1
              !      levUse = k+1
              !      call mpas_log_write("is this ever called?")
              !      exit lev_id
       
                 !If make it through loop without identifying tropopause and latitude <= +- 2 degrees, set tropopause to top of column
                 !else if ((k .EQ. 1) .AND. (abs(latCell(iCell) .LE. 0.035))) then
                 !   nChanged = 1
                 !   levUse = nVertLevels+1
                 !   exit lev_id
            
              !   end if ! tropo flags

              end if ! intCounts       
            end do lev_id
          end if

         if (nChanged .GT. 0) then      !found troposphere's highest level
           iLev_DT(iCell) = levUse      !level above troposphere (>nVertLevels if whole column below 2pvu; e.g., tropics)
         else                           !whole column above DT (e.g., arctic pv tower)
           iLev_DT(iCell) = 0
         end if
 
         ! Add in correction for locations near the equator with iLev_DT assigned to 0
         if ((iLev_DT(iCell) .EQ. 0) .AND. (latCell(iCell) .LE. 0.0436)) then ! corresponds to 2.5 deg latitude
               iLev_DT(iCell) = nVertLevels+1
         end if
       end do !iCell

       ! Do correction pass to change cells that differ significantly from surrounding cells
       !call mpas_log_write("beginning correction loop")
       call mpas_pool_get_field(diag, 'iLev_DT', iLev_DT_f)
       dminfo => iLev_DT_f % block % domain % dminfo
       call mpas_dmpar_exch_halo_field(iLev_DT_f)

       cells: do iCell=1,nCells
          intCounts = 0

          nbrloop: do iNbr = 1,nEdgesOnCell(iCell)
            iCellNbr = cellsOnCell(iNbr,iCell)

             if (abs(iLev_DT(iCell)-iLev_DT(iCellNbr)) .GT. 5) then
               intCounts = intCounts+1
               cycle nbrloop
             end if
          end do nbrloop

          if (intCounts .GE. (nEdgesOnCell(iCell)-2)) then
            ! Loop through neighbors again. Need to evaluate whether neighboring DT values are in interfaceLev array.
            nbrloop2: do iNbr = 1,nEdgesOnCell(iCell)
              iCellNbr = cellsOnCell(iNbr,iCell)
              
              !if ((abs(iLev_DT(iCell)-iLev_DT(iCellNbr)) .GT. 5) .AND. (interfaceLev(iLev_DT(iCellNbr),iCell) .GT. 0)) then
              if ((abs(iLev_DT(iCell)-iLev_DT(iCellNbr)) .GT. 5) .AND. (iLev_DT(iCellNbr) .GT. 0) .AND. (iLev_DT(iCellNbr) .LE. nVertLevels)) then
                if (interfaceLev(iLev_DT(iCellNbr),iCell) .GT. 0) then       
                   iLev_DT(iCell) = iLev_DT(iCellNbr)
                end if   
              else
                cycle nbrloop2
              end if
            end do nbrloop2

          else
            cycle cells
          end if

       end do cells

       deallocate(interfaceLev)
       deallocate(oppSignPV)
       deallocate(sgn)
       deallocate(sgn2)

    end subroutine floodFill_tropo

    !*********************************************************************************************************************
    ! NS: Subroutine to compute various fields on 2-PVU surface using the calculated PV field - potential temperature, 
    !     uZonal, uMeridional, vertical vorticity
    ! MC: Modified interpolation of vorticity to cell centers procedure to be consistent with other changes  
    !     This routine should use PV field at end of time step! 
    !********************************************************************************************************************* 
  
    subroutine interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
 
       use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
       use mpas_constants, only: r_earth=>a
       
       IMPLICIT NONE
       
       type (mpas_pool_type), intent(in) :: mesh
       type (mpas_pool_type), intent(inout) :: diag
       real(kind=RKIND) :: pvuVal, missingVal
       
       integer :: iCell, k
       integer, pointer :: nCells, nVertLevels
       integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
       integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex 
                                           
       real(kind=RKIND), dimension(:), pointer:: areaCell, latCell, u_pv, v_pv, theta_pv, vort_pv, pres_pv, height_pv
       real(kind=RKIND), dimension(:,:), pointer:: uReconstructZonal, uReconstructMeridional, vorticity, theta, ertel_pv, &
                                                   kiteAreasOnVertex, pressure, zgrid
       real(kind=RKIND), dimension(:,:), allocatable :: vVort, zCell, zCell_geo
       
       call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
       call mpas_pool_get_dimension(mesh, 'nCells', nCells)
       call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
       call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
       call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
       call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
       call mpas_pool_get_array(mesh, 'areaCell', areaCell)
       call mpas_pool_get_array(mesh, 'latCell', latCell)
       call mpas_pool_get_array(mesh, 'zgrid', zgrid)      
 
       call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
       call mpas_pool_get_array(diag, 'theta', theta)
       call mpas_pool_get_array(diag, 'pressure', pressure)
       call mpas_pool_get_array(diag, 'vorticity', vorticity)
       call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
       call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
       call mpas_pool_get_array(diag, 'u_pv', u_pv)
       call mpas_pool_get_array(diag, 'v_pv', v_pv)
       call mpas_pool_get_array(diag, 'theta_pv', theta_pv)
       call mpas_pool_get_array(diag, 'vort_pv', vort_pv)
       call mpas_pool_get_array(diag, 'pres_pv', pres_pv)
       call mpas_pool_get_array(diag, 'height_pv', height_pv)
       call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)

       allocate(vVort(nVertLevels,nCells))
       allocate(zCell(nVertLevels,nCells)) 
       allocate(zCell_geo(nVertLevels,nCells)) 
     
       ! Interpolate horizontal winds to pvuVal isosurface 
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, uReconstructZonal, u_pv, missingVal, iLev_DT)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, uReconstructMeridional, v_pv, missingVal, iLev_DT)

       ! Interpolate theta
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, theta, theta_pv, missingVal, iLev_DT)

       ! Interpolate height of cell center
       call interp_wLev_thetaLev(zgrid, nCells, nVertLevels, zCell)
       
       ! convert geometric to geopotential height: 
       zCell_geo = (zCell * r_earth)/(zCell + r_earth)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, zCell_geo, height_pv, missingVal, iLev_DT)

       ! Interpolate pressure
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, pressure, pres_pv, missingVal, iLev_DT)

       ! Interpolate absolute vertical vorticity (could just use pv_vertex?) 
       call interp_absVertVort(vorticity, nCells, nEdgesOnCell, verticesOnCell, &
                               cellsOnVertex, areaCell, kiteAreasOnVertex, vVort)
               
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, ertel_pv, vVort, &
                      vort_pv, missingVal, iLev_DT)

       deallocate(vVort)
       deallocate(zCell)
       deallocate(zCell_geo)

    end subroutine interp_pv_diagnostics     

    !*********************************************************************************************************************
    ! NS: Subroutine to compute various tendency fields on 2-PVU surface using the calculated PV field
    ! MC: This routine should use PV field and dynamic tropopause from beginning of time step, so this has been modified 
    !     accordingly. 
    !********************************************************************************************************************* 
 
    subroutine interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)
 
       use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
       
       IMPLICIT NONE
       
       type (mpas_pool_type), intent(in) :: mesh
       type (mpas_pool_type), intent(inout) :: diag
       real(kind=RKIND) :: pvuVal, missingVal
       
       integer :: iCell, k
       integer, pointer :: nCells, nVertLevels
       integer, dimension(:), pointer :: iLev_DT_prev
                                           
       real(kind=RKIND), dimension(:),pointer :: latCell, depv_dt_diab_pv, depv_dt_fric_pv, depv_dt_dyn_pv
       real(kind=RKIND), dimension(:,:),pointer :: depv_dt_diab, depv_dt_fric, depv_dt_dyn, &
                                                   ertel_pv_prev                                 ! MC changed 
      
       call mpas_pool_get_dimension(mesh, 'nCells', nCells) 
       call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
       
       call mpas_pool_get_array(mesh, 'latCell', latCell)
       call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev)  
       call mpas_pool_get_array(diag, 'iLev_DT_prev', iLev_DT_prev)                               
       call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
       call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
       call mpas_pool_get_array(diag, 'depv_dt_dyn', depv_dt_dyn)
       
       call mpas_pool_get_array(diag, 'depv_dt_diab_pv', depv_dt_diab_pv)
       call mpas_pool_get_array(diag, 'depv_dt_fric_pv', depv_dt_fric_pv)
       call mpas_pool_get_array(diag, 'depv_dt_dyn_pv', depv_dt_dyn_pv)
       
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv_prev, depv_dt_diab, depv_dt_diab_pv, missingVal, iLev_DT_prev)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv_prev, depv_dt_fric, depv_dt_fric_pv, missingVal, iLev_DT_prev)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv_prev, depv_dt_dyn, depv_dt_dyn_pv, missingVal, iLev_DT_prev)

    end subroutine interp_pvBudget_diagnostics
  
    !*********************************************************************************************************************
    ! NS: Subroutine to linearly interpolate columns of field1 to where field0 is interpVal*sign(lat) using level above 
    !     tropopause already diagnosed 
    !********************************************************************************************************************* 

    subroutine interp_pv(nCells, nLevels, interpVal, latCell, field0, &
                           field1, field_interp, missingVal, iLev_DT)
        
       IMPLICIT NONE  
       
       integer :: nCells, nLevels
       integer, intent(in) :: iLev_DT(nCells)
       real(kind=RKIND) :: interpVal, missingVal
       real(kind=RKIND), intent(in) :: latCell(nCells)
       real(kind=RKIND), intent(in) :: field0(nLevels,nCells), field1(nLevels,nCells)
       real(kind=RKIND), intent(out) :: field_interp(nCells)
 
       !  local variables
       integer :: iCell, iLev, levInd, indlNbr
       real(kind=RKIND) :: valh, vall, vallNbr, sgnh, sgnl, sgnlNbr
       real(kind=RKIND) :: dv_dl, levFrac, valInterpCell, sgnHemi
 
       do iCell = 1,nCells
         !starting from top, trap val if values on opposite side
         levInd = -1 ! what should happen with missing values?
         levFrac = 0.0
         sgnHemi = sign(1.0_RKIND, latCell(iCell)) !problem at the equator...is sign(0)=0?
         if (sgnHemi .EQ. 0.0) sgnHemi = 1.0
         valInterpCell = interpVal*sgnHemi
         
         iLev = iLev_DT(iCell)          ! lowest vertical level above the tropopause
         if (iLev .GT. nLevels) then    ! if no identified tropopause in column and all column in troposphere
           levInd = -1
           sgnl = -1.0
         else if (iLev .LT. 1) then     ! if no identified tropopause in column and all column in stratosphere
           levInd = -1
           sgnl = 1.0
         else
           valh = field0(iLev,iCell)      ! value at the level just above tropopause
           vall = field0(iLev-1,iCell)    ! value at level just below tropopause

           ! MC: need to ensure that 2 PVU is actually between valh and vall before proceeding to avoid huge erroneous interpolated vals    
           !if ((abs(valInterpCell) .LE. abs(valh)) .AND. (abs(valInterpCell) .GE. abs(vall))) then
           ! below should be more robust in situations where PV changes sign across tropopause
           if (((sgnHemi .GT. 0)  .AND. ((valInterpCell .LE. valh) .AND. (valInterpCell .GE. vall))) &
            .OR. ((sgnHemi .LT. 0)  .AND. ((valInterpCell .GE. valh) .AND. (valInterpCell .LE. vall)))) then

              !sandwiched value. equal in case val0 is a vals[l].
              !get linear interpolation: val0 = vals[l]+dvals/dl * dl
              !Avoid divide by 0 by just assuming value is 
              !halfway between...
              dv_dl = valh-vall;             ! change in PV across vertical levels
              if (abs(dv_dl)<1.e-6) then     ! if difference between PV values is tiny, set levFrac = 0.5
                 levFrac = 0.5;
              else
                 levFrac = (valInterpCell-vall)/dv_dl     ! if not tiny, calculate levFrac as difference between interp reference value and PV value just below level
              end if                                         ! of tropopause / change in PV across vertical levels 
              levInd = iLev-1                                ! index is level just below tropopause

           ! MC: need to set these incorrectly identified DT points to something... 
           else
              levInd = -1
              sgnl = 0.0
           end if   ! bounding 2 PVU
         end if !iLev in column
 
         !find value of field using index we just found
         if (levInd < 0) then                                 !didn't trap value
           if (sgnl > 0.0) then                               !column above value, take value at the lowest model level 
             field_interp(iCell) = field1(1,iCell)

           else if (sgnl < 0.0) then                          !column below value, take value at highest model level
             !field(iCell) = missingVal
             field_interp(iCell) = field1(nLevels,iCell)

           else
             field_interp(iCell) = missingVal                 ! MC: set to missing if DT incorrectly identified
           end if

         else
           valh = field1(levInd+1,iCell)    ! value of field we're interpolating at level above tropopause
           vall = field1(levInd,iCell)      ! value of field at level below tropopause

           dv_dl = valh-vall                ! change in field across vertical levels
           field_interp(iCell) = vall+dv_dl*levFrac    ! interpolated value = value below tropopause + change in value across vertical level * levFrac
         end if

       end do
       
    end subroutine interp_pv
   
    !*********************************************************************************************************************
    ! MC: Subroutine to calculate the dot product between two 3D vectors 
    !*********************************************************************************************************************

    subroutine calc_dotProduct_3D(vec1, vec2, nCells, nVertLevels, dotResult)

       IMPLICIT NONE
       
       integer, intent(in) :: nVertLevels, nCells
       real(kind=RKIND), dimension(:,:,:), intent(in) :: vec1, vec2
       real(kind=RKIND), dimension(:,:),  intent(out) :: dotResult
       integer :: iCell, k

       dotResult(:,:) = 0.0_RKIND

       do iCell=1,nCells
          do k=1,nVertLevels
             dotResult(k,iCell) = vec1(k,iCell,1)*vec2(k,iCell,1) + vec1(k,iCell,2)*vec2(k,iCell,2) + vec1(k,iCell,3)*vec2(k,iCell,3)
          end do
       end do

    end subroutine calc_dotProduct_3D

    !*********************************************************************************************************************
    ! MC: Subroutine to interpolate the absolute vertical vorticity to cell centers from the absolute 
    !     vertical vorticity at the vertices (pv_vertex) 
    !     Note: pv_vertex variable was defined incorrectly in MPAS registry -- there is no density dependence
    !********************************************************************************************************************* 

    subroutine interp_absVertVort(pv_vertex, nCells, nEdgesOnCell, verticesOnCell, &
                                cellsOnVertex, areaCell, kiteAreasOnVertex, absVort)
 
       IMPLICIT NONE
       
       integer, intent(in) :: nCells
       integer, dimension(:), intent(in) :: nEdgesOnCell
       integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
       real(kind=RKIND), dimension(:), intent(in) :: areaCell
       real(kind=RKIND), dimension(:,:), intent(in) :: pv_vertex, kiteAreasOnVertex
       real(kind=RKIND), dimension(:,:), intent(out) :: absVort
       integer :: i, j, cellIndOnVertex, iVertex
 
       absVort(:,:) = 0.0_RKIND
 
       do i=1,nCells
          do j=1,nEdgesOnCell(i)
             iVertex = verticesOnCell(j,i)
             cellIndOnVertex = FINDLOC(cellsOnVertex(:,iVertex),VALUE=i,DIM=1)
             absVort(:,i) = absVort(:,i) + kiteAreasOnVertex(cellIndOnVertex,iVertex) * pv_vertex(:,iVertex)
          end do
          absVort(:,i) = absVort(:,i) / areaCell(i)
       end do
 
    end subroutine interp_absVertVort

    !*********************************************************************************************************************
    ! MW: Subroutine that can be used to recompute the absolute vorticity at cell vertices
    !*********************************************************************************************************************
    subroutine recompute_absVort_vertex(u, nVertices, nVertLevels, vertexDegree, invAreaTriangle, &
                                 dcEdge, edgesOnVertex, edgesOnVertex_sign, fVertex, vort )

       IMPLICIT NONE
       
       integer, intent(in) :: nVertices, vertexDegree, nVertLevels
       integer, dimension(:,:), intent(in) :: edgesOnVertex
       real(kind=RKIND), dimension(:), intent(in) :: invAreaTriangle, fVertex, dcEdge
       real(kind=RKIND), dimension(:,:), intent(in) :: u, edgesOnVertex_sign
       real(kind=RKIND), dimension(:,:), intent(inout) :: vort

       ! local variables
       integer :: iVertex, iEdge, i, k
       real (kind=RKIND) :: s

       do iVertex=1,nVertices
         vort(1:nVertLevels,iVertex) = 0.0
         do i=1,vertexDegree
            iEdge = edgesOnVertex(i,iVertex)
            s = edgesOnVertex_sign(i,iVertex) * dcEdge(iEdge)
            do k=1,nVertLevels
               vort(k,iVertex) = vort(k,iVertex) + s * u(k,iEdge)
            end do
         end do
         do k=1,nVertLevels
            vort(k,iVertex) = vort(k,iVertex) * invAreaTriangle(iVertex)
            vort(k,iVertex) = vort(k,iVertex) + fVertex(iVertex)
         end do
      end do

    end subroutine recompute_absVort_vertex

    !*********************************************************************************************************************    
    ! MC: Subroutine to calculate the horizontal gradient of a field on the cell edges using field values at the 
    !     adjacent cell centers as:
    !
    ! varGrad(edgeUse,kLev) = (cellVar(cellsOnEdge(edgeUse,2),kLev)-cellVar(cellsOnEdge(edgeUse,1),kLev))/dcEdge(edgeUse)
    !
    !     and then assign the correct sign based on its direction (i.e., into or out of the cell), the value of 
    !     edgesOnCell_sign, and the convention for the u winds: "Positive u (normal) velocity is always defined as 
    !     flow from cellsOnEdge(1,jEdge) to cellsOnEdge(2,jEdge) for edge iEdge" (MPAS tutorial 2019).
    !
    !     The expression for calculating the gradient on each edge comes from Eq. 22 in Ringler et al. (2010)
    !*********************************************************************************************************************    
 
    subroutine calc_gradOnEdges(cellVar, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, varGrad)
 
       IMPLICIT NONE
       
       integer, intent(in) :: nCells, nEdges, nVertLevels
       integer, dimension(:), intent(in) :: nEdgesOnCell
       integer, dimension(:,:), intent(in) :: cellsOnEdge, edgesOnCell
       real(kind=RKIND), dimension(:), intent(in) :: dcEdge
       real(kind=RKIND), dimension(:,:), intent(in) :: cellVar, edgesOnCell_sign
       real(kind=RKIND), dimension(:,:), intent(out) :: varGrad
       integer :: iCell, jEdge, kLev, edgeSign, edgeUse, index_j1, index_j2, sign_j1, sign_j2
 
       varGrad(:,:) = 0.0_RKIND
 
       cell_loop: do iCell=1,nCells
         edge_loop: do jEdge=1,nEdgesOnCell(iCell)
            lev_loop: do kLev=1,nVertLevels
 
               ! Edges and edge signs for jEdge along parent iCell
               edgeSign = edgesOnCell_sign(jEdge,iCell)
               edgeUse = edgesOnCell(jEdge,iCell) 
 
               ! The indices of edgeUse likely differ in edgesOnCell array for each
               ! cell. Need to find the correct indices and the sign of the normal
               ! vector for each edge in edgesOnCell
               !
               ! -- if sign_j1 > 0, normal vector points out of cellsOnEdge(edgeUse,1) 
               ! -- if sign_j2 > 0, normal vector points out of cellsOnEdge(edgeUse,2) 
 
               index_j1 = FINDLOC(edgesOnCell(:,cellsOnEdge(1,edgeUse)),VALUE=edgeUse, DIM=1)
               index_j2 = FINDLOC(edgesOnCell(:,cellsOnEdge(2,edgeUse)),VALUE=edgeUse, DIM=1)
               sign_j1 = edgesOnCell_sign(index_j1,cellsOnEdge(1,edgeUse))
               sign_j2 = edgesOnCell_sign(index_j2,cellsOnEdge(2,edgeUse))
 
               ! Calculate gradient of field by taking the difference of the values
               ! at the adjacent cell centers divided by the distance between the
               ! cells
 
               varGrad(kLev,edgeUse) = cellVar(kLev,cellsOnEdge(2,edgeUse)) - cellVar(kLev,cellsOnEdge(1,edgeUse))
               varGrad(kLev,edgeUse) = varGrad(kLev,edgeUse)/dcEdge(edgeUse)
 
               ! Ensure that the sign of the gradient is consistent with the
               ! convention for the u (normal winds). Note: I think the signs are
               ! correct without doing this procedure, but I will keep it here just
               ! in case. 
     
               IF (varGrad(kLev,edgeUse) .gt. 0) THEN
                  ! Gradient vector points toward cellsOnEdge(edgeUse,2) -> 
                  ! should be directed inward for cellsOnEdge(edgeUse,2)
 
                  ! What is sign of normal vector along edgeUse for each cell? 
                  IF (sign_j2 .lt. 0) THEN       ! Normal vector points inward for cellsOnEdge(edgeUse,2) 
                                                 ! and outward for cellsOnEdge(edgeUse,1)
                      varGrad(kLev,edgeUse) = ABS(varGrad(kLev,edgeUse))
                  ELSE
                      varGrad(kLev,edgeUse) = -ABS(varGrad(kLev,edgeUse))
                  END IF
 
               ELSE IF (varGrad(kLev,edgeUse) .lt. 0) THEN
                  ! Gradient vector points toward cellsOnEdge(edgeUse,1) -> 
                  ! should be directed inward for cellsOnEdge(edgeUse,1)
 
                  ! What is sign of normal vector along edgeUse for each cell? 
                  IF (sign_j1 .lt. 0) THEN      ! Normal vector points inward for cellsOnEdge(edgeUse,1) 
                                                ! and outward for cellsOnEdge(edgeUse,2)
                      varGrad(kLev,edgeUse) = ABS(varGrad(kLev,edgeUse))
                  ELSE
                      varGrad(kLev,edgeUse) = -ABS(varGrad(kLev,edgeUse))
                  END IF
 
               END IF
 
            end do lev_loop
          end do edge_loop
       end do cell_loop
 
    end subroutine calc_gradOnEdges

    !*********************************************************************************************************************
    ! MC: Subroutine takes gradient field valid on cell edges and reconstructs the horizontal gradient vectors at the cell 
    !     center in a manner analogous to the u reconstruction of mpas_reconstruct_2d in mpas_vector_reconstruction.F
    !*********************************************************************************************************************
 
    subroutine mpas_reconstruct_grad(gradEdge, latCell, lonCell, coeffs_reconstruct, nCells, nVertLevels, &
                                     edgesOnCell, nEdgesOnCell, &
                                     gradReconstructZonal, gradReconstructMeridional)
 
       IMPLICIT NONE
 
       integer, intent(in) :: nVertLevels, nCells
       integer, dimension(:,:), intent(in) :: edgesOnCell
       integer, dimension(:), intent(in) :: nEdgesOnCell
 
       real(kind=RKIND), dimension(:),   intent(in)   :: latCell, lonCell
       real(kind=RKIND), dimension(:,:), intent(in)   :: gradEdge
       real(kind=RKIND), dimension(:,:,:), intent(in) :: coeffs_reconstruct
       real(kind=RKIND), dimension(:,:), intent(out)  :: gradReconstructZonal, gradReconstructMeridional
 
       ! local variables
       integer :: iCell, jEdge, edgeUse, kLev
       real(kind=RKIND) :: clat, slat, clon, slon
       real(kind=RKIND), dimension(:,:), allocatable :: gradReconstructX, gradReconstructY, gradReconstructZ
 
       allocate(gradReconstructX(nVertLevels,nCells))
       allocate(gradReconstructY(nVertLevels,nCells))
       allocate(gradReconstructZ(nVertLevels,nCells))
 
       gradReconstructX(nVertLevels,nCells) = 0.0_RKIND
       gradReconstructY(nVertLevels,nCells) = 0.0_RKIND
       gradReconstructZ(nVertLevels,nCells) = 0.0_RKIND
       gradReconstructZonal(nVertLevels,nCells) = 0.0_RKIND
       gradReconstructMeridional(nVertLevels,nCells) = 0.0_RKIND
 
       cell_loop: do iCell=1,nCells
          edge_loop: do jEdge=1,nEdgesOnCell(iCell)
 
             edgeUse = edgesOnCell(jEdge,iCell)
 
             gradReconstructX(:,iCell) = gradReconstructX(:,iCell) &
               + coeffs_reconstruct(1,jEdge,iCell) * gradEdge(:,edgeUse)
             gradReconstructY(:,iCell) = gradReconstructY(:,iCell) &
               + coeffs_reconstruct(2,jEdge,iCell) * gradEdge(:,edgeUse)
             gradReconstructZ(:,iCell) = gradReconstructZ(:,iCell) &
               + coeffs_reconstruct(3,jEdge,iCell) * gradEdge(:,edgeUse)
 
           end do edge_loop
 
           clat = COS(latCell(iCell))
           slat = SIN(latCell(iCell))
           clon = COS(lonCell(iCell))
           slon = SIN(lonCell(iCell))
 
           gradReconstructZonal(:,iCell) = -gradReconstructX(:,iCell)*slon + &
                                           gradReconstructY(:,iCell)*clon
 
 
           gradReconstructMeridional(:,iCell) = -(gradReconstructX(:,iCell)*clon + &
                                               gradReconstructY(:,iCell)*slon)*slat + &
                                               gradReconstructZ(:,iCell)*clat
 
       end do cell_loop
 
       deallocate(gradReconstructX)
       deallocate(gradReconstructY)
       deallocate(gradReconstructZ)
 
    end subroutine mpas_reconstruct_grad

    !*********************************************************************************************************************
    ! MC: Combined subroutines to calculate the horizontal gradient of a field on the cell edges using values at 
    !     the adjacent cell centers 
    !
    ! varGrad(edgeUse,kLev) = (cellVar(cellsOnEdge(edgeUse,2),kLev) - cellVar(cellsOnEdge(edgeUse,1),kLev)) / dcEdge(edgeUse)
    !
    !     and then assign the correct sign based on its direction (i.e., into or out of the cell), the value of 
    !     edgesOnCell_sign, and the convention for the u winds: "Positive u (normal) velocity is always defined as 
    !     flow from cellsOnEdge(1,jEdge) to cellsOnEdge(2,jEdge) for edge iEdge" (MPAS tutorial 2019).
    !
    !     The expression for calculating the gradient on each edge comes from Eq. 22 in Ringler et al. (2010)
    !
    !     Following the gradient on edge calculation, the gradient is then reconstructed to the cell centers in a 
    !     manner analogous to the u reconstruction of mpas_reconstruct_2d in mpas_vector_reconstruction.F
    !
    !     The purpose of combining these into one subroutine is to reduce the number of stored intermediate variables, 
    !     which are needed for the halo communication to work properly
    !
    !    MC 02/2024 -- note this will crash when running with DEBUG=true if using nCells and not nCellsSolve in the 
    !    calc_epv and diagnostics subroutines
    !*********************************************************************************************************************
 
    subroutine calc_gradOnEdges_reconCellCenter(cellVar, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, & 
                                latCell, lonCell, gradReconstructZonal, gradReconstructMeridional)

       IMPLICIT NONE

       integer, intent(in) :: nCells, nEdges, nVertLevels
       integer, dimension(:), intent(in) :: nEdgesOnCell
       integer, dimension(:,:), intent(in) :: cellsOnEdge, edgesOnCell
 
       real(kind=RKIND), dimension(:), intent(in) :: dcEdge, latCell, lonCell
       real(kind=RKIND), dimension(:,:), intent(in) :: cellVar, edgesOnCell_sign
       real(kind=RKIND), dimension(:,:,:), intent(in) :: coeffs_reconstruct
       real(kind=RKIND), dimension(:,:), intent(out)  :: gradReconstructZonal, gradReconstructMeridional

       ! local variables
       real(kind=RKIND) :: clat, slat, clon, slon
       real(kind=RKIND), dimension(:,:), allocatable :: varGrad, gradReconstructX, gradReconstructY, gradReconstructZ
       integer :: iCell, jEdge, kLev, edgeSign, edgeUse, index_j1, index_j2, sign_j1, sign_j2
 
       allocate(varGrad(nVertLevels,nEdges))
       allocate(gradReconstructX(nVertLevels,nCells))
       allocate(gradReconstructY(nVertLevels,nCells))
       allocate(gradReconstructZ(nVertLevels,nCells))
 
       varGrad(:,:) = 0.0_RKIND
       gradReconstructX(:,:) = 0.0_RKIND
       gradReconstructY(:,:) = 0.0_RKIND
       gradReconstructZ(:,:) = 0.0_RKIND
       gradReconstructZonal(:,:) = 0.0_RKIND
       gradReconstructMeridional(:,:) = 0.0_RKIND
 
       cell_loop: do iCell=1,nCells
         edge_loop: do jEdge=1,nEdgesOnCell(iCell)

            ! Edges and edge signs for jEdge along parent iCell
            edgeSign = edgesOnCell_sign(jEdge,iCell)
            edgeUse = edgesOnCell(jEdge,iCell)

            ! The indices of edgeUse likely differ in edgesOnCell array for each
            ! cell. Need to find the correct indices and the sign of the normal
            ! vector for each edge in edgesOnCell
            !
            ! -- if sign_j1 > 0, normal vector points out of cellsOnEdge(edgeUse,1) 
            ! -- if sign_j2 > 0, normal vector points out of cellsOnEdge(edgeUse,2) 

            index_j1 = FINDLOC(edgesOnCell(:,cellsOnEdge(1,edgeUse)),VALUE=edgeUse, DIM=1)
            index_j2 = FINDLOC(edgesOnCell(:,cellsOnEdge(2,edgeUse)),VALUE=edgeUse, DIM=1)
            sign_j1 = edgesOnCell_sign(index_j1,cellsOnEdge(1,edgeUse))
            sign_j2 = edgesOnCell_sign(index_j2,cellsOnEdge(2,edgeUse))

            lev_loop: do kLev=1,nVertLevels

               ! Calculate gradient of field by taking the difference of the values
               ! at the adjacent cell centers divided by the distance between the
               ! cells
 
               varGrad(kLev,edgeUse) = cellVar(kLev,cellsOnEdge(2,edgeUse)) - cellVar(kLev,cellsOnEdge(1,edgeUse))
               varGrad(kLev,edgeUse) = varGrad(kLev,edgeUse)/dcEdge(edgeUse)
 
               ! Ensure that the sign of the gradient is consistent with the
               ! convention for the u (normal winds). Note: I think the signs are
               ! correct without doing this procedure, but I will keep it here just
               ! in case. 
 
               IF (varGrad(kLev,edgeUse) .gt. 0) THEN
                  ! Gradient vector points toward cellsOnEdge(edgeUse,2) -> 
                  ! should be directed inward for cellsOnEdge(edgeUse,2)
 
                  ! What is sign of normal vector along edgeUse for each cell? 
                  IF (sign_j2 .lt. 0) THEN       ! Normal vector points inward for cellsOnEdge(edgeUse,2) 
                                                 ! and outward for cellsOnEdge(edgeUse,1)
                      varGrad(kLev,edgeUse) = ABS(varGrad(kLev,edgeUse))
                  ELSE
                      varGrad(kLev,edgeUse) = -ABS(varGrad(kLev,edgeUse))
                  END IF
 
               ELSE IF (varGrad(kLev,edgeUse) .lt. 0) THEN
                  ! Gradient vector points toward cellsOnEdge(edgeUse,1) -> 
                  ! should be directed inward for cellsOnEdge(edgeUse,1)
 
                  ! What is sign of normal vector along edgeUse for each cell? 
                  IF (sign_j1 .lt. 0) THEN      ! Normal vector points inward for cellsOnEdge(edgeUse,1) 
                                                ! and outward for cellsOnEdge(edgeUse,2)
                      varGrad(kLev,edgeUse) = ABS(varGrad(kLev,edgeUse))
                  ELSE
                      varGrad(kLev,edgeUse) = -ABS(varGrad(kLev,edgeUse))
                  END IF
 
               END IF
 
            end do lev_loop
 
            gradReconstructX(:,iCell) = gradReconstructX(:,iCell) &
              + coeffs_reconstruct(1,jEdge,iCell) * varGrad(:,edgeUse)
            gradReconstructY(:,iCell) = gradReconstructY(:,iCell) &
              + coeffs_reconstruct(2,jEdge,iCell) * varGrad(:,edgeUse)
            gradReconstructZ(:,iCell) = gradReconstructZ(:,iCell) &
              + coeffs_reconstruct(3,jEdge,iCell) * varGrad(:,edgeUse)
 
          end do edge_loop
 
          clat = COS(latCell(iCell))
          slat = SIN(latCell(iCell))
          clon = COS(lonCell(iCell))
          slon = SIN(lonCell(iCell))
 
          gradReconstructZonal(:,iCell) = -gradReconstructX(:,iCell)*slon + &
                                           gradReconstructY(:,iCell)*clon
 
 
          gradReconstructMeridional(:,iCell) = -(gradReconstructX(:,iCell)*clon + &
                                                gradReconstructY(:,iCell)*slon)*slat + &
                                                gradReconstructZ(:,iCell)*clat
 
       end do cell_loop
 
       deallocate(gradReconstructX)
       deallocate(gradReconstructY)
       deallocate(gradReconstructZ)
       deallocate(varGrad)
 
    end subroutine calc_gradOnEdges_reconCellCenter

   !*********************************************************************************************************************
   ! NS: Adapted from computation of circulation and relative vorticity at each vertex in atm_compute_solve_diagnostics()
   !     This takes scvt face values and computes finite volume curl at scvt vertices (triangle cell centers)
   ! MC: Modified NS's original curl subroutine to include calculation over all vertical levels
   !*********************************************************************************************************************

   subroutine calc_vertical_curl(uEdge, nEdges, nVertices, dcEdge, areaTriangle, verticesOnEdge, curlVert)

      implicit none

      integer, intent(in) :: nEdges, nVertices
      integer, dimension(:,:), intent(in) :: verticesOnEdge
      real (kind=RKIND), dimension(:), intent(in) :: dcEdge, areaTriangle
      real (kind=RKIND), dimension(:,:), intent(in) :: uEdge
      real (kind=RKIND), dimension(:,:), intent(out) :: curlVert

      integer :: jEdge, iVert
     
      curlVert(:,:) = 0.0_RKIND

      do jEdge=1,nEdges
          curlVert(:,verticesOnEdge(1,jEdge)) = curlVert(:,verticesOnEdge(1,jEdge)) - dcEdge(jEdge) * uEdge(:,jEdge)
          curlVert(:,verticesOnEdge(2,jEdge)) = curlVert(:,verticesOnEdge(2,jEdge)) + dcEdge(jEdge) * uEdge(:,jEdge)
      end do

      do iVert=1,nVertices
          curlVert(:,iVert) = curlVert(:,iVert) / areaTriangle(iVert)
      end do

   end subroutine calc_vertical_curl

   !*********************************************************************************************************************
   ! MC: Subroutine combining NS's original functions for calculating vertical derivatives, which finds values at adjacent 
   !     theta/mass levels and then calculates one-sided difference between center level and the levels above and below. 
   !     For all levels except k=1 and k=nVertLevels, these differences are then averaged to give the center difference at 
   !     the center level. Else, the one-sided differences are used. 
   ! 03/20/24: Fix this routine by using difference in zgrid_cell rather than dzu (difference in zeta, not z)
   !*********************************************************************************************************************

   subroutine calc_vertDeriv(var, nCells, nVertLevels, zCell, dvar_dz)

      IMPLICIT NONE

      integer, intent(in) :: nCells, nVertLevels
      real(kind=RKIND), dimension(:,:), intent(in) :: var, zCell
      real(kind=RKIND), dimension(:,:), intent(out) :: dvar_dz
      integer :: iCell, k
      real(kind=RKIND) :: dvar_dz_top, dvar_dz_bot

      dvar_dz(:,:) = 0.0_RKIND

      do iCell=1,nCells
         ! one-sided differences at top and bottom levels
         dvar_dz(1,iCell) = (var(2,iCell) - var(1,iCell)) / (zCell(2,iCell) - zCell(1,iCell))
         dvar_dz(nVertLevels,iCell) = (var(nVertLevels,iCell) - var(nVertLevels-1,iCell)) / &
                                       (zCell(nVertLevels,iCell) - zCell(nVertLevels-1,iCell))
         do k=2,nVertLevels-1
            dvar_dz_top = (var(k+1,iCell) - var(k,iCell)) / (zCell(k+1,iCell) - zCell(k,iCell))
            dvar_dz_bot = (var(k,iCell) - var(k-1,iCell)) / (zCell(k,iCell) - zCell(k-1,iCell))
            ! Currently top and bottom gradients are weighted equally by taking simple average 
            dvar_dz(k,iCell) = 0.5 * (dvar_dz_top + dvar_dz_bot)
         end do
      end do

   end subroutine calc_vertDeriv

   !*********************************************************************************************************************
   ! MC: Alternative method of calculating the vertical derivatives on mass levels, which calculates the vertical gradient
   !     at the lowest mass level by first extrapolating fields to the underlying w level and interpolating to the overlying 
   !     w level, and then calculating the center difference. A one-sided difference is still used at the top model level. 
   !     For all other mass levels, a weighted average of the one-sided differences is used to 
   !     calculate the center differences. 
   !*********************************************************************************************************************

   subroutine calc_vertDeriv_alt(var, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dvar_dz)

      IMPLICIT NONE

      integer, intent(in) :: nCells, nVertLevels
      real(kind=RKIND), intent(in) :: cf1, cf2, cf3

      real(kind=RKIND), dimension(:), intent(in) ::  dzu, rdzw
      real(kind=RKIND), dimension(:,:), intent(in) :: var, zgrid, zCell
      real(kind=RKIND), dimension(:,:), intent(out) :: dvar_dz

      integer :: iCell, k
      real(kind=RKIND) :: dvar_dz_top, dvar_dz_bot, var_w2, var_w1
      real(kind=RKIND), dimension(:), allocatable :: dzw
      
      allocate(dzw(nVertLevels+1))

      dzw(:) = 1./rdzw(:)
      dvar_dz(:,:) = 0.0_RKIND

      do iCell=1,nCells
        ! for bottom mass level, extrapolate to w level below,
        ! interpolate to w level above, and then take the center diff. 
         var_w1 = cf1 * var(1,iCell) + cf2 * var(2,iCell) + cf3 * var(3,iCell)
         var_w2 = (0.5/dzu(2)) * (dzw(2)*var(1,iCell) + dzw(1)*var(2,iCell)) 
         dvar_dz(1,iCell) = (var_w2 - var_w1)/(zgrid(2,iCell) - zgrid(1,iCell))

        ! one-sided differences at top level
         dvar_dz(nVertLevels,iCell) = (var(nVertLevels,iCell) - var(nVertLevels-1,iCell)) / &
                                       (zCell(nVertLevels,iCell) - zCell(nVertLevels-1,iCell))
         do k=2,nVertLevels-1
            dvar_dz_top = (var(k+1,iCell) - var(k,iCell)) / (zCell(k+1,iCell) - zCell(k,iCell))
            dvar_dz_bot = (var(k,iCell) - var(k-1,iCell)) / (zCell(k,iCell) - zCell(k-1,iCell))

            ! Alter weighting to weight bottom derivative more than top since levels closer together
            dvar_dz(k,iCell) = dzu(k+1)/(dzu(k) + dzu(k+1))*dvar_dz_bot + dzu(k)/(dzu(k) + dzu(k+1))*dvar_dz_top
         end do
      end do
      
      deallocate(dzw)

   end subroutine calc_vertDeriv_alt


   !*********************************************************************************************************************   
   ! MC: Subroutine to interpolate variable from w levels (vertical cell faces) to theta levels (cell centers)
   !*********************************************************************************************************************

   subroutine interp_wLev_thetaLev(w, nCells, nVertLevels, wCell)

      IMPLICIT NONE

      integer, intent(in) :: nCells, nVertLevels
      real(kind=RKIND), dimension(:,:), intent(in) :: w
      real(kind=RKIND), dimension(:,:), intent(out) :: wCell
      integer :: iCell, k

      do iCell=1,nCells
        do k=1,nVertLevels
           wCell(k,iCell) = 0.5*(w(k+1,iCell) + w(k,iCell))
        end do
      end do

   end subroutine interp_wLev_thetaLev

   !*********************************************************************************************************************   
   ! MC: Subroutine to store variables from the beginning of the time step to use in next timestep tendency calculations 
   !*********************************************************************************************************************

   subroutine store_previous_vars(mesh, time_lev, state, diag)

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
 
      IMPLICIT NONE

      type (mpas_pool_type), intent(in) :: mesh
      integer, intent(in) :: time_lev                                ! which time level to use from state: 1 = beginning of time step; 2 = end of time step
      type (mpas_pool_type), intent(in) :: state
      type (mpas_pool_type), intent(inout) :: diag

      integer, pointer :: nCells, nVertLevels, nEdges, nVertices

      real(kind=RKIND), dimension(:,:), pointer :: uReconstructZonal, uReconstructMeridional, wCell, theta, rho, &
                                                   pv_vertex, ertel_pv
      real(kind=RKIND), dimension(:,:), pointer :: uReconstructZonal_prev, uReconstructMeridional_prev, wCell_prev, &
                                                   theta_prev, qv_prev, rho_prev, pv_vertex_prev, ertel_pv_prev

      integer, pointer :: index_qv
      real(kind=RKIND), dimension(:,:,:), pointer :: scalars

      integer, dimension(:), pointer :: iLev_DT, iLev_DT_prev

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)

      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'wCell', wCell)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(state, 'scalars', scalars, 1)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)

      call mpas_pool_get_array(diag, 'uReconstructZonal_prev', uReconstructZonal_prev)
      call mpas_pool_get_array(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev)
      call mpas_pool_get_array(diag, 'wCell_prev', wCell_prev)
      call mpas_pool_get_array(diag, 'theta_prev', theta_prev)
      call mpas_pool_get_array(diag, 'qv_prev', qv_prev)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(diag, 'rho_prev', rho_prev)      
      call mpas_pool_get_array(diag, 'pv_vertex_prev', pv_vertex_prev)
      call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev)
      call mpas_pool_get_array(diag, 'iLev_DT_prev', iLev_DT_prev)

      uReconstructZonal_prev(:,:) = uReconstructZonal(:,:)
      uReconstructMeridional_prev(:,:) = uReconstructMeridional(:,:)
      wCell_prev(:,:) = wCell(:,:)
      theta_prev(:,:) = theta(:,:)
      qv_prev(:,:) = scalars(index_qv,:,:)
      rho_prev(:,:) = rho(:,:)
      pv_vertex_prev(:,:) = pv_vertex(:,:)
      ertel_pv_prev(:,:) = ertel_pv(:,:)
      iLev_DT_prev(:) = iLev_DT(:)

      !call mpas_log_write("maxval ertel_pv_prev after resetting is $r", realArgs=(/real(MAXVAL(ertel_pv_prev), kind=RKIND)/))

   end subroutine store_previous_vars

   !*********************************************************************************************************************   
   ! MW: Calculate density tendency term as part of the EPV dynamics tendency
   !*********************************************************************************************************************

   subroutine calc_density_term(rho, rho_prev, ertel_pv_prev, nCells, nVertLevels, dt, drho_dt)

      IMPLICIT NONE
      
      integer, intent(in) :: nCells, nVertLevels
      real(kind=RKIND), intent(in) :: dt
      real(kind=RKIND), dimension(:,:), intent(in) :: rho, rho_prev, ertel_pv_prev
      real(kind=RKIND), dimension(:,:), intent(out) :: drho_dt

      integer :: k, iCell

      do iCell=1,nCells
         do k=1,nVertLevels
           drho_dt(k,iCell) = (rho(k,iCell) - rho_prev(k,iCell)) / (rho(k,iCell)*dt)
           drho_dt(k,iCell) = ertel_pv_prev(k,iCell) * drho_dt(k,iCell)
         end do
      end do 

   end subroutine calc_density_term
   
   !*********************************************************************************************************************
   ! MC: Modified subroutine to calculate Ertel's potential vorticity
   !     PV = 1/density * [curl(wind) . grad(theta)] 
   !*********************************************************************************************************************

   subroutine calc_epv(mesh, time_lev, state, diag)

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

      IMPLICIT NONE

      type (mpas_pool_type), intent(in) :: mesh  
      integer, intent(in) :: time_lev                     ! MC: subroutines called after state/diag vars have been updated, so time_lev spec doesn't matter
      type (mpas_pool_type), intent(in) :: state
      type (mpas_pool_type), intent(inout) :: diag         

      ! Input variables -- mesh
      integer, pointer :: nCells, nVertLevels, nEdges, R3
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex
      real(kind=RKIND), dimension(:), pointer :: dzu, areaCell, latCell, lonCell, dcEdge
      real(kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex, edgesOnCell_sign, zgrid, zCell
      real(kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct 

      ! Input variables -- state/diagnostic
      real(kind=RKIND), dimension(:,:), pointer :: w, wCell, rho, theta, pv_vertex, uReconstructZonal, &
                                                   uReconstructMeridional, ertel_pv, dTheta_dz
                                                   
      ! Local variables
      real(kind=RKIND), dimension(:,:), allocatable :: duZonal_dz, duMerid_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dTheta_dxZonal, dTheta_dyMerid
      real(kind=RKIND), dimension(:,:), allocatable :: dW_dxZonal, dW_dyMerid 
      real(kind=RKIND), dimension(:,:), allocatable :: absVort
      real(kind=RKIND), dimension(:,:,:), allocatable :: absVort3D, gradTheta

      ! Uncomment if using calc_vertDeriv_alt 
      !real(kind=RKIND), pointer :: cf1, cf2, cf3
      !real(kind=RKIND), dimension(:), pointer :: rdzw

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'R3', R3)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)  
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'coeffs_reconstruct', coeffs_reconstruct)

      call mpas_pool_get_array(diag, 'zgrid_cell', zCell) 
      call mpas_pool_get_array(state, 'w', w, time_lev)  
      call mpas_pool_get_array(diag, 'wCell', wCell)    
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)  
      call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'dtheta_dz', dTheta_dz)

      ! Uncomment if using calc_vertDeriv_alt       
      !call mpas_pool_get_array(mesh, 'cf1', cf1)
      !call mpas_pool_get_array(mesh, 'cf2', cf2)
      !call mpas_pool_get_array(mesh, 'cf3', cf3)
      !call mpas_pool_get_array(mesh, 'rdzw', rdzw)
     
      ! Allocate memory to intermediate vars 
      allocate(absVort(nVertLevels,nCells+1)) 
      allocate(duZonal_dz(nVertLevels,nCells+1))
      allocate(duMerid_dz(nVertLevels,nCells+1))
      allocate(dTheta_dxZonal(nVertLevels,nCells+1))
      allocate(dTheta_dyMerid(nVertLevels,nCells+1))
      allocate(dW_dxZonal(nVertLevels,nCells+1))
      allocate(dW_dyMerid(nVertLevels,nCells+1))
      allocate(absVort3D(nVertLevels,nCells+1,3))
      allocate(gradTheta(nVertLevels,nCells+1,3))

      ertel_pv(:,:) = 0.0_RKIND
      gradTheta(:,:,:) = 0.0_RKIND
      absVort3D(:,:,:) = 0.0_RKIND

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D potential temperature gradient 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate and reconstruct horizontal potential temperature gradient to get zonal and meridional
      !     gradients at cell centers: dth_dx, dth_dy

      call calc_gradOnEdges_reconCellCenter(theta, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dTheta_dxZonal, dTheta_dyMerid)

      ! (2) Calculate the vertical potential temperature gradient: dth_dz

      call calc_vertDeriv(theta, nCells, nVertLevels, zCell, dTheta_dz)

      ! For alternative method, comment out above and uncomment below (and in all locations where vertical 
      ! derivatve is calculated). Note: the PV budget residual is lower when using the default method. 
      !
      ! call calc_vertDeriv_alt(theta, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dTheta_dz)

      ! (3) Combine theta derivatives into 3D vector

      gradTheta(:,:,1) = dTheta_dxZonal 
      gradTheta(:,:,2) = dTheta_dyMerid 
      gradTheta(:,:,3) = dTheta_dz 

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D absolute vorticity vector  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate the vertical shear of uReconstructZonal and uReconstructMeridional: du_dz and dv_dz

      call calc_vertDeriv(uReconstructZonal, nCells, nVertLevels, zCell, duZonal_dz)
      call calc_vertDeriv(uReconstructMeridional, nCells, nVertLevels, zCell, duMerid_dz)

      ! For alternative method, comment out above and uncomment below (and in all locations where vertical 
      ! derivatve is calculated). Note: the PV budget residual is lower when using the default method. 
      !
      ! call calc_vertDeriv_alt(uReconstructZonal, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duZonal_dz)
      ! call calc_vertDeriv_alt(uReconstructMeridional, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duMerid_dz)

      ! (2) Interpolate w to cell centers, calculate gradient of w on edges, and then reconstruct to get 
      !     zonal and meridional gradients at cell centers: dw_dx, dw_dy 

      call interp_wLev_thetaLev(w, nCells, nVertLevels, wCell)

      call calc_gradOnEdges_reconCellCenter(wCell, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dW_dxZonal, dW_dyMerid)

      ! (3) Reconstruct absolute vertical vorticity at vertices pv_vertex to cell centers

      !     Note: currently, pv_vertex is the absolute vertical vorticity on the cell vertices. If this 
      !     variable changes at some point, then the absolute vertical vorticity on the vertices needs 
      !     to be computed as follows: 
      !     do iVert=1,nVertices
      !        vorticity(:,iVert) = vorticity(:,iVert) + fVertex(iVert)
      !     end do

      call interp_absVertVort(pv_vertex, nCells, nEdgesOnCell, verticesOnCell, &
                                  cellsOnVertex, areaCell, kiteAreasOnVertex, absVort)

      ! (4) Combine three components into vorticity vector

      absVort3D(:,:,1) = dW_dyMerid - duMerid_dz                   ! dw/dy - dv/dz
      absVort3D(:,:,2) = duZonal_dz - dW_dxZonal                   ! du/dz - dw/dx
      absVort3D(:,:,3) = absVort                                   ! dv/dy - du/dx + f

      ! (5) Take dot product between 3D theta vector and absolute vorticity / density to compute 
      !     Ertel's PV
 
      call calc_dotProduct_3D(gradTheta, absVort3D, nCells, nVertLevels, ertel_pv)

      ertel_pv = ertel_pv / rho * 1.0e6                            !SI to PVUs

      !call mpas_log_write("minval ertel_pv: $r", realArgs=(/real(MINVAL(ertel_pv), kind=RKIND)/))
      !call mpas_log_write("maxval ertel_pv: $r", realArgs=(/real(MAXVAL(ertel_pv), kind=RKIND)/))

      deallocate(absVort)
      deallocate(absVort3D)
      deallocate(gradTheta)
      deallocate(duZonal_dz)
      deallocate(duMerid_dz)
      deallocate(dTheta_dxZonal)
      deallocate(dTheta_dyMerid)
      deallocate(dW_dxZonal)
      deallocate(dW_dyMerid)

   end subroutine calc_epv

   !*********************************************************************************************************************
   ! MC: Modified subroutine to call PV and PV advection calculations and interpolation of diagnostic fields to dynamic 
   !     tropopause
   !********************************************************************************************************************* 

   subroutine atm_compute_pv_diagnostics(configs, state, time_lev, diag, mesh) 
   
      use mpas_constants
      use mpas_derived_types, only : field2DReal
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
      use mpas_dmpar, only : mpas_dmpar_exch_halo_field
 
      implicit none

      type (mpas_pool_type), intent(in) :: configs 
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev                                           
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal, stratoPV

      type (field2DReal), pointer :: theta_f, uReconstructZonal_f, uReconstructMeridional_f, w_f, wCell_f, rho_f, &
                                     pv_vertex_f, epv_f

      call mpas_pool_get_field(diag, 'theta', theta_f)         
      call mpas_pool_get_field(diag, 'uReconstructZonal', uReconstructZonal_f)                     ! MC added
      call mpas_pool_get_field(diag, 'uReconstructMeridional', uReconstructMeridional_f)           ! MC added
      call mpas_pool_get_field(state, 'w', w_f, 1)                                                 ! MC added
      call mpas_pool_get_field(diag, 'wCell', wCell_f)                                             ! MC added 
      call mpas_pool_get_field(diag, 'rho', rho_f)                                                 ! MC added (MC_TODO: not sure if rho needed)
      call mpas_pool_get_field(diag, 'pv_vertex', pv_vertex_f)                                     ! MC added
                                         
      ! Need halo cells for everything w/ horizontal derivative 
      call mpas_dmpar_exch_halo_field(theta_f)
      call mpas_dmpar_exch_halo_field(uReconstructZonal_f)                               
      call mpas_dmpar_exch_halo_field(uReconstructMeridional_f)                         
      call mpas_dmpar_exch_halo_field(w_f)
      call mpas_dmpar_exch_halo_field(wCell_f)                                    
      call mpas_dmpar_exch_halo_field(rho_f)                                             ! MC_TODO: not sure if needed
      call mpas_dmpar_exch_halo_field(pv_vertex_f)                  
 
      ! Call calc_epv subroutine 
      call calc_epv(mesh, time_lev, state, diag)                  

      ! Halo cells need to be valid for flood fill routines called below
      call mpas_pool_get_field(diag, 'ertel_pv', epv_f)                       
      call mpas_dmpar_exch_halo_field(epv_f)
    
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      stratoPV = 10.0_RKIND

      !***********************************************************************************************
      ! Uncomment one of the following to either call floodFill_strato or floodFill_tropo to estimate 
      ! the level of the dynamic tropopause. 
      ! MC note: the two should not longer produce equivalent estimates of iLev_DT, and floodFill_tropo
      ! has been motified to mitigate spurious identification of the dynamic tropopause, which then
      ! leads to erroneous values of interpolated values. Thus, floodFill_tropo should be used by default.
      !***********************************************************************************************
     
      !call floodFill_strato(mesh, diag, pvuVal, stratoPV)    
      call floodFill_tropo(mesh, diag, pvuVal)

      ! Call interp_pv_diagnostics subroutine to interpolate fields to dynamic tropopause
      call interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
   
   end subroutine atm_compute_pv_diagnostics

   !*********************************************************************************************************************
   ! MC: Modified subroutine to calculate Ertel's potential vorticity tendency using state and diagnostic variables at 
   !     the BEGINNING of the previous time step (i.e., before they are updated by the tendencies from dynamics and physics) 
   !     and the tendencies from dynamics and physics over the previous time step. This is necessary via product rule: 
   ! 
   !     dPV/dt = 1/density * [curl(wind) . grad(theta_tendency) + curl(wind_tendency) . grad(theta)] 
   !
   !     and differs from the original formulation, which incorrectly used the updated state and diagnostic variables 
   !     at the end of the time step and the tendencies responsible for updating them!
   !*********************************************************************************************************************

   subroutine calc_pvBudget(configs, state, time_lev, diag, mesh, tend, tend_physics, diag_physics)
 
      use mpas_vector_reconstruction
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config
 
      implicit none
      
      type (mpas_pool_type), intent(in) :: configs 
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev                                            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: tend_physics
      type (mpas_pool_type), intent(in) :: diag_physics
      type (mpas_pool_type), intent(inout) :: tend                             

      ! mesh / configuration variables
      real(kind=RKIND), pointer :: config_dt
      logical, pointer :: config_pv_microphys
      
      integer, pointer :: nCells, nVertLevels, nEdges, R3, nVertices
      integer, dimension(:),   pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex, verticesOnEdge 
                                          
      real(kind=RKIND), dimension(:),     pointer :: dzu, areaCell, latCell, lonCell, dcEdge, areaTriangle                                      
      real(kind=RKIND), dimension(:,:),   pointer :: zgrid, zCell, kiteAreasOnVertex, edgesOnCell_sign
      real(kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct
      
      ! need to call stored variables from beginning of previous time step 
      real(kind=RKIND), dimension(:,:), pointer :: rho_prev, pv_vertex_prev, uReconstructZonal_prev, & 
                                                   uReconstructMeridional_prev, wCell_prev, ertel_pv_prev, qv_prev
      ! t+dt variables
      real(kind=RKIND), dimension(:,:), pointer :: rho, theta
      
      ! diabatic PV tendencies
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, &
                                                   depv_dt_mix, depv_dt_diab
      ! friction PV tendencies
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_fric, depv_dt_fric_bl, depv_dt_fric_mix, depv_dt_fric_cu

      ! dynamics PV tendency
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_dyn

      ! specific microphysics PV tendencies
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_evap_cw, depv_dt_mp_evap_rw, depv_dt_mp_depo_ice, &     
                                                   depv_dt_mp_melt_ice, depv_dt_mp_frez_ice, depv_dt_mp_allproc
      
      ! process tendencies -- friction
      real(kind=RKIND), dimension(:,:), pointer :: u_tend_diff, w_tend_diff, tend_wCell_diff, tend_u_pbl, tend_u_cu
      real(kind=RKIND), dimension(:,:), pointer :: uTend_curl_diff, uTend_curl_pbl, uTend_curl_cu                         ! Could these be allocated? MC_TODO; is allocation of curl for mixing the source of the nCellsSolve issue?

      ! process tendencies -- diabatic 
      real(kind=RKIND), dimension(:,:), pointer :: dtheta_dt_mp, dtheta_dt_mix, dtheta_dt_pbl, dtheta_dt_cu, dtheta_dt_sw, &  ! MC: calculated in atm_compute_pvBudget_diagnostics
                                                   dtheta_dt_lw
      real(kind=RKIND), dimension(:,:), pointer :: tend_theta_mp_evap_cw, tend_theta_mp_evap_rw, tend_theta_mp_depo_ice, &  
                                                   tend_theta_mp_melt_ice, tend_theta_mp_frez_ice

      ! process tendencies -- dynamics 
      real(kind=RKIND), dimension(:,:), pointer :: dtheta_dt_dyn, du_dt_dyn, dw_dt_dyn, tenddyn_wCell
      real(kind=RKIND), dimension(:,:), pointer :: uTend_curl_dyn
     
      ! local static variables 
      real(kind=RKIND), dimension(:,:), allocatable :: duZonal_dz, duMerid_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dTheta_dxZonal, dTheta_dyMerid, dTheta_dz           ! note: dtheta_dz is stored. make consistent
      real(kind=RKIND), dimension(:,:), allocatable :: dW_dxZonal, dW_dyMerid
      real(kind=RKIND), dimension(:,:), allocatable :: absVort
      real(kind=RKIND), dimension(:,:,:), allocatable :: absVort3D, gradTheta

      ! diabatic tendency variables
      real(kind=RKIND), dimension(:,:), allocatable :: dLWtend_dxZonal, dLWtend_dyMerid, dLWtend_dz      ! Gradients of theta tendencies from LW radiation
      real(kind=RKIND), dimension(:,:), allocatable :: dSWtend_dxZonal, dSWtend_dyMerid, dSWtend_dz      ! Gradients of theta tendencies from SW radiation
      real(kind=RKIND), dimension(:,:), allocatable :: dBLtend_dxZonal, dBLtend_dyMerid, dBLtend_dz      ! Gradients of theta tendencies from PBL scheme
      real(kind=RKIND), dimension(:,:), allocatable :: dCUtend_dxZonal, dCUtend_dyMerid, dCUtend_dz      ! Gradients of theta tendencies from cumulus scheme
      real(kind=RKIND), dimension(:,:), allocatable :: dMPtend_dxZonal, dMPtend_dyMerid, dMPtend_dz      ! Gradients of theta tendencies from microphysics
      real(kind=RKIND), dimension(:,:), allocatable :: dMXtend_dxZonal, dMXtend_dyMerid, dMXtend_dz      ! Gradients of theta tendencies from explicit horiz mixing 
      real(kind=RKIND), dimension(:,:,:), allocatable :: grad_diabatic_LW, grad_diabatic_SW, grad_diabatic_BL, &
                                                         grad_diabatic_CU, grad_diabatic_MP, grad_diabatic_MX

      real(kind=RKIND), dimension(:,:), allocatable :: tend_theta_mp_sum
      real(kind=RKIND), dimension(:,:), allocatable :: dMPevapcwtend_dxZonal, dMPevapcwtend_dyMerid, &   ! Gradients of theta tendencies from specific microphysical processes
                                                       dMPevapcwtend_dz, &
                                                       dMPevaprwtend_dxZonal, dMPevaprwtend_dyMerid, &
                                                       dMPevaprwtend_dz, &
                                                       dMPdepotend_dxZonal, dMPdepotend_dyMerid, &
                                                       dMPdepotend_dz, &
                                                       dMPmelttend_dxZonal, dMPmelttend_dyMerid, &
                                                       dMPmelttend_dz, &
                                                       dMPfreztend_dxZonal, dMPfreztend_dyMerid, &
                                                       dMPfreztend_dz, &
                                                       dMPsumtend_dxZonal, dMPsumtend_dyMerid, &
                                                       dMPsumtend_dz
      real(kind=RKIND), dimension(:,:,:), allocatable :: grad_diabatic_MP_evap_cw, grad_diabatic_MP_evap_rw, grad_diabatic_MP_depo, &
                                                         grad_diabatic_MP_melt, grad_diabatic_MP_frez, grad_diabatic_MP_sum

      ! friction tendency vars 
      real(kind=RKIND), dimension(:,:), allocatable :: dWtend_dxZonal, dWtend_dyMerid, duZonalTend_dz_mix, duMeridTend_dz_mix, &
                                                       vertVortTend_mix, tenduX_mix, tenduY_mix, tenduZ_mix, tend_uZonal_mix, tend_uMerid_mix
      real(kind=RKIND), dimension(:,:), allocatable :: duZonalTend_dz_pbl, duMeridTend_dz_pbl, &
                                                       vertVortTend_pbl, tenduX_pbl, tenduY_pbl, tenduZ_pbl, tend_uZonal_pbl, tend_uMerid_pbl       
      real(kind=RKIND), dimension(:,:), allocatable :: duZonalTend_dz_cu, duMeridTend_dz_cu, &
                                                       vertVortTend_cu, tenduX_cu, tenduY_cu, tenduZ_cu, tend_uZonal_cu, tend_uMerid_cu
      real(kind=RKIND), dimension(:,:), allocatable :: dWtend_dxZonal_phys, dWtend_dyMerid_phys
      real(kind=RKIND), dimension(:,:,:), allocatable :: vortTend3D_mix, vortTend3D_pbl, vortTend3D_cu

      ! dynamics vars
      real(kind=RKIND), dimension(:,:), allocatable :: depv_dt_graddyn, depv_dt_vortdyn, drho_dt_term, tenduX_dyn, tenduY_dyn, tenduZ_dyn
      real(kind=RKIND), dimension(:,:), allocatable :: tend_uZonal_dyn, tend_uMerid_dyn, dDYNtend_dxZonal, dDYNtend_dyMerid, &
                                                       dDYNtend_dz, vertVortTend_dyn, duZonalTend_dz, duMeridTend_dz
      real(kind=RKIND), dimension(:,:,:), allocatable :: grad_DYN, vortTend3D_DYN      

      ! needed for alternative vertical derivative calculation (not recommended)
      !real(kind=RKIND), pointer :: cf1, cf2, cf3
      !real(kind=RKIND), dimension(:), pointer :: rdzw

      call mpas_log_write("In PV budget calculation subroutine.")
      
      ! mesh / config vars
      call mpas_pool_get_config(configs,'config_dt',config_dt)
      call mpas_pool_get_config(configs, 'config_pv_microphys', config_pv_microphys)
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'R3', R3)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'coeffs_reconstruct', coeffs_reconstruct)
      call mpas_pool_get_array(diag, 'zgrid_cell', zCell)

      ! time-level t variables
      call mpas_pool_get_array(diag, 'rho_prev', rho_prev)
      call mpas_pool_get_array(diag, 'pv_vertex_prev', pv_vertex_prev)
      call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev)
      call mpas_pool_get_array(diag, 'uReconstructZonal_prev', uReconstructZonal_prev)
      call mpas_pool_get_array(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev)
      call mpas_pool_get_array(diag, 'wCell_prev', wCell_prev)
      
      ! time-level t+dt variables
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta', theta)

      ! diabatic PV tendencies
      call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)                            ! diabatic PV tendency from longwave radiation
      call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)                            ! diabatic PV tendency from shortwave radiation
      call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)                            ! diabatic PV tendency from PBL scheme
      call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)                            ! diabatic PV tendency from cumulus scheme 
      call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)                            ! diabatic PV tendency from microphysics scheme
      call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)                          ! diabatic PV tendency from explict mixing 
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)                        ! total diabatic PV tendency 
      
      ! frictional PV tendencies
      call mpas_pool_get_array(diag, 'depv_dt_fric_bl', depv_dt_fric_bl)                  ! frictional PV tendency from PBL + GWD schemes 
      call mpas_pool_get_array(diag, 'depv_dt_fric_cu', depv_dt_fric_cu)                  ! frictional PV tendency from cumulus scheme (only nonzero if scheme modifies winds)
      call mpas_pool_get_array(diag, 'depv_dt_fric_mix', depv_dt_fric_mix)                ! frictional PV tendency from explicit mixing 
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)                        ! total frictional PV tendency 

      ! dynamics PV tendency
      call mpas_pool_get_array(diag,'depv_dt_dyn',depv_dt_dyn)                            ! total PV tendency from dynamics (includes transport, decoupling, rho tendency) 
      
      ! specific microphysics PV tendencies (Thompson only)
      ! MC_TODO: these should only be activated when config_pv_microphys is true. 
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_cw', depv_dt_mp_evap_cw)            ! diabatic PV tendency from net condensation/evaporation of cloud water
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_rw', depv_dt_mp_evap_rw)            ! diabatic PV tendency from evaporation of rain water
      call mpas_pool_get_array(diag, 'depv_dt_mp_depo_ice', depv_dt_mp_depo_ice)          ! diabatic PV tendency from net deposition/sublimation of all ice hydrometeors
      call mpas_pool_get_array(diag, 'depv_dt_mp_melt_ice', depv_dt_mp_melt_ice)          ! diabatic PV tendency from melting of all ice hydrometeors
      call mpas_pool_get_array(diag, 'depv_dt_mp_frez_ice', depv_dt_mp_frez_ice)          ! diabatic PV tendency from freezing of all ice hydrometeors
      call mpas_pool_get_array(diag, 'depv_dt_mp_allproc', depv_dt_mp_allproc)            ! total diabatic PV tendency from all specific processes 

      
      ! process tendencies -- friction:
      call mpas_pool_get_array(diag, 'u_tend_diff', u_tend_diff)                          ! Normal wind tendency from explicit mixing on cell edges
      call mpas_pool_get_array(diag, 'w_tend_diff', w_tend_diff)                          ! Vertical wind tendency from explicit mixing 
      call mpas_pool_get_array(diag, 'tend_wCell_diff', tend_wCell_diff)                  ! w_tend_diff interpolated to mass levels 
      call mpas_pool_get_array(diag, 'tend_u_pbl', tend_u_pbl)                            ! Normal wind tendency from PBL + GWD on cell edges
      call mpas_pool_get_array(diag, 'tend_u_cu', tend_u_cu)                              ! Normal wind tendency from cumulus scheme (only nonzero if scheme modifies winds)
      call mpas_pool_get_array(diag, 'uTend_curl_diff', uTend_curl_diff)                  ! Vertical curl of u_tend_diff at cell vertices 
      call mpas_pool_get_array(diag, 'uTend_curl_pbl', uTend_curl_pbl)                    ! Vertical curl of tend_u_pbl at cell vertices 
      call mpas_pool_get_array(diag, 'uTend_curl_cu', uTend_curl_cu)                      ! Vertical curl of tend_u_cu at cell vertices 

      ! process tendencies -- diabatic:
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)                      ! Derived potential temperature tendency from explicit horizontal mixing
      call mpas_pool_get_array(diag, 'dtheta_dt_cu', dtheta_dt_cu)                        ! Derived potential temperature tendency from cumulus
      call mpas_pool_get_array(diag, 'dtheta_dt_pbl', dtheta_dt_pbl)                      ! Derived potential temperature tendency from PBL 
      call mpas_pool_get_array(diag, 'dtheta_dt_sw', dtheta_dt_sw)                        ! Derived potential temperature tendency from SW radiation
      call mpas_pool_get_array(diag, 'dtheta_dt_lw', dtheta_dt_lw)                        ! Derived potential temperature tendency from LW radiation
      call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)                        ! Derived potential temperature tendency from microphysics 
      ! MC: theta tendencies directly output from schemes are no longer used 
      !call mpas_pool_get_array(tend_physics, 'rthblten', dtheta_dt_pbl)
      !call mpas_pool_get_array(tend_physics, 'rthcuten', dtheta_dt_cu)
      !call mpas_pool_get_array(tend_physics, 'rthratenlw', dtheta_dt_lw)
      !call mpas_pool_get_array(tend_physics, 'rthratensw', dtheta_dt_sw)

      ! specific microphysics process potential temperature tendencies
      ! MC note: unlike dtheta_dt_mp, which is derived from theta_m tendency and is more precise,
      !          these tendencies are directly obtained from the theta rates in the Thompson
      !          scheme. thus, the sum of these will not exactly equal dtheta_dt_mp.
      ! MC_TODO: these should only be activated when config_pv_microphys is true. 
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_evap_cw', tend_theta_mp_evap_cw)     
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_evap_rw', tend_theta_mp_evap_rw)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_depo_ice', tend_theta_mp_depo_ice)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_melt_ice', tend_theta_mp_melt_ice)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_frez_ice', tend_theta_mp_frez_ice)

      ! process tendencies -- dynamics:
      call mpas_pool_get_array(diag, 'dtheta_dt_dyn', dtheta_dt_dyn)                       ! Derived potential temperature tendency from dynamics 
      call mpas_pool_get_array(diag, 'du_dt_dyn', du_dt_dyn)                               ! Normal wind tendency from dynamics on cell edges
      call mpas_pool_get_array(diag, 'dw_dt_dyn', dw_dt_dyn)                               ! Vertical wind tendency from dynamics
      call mpas_pool_get_array(diag, 'tenddyn_wCell', tenddyn_wCell)                       ! dw_dt_dyn interpolated to mass levels 
      call mpas_pool_get_array(diag, 'uTend_curl_dyn', uTend_curl_dyn)                     ! Vertical curl of du_dt_dyn at cell vertices 
       
      ! needed for alternative vertical derivative calculation (not recommended)
      !call mpas_pool_get_array(mesh, 'cf1', cf1)
      !call mpas_pool_get_array(mesh, 'cf2', cf2)
      !call mpas_pool_get_array(mesh, 'cf3', cf3)
      !call mpas_pool_get_array(mesh, 'rdzw', rdzw)

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Allocate Variables
      
      ! local static vars
      allocate(duZonal_dz(nVertLevels,nCells+1))
      allocate(duMerid_dz(nVertLevels,nCells+1))
      allocate(dTheta_dxZonal(nVertLevels,nCells+1))
      allocate(dTheta_dyMerid(nVertLevels,nCells+1))
      allocate(dTheta_dz(nVertLevels,nCells+1))
      allocate(dW_dxZonal(nVertLevels,nCells+1))
      allocate(dW_dyMerid(nVertLevels,nCells+1))
      allocate(absVort(nVertLevels,nCells+1))
      ! 3D static vectors
      allocate(absVort3D(nVertLevels,nCells+1,R3))
      allocate(gradTheta(nVertLevels,nCells+1,R3))

      ! allocate diabatic tendency variables
      allocate(dLWtend_dxZonal(nVertLevels,nCells+1))
      allocate(dLWtend_dyMerid(nVertLevels,nCells+1))
      allocate(dLWtend_dz(nVertLevels,nCells+1))
      allocate(dSWtend_dxZonal(nVertLevels,nCells+1))
      allocate(dSWtend_dyMerid(nVertLevels,nCells+1))
      allocate(dSWtend_dz(nVertLevels,nCells+1))
      allocate(dBLtend_dxZonal(nVertLevels,nCells+1))
      allocate(dBLtend_dyMerid(nVertLevels,nCells+1))
      allocate(dBLtend_dz(nVertLevels,nCells+1))
      allocate(dCUtend_dxZonal(nVertLevels,nCells+1))
      allocate(dCUtend_dyMerid(nVertLevels,nCells+1))
      allocate(dCUtend_dz(nVertLevels,nCells+1))
      allocate(dMPtend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPtend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPtend_dz(nVertLevels,nCells+1))
      allocate(dMXtend_dxZonal(nVertLevels,nCells+1))
      allocate(dMXtend_dyMerid(nVertLevels,nCells+1))
      allocate(dMXtend_dz(nVertLevels,nCells+1))
      ! 3D tendency vectors
      allocate(grad_diabatic_LW(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_SW(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_BL(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_CU(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MX(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MP(nVertLevels,nCells+1,R3))

      ! allocate diabatic tendency variables from specific microphys processes
      allocate(dMPevapcwtend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPevapcwtend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPevapcwtend_dz(nVertLevels,nCells+1))
      allocate(dMPevaprwtend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPevaprwtend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPevaprwtend_dz(nVertLevels,nCells+1))
      allocate(dMPdepotend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPdepotend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPdepotend_dz(nVertLevels,nCells+1))
      allocate(dMPmelttend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPmelttend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPmelttend_dz(nVertLevels,nCells+1))
      allocate(dMPfreztend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPfreztend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPfreztend_dz(nVertLevels,nCells+1))
      allocate(dMPsumtend_dxZonal(nVertLevels,nCells+1))
      allocate(dMPsumtend_dyMerid(nVertLevels,nCells+1))
      allocate(dMPsumtend_dz(nVertLevels,nCells+1))
      allocate(tend_theta_mp_sum(nVertLevels,nCells+1))
      ! 3D tendency vectors
      allocate(grad_diabatic_MP_evap_cw(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MP_evap_rw(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MP_depo(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MP_melt(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MP_frez(nVertLevels,nCells+1,R3))
      allocate(grad_diabatic_MP_sum(nVertLevels,nCells+1,R3))

      ! allocate friction tendency variables
      ! mixing
      allocate(dWtend_dxZonal(nVertLevels,nCells+1))     ! also used for dynamics
      allocate(dWtend_dyMerid(nVertLevels,nCells+1))     ! also used for dynamics
      allocate(duZonalTend_dz_mix(nVertLevels,nCells+1))
      allocate(duMeridTend_dz_mix(nVertLevels,nCells+1))
      allocate(tend_uZonal_mix(nVertLevels,nCells+1))    ! reconstructing in tend subroutines
      allocate(tend_uMerid_mix(nVertLevels,nCells+1))    ! reconstructing in tend subroutines
      allocate(vertVortTend_mix(nVertLevels,nCells+1))
      allocate(tenduX_mix(nVertLevels,nCells+1))
      allocate(tenduY_mix(nVertLevels,nCells+1))
      allocate(tenduZ_mix(nVertLevels,nCells+1))
      ! PBL
      allocate(duZonalTend_dz_pbl(nVertLevels,nCells+1))
      allocate(duMeridTend_dz_pbl(nVertLevels,nCells+1))
      allocate(tend_uZonal_pbl(nVertLevels,nCells+1))    ! reconstructing in tend subroutines
      allocate(tend_uMerid_pbl(nVertLevels,nCells+1))    ! reconstructing in tend subroutines
      allocate(vertVortTend_pbl(nVertLevels,nCells+1))
      allocate(tenduX_pbl(nVertLevels,nCells+1))
      allocate(tenduY_pbl(nVertLevels,nCells+1))
      allocate(tenduZ_pbl(nVertLevels,nCells+1))
      ! cumulus
      allocate(duZonalTend_dz_cu(nVertLevels,nCells+1))
      allocate(duMeridTend_dz_cu(nVertLevels,nCells+1))
      allocate(tend_uZonal_cu(nVertLevels,nCells+1))    ! reconstructing in tend subroutines
      allocate(tend_uMerid_cu(nVertLevels,nCells+1))    ! reconstructing in tend subroutines
      allocate(vertVortTend_cu(nVertLevels,nCells+1))
      allocate(tenduX_cu(nVertLevels,nCells+1))
      allocate(tenduY_cu(nVertLevels,nCells+1))
      allocate(tenduZ_cu(nVertLevels,nCells+1))
      ! zeroed variables for w tendency from phys
      allocate(dWtend_dxZonal_phys(nVertLevels,nCells+1))
      allocate(dWtend_dyMerid_phys(nVertLevels,nCells+1))
      ! 3D tendency vectors
      allocate(vortTend3D_mix(nVertLevels,nCells+1,R3))
      allocate(vortTend3D_pbl(nVertLevels,nCells+1,R3))
      allocate(vortTend3D_cu(nVertLevels,nCells+1,R3))

      ! allocate dynamics tendency variables
      allocate(duZonalTend_dz(nVertLevels,nCells+1))
      allocate(duMeridTend_dz(nVertLevels,nCells+1))
      allocate(tend_uZonal_dyn(nVertLevels,nCells+1))  
      allocate(tend_uMerid_dyn(nVertLevels,nCells+1))
      allocate(vertVortTend_dyn(nVertLevels,nCells+1))
      allocate(tenduX_dyn(nVertLevels,nCells+1))
      allocate(tenduY_dyn(nVertLevels,nCells+1))
      allocate(tenduZ_dyn(nVertLevels,nCells+1))
      allocate(vortTend3D_DYN(nVertLevels,nCells+1,R3))
      allocate(dDYNtend_dxZonal(nVertLevels,nCells+1))   ! used for dyn theta tendency
      allocate(dDYNtend_dyMerid(nVertLevels,nCells+1))   ! used for dyn theta tendency
      allocate(dDYNtend_dz(nVertLevels,nCells+1))        ! used for dyn theta tendency
      allocate(grad_DYN(nVertLevels,nCells+1,R3))
      allocate(depv_dt_graddyn(nVertLevels,nCells+1))
      allocate(depv_dt_vortdyn(nVertLevels,nCells+1))
      allocate(drho_dt_term(nVertLevels,nCells+1))
      
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Initialize vars
      depv_dt_lw(:,:) = 0.0_RKIND
      depv_dt_sw(:,:) = 0.0_RKIND
      depv_dt_bl(:,:) = 0.0_RKIND
      depv_dt_cu(:,:) = 0.0_RKIND
      depv_dt_mp(:,:) = 0.0_RKIND 
      depv_dt_mix(:,:) = 0.0_RKIND
      depv_dt_diab(:,:) = 0.0_RKIND
      depv_dt_fric(:,:) = 0.0_RKIND
      depv_dt_fric_bl(:,:) = 0.0_RKIND
      depv_dt_fric_mix(:,:) = 0.0_RKIND 
      depv_dt_fric_cu(:,:) = 0.0_RKIND
      depv_dt_mp_evap_cw(:,:) =  0.0_RKIND
      depv_dt_mp_evap_rw(:,:) =  0.0_RKIND
      depv_dt_mp_depo_ice(:,:) =  0.0_RKIND
      depv_dt_mp_melt_ice(:,:) =  0.0_RKIND
      depv_dt_mp_frez_ice(:,:) =  0.0_RKIND
      tend_theta_mp_sum(:,:) =  0.0_RKIND
      depv_dt_mp_allproc(:,:) =  0.0_RKIND
      depv_dt_dyn(:,:) = 0.0_RKIND
      depv_dt_graddyn(:,:) = 0.0_RKIND
      depv_dt_vortdyn(:,:) = 0.0_RKIND
      drho_dt_term(:,:) = 0.0_RKIND

      ! Gradient of w tendency from phys (remains 0)
      dWtend_dxZonal_phys(:,:) = 0.0_RKIND
      dWtend_dyMerid_phys(:,:) = 0.0_RKIND


      !***********************************************************************************************
      ! Calculate terms needed for PV tendency equation
      !***********************************************************************************************

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D potential temperature gradient using theta at end of time step  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate and reconstruct horizontal potential temperature gradient to get zonal and meridional
      !     gradients at cell centers: dth_dx, dth_dy
      call calc_gradOnEdges_reconCellCenter(theta, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dTheta_dxZonal, dTheta_dyMerid)

      ! (2) Calculate the vertical potential temperature gradient: dth_dz
      call calc_vertDeriv(theta, nCells, nVertLevels, zCell, dTheta_dz)

      ! For alternative method, comment out above and uncomment below. 
      ! call calc_vertDeriv_alt(theta, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dTheta_dz)

      ! (3) Combine theta derivatives into 3D vector
      gradTheta(:,:,1) = dTheta_dxZonal
      gradTheta(:,:,2) = dTheta_dyMerid
      gradTheta(:,:,3) = dTheta_dz

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D absolute vorticity vector using winds at beginning of time step  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate the vertical shear of uReconstructZonal and uReconstructMeridional: du_dz and dv_dz
      call calc_vertDeriv(uReconstructZonal_prev, nCells, nVertLevels, zCell, duZonal_dz)
      call calc_vertDeriv(uReconstructMeridional_prev, nCells, nVertLevels, zCell, duMerid_dz)

      ! For alternative method, comment out above and uncomment below. 
      ! call calc_vertDeriv_alt(uReconstructZonal_prev, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duZonal_dz)
      ! call calc_vertDeriv_alt(uReconstructMeridional_prev, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duMerid_dz)

      ! (2) Interpolate w to cell centers, calculate gradient of w on edges, and then reconstruct to get 
      !     zonal and meridional gradients at cell centers: dw_dx, dw_dy 
      call calc_gradOnEdges_reconCellCenter(wCell_prev, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dW_dxZonal, dW_dyMerid)

      ! (3) Reconstruct absolute vertical vorticity at vertices pv_vertex to cell centers
      !     Note: currently, pv_vertex is the absolute vertical vorticity on the cell vertices. If this 
      !     variable changes at some point, then the absolute vertical vorticity on the vertices needs 
      !     to be computed as follows: 
      !     do iVert=1,nVertices
      !        vorticity(:,iVert) = vorticity(:,iVert) + fVertex(iVert)
      !     end do
      call interp_absVertVort(pv_vertex_prev, nCells, nEdgesOnCell, verticesOnCell, &
                                cellsOnVertex, areaCell, kiteAreasOnVertex, absVort)

      ! (4) Combine three components into vorticity vector
      absVort3D(:,:,1) = dW_dyMerid - duMerid_dz        ! dw/dy - dv/dz
      absVort3D(:,:,2) = duZonal_dz - dW_dxZonal        ! du/dz - dw/dx
      absVort3D(:,:,3) = absVort                        ! dv/dy - du/dx + f

      !***********************************************************************************************
      ! Calculate diabatic PV tendency terms:
      !***********************************************************************************************

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Longwave radiation tendency: depv_dt_lw
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_lw)) then
      
         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
         call calc_gradOnEdges_reconCellCenter(dtheta_dt_lw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                  edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                  latCell, lonCell, dLWtend_dxZonal, &
                                  dLWtend_dyMerid)

         ! (2) Calculate vertical gradient of theta tendency
         call calc_vertDeriv(dtheta_dt_lw, nCells, nVertLevels, zCell, dLWtend_dz)

         ! For alternative method, comment out above and uncomment below. 
         ! call calc_vertDeriv_alt(dtheta_dt_lw, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dLWtend_dz)

         ! (3) Combine into 3D theta tendency gradient vector 
         grad_diabatic_LW(:,:,1) = dLWtend_dxZonal
         grad_diabatic_LW(:,:,2) = dLWtend_dyMerid
         grad_diabatic_LW(:,:,3) = dLWtend_dz

         ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
         call calc_dotProduct_3D(grad_diabatic_LW, absVort3D, nCells, nVertLevels, depv_dt_lw)
         
         depv_dt_lw = depv_dt_lw / rho * 1.0e6                             
      else
         depv_dt_lw = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Shortwave radiation tendency: depv_dt_sw
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_sw)) then         
      
         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_sw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dSWtend_dxZonal, &
                                   dSWtend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_sw, nCells, nVertLevels, zCell, dSWtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(dtheta_dt_sw, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dSWtend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_SW(:,:,1) = dSWtend_dxZonal
          grad_diabatic_SW(:,:,2) = dSWtend_dyMerid
          grad_diabatic_SW(:,:,3) = dSWtend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_SW, absVort3D, nCells, nVertLevels, depv_dt_sw)
 
          depv_dt_sw = depv_dt_sw / rho * 1.0e6
      else
          depv_dt_sw = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PBL diabatic tendency: depv_dt_bl
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_pbl)) then 
      
         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_pbl, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dBLtend_dxZonal, &
                                   dBLtend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_pbl, nCells, nVertLevels, zCell, dBLtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(dtheta_dt_pbl, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dBLtend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_BL(:,:,1) = dBLtend_dxZonal
          grad_diabatic_BL(:,:,2) = dBLtend_dyMerid
          grad_diabatic_BL(:,:,3) = dBLtend_dz
          
         ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_BL, absVort3D, nCells, nVertLevels, depv_dt_bl)
 
          depv_dt_bl = depv_dt_bl / rho * 1.0e6
      else
          depv_dt_bl = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Cumulus diabatic tendency: depv_dt_cu
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_cu)) then
      
         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_cu, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dCUtend_dxZonal, &
                                   dCUtend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_cu, nCells, nVertLevels, zCell, dCUtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(dtheta_dt_cu, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dCUtend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_CU(:,:,1) = dCUtend_dxZonal
          grad_diabatic_CU(:,:,2) = dCUtend_dyMerid
          grad_diabatic_CU(:,:,3) = dCUtend_dz

         ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_CU, absVort3D, nCells, nVertLevels, depv_dt_cu)
 
          depv_dt_cu = depv_dt_cu / rho * 1.0e6
      else
          depv_dt_cu = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Microphysics diabatic tendency: depv_dt_mp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_mp)) then
      
         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_mp, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPtend_dxZonal, &
                                   dMPtend_dyMerid)
 
          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_mp, nCells, nVertLevels, zCell, dMPtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(dtheta_dt_mp, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPtend_dz)
 
          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP(:,:,1) = dMPtend_dxZonal
          grad_diabatic_MP(:,:,2) = dMPtend_dyMerid
          grad_diabatic_MP(:,:,3) = dMPtend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP, absVort3D, nCells, nVertLevels, depv_dt_mp)
 
          depv_dt_mp = depv_dt_mp / rho * 1.0e6
      else
          depv_dt_mp = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Diabatic tendency from explicit mixing: depv_dt_mix
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_mix)) then
      
         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_mix, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMXtend_dxZonal, &
                                   dMXtend_dyMerid)
 
          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_mix, nCells, nVertLevels, zCell, dMXtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(dtheta_dt_mix, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMXtend_dz)
 
          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MX(:,:,1) = dMXtend_dxZonal
          grad_diabatic_MX(:,:,2) = dMXtend_dyMerid
          grad_diabatic_MX(:,:,3) = dMXtend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MX, absVort3D, nCells, nVertLevels, depv_dt_mix)
 
          depv_dt_mix = depv_dt_mix / rho * 1.0e6
      else
          depv_dt_mix = 0.0_RKIND
      end if
 
      ! Sum of all diabatic contributions to PV through potential temperature tendencies   
      depv_dt_diab = depv_dt_mix + depv_dt_lw + depv_dt_sw + depv_dt_bl + depv_dt_cu + depv_dt_mp 

 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Individual diabatic tendencies from specific microphysical processes (not included in budget)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !call mpas_log_write('config_pv_microphys is: $l', logicArgs=(/config_pv_microphys/))
      
      ! Net cloud water condensation and evaporation
      if ((config_pv_microphys) .and. (associated(tend_theta_mp_evap_cw))) then

          ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_evap_cw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPevapcwtend_dxZonal, &
                                   dMPevapcwtend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_evap_cw, nCells, nVertLevels, zCell, dMPevapcwtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_theta_mp_evap_cw, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPevapcwtend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_evap_cw(:,:,1) = dMPevapcwtend_dxZonal
          grad_diabatic_MP_evap_cw(:,:,2) = dMPevapcwtend_dyMerid
          grad_diabatic_MP_evap_cw(:,:,3) = dMPevapcwtend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP_evap_cw, absVort3D, nCells, nVertLevels, depv_dt_mp_evap_cw)

          depv_dt_mp_evap_cw = depv_dt_mp_evap_cw / rho * 1.0e6
      else
          depv_dt_mp_evap_cw = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Rain water evaporation
      if ((config_pv_microphys) .and. (associated(tend_theta_mp_evap_rw))) then

          ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_evap_rw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPevaprwtend_dxZonal, &
                                   dMPevaprwtend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_evap_rw, nCells, nVertLevels, zCell, dMPevaprwtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_theta_mp_evap_rw, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPevaprwtend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_evap_rw(:,:,1) = dMPevaprwtend_dxZonal
          grad_diabatic_MP_evap_rw(:,:,2) = dMPevaprwtend_dyMerid
          grad_diabatic_MP_evap_rw(:,:,3) = dMPevaprwtend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP_evap_rw, absVort3D, nCells, nVertLevels, depv_dt_mp_evap_rw)

          depv_dt_mp_evap_rw = depv_dt_mp_evap_rw / rho * 1.0e6
      else
          depv_dt_mp_evap_rw = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Net sublimation/deposition
      if ((config_pv_microphys) .and. (associated(tend_theta_mp_depo_ice))) then

          ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_depo_ice, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPdepotend_dxZonal, &
                                   dMPdepotend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_depo_ice, nCells, nVertLevels, zCell, dMPdepotend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_theta_mp_depo_ice, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPdepotend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_depo(:,:,1) = dMPdepotend_dxZonal
          grad_diabatic_MP_depo(:,:,2) = dMPdepotend_dyMerid
          grad_diabatic_MP_depo(:,:,3) = dMPdepotend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP_depo, absVort3D, nCells, nVertLevels, depv_dt_mp_depo_ice)

          depv_dt_mp_depo_ice = depv_dt_mp_depo_ice / rho * 1.0e6

      else
          depv_dt_mp_depo_ice = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Melting
      if ((config_pv_microphys) .and. (associated(tend_theta_mp_melt_ice))) then

          ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_melt_ice, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPmelttend_dxZonal, &
                                   dMPmelttend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_melt_ice, nCells, nVertLevels, zCell, dMPmelttend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_theta_mp_melt_ice, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPmelttend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_melt(:,:,1) = dMPmelttend_dxZonal
          grad_diabatic_MP_melt(:,:,2) = dMPmelttend_dyMerid
          grad_diabatic_MP_melt(:,:,3) = dMPmelttend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP_melt, absVort3D, nCells, nVertLevels, depv_dt_mp_melt_ice)

          depv_dt_mp_melt_ice = depv_dt_mp_melt_ice / rho * 1.0e6
      else
          depv_dt_mp_melt_ice = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Freezing
      if ((config_pv_microphys) .and. (associated(tend_theta_mp_frez_ice))) then

          ! (1) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_frez_ice, nCells, nEdges, nVertLevels, nEdgesOnCell, &   
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPfreztend_dxZonal, &
                                   dMPfreztend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_frez_ice, nCells, nVertLevels, zCell, dMPfreztend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_theta_mp_frez_ice, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPfreztend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_frez(:,:,1) = dMPfreztend_dxZonal
          grad_diabatic_MP_frez(:,:,2) = dMPfreztend_dyMerid
          grad_diabatic_MP_frez(:,:,3) = dMPfreztend_dz

          ! (4) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP_frez, absVort3D, nCells, nVertLevels, depv_dt_mp_frez_ice)

          depv_dt_mp_frez_ice = depv_dt_mp_frez_ice / rho * 1.0e6
      else
          depv_dt_mp_frez_ice = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Sum all processes together to find combined PV tendency from microphysics. Compare to depv_dt_mp
      if ((config_pv_microphys) .and. (associated(tend_theta_mp_frez_ice))) then

          ! (1) Sum individual process tendencies 
          tend_theta_mp_sum = tend_theta_mp_frez_ice + tend_theta_mp_melt_ice + tend_theta_mp_depo_ice + &
                             tend_theta_mp_evap_rw + tend_theta_mp_evap_cw

          ! (2) Calculate and reconstruct horizontal gradients of theta tendency
          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_sum, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPsumtend_dxZonal, &
                                   dMPsumtend_dyMerid)

          ! (3) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_sum, nCells, nVertLevels, zCell, dMPsumtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_theta_mp_sum, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dMPsumtend_dz)

          ! (4) Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_sum(:,:,1) = dMPsumtend_dxZonal
          grad_diabatic_MP_sum(:,:,2) = dMPsumtend_dyMerid
          grad_diabatic_MP_sum(:,:,3) = dMPsumtend_dz

          ! (5) Take dot product between 3D theta tendency vector and absolute vorticity / density 
          call calc_dotProduct_3D(grad_diabatic_MP_sum, absVort3D, nCells, nVertLevels, depv_dt_mp_allproc)

          depv_dt_mp_allproc = depv_dt_mp_allproc / rho * 1.0e6
      else
          depv_dt_mp_allproc = 0.0_RKIND
      end if


      !***********************************************************************************************
      ! Calculate frictional tendency terms:
      !
      ! Friction terms are essentially the vorticity tendency due to friction. Need to use the u, v, w
      ! tendencies
      !***********************************************************************************************
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Frictional tendency from explicit mixing: depv_dt_fric_mix
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if ((associated(u_tend_diff)) .and. (associated(w_tend_diff))) then

         ! (1) Reconstruct u tendency from diffusion to cell center 
         call mpas_reconstruct(mesh, u_tend_diff, tenduX_mix, tenduY_mix, tenduZ_mix, &
                                   tend_uZonal_mix, tend_uMerid_mix)

         ! (2) Calculate vertical derivative of tend_uZonal_mix, tend_uMerid_mix
         call calc_vertDeriv(tend_uZonal_mix, nCells, nVertLevels, zCell, duZonalTend_dz_mix)
         call calc_vertDeriv(tend_uMerid_mix, nCells, nVertLevels, zCell, duMeridTend_dz_mix)

         ! For alternative method, comment out above and uncomment below. 
         ! call calc_vertDeriv_alt(tend_uZonal_mix, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duZonalTend_dz_mix)
         ! call calc_vertDeriv_alt(tend_uMerid_mix, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duMeridTend_dz_mix)

         ! (3) Interpolate w tendency from diffusion to cell center, calclulate gradient of tend_wCell_diff on edges,
         !     and then reconstruct to get zonal and meridional gradients at cell center:
         !     dWtend_dxZonal, dWtend_dyMerid
         call interp_wLev_thetaLev(w_tend_diff, nCells, nVertLevels, tend_wCell_diff)

         call calc_gradOnEdges_reconCellCenter(tend_wCell_diff, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dWtend_dxZonal, dWtend_dyMerid)

         ! (4) Calculate vertical vorticity tendency on vertices and reconstruct to cell center
         !     Note: the analogy of this procedure in the ertel_pv calculation is missing because vertical vorticity
         !           is already computed and output in MPAS 
         call calc_vertical_curl(u_tend_diff, nEdges, nVertices, dcEdge, areaTriangle, verticesOnEdge, uTend_curl_diff)

         call interp_absVertVort(uTend_curl_diff, nCells, nEdgesOnCell, verticesOnCell, &
                                  cellsOnVertex, areaCell, kiteAreasOnVertex, vertVortTend_mix)

         ! (5) Combine three components into vorticity tendency vector
         vortTend3D_mix(:,:,1) = dWtend_dyMerid - duMeridTend_dz_mix        ! dFz/dy - dFy/dz
         vortTend3D_mix(:,:,2) = duZonalTend_dz_mix - dWtend_dxZonal        ! dFx/dz - dFz/dx
         vortTend3D_mix(:,:,3) = vertVortTend_mix                           ! dFy/dy - dFx/dx 

         ! (6) Take dot product between 3D theta gradient and absolute vorticity tendency vector / density 
         call calc_dotProduct_3D(gradTheta, vortTend3D_mix, nCells, nVertLevels, depv_dt_fric_mix)

         depv_dt_fric_mix = depv_dt_fric_mix / rho * 1.0e6
      else
          depv_dt_fric_mix = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Frictional tendency from PBL and GWD schemes: depv_dt_fric_bl
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(tend_u_pbl)) then 

         ! (1) Reconstruct u tendency from PBL to cell center 
         call mpas_reconstruct(mesh, tend_u_pbl, tenduX_pbl, tenduY_pbl, tenduZ_pbl, tend_uZonal_pbl, tend_uMerid_pbl)

         ! (2) Calculate vertical derivative of tend_uZonal_pbl, tend_uMerid_pbl
         call calc_vertDeriv(tend_uZonal_pbl, nCells, nVertLevels, zCell, duZonalTend_dz_pbl)
         call calc_vertDeriv(tend_uMerid_pbl, nCells, nVertLevels, zCell, duMeridTend_dz_pbl)

         ! For alternative method, comment out above and uncomment below. 
         ! call calc_vertDeriv_alt(tend_uZonal_pbl, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duZonalTend_dz_pbl)
         ! call calc_vertDeriv_alt(tend_uMerid_pbl, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duMeridTend_dz_pbl)

         ! (3) Since w tendency from PBL scheme is zero, just set dWtend_dxZonal_phys, dWtend_dyMerid_phys = 0.0 (done above) 
         !     If model is ever updated to include a cumulus scheme that alters w, will need to modify this.

         ! (4) Calculate vertical vorticity tendency on vertices and reconstruct to cell center
         !     Note: the analogy of this procedure in the ertel_pv calculation is missing because vertical vorticity
         !           is already computed and output in MPAS 
         call calc_vertical_curl(tend_u_pbl, nEdges, nVertices, dcEdge, areaTriangle, verticesOnEdge, uTend_curl_pbl)

         call interp_absVertVort(uTend_curl_pbl, nCells, nEdgesOnCell, verticesOnCell, &
                               cellsOnVertex, areaCell, kiteAreasOnVertex, vertVortTend_pbl)

         ! (5) Combine three components into vorticity tendency vector
         vortTend3D_pbl(:,:,1) = dWtend_dyMerid_phys - duMeridTend_dz_pbl        ! dFz/dy - dFy/dz
         vortTend3D_pbl(:,:,2) = duZonalTend_dz_pbl - dWtend_dxZonal_phys        ! dFx/dz - dFz/dx
         vortTend3D_pbl(:,:,3) = vertVortTend_pbl                                ! dFy/dy - dFx/dx 

         ! (6) Take dot product between 3D theta gradient and absolute vorticity tendency vector / density 
         call calc_dotProduct_3D(gradTheta, vortTend3D_pbl, nCells, nVertLevels, depv_dt_fric_bl)

         depv_dt_fric_bl = depv_dt_fric_bl / rho * 1.0e6
      else
          depv_dt_fric_bl = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Cumulus scheme 
      if (associated(tend_u_cu)) then
 
         ! (1) Reconstruct u tendency from cumulus to cell center and calculate vertical derivative of 
         !     tend_uZonal, tend_uMerid
         call mpas_reconstruct(mesh, tend_u_cu, tenduX_cu, tenduY_cu, tenduZ_cu, tend_uZonal_cu, tend_uMerid_cu)

         ! (2) Calculate vertical derivative of tend_uZonal_pbl, tend_uMerid_pbl
         call calc_vertDeriv(tend_uZonal_cu, nCells, nVertLevels, zCell, duZonalTend_dz_cu)
         call calc_vertDeriv(tend_uMerid_cu, nCells, nVertLevels, zCell, duMeridTend_dz_cu)

         ! For alternative method, comment out above and uncomment below. 
         ! call calc_vertDeriv_alt(tend_uZonal_cu, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duZonalTend_dz_cu)
         ! call calc_vertDeriv_alt(tend_uMerid_cu, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duMeridTend_dz_cu)
 
         ! (3) Since w tendency from PBL scheme is zero, just set dWtend_dxZonal_phys, dWtend_dyMerid_phys = 0.0 (done above) 
         !     If model is ever updated to include a cumulus scheme that alters w, will need to modify this.

         ! (4) Calculate vertical vorticity tendency on vertices and reconstruct to cell center
         !     Note: the analogy of this procedure in the ertel_pv calculation is missing because vertical vorticity
         !           is already computed and output in MPAS 
         call calc_vertical_curl(tend_u_cu, nEdges, nVertices, dcEdge, areaTriangle, verticesOnEdge, uTend_curl_cu)

         call interp_absVertVort(uTend_curl_cu, nCells, nEdgesOnCell, verticesOnCell, &
                               cellsOnVertex, areaCell, kiteAreasOnVertex, vertVortTend_cu)

         ! (5) Combine three components into vorticity tendency vector
         vortTend3D_cu(:,:,1) = dWtend_dyMerid_phys - duMeridTend_dz_cu        ! dFz/dy - dFy/dz
         vortTend3D_cu(:,:,2) = duZonalTend_dz_cu - dWtend_dxZonal_phys        ! dFx/dz - dFz/dx
         vortTend3D_cu(:,:,3) = vertVortTend_cu                                ! dFy/dy - dFx/dx 

         ! (6) Take dot product between 3D theta gradient and absolute vorticity tendency vector / density 
         call calc_dotProduct_3D(gradTheta, vortTend3D_cu, nCells, nVertLevels, depv_dt_fric_cu)

         depv_dt_fric_cu = depv_dt_fric_cu / rho * 1.0e6
      else
          depv_dt_fric_cu = 0.0_RKIND
      end if


      ! Sum of all frictional contributions to PV through momentum tendencies   
       depv_dt_fric = depv_dt_fric_mix + depv_dt_fric_bl + depv_dt_fric_cu 


      !***********************************************************************************************
      ! Calculate dynamics tendency term : depv_dt_dyn 
      !***********************************************************************************************
      ! --------------------------------------------
      ! The theta gradient dynamics tendency piece
      ! --------------------------------------------
      if (associated(dtheta_dt_dyn)) then 

         ! (1) Calculate and reconstruct horizontal gradients of theta tendency
         call calc_gradOnEdges_reconCellCenter(dtheta_dt_dyn, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dDYNtend_dxZonal, &
                                   dDYNtend_dyMerid)

          ! (2) Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_dyn, nCells, nVertLevels, zCell, dDYNtend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(dtheta_dt_dyn, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, dDYNtend_dz)

          ! (3) Combine into 3D theta tendency gradient vector 
          grad_DYN(:,:,1) = dDYNtend_dxZonal
          grad_DYN(:,:,2) = dDYNtend_dyMerid
          grad_DYN(:,:,3) = dDYNtend_dz
          
          ! (4) Take dot product between 3D theta tendency and absolute vorticity vector / density 
          call calc_dotProduct_3D(grad_DYN, absVort3D, nCells, nVertLevels, depv_dt_graddyn)

          depv_dt_graddyn = depv_dt_graddyn / rho * 1.0e6 
      else
          depv_dt_graddyn = 0.0_RKIND
      end if
      
      ! --------------------------------------------
      ! The vorticity dynamics tendency piece
      ! --------------------------------------------
      if (associated(du_dt_dyn)) then

          ! (1) Reconstruct du_dt_dyn to cell center and 
          call mpas_reconstruct(mesh, du_dt_dyn, tenduX_dyn, tenduY_dyn, tenduZ_dyn, tend_uZonal_dyn, tend_uMerid_dyn)

          ! (2) Calculate vertical derivative of tend_uZonal_dyn and tend_uMerid_dyn
          call calc_vertDeriv(tend_uZonal_dyn, nCells, nVertLevels, zCell, duZonalTend_dz)
          call calc_vertDeriv(tend_uMerid_dyn, nCells, nVertLevels, zCell, duMeridTend_dz)

          ! For alternative method, comment out above and uncomment below. 
          ! call calc_vertDeriv_alt(tend_uZonal_dyn, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duZonalTend_dz)
          ! call calc_vertDeriv_alt(tend_uMerid_dyn, nCells, nVertLevels, dzu, cf1, cf2, cf3, rdzw, zgrid, zCell, duMeridTend_dz)
          
          ! (3) Interpolate w tendency from dynamics to cell center, uncouple from density, calclulate gradient 
          !     of tend_wCell on edges, and then reconstruct to get zonal and meridional gradients at cell center:
          !     dWtend_dxZonal, dWtend_dyMerid
          !
          !     Note: NS's procedure interpolated density to w levels before uncoupling, but this method 
          !           is consistent with the procedure in the vorticity calculation 
          call interp_wLev_thetaLev(dw_dt_dyn, nCells, nVertLevels, tenddyn_wCell)

          call calc_gradOnEdges_reconCellCenter(tenddyn_wCell, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dWtend_dxZonal, dWtend_dyMerid)

          ! (4) Calculate vertical vorticity tendency on vertices and reconstruct to cell center
          !     Note: the analogy of this procedure in the ertel_pv calculation is missing because vertical vorticity
          !           is already computed and output in MPAS 
          call calc_vertical_curl(du_dt_dyn, nEdges, nVertices, dcEdge, areaTriangle, verticesOnEdge, uTend_curl_dyn)

          call interp_absVertVort(uTend_curl_dyn, nCells, nEdgesOnCell, verticesOnCell, &
                               cellsOnVertex, areaCell, kiteAreasOnVertex, vertVortTend_dyn)

          ! (5) Combine three components into vorticity tendency vector
          vortTend3D_DYN(:,:,1)= dWtend_dyMerid - duMeridTend_dz       ! dFz/dy - dFy/dz
          vortTend3D_DYN(:,:,2)= duZonalTend_dz - dWtend_dxZonal       ! dFx/dz - dFz/dx 
          vortTend3D_DYN(:,:,3)= vertVortTend_dyn                      ! dFy/dy - dFx/dx 
                    
          ! (6) Take dot product between 3D theta gradient and absolute vorticity tendency vector / density 
          call calc_dotProduct_3D(gradTheta, vortTend3D_DYN, nCells, nVertLevels, depv_dt_vortdyn)
 
          depv_dt_vortdyn = depv_dt_vortdyn / rho * 1.0e6
      else
          depv_dt_vortdyn = 0.0_RKIND
      end if

      ! --------------------------------------------
      ! The density tendency piece 
      ! --------------------------------------------

      call calc_density_term(rho, rho_prev, ertel_pv_prev, nCells, nVertLevels, config_dt, drho_dt_term) 

      ! --------------------------------------------
      ! Combine all pieces into full dynamics tendency 
      ! --------------------------------------------
      
      depv_dt_dyn = depv_dt_graddyn + depv_dt_vortdyn - drho_dt_term

      
      ! --------------------------------------------------------------------------------------------------
      ! deallocate local static variables
      deallocate(duZonal_dz)
      deallocate(duMerid_dz)
      deallocate(dTheta_dxZonal)
      deallocate(dTheta_dyMerid)
      deallocate(dTheta_dz)
      deallocate(dW_dxZonal)
      deallocate(dW_dyMerid)
      deallocate(absVort)
      deallocate(absVort3D)
      deallocate(gradTheta)

      ! deallocate diabatic tendency variables
      deallocate(dLWtend_dxZonal)
      deallocate(dLWtend_dyMerid)
      deallocate(dLWtend_dz)
      deallocate(dSWtend_dxZonal)
      deallocate(dSWtend_dyMerid)
      deallocate(dSWtend_dz)
      deallocate(dBLtend_dxZonal)
      deallocate(dBLtend_dyMerid)
      deallocate(dBLtend_dz)
      deallocate(dCUtend_dxZonal)
      deallocate(dCUtend_dyMerid)
      deallocate(dCUtend_dz)
      deallocate(dMPtend_dxZonal)
      deallocate(dMPtend_dyMerid)
      deallocate(dMPtend_dz)
      deallocate(dMXtend_dxZonal)
      deallocate(dMXtend_dyMerid)
      deallocate(dMXtend_dz)
      deallocate(grad_diabatic_LW)
      deallocate(grad_diabatic_SW)
      deallocate(grad_diabatic_BL)
      deallocate(grad_diabatic_CU)
      deallocate(grad_diabatic_MX)
      deallocate(grad_diabatic_MP)
     
     ! deallocate diabatic tendency variables from specific microphys processes
      deallocate(dMPevapcwtend_dxZonal)
      deallocate(dMPevapcwtend_dyMerid)
      deallocate(dMPevapcwtend_dz)
      deallocate(dMPevaprwtend_dxZonal)
      deallocate(dMPevaprwtend_dyMerid)
      deallocate(dMPevaprwtend_dz)
      deallocate(dMPdepotend_dxZonal)
      deallocate(dMPdepotend_dyMerid)
      deallocate(dMPdepotend_dz)
      deallocate(dMPmelttend_dxZonal)
      deallocate(dMPmelttend_dyMerid)
      deallocate(dMPmelttend_dz)
      deallocate(dMPfreztend_dxZonal)
      deallocate(dMPfreztend_dyMerid)
      deallocate(dMPfreztend_dz)
      deallocate(dMPsumtend_dxZonal)
      deallocate(dMPsumtend_dyMerid)
      deallocate(dMPsumtend_dz)
      deallocate(tend_theta_mp_sum)
      deallocate(grad_diabatic_MP_evap_cw)
      deallocate(grad_diabatic_MP_evap_rw)
      deallocate(grad_diabatic_MP_depo)
      deallocate(grad_diabatic_MP_melt)
      deallocate(grad_diabatic_MP_frez)
      deallocate(grad_diabatic_MP_sum)

      ! deallocate friction tendency variables
      ! mixing
      deallocate(dWtend_dxZonal)
      deallocate(dWtend_dyMerid)
      deallocate(duZonalTend_dz_mix)
      deallocate(duMeridTend_dz_mix)
      deallocate(tend_uZonal_mix)
      deallocate(tend_uMerid_mix)
      deallocate(vertVortTend_mix)
      deallocate(tenduX_mix)
      deallocate(tenduY_mix)
      deallocate(tenduZ_mix)
      ! pbl
      deallocate(duZonalTend_dz_pbl)
      deallocate(duMeridTend_dz_pbl)
      deallocate(vertVortTend_pbl)
      deallocate(tenduX_pbl)
      deallocate(tenduY_pbl)
      deallocate(tenduZ_pbl)
      deallocate(tend_uZonal_pbl)  
      deallocate(tend_uMerid_pbl)
      !cumulus
      deallocate(duZonalTend_dz_cu)
      deallocate(duMeridTend_dz_cu)
      deallocate(vertVortTend_cu)
      deallocate(tenduX_cu)
      deallocate(tenduY_cu)
      deallocate(tenduZ_cu)
      deallocate(tend_uZonal_cu)  
      deallocate(tend_uMerid_cu)
      ! zeroed variables for w tendency from phys
      deallocate(dWtend_dxZonal_phys)
      deallocate(dWtend_dyMerid_phys)
      ! 3D tendency vectors
      deallocate(vortTend3D_mix)
      deallocate(vortTend3D_pbl)
      deallocate(vortTend3D_cu)

     ! deallocate dynamics tendency variables
      deallocate(duZonalTend_dz)
      deallocate(duMeridTend_dz)
      deallocate(tend_uZonal_dyn)  
      deallocate(tend_uMerid_dyn)
      deallocate(vertVortTend_dyn)
      deallocate(tenduX_dyn)
      deallocate(tenduY_dyn)
      deallocate(tenduZ_dyn)
      deallocate(vortTend3D_DYN)
      deallocate(dDYNtend_dxZonal)
      deallocate(dDYNtend_dyMerid)
      deallocate(dDYNtend_dz)
      deallocate(grad_DYN)
      deallocate(depv_dt_graddyn)
      deallocate(depv_dt_vortdyn)
      deallocate(drho_dt_term)
      
   end subroutine calc_pvBudget


   !*********************************************************************************************************************
   ! MC: Modified subroutine to call variables at the correct time levels for PV tendency calculations and interpolate
   !     fields to the DT identified at the beginning of the time step. Additionally, this subroutine now contains
   !     calculations of the derived theta tendencies from physical processes and mixing from the corresponding tendencies 
   !     for theta_m and qv (if applicable), consistent with the discretized equations for these variables. Doing this
   !     calculation rather than using the theta tendencies from the physics schemes directly enables us to close the 
   !     theta and PV budgets. 
   !  MC_TODO: CLEAN UP!
   !********************************************************************************************************************* 
   
   subroutine atm_compute_pvBudget_diagnostics(configs, state, time_lev, diag, mesh, tend, tend_physics, diag_physics)
      
      use mpas_constants
      use mpas_derived_types, only : field2DReal
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field, mpas_pool_get_config
      use mpas_dmpar, only : mpas_dmpar_exch_halo_field
      use isobaric_diagnostics, only : isobaric_diagnostics_pv
 
      implicit none

      type (mpas_pool_type), intent(in) :: configs
      type (mpas_pool_type), intent(in) :: state, mesh, tend_physics, diag_physics
      type (mpas_pool_type), intent(inout) :: diag, tend
      
      logical, pointer :: config_pv_isobaric, config_pv_microphys
      integer, intent(in) :: time_lev                                  
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal 

      ! For calculating theta tendencies from theta_m tendencies
      ! Note: this is more accurate than using theta tendencies directly from physics schemes, which 
      ! will lead to residuals in the theta and PV budgets owing to coupling between theta and qv via theta_m
      real (kind=RKIND), dimension(:,:), pointer :: theta, qv_prev
      real (kind=RKIND), dimension(:,:), pointer :: dthetam_dt_dyn, dtheta_dt_dyn, dqv_dt_dyn 
      real (kind=RKIND), dimension(:,:), pointer :: dthetam_dt_mix, dtheta_dt_mix
      real (kind=RKIND), dimension(:,:), pointer :: thmblten, qvblten, dtheta_dt_pbl
      real (kind=RKIND), dimension(:,:), pointer :: thmcuten, qvcuten, dtheta_dt_cu
      real (kind=RKIND), dimension(:,:), pointer :: thmswten, dtheta_dt_sw, thmlwten, dtheta_dt_lw
      real (kind=RKIND), dimension(:,:), pointer :: thmmpten, qvmpten, dtheta_dt_mp

      ! For halo exchange 
      ! static vars  
      type (field2DReal), pointer :: theta_prev_f, uReconstructZonal_prev_f, uReconstructMeridional_prev_f, &
                                     wCell_prev_f, pv_vertex_prev_f, qv_prev_f, theta_f, rho_f, rho_prev_f, &
                                     ertel_pv_prev_f

      ! theta tendencies -- diabatic and dynamics 
      type (field2DReal), pointer :: dtheta_dt_dyn_f, dtheta_dt_mix_f, dtheta_dt_pbl_f, dtheta_dt_cu_f, &
                                     dtheta_dt_sw_f, dtheta_dt_lw_f, dtheta_dt_mp_f

      ! theta tendencies -- specific microphysical processes 
      type (field2DReal), pointer :: tend_theta_mp_evap_cw_f, tend_theta_mp_evap_rw_f, tend_theta_mp_depo_ice_f, & 
                                     tend_theta_mp_frez_ice_f, tend_theta_mp_melt_ice_f

      ! momentum tendencies -- diffusion, physics, dynamics 
      type (field2DReal), pointer :: u_tend_diff_f, w_tend_diff_f, uTend_curl_diff_f, tend_wCell_diff_f 
      type (field2DReal), pointer :: tend_u_pbl_f, uTend_curl_pbl_f, tend_u_cu_f, uTend_curl_cu_f
      type (field2DReal), pointer :: du_dt_dyn_f, dw_dt_dyn_f, tenddyn_wCell_f, uTend_curl_dyn_f
   

      ! -----------------------------------------------------------
      ! Calculate theta tendencies from theta_m tendencies
      ! -----------------------------------------------------------
      
      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      ! Dynamics
      call mpas_pool_get_array(diag, 'qv_prev', qv_prev) 
      call mpas_pool_get_array(diag, 'dqv_dt_dyn', dqv_dt_dyn)
      call mpas_pool_get_array(diag, 'dthetam_dt_dyn', dthetam_dt_dyn)
      call mpas_pool_get_array(diag, 'dtheta_dt_dyn', dtheta_dt_dyn)
      
      ! Horizontal mixing
      call mpas_pool_get_array(diag, 'dthetam_dt_mix', dthetam_dt_mix)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)

      ! PBL heating
      call mpas_pool_get_array(diag, 'thmblten', thmblten)
      call mpas_pool_get_array(diag, 'qvblten', qvblten)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'dtheta_dt_pbl', dtheta_dt_pbl)

      ! cumulus heating
      call mpas_pool_get_array(diag, 'thmcuten', thmcuten)
      call mpas_pool_get_array(diag, 'qvcuten', qvcuten)
      call mpas_pool_get_array(diag, 'dtheta_dt_cu', dtheta_dt_cu)

      ! radiation
      call mpas_pool_get_array(diag, 'thmswten', thmswten)
      call mpas_pool_get_array(diag, 'thmlwten', thmlwten)
      call mpas_pool_get_array(diag, 'dtheta_dt_sw', dtheta_dt_sw)
      call mpas_pool_get_array(diag, 'dtheta_dt_lw', dtheta_dt_lw)

      ! microphysics
      call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)
      call mpas_pool_get_array(diag, 'qvmpten', qvmpten)
      call mpas_pool_get_array(diag, 'thmmpten', thmmpten)

      dtheta_dt_dyn(:,:) = 0.0_RKIND
      dtheta_dt_mix(:,:) = 0.0_RKIND
      dtheta_dt_pbl(:,:) = 0.0_RKIND
      dtheta_dt_cu(:,:) = 0.0_RKIND
      dtheta_dt_lw(:,:) = 0.0_RKIND
      dtheta_dt_sw(:,:) = 0.0_RKIND
      dtheta_dt_mp(:,:) = 0.0_RKIND

      ! MC: Modified the diabatic tendency terms used in the PV diagnostics tendency calculations, which are derived here from the theta_m
      !     diabatic process tendencies. These modifications are necessary to conserve theta and PV
      do iCell = 1,nCells
         do k = 1,nVertLevels
           ! processes that also modify qv
           dtheta_dt_dyn(k,iCell) = (dthetam_dt_dyn(k,iCell) - rvord*theta(k,iCell)*dqv_dt_dyn(k,iCell)) / (1._RKIND + rvord*qv_prev(k,iCell))
           dtheta_dt_pbl(k,iCell) = (thmblten(k,iCell) - rvord*theta(k,iCell)*qvblten(k,iCell)) / (1._RKIND + rvord*qv_prev(k,iCell))
           dtheta_dt_cu(k,iCell) = (thmcuten(k,iCell) - rvord*theta(k,iCell)*qvcuten(k,iCell)) / (1._RKIND + rvord*qv_prev(k,iCell))
           dtheta_dt_mp(k,iCell) = (thmmpten(k,iCell) - rvord*theta(k,iCell)*qvmpten(k,iCell)) / (1._RKIND + rvord*qv_prev(k,iCell)) 

           ! processes that do not modify qv
           dtheta_dt_mix(k,iCell) = dthetam_dt_mix(k,iCell) / (1._RKIND + rvord * qv_prev(k,iCell)) 
           dtheta_dt_lw(k,iCell) = thmlwten(k,iCell) / (1._RKIND + rvord * qv_prev(k,iCell)) 
           dtheta_dt_sw(k,iCell) = thmswten(k,iCell) / (1._RKIND + rvord * qv_prev(k,iCell)) 
        end do
      end do

      ! -----------------------------------------------------------
      ! Do halo communication
      ! -----------------------------------------------------------

      call mpas_pool_get_config(configs, 'config_pv_microphys', config_pv_microphys)

      ! State/diag vars
      call mpas_pool_get_field(diag, 'theta_prev', theta_prev_f)
      call mpas_pool_get_field(diag, 'uReconstructZonal_prev', uReconstructZonal_prev_f)           
      call mpas_pool_get_field(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev_f) 
      call mpas_pool_get_field(diag, 'wCell_prev', wCell_prev_f)                                   
      call mpas_pool_get_field(diag, 'pv_vertex_prev', pv_vertex_prev_f)    
      call mpas_pool_get_field(diag, 'qv_prev', qv_prev_f)
      call mpas_pool_get_field(diag, 'theta', theta_f)
      call mpas_pool_get_field(diag, 'rho', rho_f)
      call mpas_pool_get_field(diag, 'rho_prev', rho_prev_f)
      call mpas_pool_get_field(diag, 'ertel_pv_prev', ertel_pv_prev_f)

      call mpas_dmpar_exch_halo_field(theta_prev_f)
      call mpas_dmpar_exch_halo_field(uReconstructZonal_prev_f)                              
      call mpas_dmpar_exch_halo_field(uReconstructMeridional_prev_f)     
      call mpas_dmpar_exch_halo_field(wCell_prev_f)                                           
      call mpas_dmpar_exch_halo_field(pv_vertex_prev_f)               
      call mpas_dmpar_exch_halo_field(qv_prev_f)
      call mpas_dmpar_exch_halo_field(theta_f)
      call mpas_dmpar_exch_halo_field(rho_f)
      call mpas_dmpar_exch_halo_field(rho_prev_f)
      call mpas_dmpar_exch_halo_field(ertel_pv_prev_f)

      ! Diabatic tendencies 
      call mpas_pool_get_field(diag, 'dtheta_dt_dyn', dtheta_dt_dyn_f) 
      call mpas_pool_get_field(diag, 'dtheta_dt_mix', dtheta_dt_mix_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_pbl', dtheta_dt_pbl_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_cu', dtheta_dt_cu_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_sw', dtheta_dt_sw_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_lw', dtheta_dt_lw_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mp', dtheta_dt_mp_f)

      call mpas_dmpar_exch_halo_field(dtheta_dt_dyn_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mix_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_pbl_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_cu_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_sw_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_lw_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mp_f)

      ! Microphysics process tendencies
      if (config_pv_microphys) then 
          call mpas_pool_get_field(diag_physics, 'tend_theta_mp_evap_cw', tend_theta_mp_evap_cw_f)
          call mpas_pool_get_field(diag_physics, 'tend_theta_mp_evap_rw', tend_theta_mp_evap_rw_f)
          call mpas_pool_get_field(diag_physics, 'tend_theta_mp_depo_ice', tend_theta_mp_depo_ice_f)
          call mpas_pool_get_field(diag_physics, 'tend_theta_mp_melt_ice', tend_theta_mp_melt_ice_f)
          call mpas_pool_get_field(diag_physics, 'tend_theta_mp_frez_ice', tend_theta_mp_frez_ice_f)

          call mpas_dmpar_exch_halo_field(tend_theta_mp_evap_cw_f)
          call mpas_dmpar_exch_halo_field(tend_theta_mp_evap_rw_f)
          call mpas_dmpar_exch_halo_field(tend_theta_mp_depo_ice_f)
          call mpas_dmpar_exch_halo_field(tend_theta_mp_melt_ice_f)
          call mpas_dmpar_exch_halo_field(tend_theta_mp_frez_ice_f)
      end if 

      ! momentum tendencies -- friction and dynamics
      call mpas_pool_get_field(diag, 'u_tend_diff', u_tend_diff_f)        ! MW: using ru explicit horiz mixing term from tendencies
      call mpas_pool_get_field(diag, 'w_tend_diff', w_tend_diff_f)        ! MW: using rw explicit horiz mixing term from tendencies
      call mpas_pool_get_field(diag, 'uTend_curl_diff', uTend_curl_diff_f)  
      call mpas_pool_get_field(diag, 'tend_wCell_diff', tend_wCell_diff_f)
      call mpas_pool_get_field(diag, 'tend_u_pbl', tend_u_pbl_f)
      call mpas_pool_get_field(diag, 'tend_u_cu', tend_u_cu_f)
      call mpas_pool_get_field(diag, 'uTend_curl_pbl', uTend_curl_pbl_f)  
      call mpas_pool_get_field(diag, 'uTend_curl_cu', uTend_curl_cu_f)   
      call mpas_pool_get_field(diag, 'du_dt_dyn', du_dt_dyn_f) 
      call mpas_pool_get_field(diag, 'dw_dt_dyn', dw_dt_dyn_f) 
      call mpas_pool_get_field(diag, 'tenddyn_wCell', tenddyn_wCell_f) 
      call mpas_pool_get_field(diag, 'uTend_curl_dyn', uTend_curl_dyn_f)            

      call mpas_dmpar_exch_halo_field(u_tend_diff_f)
      call mpas_dmpar_exch_halo_field(w_tend_diff_f)
      call mpas_dmpar_exch_halo_field(uTend_curl_diff_f)
      call mpas_dmpar_exch_halo_field(tend_wCell_diff_f)
      call mpas_dmpar_exch_halo_field(tend_u_pbl_f)
      call mpas_dmpar_exch_halo_field(tend_u_cu_f)
      call mpas_dmpar_exch_halo_field(uTend_curl_pbl_f)                  
      call mpas_dmpar_exch_halo_field(uTend_curl_cu_f)                        
      call mpas_dmpar_exch_halo_field(du_dt_dyn_f)
      call mpas_dmpar_exch_halo_field(dw_dt_dyn_f)
      call mpas_dmpar_exch_halo_field(tenddyn_wCell_f)
      call mpas_dmpar_exch_halo_field(uTend_curl_dyn_f)                   
                      
      

      ! -----------------------------------------------------------
      ! Call subroutines:
      ! -----------------------------------------------------------

      ! Calculating PV budget:
      call mpas_log_write("Calling calc_pvBudget:") 
      call calc_pvBudget(configs, state, time_lev, diag, mesh, tend, tend_physics, diag_physics)

      ! Interpolate fields to DT:
      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND

      call mpas_log_write("Calling interp_pvBudget_diagnostics")
      call interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)

      ! Interpolate fields to isobaric levels if desired: 
      call mpas_pool_get_config(configs, 'config_pv_isobaric', config_pv_isobaric)

      if (config_pv_isobaric) then
          call mpas_log_write("Calling isobaric_diagnostics_pv:")
          call isobaric_diagnostics_pv()
      end if 

      ! Accumulate tendencies over multiple time steps: 
      call mpas_log_write("Calling acc_pvBudget:")
      call acc_pvBudget(mesh, diag, tend_physics)
   
   end subroutine atm_compute_pvBudget_diagnostics



   !*********************************************************************************************************************
   ! MW: Adding a subroutine to accumulate the PV budget tendency terms at each time step. Accumulated values are output 
   !     at the user-specific diagnostic output interval, which determines the time-averaging window of the tendencies.
   !********************************************************************************************************************* 

   ! MC_TODO: add friction terms
   ! MC_TODO: do something about restarts
   
   subroutine acc_pvBudget(mesh, diag, tend_physics)
 
      use mpas_vector_reconstruction
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      use mpas_pool_routines, only: mpas_pool_get_config

      implicit none

      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: tend_physics 
      
      logical, pointer :: config_pv_tend, config_pv_isobaric, config_pv_microphys

      integer, pointer :: nCells, nVertLevels
      integer :: iCell, k
      
      ! Instantaneous and accumulated PV tendencies 
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_lw, acc_depv_dt_sw, acc_depv_dt_bl, acc_depv_dt_cu, acc_depv_dt_mp, acc_depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_fric_mix, depv_dt_fric_bl, depv_dt_fric_cu
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_fric_mix, acc_depv_dt_fric_bl, acc_depv_dt_fric_cu
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_diab, depv_dt_fric, depv_dt_dyn
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_diab, acc_depv_dt_fric, acc_depv_dt_dyn
      real(kind=RKIND), dimension(:),   pointer :: depv_dt_diab_pv, depv_dt_fric_pv, depv_dt_dyn_pv
      real(kind=RKIND), dimension(:),   pointer :: acc_depv_dt_diab_pv, acc_depv_dt_fric_pv, acc_depv_dt_dyn_pv

      ! Specific microphysics process tendencies 
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_evap_cw, depv_dt_mp_evap_rw, depv_dt_mp_depo_ice, depv_dt_mp_melt_ice
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_mp_evap_cw, acc_depv_dt_mp_evap_rw, acc_depv_dt_mp_depo_ice, acc_depv_dt_mp_melt_ice
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_frez_ice, depv_dt_mp_allproc 
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_mp_frez_ice, acc_depv_dt_mp_allproc 

      ! Latent heating tendencies 
      real(kind=RKIND), dimension(:,:), pointer :: dtheta_dt_cu, dtheta_dt_mp 
      real(kind=RKIND), dimension(:,:), pointer :: acc_dtheta_dt_cu, acc_dtheta_dt_mp

      ! Tendencies interpolated to isobaric levels 
      ! Note: all PV isobaric variables are defined in Registry_isobaric.xml
      ! Isobaric interpolation occurs in isobaric_diagnostics.F
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw_isobaric, depv_dt_sw_isobaric, depv_dt_bl_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_cu_isobaric, depv_dt_mp_isobaric, depv_dt_mix_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_fric_cu_isobaric, depv_dt_fric_bl_isobaric, depv_dt_fric_mix_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_diab_isobaric, depv_dt_fric_isobaric, depv_dt_dyn_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_evap_rw_isobaric, depv_dt_mp_evap_cw_isobaric, &
                                                   depv_dt_mp_depo_ice_isobaric, depv_dt_mp_melt_ice_isobaric, &
                                                   depv_dt_mp_frez_ice_isobaric, depv_dt_mp_allproc_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: dtheta_dt_cu_isobaric, dtheta_dt_mp_isobaric 
                                                                              
      ! Accumulated tendencies on isobaric levels 
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_lw_isobaric, acc_depv_dt_sw_isobaric, acc_depv_dt_bl_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_cu_isobaric, acc_depv_dt_mp_isobaric, acc_depv_dt_mix_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_fric_bl_isobaric, acc_depv_dt_fric_cu_isobaric, acc_depv_dt_fric_mix_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_diab_isobaric, acc_depv_dt_fric_isobaric, acc_depv_dt_dyn_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_mp_evap_rw_isobaric, acc_depv_dt_mp_evap_cw_isobaric, &
                                                   acc_depv_dt_mp_depo_ice_isobaric, acc_depv_dt_mp_melt_ice_isobaric, &
                                                   acc_depv_dt_mp_frez_ice_isobaric, acc_depv_dt_mp_allproc_isobaric
      real(kind=RKIND), dimension(:,:), pointer :: acc_dtheta_dt_cu_isobaric, acc_dtheta_dt_mp_isobaric 

      call mpas_pool_get_config(configs, 'config_pv_tend', config_pv_tend)
      call mpas_pool_get_config(configs, 'config_pv_microphys', config_pv_microphys)
      call mpas_pool_get_config(configs, 'config_pv_isobaric', config_pv_isobaric)

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      ! Instantaneous and accumulated PV tendencies 
      call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)
      call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)
      call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)
      call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)
      call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)
      call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)
      call mpas_pool_get_array(diag, 'depv_dt_fric_mix', depv_dt_fric_mix)
      call mpas_pool_get_array(diag, 'depv_dt_fric_bl', depv_dt_fric_bl)
      call mpas_pool_get_array(diag, 'depv_dt_fric_cu', depv_dt_fric_cu)
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      call mpas_pool_get_array(diag, 'depv_dt_dyn', depv_dt_dyn) 
      call mpas_pool_get_array(diag, 'depv_dt_diab_pv', depv_dt_diab_pv)
      call mpas_pool_get_array(diag, 'depv_dt_fric_pv', depv_dt_fric_pv)
      call mpas_pool_get_array(diag, 'depv_dt_dyn_pv', depv_dt_dyn_pv)

      call mpas_pool_get_array(diag, 'acc_depv_dt_lw', acc_depv_dt_lw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_sw', acc_depv_dt_sw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_bl', acc_depv_dt_bl)
      call mpas_pool_get_array(diag, 'acc_depv_dt_cu', acc_depv_dt_cu)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp', acc_depv_dt_mp)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mix', acc_depv_dt_mix)
      call mpas_pool_get_array(diag, 'acc_depv_dt_diab', acc_depv_dt_diab)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_mix', acc_depv_dt_fric_mix)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_bl', acc_depv_dt_fric_bl)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_cu', acc_depv_dt_fric_cu)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric', acc_depv_dt_fric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_dyn', acc_depv_dt_dyn)
      call mpas_pool_get_array(diag, 'acc_depv_dt_diab_pv', acc_depv_dt_diab_pv)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_pv', acc_depv_dt_fric_pv)
      call mpas_pool_get_array(diag, 'acc_depv_dt_dyn_pv', acc_depv_dt_dyn_pv)
      
      ! Specific microphysics process tendencies 
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_cw', depv_dt_mp_evap_cw)
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_rw', depv_dt_mp_evap_rw)
      call mpas_pool_get_array(diag, 'depv_dt_mp_depo_ice', depv_dt_mp_depo_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_melt_ice', depv_dt_mp_melt_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_frez_ice', depv_dt_mp_frez_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_allproc', depv_dt_mp_allproc)

      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_evap_cw', acc_depv_dt_mp_evap_cw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_evap_rw', acc_depv_dt_mp_evap_rw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_depo_ice', acc_depv_dt_mp_depo_ice)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_melt_ice', acc_depv_dt_mp_melt_ice)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_frez_ice', acc_depv_dt_mp_frez_ice)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_allproc', acc_depv_dt_mp_allproc)

      ! Latent heating tendencies 
      call mpas_pool_get_array(diag, 'dtheta_dt_cu', dtheta_dt_cu)
      call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)
      
      call mpas_pool_get_array(diag, 'acc_dtheta_dt_cu', acc_dtheta_dt_cu)
      call mpas_pool_get_array(diag, 'acc_dtheta_dt_mp', acc_dtheta_dt_mp)

      ! Isobaric tendencies -- interpolate prior to accumulating
      call mpas_pool_get_array(diag, 'depv_dt_lw_isobaric', depv_dt_lw_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_sw_isobaric', depv_dt_sw_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_bl_isobaric', depv_dt_bl_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_cu_isobaric', depv_dt_cu_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_isobaric', depv_dt_mp_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mix_isobaric', depv_dt_mix_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_fric_cu_isobaric', depv_dt_fric_cu_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_fric_bl_isobaric', depv_dt_fric_bl_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_fric_mix_isobaric', depv_dt_fric_mix_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_diab_isobaric', depv_dt_diab_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_fric_isobaric', depv_dt_fric_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_dyn_isobaric', depv_dt_dyn_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_rw_isobaric', depv_dt_mp_evap_rw_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_cw_isobaric', depv_dt_mp_evap_cw_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_depo_ice_isobaric', depv_dt_mp_depo_ice_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_melt_ice_isobaric', depv_dt_mp_melt_ice_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_frez_ice_isobaric', depv_dt_mp_frez_ice_isobaric)
      call mpas_pool_get_array(diag, 'depv_dt_mp_allproc_isobaric', depv_dt_mp_allproc_isobaric)
      call mpas_pool_get_array(diag, 'dtheta_dt_cu_isobaric', dtheta_dt_cu_isobaric)
      call mpas_pool_get_array(diag, 'dtheta_dt_mp_isobaric', dtheta_dt_mp_isobaric)
      
      ! Accumulated isobaric tendencies
      call mpas_pool_get_array(diag, 'acc_depv_dt_lw_isobaric', acc_depv_dt_lw_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_sw_isobaric', acc_depv_dt_sw_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_bl_isobaric', acc_depv_dt_bl_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_cu_isobaric', acc_depv_dt_cu_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_isobaric', acc_depv_dt_mp_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mix_isobaric', acc_depv_dt_mix_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_cu_isobaric', acc_depv_dt_fric_cu_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_bl_isobaric', acc_depv_dt_fric_bl_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_mix_isobaric', acc_depv_dt_fric_mix_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_diab_isobaric', acc_depv_dt_diab_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_isobaric', acc_depv_dt_fric_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_dyn_isobaric', acc_depv_dt_dyn_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_evap_rw_isobaric', acc_depv_dt_mp_evap_rw_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_evap_cw_isobaric', acc_depv_dt_mp_evap_cw_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_depo_ice_isobaric', acc_depv_dt_mp_depo_ice_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_melt_ice_isobaric', acc_depv_dt_mp_melt_ice_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_frez_ice_isobaric', acc_depv_dt_mp_frez_ice_isobaric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_allproc_isobaric', acc_depv_dt_mp_allproc_isobaric)
      call mpas_pool_get_array(diag, 'acc_dtheta_dt_cu_isobaric', acc_dtheta_dt_cu_isobaric)
      call mpas_pool_get_array(diag, 'acc_dtheta_dt_mp_isobaric', acc_dtheta_dt_mp_isobaric)

      ! Accumulate tendencies 
      acc_depv_dt_lw(:,:) = acc_depv_dt_lw(:,:) + depv_dt_lw(:,:)
      acc_depv_dt_sw(:,:) = acc_depv_dt_sw(:,:) + depv_dt_sw(:,:)
      acc_depv_dt_bl(:,:) = acc_depv_dt_bl(:,:) + depv_dt_bl(:,:)
      acc_depv_dt_cu(:,:) = acc_depv_dt_cu(:,:) + depv_dt_cu(:,:)
      acc_depv_dt_mp(:,:) = acc_depv_dt_mp(:,:) + depv_dt_mp(:,:)
      acc_depv_dt_mix(:,:) = acc_depv_dt_mix(:,:) + depv_dt_mix(:,:)
      acc_depv_dt_fric_bl(:,:) = acc_depv_dt_fric_bl(:,:) + depv_dt_fric_bl(:,:)
      acc_depv_dt_fric_cu(:,:) = acc_depv_dt_fric_cu(:,:) + depv_dt_fric_cu(:,:)
      acc_depv_dt_fric_mix(:,:) = acc_depv_dt_fric_mix(:,:) + depv_dt_fric_mix(:,:)
      acc_depv_dt_diab(:,:) = acc_depv_dt_diab(:,:) + depv_dt_diab(:,:)
      acc_depv_dt_fric(:,:) = acc_depv_dt_fric(:,:) + depv_dt_fric(:,:)
      acc_depv_dt_dyn(:,:) = acc_depv_dt_dyn(:,:) + depv_dt_dyn(:,:) 
      acc_depv_dt_diab_pv(:) = acc_depv_dt_diab_pv(:) + depv_dt_diab_pv(:)
      acc_depv_dt_fric_pv(:) = acc_depv_dt_fric_pv(:) + depv_dt_fric_pv(:)
      acc_depv_dt_dyn_pv(:)  = acc_depv_dt_dyn_pv(:)  + depv_dt_dyn_pv(:)
      
      ! Accumulate heating tendencies from microphysics and cumulus schemes 
      if (associated(dtheta_dt_mp)) then 
          acc_dtheta_dt_mp(:,:) = acc_dtheta_dt_mp(:,:) + dtheta_dt_mp(:,:)
      end if
      
      if (associated(dtheta_dt_cu)) then
          acc_dtheta_dt_cu(:,:) = acc_dtheta_dt_cu(:,:) + dtheta_dt_cu(:,:)
      end if 

      ! If specific microphysics process tendencies are requested, accumulate
      if (config_pv_microphys) then
          acc_depv_dt_mp_evap_cw(:,:) = acc_depv_dt_mp_evap_cw(:,:) + depv_dt_mp_evap_cw(:,:) 
          acc_depv_dt_mp_evap_rw(:,:) = acc_depv_dt_mp_evap_rw(:,:) + depv_dt_mp_evap_rw(:,:) 
          acc_depv_dt_mp_depo_ice(:,:) = acc_depv_dt_mp_depo_ice(:,:) + depv_dt_mp_depo_ice(:,:) 
          acc_depv_dt_mp_melt_ice(:,:) = acc_depv_dt_mp_melt_ice(:,:) + depv_dt_mp_melt_ice(:,:) 
          acc_depv_dt_mp_frez_ice(:,:) = acc_depv_dt_mp_frez_ice(:,:) + depv_dt_mp_frez_ice(:,:) 
          acc_depv_dt_mp_allproc(:,:) = acc_depv_dt_mp_allproc(:,:) + depv_dt_mp_allproc(:,:) 
      end if 

      ! If PV variables on isobaric levels are requested, accumulate the isobaric tendencies 
      ! NOTE: This requires specifications in ...... 
      
      if (config_pv_isobaric) then
          acc_depv_dt_lw_isobaric(:,:) = acc_depv_dt_lw_isobaric(:,:) + depv_dt_lw_isobaric(:,:)
          acc_depv_dt_sw_isobaric(:,:) = acc_depv_dt_sw_isobaric(:,:) + depv_dt_sw_isobaric(:,:)
          acc_depv_dt_bl_isobaric(:,:) = acc_depv_dt_bl_isobaric(:,:) + depv_dt_bl_isobaric(:,:)
          acc_depv_dt_cu_isobaric(:,:) = acc_depv_dt_cu_isobaric(:,:) + depv_dt_cu_isobaric(:,:)
          acc_depv_dt_mp_isobaric(:,:) = acc_depv_dt_mp_isobaric(:,:) + depv_dt_mp_isobaric(:,:)
          acc_depv_dt_mix_isobaric(:,:) = acc_depv_dt_mix_isobaric(:,:) + depv_dt_mix_isobaric(:,:)
          acc_depv_dt_fric_bl_isobaric(:,:) = acc_depv_dt_fric_bl_isobaric(:,:) + depv_dt_fric_bl_isobaric(:,:)
          acc_depv_dt_fric_cu_isobaric(:,:) = acc_depv_dt_fric_cu_isobaric(:,:) + depv_dt_fric_cu_isobaric(:,:)
          acc_depv_dt_fric_mix_isobaric(:,:) = acc_depv_dt_fric_mix_isobaric(:,:) + depv_dt_fric_mix_isobaric(:,:)
          acc_depv_dt_diab_isobaric(:,:) = acc_depv_dt_diab_isobaric(:,:) + depv_dt_diab_isobaric(:,:)
          acc_depv_dt_fric_isobaric(:,:) = acc_depv_dt_fric_isobaric(:,:) + depv_dt_fric_isobaric(:,:)
          acc_depv_dt_dyn_isobaric(:,:) = acc_depv_dt_dyn_isobaric(:,:) + depv_dt_dyn_isobaric(:,:)

          if (associated(dtheta_dt_cu_isobaric)) then
              acc_dtheta_dt_cu_isobaric(:,:) = acc_dtheta_dt_cu_isobaric(:,:) + dtheta_dt_cu_isobaric(:,:) 
          end if
          if (associated(dtheta_dt_mp_isobaric)) then
              acc_dtheta_dt_mp_isobaric(:,:) = acc_dtheta_dt_mp_isobaric(:,:) + dtheta_dt_mp_isobaric(:,:) 
          end if 

          ! If specific microphysics process tendencies are requested, accumulate the isobaric tendencies 
          if (config_pv_microphys) then
               acc_depv_dt_mp_evap_cw_isobaric(:,:) = acc_depv_dt_mp_evap_cw_isobaric(:,:) + depv_dt_mp_evap_cw_isobaric(:,:) 
               acc_depv_dt_mp_evap_rw_isobaric(:,:) = acc_depv_dt_mp_evap_rw_isobaric(:,:) + depv_dt_mp_evap_rw_isobaric(:,:) 
               acc_depv_dt_mp_depo_ice_isobaric(:,:) = acc_depv_dt_mp_depo_ice_isobaric(:,:) + depv_dt_mp_depo_ice_isobaric(:,:) 
               acc_depv_dt_mp_melt_ice_isobaric(:,:) = acc_depv_dt_mp_melt_ice_isobaric(:,:) + depv_dt_mp_melt_ice_isobaric(:,:) 
               acc_depv_dt_mp_frez_ice_isobaric(:,:) = acc_depv_dt_mp_frez_ice_isobaric(:,:) + depv_dt_mp_frez_ice_isobaric(:,:) 
               acc_depv_dt_mp_allproc_isobaric(:,:) = acc_depv_dt_mp_allproc_isobaric(:,:) + depv_dt_mp_allproc_isobaric(:,:) 
          end if  
      end if 

   end subroutine acc_pvBudget




!=================================================================================================================
end module pv_diagnostics
!=================================================================================================================
