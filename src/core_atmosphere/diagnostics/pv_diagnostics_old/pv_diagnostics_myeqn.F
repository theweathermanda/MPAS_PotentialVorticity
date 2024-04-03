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
    use mpas_kind_types, only : RKIND

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: tend
    type (MPAS_pool_type), pointer :: tend_physics
    type (MPAS_pool_type), pointer :: diag_physics
    type (MPAS_pool_type), pointer :: configs 

    type (MPAS_clock_type), pointer :: clock

    public :: pv_diagnostics_setup, &                          ! Order diagnostics subroutines are called in mpas_atm_core.F is setup, reset, update, compute, reset, etc. 
              pv_diagnostics_compute, &
              pv_diagnostics_reset, &
              pv_diagnostics_update , &                         
              pv_diagnostics_init                       ! MC: added to initialize ertel_pv_scalar as initial PV field and initialize dyn terms 

    private

    logical :: need_ertel_pv, need_u_pv, need_v_pv, need_theta_pv, need_vort_pv, need_iLev_DT, &
               need_advection_zonal, need_advection_merid, need_advection_vertical, need_wcell, &                  
               need_tend_lw, need_tend_sw, need_tend_bl, need_tend_cu, need_tend_mix, need_tend_mp, &
               need_tend_adv, need_tend_rho, &
               need_tend_diab, need_tend_fric, need_tend_diab_pv, need_tend_fric_pv, need_dtheta_mp, & 
               need_uzonal_prev, need_umerid_prev, need_wcell_prev, need_rho_prev, need_theta_prev, &                   
               need_rho_zz_prev, & ! MW added for decoupling friction term 
               need_pvvert_prev, need_epv_prev, need_iLev_prev, need_epv_scalar, need_tends_mp_proc, &
               need_acc_tend_lw, need_acc_tend_sw, need_acc_tend_bl, need_acc_tend_cu, need_acc_tend_mix, &
               need_acc_tend_mp, need_Acc_tend_diab, need_acc_tend_fric, need_acc_tend_diab_pv, need_acc_tend_fric_pv, &
               need_acc_tend_adv

    !=================================================================================================================
    ! MPAS code to compute Ertel's potential vorticity and tendency contributions
    ! Original diagnostics code was written by Nick Szapiro - 2016
    ! Significant changes to the code were made by Manda Chasteen (chasteen@ucar.edu) - 2022
    !
    ! subroutines contained in pv_diagnostics.F:
    ! ------------------------------------------
    ! pv_diagnostics_setup    :
    ! pv_diagnostics_update   :  
    ! pv_diagnostics_compute  :
    ! pv_diagnostics_reset    :


    ! Changes made by Manda Chasteen include:
    ! * Different formulation for calculation of horizontal gradients on native MPAS grid. The updated method is based on 
    !   Eq. 22 in Ringler et al. (2010) and is more robust than the previous method implemented by NS 
    ! * Reconstruction of horizontal gradients on each cell's edges to the cell center following the same method as the 
    !   horizontal wind reconstruction in mpas_vector_reconstruction.F
    ! * Incorporation of simple zonal, meridional, and vertical PV advection calculatons using the reconstructed PV gradients 
    ! * Changes to the calculation of the PV tendency terms to ensure that the potential temperature, wind, and density 
    !   fields from the beginning of the time step are used alongside the theta and wind tendencies in the calculations.
    !   In other words, the 3D gradient of the theta tendencies from a given time step should be dotted with the curl of 
    !   the UNMODIFIED 3D wind field at the beginning of the time step before any wind tendencies were added, and the curl 
    !   of the wind tendencies from a given time step should be dotted with the 3D gradient of the UNMODIFIED theta field
    !   at the beginning of the time step. This follows from the product rule.  
    !   This important change requires storing fields from the beginning of the time step to be used in the PV tendency 
    !   calculations because the model state and diagnostic fields are updated and assigned to time level 1 before the PV 
    !   diagnostics are called at the end of the time step in mpas_atm_core.F. Thus, before this change was implemented, 
    !   the updated variables from the end of the time step were incorrectly used alongside these tendencies.  
    ! * Added tend_u_phys(:,:) = 0._RKIND to mpas_atmphys_todynamics.F to initialize and update this variable as 0 for each
    !   time step. Without this modification, the tend_u_phys values would incorrectly accumulate over successive time steps.
    ! * Added PV tendencies for specific microphysical processes - net condensation/evaporation of cloud water, evaporation of 
    !   rain water, net deposition/sublimation, melting, and freezing
    !
    ! Modifications/additions needed:
    ! * Sum up tendencies over all over timesteps in outfile to compute averages over period? 
    ! * Incorporate PV into model integration so that it is advected like other state variables and/or make calculation
    !   of advection terms more sophisticated so the budget will balance?  


    contains

    !*********************************************************************************************************************
    ! pv_diagnostics_setup: initialize the PV diagnostics when called in mpas_atm_diagnostics_manager.F
    !*********************************************************************************************************************

    subroutine pv_diagnostics_setup(configs_in, all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type, MPAS_STREAM_OUTPUT, MPAS_STREAM_INPUT, &
                                       MPAS_STREAM_INPUT_OUTPUT
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none

        type (MPAS_pool_type), pointer :: configs_in
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'tend', tend)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)

        clock => simulation_clock

        configs => configs_in
  
    end subroutine pv_diagnostics_setup

    !*********************************************************************************************************************
    ! pv_diagnostics_init: 
    !                             A subroutine was created in mpas_atm_diagnostics_manager.F that calls this subroutine and 
    !                             is called explicitly in mpas_atm_core.F during the first time step so that values are 
    !                             updated after ertel_pv has been computed initially. 
    !                           - initializing the pv field to be advected by the model dynamics 
    !                           - initializing the *_dyn fields before the first time step, so ertel_pv_dyn is not NaN at time zero
    !*********************************************************************************************************************

    subroutine pv_diagnostics_init()
        use mpas_log, only : mpas_log_write

        implicit none

        call mpas_log_write('Doing pv_diagnostics_init()')
        ! Calculate PV diagnostics at each time step
        call mpas_log_write("Calling pv_diagnostics_scalar_init()")
        call pv_diagnostics_scalar_init()

        call mpas_log_write("Calling pv_diagnostics_dyn_init()")
        call pv_diagnostics_dyn_init() 

    end subroutine pv_diagnostics_init

    !*********************************************************************************************************************
    ! pv_diagnostics_scalar_init: initialize the pv field to be advected by the model dynamics 
    !                          - first called in pv_diagnostics_init at initial time
    !                          - subsequently called in pv_diagnostics_reset
    !*********************************************************************************************************************

    subroutine pv_diagnostics_scalar_init()
        use mpas_log, only : mpas_log_write
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none 

        integer, pointer :: index_ertel_pv_scalar
        real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
        real(kind=RKIND), dimension(:,:,:), pointer :: scalars

        call mpas_pool_get_dimension(state, 'index_pv_scalar', index_ertel_pv_scalar)  ! MW changed
        call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
        call mpas_pool_get_array(state, 'pv_scalars', scalars)  ! MW changed 

        call calc_epv(mesh, 1, state, diag) 

        scalars(index_ertel_pv_scalar,:,:) = ertel_pv(:,:)

        !call mpas_log_write("minval scalars(index_ertel_pv_scalar,:,:) is $r", realArgs=(/real(MINVAL(scalars(index_ertel_pv_scalar,:,:)))/))
        !call mpas_log_write("minval ertel_pv is $r", realArgs=(/real(MINVAL(ertel_pv))/)) 
        !print *, "minval ertel_pv", MINVAL(ertel_pv)


    end subroutine pv_diagnostics_scalar_init

    !*********************************************************************************************************************
    ! pv_diagnostics_dyn_init: initialize the _dyn fields that are used to compute ertel_pv_dyn 
    !                          - only called in pv_diagnostics_init at initial time
    !*********************************************************************************************************************

    subroutine pv_diagnostics_dyn_init()
        use mpas_log, only : mpas_log_write
        use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

        implicit none 

        real(kind=RKIND), dimension(:,:), pointer :: w_dyn, theta_dyn, u_dyn, rho_dyn, uReconstructZonal_dyn, uReconstructMeridional_dyn
        real(kind=RKIND), dimension(:,:), pointer :: w_1, theta_1, u_1, rho_1, uReconstructZonal_1, uReconstructMeridional_1

      call mpas_pool_get_array(diag, 'w_dyn', w_dyn)
      call mpas_pool_get_array(diag, 'theta_dyn', theta_dyn)
      call mpas_pool_get_array(diag, 'u_dyn', u_dyn) 
      call mpas_pool_get_array(diag, 'rho_dyn', rho_dyn) 
      call mpas_pool_get_array(diag, 'uReconstructZonal_dyn', uReconstructZonal_dyn)
      call mpas_pool_get_array(diag, 'uReconstructMeridional_dyn', uReconstructMeridional_dyn)

      call mpas_pool_get_array(state, 'w', w_1, 1)
      call mpas_pool_get_array(diag, 'theta', theta_1)
      call mpas_pool_get_array(state, 'u', u_1, 1) 
      call mpas_pool_get_array(diag, 'rho', rho_1) 
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal_1)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional_1)

      call mpas_log_write('Doing dyn init')

      w_dyn = w_1
      theta_dyn = theta_1 
      u_dyn = u_1
      rho_dyn = rho_1
      uReconstructZonal_dyn = uReconstructZonal_1
      uReconstructMeridional_dyn = uReconstructMeridional_1

      call mpas_log_write('Done with dyn init')
    end subroutine pv_diagnostics_dyn_init


    !*********************************************************************************************************************
    ! pv_diagnostics_update: compute the PV diagnostics over each time step when called in mpas_atm_diagnostics_manager.F
    !*********************************************************************************************************************

    subroutine pv_diagnostics_update()
        use mpas_log, only : mpas_log_write

        implicit none

        ! Calculate PV diagnostics at each time step
        call mpas_log_write("Computing Ertel's PV and PV advection")
        call atm_compute_pv_diagnostics(configs, state, 1, diag, mesh)
 
        call mpas_log_write('Computing PV tendency diagnostics.')
        call atm_compute_pvBudget_diagnostics(configs, state, 1, diag, mesh, tend, tend_physics, diag_physics)

    end subroutine pv_diagnostics_update

    !*********************************************************************************************************************
    ! pv_diagnostics_compute: compute PV diagnostic fields at the end of the timestep and/or write PV diagnostic fields to 
    !                         model output when called in mpas_atm_diagnostics_manager.F
    !
    ! MC - since PV diagnostics are now computed at each time step, subroutine was changed to only write fields to outfiles                            
    !*********************************************************************************************************************

    subroutine pv_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_log, only : mpas_log_write
     
        implicit none

        logical :: need_any_diags, need_any_budget, need_prev_fields

        need_any_diags = .false.
        need_any_budget = .false.

        ! Ertel's PV and fields interpolated to dynamic tropopause
        need_ertel_pv = MPAS_field_will_be_written('ertel_pv')
        need_any_diags = need_any_diags .or. need_ertel_pv
        need_u_pv = MPAS_field_will_be_written('u_pv')
        need_any_diags = need_any_diags .or. need_u_pv
        need_v_pv = MPAS_field_will_be_written('v_pv')
        need_any_diags = need_any_diags .or. need_v_pv
        need_theta_pv = MPAS_field_will_be_written('theta_pv')
        need_any_diags = need_any_diags .or. need_theta_pv
        need_vort_pv = MPAS_field_will_be_written('vort_pv')
        need_any_diags = need_any_diags .or. need_vort_pv
        need_iLev_DT = MPAS_field_will_be_written('iLev_DT')
        need_any_diags = need_any_diags .or. need_iLev_DT
        need_wCell = MPAS_field_will_be_written('wCell')
        need_any_diags = need_any_diags .or. need_wCell

        ! PV advection terms 
        need_advection_zonal = MPAS_field_will_be_written('advectionPV_zonal')                 ! MC added
        need_any_diags = need_any_diags .or. need_advection_zonal                              ! MC added
        need_advection_merid = MPAS_field_will_be_written('advectionPV_meridional')            ! MC added
        need_any_diags = need_any_diags .or. need_advection_merid                              ! MC added
        need_advection_vertical = MPAS_field_will_be_written('advectionPV_vertical')           ! MC added
        need_any_diags = need_any_diags .or. need_advection_vertical

        ! PV tendency terms 
        need_tend_lw = MPAS_field_will_be_written('depv_dt_lw')
        need_any_budget = need_any_budget .or. need_tend_lw
        need_tend_sw = MPAS_field_will_be_written('depv_dt_sw')
        need_any_budget = need_any_budget .or. need_tend_sw
        need_tend_bl = MPAS_field_will_be_written('depv_dt_bl')
        need_any_budget = need_any_budget .or. need_tend_bl
        need_tend_cu = MPAS_field_will_be_written('depv_dt_cu')
        need_any_budget = need_any_budget .or. need_tend_cu
        need_tend_mix = MPAS_field_will_be_written('depv_dt_mix')
        need_any_budget = need_any_budget .or. need_tend_mix
        need_dtheta_mp = MPAS_field_will_be_written('dtheta_dt_mp')
        need_any_budget = need_any_budget .or. need_dtheta_mp
        need_tend_mp = MPAS_field_will_be_written('depv_dt_mp')
        need_any_budget = need_any_budget .or. need_tend_mp
        need_tend_adv = MPAS_field_will_be_written('depv_dt_adv')
        need_any_budget = need_any_budget .or. need_tend_adv
        need_tend_diab = MPAS_field_will_be_written('depv_dt_diab')
        need_any_budget = need_any_budget .or. need_tend_diab
        need_tend_fric = MPAS_field_will_be_written('depv_dt_fric')
        need_any_budget = need_any_budget .or. need_tend_fric
        need_tend_rho = MPAS_field_will_be_written('depv_dt_rho')

        call mpas_log_write("wrote depv_dt_rho") 
        need_any_budget = need_any_budget .or. need_tend_rho
        need_tend_diab_pv = MPAS_field_will_be_written('depv_dt_diab_pv')
        need_any_budget = need_any_budget .or. need_tend_diab_pv
        need_tend_fric_pv = MPAS_field_will_be_written('depv_dt_fric_pv')
        need_any_budget = need_any_budget .or. need_tend_fric_pv


        ! Accumulated PV tendency terms
        need_acc_tend_lw = MPAS_field_will_be_written('acc_depv_dt_lw')
        need_any_budget = need_any_budget.or. need_acc_tend_lw
        need_acc_tend_sw = MPAS_field_will_be_written('acc_depv_dt_sw')
        need_any_budget = need_any_budget.or. need_acc_tend_sw
        need_acc_tend_bl = MPAS_field_will_be_written('acc_depv_dt_bl')
        need_any_budget = need_any_budget.or. need_acc_tend_bl
        need_acc_tend_cu = MPAS_field_will_be_written('acc_depv_dt_cu')
        need_any_budget = need_any_budget.or. need_acc_tend_cu
        need_acc_tend_mp = MPAS_field_will_be_written('acc_depv_dt_mp')
        need_any_budget = need_any_budget.or. need_acc_tend_mp
        need_acc_tend_mix = MPAS_field_will_be_written('acc_depv_dt_mix')
        need_any_budget = need_any_budget.or. need_acc_tend_mix
        need_acc_tend_diab = MPAS_field_will_be_written('acc_depv_dt_diab')
        need_any_budget = need_any_budget.or. need_acc_tend_diab
        need_acc_tend_fric = MPAS_field_will_be_written('acc_depv_dt_fric')
        need_any_budget = need_any_budget.or. need_acc_tend_fric
        need_acc_tend_adv = MPAS_field_will_be_written('acc_depv_dt_adv')
        need_any_budget = need_any_budget.or. need_acc_tend_adv
        need_acc_tend_diab_pv = MPAS_field_will_be_written('acc_depv_dt_diab_pv')
        need_any_budget = need_any_budget.or. need_acc_tend_diab_pv
        need_acc_tend_fric_pv = MPAS_field_will_be_written('acc_depv_dt_fric_pv')
        need_any_budget = need_any_budget.or. need_acc_tend_fric_pv


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

        need_any_budget = need_any_budget .or. need_tends_mp_proc

        ! MC - If fields do not need to be computed each time step, then uncomment the lines below and 
        ! comment the same lines in pv_diagnostics_update above
        ! Caution: only make this change if the PV tendency variables are not needed because they require
        ! stored fields from the previous time step!         
      
        ! Calculate PV diagnostics before outfiles are written at end of each time step
        !if (need_any_diags) then
            !call mpas_log_write("Computing Ertel's PV and PV advection")
            !call atm_compute_pv_diagnostics(state, 1, diag, mesh)
        !end if
        !if (need_any_budget) then
            !call mpas_log_write('Computing PV tendency diagnostics.')
            !call atm_compute_pvBudget_diagnostics(state, 1, diag, mesh, tend, tend_physics)
        !end if

        ! Save fields from beginning of time step to outfile before updating them 
        call mpas_log_write('Saving beginning of time step variables to outfile.')

        need_uzonal_prev = MPAS_field_will_be_written('uReconstructZonal_prev')
        need_umerid_prev = MPAS_field_will_be_written('uReconstructMeridional_prev')
        need_wcell_prev = MPAS_field_will_be_written('wCell_prev')
        need_theta_prev = MPAS_field_will_be_written('theta_prev')
        need_rho_prev = MPAS_field_will_be_written('rho_prev')
        need_rho_zz_prev = MPAS_field_will_be_written('rho_zz_prev')  ! MW added for decoupling friction term
        need_pvvert_prev = MPAS_field_will_be_written('pv_vertex_prev')
        need_epv_prev = MPAS_field_will_be_written('ertel_pv_prev')
        need_iLev_prev = MPAS_field_will_be_written('iLev_DT_prev')

        ! Write ertel_pv_scalar 
        need_epv_scalar = MPAS_field_will_be_written('ertel_pv_scalar')

    end subroutine pv_diagnostics_compute

    !*********************************************************************************************************************
    ! pv_diagnostics_reset: Update beginning of time step fields after they have been written to outfile for use in next 
    !                       time step PV tendency calculations when called in mpas_atm_diagnostics_manager.F
    !*********************************************************************************************************************

    subroutine pv_diagnostics_reset()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_log, only : mpas_log_write

        implicit none

        ! Called immediately after diagnostics have been written
        ! Update previous variables for next timestep calculation
        call mpas_log_write('Doing pv_diagnostics_reset()')

        call mpas_log_write('Updating previous fields.')
        call store_previous_vars(mesh, 2, state, diag)

        ! Update ertel_pv_scalar to be calculated ertel_pv for next advection cycle.
!        print *, "calling init to reset pv scalar variable"
        call pv_diagnostics_scalar_init() 

        ! Populate the fields used to calculate ertel_pv_dyn with default values from time-level 1 
        call pv_diagnostics_dyn_init()

    end subroutine pv_diagnostics_reset
     
    !*********************************************************************************************************************
    ! NS: Below are two subroutines (floodFill_strato and floodFill_tropo), designed to determine the first model level 
    ! above the dynamic tropopause, iLev_DT, which is designated as the 2-PVU isosurface. Only one of these subroutines 
    ! is used (toggled with "call floodFill_strato(mesh, diag, pvuVal, stratoPV)" and "call floodFill_tropo(mesh,diag,pvuVal)"
    ! in the atm_compute_pv_diagnostics subroutine below. The routines *should* produce equivalent estimates for iLev_DT. 
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
 
       use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
       use mpas_dmpar, only : mpas_dmpar_max_int,mpas_dmpar_exch_halo_field
       use mpas_derived_types, only : dm_info, field2DInteger
      
       implicit none
       
       type (mpas_pool_type), intent(in) :: mesh
       type (mpas_pool_type), intent(inout) :: diag
       real(kind=RKIND), intent(in) :: pvuVal
 
       integer :: iCell, k, nChanged, iNbr, iCellNbr, levInd, haloChanged, global_haloChanged
       integer, pointer :: nCells, nVertLevels
       integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
       integer, dimension(:,:), pointer :: cellsOnCell, inTropo
 
       type (field2DInteger), pointer :: inTropo_f
 
       real(kind=RKIND) :: sgnHemi, sgn
       real(kind=RKIND),dimension(:),pointer:: latCell
       real(kind=RKIND), dimension(:,:), pointer :: ertel_pv
       
       type (dm_info), pointer :: dminfo
 
       integer, dimension(:,:), allocatable :: candInTropo !whether in troposphere
       
       call mpas_pool_get_dimension(mesh, 'nCells', nCells)
       call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
       call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
       call mpas_pool_get_array(mesh, 'cellsOnCell', cellsOnCell)
       call mpas_pool_get_array(mesh, 'latCell', latCell)
 
       call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
       !call mpas_pool_get_array(diag, 'iLev_DT_trop', iLev_DT)
       call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
       call mpas_pool_get_array(diag, 'inTropo', inTropo)
       
       allocate(candInTropo(nVertLevels, nCells+1))
       candInTropo(:,:) = 0
       inTropo(:,:) = 0
       !store whether each level above DT to avoid repeating logic. we'll
       !use cand as a isVisited marker further below.
       do iCell=1,nCells
          sgnHemi = sign(1.0_RKIND, latCell(iCell)) !at the equator, sign(0)=0
          if (sgnHemi .EQ. 0.0) sgnHemi = 1.0_RKIND
          do k=1,nVertLevels
             sgn = ertel_pv(k,iCell)*sgnHemi-pvuVal
             if (sgn .LT. 0) candInTropo(k,iCell) = 1
          end do
       end do
       
       !seed flood fill with near surface that's below DT (can have
       !surface above 2pvu from pv anoms).
       !Note that this would be wrong if low PV "stratospheric" blobs are
       !right above the surface
       nChanged = 0
       levInd = min(nVertLevels, 3)
       do iCell=1,nCells
          do k=1,levInd
             if (candInTropo(k,iCell) .GT. 0) then
                inTropo(k,iCell) = 1
                !candInTropo(k,iCell) = 0
                nChanged = nChanged+1
             end if
          end do
       end do
       
       !flood fill from the given seeds. since I don't know enough
       !fortran,
       !we'll just brute force a continuing loop rather than queue.
       call mpas_pool_get_field(diag, 'inTropo', inTropo_f)
       dminfo => inTropo_f % block % domain % dminfo
       global_haloChanged = 1
       do while(global_haloChanged .GT. 0) !any cell in a halo has changed, to propagate to other domains
         global_haloChanged = 0 !aggregate the number of changed cells w/in the loop below
         do while(nChanged .GT. 0)
           nChanged = 0
           do iCell=1,nCells !should we look for neighbors of halo cells?
             do k=1,nVertLevels
                !update if candidate and neighbor in troposphere
                if ((candInTropo(k,iCell) .GT. 0) .AND. (inTropo(k,iCell).LT.1) ) then
                   !nbr below
                   if (k .GT. 1) then
                     if (inTropo(k-1,iCell) .GT. 0) then
                       inTropo(k,iCell) = 1
                       !candInTropo(k,iCell) = 0
                       nChanged = nChanged+1
                       cycle
                     end if
                   end if
 
                   !side nbrs
                   do iNbr = 1, nEdgesOnCell(iCell)
                     iCellNbr = cellsOnCell(iNbr,iCell)
                     if (inTropo(k,iCellNbr) .GT. 0) then
                       inTropo(k,iCell) = 1
                       !candInTropo(k,iCell) = 0
                       nChanged = nChanged+1
                       exit
                     end if
                   end do
 
                   !nbr above
                   if (k .LT. nVertLevels) then
                     if (inTropo(k+1,iCell) .GT. 0) then
                       inTropo(k,iCell) = 1
                       !candInTropo(k,iCell) = 0
                       nChanged = nChanged+1
                       cycle
                     end if
                   end if
 
                end if !candIn
             end do !levels
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
       deallocate(candInTropo)
       
       !Fill iLev_DT with the lowest level above the tropopause (If DT
       !above column, iLev>nVertLevels. If DT below column, iLev=0.
       do iCell=1,nCells
         nChanged = 0
         do k=nVertLevels,1,-1
           if (inTropo(k,iCell) .GT. 0) then
             nChanged = 1
             exit
           end if
         end do !k
         
         if (nChanged .GT. 0) then !found troposphere's highest level
           iLev_DT(iCell) = k+1 !level above troposphere (>nVertLevels if whole column below 2pvu; e.g., tropics)
         else !whole column above DT (e.g., arctic pv tower)
           iLev_DT(iCell) = 0
         end if
       end do !iCell
      
    end subroutine floodFill_tropo

    !*********************************************************************************************************************
    ! NS: Subroutine to compute various fields on 2-PVU surface using the calculated PV field - potential temperature, 
    !     uZonal, uMeridional, vertical vorticity
    ! MC: Modified interpolation of vorticity to cell centers procedure to be consistent with other changes  
    !     This routine should use PV field at end of time step! 
    !********************************************************************************************************************* 
  
    subroutine interp_pv_diagnostics(mesh, diag, pvuVal, missingVal)
 
       use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
       
       IMPLICIT NONE
       
       type (mpas_pool_type), intent(in) :: mesh
       type (mpas_pool_type), intent(inout) :: diag
       real(kind=RKIND) :: pvuVal, missingVal
       
       integer :: iCell, k
       integer, pointer :: nCells, nVertLevels
       integer, dimension(:), pointer :: nEdgesOnCell, iLev_DT
       integer, dimension(:,:), pointer :: cellsOnCell, cellsOnEdge, verticesOnCell, cellsOnVertex 
                                           
       real(kind=RKIND), dimension(:), pointer:: areaCell, latCell, u_pv, v_pv, theta_pv, vort_pv
       real(kind=RKIND), dimension(:,:), pointer:: uReconstructZonal, uReconstructMeridional, vorticity, theta, ertel_pv, &
                                                   kiteAreasOnVertex
       real(kind=RKIND), dimension(:,:), allocatable :: vVort
       
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
       
       call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
       call mpas_pool_get_array(diag, 'theta', theta)
       call mpas_pool_get_array(diag, 'vorticity', vorticity)
       call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
       call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
       call mpas_pool_get_array(diag, 'u_pv', u_pv)
       call mpas_pool_get_array(diag, 'v_pv', v_pv)
       call mpas_pool_get_array(diag, 'theta_pv', theta_pv)
       call mpas_pool_get_array(diag, 'vort_pv', vort_pv)
       call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)

       allocate(vVort(nVertLevels,nCells))                         ! MC: Nick's code and other parts of model allocate variables 
                                                                   ! with dims of nCells+1... not sure why?
      
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, uReconstructZonal, u_pv, missingVal, iLev_DT)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, uReconstructMeridional, v_pv, missingVal, iLev_DT)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv, theta, theta_pv, missingVal, iLev_DT)

       call interp_absVertVort(vorticity, nCells, nEdgesOnCell, verticesOnCell, &
                               cellsOnVertex, areaCell, kiteAreasOnVertex, vVort)
               
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, ertel_pv, vVort, &
                      vort_pv, missingVal, iLev_DT)

       deallocate(vVort)

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
                                           
       real(kind=RKIND), dimension(:),pointer :: latCell, depv_dt_diab_pv, depv_dt_fric_pv
       real(kind=RKIND), dimension(:,:),pointer :: depv_dt_diab, depv_dt_fric, ertel_pv_prev      ! MC changed
      
       call mpas_pool_get_dimension(mesh, 'nCells', nCells) 
       call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
       
       call mpas_pool_get_array(mesh, 'latCell', latCell)
       call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev)                             ! MC changed
       call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
       call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
       call mpas_pool_get_array(diag, 'depv_dt_diab_pv', depv_dt_diab_pv)
       call mpas_pool_get_array(diag, 'depv_dt_fric_pv', depv_dt_fric_pv)
       call mpas_pool_get_array(diag, 'iLev_DT_prev', iLev_DT_prev)                               ! MC changed
       
       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv_prev, depv_dt_diab, depv_dt_diab_pv, missingVal, iLev_DT_prev)

       call interp_pv(nCells, nVertLevels, pvuVal, latCell, &
                      ertel_pv_prev, depv_dt_fric, depv_dt_fric_pv, missingVal, iLev_DT_prev)

    end subroutine interp_pvBudget_diagnostics
   
    !*********************************************************************************************************************
    ! NS: Subroutine to linearly interpolate columns of field1 to where field0 is interpVal*sign(lat) using level above 
    !     tropopause already diagnosed 
    !********************************************************************************************************************* 

    subroutine interp_pv(nCells, nLevels, interpVal, latCell, field0, &
                           field1, field, missingVal, iLev_DT)
        
       IMPLICIT NONE  
 
       integer :: nCells, nLevels
       integer, intent(in) :: iLev_DT(nCells)
       real(kind=RKIND) :: interpVal, missingVal
       real(kind=RKIND), intent(in) :: latCell(nCells)
       real(kind=RKIND), intent(in) :: field0(nLevels,nCells), field1(nLevels,nCells)
       real(kind=RKIND), intent(out) :: field(nCells)
 
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
         
         iLev = iLev_DT(iCell)
         if (iLev .GT. nLevels) then
           levInd = -1
           sgnl = -1.0
         else if (iLev .LT. 1) then
           levInd = -1
           sgnl = 1.0
         else
           valh = field0(iLev,iCell)
           vall = field0(iLev-1,iCell)
           !sandwiched value. equal in case val0 is a vals[l].
           !get linear interpolation: val0 = vals[l]+dvals/dl * dl
           !Avoid divide by 0 by just assuming value is 
           !halfway between...
    
           dv_dl = valh-vall;
           if (abs(dv_dl)<1.e-6) then
             levFrac = 0.5;
           else
             levFrac = (valInterpCell-vall)/dv_dl
           end if
           
           levInd = iLev-1
         end if !iLev in column
 
         !find value of field using index we just found
         if (levInd<0) then !didn't trap value
           if (sgnl>0.0) then !column above value, take surface
             field(iCell) = field1(1,iCell)
           else !column below value, take top
             !field(iCell) = missingVal
             field(iCell) = field1(nLevels,iCell)
           end if
         else
           valh = field1(levInd+1,iCell)
           vall = field1(levInd,iCell)
         
           dv_dl = valh-vall
           field(iCell) = vall+dv_dl*levFrac
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
    ! MW: Test code - adding this to recompute pv_vertex using the velocities after dynamics
    !********************************************************************************************************************* 
    subroutine recompute_pv_vertex(u, nVertices, nVertLevels, vertexDegree, invAreaTriangle, &
                                 dcEdge, edgesOnVertex, edgesOnVertex_sign, fVertex, vort )
 
       IMPLICIT NONE
 
       integer, intent(in) :: nVertices, vertexDegree, nVertLevels
       integer, dimension(:,:), intent(in) :: edgesOnVertex
       real(kind=RKIND), dimension(:,:), intent(in) :: u, edgesOnVertex_sign 
       real(kind=RKIND), dimension(:), intent(in) :: invAreaTriangle, fVertex, dcEdge
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
 
    end subroutine recompute_pv_vertex 


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
       integer, dimension(:,:), intent(in) :: cellsOnEdge, edgesOnCell
       integer, dimension(:), intent(in) :: nEdgesOnCell
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
   ! MC: Modified original curl subroutine to include calculation over all vertical levels
   !
   ! NS: Adapted from computation of circulation and relative vorticity at each vertex in atm_compute_solve_diagnostics()
   !     This takes scvt face values and computes finite volume curl at scvt vertices (triangle cell centers)
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
   !*********************************************************************************************************************

   subroutine calc_vertDeriv(var, nCells, nVertLevels, dzu, dvar_dz)

      IMPLICIT NONE

      integer, intent(in) :: nCells, nVertLevels
      real(kind=RKIND), dimension(:), intent(in) ::  dzu
      real(kind=RKIND), dimension(:,:), intent(in) :: var
      real(kind=RKIND), dimension(:,:), intent(out) :: dvar_dz
      integer :: iCell, k
      real(kind=RKIND) :: dvar_dz_top, dvar_dz_bot

      dvar_dz(:,:) = 0.0_RKIND

      do iCell=1,nCells
         dvar_dz(1,iCell) = (var(2,iCell) - var(1,iCell)) / dzu(2)
         dvar_dz(nVertLevels,iCell) = (var(nVertLevels,iCell) - var(nVertLevels-1,iCell)) / &
                                       dzu(nVertLevels)
         do k=2,nVertLevels-1
            dvar_dz_top = (var(k+1,iCell) - var(k,iCell)) / dzu(k+1)
            dvar_dz_bot = (var(k,iCell) - var(k-1,iCell)) / dzu(k)
            dvar_dz(k,iCell) = 0.5 * (dvar_dz_top + dvar_dz_bot)
         end do
      end do
      
   end subroutine calc_vertDeriv

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
   ! MC: Subroutine to calculate density tendency using density difference between current and previous time steps
   !*********************************************************************************************************************

   subroutine calc_rho_tend(rho, rho_prev, config_dt, rho_tend)

      IMPLICIT NONE

      real(kind=RKIND), dimension(:,:), intent(in) :: rho, rho_prev
      real(kind=RKIND), intent(in) :: config_dt
      real(kind=RKIND), dimension(:,:), intent(out) :: rho_tend

      rho_tend(:,:) = 0.0_RKIND
      rho_tend = (rho - rho_prev) / config_dt

   end subroutine calc_rho_tend




   !*********************************************************************************************************************   
   ! MC: Subroutine to store variables from the beginning of the time step to use in next timestep tendency calculations 
   !*********************************************************************************************************************

   subroutine store_previous_vars(mesh, time_lev, state, diag)

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      use mpas_log, only : mpas_log_write
 
      IMPLICIT NONE
   
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev                                ! which time level to use from state: 1 = beginning of time step; 2 = end of time step
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh

      integer, pointer :: nCells, nVertLevels, nEdges, nVertices

      real(kind=RKIND), dimension(:,:), pointer :: uReconstructZonal, uReconstructMeridional, wCell, theta, rho, &
                                                   pv_vertex, ertel_pv, rho_edge
      real(kind=RKIND), dimension(:,:), pointer :: uReconstructZonal_prev, uReconstructMeridional_prev, wCell_prev, &
                                                   theta_prev, rho_prev, pv_vertex_prev, ertel_pv_prev, &
                                                   rho_edge_prev
      real(kind=RKIND), dimension(:,:), pointer :: rho_zz, rho_zz_prev ! MW: added for decoupling friction term

      integer, dimension(:), pointer :: iLev_DT, iLev_DT_prev

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)

      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'wCell', wCell)
      call mpas_pool_get_array(diag, 'theta', theta)
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(state, 'rho_zz', rho_zz, 1)   ! MW added for decoupling friction term
      call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'iLev_DT', iLev_DT)
      call mpas_pool_get_array(diag, 'rho_edge', rho_edge)  ! MW: this is rho_zz at edge

      call mpas_pool_get_array(diag, 'uReconstructZonal_prev', uReconstructZonal_prev)
      call mpas_pool_get_array(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev)
      call mpas_pool_get_array(diag, 'wCell_prev', wCell_prev)
      call mpas_pool_get_array(diag, 'theta_prev', theta_prev)
      call mpas_pool_get_array(diag, 'rho_prev', rho_prev)      
      call mpas_pool_get_array(diag, 'rho_zz_prev', rho_zz_prev)  ! MW added for decoupling friction term 
      call mpas_pool_get_array(diag, 'pv_vertex_prev', pv_vertex_prev)
      call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev)
      call mpas_pool_get_array(diag, 'iLev_DT_prev', iLev_DT_prev)
      call mpas_pool_get_array(diag, 'rho_edge_prev', rho_edge_prev)
      !call mpas_pool_get_array(diag, 'uReconstructZonal_prev', uReconstructZonal_prev, time_lev)
      !call mpas_pool_get_array(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev, time_lev)
      !call mpas_pool_get_array(diag, 'wCell_prev', wCell_prev, time_lev)
      !call mpas_pool_get_array(diag, 'theta_prev', theta_prev, time_lev)
      !call mpas_pool_get_array(diag, 'rho_prev', rho_prev, time_lev)      
      !call mpas_pool_get_array(diag, 'rho_zz_prev', rho_zz_prev, time_lev)  ! MW added for decoupling friction term 
      !call mpas_pool_get_array(diag, 'pv_vertex_prev', pv_vertex_prev, time_lev)
      !call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev, time_lev)
      !call mpas_pool_get_array(diag, 'iLev_DT_prev', iLev_DT_prev, time_lev)
      !call mpas_pool_get_array(diag, 'rho_edge_prev', rho_edge_prev, time_lev)

      !call mpas_log_write("maxval ertel_pv_prev before resetting is $r", realArgs=(/real(MAXVAL(ertel_pv_prev), kind=RKIND)/))

      uReconstructZonal_prev(:,:) = uReconstructZonal(:,:)
      uReconstructMeridional_prev(:,:) = uReconstructMeridional(:,:)
      wCell_prev(:,:) = wCell(:,:)
      theta_prev(:,:) = theta(:,:)
      !call mpas_log_write("in setting _prev variables, theta_prev = $r, theta = $r", realArgs=(/ theta_prev(55,100), theta(55,100)/))
      !call mpas_log_write("in setting _prev variables, wCell_prev = $r, wCell = $r", realArgs=(/ wCell_prev(55,100), wCell(55,100)/))
      rho_prev(:,:) = rho(:,:)
      rho_zz_prev(:,:) = rho_zz(:,:) ! MW: added for decoupling friction term
      pv_vertex_prev(:,:) = pv_vertex(:,:)
      !call mpas_log_write("in setting _prev variables, ertel_pv_prev = $r, ertel_pv = $r", realArgs=(/ ertel_pv_prev(55,100), ertel_pv(55,100)/))
      ertel_pv_prev(:,:) = ertel_pv(:,:)
      !call mpas_log_write("in setting _prev variables, ertel_pv_prev = $r, ertel_pv = $r", realArgs=(/ ertel_pv_prev(55,100), ertel_pv(55,100)/))
      iLev_DT_prev(:) = iLev_DT(:)
      rho_edge_prev(:,:) = rho_edge(:,:)

      !call mpas_log_write("maxval ertel_pv_prev after resetting is $r", realArgs=(/real(MAXVAL(ertel_pv_prev), kind=RKIND)/))

   end subroutine store_previous_vars

   !*********************************************************************************************************************
   ! MC: Modified subroutine to calculate Ertel's potential vorticity
   !     PV = 1/density * [curl(wind) . grad(theta)] 
   !*********************************************************************************************************************

   subroutine calc_epv(mesh, time_lev, state, diag)
      use mpas_log, only : mpas_log_write ! MW debug

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

      IMPLICIT NONE

      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev                     ! which time level to use from state: 1 = beginning of time step; 2 = end of time step
      type (mpas_pool_type), intent(inout) :: diag        ! MC: since these subroutines are called after state/diag vars have been updated, time_lev 
      type (mpas_pool_type), intent(in) :: mesh           !     specification here doesn't really matter

      ! input vars
      integer, pointer :: nCells, nVertLevels, nEdges, R3
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex

      real(kind=RKIND), dimension(:), pointer :: dzu, areaCell, latCell, lonCell, dcEdge
      real(kind=RKIND), dimension(:,:), pointer :: w, rho, theta, pv_vertex, zgrid, kiteAreasOnVertex, &
                                                   uReconstructZonal, uReconstructMeridional, ertel_pv, & 
                                                   edgesOnCell_sign, wCell
      real(kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct 

      ! local vars
      real(kind=RKIND), dimension(:,:), allocatable :: duZonal_dz, duMerid_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dTheta_dxZonal, dTheta_dyMerid, dTheta_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dW_dxZonal, dW_dyMerid 
      real(kind=RKIND), dimension(:,:), allocatable :: absVort
      real(kind=RKIND), dimension(:,:,:), allocatable :: absVort3D, gradTheta

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'R3', R3)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'coeffs_reconstruct', coeffs_reconstruct)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(state, 'w', w, time_lev)     
      call mpas_pool_get_array(diag, 'theta', theta)  
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
      call mpas_pool_get_array(diag, 'wCell', wCell)    

      ! Allocate memory to intermediate vars 
      allocate(absVort(nVertLevels,nCells)) 
      allocate(duZonal_dz(nVertLevels,nCells))
      allocate(duMerid_dz(nVertLevels,nCells))
      allocate(dTheta_dxZonal(nVertLevels,nCells))
      allocate(dTheta_dyMerid(nVertLevels,nCells))
      allocate(dTheta_dz(nVertLevels,nCells))
      allocate(dW_dxZonal(nVertLevels,nCells))
      allocate(dW_dyMerid(nVertLevels,nCells))
      allocate(absVort3D(nVertLevels,nCells,3))
      allocate(gradTheta(nVertLevels,nCells,3))

      ertel_pv(:,:) = 0.0_RKIND
      gradTheta(:,:,:) = 0.0_RKIND
      absVort3D(:,:,:) = 0.0_RKIND

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D potential temperature gradient 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !call mpas_log_write('calc_epv: theta(k=55,iCell=100) = $r',realArgs=(/ theta(55,100) /) )
      !call mpas_log_write('calc_epv: uReconstructZonal(k=55,iCell=100) = $r',realArgs=(/ uReconstructZonal(55,100) /) )
      !call mpas_log_write('calc_epv: uReconstructMeridional(k=55,iCell=100) = $r',realArgs=(/ uReconstructMeridional(55,100) /) )


      ! (1) Calculate and reconstruct horizontal potential temperature gradient to get zonal and meridional
      !     gradients at cell centers: dth_dx, dth_dy

      call calc_gradOnEdges_reconCellCenter(theta, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dTheta_dxZonal, dTheta_dyMerid)


      ! (2) Calculate the vertical potential temperature gradient: dth_dz

      call calc_vertDeriv(theta, nCells, nVertLevels, dzu, dTheta_dz)

      ! (3) Combine theta derivatives into 3D vector

      gradTheta(:,:,1) = dTheta_dxZonal
      gradTheta(:,:,2) = dTheta_dyMerid
      gradTheta(:,:,3) = dTheta_dz

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D absolute vorticity vector  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate the vertical shear of uReconstructZonal and uReconstructMeridional: du_dz and dv_dz

      call calc_vertDeriv(uReconstructZonal, nCells, nVertLevels, dzu, duZonal_dz)
      call calc_vertDeriv(uReconstructMeridional, nCells, nVertLevels, dzu, duMerid_dz)
      !call mpas_log_write('calc_epv: duZonal_dz(k=55,iCell=100) = $r',realArgs=(/ duZonal_dz(55,100) /) )
      !call mpas_log_write('calc_epv: duMerid_dz(k=55,iCell=100) = $r',realArgs=(/ duMerid_dz(55,100) /) )

      ! (2) Interpolate w to cell centers, calculate gradient of w on edges, and then reconstruct to get 
      !     zonal and meridional gradients at cell centers: dw_dx, dw_dy 

      call interp_wLev_thetaLev(w, nCells, nVertLevels, wCell)
      !call mpas_log_write('calc_epv: wCell(k=55,iCell=100) = $r',realArgs=(/ wCell(55,100) /) )

      call calc_gradOnEdges_reconCellCenter(wCell, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dW_dxZonal, dW_dyMerid)

      ! (3) Recconstruct absolute vertical vorticity at vertices pv_vertex to cell centers

      !     Note: currently, pv_vertex is the absolute vertical vorticity on the cell vertices. If this 
      !     variable changes at some point, then the absolute vertical vorticity on the vertices needs 
      !     to be computed as follows: 
      !     do iVert=1,nVertices
      !        vorticity(:,iVert) = vorticity(:,iVert) + fVertex(iVert)
      !     end do

      call interp_absVertVort(pv_vertex, nCells, nEdgesOnCell, verticesOnCell, &
                                  cellsOnVertex, areaCell, kiteAreasOnVertex, absVort)
      !call mpas_log_write('calc_epv: pv_vertex(k=55,iCell=100) = $r',realArgs=(/ pv_vertex(55,100) /) )

      ! (4) Combine three components into vorticity vector

      absVort3D(:,:,1) = dW_dyMerid - duMerid_dz                   ! dw/dy - dv/dz
      absVort3D(:,:,2) = duZonal_dz - dW_dxZonal                   ! du/dz - dw/dx
      absVort3D(:,:,3) = absVort                                   ! dv/dy - du/dx + f

      ! (5) Take dot product between 3D theta vector and absolute vorticity / density to compute 
      !     Ertel's PV
  

      !call mpas_log_write('calc_epv: absVort3D components (k=55,iCell=100) = $r $r $r',realArgs=(/ absVort3D(55,100,1), absVort3D(55,100,2), absVort3D(55,100,3) /) )
      !call mpas_log_write('calc_epv: gradTheta components (k=55,iCell=100) = $r $r $r',realArgs=(/ gradTheta(55,100,1), gradTheta(55,100,2), gradTheta(55,100,3) /) )
      call calc_dotProduct_3D(gradTheta, absVort3D, nCells, nVertLevels, ertel_pv)
 
      !call mpas_log_write('calc_epv: ertel_pv(k=55,iCell=100) = $r',realArgs=(/ ertel_pv(55,100) /) )
      !call mpas_log_write('calc_epv: rho(k=55,iCell=100) = $r',realArgs=(/ rho(55,100) /) )
      ertel_pv = ertel_pv / rho * 1.0e6                            !SI to PVUs
      !call mpas_log_write('calc_epv: ertel_pv(k=55,iCell=100) = $r',realArgs=(/ ertel_pv(55,100) /) )
      !call mpas_log_write('calc_epv: rho(k=55,iCell=100) = $r',realArgs=(/ rho(55,100) /) )

      deallocate(absVort)
      deallocate(absVort3D)
      deallocate(gradTheta)
      deallocate(duZonal_dz)
      deallocate(duMerid_dz)
      deallocate(dTheta_dz)
      deallocate(dTheta_dxZonal)
      deallocate(dTheta_dyMerid)
      deallocate(dW_dxZonal)
      deallocate(dW_dyMerid)

   end subroutine calc_epv

   !*********************************************************************************************************************
   ! MC: Modified subroutine to calculate Ertel's potential vorticity
   !     PV = 1/density * [curl(wind) . grad(theta)] 
   !*********************************************************************************************************************

   subroutine calc_epv_dyn(configs, mesh, time_lev, state, diag)
      use mpas_log, only : mpas_log_write ! MW debug

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_config

      IMPLICIT NONE

      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev                     ! which time level to use from state: 1 = beginning of time step; 2 = end of time step
      type (mpas_pool_type), intent(inout) :: diag        ! MC: since these subroutines are called after state/diag vars have been updated, time_lev 
      type (mpas_pool_type), intent(in) :: mesh           !     specification here doesn't really matter
      type (mpas_pool_type), intent(in) :: configs 

      ! input vars
      integer, pointer :: nCells, nVertLevels, nEdges, R3
      integer, pointer :: nVertices, vertexDegree ! MW: adding this to recompute pv_vertex 
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex

      real(kind=RKIND), dimension(:), pointer :: dzu, areaCell, latCell, lonCell, dcEdge
      real(kind=RKIND), dimension(:), pointer :: invAreaTriangle ! MW: adding this to recompute pv_vertex
      integer, dimension(:,:), pointer :: edgesOnVertex ! MW: adding this to recompute pv_vertex
      real(kind=RKIND), dimension(:,:), pointer ::  edgesOnVertex_sign ! MW: adding this to recompute pv_vertex
      real(kind=RKIND), dimension(:),   pointer ::  fVertex ! MW: adding this to recompute pv_vertex
      real(kind=RKIND), dimension(:,:), pointer :: w, rho, theta, pv_vertex, zgrid, kiteAreasOnVertex, &
                                                   uReconstructZonal, uReconstructMeridional, ertel_pv, & 
                                                   edgesOnCell_sign, wCell
      real(kind=RKIND), dimension(:,:), pointer :: u ! MW: adding u to recompute pv_vertex
      real(kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct 

      ! local vars
      real(kind=RKIND), dimension(:,:), allocatable :: duZonal_dz, duMerid_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dTheta_dxZonal, dTheta_dyMerid, dTheta_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dW_dxZonal, dW_dyMerid 
      real(kind=RKIND), dimension(:,:), allocatable :: absVort
      real(kind=RKIND), dimension(:,:,:), allocatable :: absVort3D, gradTheta
      real(kind=RKIND), dimension(:,:), pointer :: pv_vertex_dyn ! MW: adding this to recompute pv_vertex

      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_adv, ertel_pv_dyn 
      real(kind=RKIND), pointer:: config_dt

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)  ! MW: adding this to recompute pv_vertex
      call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)  ! MW: adding this to recompute pv_vertex
      call mpas_pool_get_array(mesh, 'edgesOnVertex', edgesOnVertex)  ! MW: adding this to recompute pv_vertex
      call mpas_pool_get_array(mesh, 'edgesOnVertex_sign', edgesOnVertex_sign)  ! MW: adding this to recompute pv_vertex
      call mpas_pool_get_array(mesh, 'invAreaTriangle', invAreaTriangle)  ! MW: adding this to recompute pv_vertex
      call mpas_pool_get_array(mesh, 'fVertex', fVertex)  ! MW: adding this to recompute pv_vertex

      call mpas_pool_get_dimension(mesh, 'R3', R3)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'coeffs_reconstruct', coeffs_reconstruct)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(diag, 'w_dyn', w)     
      call mpas_pool_get_array(diag, 'theta_dyn', theta) 
      call mpas_pool_get_array(diag, 'u_dyn', u)  ! MW: adding u_dyn to recompute relative vorticyt, pv_vertex 
      !MW call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'rho_dyn', rho)  ! MW: when testing with no dynamics update - need to use rho_prev
      !call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)  ! MW: we will recompute this field using u_dyn
      call mpas_pool_get_array(diag, 'uReconstructZonal_dyn', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional_dyn', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'ertel_pv_dyn', ertel_pv_dyn)
      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)  ! this should be valid at the beginning of the time step 
      call mpas_pool_get_array(diag, 'wCell_dyn', wCell)    
      call mpas_pool_get_array(diag, 'pv_vertex_dyn', pv_vertex_dyn)
      call mpas_pool_get_array(diag, 'depv_dt_adv', depv_dt_adv)

      call mpas_pool_get_config(configs,'config_dt',config_dt)
      !call mpas_log_write('calc_epv_dyn: config_dt = $r', realArgs=(/config_dt /))

      ! Allocate memory to intermediate vars 
      allocate(absVort(nVertLevels,nCells)) 
      allocate(duZonal_dz(nVertLevels,nCells))
      allocate(duMerid_dz(nVertLevels,nCells))
      allocate(dTheta_dxZonal(nVertLevels,nCells))
      allocate(dTheta_dyMerid(nVertLevels,nCells))
      allocate(dTheta_dz(nVertLevels,nCells))
      allocate(dW_dxZonal(nVertLevels,nCells))
      allocate(dW_dyMerid(nVertLevels,nCells))
      allocate(absVort3D(nVertLevels,nCells,3))
      allocate(gradTheta(nVertLevels,nCells,3))

      ertel_pv_dyn(:,:) = 0.0_RKIND
      gradTheta(:,:,:) = 0.0_RKIND
      absVort3D(:,:,:) = 0.0_RKIND
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D potential temperature gradient 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !call mpas_log_write('calc_epv_dyn: theta(k=55,iCell=100) = $r',realArgs=(/ theta(55,100) /) )
      !call mpas_log_write('calc_epv_dyn: uReconstructZonal(k=55,iCell=100) = $r',realArgs=(/ uReconstructZonal(55,100) /) )
      !call mpas_log_write('calc_epv_dyn: uReconstructMeridional(k=55,iCell=100) = $r',realArgs=(/ uReconstructMeridional(55,100) /) )

      ! (1) Calculate and reconstruct horizontal potential temperature gradient to get zonal and meridional
      !     gradients at cell centers: dth_dx, dth_dy

      call calc_gradOnEdges_reconCellCenter(theta, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dTheta_dxZonal, dTheta_dyMerid)


      ! (2) Calculate the vertical potential temperature gradient: dth_dz

      call calc_vertDeriv(theta, nCells, nVertLevels, dzu, dTheta_dz)

      ! (3) Combine theta derivatives into 3D vector

      gradTheta(:,:,1) = dTheta_dxZonal
      gradTheta(:,:,2) = dTheta_dyMerid
      gradTheta(:,:,3) = dTheta_dz

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D absolute vorticity vector  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate the vertical shear of uReconstructZonal and uReconstructMeridional: du_dz and dv_dz

      call calc_vertDeriv(uReconstructZonal, nCells, nVertLevels, dzu, duZonal_dz)
      call calc_vertDeriv(uReconstructMeridional, nCells, nVertLevels, dzu, duMerid_dz)
      !call mpas_log_write('calc_epv_dyn: duZonal_dz(k=55,iCell=100) = $r',realArgs=(/ duZonal_dz(55,100) /) )
      !call mpas_log_write('calc_epv_dyn: duMerid_dz(k=55,iCell=100) = $r',realArgs=(/ duMerid_dz(55,100) /) )


      ! (2) Interpolate w to cell centers, calculate gradient of w on edges, and then reconstruct to get 
      !     zonal and meridional gradients at cell centers: dw_dx, dw_dy 

      call interp_wLev_thetaLev(w, nCells, nVertLevels, wCell)
      !call mpas_log_write('calc_epv_dyn: wCell(k=55,iCell=100) = $r',realArgs=(/ wCell(55,100) /) )

      call calc_gradOnEdges_reconCellCenter(wCell, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dW_dxZonal, dW_dyMerid)

      ! (3) Recconstruct absolute vertical vorticity at vertices pv_vertex to cell centers

      !     Note: currently, pv_vertex is the absolute vertical vorticity on the cell vertices. If this 
      !     variable changes at some point, then the absolute vertical vorticity on the vertices needs 
      !     to be computed as follows: 
      !     do iVert=1,nVertices
      !        vorticity(:,iVert) = vorticity(:,iVert) + fVertex(iVert)
      !     end do

      ! MW: recompute pv_vertex based on velocities after dynamics
      pv_vertex_dyn(:,:) = 0.0_RKIND
      call recompute_pv_vertex(u, nVertices, nVertLevels, vertexDegree, invAreaTriangle, &
                                 dcEdge, edgesOnVertex, edgesOnVertex_sign, fVertex, pv_vertex_dyn ) 
      !call mpas_log_write('calc_epv_dyn: pv_vertex(k=55,iCell=100) = $r',realArgs=(/ pv_vertex_dyn(55,100) /) )

      call interp_absVertVort(pv_vertex_dyn, nCells, nEdgesOnCell, verticesOnCell, &
                                  cellsOnVertex, areaCell, kiteAreasOnVertex, absVort)

      ! (4) Combine three components into vorticity vector

      absVort3D(:,:,1) = dW_dyMerid - duMerid_dz                   ! dw/dy - dv/dz
      absVort3D(:,:,2) = duZonal_dz - dW_dxZonal                   ! du/dz - dw/dx
      absVort3D(:,:,3) = absVort                                   ! dv/dy - du/dx + f

      ! (5) Take dot product between 3D theta vector and absolute vorticity / density to compute 
      !     Ertel's PV
      !call mpas_log_write('calc_epv_dyn: absVort3D components (k=55,iCell=100) = $r $r $r',realArgs=(/ absVort3D(55,100,1), absVort3D(55,100,2), absVort3D(55,100,3) /) )
      !call mpas_log_write('calc_epv_dyn: gradTheta components (k=55,iCell=100) = $r $r $r',realArgs=(/ gradTheta(55,100,1), gradTheta(55,100,2), gradTheta(55,100,3) /) )
      call calc_dotProduct_3D(gradTheta, absVort3D, nCells, nVertLevels, ertel_pv_dyn)
 
      !call mpas_log_write('calc_epv_dyn: ertel_pv_dyn(k=55,iCell=100) = $r',realArgs=(/ ertel_pv_dyn(55,100) /) )
      !call mpas_log_write('calc_epv_dyn: rho(k=55,iCell=100) = $r',realArgs=(/ rho(55,100) /) )
      ertel_pv_dyn = ertel_pv_dyn / rho * 1.0e6                            !SI to PVUs
      !call mpas_log_write('calc_epv_dyn: ertel_pv_dyn(k=55,iCell=100) = $r',realArgs=(/ ertel_pv_dyn(55,100) /) )

      depv_dt_adv(:,:) = 0.0_RKIND
      depv_dt_adv = ( ertel_pv_dyn - ertel_pv ) / config_dt


      deallocate(absVort)
      deallocate(absVort3D)
      deallocate(gradTheta)
      deallocate(duZonal_dz)
      deallocate(duMerid_dz)
      deallocate(dTheta_dz)
      deallocate(dTheta_dxZonal)
      deallocate(dTheta_dyMerid)
      deallocate(dW_dxZonal)
      deallocate(dW_dyMerid)

      call mpas_log_write('done with advection')

   end subroutine calc_epv_dyn

   !***********************************************************************************************
   ! MC: Subroutine to calculate Ertel's potential vorticity advection terms in zonal, meridional, 
   !     vertical directions  
   !***********************************************************************************************

   subroutine calc_pvAdvection(state, time_lev, diag, mesh)

      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array

      implicit none

      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev                        
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh

      integer, pointer :: nCells, nVertLevels, nEdges
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell

      real(kind=RKIND), dimension(:),   pointer :: dzu, latCell, lonCell, dcEdge

      real(kind=RKIND), dimension(:,:), pointer :: uReconstructZonal, uReconstructMeridional, edgesOnCell_sign, &
                                                   wCell, ertel_pv
      real(kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

      ! tendency variables - PV advection
      real(kind=RKIND), dimension(:,:), pointer :: dEPV_dxZonal, dEPV_dyMerid, dEPV_dz, &
                                                   advectionPV_zonal, advectionPV_meridional, advectionPV_vertical

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'coeffs_reconstruct', coeffs_reconstruct)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'wCell', wCell)

      call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)                
      call mpas_pool_get_array(diag, 'dEPV_dxZonal', dEPV_dxZonal)
      call mpas_pool_get_array(diag, 'dEPV_dyMerid', dEPV_dyMerid)
      call mpas_pool_get_array(diag, 'dEPV_dz', dEPV_dz)

      call mpas_pool_get_array(diag, 'advectionPV_zonal', advectionPV_zonal)   
      call mpas_pool_get_array(diag, 'advectionPV_meridional', advectionPV_meridional)
      call mpas_pool_get_array(diag, 'advectionPV_vertical', advectionPV_vertical)

      dEPV_dxZonal(:,:) = 0.0_RKIND
      dEPV_dyMerid(:,:) = 0.0_RKIND
      dEPV_dz(:,:) = 0.0_RKIND
      advectionPV_zonal(:,:) = 0.0_RKIND
      advectionPV_meridional(:,:) = 0.0_RKIND
      advectionPV_vertical(:,:) = 0.0_RKIND

      !***********************************************************************************************
      ! Calculate the PV advection terms
      ! Note: these are computed using the reconstructed edge gradients and zonal/meridional wind 
      ! components, which is much different and simpler than how advection is computed in the MPAS 
      ! dynamics code.
      !***********************************************************************************************

      ! (1) Calculate and reconstruct horizontal potential vorticity gradient to get zonal and meridional
      !     gradients at cell centers: dPV_dx, dPV_dy

      call calc_gradOnEdges_reconCellCenter(ertel_pv, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dEPV_dxZonal, dEPV_dyMerid)

      ! (2) Calculate the vertical potential vorticity gradient: dEPV_dz

      call calc_vertDeriv(ertel_pv, nCells, nVertLevels, dzu, dEPV_dz)

      ! (3) Calculate advection terms: u*dPV/dx, v*dPV/dy, w*dPV/dz

      advectionPV_zonal = uReconstructZonal * dEPV_dxZonal
      advectionPV_meridional = uReconstructMeridional * dEPV_dyMerid
      advectionPV_vertical = wCell * dEPV_dz

   end subroutine calc_pvAdvection

   !*********************************************************************************************************************
   ! MC: Modified subroutine to call PV and PV advection calculations and interpolation of diagnostic fields to dynamic 
   !     tropopause
   !********************************************************************************************************************* 

   subroutine atm_compute_pv_diagnostics(configs, state, time_lev, diag, mesh) 
   
      use mpas_constants
      use mpas_derived_types, only : field2DReal
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
      use mpas_dmpar, only : mpas_dmpar_exch_halo_field
      use mpas_log, only : mpas_log_write
 
      implicit none
   
      type (mpas_pool_type), intent(inout) :: state
      integer, intent(in) :: time_lev                                           
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: configs 
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal, stratoPV

      type (field2DReal), pointer :: theta_f, uReconstructZonal_f, uReconstructMeridional_f, w_f, epv_f, pv_vertex_f, &
                                     wCell_f !, vorticity_f
      type (field2DReal), pointer :: u_f ! MW: added for recomputing pv_vertex 
      type (field2DReal), pointer :: rho_f ! MW: added for rho_dyn 
      type (field2DReal), pointer :: depv_dt_adv_f

      !!!!!!!!!!!!!!!!!!!!
      ! MW: call calc_epv_dyn subroutine
      !     - call calc_epv_dyn before calc_epv to get the change in epv due to advection
      !           = ertel_pv_dyn - ertel_pv at the beginning of the time step
      call mpas_pool_get_field(diag, 'u_dyn', u_f)
      call mpas_pool_get_field(diag, 'w_dyn', w_f)
      call mpas_pool_get_field(diag, 'wCell_dyn', wCell_f)
      call mpas_pool_get_field(diag, 'uReconstructZonal_dyn', uReconstructZonal_f)
      call mpas_pool_get_field(diag, 'uReconstructMeridional_dyn', uReconstructMeridional_f)
      call mpas_pool_get_field(diag, 'theta_dyn', theta_f)
      call mpas_pool_get_field(diag, 'rho_dyn', rho_f)
      call mpas_pool_get_field(diag, 'ertel_pv', epv_f )  ! this should be values at the end of last time step


      call mpas_dmpar_exch_halo_field(u_f)
      call mpas_dmpar_exch_halo_field(w_f)
      call mpas_dmpar_exch_halo_field(wCell_f)
      call mpas_dmpar_exch_halo_field(uReconstructZonal_f)
      call mpas_dmpar_exch_halo_field(uReconstructMeridional_f)
      call mpas_dmpar_exch_halo_field(theta_f)
      call mpas_dmpar_exch_halo_field(rho_f)
      call mpas_dmpar_exch_halo_field(epv_f)

      call calc_epv_dyn(configs, mesh, time_lev, state, diag)

      call mpas_pool_get_field(diag, 'ertel_pv_dyn', epv_f)
      call mpas_dmpar_exch_halo_field(epv_f)
      call mpas_pool_get_field(diag, 'depv_dt_adv', depv_dt_adv_f)
      call mpas_dmpar_exch_halo_field(depv_dt_adv_f)

      ! 




      call mpas_pool_get_field(state, 'w', w_f, time_lev)                                          
      call mpas_pool_get_field(diag, 'uReconstructZonal', uReconstructZonal_f)                     ! MC added
      call mpas_pool_get_field(diag, 'uReconstructMeridional', uReconstructMeridional_f)           ! MC added
      call mpas_pool_get_field(diag, 'wCell', wCell_f)                                             ! MC added 
      call mpas_pool_get_field(diag, 'theta', theta_f)
      call mpas_pool_get_field(diag, 'pv_vertex', pv_vertex_f)                                     ! MC added
      call mpas_pool_get_field(diag, 'rho', rho_f)    !!! NOTE: VAR RHO_F CALLED TWICE                                 ! MW testing 
      !call mpas_pool_get_field(diag, 'vorticity', vorticity_f) 

      !need halo cells for everything w/ horizontal derivative 
      call mpas_dmpar_exch_halo_field(w_f)
      call mpas_dmpar_exch_halo_field(uReconstructZonal_f)                               ! MC added
      call mpas_dmpar_exch_halo_field(uReconstructMeridional_f)                          ! MC added
      call mpas_dmpar_exch_halo_field(wCell_f)                                           ! MC added
      call mpas_dmpar_exch_halo_field(theta_f)
      call mpas_dmpar_exch_halo_field(pv_vertex_f)                                       ! MC added
      call mpas_dmpar_exch_halo_field(rho_f)                                       ! MW testing
      !call mpas_dmpar_exch_halo_field(vorticity_f)
 
      ! call calc_epv subroutine 
      call calc_epv(mesh, time_lev, state, diag)                  

      ! halo cells need to be valid for flood fill routines called below
      call mpas_pool_get_field(diag, 'ertel_pv', epv_f)                       
      call mpas_dmpar_exch_halo_field(epv_f)

    
      ! call calc_pvAdvection subroutine
      call calc_pvAdvection(state, time_lev, diag, mesh)

      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND
      stratoPV = 10.0_RKIND

      !***********************************************************************************************
      ! Uncomment one of the following to either call floodFill_strato or floodFill_tropo
      ! to estimate the level of the dynamic tropopause.
      !***********************************************************************************************
     
      !call floodFill_strato(mesh, diag, pvuVal, stratoPV)    
      call floodFill_tropo(mesh,diag,pvuVal)

      ! call interp_pv_diagnostics subroutine to interpolate fields to dynamic tropopause
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
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      use mpas_log, only : mpas_log_write
 
      implicit none
      
      type (mpas_pool_type), intent(in) :: configs
      type (mpas_pool_type), intent(in) :: state
      integer, intent(in) :: time_lev                                            ! which time level to use from state
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh
      type (mpas_pool_type), intent(in) :: tend_physics
      type (mpas_pool_type), intent(in) :: diag_physics
      type (mpas_pool_type), intent(inout) :: tend                               ! modify tend_w_euler to uncouple with density

      ! input vars
      integer, pointer :: nCells, nVertLevels, nEdges, R3, nVertices
      integer, dimension(:), pointer :: nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, verticesOnCell, &
                                          cellsOnVertex, verticesOnEdge   

      real(kind=RKIND), pointer:: config_dt
      real(kind=RKIND), dimension(:),   pointer :: dzu, areaCell, latCell, lonCell, dcEdge, areaTriangle

      ! call variables from end of time step
      real(kind=RKIND), dimension(:,:), pointer :: rho, pv_vertex, uReconstructZonal, uReconstructMeridional, wCell

      ! need to call stored variables from beginning of previous time step 
      real(kind=RKIND), dimension(:,:), pointer :: rho_prev, theta_prev, zgrid, kiteAreasOnVertex, &
                                                   edgesOnCell_sign, &
                                                   ertel_pv_prev
                                                   ! pv_vertex_prev, rho_edge_prev, uReconstructZonal_prev, uReconstructMeridional_prev, wCell_prev 
      !real(kind=RKIND), dimension(:,:), pointer :: rho_zz_prev ! MW: added for decoupling friction term

      real(kind=RKIND), dimension(:,:,:), pointer :: coeffs_reconstruct

      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_diab, depv_dt_fric
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_rho

      real(kind=RKIND), dimension(:,:), pointer :: tend_u_phys, tend_u_euler, rho_edge, tend_w_euler
      real(kind=RKIND), dimension(:,:), pointer :: rthblten, rthcuten, rthratenlw, rthratensw, &
                                                   dtheta_dt_mp, dtheta_dt_mix

      real(kind=RKIND), dimension(:,:), pointer :: tend_theta_mp_evap_cw, tend_theta_mp_evap_rw, tend_theta_mp_depo_ice, &  ! Potential temperature tendency contributions from specific microphysical processes
                                                   tend_theta_mp_melt_ice, tend_theta_mp_frez_ice

      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_evap_cw, depv_dt_mp_evap_rw, depv_dt_mp_depo_ice, &      ! PV tendency contrubutions from specific microphysical proceses 
                                                   depv_dt_mp_melt_ice, depv_dt_mp_frez_ice, depv_dt_mp_allproc

      ! local static variables 
      real(kind=RKIND), dimension(:,:), allocatable :: duZonal_dz, duMerid_dz
      real(kind=RKIND), dimension(:,:), allocatable :: dTheta_dxZonal, dTheta_dyMerid, dTheta_dz 
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

      real(kind=RKIND), dimension(:,:,:), allocatable :: grad_diabatic_LW, grad_diabatic_SW, grad_diabatic_BL, &
                                                         grad_diabatic_CU, grad_diabatic_MP, grad_diabatic_MX

      real(kind=RKIND), dimension(:,:,:), allocatable :: grad_diabatic_MP_evap_cw, grad_diabatic_MP_evap_rw, grad_diabatic_MP_depo, &
                                                         grad_diabatic_MP_melt, grad_diabatic_MP_frez, grad_diabatic_MP_sum

      real(kind=RKIND), dimension(:,:), allocatable :: tend_theta_mp_sum

      ! friction tendency vars 
      real(kind=RKIND), dimension(:,:), allocatable :: dWtend_dxZonal, dWtend_dyMerid, duZonalTend_dz, duMeridTend_dz, &
                                                       vertVortTend, tenduX, tenduY, tenduZ, tend_uZonal, tend_uMerid
      real(kind=RKIND), dimension(:,:,:), allocatable :: vortTend3D
      real(kind=RKIND), dimension(:,:), pointer :: tend_wCell, uTend_curl

      ! density tendency vars
      real(kind=RKIND), dimension(:,:), allocatable :: rho_tend

      ! MW : testing advection tendency
      real(kind=RKIND), dimension(:,:), pointer :: w_dyn, theta_dyn
      real(kind=RKIND), dimension(:,:), pointer :: uReconstructZonal_dyn, uReconstructMeridional_dyn

      call mpas_log_write('check 3')

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'R3', R3)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(mesh, 'nEdges', nEdges)
      call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
      call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(mesh, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(mesh, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
      call mpas_pool_get_array(mesh, 'verticesOnEdge', verticesOnEdge)
      call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
      call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
      call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
      call mpas_pool_get_array(mesh, 'edgesOnCell_sign', edgesOnCell_sign)
      call mpas_pool_get_array(mesh, 'latCell', latCell)
      call mpas_pool_get_array(mesh, 'lonCell', lonCell)
      call mpas_pool_get_array(mesh, 'dcEdge', dcEdge)
      call mpas_pool_get_array(mesh, 'areaCell', areaCell)
      call mpas_pool_get_array(mesh, 'coeffs_reconstruct', coeffs_reconstruct)
      call mpas_pool_get_array(mesh, 'zgrid', zgrid)
      call mpas_pool_get_array(mesh, 'dzu', dzu)
      call mpas_pool_get_config(configs,'config_dt',config_dt)

      call mpas_log_write('check 4')
      call mpas_pool_get_array(diag, 'rho', rho)
      call mpas_pool_get_array(diag, 'theta_prev', theta_prev)    
      call mpas_pool_get_array(diag, 'rho_prev', rho_prev)
      call mpas_pool_get_array(diag, 'pv_vertex', pv_vertex)                      ! CHANGING: modified time level
      call mpas_pool_get_array(diag, 'uReconstructZonal', uReconstructZonal)
      call mpas_pool_get_array(diag, 'uReconstructMeridional', uReconstructMeridional)
      call mpas_pool_get_array(diag, 'wCell', wCell)
      call mpas_pool_get_array(diag, 'ertel_pv_prev', ertel_pv_prev)

!      call mpas_pool_get_array(diag, 'pv_vertex_prev', pv_vertex_prev)
!      call mpas_pool_get_array(diag, 'uReconstructZonal_prev', uReconstructZonal_prev)
!      call mpas_pool_get_array(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev)
!      call mpas_pool_get_array(diag, 'wCell_prev', wCell_prev)
!      call mpas_pool_get_array(diag, 'rho_zz_prev', rho_zz_prev)  ! MW: added for decoupling friction terms
      
      call mpas_pool_get_array(tend_physics, 'rthblten', rthblten)
      call mpas_pool_get_array(tend_physics, 'rthcuten', rthcuten)
      call mpas_pool_get_array(tend_physics, 'rthratenlw', rthratenlw)
      call mpas_pool_get_array(tend_physics, 'rthratensw', rthratensw)
      call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)                      ! Potential temperature tendency from explicit horizontal mixing
      
      call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)
      call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)
      call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)
      call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)
      call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)
      call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      call mpas_pool_get_array(diag, 'depv_dt_rho', depv_dt_rho)

      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_evap_cw', tend_theta_mp_evap_cw)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_evap_rw', tend_theta_mp_evap_rw)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_depo_ice', tend_theta_mp_depo_ice)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_melt_ice', tend_theta_mp_melt_ice)
      call mpas_pool_get_array(diag_physics, 'tend_theta_mp_frez_ice', tend_theta_mp_frez_ice)

      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_cw', depv_dt_mp_evap_cw)
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_rw', depv_dt_mp_evap_rw)
      call mpas_pool_get_array(diag, 'depv_dt_mp_depo_ice', depv_dt_mp_depo_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_melt_ice', depv_dt_mp_melt_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_frez_ice', depv_dt_mp_frez_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_allproc', depv_dt_mp_allproc)

      call mpas_pool_get_array(diag, 'tend_u_phys', tend_u_phys)                          ! Normal wind tendencies from physics (defined in Registry_pv.xml)
      call mpas_pool_get_array(diag, 'tend_wCell', tend_wCell)
      !call mpas_pool_get_array(diag, 'rho_edge_prev', rho_edge_prev)
      !MW call mpas_pool_get_array(tend, 'u_euler', tend_u_euler)                             ! Tendency of u from explicit horizontal mixing 
      !MW call mpas_pool_get_array(tend, 'w_euler', tend_w_euler)                             ! Tendency of w from explicit horizontal mixing 
      call mpas_pool_get_array(diag,'u_tend_diff',tend_u_euler)   ! MW: using ru explicit horiz mixing from tendencies
      call mpas_pool_get_array(diag,'w_tend_diff',tend_w_euler)   ! MW: using rw explicit horiz mixing from tendencies 
      call mpas_pool_get_array(diag, 'uTend_curl', uTend_curl)

      ! MW: test - advection tendency term
      call mpas_pool_get_array(diag, 'w_dyn', w_dyn )         ! horizontal winds after dynamics 
      call mpas_pool_get_array(diag, 'theta_dyn', theta_dyn ) ! potential temperature after advection
      call mpas_pool_get_array(diag,'uReconstructZonal_dyn', uReconstructZonal_dyn)
      call mpas_pool_get_array(diag,'uReconstructMeridional_dyn', uReconstructMeridional_dyn)

      call mpas_log_write('check 5')

      ! allocate local static variables
      allocate(duZonal_dz(nVertLevels,nCells))
      allocate(duMerid_dz(nVertLevels,nCells))
      allocate(dTheta_dxZonal(nVertLevels,nCells))
      allocate(dTheta_dyMerid(nVertLevels,nCells))
      allocate(dTheta_dz(nVertLevels,nCells))
      allocate(dW_dxZonal(nVertLevels,nCells))
      allocate(dW_dyMerid(nVertLevels,nCells))
      allocate(absVort(nVertLevels,nCells))
      allocate(absVort3D(nVertLevels,nCells,3))
      allocate(gradTheta(nVertLevels,nCells,3))
     
      ! allocate diabatic tendency variables
      allocate(dLWtend_dxZonal(nVertLevels,nCells))
      allocate(dLWtend_dyMerid(nVertLevels,nCells))
      allocate(dLWtend_dz(nVertLevels,nCells))
      allocate(dSWtend_dxZonal(nVertLevels,nCells))
      allocate(dSWtend_dyMerid(nVertLevels,nCells))
      allocate(dSWtend_dz(nVertLevels,nCells))
      allocate(dBLtend_dxZonal(nVertLevels,nCells))
      allocate(dBLtend_dyMerid(nVertLevels,nCells))
      allocate(dBLtend_dz(nVertLevels,nCells))
      allocate(dCUtend_dxZonal(nVertLevels,nCells))
      allocate(dCUtend_dyMerid(nVertLevels,nCells))
      allocate(dCUtend_dz(nVertLevels,nCells))
      allocate(dMPtend_dxZonal(nVertLevels,nCells))
      allocate(dMPtend_dyMerid(nVertLevels,nCells))
      allocate(dMPtend_dz(nVertLevels,nCells))
      allocate(dMXtend_dxZonal(nVertLevels,nCells))
      allocate(dMXtend_dyMerid(nVertLevels,nCells))
      allocate(dMXtend_dz(nVertLevels,nCells))
      allocate(grad_diabatic_LW(nVertLevels,nCells,R3))
      allocate(grad_diabatic_SW(nVertLevels,nCells,R3))
      allocate(grad_diabatic_BL(nVertLevels,nCells,R3))
      allocate(grad_diabatic_CU(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MX(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MP(nVertLevels,nCells,R3))

      ! allocate diabatic tendency variables from specific microphys processes
      allocate(dMPevapcwtend_dxZonal(nVertLevels,nCells))
      allocate(dMPevapcwtend_dyMerid(nVertLevels,nCells))
      allocate(dMPevapcwtend_dz(nVertLevels,nCells))
      allocate(dMPevaprwtend_dxZonal(nVertLevels,nCells))
      allocate(dMPevaprwtend_dyMerid(nVertLevels,nCells))
      allocate(dMPevaprwtend_dz(nVertLevels,nCells))
      allocate(dMPdepotend_dxZonal(nVertLevels,nCells))
      allocate(dMPdepotend_dyMerid(nVertLevels,nCells))
      allocate(dMPdepotend_dz(nVertLevels,nCells))
      allocate(dMPmelttend_dxZonal(nVertLevels,nCells))
      allocate(dMPmelttend_dyMerid(nVertLevels,nCells))
      allocate(dMPmelttend_dz(nVertLevels,nCells))
      allocate(dMPfreztend_dxZonal(nVertLevels,nCells))
      allocate(dMPfreztend_dyMerid(nVertLevels,nCells))
      allocate(dMPfreztend_dz(nVertLevels,nCells))
      allocate(dMPsumtend_dxZonal(nVertLevels,nCells))
      allocate(dMPsumtend_dyMerid(nVertLevels,nCells))
      allocate(dMPsumtend_dz(nVertLevels,nCells))
      allocate(tend_theta_mp_sum(nVertLevels,nCells))

      allocate(grad_diabatic_MP_evap_cw(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MP_evap_rw(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MP_depo(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MP_melt(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MP_frez(nVertLevels,nCells,R3))
      allocate(grad_diabatic_MP_sum(nVertLevels,nCells,R3))

      ! allocate friction tendency variables
      allocate(dWtend_dxZonal(nVertLevels,nCells))
      allocate(dWtend_dyMerid(nVertLevels,nCells))
      allocate(duZonalTend_dz(nVertLevels,nCells))
      allocate(duMeridTend_dz(nVertLevels,nCells))
      allocate(vertVortTend(nVertLevels,nCells))
      allocate(tenduX(nVertLevels,nCells+1))
      allocate(tenduY(nVertLevels,nCells+1))
      allocate(tenduZ(nVertLevels,nCells+1))
      allocate(tend_uZonal(nVertLevels,nCells+1))
      allocate(tend_uMerid(nVertLevels,nCells+1))
      allocate(vortTend3D(nVertLevels,nCells,3))

      ! allocate rho tendency variables
      allocate(rho_tend(nVertLevels,nCells))

      depv_dt_lw(:,:) = 0.0_RKIND
      depv_dt_sw(:,:) = 0.0_RKIND
      depv_dt_bl(:,:) = 0.0_RKIND
      depv_dt_cu(:,:) = 0.0_RKIND
      depv_dt_mp(:,:) = 0.0_RKIND 
      depv_dt_mix(:,:) = 0.0_RKIND
      depv_dt_diab(:,:) = 0.0_RKIND
      depv_dt_fric(:,:) = 0.0_RKIND
      depv_dt_rho(:,:) = 0.0_RKIND

      dWtend_dxZonal(:,:) =  0.0_RKIND
      dWtend_dyMerid(:,:) =  0.0_RKIND
      depv_dt_mp_evap_cw(:,:) =  0.0_RKIND
      depv_dt_mp_evap_rw(:,:) =  0.0_RKIND
      depv_dt_mp_depo_ice(:,:) =  0.0_RKIND
      depv_dt_mp_melt_ice(:,:) =  0.0_RKIND
      depv_dt_mp_frez_ice(:,:) =  0.0_RKIND
      tend_theta_mp_sum(:,:) =  0.0_RKIND
      depv_dt_mp_allproc(:,:) =  0.0_RKIND

      call mpas_log_write('Beginning calculating terms for PV tendencies')

      !***********************************************************************************************
      ! Calculate terms needed for PV tendency equation
      !***********************************************************************************************

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D potential temperature gradient using theta at beginning of time step  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate and reconstruct horizontal potential temperature gradient to get zonal and meridional
      !     gradients at cell centers: dth_dx, dth_dy

      call calc_gradOnEdges_reconCellCenter(theta_prev, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dTheta_dxZonal, dTheta_dyMerid)

      ! (2) Calculate the vertical potential temperature gradient: dth_dz

      call calc_vertDeriv(theta_prev, nCells, nVertLevels, dzu, dTheta_dz)

      ! (3) Combine theta derivatives into 3D vector

      gradTheta(:,:,1) = dTheta_dxZonal
      gradTheta(:,:,2) = dTheta_dyMerid
      gradTheta(:,:,3) = dTheta_dz


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! MC: TIME LEVEL MOD. Calculate the 3D absolute vorticity vector using winds at end of time step  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate the vertical shear of uReconstructZonal and uReconstructMeridional: du_dz and dv_dz

      call calc_vertDeriv(uReconstructZonal, nCells, nVertLevels, dzu, duZonal_dz)
      call calc_vertDeriv(uReconstructMeridional, nCells, nVertLevels, dzu, duMerid_dz)

      ! (2) Interpolate w to cell centers, calculate gradient of w on edges, and then reconstruct to get 
      !     zonal and meridional gradients at cell centers: dw_dx, dw_dy 

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

      absVort3D(:,:,1) = dW_dyMerid - duMerid_dz        ! dw/dy - dv/dz
      absVort3D(:,:,2) = duZonal_dz - dW_dxZonal        ! du/dz - dw/dx
      absVort3D(:,:,3) = absVort                        ! dv/dy - du/dx + f



      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Calculate the 3D absolute vorticity vector using winds at beginning of time step  
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      ! (1) Calculate the vertical shear of uReconstructZonal and uReconstructMeridional: du_dz and dv_dz

!      call calc_vertDeriv(uReconstructZonal_prev, nCells, nVertLevels, dzu, duZonal_dz)
!      call calc_vertDeriv(uReconstructMeridional_prev, nCells, nVertLevels, dzu, duMerid_dz)

      ! (2) Interpolate w to cell centers, calculate gradient of w on edges, and then reconstruct to get 
      !     zonal and meridional gradients at cell centers: dw_dx, dw_dy 

!      call calc_gradOnEdges_reconCellCenter(wCell_prev, nCells, nEdges, nVertLevels, nEdgesOnCell, &
!                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
!                                latCell, lonCell, dW_dxZonal, dW_dyMerid)

      ! (3) Reconstruct absolute vertical vorticity at vertices pv_vertex to cell centers
   
      !     Note: currently, pv_vertex is the absolute vertical vorticity on the cell vertices. If this 
      !     variable changes at some point, then the absolute vertical vorticity on the vertices needs 
      !     to be computed as follows: 
      !     do iVert=1,nVertices
      !        vorticity(:,iVert) = vorticity(:,iVert) + fVertex(iVert)
      !     end do

!      call interp_absVertVort(pv_vertex_prev, nCells, nEdgesOnCell, verticesOnCell, &
!                                cellsOnVertex, areaCell, kiteAreasOnVertex, absVort)

      ! (4) Combine three components into vorticity vector

!      absVort3D(:,:,1) = dW_dyMerid - duMerid_dz        ! dw/dy - dv/dz
!      absVort3D(:,:,2) = duZonal_dz - dW_dxZonal        ! du/dz - dw/dx
!      absVort3D(:,:,3) = absVort                        ! dv/dy - du/dx + f

      !***********************************************************************************************
      ! Calculate diabatic PV tendency terms:
      !***********************************************************************************************

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Longwave radiation tendency: depv_dt_lw
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(rthratenlw)) then
         ! Calculate and reconstruct horizontal gradients of theta tendency

         call calc_gradOnEdges_reconCellCenter(rthratenlw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                  edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                  latCell, lonCell, dLWtend_dxZonal, & 
                                  dLWtend_dyMerid)

         ! Calculate vertical gradient of theta tendency
         call calc_vertDeriv(rthratenlw, nCells, nVertLevels, dzu, dLWtend_dz)
          
         ! Combine into 3D theta tendency gradient vector 
         grad_diabatic_LW(:,:,1) = dLWtend_dxZonal
         grad_diabatic_LW(:,:,2) = dLWtend_dyMerid
         grad_diabatic_LW(:,:,3) = dLWtend_dz
         
         call calc_dotProduct_3D(grad_diabatic_LW, absVort3D, nCells, nVertLevels, depv_dt_lw)
         
         !depv_dt_lw = depv_dt_lw / rho_prev * 1.0e6                             ! MC: using density from beginning of time step here... 
         depv_dt_lw = depv_dt_lw / rho * 1.0e6 
      else
         depv_dt_lw = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Shortwave radiation tendency: depv_dt_sw
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(rthratensw)) then
         ! Calculate and reconstruct horizontal gradients of theta tendency
 
          call calc_gradOnEdges_reconCellCenter(rthratensw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dSWtend_dxZonal, &
                                   dSWtend_dyMerid)
 
          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(rthratensw, nCells, nVertLevels, dzu, dSWtend_dz)
 
          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_SW(:,:,1) = dSWtend_dxZonal
          grad_diabatic_SW(:,:,2) = dSWtend_dyMerid
          grad_diabatic_SW(:,:,3) = dSWtend_dz
 
          call calc_dotProduct_3D(grad_diabatic_SW, absVort3D, nCells, nVertLevels, depv_dt_sw)
 
          !depv_dt_sw = depv_dt_sw / rho_prev * 1.0e6
          depv_dt_sw = depv_dt_sw / rho * 1.0e6
      else
          depv_dt_sw = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! PBL diabatic tendency: depv_dt_bl
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(rthblten)) then
         ! Calculate and reconstruct horizontal gradients of theta tendency
 
          call calc_gradOnEdges_reconCellCenter(rthblten, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dBLtend_dxZonal, &
                                   dBLtend_dyMerid)
 
          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(rthblten, nCells, nVertLevels, dzu, dBLtend_dz)
 
          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_BL(:,:,1) = dBLtend_dxZonal
          grad_diabatic_BL(:,:,2) = dBLtend_dyMerid
          grad_diabatic_BL(:,:,3) = dBLtend_dz
 
          call calc_dotProduct_3D(grad_diabatic_BL, absVort3D, nCells, nVertLevels, depv_dt_bl)
 
          !depv_dt_bl = depv_dt_bl / rho_prev * 1.0e6
          depv_dt_bl = depv_dt_bl / rho * 1.0e6
      else
          depv_dt_bl = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Cumulus diabatic tendency: depv_dt_cu
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(rthcuten)) then
         ! Calculate and reconstruct horizontal gradients of theta tendency
 
          call calc_gradOnEdges_reconCellCenter(rthcuten, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dCUtend_dxZonal, &
                                   dCUtend_dyMerid)
  
          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(rthcuten, nCells, nVertLevels, dzu, dCUtend_dz)
 
          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_CU(:,:,1) = dCUtend_dxZonal
          grad_diabatic_CU(:,:,2) = dCUtend_dyMerid
          grad_diabatic_CU(:,:,3) = dCUtend_dz
 
          call calc_dotProduct_3D(grad_diabatic_CU, absVort3D, nCells, nVertLevels, depv_dt_cu)

          !depv_dt_cu = depv_dt_cu / rho_prev * 1.0e6 
          depv_dt_cu = depv_dt_cu / rho * 1.0e6
      else
          depv_dt_cu = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Microphysics diabatic tendency: depv_dt_mp
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_mp)) then
         ! Calculate and reconstruct horizontal gradients of theta tendency
 
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_mp, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPtend_dxZonal, &
                                   dMPtend_dyMerid)
 
          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_mp, nCells, nVertLevels, dzu, dMPtend_dz)
 
          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP(:,:,1) = dMPtend_dxZonal
          grad_diabatic_MP(:,:,2) = dMPtend_dyMerid
          grad_diabatic_MP(:,:,3) = dMPtend_dz
 
          call calc_dotProduct_3D(grad_diabatic_MP, absVort3D, nCells, nVertLevels, depv_dt_mp)
 
          !depv_dt_mp = depv_dt_mp / rho_prev * 1.0e6
          depv_dt_mp = depv_dt_mp / rho * 1.0e6
      else
          depv_dt_mp = 0.0_RKIND
      end if
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Diabatic tendency from explicit mixing: depv_dt_mix
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (associated(dtheta_dt_mix)) then
         ! Calculate and reconstruct horizontal gradients of theta tendency
  
          call calc_gradOnEdges_reconCellCenter(dtheta_dt_mix, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMXtend_dxZonal, &
                                   dMXtend_dyMerid)
 
          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(dtheta_dt_mix, nCells, nVertLevels, dzu, dMXtend_dz)
 
          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MX(:,:,1) = dMXtend_dxZonal
          grad_diabatic_MX(:,:,2) = dMXtend_dyMerid
          grad_diabatic_MX(:,:,3) = dMXtend_dz
 
          call calc_dotProduct_3D(grad_diabatic_MX, absVort3D, nCells, nVertLevels, depv_dt_mix)
 
          !depv_dt_mix = depv_dt_mix / rho_prev * 1.0e6
          depv_dt_mix = depv_dt_mix / rho * 1.0e6
      else
          depv_dt_mix = 0.0_RKIND
      end if
 
      ! Sum of all diabatic contributions to PV through potential temperature tendencies   
      depv_dt_diab = depv_dt_lw + depv_dt_sw + depv_dt_bl + depv_dt_cu + depv_dt_mp + depv_dt_mix

 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Individual diabatic tendencies from specific microphysical processes (not included in budget)
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Net cloud water condensation and evaporatiom
      if (associated(tend_theta_mp_evap_cw)) then

          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_evap_cw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPevapcwtend_dxZonal, &
                                   dMPevapcwtend_dyMerid)

          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_evap_cw, nCells, nVertLevels, dzu, dMPevapcwtend_dz)

          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_evap_cw(:,:,1) = dMPevapcwtend_dxZonal
          grad_diabatic_MP_evap_cw(:,:,2) = dMPevapcwtend_dyMerid
          grad_diabatic_MP_evap_cw(:,:,3) = dMPevapcwtend_dz

          call calc_dotProduct_3D(grad_diabatic_MP_evap_cw, absVort3D, nCells, nVertLevels, depv_dt_mp_evap_cw)

          !depv_dt_mp_evap_cw = depv_dt_mp_evap_cw / rho_prev * 1.0e6
          depv_dt_mp_evap_cw = depv_dt_mp_evap_cw / rho * 1.0e6
      else
          depv_dt_mp_evap_cw = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Rain water evaporation
      if (associated(tend_theta_mp_evap_rw)) then

          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_evap_rw, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPevaprwtend_dxZonal, &
                                   dMPevaprwtend_dyMerid)

          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_evap_rw, nCells, nVertLevels, dzu, dMPevaprwtend_dz)

          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_evap_rw(:,:,1) = dMPevaprwtend_dxZonal
          grad_diabatic_MP_evap_rw(:,:,2) = dMPevaprwtend_dyMerid
          grad_diabatic_MP_evap_rw(:,:,3) = dMPevaprwtend_dz

          call calc_dotProduct_3D(grad_diabatic_MP_evap_rw, absVort3D, nCells, nVertLevels, depv_dt_mp_evap_rw)

          !depv_dt_mp_evap_rw = depv_dt_mp_evap_rw / rho_prev * 1.0e6
          depv_dt_mp_evap_rw = depv_dt_mp_evap_rw / rho * 1.0e6
      else
          depv_dt_mp_evap_rw = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Net sublimation/deposition
      if (associated(tend_theta_mp_depo_ice)) then

          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_depo_ice, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPdepotend_dxZonal, &
                                   dMPdepotend_dyMerid)

          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_depo_ice, nCells, nVertLevels, dzu, dMPdepotend_dz)

          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_depo(:,:,1) = dMPdepotend_dxZonal
          grad_diabatic_MP_depo(:,:,2) = dMPdepotend_dyMerid
          grad_diabatic_MP_depo(:,:,3) = dMPdepotend_dz

          call calc_dotProduct_3D(grad_diabatic_MP_depo, absVort3D, nCells, nVertLevels, depv_dt_mp_depo_ice)

          !depv_dt_mp_depo_ice = depv_dt_mp_depo_ice / rho_prev * 1.0e6
          depv_dt_mp_depo_ice = depv_dt_mp_depo_ice / rho * 1.0e6

      else
          depv_dt_mp_depo_ice = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Melting
      if (associated(tend_theta_mp_melt_ice)) then

          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_melt_ice, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPmelttend_dxZonal, &
                                   dMPmelttend_dyMerid)

          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_melt_ice, nCells, nVertLevels, dzu, dMPmelttend_dz)

          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_melt(:,:,1) = dMPmelttend_dxZonal
          grad_diabatic_MP_melt(:,:,2) = dMPmelttend_dyMerid
          grad_diabatic_MP_melt(:,:,3) = dMPmelttend_dz

          call calc_dotProduct_3D(grad_diabatic_MP_melt, absVort3D, nCells, nVertLevels, depv_dt_mp_melt_ice)

          !depv_dt_mp_melt_ice = depv_dt_mp_melt_ice / rho_prev * 1.0e6
          depv_dt_mp_melt_ice = depv_dt_mp_melt_ice / rho * 1.0e6
      else
          depv_dt_mp_melt_ice = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Freezing
      if (associated(tend_theta_mp_frez_ice)) then

          call calc_gradOnEdges_reconCellCenter(tend_theta_mp_frez_ice, nCells, nEdges, nVertLevels, nEdgesOnCell, &   
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPfreztend_dxZonal, &
                                   dMPfreztend_dyMerid)

          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_frez_ice, nCells, nVertLevels, dzu, dMPfreztend_dz)

          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_frez(:,:,1) = dMPfreztend_dxZonal
          grad_diabatic_MP_frez(:,:,2) = dMPfreztend_dyMerid
          grad_diabatic_MP_frez(:,:,3) = dMPfreztend_dz

          call calc_dotProduct_3D(grad_diabatic_MP_frez, absVort3D, nCells, nVertLevels, depv_dt_mp_frez_ice)

          !depv_dt_mp_frez_ice = depv_dt_mp_frez_ice / rho_prev * 1.0e6
          depv_dt_mp_frez_ice = depv_dt_mp_frez_ice / rho * 1.0e6
      else
          depv_dt_mp_frez_ice = 0.0_RKIND
      end if

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
      ! Sum all processes together to find combined PV tendency from microphysics. Compare to depv_dt_mp
      if (associated(tend_theta_mp_frez_ice)) then
      
         tend_theta_mp_sum = tend_theta_mp_frez_ice + tend_theta_mp_melt_ice + tend_theta_mp_depo_ice + &
                             tend_theta_mp_evap_rw + tend_theta_mp_evap_cw

         call calc_gradOnEdges_reconCellCenter(tend_theta_mp_sum, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                   edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                   latCell, lonCell, dMPsumtend_dxZonal, &
                                   dMPsumtend_dyMerid)

          ! Calculate vertical gradient of theta tendency
          call calc_vertDeriv(tend_theta_mp_sum, nCells, nVertLevels, dzu, dMPsumtend_dz)

          ! Combine into 3D theta tendency gradient vector 
          grad_diabatic_MP_sum(:,:,1) = dMPsumtend_dxZonal
          grad_diabatic_MP_sum(:,:,2) = dMPsumtend_dyMerid
          grad_diabatic_MP_sum(:,:,3) = dMPsumtend_dz

          call calc_dotProduct_3D(grad_diabatic_MP_sum, absVort3D, nCells, nVertLevels, depv_dt_mp_allproc)

          !depv_dt_mp_allproc = depv_dt_mp_allproc / rho_prev * 1.0e6
          depv_dt_mp_allproc = depv_dt_mp_allproc / rho * 1.0e6
      else
          depv_dt_mp_allproc = 0.0_RKIND
      end if


      !***********************************************************************************************
      ! Calculate friction tendency terms:
      !***********************************************************************************************
 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Friction terms are essentially the vorticity tendency due to friction. Need to use the u, v, w
      ! tendencies 
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      ! (1) Uncouple the u tendency on cell edges from dynamics and add to the u tendency from physics   ! MW: Why add this to physics?
      !       MW note: tend_u_euler already decoupled from mass
      tend_u_phys(:,:) = tend_u_phys(:,:) + tend_u_euler(:,:) 
 
      ! (2) Reconstruct tend_u_phys to cell center and calculate vertical derivative of 
      !     tend_uZonal and tend_uMerid
 
      call mpas_reconstruct(mesh, tend_u_phys, tenduX, tenduY, tenduZ, tend_uZonal, tend_uMerid)
 
      call calc_vertDeriv(tend_uZonal, nCells, nVertLevels, dzu, duZonalTend_dz)
      call calc_vertDeriv(tend_uMerid, nCells, nVertLevels, dzu, duMeridTend_dz)
 
      ! (3) Interpolate w tendency from dynamics to cell center, uncouple from density, calclulate gradient 
      !     of tend_wCell on edges, and then reconstruct to get zonal and meridional gradients at cell center:
      !     dwTend_dx, dwTend_dy  
      !
      !     Note: Nick's procedure interpolated density to w levels before uncoupling, but this method 
      !           is consistent with the procedure in the vorticity calculation 
 
      call interp_wLev_thetaLev(tend_w_euler, nCells, nVertLevels, tend_wCell)
    
      ! MW note: tend_w_euler from w_tend_diff is now decoupled before vertical interpolation 
      !tend_wCell(:,:) = tend_wCell(:,:) / rho_prev(:,:)   ! MW: I think this needs a rho_zz term
      
      !call mpas_log_write("maxval tend_wCell is $r", realArgs=(/real(MAXVAL(tend_wCell), kind=RKIND)/))
 
      call calc_gradOnEdges_reconCellCenter(tend_wCell, nCells, nEdges, nVertLevels, nEdgesOnCell, &
                                edgesOnCell, edgesOnCell_sign, cellsOnEdge, dcEdge, coeffs_reconstruct, &
                                latCell, lonCell, dWtend_dxZonal, dWtend_dyMerid)
 
      ! (4) Calculate vertical vorticity tendency on vertices and reconstruct to cell center
      !     Note: the analogy of this procedure in the ertel_pv calculation is missing because vertical vorticity
      !           is already computed and output in MPAS 
       
      call calc_vertical_curl(tend_u_phys, nEdges, nVertices, dcEdge, areaTriangle, verticesOnEdge, uTend_curl)
 
      call interp_absVertVort(uTend_curl, nCells, nEdgesOnCell, verticesOnCell, &
                               cellsOnVertex, areaCell, kiteAreasOnVertex, vertVortTend)
 
      ! (5) Combine three components into vorticity tendency vector
 
      vortTend3D(:,:,1) = dWtend_dyMerid - duMeridTend_dz        ! dFz/dy - dFy/dz
      vortTend3D(:,:,2) = duZonalTend_dz - dWtend_dxZonal        ! dFx/dz - dFz/dx
      vortTend3D(:,:,3) = vertVortTend                           ! dFy/dy - dFx/dx 
 
      call calc_dotProduct_3D(gradTheta, vortTend3D, nCells, nVertLevels, depv_dt_fric)
 
      !depv_dt_fric = depv_dt_fric / rho_prev * 1.0e6
      depv_dt_fric = depv_dt_fric / rho * 1.0e6

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Density tendency contributions to PV budget from discretized PV tendency equation
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      call mpas_log_write("through everything except calc_rho_tend")

      call calc_rho_tend(rho, rho_prev, config_dt, rho_tend)
      call mpas_log_write("through calc_rho_tend")

      depv_dt_rho(:,:) = -1 * ertel_pv_prev(:,:)/rho(:,:) * rho_tend(:,:)  
      call mpas_log_write("through depv_dt_rho")


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
      deallocate(dWtend_dxZonal)
      deallocate(dWtend_dyMerid)
      deallocate(duZonalTend_dz)
      deallocate(duMeridTend_dz)
      deallocate(vertVortTend)
      deallocate(tenduX)
      deallocate(tenduY)
      deallocate(tenduZ)
      deallocate(tend_uZonal)
      deallocate(tend_uMerid)
      deallocate(vortTend3D)
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
      deallocate(rho_tend)

   end subroutine calc_pvBudget


   !*********************************************************************************************************************
   ! MW: Adding a subroutine to accumulate the PV budget tendency terms at each time step. Accumulated values are output 
   !     at the user-specific diagnostic output interval, which determines the time-averaging window of the tendencies.
   !********************************************************************************************************************* 

   subroutine acc_pvBudget(mesh, diag)
 
      use mpas_vector_reconstruction
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array
      use mpas_log, only : mpas_log_write
 
      implicit none
      
      type (mpas_pool_type), intent(inout) :: diag
      type (mpas_pool_type), intent(in) :: mesh

      real(kind=RKIND), dimension(:),   pointer :: depv_dt_diab_pv, depv_dt_fric_pv
      real(kind=RKIND), dimension(:),   pointer :: acc_depv_dt_diab_pv, acc_depv_dt_fric_pv
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_lw, depv_dt_sw, depv_dt_bl, depv_dt_cu, depv_dt_mp, depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_diab, depv_dt_fric
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_lw, acc_depv_dt_sw, acc_depv_dt_bl, acc_depv_dt_cu, acc_depv_dt_mp, acc_depv_dt_mix
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_diab, acc_depv_dt_fric
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_adv, acc_depv_dt_adv 

      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_mp_evap_cw, acc_depv_dt_mp_evap_rw, acc_depv_dt_mp_depo_ice, acc_depv_dt_mp_melt_ice
      real(kind=RKIND), dimension(:,:), pointer :: acc_depv_dt_mp_frez_ice, acc_depv_dt_mp_allproc 
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_evap_cw, depv_dt_mp_evap_rw, depv_dt_mp_depo_ice, depv_dt_mp_melt_ice
      real(kind=RKIND), dimension(:,:), pointer :: depv_dt_mp_frez_ice, depv_dt_mp_allproc 


      integer, pointer :: nCells, nVertLevels
      integer :: iCell, k

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(diag, 'depv_dt_diab_pv', depv_dt_diab_pv)
      call mpas_pool_get_array(diag, 'depv_dt_fric_pv', depv_dt_fric_pv)

      call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)
      call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)
      call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)
      call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)
      call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)
      call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)
      call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
      call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
      call mpas_pool_get_array(diag, 'depv_dt_adv', depv_dt_adv)

      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_cw', depv_dt_mp_evap_cw)
      call mpas_pool_get_array(diag, 'depv_dt_mp_evap_rw', depv_dt_mp_evap_rw)
      call mpas_pool_get_array(diag, 'depv_dt_mp_depo_ice', depv_dt_mp_depo_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_melt_ice', depv_dt_mp_melt_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_frez_ice', depv_dt_mp_frez_ice)
      call mpas_pool_get_array(diag, 'depv_dt_mp_allproc', depv_dt_mp_allproc)

      call mpas_pool_get_array(diag, 'acc_depv_dt_diab_pv', acc_depv_dt_diab_pv)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric_pv', acc_depv_dt_fric_pv)

      call mpas_pool_get_array(diag, 'acc_depv_dt_lw', acc_depv_dt_lw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_sw', acc_depv_dt_sw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_bl', acc_depv_dt_bl)
      call mpas_pool_get_array(diag, 'acc_depv_dt_cu', acc_depv_dt_cu)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp', acc_depv_dt_mp)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mix', acc_depv_dt_mix)
      call mpas_pool_get_array(diag, 'acc_depv_dt_diab', acc_depv_dt_diab)
      call mpas_pool_get_array(diag, 'acc_depv_dt_fric', acc_depv_dt_fric)
      call mpas_pool_get_array(diag, 'acc_depv_dt_adv', acc_depv_dt_adv)

      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_evap_cw', acc_depv_dt_mp_evap_cw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_evap_rw', acc_depv_dt_mp_evap_rw)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_depo_ice', acc_depv_dt_mp_depo_ice)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_melt_ice', acc_depv_dt_mp_melt_ice)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_frez_ice', acc_depv_dt_mp_frez_ice)
      call mpas_pool_get_array(diag, 'acc_depv_dt_mp_allproc', acc_depv_dt_mp_allproc)

      do iCell=1,nCells
         acc_depv_dt_diab_pv(iCell) = acc_depv_dt_diab_pv(iCell) + depv_dt_diab_pv(iCell)
         acc_depv_dt_fric_pv(iCell) = acc_depv_dt_fric_pv(iCell) + depv_dt_fric_pv(iCell)

         do  k=1,nVertLevels
              acc_depv_dt_lw(k,iCell) = acc_depv_dt_lw(k,iCell) + depv_dt_lw(k,iCell)
              acc_depv_dt_sw(k,iCell) = acc_depv_dt_sw(k,iCell) + depv_dt_sw(k,iCell)
              acc_depv_dt_bl(k,iCell) = acc_depv_dt_bl(k,iCell) + depv_dt_bl(k,iCell)
              acc_depv_dt_cu(k,iCell) = acc_depv_dt_cu(k,iCell) + depv_dt_cu(k,iCell)
              acc_depv_dt_mp(k,iCell) = acc_depv_dt_mp(k,iCell) + depv_dt_mp(k,iCell)
              acc_depv_dt_mix(k,iCell) = acc_depv_dt_mix(k,iCell) + depv_dt_mix(k,iCell)
              acc_depv_dt_diab(k,iCell) = acc_depv_dt_diab(k,iCell) + depv_dt_diab(k,iCell)
              acc_depv_dt_fric(k,iCell) = acc_depv_dt_fric(k,iCell) + depv_dt_fric(k,iCell)
              acc_depv_dt_adv(k,iCell) = acc_depv_dt_adv(k,iCell) + depv_dt_adv(k,iCell) 
     
              acc_depv_dt_mp_evap_cw(k,iCell) = acc_depv_dt_mp_evap_cw(k,iCell) + depv_dt_mp_evap_cw(k,iCell) 
              acc_depv_dt_mp_evap_rw(k,iCell) = acc_depv_dt_mp_evap_rw(k,iCell) + depv_dt_mp_evap_rw(k,iCell) 
              acc_depv_dt_mp_depo_ice(k,iCell) = acc_depv_dt_mp_depo_ice(k,iCell) + depv_dt_mp_depo_ice(k,iCell) 
              acc_depv_dt_mp_melt_ice(k,iCell) = acc_depv_dt_mp_melt_ice(k,iCell) + depv_dt_mp_melt_ice(k,iCell) 
              acc_depv_dt_mp_frez_ice(k,iCell) = acc_depv_dt_mp_frez_ice(k,iCell) + depv_dt_mp_frez_ice(k,iCell) 
              acc_depv_dt_mp_allproc(k,iCell) = acc_depv_dt_mp_allproc(k,iCell) + depv_dt_mp_allproc(k,iCell) 
         end do
      end do 

   end subroutine acc_pvBudget

   !*********************************************************************************************************************
   ! MC: Modified subroutine to call PV tendency calculations and interpolation of tendency fields to dynamic 
   !     tropopause at beginning of time step 
   !********************************************************************************************************************* 
   
   subroutine atm_compute_pvBudget_diagnostics(configs, state, time_lev, diag, mesh, tend, tend_physics, diag_physics)
      
      use mpas_constants
      use mpas_derived_types, only : field2DReal
      use mpas_pool_routines, only : mpas_pool_get_dimension, mpas_pool_get_array, mpas_pool_get_field
      use mpas_dmpar, only : mpas_dmpar_exch_halo_field
      use mpas_log, only : mpas_log_write
 
      implicit none
      
      type (mpas_pool_type), intent(in) :: configs
      type (mpas_pool_type), intent(inout) :: diag, tend
      type (mpas_pool_type), intent(in) :: state, mesh, tend_physics, diag_physics
      integer, intent(in) :: time_lev                                  
   
      integer :: iCell, k
      integer, pointer :: nCells, nVertLevels, index_qv
      real (kind=RKIND) :: pvuVal, missingVal 
      real (kind=RKIND), dimension(:,:), pointer :: dtheta_dt_mix, tend_theta_euler 
      type (field2DReal), pointer :: rthratenlw_f, rthratensw_f, rthcuten_f, rthblten_f, dtheta_dt_mp_f, dtheta_dt_mix_f    
      type (field2DReal), pointer :: tend_u_phys_f, tend_u_euler_f, tend_w_euler_f, tend_wCell_f, uTend_curl_f

      type (field2DReal), pointer ::  tend_theta_mp_evap_cw_f, tend_theta_mp_evap_rw_f,  tend_theta_mp_depo_ice_f, & 
                                      tend_theta_mp_frez_ice_f,  tend_theta_mp_melt_ice_f

      real (kind=RKIND), dimension(:,:,:), pointer :: scalars

      ! need halo communication for beginning of time step variables for PV tendency calculations 
      type (field2DReal), pointer :: theta_prev_f !, uReconstructZonal_prev_f, uReconstructMeridional_prev_f, wCell_prev_f, pv_vertex_prev_f

      ! also need variables from end of time step from PV tendency calcs...
      type (field2DReal), pointer :: uReconstructZonal_f, uReconstructMeridional_f, wCell_f, pv_vertex_f
      !!!!!!!
      type (field2DReal), pointer :: theta_dyn_f, w_dyn_f, uReconstructZonal_dyn_f, uReconstructMeridional_dyn_f 

      type (field2DReal), pointer :: depv_dt_adv_f

      call mpas_pool_get_dimension(mesh, 'nCells', nCells)
      call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(state, 'index_qv', index_qv)
      call mpas_pool_get_array(state, 'scalars', scalars, time_lev)

      !need halo cells for everything w/ horizontal derivative             
      !Dtheta/Dt
      call mpas_pool_get_array(tend, 'theta_euler', tend_theta_euler)
      call mpas_pool_get_array(diag, 'dtheta_dt_mix', dtheta_dt_mix)       

      dtheta_dt_mix(:,:) = 0.0_RKIND

      ! NS: With modified moist potential temperature being the model state variable being mixed,
      !     assume qv field is not mixed and so there's no tend_qv to consider
      ! MC: This variable still needs to be checked as of 01/04/23, but horizontal mixing contribution to PV tendency is 
      !     much smaller than other terms, so it's likely fine regardless of mixing of moisture 
      !     Are the theta tendencies and scalars valid at the same time here? 
      !
      !     03/01/2023 - Nick converted theta_m to theta but still needs to be decoupled from rho/zz!

      do iCell = 1,nCells
         do k = 1,nVertLevels
            dtheta_dt_mix(k,iCell) = tend_theta_euler(k,iCell) / (1._RKIND + rvord * scalars(index_qv,k,iCell))
         end do
      end do

      !call mpas_pool_get_field(diag, 'uReconstructZonal_prev', uReconstructZonal_prev_f)           ! MC added
      !call mpas_pool_get_field(diag, 'uReconstructMeridional_prev', uReconstructMeridional_prev_f) ! MC added
      call mpas_pool_get_field(diag, 'theta_prev', theta_prev_f)
      !call mpas_pool_get_field(diag, 'wCell_prev', wCell_prev_f)                                   ! MC added 
      !call mpas_pool_get_field(diag, 'pv_vertex_prev', pv_vertex_prev_f)                           ! MC added

      call mpas_pool_get_field(diag, 'uReconstructZonal', uReconstructZonal_f)           ! MC added
      call mpas_pool_get_field(diag, 'uReconstructMeridional', uReconstructMeridional_f) ! MC added
      call mpas_pool_get_field(diag, 'wCell', wCell_f)                                   ! MC added 
      call mpas_pool_get_field(diag, 'pv_vertex', pv_vertex_f)                           ! MC added    
 
      call mpas_dmpar_exch_halo_field(theta_prev_f)
      !call mpas_dmpar_exch_halo_field(uReconstructZonal_prev_f)                               ! MC added
      !call mpas_dmpar_exch_halo_field(uReconstructMeridional_prev_f)                          ! MC added
      !call mpas_dmpar_exch_halo_field(wCell_prev_f)                                           ! MC added
      !call mpas_dmpar_exch_halo_field(pv_vertex_prev_f)                                       ! MC added

      call mpas_dmpar_exch_halo_field(uReconstructZonal_f)                               ! MC added
      call mpas_dmpar_exch_halo_field(uReconstructMeridional_f)                          ! MC added
      call mpas_dmpar_exch_halo_field(wCell_f)                                           ! MC added
      call mpas_dmpar_exch_halo_field(pv_vertex_f)                                       ! MC added


      call mpas_pool_get_field(tend_physics, 'rthratenlw', rthratenlw_f)
      call mpas_pool_get_field(tend_physics, 'rthratensw', rthratensw_f)
      call mpas_pool_get_field(tend_physics, 'rthcuten', rthcuten_f)
      call mpas_pool_get_field(tend_physics, 'rthblten', rthblten_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mp', dtheta_dt_mp_f)
      call mpas_pool_get_field(diag, 'dtheta_dt_mix', dtheta_dt_mix_f)
 
      call mpas_dmpar_exch_halo_field(rthratenlw_f)
      call mpas_dmpar_exch_halo_field(rthratensw_f)
      call mpas_dmpar_exch_halo_field(rthcuten_f)
      call mpas_dmpar_exch_halo_field(rthblten_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mp_f)
      call mpas_dmpar_exch_halo_field(dtheta_dt_mix_f)
  
      !specific microphysics processes
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

      !friction
      call mpas_pool_get_field(diag, 'tend_u_phys', tend_u_phys_f)
      !MW call mpas_pool_get_field(tend, 'u_euler', tend_u_euler_f)
      !MW call mpas_pool_get_field(tend, 'w_euler', tend_w_euler_f)
      call mpas_pool_get_field(diag, 'u_tend_diff', tend_u_euler_f) ! MW: using ru explicit horiz mixing term from tendencies
      call mpas_pool_get_field(diag, 'w_tend_diff', tend_w_euler_f) ! MW: using rw explicit horiz mixing term from tendencies
      call mpas_pool_get_field(diag, 'tend_wCell', tend_wCell_f)
      call mpas_pool_get_field(diag, 'uTend_curl', uTend_curl_f)

      call mpas_dmpar_exch_halo_field(tend_u_phys_f)
      call mpas_dmpar_exch_halo_field(tend_u_euler_f)
      call mpas_dmpar_exch_halo_field(tend_w_euler_f)
      call mpas_dmpar_exch_halo_field(tend_wCell_f)
      call mpas_dmpar_exch_halo_field(uTend_curl_f)

      call mpas_log_write("beginning calc pv budget")
      call calc_pvBudget(configs, state, time_lev, diag, mesh, tend, tend_physics, diag_physics)
      call mpas_log_write("done with calc pv budget")

      pvuVal = 2.0_RKIND
      missingVal = -99999.0_RKIND

      ! need to interpolate PV tendencies to dynamic tropopause at BEGINNING of time step!
      call interp_pvBudget_diagnostics(mesh, diag, pvuVal, missingVal)
      call mpas_log_write("done with interp pv budget")

      ! MW: added subroutine to accumulate the tendencies for output purposes
      !call mpas_pool_get_field(diag, 'depv_dt_adv', depv_dt_adv_f)
      !call mpas_dmpar_exch_halo_field(depv_dt_adv_f)
      !call mpas_log_write("called and exchanged advection tend?")
      call acc_pvBudget(mesh, diag)
      call mpas_log_write("done with acc pv budget")
   
   end subroutine atm_compute_pvBudget_diagnostics


!=================================================================================================================
end module pv_diagnostics
!=================================================================================================================
