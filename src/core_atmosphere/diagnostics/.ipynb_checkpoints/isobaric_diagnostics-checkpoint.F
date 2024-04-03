! Copyright (c) 2016,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module isobaric_diagnostics

    use mpas_dmpar
    use mpas_kind_types
    use mpas_derived_types
    use mpas_pool_routines
    use mpas_constants, only: rvord, r_earth=>a
    use mpas_log, only : mpas_log_write

    type (MPAS_pool_type), pointer :: mesh
    type (MPAS_pool_type), pointer :: state
    type (MPAS_pool_type), pointer :: diag
    type (MPAS_pool_type), pointer :: diag_physics
    type (MPAS_pool_type), pointer :: tend_physics
    type (MPAS_pool_type), pointer :: configs 

    type (MPAS_clock_type), pointer :: clock

    public :: isobaric_diagnostics_setup, &
              isobaric_diagnostics_compute, &
              isobaric_diagnostics_pv

    private

    logical :: need_mslp, &
               need_relhum, need_dewpoint, need_qvapor, need_refl10cm, need_temperature, need_w, need_uzonal, need_umeridional, &
               need_vorticity, need_height, &
               need_geoheight, need_isobaric, &
               ! PV variables
               need_ertelpv, need_lw, need_sw, need_bl, need_cu, need_mp, need_mix, &
               need_fric_mix, need_fric_bl, need_fric_cu, &
               need_diab, need_fric, need_dyn, &
               need_evap_rw, need_evap_cw, need_depo, need_melt, need_frez, need_mp_all, &
               need_pv_diags, need_pvtend_isobaric, &
               ! Latent heating rates 
               need_dtheta_dt_cu, need_dtheta_dt_mp, need_thtend_isobaric  

    contains


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_setup
    !
    !> \brief Set up the isobaric diagnostics module
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  This routine sets up the isobaric diagnostics module, principally by
    !>  saving pointers to pools that are used in the computation of diagnostics.
    !>  MC: added configs input
    !-----------------------------------------------------------------------
    subroutine isobaric_diagnostics_setup(configs_in, all_pools, simulation_clock)

        use mpas_derived_types, only : MPAS_pool_type, MPAS_clock_type
        use mpas_pool_routines, only : mpas_pool_get_subpool

        implicit none
      
        type (MPAS_pool_type), pointer :: configs_in
        type (MPAS_pool_type), pointer :: all_pools
        type (MPAS_clock_type), pointer :: simulation_clock

        clock => simulation_clock

        call mpas_pool_get_subpool(all_pools, 'mesh', mesh)
        call mpas_pool_get_subpool(all_pools, 'state', state)
        call mpas_pool_get_subpool(all_pools, 'diag', diag)
        call mpas_pool_get_subpool(all_pools, 'diag_physics', diag_physics)
        call mpas_pool_get_subpool(all_pools, 'tend_physics', tend_physics)
   
        configs => configs_in

    end subroutine isobaric_diagnostics_setup


    !-----------------------------------------------------------------------
    !  routine isobaric_diagnostics_compute
    !
    !> \brief Compute isobaric diagnostic before model output is written
    !> \author Michael Duda
    !> \date   21 October 2016
    !> \details
    !>  Compute isobaric diagnostic before model output is written. Code called
    !>  from here was previously in mpas_atm_interp_diagnostics.F.
    !
    !-----------------------------------------------------------------------
    subroutine isobaric_diagnostics_compute()

        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written

        implicit none

        logical :: need_any_diags

        need_any_diags = .false.
        
        need_isobaric = .false.

        need_mslp = MPAS_field_will_be_written('mslp')
        need_any_diags = need_any_diags .or. need_mslp
        
        need_relhum =  MPAS_field_will_be_written('relhum_isobaric')
        need_isobaric = need_isobaric .or. need_relhum
        need_dewpoint =  MPAS_field_will_be_written('dewpoint_isobaric')
        need_isobaric = need_isobaric .or. need_dewpoint
        need_qvapor =  MPAS_field_will_be_written('qvapor_isobaric')
        need_isobaric = need_isobaric .or. need_qvapor
        need_refl10cm =  MPAS_field_will_be_written('refl10cm_isobaric')
        need_isobaric = need_isobaric .or. need_refl10cm
        need_temperature = MPAS_field_will_be_written('temperature_isobaric')
        need_isobaric = need_isobaric .or. need_temperature
        need_height = MPAS_field_will_be_written('height_isobaric')
        need_isobaric = need_isobaric .or. need_height
        need_geoheight = MPAS_field_will_be_written('geoheight_isobaric')
        need_isobaric = need_isobaric .or. need_geoheight
        need_uzonal = MPAS_field_will_be_written('uzonal_isobaric')
        need_isobaric = need_isobaric .or. need_uzonal            
        need_umeridional = MPAS_field_will_be_written('umeridional_isobaric')
        need_isobaric = need_isobaric .or. need_umeridional
        need_w = MPAS_field_will_be_written('w_isobaric')
        need_isobaric = need_isobaric .or. need_w
        need_vorticity = MPAS_field_will_be_written('vorticity_isobaric')
        need_isobaric = need_isobaric .or. need_vorticity
        need_any_diags = need_any_diags .or. need_isobaric


        if (need_any_diags) then
            call interp_diagnostics(mesh, state, 1, diag, diag_physics)
        end if
   
    end subroutine isobaric_diagnostics_compute

   !==================================================================================================
    subroutine isobaric_diagnostics_pv()
    ! MC_TODO: clean up
   !==================================================================================================
   ! MC: add subroutine to be called by pv_diagnostics.F to interpolate PV variables and tendencies
   !     to isobaric levels prior to accumulating 
   
        use mpas_atm_diagnostics_utils, only : MPAS_field_will_be_written
        use mpas_pool_routines, only: mpas_pool_get_config

        implicit none

        logical, pointer :: config_pv_isobaric, config_pv_microphys
    
        need_pv_diags = .false.           ! PV field only
        need_pvtend_isobaric = .false.    ! PV tendencies 
        need_thtend_isobaric = .false.    ! Latent heating tendencies 

        call mpas_log_write('In isobaric_diagnostics_pv routine.')

        call mpas_pool_get_config(configs, 'config_pv_isobaric', config_pv_isobaric)
        call mpas_pool_get_config(configs, 'config_pv_microphys', config_pv_microphys)

        call mpas_log_write('config_pv_isobaric is: $l', logicArgs=(/config_pv_isobaric/))
        call mpas_log_write('config_pv_microphys is: $l', logicArgs=(/config_pv_microphys/))

        ! MC_TODO: do these flags do anything regarding writing vars? They are evaluated below...
        if (config_pv_isobaric) then 
        
            ! PV field 
            need_ertelpv =  MPAS_field_will_be_written('ertel_pv_isobaric')
            need_pv_diags = need_pv_diags .or. need_ertelpv

            ! tendencies 
            need_lw = MPAS_field_will_be_written('depv_dt_lw_isobaric') 
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_lw
            need_sw = MPAS_field_will_be_written('depv_dt_sw_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_sw
            need_bl = MPAS_field_will_be_written('depv_dt_bl_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_bl
            need_cu = MPAS_field_will_be_written('depv_dt_cu_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_cu
            need_mp = MPAS_field_will_be_written('depv_dt_mp_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_mp
            need_mix = MPAS_field_will_be_written('depv_dt_mix_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_mix
            need_fric_mix = MPAS_field_will_be_written('depv_dt_fric_mix_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_fric_mix
            need_fric_bl = MPAS_field_will_be_written('depv_dt_fric_bl_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_fric_bl
            need_fric_cu = MPAS_field_will_be_written('depv_dt_fric_cu_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_fric_cu
            need_diab = MPAS_field_will_be_written('depv_dt_diab_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_diab
            need_fric = MPAS_field_will_be_written('depv_dt_fric_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_fric
            need_dyn = MPAS_field_will_be_written('depv_dt_dyn_isobaric')
            need_pvtend_isobaric = need_pvtend_isobaric .or. need_dyn


            if (config_pv_microphys) then
                need_evap_rw = MPAS_field_will_be_written('depv_dt_mp_evap_rw_isobaric')
                need_pvtend_isobaric = need_pvtend_isobaric .or. need_evap_rw
                need_evap_cw = MPAS_field_will_be_written('depv_dt_mp_evap_cw_isobaric')
                need_pvtend_isobaric = need_pvtend_isobaric .or. need_evap_cw
                need_depo = MPAS_field_will_be_written('depv_dt_mp_depo_ice_isobaric')
                need_pvtend_isobaric = need_pvtend_isobaric .or. need_depo
                need_melt = MPAS_field_will_be_written('depv_dt_mp_melt_ice_isobaric')
                need_pvtend_isobaric = need_pvtend_isobaric .or. need_melt
                need_frez = MPAS_field_will_be_written('depv_dt_mp_frez_ice_isobaric')
                need_pvtend_isobaric = need_pvtend_isobaric .or. need_frez
                need_mp_all = MPAS_field_will_be_written('depv_dt_mp_allproc_isobaric')
                need_pvtend_isobaric = need_pvtend_isobaric .or. need_mp_all
            end if 

                
            need_dtheta_dt_cu = MPAS_field_will_be_written('dtheta_dt_cu_isobaric')
            need_thtend_isobaric = need_thtend_isobaric .or. need_dtheta_dt_cu
            need_dtheta_dt_mp = MPAS_field_will_be_written('acc_dtheta_dt_mp_isobaric')
            need_thtend_isobaric = need_thtend_isobaric .or. need_dtheta_dt_mp

            need_pv_diags = need_pv_diags .or. need_pvtend_isobaric .or. need_thtend_isobaric

            ! MC: removed if condition here because statement didn't seem to get activated... not sure why
            !if (need_pv_diags) then
            call mpas_log_write('Calling interp_pv_diag_isobaric.')
            call interp_pv_diag_isobaric(mesh, state, 1, diag, tend_physics)
        end if

     end subroutine isobaric_diagnostics_pv


   !==================================================================================================
    subroutine interp_diagnostics(mesh, state, time_lev, diag, diag_physics)
   !==================================================================================================

        implicit none 
        
       !input arguments:
        type (mpas_pool_type), intent(in)  :: mesh
        type (mpas_pool_type), intent(in) :: state
        type (mpas_pool_type), intent(in) :: diag_physics 
        integer, intent(in) :: time_lev              ! which time level to use from state
       
       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag
       
       !local variables:
        integer :: iCell, iVert, iVertD, k, kk
        integer, pointer :: nCells, nVertLevels, nVertices, vertexDegree, nIsoLevels, nVertLevelsP1
        integer, pointer :: index_qv, num_scalars 
       
        type (field2DReal), pointer:: pressure_p_field, vorticity_f
       
        integer, dimension(:), pointer :: nEdgesOnCell
        integer, dimension(:,:), pointer :: verticesOnCell, cellsOnVertex
        real (kind=RKIND), dimension(:), pointer :: areaTriangle, areaCell
        real (kind=RKIND), dimension(:,:), pointer :: kiteAreasOnVertex
        
        real (kind=RKIND), dimension(:,:), pointer :: exner, height
        real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
        real (kind=RKIND), dimension(:,:), pointer :: relhum, theta_m, vorticity
        real (kind=RKIND), dimension(:,:), pointer :: qvapor, refl10cm 
        real (kind=RKIND), dimension(:,:), pointer :: umeridional, uzonal, vvel
        real (kind=RKIND), dimension(:,:,:), pointer :: scalars
       
        real (kind=RKIND), dimension(:), pointer :: iso_levels
        
        real (kind=RKIND), dimension(:,:), pointer :: temperature_isobaric, relhum_isobaric, &
                                                      dewpoint_isobaric, uzonal_isobaric, &
                                                      umeridional_isobaric, height_isobaric, &
                                                      geoheight_isobaric, w_isobaric, &
                                                      vorticity_isobaric, qvapor_isobaric, refl10cm_isobaric 
                                                      
        real (kind=RKIND) :: evp
       
       !--------------------
       
        real (kind=RKIND), dimension(:), pointer :: mslp
       
        real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1, pressure2, pressure_v, temperature
        real (kind=RKIND), dimension(:,:), allocatable :: dewpoint, vorticity_cell
       
       !local interpolated fields:
        real (kind=RKIND) :: w1,w2,z0,z1,z2
        real (kind=RKIND), dimension(:,:), allocatable :: field_in, press_in, press_in2
        real (kind=RKIND), dimension(:,:), allocatable :: field_interp, press_interp
 
       !--------------------------------------------------------------------------------------------------
       
       ! call mpas_log_write('')
       ! call mpas_log_write('--- enter subroutine interp_diagnostics:')
       
        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nVertLevelsP1', nVertLevelsP1)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nVertices', nVertices)
        call mpas_pool_get_dimension(mesh, 'vertexDegree', vertexDegree)
        call mpas_pool_get_dimension(mesh, 'nIsoLevels', nIsoLevels)
        call mpas_pool_get_dimension(state, 'index_qv', index_qv)
        call mpas_pool_get_dimension(state, 'num_scalars', num_scalars)
       
        call mpas_pool_get_array(mesh, 'nEdgesOnCell', nEdgesOnCell)
        call mpas_pool_get_array(mesh, 'verticesOnCell', verticesOnCell)
        call mpas_pool_get_array(mesh, 'cellsOnVertex', cellsOnVertex)
        call mpas_pool_get_array(mesh, 'areaTriangle', areaTriangle)
        call mpas_pool_get_array(mesh, 'areaCell', areaCell)
        call mpas_pool_get_array(mesh, 'kiteAreasOnVertex', kiteAreasOnVertex)
       
        call mpas_pool_get_array(diag, 'iso_levels', iso_levels)
        
        call mpas_pool_get_array(mesh, 'zgrid', height)
        call mpas_pool_get_array(state, 'w', vvel, time_lev)
        call mpas_pool_get_array(state, 'theta_m', theta_m, time_lev)
        call mpas_pool_get_array(state, 'scalars', scalars, time_lev)
       
        call mpas_pool_get_array(diag, 'exner', exner)
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
        call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
        call mpas_pool_get_array(diag, 'vorticity', vorticity)
        call mpas_pool_get_array(diag, 'uReconstructMeridional', umeridional)
        call mpas_pool_get_array(diag, 'uReconstructZonal', uzonal)
        call mpas_pool_get_array(diag, 'relhum', relhum)
        call mpas_pool_get_array(diag, 'qvapor', qvapor)
        call mpas_pool_get_array(diag_physics, 'refl10cm_3D', refl10cm)
               
        call mpas_pool_get_array(diag, 'mslp', mslp)
        
        call mpas_pool_get_field(diag, 'pressure_p', pressure_p_field)
        call mpas_dmpar_exch_halo_field(pressure_p_field)
        
        call mpas_pool_get_field(diag, 'vorticity', vorticity_f)                                     
        call mpas_dmpar_exch_halo_field(vorticity_f)
        
        ! Variables to interpolate as common arrays:
        call mpas_pool_get_array(diag, 'temperature_isobaric', temperature_isobaric)
        call mpas_pool_get_array(diag, 'relhum_isobaric', relhum_isobaric)
        call mpas_pool_get_array(diag, 'refl10cm_isobaric', refl10cm_isobaric)
        call mpas_pool_get_array(diag, 'dewpoint_isobaric', dewpoint_isobaric)
        call mpas_pool_get_array(diag, 'qvapor_isobaric', qvapor_isobaric)
        call mpas_pool_get_array(diag, 'uzonal_isobaric', uzonal_isobaric)
        call mpas_pool_get_array(diag, 'umeridional_isobaric', umeridional_isobaric)
        call mpas_pool_get_array(diag, 'height_isobaric', height_isobaric)
        call mpas_pool_get_array(diag, 'geoheight_isobaric', geoheight_isobaric)
        call mpas_pool_get_array(diag, 'w_isobaric', w_isobaric)
        call mpas_pool_get_array(diag, 'vorticity_isobaric', vorticity_isobaric)
 
        if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
        if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )
        if(.not.allocated(pressure2)   ) allocate(pressure2(nVertLevelsP1,nCells)   )
        if(.not.allocated(pressure_v)  ) allocate(pressure_v(nVertLevels,nVertices) )
        if(.not.allocated(temperature) ) allocate(temperature(nVertLevels,nCells)   )
        if(.not.allocated(dewpoint) )    allocate(dewpoint(nVertLevels,nCells)   )
        
        temperature(:,:) = 0.0
        dewpoint(:,:) = 0.0

        ! Definition of isobaric levels.
        iso_levels(:) = (/10000.0, 12500.0, 15000.0, 17500.0, 20000.0, 22500.0, 25000.0, 27500.0, 30000.0, &
                          32500.0, 35000.0, 40000.0, 45000.0, 50000.0, 55000.0, 60000.0, 65000.0, 70000.0, &
                          75000.0, 77500.0, 80000.0, 82500.0, 85000.0, 87500.0, 90000.0, 92500.0, 95000.0, 100000.0/)
       
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !calculation of total pressure at cell centers (at mass points):
        do iCell = 1, nCells
            do k = 1, nVertLevels
               pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
               pressureCp1(k,iCell) = pressure(k,iCell)
            end do
        end do
        do iCell = nCells+1,nCells+1
            do k =1,nVertLevels
               pressureCp1(k,iCell) = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
            end do
        end do
       
       !calculation of total pressure at cell centers (at vertical velocity points):
        k = nVertLevelsP1
        do iCell=1,nCells
           z0 = height(k,iCell)
           z1 = 0.5*(height(k,iCell)+height(k-1,iCell)) 
           z2 = 0.5*(height(k-1,iCell)+height(k-2,iCell))
           w1 = (z0-z2)/(z1-z2)
           w2 = 1.-w1
           !use log of pressure to avoid occurrences of negative top-of-the-model pressure.
           pressure2(k,iCell) = exp(w1*log(pressure(k-1,iCell))+w2*log(pressure(k-2,iCell)))
        end do
        
        do k=2,nVertLevels
            do iCell=1,nCells
               w1 = (height(k,iCell)-height(k-1,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
               w2 = (height(k+1,iCell)-height(k,iCell)) / (height(k+1,iCell)-height(k-1,iCell))
               pressure2(k,iCell) = exp(w1*log(pressure(k,iCell)) + w2*log(pressure(k-1,iCell)))
            end do
        end do
        
        k = 1
        do iCell=1,nCells
           z0 = height(k,iCell)
           z1 = 0.5*(height(k,iCell)+height(k+1,iCell)) 
           z2 = 0.5*(height(k+1,iCell)+height(k+2,iCell))
           w1 = (z0-z2)/(z1-z2)
           w2 = 1.-w1
           pressure2(k,iCell) = exp(w1*log(pressure(k,iCell))+w2*log(pressure(k+1,iCell)))
        end do
       
       !calculation of total pressure at cell vertices (at mass points):
        do iVert = 1, nVertices
           pressure_v(:,iVert) = 0._RKIND
       
           do k=1,nVertLevels
               do iVertD = 1, vertexDegree
                  pressure_v(k,iVert) = pressure_v(k,iVert) &
                      + kiteAreasOnVertex(iVertD,iVert)*pressureCp1(k,cellsOnVertex(iVertD,iVert))
               end do
               pressure_v(k,iVert) = pressure_v(k,iVert) / areaTriangle(iVert)
           end do
        end do
       
       
        if (need_isobaric .or. need_relhum .or. need_dewpoint .or. need_temperature .or. need_mslp) then
           !calculation of temperature at cell centers:
            do iCell=1,nCells
               do k=1,nVertLevels
                  temperature(k,iCell) = (theta_m(k,iCell)/(1._RKIND+rvord*scalars(index_qv,k,iCell)))*exner(k,iCell) 

                  ! Vapor pressure (NB: pressure here is already in hPa)
                  evp = pressure(k,iCell) * scalars(index_qv,k,iCell) / (scalars(index_qv,k,iCell) + 0.622_RKIND)
                  evp = max(evp, 1.0e-8_RKIND)

                  ! Dewpoint temperature following Bolton (1980)
                  dewpoint(k,iCell) = (243.5_RKIND * log(evp/6.112_RKIND)) / (17.67_RKIND - log(evp/6.112_RKIND))
                  dewpoint(k,iCell) = dewpoint(k,iCell) + 273.15
              end do
            end do
        end if
       
       
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!! Interpolate fields to array of pressure levels !!!!!!!!!!!
        
        if (need_isobaric) then
     
            !if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
            !if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
            !if(.not.allocated(press_in2)) allocate(press_in2(nCells,nVertLevelsP1))
            !if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
            if(.not.allocated(press_interp)) allocate(press_interp(nCells, nIsoLevels))
            
            ! populate array with pressure levels for interpolation [in Pa]
            do k=1,nIsoLevels
               press_interp(:,k) = iso_levels(k)
            end do
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! reverse the vertical axis of pressure array and put back into Pa
         !   do iCell=1,nCells
         !       ! for mass levels
         !       do k=1,nVertLevels
         !           kk=nVertLevels+1-k
         !           press_in(iCell,kk) = pressure(k,iCell) * 100.
         !       end do
                
         !       ! for w levels
         !       do k=1,nVertLevelsP1
         !          kk = nVertLevelsP1+1-k
         !          press_in2(iCell,kk) = pressure2(k,iCell) 
         !       end do
                
         !   end do
         
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate temperature:
            if (need_temperature) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = temperature(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    temperature_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
           end if
           
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate dewpoint:
            if (need_dewpoint) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = dewpoint(k,iCell)
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    dewpoint_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if
            
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate relative humidity:
            if (need_relhum) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = relhum(k,iCell)
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                    end do
                end do
                
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    relhum_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if 
            
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate qv (water vapor mixing ratio):
            if (need_qvapor) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = qvapor(k,iCell)
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                    end do
                end do
                
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    qvapor_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if 
            
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate 10 cm radar reflectivity:
            if (need_refl10cm) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = refl10cm(k,iCell)
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                    end do
                end do
                
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    refl10cm_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if 

           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate uzonal:
            if (need_uzonal) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = uzonal(k,iCell)
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    uzonal_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
             end if
             
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate umeridional:
            if (need_umeridional) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = umeridional(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    umeridional_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if 
                
           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Interpolate vertical velocity:
            if (need_w) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevelsP1
                        kk = nVertLevelsP1+1-k
                        field_in(iCell,kk) = vvel(k,iCell)
                        press_in(iCell,kk) = pressure2(k,iCell) * 100.
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevelsP1, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    w_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
      
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if 
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate vertical vorticity 
            if (need_vorticity) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                if(.not.allocated(vorticity_cell)) allocate(vorticity_cell(nVertLevels,nCells))
 
                vorticity_cell(:,:) = 0.0

               ! first, reconstruct vorticity to cell center (decreases number of points by roughly half) 
                call interp_absVertVort(vorticity, nCells, nEdgesOnCell, verticesOnCell, &
                                cellsOnVertex, areaCell, kiteAreasOnVertex, vorticity_cell)
                                
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        field_in(iCell,kk) = vorticity_cell(k,iCell)
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                    end do
                end do
                
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    vorticity_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                if(allocated(vorticity_cell)) deallocate(vorticity_cell)
            end if 
        
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate geometric height and convert to geopotential height

            if (need_height .or. need_geoheight) then
                !reverse the vertical axis of array
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevelsP1))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevelsP1))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))

                do iCell=1,nCells
                    do k=1,nVertLevelsP1
                        kk = nVertLevelsP1+1-k
                        field_in(iCell,kk) = height(k,iCell) 
                        press_in(iCell,kk) = pressure2(k,iCell) * 100.
                        !field_in(iCell,kk) = 0.5*(height(k,iCell)+height(k+1,iCell))  ! if only looping over nVertLevels
                    end do
                end do
                call interp_tofixed_pressure(nCells, nVertLevelsP1, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    height_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                
                if (need_geoheight) then
                    geoheight_isobaric(:,1:nCells) = (r_earth * height_isobaric(:,1:nCells)) / (r_earth + height_isobaric(:,1:nCells))
                end if     

                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
            end if 
       
        end if 
        if(allocated(press_interp)) deallocate(press_interp)
        if(allocated(pressureCp1)) deallocate(pressureCp1)
        if(allocated(pressure_v)) deallocate(pressure_v)

        if (need_mslp) then
            !... compute SLP (requires temp, height, pressure, qvapor)
             call compute_slp(nCells, nVertLevels, num_scalars, temperature, height, pressure, index_qv, scalars, mslp)
             mslp(:) = mslp(:) * 100.0   ! Convert from hPa to Pa
            !... alternative way
            !do iCell = 1, nCells
            !   mslp(iCell) = diag % surface_pressure % array(iCell) + 11.38*height(1,iCell)
            !   mslp(iCell) = mslp(iCell)/100.
            !enddo
        end if
   
        if(allocated(temperature)) deallocate(temperature)
        if(allocated(pressure2)) deallocate(pressure2)
        if(allocated(pressure)) deallocate(pressure)
        if(allocated(dewpoint)) deallocate(dewpoint)
       
    end subroutine interp_diagnostics


   !==================================================================================================
    subroutine interp_pv_diag_isobaric(mesh, state, time_lev, diag, tend_physics)
   !==================================================================================================

        use mpas_pool_routines, only: mpas_pool_get_config

        implicit none 
        
        logical, pointer :: config_pv_microphys

       !input arguments:
        type (mpas_pool_type), intent(in)  :: mesh
        type (mpas_pool_type), intent(in)  :: tend_physics 
        type (mpas_pool_type), intent(in) :: state
        integer, intent(in) :: time_lev              

       !inout arguments:
        type (mpas_pool_type), intent(inout) :: diag

       !local variables:
        integer :: iCell, k, kk
        integer, pointer :: nCells, nVertLevels, nIsoLevels

        type (field2DReal), pointer:: pressure_p_field !, depv_dt_dyn_f      ! MC: commented _dyn. not sure why needed...

        real (kind=RKIND), dimension(:,:), pointer :: pressure_b, pressure_p 
        
        real (kind=RKIND), dimension(:,:), pointer :: ertel_pv, depv_dt_lw, depv_dt_sw, &
                                                      depv_dt_bl, depv_dt_cu, depv_dt_mp, &
                                                      depv_dt_mix,  &
                                                      depv_dt_fric_bl, depv_dt_fric_cu, depv_dt_fric_mix, &
                                                      depv_dt_diab, depv_dt_fric, depv_dt_dyn, &
                                                      depv_dt_mp_evap_rw, depv_dt_mp_evap_cw, &
                                                      depv_dt_mp_depo_ice, depv_dt_mp_melt_ice, &
                                                      depv_dt_mp_frez_ice, depv_dt_mp_allproc, &
                                                      dtheta_dt_cu, dtheta_dt_mp 

        real (kind=RKIND), dimension(:), pointer :: iso_levels

        real (kind=RKIND), dimension(:,:), pointer :: ertel_pv_isobaric, depv_dt_lw_isobaric, &
                                                      depv_dt_sw_isobaric, depv_dt_bl_isobaric, &
                                                      depv_dt_cu_isobaric, depv_dt_mp_isobaric, &
                                                      depv_dt_mix_isobaric, &
                                                      depv_dt_fric_bl_isobaric, depv_dt_fric_cu_isobaric, &
                                                      depv_dt_fric_mix_isobaric, &
                                                      depv_dt_diab_isobaric, depv_dt_fric_isobaric, &
                                                      depv_dt_dyn_isobaric, &
                                                      depv_dt_mp_evap_rw_isobaric, depv_dt_mp_evap_cw_isobaric, &
                                                      depv_dt_mp_depo_ice_isobaric, depv_dt_mp_melt_ice_isobaric, &
                                                      depv_dt_mp_frez_ice_isobaric, depv_dt_mp_allproc_isobaric, &
                                                      dtheta_dt_cu_isobaric, dtheta_dt_mp_isobaric 
       !--------------------

        real (kind=RKIND), dimension(:,:), allocatable :: pressure, pressureCp1

       !local interpolated fields:
        real (kind=RKIND) :: w1,w2,z0,z1,z2
        real (kind=RKIND), dimension(:,:), allocatable :: field_in, press_in
        real (kind=RKIND), dimension(:,:), allocatable :: field_interp, press_interp

       !--------------------------------------------------------------------------------------------------

        call mpas_pool_get_dimension(mesh, 'nCells', nCells)
        call mpas_pool_get_dimension(mesh, 'nVertLevels', nVertLevels)
        call mpas_pool_get_dimension(mesh, 'nIsoLevels', nIsoLevels)
        call mpas_pool_get_array(diag, 'iso_levels', iso_levels)
        
        call mpas_pool_get_array(diag, 'pressure_base', pressure_b)
        call mpas_pool_get_array(diag, 'pressure_p', pressure_p)
        
        ! PV variables to interpolate
        call mpas_pool_get_array(diag, 'ertel_pv', ertel_pv)
        call mpas_pool_get_array(diag, 'depv_dt_lw', depv_dt_lw)
        call mpas_pool_get_array(diag, 'depv_dt_sw', depv_dt_sw)
        call mpas_pool_get_array(diag, 'depv_dt_bl', depv_dt_bl)
        call mpas_pool_get_array(diag, 'depv_dt_cu', depv_dt_cu)
        call mpas_pool_get_array(diag, 'depv_dt_mp', depv_dt_mp)
        call mpas_pool_get_array(diag, 'depv_dt_mix', depv_dt_mix)
        call mpas_pool_get_array(diag, 'depv_dt_fric_bl', depv_dt_fric_bl)
        call mpas_pool_get_array(diag, 'depv_dt_fric_cu', depv_dt_fric_cu)
        call mpas_pool_get_array(diag, 'depv_dt_fric_mix', depv_dt_fric_mix)
        call mpas_pool_get_array(diag, 'depv_dt_diab', depv_dt_diab)
        call mpas_pool_get_array(diag, 'depv_dt_fric', depv_dt_fric)
        call mpas_pool_get_array(diag, 'depv_dt_dyn', depv_dt_dyn)

        ! microphysics
        call mpas_pool_get_config(configs, 'config_pv_microphys', config_pv_microphys)
        call mpas_pool_get_array(diag, 'depv_dt_mp_evap_rw', depv_dt_mp_evap_rw)
        call mpas_pool_get_array(diag, 'depv_dt_mp_evap_cw', depv_dt_mp_evap_cw)
        call mpas_pool_get_array(diag, 'depv_dt_mp_depo_ice', depv_dt_mp_depo_ice)
        call mpas_pool_get_array(diag, 'depv_dt_mp_melt_ice', depv_dt_mp_melt_ice)
        call mpas_pool_get_array(diag, 'depv_dt_mp_frez_ice', depv_dt_mp_frez_ice)
        call mpas_pool_get_array(diag, 'depv_dt_mp_allproc', depv_dt_mp_allproc)
        call mpas_pool_get_array(diag, 'dtheta_dt_cu', dtheta_dt_cu)            ! MC: variable now defined in pv_diagnostics.F
        call mpas_pool_get_array(diag, 'dtheta_dt_mp', dtheta_dt_mp)            ! MC: variable now defined in pv_diagnostics.F
        
        ! Interpolated variables
        call mpas_pool_get_array(diag, 'ertel_pv_isobaric', ertel_pv_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_lw_isobaric', depv_dt_lw_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_sw_isobaric', depv_dt_sw_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_bl_isobaric', depv_dt_bl_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_cu_isobaric', depv_dt_cu_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_mp_isobaric', depv_dt_mp_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_mix_isobaric', depv_dt_mix_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_fric_bl_isobaric', depv_dt_fric_bl_isobaric)
        call mpas_pool_get_array(diag, 'depv_dt_fric_cu_isobaric', depv_dt_fric_cu_isobaric)
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

        call mpas_pool_get_field(diag, 'pressure_p', pressure_p_field)
        call mpas_dmpar_exch_halo_field(pressure_p_field)

        ! MC -- commenting below... not sure why needed
        !call mpas_pool_get_field(diag, 'depv_dt_dyn', depv_dt_dyn_f)
        !call mpas_dmpar_exch_halo_field(depv_dt_dyn_f)


        ! ---------------------------------------------------------------------------
        ! MC: Explicitly adding these declarations here because they aren't being activated otherwise. 
        ! Needs to be fixed...
        ! MC_TODO: could be moved to top of script...

        need_pv_diags = .true.
        need_pvtend_isobaric = .true.
        need_ertelpv = .true.
        need_lw = .true.
        need_sw = .true.
        need_bl = .true.
        need_cu = .true.
        need_mp = .true.
        need_mix = .true.
        need_fric_mix = .true.
        need_fric_bl = .true.
        need_fric_cu = .true.
        need_diab = .true.
        need_fric = .true.
        need_dyn = .true.
        need_evap_rw = .true.
        need_evap_cw = .true.
        need_depo = .true.
        need_melt = .true.
        need_frez = .true.
        need_mp_all = .true.
        need_dtheta_dt_cu = .true.
        need_dtheta_dt_mp = .true.
        need_thtend_isobaric = .true.
        
        ! ---------------------------------------------------------------------------
        

        if(.not.allocated(pressure)    ) allocate(pressure(nVertLevels,nCells)      )
        if(.not.allocated(pressureCp1) ) allocate(pressureCp1(nVertLevels,nCells+1) )

        call mpas_log_write('Declaring iso_levels for PV variable interpolation.')
        
        iso_levels(:) = (/10000.0, 12500.0, 15000.0, 17500.0, 20000.0, 22500.0, 25000.0, 27500.0, 30000.0, &
                        32500.0, 35000.0, 40000.0, 45000.0, 50000.0, 55000.0, 60000.0, 65000.0, 70000.0, &
                        75000.0, 77500.0, 80000.0, 82500.0, 85000.0, 87500.0, 90000.0, 92500.0, 95000.0, 100000.0/)
                
       !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !calculation of total pressure at cell centers (at mass points):
        do iCell = 1,nCells
            do k = 1,nVertLevels
               pressure(k,iCell)    = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
               pressureCp1(k,iCell) = pressure(k,iCell)
            end do
        end do
        do iCell = nCells+1,nCells+1
            do k =1,nVertLevels
               pressureCp1(k,iCell) = (pressure_p(k,iCell) + pressure_b(k,iCell)) / 100._RKIND
            end do
        end do
       
        ! call mpas_log_write('did pressure calcualtion')
        
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!! Interpolate fields to array of pressure levels !!!!!!!!!!!
        
        call mpas_log_write('Beginning PV interpolation:')
        ! None of these flags are being activated... 
        if (need_pv_diags) then
            call mpas_log_write('needed pvtend_isobaric')
            !if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
            !if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
            !if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
            if(.not.allocated(press_interp)) allocate(press_interp(nCells, nIsoLevels))
     
            ! populate array with pressure levels for interpolation [in Pa]
            do k=1,nIsoLevels
               press_interp(:,k) = iso_levels(k)
            end do
   
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       
            !if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
            !if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
            
            ! reverse the vertical axis of pressure array and put back into Pa
            !do iCell=1,nCells
            !    ! for mass levels
            !    do k=1,nVertLevels
            !        kk=nVertLevels+1-k
            !        press_in(iCell,kk) = pressure(k,iCell) * 100.
            !    end do
                
            !end do
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate PV    
            if (need_ertelpv .and. associated(ertel_pv)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = ertel_pv(k,iCell)
                    end do
                end do
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    ertel_pv_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of ertel_pv.')
            end if
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate LW PV tendency
            
            if (need_lw .and. associated(depv_dt_lw)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_lw(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)
                
                do k=1,nIsoLevels
                    depv_dt_lw_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of LW tendency.')
            end if
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate SW PV tendency
            
            if (need_sw .and. associated(depv_dt_sw)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_sw(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)
                
                do k=1,nIsoLevels
                    depv_dt_sw_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of SW tendency.')
            end if
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate BL PV tendency
            
            if (need_bl .and. associated(depv_dt_bl)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_bl(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)
                
                do k=1,nIsoLevels
                    depv_dt_bl_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of BL tendency.')
            end if
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate CU PV tendency
            
            if (need_cu .and. associated(depv_dt_cu)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_cu(k,iCell)
                    end do
                end do
            
                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)
                
                do k=1,nIsoLevels
                    depv_dt_cu_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of CU tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP PV tendency
            
            if (need_mp .and. associated(depv_dt_mp)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of MP tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate Mixing PV tendency
            
            if (need_mix .and. associated(depv_dt_mix)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mix(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mix_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of mixing tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate Mixing friction PV tendency
            
            if (need_fric_mix .and. associated(depv_dt_fric_mix)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_fric_mix(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_fric_mix_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of friction mixing tendency.')
            end if
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate PBL friction PV tendency
            
            if (need_fric_bl .and. associated(depv_dt_fric_bl)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_fric_bl(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_fric_bl_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of friction PBL tendency.')
            end if


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate cumulus friction PV tendency
            
            if (need_fric_cu .and. associated(depv_dt_fric_cu)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_fric_cu(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_fric_cu_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of friction CU tendency.')
            end if

            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate Diabatic PV tendency
            
            if (need_diab .and. associated(depv_dt_diab)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_diab(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_diab_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of diabatic PV tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate Friction PV tendency
            
            if (need_fric .and. associated(depv_dt_fric)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_fric(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_fric_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of friction PV tendency.')
            end if
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate Dynamics PV tendency
            
            if (need_dyn .and. associated(depv_dt_dyn)) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_dyn(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_dyn_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of dynamics PV tendency.')
            end if


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP - evap_rw PV tendency
            
            if (config_pv_microphys .and. associated(depv_dt_mp_evap_rw)) then 
            !if (need_evap_rw) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp_evap_rw(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_evap_rw_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of evaporation RW PV tendency.')
            end if
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP - evap_cw PV tendency
           
            if (config_pv_microphys .and. associated(depv_dt_mp_evap_cw)) then 
            !if (need_evap_cw) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp_evap_cw(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_evap_cw_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of evaporation CW PV tendency.')
            end if


            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP - depo_ice PV tendency
            
            if (config_pv_microphys .and. associated(depv_dt_mp_depo_ice)) then 
            !if (need_depo) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp_depo_ice(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_depo_ice_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of deposition PV tendency.')
            end if
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP - melt_ice PV tendency
            
            if (config_pv_microphys .and. associated(depv_dt_mp_melt_ice)) then 
            !if (need_melt) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp_melt_ice(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_melt_ice_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of melting PV tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP - frez_ice PV tendency
            
            if (config_pv_microphys .and. associated(depv_dt_mp_frez_ice)) then 
            !if (need_frez) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp_frez_ice(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_frez_ice_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of freezing PV tendency.')
            end if
            
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate MP - all PV tendency
            
            if (config_pv_microphys .and. associated(depv_dt_mp_allproc)) then 
            !if (need_mp_all) then
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = depv_dt_mp_allproc(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    depv_dt_mp_allproc_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of all MP processes PV tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate dtheta_dt_cu
            
            if (need_dtheta_dt_cu .and. associated(dtheta_dt_cu)) then 
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = dtheta_dt_cu(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    dtheta_dt_cu_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of cumulus heating tendency.')
            end if

            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! Interpolate dtheta_dt_mp

            if (need_dtheta_dt_mp .and. associated(dtheta_dt_mp)) then 
                if(.not.allocated(field_in)) allocate(field_in(nCells,nVertLevels))
                if(.not.allocated(press_in)) allocate(press_in(nCells,nVertLevels))
                if(.not.allocated(field_interp)) allocate(field_interp(nCells, nIsoLevels))
                !reverse the vertical axis of array
                do iCell=1,nCells
                    do k=1,nVertLevels
                        kk = nVertLevels+1-k
                        press_in(iCell,kk) = pressure(k,iCell) * 100.
                        field_in(iCell,kk) = dtheta_dt_mp(k,iCell)
                    end do
                end do

                call interp_tofixed_pressure(nCells, nVertLevels, nIsoLevels, press_in, field_in, press_interp, field_interp)

                do k=1,nIsoLevels
                    dtheta_dt_mp_isobaric(k,1:nCells) = field_interp(1:nCells,k)
                end do
                if(allocated(press_in)) deallocate(press_in)
                if(allocated(field_in)) deallocate(field_in)
                if(allocated(field_interp)) deallocate(field_interp)
                call mpas_log_write('Did isobaric interp of microphysics heating tendency.')
            end if

        end if
       
        if(allocated(press_interp)) deallocate(press_interp)
        if(allocated(pressureCp1)) deallocate(pressureCp1)
        call mpas_log_write('Finished isobaric interpolation of PV vars.')
        
    end subroutine interp_pv_diag_isobaric

   !==================================================================================================
    subroutine interp_tofixed_pressure(ncol,nlev_in,nlev_out,pres_in,field_in,pres_out,field_out)
   !==================================================================================================
   
   !input arguments:
    integer,intent(in):: ncol,nlev_in,nlev_out
   
    real(kind=RKIND),intent(in),dimension(ncol,nlev_in) :: pres_in,field_in
    real(kind=RKIND),intent(in),dimension(ncol,nlev_out):: pres_out
   
   !output arguments:
    real(kind=RKIND),intent(out),dimension(ncol,nlev_out):: field_out
   
   !local variables:
   ! integer:: i1,i2,icol,k,kk
    integer:: icol,k,kk
    integer:: kkstart,kount
    integer,dimension(ncol):: kupper
   
    real(kind=RKIND):: dpl,dpu
   
   !--------------------------------------------------------------------------------------------------
   
   !call mpas_log_write('')
   !call mpas_log_write('--- enter subroutine interp_tofixed_pressure:')
   !call mpas_log_write('... ncol     = $i',intArgs=(/ncol/))
   !call mpas_log_write('... nlev_in  = $i',intArgs=(/nlev_in/))
   !call mpas_log_write('... nlev_out = $i',intArgs=(/nlev_out/))
   !i1=1 ; i2=ncol
   !do k = 1, nlev_in
   !   call mpas_log_write('$i $r $r $r $r', intArgs=(/k/), realArgs=(/pres_in(i1,k),field_in(i1,k),pres_in(i2,k),field_in(i2,k)/))
   !enddo
   !call mpas_log_write('')
   
    do icol = 1, ncol
       kupper(icol) = 1
    enddo
   
    do k = 1, nlev_out 
   
       kkstart = nlev_in
       do icol = 1, ncol
          kkstart = min0(kkstart,kupper(icol))
       enddo
       kount = 0
   
       do kk = kkstart, nlev_in-1
          do icol = 1, ncol
             if(pres_out(icol,k).gt.pres_in(icol,kk).and.pres_out(icol,k).le.pres_in(icol,kk+1)) then
                kupper(icol) = kk
                kount = kount + 1
   !            call mpas_log_write('$i $i $r $r $r', intArgs=(/k,kupper(icol)/), realArgs=(/pres_out(icol,k),pres_in(icol,kk),pres_in(icol,kk+1)/))
             endif
          enddo
   
          if(kount.eq.ncol) then
             do icol = 1, ncol
                dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
                dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
                field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                                  + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
             end do
             goto 35
           end if
       enddo
   
       do icol = 1, ncol
          if(pres_out(icol,k) .lt. pres_in(icol,1)) then
             field_out(icol,k) = field_in(icol,1)*pres_out(icol,k)/pres_in(icol,1)
          elseif(pres_out(icol,k) .gt. pres_in(icol,nlev_in)) then
             field_out(icol,k) = field_in(icol,nlev_in)
          else
             dpu = pres_out(icol,k) - pres_in(icol,kupper(icol))
             dpl = pres_in(icol,kupper(icol)+1) - pres_out(icol,k)
             field_out(icol,k) = (field_in(icol,kupper(icol))*dpl &
                               + field_in(icol,kupper(icol)+1)*dpu)/(dpl + dpu)
          endif
       enddo
   
    35 continue
   !   call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/kupper(i1)/), &
   !                       realArgs=(/pres_out(i1,k),pres_in(i1,kupper(i1)),pres_in(i1,kupper(i1)+1),field_out(i1,k),field_in(i1,kupper(i1)),field_in(i1,kupper(i1)+1)/))
   !   call mpas_log_write('$i $r $r $r $r $r $r', intArgs=(/kupper(i2)/), &
   !                       realArgs=(/pres_out(i2,k),pres_in(i2,kupper(i2)),pres_in(i2,kupper(i2)+1),field_out(i2,k),field_in(i2,kupper(i2)),field_in(i2,kupper(i2)+1)/))
   
    enddo
   
    end subroutine interp_tofixed_pressure
   

   !==================================================================================================
    subroutine compute_slp(ncol,nlev_in,nscalars,t,height,p,index_qv,scalars,slp)
   !==================================================================================================
    
       implicit none
   
      !input arguments:
       integer, intent(in) :: ncol, nlev_in, nscalars
      
      !p: in mb
      !t: in K
      !scalars: in kg/kg
      !height: in m
       real(kind=RKIND), intent(in), dimension(nlev_in,ncol) :: p,t
       real(kind=RKIND), intent(in), dimension(nlev_in+1,ncol) :: height
       integer, intent(in) :: index_qv
       real(kind=RKIND), intent(in), dimension(nscalars,nlev_in,ncol) :: scalars
      
      !output arguments:
       real(kind=RKIND), intent(out), dimension(ncol) :: slp
      
      !local variables:
       integer :: icol, k, kcount
       integer :: klo, khi
      
       real(kind=RKIND) :: gamma, rr, grav
       parameter (rr=287.0, grav=9.80616, gamma=0.0065)
      
       real(kind=RKIND) :: tc, pconst
       parameter (tc=273.16+17.5, pconst=100.)
      
       logical mm5_test
       parameter (mm5_test=.true.)
      
       integer, dimension(:), allocatable :: level
       real(kind=RKIND), dimension(:), allocatable :: t_surf, t_msl
       real(kind=RKIND) :: plo , phi , tlo, thi , zlo , zhi
       real(kind=RKIND) :: p_at_pconst , t_at_pconst , z_at_pconst, z_half_lowest
      
       logical :: l1, l2, l3, found
      
      ! Find least zeta level that is PCONST Pa above the surface.  We later use this
      ! level to extrapolate a surface pressure and temperature, which is supposed
      ! to reduce the effect of the diurnal heating cycle in the pressure field.
      
       if (.not.allocated(level))  allocate(level(ncol))
       if (.not.allocated(t_surf)) allocate(t_surf(ncol))
       if (.not.allocated(t_msl))  allocate(t_msl(ncol))
      
       do icol = 1 , ncol
          level(icol) = -1
      
          k = 1
          found = .false.
          do while ( (.not. found) .and. (k.le.nlev_in))
                if ( p(k,icol) .lt. p(1,icol)-pconst ) then
                   level(icol) = k
                   found = .true.
                end if
                k = k+1
          end do
      
          if ( level(icol) .eq. -1 ) then
             call mpas_log_write('Troubles finding level $r above ground.', realArgs=(/pconst/))
             call mpas_log_write('Problems first occur at ($i)', intArgs=(/icol/))
             call mpas_log_write('Surface pressure = $r hPa.', realArgs=(/p(1,icol)/))
             call mpas_log_write('*** MSLP field will not be computed')
             slp(:) = 0.0
             return
          end if
      
       end do
      
      ! Get temperature PCONST hPa above surface.  Use this to extrapolate
      ! the temperature at the surface and down to sea level.
      
       do icol = 1 , ncol
      
          klo = max ( level(icol) - 1 , 1      )
          khi = min ( klo + 1        , nlev_in - 1 )
      
          if ( klo .eq. khi ) then
             call mpas_log_write('Trapping levels are weird.')
             call mpas_log_write('icol = $i', intArgs=(/icol/))
             call mpas_log_write('klo = $i, khi = $i: and they should not be equal.', intArgs=(/klo,khi/))
             call mpas_log_write('*** MSLP field will not be computed')
             slp(:) = 0.0
             return
          end if
      
          plo = p(klo,icol)
          phi = p(khi,icol)
          tlo = t(klo,icol) * (1. + 0.608 * scalars(index_qv,klo,icol))
          thi = t(khi,icol) * (1. + 0.608 * scalars(index_qv,khi,icol))
          zlo = 0.5*(height(klo,icol)+height(klo+1,icol))
          zhi = 0.5*(height(khi,icol)+height(khi+1,icol))
      
          p_at_pconst = p(1,icol) - pconst
          t_at_pconst = thi-(thi-tlo)*log(p_at_pconst/phi)*log(plo/phi)
          z_at_pconst = zhi-(zhi-zlo)*log(p_at_pconst/phi)*log(plo/phi)
      
          t_surf(icol) = t_at_pconst*(p(1,icol)/p_at_pconst)**(gamma*rr/grav)
          t_msl(icol) = t_at_pconst+gamma*z_at_pconst
      
       end do
      
      ! If we follow a traditional computation, there is a correction to the sea level
      ! temperature if both the surface and sea level temnperatures are *too* hot.
      
       if ( mm5_test ) then
          kcount = 0
          do icol = 1 , ncol
                l1 = t_msl(icol) .lt. tc
                l2 = t_surf(icol) .le. tc
                l3 = .not. l1
                if ( l2 .and. l3 ) then
                   t_msl(icol) = tc
                else
                   t_msl(icol) = tc - 0.005*(t_surf(icol)-tc)**2
                   kcount = kcount+1
                end if
          end do
      !   call mpas_log_write('These number of points had t_msl adjusted $i', intArgs=(/kcount/))
       end if
      
       do icol = 1 , ncol
          z_half_lowest=0.5*(height(1,icol)+height(2,icol))
          slp(icol) = p(1,icol) * exp((2.*grav*z_half_lowest)/ &
                                    (rr*(t_msl(icol)+t_surf(icol))))
       end do
      
       if (allocated(level))  deallocate(level)
       if (allocated(t_surf)) deallocate(t_surf)
       if (allocated(t_msl))  deallocate(t_msl)
   
    end subroutine compute_slp


    !*********************************************************************************************************************
    ! MC: Subroutine to interpolate vertical vorticity to cell centers from the vertical vorticity at vertices
    !********************************************************************************************************************* 
   !==================================================================================================
    subroutine interp_absVertVort(vorticity_vertex, nCells, nEdgesOnCell, verticesOnCell, &
                                cellsOnVertex, areaCell, kiteAreasOnVertex, vorticity_cell)
   !==================================================================================================

       IMPLICIT NONE

       integer, intent(in) :: nCells
       integer, dimension(:), intent(in) :: nEdgesOnCell
       integer, dimension(:,:), intent(in) :: verticesOnCell, cellsOnVertex
       real(kind=RKIND), dimension(:), intent(in) :: areaCell
       real(kind=RKIND), dimension(:,:), intent(in) :: vorticity_vertex, kiteAreasOnVertex
       real(kind=RKIND), dimension(:,:), intent(out) :: vorticity_cell
       integer :: i, j, cellIndOnVertex, iVertex

       vorticity_cell(:,:) = 0.0_RKIND

       do i=1,nCells
          do j=1,nEdgesOnCell(i)
             iVertex = verticesOnCell(j,i)
             cellIndOnVertex = FINDLOC(cellsOnVertex(:,iVertex),VALUE=i,DIM=1)
             vorticity_cell(:,i) = vorticity_cell(:,i) + kiteAreasOnVertex(cellIndOnVertex,iVertex) * vorticity_vertex(:,iVertex)
          end do
          vorticity_cell(:,i) = vorticity_cell(:,i) / areaCell(i)
       end do

    end subroutine interp_absVertVort

   !***********************************************************************
   !
   !  routine compute_layer_mean
   !
   !> \brief   Computes the mean of a field in the specified layer.
   !> \author  Michael Duda
   !> \date    3 July 2014
   !> \details
   !>  Given a 3d pressure field, press_in(nCells,nVertLevels), with pressure 
   !>  increasing with vertical index, and a 3d field, 
   !>  field_in(nCells,nVertLevels) with levels in the same order, this routine
   !>  will compute the mean of the field for each column between pressures
   !>  p1 and p2.
   !
   !----------------------------------------------------------------------- 
    subroutine compute_layer_mean(layerMean, p1, p2, field_in, press_in)
   
       implicit none
   
       real(kind=RKIND), dimension(:), intent(out) :: layerMean
       real(kind=RKIND), intent(in) :: p1, p2
       real(kind=RKIND), dimension(:,:), intent(in) :: field_in
       real(kind=RKIND), dimension(:,:), intent(in) :: press_in
   
       integer :: nCells, nVertLevels
       integer :: iCell, k
       integer :: k_bot, k_top
       real(kind=RKIND) :: p_bot, p_top
       real(kind=RKIND) :: wtop_p, wtop_m
       real(kind=RKIND) :: wbot_p, wbot_m
       real(kind=RKIND) :: wtotal, w
       real(kind=RKIND) :: temp
   
   
       !
       ! Get dimensions of input arrays
       !
       nCells = size(field_in, 1) 
       nVertLevels = size(field_in, 2) 
   
   
       !
       ! Check that pressure is increasing with index
       !
       if (press_in(1,1) > press_in(1,nVertLevels)) then
           call mpas_log_write('Error in compute_layer_mean: pressure should increase with index', messageType=MPAS_LOG_ERR)
           layerMean(:) = 0.0
           return
       end if
   
       
       !
       ! Set the pressure at the top and bottom of the layer
       !
       if (p1 < p2) then
          p_top = p1
          p_bot = p2
       else
          p_top = p2
          p_bot = p1
       end if
   
   
       !
       ! For each column, compute the mean value of the field between p_bot and
       ! p_top, with the field weighted by delta-p in each layer
       !
       do iCell=1,nCells
          k_bot = -1
          k_top = -1
   
          ! Search for trapping levels: k_top is the index just above (or equal to)
          ! p_top, and k_bot is the index just above (or equal to) p_bot.
          do k=1,nVertLevels-1
             if (press_in(iCell,k) <= p_top .and. press_in(iCell,k+1) > p_top) then
                k_top = k
                wtop_p = (p_top - press_in(iCell,k)) / (press_in(iCell,k+1) - press_in(iCell,k))
                wtop_m = (press_in(iCell,k+1) - p_top) / (press_in(iCell,k+1) - press_in(iCell,k))
             end if
             if (press_in(iCell,k) <= p_bot .and. press_in(iCell,k+1) > p_bot) then
                k_bot = k
                wbot_m = (p_bot - press_in(iCell,k)) / (press_in(iCell,k+1) - press_in(iCell,k))
                wbot_p = (press_in(iCell,k+1) - p_bot) / (press_in(iCell,k+1) - press_in(iCell,k))
             end if
          end do
   
          if (k_top == -1 .or. k_bot == -1) then      ! Layer intersects top or bottom boundary
   
             layerMean(iCell) = 0.0  
   
          else if (k_top == k_bot) then               ! Layer lies entirely within a single model layer
   
             layerMean(iCell) = wtop_m * field_in(iCell,k_top) + wtop_p * field_in(iCell,k_top+1)
             layerMean(iCell) = layerMean(iCell) + wbot_m * field_in(iCell,k_bot) + wbot_p * field_in(iCell,k_bot+1)
             layerMean(iCell) = 0.5 * layerMean(iCell)
   
          else
   
             ! First layer: from p_top down to press_in(iCell,k_top+1)
             wtotal = press_in(iCell,k_top+1) - p_top  
             temp = wtop_m * field_in(iCell,k_top) + wtop_p * field_in(iCell,k_top+1)
             layerMean(iCell) = wtotal * 0.5 * (field_in(iCell,k_top+1) + temp)
   
             ! Middle layers
             do k=k_top+1,k_bot-1
                w = press_in(iCell,k+1) - press_in(iCell,k)
                wtotal = wtotal + w
                layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell,k) + field_in(iCell,k+1))
             end do
   
             ! Last layer: from press_in(iCell,k_bot) down to p_bot
             w = p_bot - press_in(iCell,k_bot)
             wtotal = wtotal + w
             temp = wbot_m * field_in(iCell,k_bot) + wbot_p * field_in(iCell,k_bot+1)
             layerMean(iCell) = layerMean(iCell) + w * 0.5 * (field_in(iCell,k_bot) + temp)
   
             layerMean(iCell) = layerMean(iCell) / wtotal
          end if
   
       end do
   
    end subroutine compute_layer_mean
   
end module isobaric_diagnostics
