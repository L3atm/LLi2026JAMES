module microp_aero

!---------------------------------------------------------------------------------
! Purpose:
!   CAM driver layer for aerosol activation processes.
!
! ***N.B.*** This module is currently hardcoded to recognize only the aerosols/modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!
! Author: Andrew Gettelman
! Based on code from: Hugh Morrison, Xiaohong Liu and Steve Ghan
! May 2010
! Description in: Morrison and Gettelman, 2008. J. Climate (MG2008)
!                 Gettelman et al., 2010 J. Geophys. Res. - Atmospheres (G2010)         
! for questions contact Andrew Gettelman  (andrew@ucar.edu)
! Modifications: A. Gettelman Nov 2010  - changed to support separation of 
!                  microphysics and macrophysics and concentrate aerosol information here
!                B. Eaton, Sep 2014 - Refactored to move CAM interface code into the CAM
!                  interface modules and preserve just the driver layer functionality here.
!
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8=>shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use ref_pres,         only: top_lev => trop_cloud_top_lev
use physconst,        only: rair
use constituents,     only: cnst_get_ind
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init, physics_ptend_sum, &
                            physics_state_copy, physics_update
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_old_tim_idx, pbuf_get_field
use phys_control,     only: phys_getopts, use_hetfrz_classnuc
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_aer_mmr, rad_cnst_get_aer_props, &
                            rad_cnst_get_mode_num

use nucleate_ice_cam, only: use_preexisting_ice, nucleate_ice_cam_readnl, nucleate_ice_cam_register, &
                            nucleate_ice_cam_init, nucleate_ice_cam_calc

use ndrop,            only: ndrop_init, dropmixnuc
use ndrop_bam,        only: ndrop_bam_init, ndrop_bam_run, ndrop_bam_ccn

use hetfrz_classnuc_cam, only: hetfrz_classnuc_cam_readnl, hetfrz_classnuc_cam_register, hetfrz_classnuc_cam_init, &
                               hetfrz_classnuc_cam_save_cbaero, hetfrz_classnuc_cam_calc

use cam_history,      only: addfld, add_default, outfld
use cam_logfile,      only: iulog
use cam_abortutils,   only: endrun

implicit none
private
save

public :: microp_aero_init, microp_aero_run, microp_aero_readnl, microp_aero_register

! Private module data

character(len=16)   :: eddy_scheme

! contact freezing due to dust
! dust number mean radius (m), Zender et al JGR 2003 assuming number mode radius of 0.6 micron, sigma=2
real(r8), parameter :: rn_dst1 = 0.258e-6_r8
real(r8), parameter :: rn_dst2 = 0.717e-6_r8
real(r8), parameter :: rn_dst3 = 1.576e-6_r8
real(r8), parameter :: rn_dst4 = 3.026e-6_r8

real(r8) :: bulk_scale    ! prescribed aerosol bulk sulfur scale factor

! smallest mixing ratio considered in microphysics
real(r8), parameter :: qsmall = 1.e-18_r8

! minimum allowed cloud fraction
real(r8), parameter :: mincld = 0.0001_r8

! indices in state%q and pbuf structures
integer :: cldliq_idx = -1
integer :: cldice_idx = -1
integer :: numliq_idx = -1
integer :: numice_idx = -1
integer :: kvh_idx = -1
integer :: tke_idx = -1
integer :: wp2_idx = -1
integer :: ast_idx = -1
integer :: cldo_idx = -1
integer :: dgnumwet_idx = -1

! Bulk aerosols
character(len=20), allocatable :: aername(:)
real(r8), allocatable :: num_to_mass_aer(:)

integer :: naer_all      ! number of aerosols affecting climate
integer :: idxsul   = -1 ! index in aerosol list for sulfate
integer :: idxdst2  = -1 ! index in aerosol list for dust2
integer :: idxdst3  = -1 ! index in aerosol list for dust3
integer :: idxdst4  = -1 ! index in aerosol list for dust4

! modal aerosols
logical :: clim_modal_aero

integer :: mode_accum_idx  = -1  ! index of accumulation mode
integer :: mode_aitken_idx = -1  ! index of aitken mode
#if ( defined MODAL_AERO_10MODE )
integer :: mode_coarse6_idx = -1  ! index of coarse mode
integer :: mode_coarse7_idx = -1  ! index of coarse mode
integer :: mode_coarse8_idx = -1  ! index of coarse mode
integer :: mode_coarse9_idx = -1  ! index of coarse mode
integer :: mode_coarse10_idx = -1  ! index of coarse mode
!#if ( defined MODAL_AERO_4MODE )
integer :: mode_coarse_idx = -1  ! index of coarse mode
#endif
integer :: mode_coarse_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse7_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse8_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse9_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse10_dst_idx = -1  ! index of coarse dust mode
integer :: mode_coarse_slt_idx = -1  ! index of coarse sea salt mode
#if ( defined DUSTTYPES )
#if ( defined MODAL_AERO_4MODE )
integer :: coarse_dust1_idx = -1
integer :: coarse_dust2_idx = -1
integer :: coarse_dust3_idx = -1
integer :: coarse_dust4_idx = -1
integer :: coarse_dust5_idx = -1
integer :: coarse_dust6_idx = -1
integer :: coarse_dust7_idx = -1
integer :: coarse_dust8_idx = -1
#elif ( defined MODAL_AERO_10MODE )
! ++ Longlei Li MAM10 ++ 
integer :: coarse7_dust1_idx = -1
integer :: coarse7_dust2_idx = -1
integer :: coarse7_dust3_idx = -1
integer :: coarse7_dust4_idx = -1
integer :: coarse7_dust5_idx = -1
integer :: coarse7_dust6_idx = -1
integer :: coarse7_dust7_idx = -1
integer :: coarse7_dust8_idx = -1

integer :: coarse8_dust1_idx = -1
integer :: coarse8_dust2_idx = -1
integer :: coarse8_dust3_idx = -1
integer :: coarse8_dust4_idx = -1
integer :: coarse8_dust5_idx = -1
integer :: coarse8_dust6_idx = -1
integer :: coarse8_dust7_idx = -1
integer :: coarse8_dust8_idx = -1

integer :: coarse9_dust1_idx = -1
integer :: coarse9_dust2_idx = -1
integer :: coarse9_dust3_idx = -1
integer :: coarse9_dust4_idx = -1
integer :: coarse9_dust5_idx = -1
integer :: coarse9_dust6_idx = -1
integer :: coarse9_dust7_idx = -1
integer :: coarse9_dust8_idx = -1

integer :: coarse10_dust1_idx = -1
integer :: coarse10_dust2_idx = -1
integer :: coarse10_dust3_idx = -1
integer :: coarse10_dust4_idx = -1
integer :: coarse10_dust5_idx = -1
integer :: coarse10_dust6_idx = -1
integer :: coarse10_dust7_idx = -1
integer :: coarse10_dust8_idx = -1

! -- Longlei Li MAM10 -- 
#endif
#else
integer :: coarse_dust_idx = -1  ! index of dust in coarse mode
#endif
integer :: coarse_nacl_idx = -1  ! index of nacl in coarse mode
#if (defined MODAL_AERO_10MODE)
integer :: coarse6_so4_idx = -1  ! index of sulfate in coarse mode
integer :: coarse7_so4_idx = -1  ! index of sulfate in coarse mode
integer :: coarse8_so4_idx = -1  ! index of sulfate in coarse mode
integer :: coarse9_so4_idx = -1  ! index of sulfate in coarse mode
integer :: coarse10_so4_idx = -1  ! index of sulfate in coarse mode
!#elif (defined MODAL_AERO_4MODE)
integer :: coarse_so4_idx = -1  ! index of sulfate in coarse mode
#endif

integer :: npccn_idx, rndst_idx, nacon_idx

!logical  :: separate_dust = .false.
! Longlei Li
logical  :: separate_dust = .false.

!=========================================================================================
contains
!=========================================================================================

subroutine microp_aero_register
   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Register pbuf fields for aerosols needed by microphysics
   ! 
   ! Author: Cheryl Craig October 2012
   ! 
   !-----------------------------------------------------------------------
   use ppgrid,         only: pcols
   use physics_buffer, only: pbuf_add_field, dtype_r8

   call pbuf_add_field('NPCCN',      'physpkg',dtype_r8,(/pcols,pver/), npccn_idx)

   call pbuf_add_field('RNDST',      'physpkg',dtype_r8,(/pcols,pver,4/), rndst_idx)
   call pbuf_add_field('NACON',      'physpkg',dtype_r8,(/pcols,pver,4/), nacon_idx)
 
   call nucleate_ice_cam_register()
   call hetfrz_classnuc_cam_register()

end subroutine microp_aero_register

!=========================================================================================

subroutine microp_aero_init

   !----------------------------------------------------------------------- 
   ! 
   ! Purpose: 
   ! Initialize constants for aerosols needed by microphysics
   ! 
   ! Author: Andrew Gettelman May 2010
   ! 
   !-----------------------------------------------------------------------

   use cam_logfile,       only: iulog

   ! local variables
   integer  :: iaer, ierr
   integer  :: m, n, nmodes, nspec

   character(len=32) :: str32
   character(len=*), parameter :: routine = 'microp_aero_init'
   logical :: history_amwg
   !-----------------------------------------------------------------------

   ! Query the PBL eddy scheme
   call phys_getopts(eddy_scheme_out          = eddy_scheme,  &
                     history_amwg_out         = history_amwg )

   ! Access the physical properties of the aerosols that are affecting the climate
   ! by using routines from the rad_constituents module.

   ! get indices into state and pbuf structures
   call cnst_get_ind('CLDLIQ', cldliq_idx)
   call cnst_get_ind('CLDICE', cldice_idx)
   call cnst_get_ind('NUMLIQ', numliq_idx)
   call cnst_get_ind('NUMICE', numice_idx)
   select case(trim(eddy_scheme))
   case ('diag_TKE')
      tke_idx      = pbuf_get_index('tke')   
   case ('CLUBB_SGS')
      wp2_idx = pbuf_get_index('WP2_nadv')
   case default
      kvh_idx      = pbuf_get_index('kvh')
   end select

   ! clim_modal_aero determines whether modal aerosols are used in the climate calculation.
   ! The modal aerosols can be either prognostic or prescribed.
   call rad_cnst_get_info(0, nmodes=nmodes)
   clim_modal_aero = (nmodes > 0)

   ast_idx = pbuf_get_index('AST')

   if (clim_modal_aero) then

      cldo_idx     = pbuf_get_index('CLDO')
      dgnumwet_idx = pbuf_get_index('DGNUMWET')

      call ndrop_init()

      ! Init indices for specific modes/species

      ! mode index for specified mode types
      do m = 1, nmodes
         call rad_cnst_get_info(0, m, mode_type=str32)
         select case (trim(str32))
         case ('accum')
            mode_accum_idx = m
         case ('aitken')
            mode_aitken_idx = m
        ! Longlei Li commented out for MAM9
        ! case ('coarse')
        !    mode_coarse_idx = m
        ! case ('coarse_dust')
        !    mode_coarse_dst_idx = m
         case ('coarse_seasalt')
            mode_coarse_slt_idx = m
         case ('coarse7_dust')
            mode_coarse7_dst_idx = m
            mode_coarse_idx = m
            mode_coarse_dst_idx = m
         case ('coarse8_dust')
            mode_coarse8_dst_idx = m
         case ('coarse9_dust')
            mode_coarse9_dst_idx = m
         case ('coarse10_dust')
            mode_coarse10_dst_idx = m
         end select
      end do

      !! check if coarse dust is in separate mode
      !separate_dust = mode_coarse_dst_idx > 0

      !! for 3-mode 
      !if ( mode_coarse_dst_idx<0 ) mode_coarse_dst_idx = mode_coarse_idx
      !if ( mode_coarse_slt_idx<0 ) mode_coarse_slt_idx = mode_coarse_idx

    !write(iulog,*)'  Longlei debugging (microp_aero.F90: test01)'

      ! Check that required mode types were found
      if (mode_accum_idx == -1 .or. mode_aitken_idx == -1 .or. &
          mode_coarse_slt_idx == -1 .or. mode_coarse9_dst_idx == -1) then
         write(iulog,*) routine//': ERROR required mode type not found - mode idx:', &
            mode_accum_idx, mode_aitken_idx, mode_coarse_dst_idx, mode_coarse_slt_idx
            !mode_accum_idx, mode_aitken_idx, mode_coarse8_dst_idx == -1, mode_coarse9_dst_idx == -1, mode_coarse_slt_idx
         call endrun(routine//': ERROR required mode type not found')
      end if

    !write(iulog,*)'  Longlei debugging (microp_aero.F90: test02)'


      ! species indices for specified types
      ! find indices for the dust and seasalt species in the coarse mode
      call rad_cnst_get_info(0, mode_coarse_dst_idx, nspec=nspec)
      ! ++ Longlei Li MAM10 ++
      call rad_cnst_get_info(0, mode_coarse7_dst_idx, nspec=nspec)
      call rad_cnst_get_info(0, mode_coarse8_dst_idx, nspec=nspec)
      call rad_cnst_get_info(0, mode_coarse9_dst_idx, nspec=nspec)
      call rad_cnst_get_info(0, mode_coarse10_dst_idx, nspec=nspec)

      do n = 1, nspec
         !call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse7_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse8_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse9_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse10_dst_idx, n, spec_type=str32)

         select case (trim(str32))
#if ( defined DUSTTYPES )
#if ( defined MODAL_AERO_4MODE )
         case ('dust-type1')
            coarse_dust1_idx = n
         case ('dust-type2')
            coarse_dust2_idx = n
         case ('dust-type3')
            coarse_dust3_idx = n
         case ('dust-type4')
            coarse_dust4_idx = n
         case ('dust-type5')
            coarse_dust5_idx = n
         case ('dust-type6')
            coarse_dust6_idx = n
         case ('dust-type7')
            coarse_dust7_idx = n
         case ('dust-type8')
            coarse_dust8_idx = n
#elif ( defined MODAL_AERO_10MODE )
         case ('dust-type1')
            coarse7_dust1_idx = n
            coarse8_dust1_idx = n
            coarse9_dust1_idx = n
            coarse10_dust1_idx = n
         case ('dust-type2')
            coarse7_dust2_idx = n
            coarse8_dust2_idx = n
            coarse9_dust2_idx = n
            coarse10_dust2_idx = n
         case ('dust-type3')
            coarse7_dust3_idx = n
            coarse8_dust3_idx = n
            coarse9_dust3_idx = n
            coarse10_dust3_idx = n
         case ('dust-type4')
            coarse7_dust4_idx = n
            coarse8_dust4_idx = n
            coarse9_dust4_idx = n
            coarse10_dust4_idx = n
         case ('dust-type5')
            coarse7_dust5_idx = n
            coarse8_dust5_idx = n
            coarse9_dust5_idx = n
            coarse10_dust5_idx = n
         case ('dust-type6')
            coarse7_dust6_idx = n
            coarse8_dust6_idx = n
            coarse9_dust6_idx = n
            coarse10_dust6_idx = n
         case ('dust-type7')
            coarse7_dust7_idx = n
            coarse8_dust7_idx = n
            coarse9_dust7_idx = n
            coarse10_dust7_idx = n
         case ('dust-type8')
            coarse7_dust8_idx = n
            coarse8_dust8_idx = n
            coarse9_dust8_idx = n
            coarse10_dust8_idx = n
#endif
!#else
!         case ('dust')
!            coarse_dust_idx = n
#endif
         end select
      end do

      do n = 1, nspec
         !call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse7_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust-type1')
            coarse7_dust1_idx = n
         case ('dust-type2')
            coarse7_dust2_idx = n
         case ('dust-type3')
            coarse7_dust3_idx = n
         case ('dust-type4')
            coarse7_dust4_idx = n
         case ('dust-type5')
            coarse7_dust5_idx = n
         case ('dust-type6')
            coarse7_dust6_idx = n
         case ('dust-type7')
            coarse7_dust7_idx = n
         case ('dust-type8')
            coarse7_dust8_idx = n
         end select
      end do

      do n = 1, nspec
         !call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse8_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust-type1')
            coarse8_dust1_idx = n
         case ('dust-type2')
            coarse8_dust2_idx = n
         case ('dust-type3')
            coarse8_dust3_idx = n
         case ('dust-type4')
            coarse8_dust4_idx = n
         case ('dust-type5')
            coarse8_dust5_idx = n
         case ('dust-type6')
            coarse8_dust6_idx = n
         case ('dust-type7')
            coarse8_dust7_idx = n
         case ('dust-type8')
            coarse8_dust8_idx = n
         end select
      end do

      do n = 1, nspec
         !call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse9_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust-type1')
            coarse9_dust1_idx = n
         case ('dust-type2')
            coarse9_dust2_idx = n
         case ('dust-type3')
            coarse9_dust3_idx = n
         case ('dust-type4')
            coarse9_dust4_idx = n
         case ('dust-type5')
            coarse9_dust5_idx = n
         case ('dust-type6')
            coarse9_dust6_idx = n
         case ('dust-type7')
            coarse9_dust7_idx = n
         case ('dust-type8')
            coarse9_dust8_idx = n
         end select
      end do

      do n = 1, nspec
         !call rad_cnst_get_info(0, mode_coarse_dst_idx, n, spec_type=str32)
         call rad_cnst_get_info(0, mode_coarse10_dst_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('dust-type1')
            coarse10_dust1_idx = n
         case ('dust-type2')
            coarse10_dust2_idx = n
         case ('dust-type3')
            coarse10_dust3_idx = n
         case ('dust-type4')
            coarse10_dust4_idx = n
         case ('dust-type5')
            coarse10_dust5_idx = n
         case ('dust-type6')
            coarse10_dust6_idx = n
         case ('dust-type7')
            coarse10_dust7_idx = n
         case ('dust-type8')
            coarse10_dust8_idx = n
         end select
      end do

      call rad_cnst_get_info(0, mode_coarse_slt_idx, nspec=nspec)
      do n = 1, nspec
         call rad_cnst_get_info(0, mode_coarse_slt_idx, n, spec_type=str32)
         select case (trim(str32))
         case ('seasalt')
            coarse_nacl_idx = n
         end select
      end do

#if ( defined MODAL_AERO_10MODE )
      if (mode_coarse7_dst_idx>0) then
         call rad_cnst_get_info(0, mode_coarse7_dst_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse7_dst_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse7_so4_idx = n
            end select
         end do
      endif

      if (mode_coarse8_dst_idx>0) then
         call rad_cnst_get_info(0, mode_coarse8_dst_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse8_dst_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse8_so4_idx = n
            end select
         end do
      endif

      if (mode_coarse9_dst_idx>0) then
         call rad_cnst_get_info(0, mode_coarse9_dst_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse9_dst_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse9_so4_idx = n
            end select
         end do
      endif

      if (mode_coarse10_dst_idx>0) then
         call rad_cnst_get_info(0, mode_coarse10_dst_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse10_dst_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse10_so4_idx = n
            end select
         end do
      endif

#elif
      if (mode_coarse_idx>0) then
         call rad_cnst_get_info(0, mode_coarse_idx, nspec=nspec)
         do n = 1, nspec
            call rad_cnst_get_info(0, mode_coarse_idx, n, spec_type=str32)
            select case (trim(str32))
            case ('sulfate')
               coarse_so4_idx = n
            end select
         end do
      endif
#endif

      ! Check that required mode specie types were found
#if ( defined DUSTTYPES )
#if ( defined MODAL_AERO_4MODE )
      if ( coarse_dust1_idx == -1 .or. coarse_dust2_idx == -1 .or. &
           coarse_dust3_idx == -1 .or. coarse_dust4_idx == -1 .or. &
           coarse_dust5_idx == -1 .or. coarse_dust6_idx == -1 .or. &
           coarse_dust7_idx == -1 .or. coarse_dust8_idx == -1 .or. &
           coarse_nacl_idx == -1) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust1_idx, coarse_dust2_idx, coarse_dust3_idx, coarse_dust4_idx, &
            coarse_dust5_idx, coarse_dust6_idx, coarse_dust7_idx, coarse_dust8_idx, &
            coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
#elif ( defined MODAL_AERO_10MODE )
      if ( coarse7_dust1_idx == -1 .or. coarse7_dust2_idx == -1 .or. &
           coarse7_dust3_idx == -1 .or. coarse7_dust4_idx == -1 .or. &
           coarse7_dust5_idx == -1 .or. coarse7_dust6_idx == -1 .or. &
           coarse7_dust7_idx == -1 .or. coarse7_dust8_idx == -1 .or. &
           coarse_nacl_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse7_dust1_idx, coarse7_dust2_idx, coarse7_dust3_idx, coarse7_dust4_idx, &
            coarse7_dust5_idx, coarse7_dust6_idx, coarse7_dust7_idx, coarse7_dust8_idx, &
            coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      if ( coarse8_dust1_idx == -1 .or. coarse8_dust2_idx == -1 .or. &
           coarse8_dust3_idx == -1 .or. coarse8_dust4_idx == -1 .or. &
           coarse8_dust5_idx == -1 .or. coarse8_dust6_idx == -1 .or. &
           coarse8_dust7_idx == -1 .or. coarse8_dust8_idx == -1 .or. &
           coarse_nacl_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse8_dust1_idx, coarse8_dust2_idx, coarse8_dust3_idx, coarse8_dust4_idx, &
            coarse8_dust5_idx, coarse8_dust6_idx, coarse8_dust7_idx, coarse8_dust8_idx, &
            coarse_nacl_idx 
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      if ( coarse9_dust1_idx == -1 .or. coarse9_dust2_idx == -1 .or. &
           coarse9_dust3_idx == -1 .or. coarse9_dust4_idx == -1 .or. &
           coarse9_dust5_idx == -1 .or. coarse9_dust6_idx == -1 .or. &
           coarse9_dust7_idx == -1 .or. coarse9_dust8_idx == -1 .or. &
           coarse_nacl_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse9_dust1_idx, coarse9_dust2_idx, coarse9_dust3_idx, coarse9_dust4_idx, &
            coarse9_dust5_idx, coarse9_dust6_idx, coarse9_dust7_idx, coarse9_dust8_idx, &
            coarse_nacl_idx 
         call endrun(routine//': ERROR required mode-species type not found')
      end if

      if ( coarse10_dust1_idx == -1 .or. coarse10_dust2_idx == -1 .or. &
           coarse10_dust3_idx == -1 .or. coarse10_dust4_idx == -1 .or. &
           coarse10_dust5_idx == -1 .or. coarse10_dust6_idx == -1 .or. &
           coarse10_dust7_idx == -1 .or. coarse10_dust8_idx == -1 .or. &
           coarse_nacl_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse10_dust1_idx, coarse10_dust2_idx, coarse10_dust3_idx, coarse10_dust4_idx, &
            coarse10_dust5_idx, coarse10_dust6_idx, coarse10_dust7_idx, coarse10_dust8_idx, &
            coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if


#endif
#else
      if ( coarse_dust_idx == -1 .or. coarse_nacl_idx == -1 ) then
         write(iulog,*) routine//': ERROR required mode-species type not found - indicies:', &
            coarse_dust_idx, coarse_nacl_idx
         call endrun(routine//': ERROR required mode-species type not found')
      end if
#endif
   else

      ! Props needed for BAM number concentration calcs.

      call rad_cnst_get_info(0, naero=naer_all)
      allocate( &
         aername(naer_all),        &
         num_to_mass_aer(naer_all) )

      do iaer = 1, naer_all
         call rad_cnst_get_aer_props(0, iaer, &
            aername         = aername(iaer), &
            num_to_mass_aer = num_to_mass_aer(iaer) )

         ! Look for sulfate, dust, and soot in this list (Bulk aerosol only)
         if (trim(aername(iaer)) == 'SULFATE') idxsul = iaer
         if (trim(aername(iaer)) == 'DUST2') idxdst2 = iaer
         if (trim(aername(iaer)) == 'DUST3') idxdst3 = iaer
         if (trim(aername(iaer)) == 'DUST4') idxdst4 = iaer
      end do

      call ndrop_bam_init()

   end if

   call addfld('LCLOUD', (/ 'lev' /), 'A', ' ',   'Liquid cloud fraction used in stratus activation')

   call addfld('WSUB',   (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity'                   )
   call addfld('WSUBI',  (/ 'lev' /), 'A', 'm/s', 'Diagnostic sub-grid vertical velocity for ice'           )

   if (history_amwg) then
      call add_default ('WSUB     ', 1, ' ')
   end if

   call nucleate_ice_cam_init(mincld, bulk_scale)
   call hetfrz_classnuc_cam_init(mincld)
end subroutine microp_aero_init

!=========================================================================================

subroutine microp_aero_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   use cam_logfile,       only: iulog
   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Namelist variables
   real(r8) :: microp_aero_bulk_scale = 2._r8  ! prescribed aerosol bulk sulfur scale factor
 
   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'microp_aero_readnl'

   namelist /microp_aero_nl/ microp_aero_bulk_scale
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'microp_aero_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, microp_aero_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variable
   call mpibcast(microp_aero_bulk_scale, 1, mpir8, 0, mpicom)
#endif

   ! set local variables
   bulk_scale = microp_aero_bulk_scale

   call nucleate_ice_cam_readnl(nlfile)
   call hetfrz_classnuc_cam_readnl(nlfile)
end subroutine microp_aero_readnl

!=========================================================================================

subroutine microp_aero_run ( &
   state, ptend_all, deltatin, pbuf)

   use cam_logfile,       only: iulog

   ! input arguments
   type(physics_state),         intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend_all
   real(r8),                    intent(in)    :: deltatin     ! time step (s)
   type(physics_buffer_desc),   pointer       :: pbuf(:)

   ! local workspace
   ! all units mks unless otherwise stated

   integer :: i, k, m
   integer :: itim_old
   integer :: nmodes

   type(physics_state) :: state1                ! Local copy of state variable
   type(physics_ptend) :: ptend_loc

   real(r8), pointer :: ast(:,:)        

   real(r8), pointer :: npccn(:,:)      ! number of CCN (liquid activated)

   real(r8), pointer :: rndst(:,:,:)    ! radius of 4 dust bins for contact freezing
   real(r8), pointer :: nacon(:,:,:)    ! number in 4 dust bins for contact freezing

#if ( defined MODAL_AERO_10MODE )
   real(r8), pointer :: num_coarse7(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: num_coarse8(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: num_coarse9(:,:) ! number m.r. of coarse mode
   real(r8), pointer :: num_coarse10(:,:) ! number m.r. of coarse mode
#else
   real(r8), pointer :: num_coarse(:,:) ! number m.r. of coarse mode
#endif

#if ( defined DUSTTYPES )
   real(r8), pointer :: coarse_dust1(:,:)
   real(r8), pointer :: coarse_dust2(:,:)
   real(r8), pointer :: coarse_dust3(:,:)
   real(r8), pointer :: coarse_dust4(:,:)
   real(r8), pointer :: coarse_dust5(:,:)
   real(r8), pointer :: coarse_dust6(:,:)
   real(r8), pointer :: coarse_dust7(:,:)
   real(r8), pointer :: coarse_dust8(:,:)

   real(r8), pointer :: coarse7_dust1(:,:)
   real(r8), pointer :: coarse7_dust2(:,:)
   real(r8), pointer :: coarse7_dust3(:,:)
   real(r8), pointer :: coarse7_dust4(:,:)
   real(r8), pointer :: coarse7_dust5(:,:)
   real(r8), pointer :: coarse7_dust6(:,:)
   real(r8), pointer :: coarse7_dust7(:,:)
   real(r8), pointer :: coarse7_dust8(:,:)

   real(r8), pointer :: coarse8_dust1(:,:)
   real(r8), pointer :: coarse8_dust2(:,:)
   real(r8), pointer :: coarse8_dust3(:,:)
   real(r8), pointer :: coarse8_dust4(:,:)
   real(r8), pointer :: coarse8_dust5(:,:)
   real(r8), pointer :: coarse8_dust6(:,:)
   real(r8), pointer :: coarse8_dust7(:,:)
   real(r8), pointer :: coarse8_dust8(:,:)

   real(r8), pointer :: coarse9_dust1(:,:)
   real(r8), pointer :: coarse9_dust2(:,:)
   real(r8), pointer :: coarse9_dust3(:,:)
   real(r8), pointer :: coarse9_dust4(:,:)
   real(r8), pointer :: coarse9_dust5(:,:)
   real(r8), pointer :: coarse9_dust6(:,:)
   real(r8), pointer :: coarse9_dust7(:,:)
   real(r8), pointer :: coarse9_dust8(:,:)

   real(r8), pointer :: coarse10_dust1(:,:)
   real(r8), pointer :: coarse10_dust2(:,:)
   real(r8), pointer :: coarse10_dust3(:,:)
   real(r8), pointer :: coarse10_dust4(:,:)
   real(r8), pointer :: coarse10_dust5(:,:)
   real(r8), pointer :: coarse10_dust6(:,:)
   real(r8), pointer :: coarse10_dust7(:,:)
   real(r8), pointer :: coarse10_dust8(:,:)
#else
   real(r8), pointer :: coarse_dust(:,:) ! mass m.r. of coarse dust
#endif
   real(r8), pointer :: coarse_nacl(:,:) ! mass m.r. of coarse nacl
#if ( defined MODAL_AERO_10MODE )
   real(r8), pointer :: coarse7_so4(:,:)  ! mass m.r. of coarse sulfate
   real(r8), pointer :: coarse8_so4(:,:)  ! mass m.r. of coarse sulfate
   real(r8), pointer :: coarse9_so4(:,:)  ! mass m.r. of coarse sulfate
   real(r8), pointer :: coarse10_so4(:,:)  ! mass m.r. of coarse sulfate
!#if ( defined MODAL_AERO_4MODE )
   real(r8), pointer :: coarse_so4(:,:)  ! mass m.r. of coarse sulfate
#endif
   real(r8), pointer :: kvh(:,:)        ! vertical eddy diff coef (m2 s-1)
   real(r8), pointer :: tke(:,:)        ! TKE from the UW PBL scheme (m2 s-2)
   real(r8), pointer :: wp2(:,:)        ! CLUBB vertical velocity variance

   real(r8), pointer :: cldn(:,:)       ! cloud fraction
   real(r8), pointer :: cldo(:,:)       ! old cloud fraction

   real(r8), pointer :: dgnumwet(:,:,:) ! aerosol mode diameter

   real(r8), pointer :: aer_mmr(:,:)    ! aerosol mass mixing ratio

   real(r8) :: rho(pcols,pver)     ! air density (kg m-3)

   real(r8) :: lcldm(pcols,pver)   ! liq cloud fraction

   real(r8) :: lcldn(pcols,pver)   ! fractional coverage of new liquid cloud
   real(r8) :: lcldo(pcols,pver)   ! fractional coverage of old liquid cloud
   real(r8) :: cldliqf(pcols,pver) ! fractional of total cloud that is liquid
   real(r8) :: qcld                ! total cloud water
   real(r8) :: nctend_mixnuc(pcols,pver)
   real(r8) :: dum, dum2           ! temporary dummy variable
   real(r8) :: dmc, ssmc, so4mc    ! variables for modal scheme.
   integer  :: dst_idx, num_idx

   ! bulk aerosol variables
   real(r8), allocatable :: naer2(:,:,:)    ! bulk aerosol number concentration (1/m3)
   real(r8), allocatable :: maerosol(:,:,:) ! bulk aerosol mass conc (kg/m3)

   real(r8) :: wsub(pcols,pver)    ! diagnosed sub-grid vertical velocity st. dev. (m/s)
   real(r8) :: wsubi(pcols,pver)   ! diagnosed sub-grid vertical velocity ice (m/s)
   real(r8) :: nucboas

   real(r8) :: wght

   integer :: lchnk, ncol

   real(r8), allocatable :: factnum(:,:,:) ! activation fraction for aerosol number
   !-------------------------------------------------------------------------------
   call physics_state_copy(state,state1)

   lchnk = state1%lchnk
   ncol  = state1%ncol

   itim_old = pbuf_old_tim_idx()
   call pbuf_get_field(pbuf, ast_idx,      ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/))

   call pbuf_get_field(pbuf, npccn_idx, npccn)

   call pbuf_get_field(pbuf, nacon_idx, nacon)
   call pbuf_get_field(pbuf, rndst_idx, rndst)

   call physics_ptend_init(ptend_all, state%psetcols, 'microp_aero')

   if (clim_modal_aero) then

      itim_old = pbuf_old_tim_idx()
      
      call pbuf_get_field(pbuf, ast_idx,  cldn, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )
      call pbuf_get_field(pbuf, cldo_idx, cldo, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

      call rad_cnst_get_info(0, nmodes=nmodes)
      call pbuf_get_field(pbuf, dgnumwet_idx, dgnumwet, start=(/1,1,1/), kount=(/pcols,pver,nmodes/) )

      allocate(factnum(pcols,pver,nmodes))

   end if

   ! initialize output
   npccn(1:ncol,1:pver)    = 0._r8  

   nacon(1:ncol,1:pver,:)  = 0._r8
   ! set default or fixed dust bins for contact freezing
   rndst(1:ncol,1:pver,1) = rn_dst1
   rndst(1:ncol,1:pver,2) = rn_dst2
   rndst(1:ncol,1:pver,3) = rn_dst3
   rndst(1:ncol,1:pver,4) = rn_dst4

   ! save copy of cloud borne aerosols for use in heterogeneous freezing
   if (use_hetfrz_classnuc) then
      call hetfrz_classnuc_cam_save_cbaero(state1, pbuf)
   end if

   ! initialize time-varying parameters
   do k = top_lev, pver
      do i = 1, ncol
         rho(i,k) = state1%pmid(i,k)/(rair*state1%t(i,k))
      end do
   end do

   if (clim_modal_aero) then
      ! mode number mixing ratios
#if ( defined MODAL_AERO_10MODE )
      call rad_cnst_get_mode_num(0, mode_coarse7_dst_idx, 'a', state1, pbuf, num_coarse7)
      call rad_cnst_get_mode_num(0, mode_coarse8_dst_idx, 'a', state1, pbuf, num_coarse8)
      call rad_cnst_get_mode_num(0, mode_coarse9_dst_idx, 'a', state1, pbuf, num_coarse9)
      call rad_cnst_get_mode_num(0, mode_coarse10_dst_idx, 'a', state1, pbuf, num_coarse10)
#else 
      call rad_cnst_get_mode_num(0, mode_coarse_dst_idx, 'a', state1, pbuf, num_coarse)
#endif

      ! mode specie mass m.r.
#if ( defined DUSTTYPES )
#if ( defined MODAL_AERO_4MODE )
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust1_idx, 'a', state1, pbuf, coarse_dust1)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust2_idx, 'a', state1, pbuf, coarse_dust2)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust3_idx, 'a', state1, pbuf, coarse_dust3)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust4_idx, 'a', state1, pbuf, coarse_dust4)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust5_idx, 'a', state1, pbuf, coarse_dust5)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust6_idx, 'a', state1, pbuf, coarse_dust6)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust7_idx, 'a', state1, pbuf, coarse_dust7)
      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust8_idx, 'a', state1, pbuf, coarse_dust8)
#elif ( defined MODAL_AERO_10MODE )
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust1_idx, 'a', state1, pbuf, coarse7_dust1)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust2_idx, 'a', state1, pbuf, coarse7_dust2)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust3_idx, 'a', state1, pbuf, coarse7_dust3)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust4_idx, 'a', state1, pbuf, coarse7_dust4)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust5_idx, 'a', state1, pbuf, coarse7_dust5)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust6_idx, 'a', state1, pbuf, coarse7_dust6)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust7_idx, 'a', state1, pbuf, coarse7_dust7)
      call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_dust8_idx, 'a', state1, pbuf, coarse7_dust8)

      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust1_idx, 'a', state1, pbuf, coarse8_dust1)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust2_idx, 'a', state1, pbuf, coarse8_dust2)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust3_idx, 'a', state1, pbuf, coarse8_dust3)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust4_idx, 'a', state1, pbuf, coarse8_dust4)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust5_idx, 'a', state1, pbuf, coarse8_dust5)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust6_idx, 'a', state1, pbuf, coarse8_dust6)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust7_idx, 'a', state1, pbuf, coarse8_dust7)
      call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_dust8_idx, 'a', state1, pbuf, coarse8_dust8)

      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust1_idx, 'a', state1, pbuf, coarse9_dust1)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust2_idx, 'a', state1, pbuf, coarse9_dust2)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust3_idx, 'a', state1, pbuf, coarse9_dust3)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust4_idx, 'a', state1, pbuf, coarse9_dust4)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust5_idx, 'a', state1, pbuf, coarse9_dust5)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust6_idx, 'a', state1, pbuf, coarse9_dust6)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust7_idx, 'a', state1, pbuf, coarse9_dust7)
      call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_dust8_idx, 'a', state1, pbuf, coarse9_dust8)

      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust1_idx, 'a', state1, pbuf, coarse10_dust1)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust2_idx, 'a', state1, pbuf, coarse10_dust2)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust3_idx, 'a', state1, pbuf, coarse10_dust3)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust4_idx, 'a', state1, pbuf, coarse10_dust4)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust5_idx, 'a', state1, pbuf, coarse10_dust5)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust6_idx, 'a', state1, pbuf, coarse10_dust6)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust7_idx, 'a', state1, pbuf, coarse10_dust7)
      call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_dust8_idx, 'a', state1, pbuf, coarse10_dust8)
#endif
!#else
!      call rad_cnst_get_aer_mmr(0, mode_coarse_dst_idx, coarse_dust_idx, 'a', state1, pbuf, coarse_dust)
#endif
      call rad_cnst_get_aer_mmr(0, mode_coarse_slt_idx, coarse_nacl_idx, 'a', state1, pbuf, coarse_nacl)
!++Longlei Li comments out
      if (mode_coarse_idx>0) then
         !call rad_cnst_get_aer_mmr(0, mode_coarse_idx, coarse_so4_idx, 'a', state1, pbuf, coarse_so4)
         !call rad_cnst_get_aer_mmr(0, mode_coarse6_idx, coarse6_so4_idx, 'a', state1, pbuf, coarse6_so4)
         call rad_cnst_get_aer_mmr(0, mode_coarse7_dst_idx, coarse7_so4_idx, 'a', state1, pbuf, coarse7_so4)
         call rad_cnst_get_aer_mmr(0, mode_coarse8_dst_idx, coarse8_so4_idx, 'a', state1, pbuf, coarse8_so4)
         call rad_cnst_get_aer_mmr(0, mode_coarse9_dst_idx, coarse9_so4_idx, 'a', state1, pbuf, coarse9_so4)
         call rad_cnst_get_aer_mmr(0, mode_coarse10_dst_idx, coarse10_so4_idx, 'a', state1, pbuf, coarse10_so4)
      endif
!--Longlei Li comments out
   else
      ! init number/mass arrays for bulk aerosols
      allocate( &
         naer2(pcols,pver,naer_all), &
         maerosol(pcols,pver,naer_all))

      do m = 1, naer_all
         call rad_cnst_get_aer_mmr(0, m, state1, pbuf, aer_mmr)
         maerosol(:ncol,:,m) = aer_mmr(:ncol,:)*rho(:ncol,:)
         
         if (m .eq. idxsul) then
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)*bulk_scale
         else
            naer2(:ncol,:,m) = maerosol(:ncol,:,m)*num_to_mass_aer(m)
         end if
      end do
   end if

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! More refined computation of sub-grid vertical velocity 
   ! Set to be zero at the surface by initialization.

   select case (trim(eddy_scheme))
   case ('diag_TKE')
      call pbuf_get_field(pbuf, tke_idx, tke)
   case ('CLUBB_SGS')
      itim_old = pbuf_old_tim_idx()
      call pbuf_get_field(pbuf, wp2_idx, wp2, start=(/1,1,itim_old/),kount=(/pcols,pverp,1/))
      allocate(tke(pcols,pverp))
      tke(:ncol,:) = (3._r8/2._r8)*wp2(:ncol,:)

   case default
      call pbuf_get_field(pbuf, kvh_idx, kvh)
   end select

   ! Set minimum values above top_lev.
   wsub(:ncol,:top_lev-1)  = 0.20_r8
   wsubi(:ncol,:top_lev-1) = 0.001_r8

   do k = top_lev, pver
      do i = 1, ncol

         select case (trim(eddy_scheme))
         case ('diag_TKE', 'CLUBB_SGS')
            wsub(i,k) = sqrt(0.5_r8*(tke(i,k) + tke(i,k+1))*(2._r8/3._r8))
            wsub(i,k) = min(wsub(i,k),10._r8)
         case default 
            ! get sub-grid vertical velocity from diff coef.
            ! following morrison et al. 2005, JAS
            ! assume mixing length of 30 m
            dum = (kvh(i,k) + kvh(i,k+1))/2._r8/30._r8
            ! use maximum sub-grid vertical vel of 10 m/s
            dum = min(dum, 10._r8)
            ! set wsub to value at current vertical level
            wsub(i,k)  = dum
         end select

         wsubi(i,k) = max(0.001_r8, wsub(i,k))
         if (.not. use_preexisting_ice) then
            wsubi(i,k) = min(wsubi(i,k), 0.2_r8)
         endif

         wsub(i,k)  = max(0.20_r8, wsub(i,k))

      end do
   end do

   call outfld('WSUB',   wsub, pcols, lchnk)
   call outfld('WSUBI', wsubi, pcols, lchnk)

   if (trim(eddy_scheme) == 'CLUBB_SGS') deallocate(tke)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !ICE Nucleation

   call nucleate_ice_cam_calc(state1, wsubi, pbuf, deltatin, ptend_loc)

   call physics_ptend_sum(ptend_loc, ptend_all, ncol)

   !write(iulog,*) 'Longlei Li marks begin in microp_aero.F90: test03'

   call physics_update(state1, ptend_loc, deltatin)

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! get liquid cloud fraction, check for minimum

   do k = top_lev, pver
      do i = 1, ncol
         lcldm(i,k) = max(ast(i,k), mincld)
      end do
   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Droplet Activation

   if (clim_modal_aero) then

      ! for modal aerosol

      ! partition cloud fraction into liquid water part
      lcldn = 0._r8
      lcldo = 0._r8
      cldliqf = 0._r8
      do k = top_lev, pver
         do i = 1, ncol
            qcld = state1%q(i,k,cldliq_idx) + state1%q(i,k,cldice_idx)
            if (qcld > qsmall) then
               lcldn(i,k)   = cldn(i,k)*state1%q(i,k,cldliq_idx)/qcld
               lcldo(i,k)   = cldo(i,k)*state1%q(i,k,cldliq_idx)/qcld
               cldliqf(i,k) = state1%q(i,k,cldliq_idx)/qcld
            end if
         end do
      end do

      call outfld('LCLOUD', lcldn, pcols, lchnk)

      ! If not using preexsiting ice, then only use cloudbourne aerosol for the
      ! liquid clouds. This is the same behavior as CAM5.
      if (use_preexisting_ice) then
         call dropmixnuc( &
            state1, ptend_loc, deltatin, pbuf, wsub, &
            cldn, cldo, cldliqf, nctend_mixnuc, factnum)
      else   
         cldliqf = 1._r8
         call dropmixnuc( &
            state1, ptend_loc, deltatin, pbuf, wsub, &
            lcldn, lcldo, cldliqf, nctend_mixnuc, factnum)
      end if

      npccn(:ncol,:) = nctend_mixnuc(:ncol,:)

   else

      ! for bulk aerosol

      ! no tendencies returned from ndrop_bam_run, so just init ptend here
      call physics_ptend_init(ptend_loc, state1%psetcols, 'none')

      do k = top_lev, pver
         do i = 1, ncol

            if (state1%q(i,k,cldliq_idx) >= qsmall) then

               ! get droplet activation rate

               call ndrop_bam_run( &
                  wsub(i,k), state1%t(i,k), rho(i,k), naer2(i,k,:), naer_all, &
                  naer_all, maerosol(i,k,:),  &
                  dum2)
               dum = dum2
            else
               dum = 0._r8
            end if

            npccn(i,k) = (dum*lcldm(i,k) - state1%q(i,k,numliq_idx))/deltatin
         end do
      end do

   end if

   call physics_ptend_sum(ptend_loc, ptend_all, ncol)

   call physics_update(state1, ptend_loc, deltatin)

    !write(iulog,*) 'Longlei marks 08 in microp_aero'
   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   ! Contact freezing  (-40<T<-3 C) (Young, 1974) with hooks into simulated dust
   ! estimate rndst and nanco for 4 dust bins here to pass to MG microphysics

   do k = top_lev, pver
      do i = 1, ncol

         if (state1%t(i,k) < 269.15_r8) then

            if (clim_modal_aero) then

               ! For modal aerosols:
               !  use size '3' for dust coarse mode...
               !  scale by dust fraction in coarse mode         
#if ( defined DUSTTYPES )
#if ( defined MODAL_AERO_10MODE )

               dmc  = coarse7_dust1(i,k)+coarse7_dust2(i,k)+coarse7_dust3(i,k)+coarse7_dust4(i,k)+ &
                      coarse7_dust5(i,k)+coarse7_dust6(i,k)+coarse7_dust7(i,k)+coarse7_dust8(i,k)+ &

                      coarse8_dust1(i,k)+coarse8_dust2(i,k)+coarse8_dust3(i,k)+coarse8_dust4(i,k)+ &
                      coarse8_dust5(i,k)+coarse8_dust6(i,k)+coarse8_dust7(i,k)+coarse8_dust8(i,k)+ &

                      coarse9_dust1(i,k)+coarse9_dust2(i,k)+coarse9_dust3(i,k)+coarse9_dust4(i,k)+ &
                      coarse9_dust5(i,k)+coarse9_dust6(i,k)+coarse9_dust7(i,k)+coarse9_dust8(i,k)+ &

                      coarse10_dust1(i,k)+coarse10_dust2(i,k)+coarse10_dust3(i,k)+coarse10_dust4(i,k)+ &
                      coarse10_dust5(i,k)+coarse10_dust6(i,k)+coarse10_dust7(i,k)+coarse10_dust8(i,k)

#else
               dmc  = coarse_dust1(i,k)+coarse_dust2(i,k)+coarse_dust3(i,k)+coarse_dust4(i,k)+ &
                      coarse_dust5(i,k)+coarse_dust6(i,k)+coarse_dust7(i,k)+coarse_dust8(i,k)
#endif

!#else
!               dmc  = coarse_dust(i,k)
#endif
               ssmc = coarse_nacl(i,k)

               if ( separate_dust ) then
                  ! 7-mode -- has separate dust and seasalt mode types and no need for weighting 
                  wght = 1._r8
               else
!++Longlei comments out
                  !so4mc = coarse_so4(i,k)
                  so4mc = coarse7_so4(i,k)+coarse8_so4(i,k)+coarse9_so4(i,k)+coarse10_so4(i,k)
                  !so4mc = 1.0e-36_r8
!--Longlei comments out
                  ! 3-mode -- needs weighting for dust since dust, seasalt, and sulfate  are combined in the "coarse" mode type
!++Longlei deleted sea salt: ssmc
                  !wght = dmc/(ssmc + dmc + so4mc)
                  wght = dmc/(dmc + so4mc)
!--Longlei deleted so4mc
               endif

               if (dmc > 0.0_r8) then
#if ( defined MODAL_AERO_10MODE )
                  nacon(i,k,3) = wght*(num_coarse7(i,k)+num_coarse8(i,k)+num_coarse9(i,k)+num_coarse10(i,k))*rho(i,k)
#else
                  nacon(i,k,3) = wght*num_coarse(i,k)*rho(i,k)
#endif
               else
                  nacon(i,k,3) = 0._r8
               end if

               !also redefine parameters based on size...

               !rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse_dst_idx)
               ! Longlei Li: just set rndst(i,k,3) to be rn_dst3 now

               rndst(i,k,3) = 0.5_r8*dgnumwet(i,k,mode_coarse8_dst_idx)
               if (rndst(i,k,3) <= 0._r8) then 
                  rndst(i,k,3) = rn_dst3
               end if

               !rndst(i,k,3) = rn_dst3

            else

               !For Bulk Aerosols: set equal to aerosol number for dust for bins 2-4 (bin 1=0)

               if (idxdst2 > 0) then 
                  nacon(i,k,2) = naer2(i,k,idxdst2)
               end if
               if (idxdst3 > 0) then 
                  nacon(i,k,3) = naer2(i,k,idxdst3)
               end if
               if (idxdst4 > 0) then 
                  nacon(i,k,4) = naer2(i,k,idxdst4)
               end if
            end if

         end if
      end do

   end do

   !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   !bulk aerosol ccn concentration (modal does it in ndrop, from dropmixnuc)

   if (.not. clim_modal_aero) then

      ! ccn concentration as diagnostic
      call ndrop_bam_ccn(lchnk, ncol, maerosol, naer2)

      deallocate( &
         naer2,    &
         maerosol)

   end if

   ! heterogeneous freezing
   if (use_hetfrz_classnuc) then

      call hetfrz_classnuc_cam_calc(state1, deltatin, factnum, pbuf)

   end if

   if (clim_modal_aero) then
      deallocate(factnum)
   end if
    !write(iulog,*) 'Longlei marks end in microp_aero'

end subroutine microp_aero_run

!=========================================================================================

end module microp_aero
