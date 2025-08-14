!===============================================================================
! Dust for Modal Aerosol Model
!===============================================================================
module dust_model 
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun
  use modal_aero_data,  only: ntot_amode, ndst=>nDust

  implicit none
  private

  public :: dust_names
  public :: dust_nbin
  public :: dust_nnum
#if ( defined DUSTTYPES )
  public :: mine_nbin ! define this parameter here -Longlei Li
#endif
  public :: dust_indices
  public :: dust_emis
  public :: dust_readnl
  public :: dust_init
  public :: dust_active

  integer, protected :: dust_nbin != 2
  integer, protected :: dust_nnum != 2
#if ( defined DUSTTYPES )
  integer, protected :: mine_nbin
#endif

#if  ( defined MODAL_AERO_3MODE )
  character(len=6), protected, allocatable :: dust_names(:)

!++ Longlei Li MAM10 ++
#elif ( defined MODAL_AERO_10MODE )
  character(len=8), protected, allocatable :: dust_names(:)
!-- Longlei Li MAM10 -- 

#endif
  real(r8), allocatable :: dust_dmt_grd(:)
  real(r8), allocatable :: dust_emis_sclfctr(:)
#if ( defined DUSTTYPES )
  real(r8), allocatable :: brit_emis_sclfctr(:)
#endif

  integer , protected, allocatable :: dust_indices(:)
  real(r8), allocatable :: dust_dmt_vwr(:)
  real(r8), allocatable :: dust_stk_crc(:)

  real(r8)          :: dust_emis_fact = -1.e36_r8        ! tuning parameter for dust emissions
  character(len=cl) :: soil_erod_file = 'soil_erod_file' ! full pathname for soil erodibility dataset

  logical :: dust_active = .false.

 contains

  !=============================================================================
  ! reads dust namelist options
  !=============================================================================
  subroutine dust_readnl(nlfile)

    use namelist_utils,  only: find_group_name
    use units,           only: getunit, freeunit
    use mpishorthand
    use cam_logfile,       only: iulog

    character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

    ! Local variables
    integer :: unitn, ierr
    character(len=*), parameter :: subname = 'dust_readnl'

    namelist /dust_nl/ dust_emis_fact, soil_erod_file

    !-----------------------------------------------------------------------------
    ! Read namelist
    if (masterproc) then
       unitn = getunit()
       open( unitn, file=trim(nlfile), status='old' )
       call find_group_name(unitn, 'dust_nl', status=ierr)
       if (ierr == 0) then
          read(unitn, dust_nl, iostat=ierr)
          if (ierr /= 0) then
             call endrun(subname // ':: ERROR reading namelist')
          end if
       end if
       close(unitn)
       call freeunit(unitn)
    end if

#ifdef SPMD
    ! Broadcast namelist variables
    call mpibcast(dust_emis_fact, 1,                   mpir8,   0, mpicom)
    call mpibcast(soil_erod_file, len(soil_erod_file), mpichar, 0, mpicom)
#endif

  end subroutine dust_readnl

  !=============================================================================
  !=============================================================================
  subroutine dust_init()
    use soil_erod_mod, only: soil_erod_init
    use constituents,  only: cnst_get_ind
    use rad_constituents, only: rad_cnst_get_info
    use dust_common,   only: dust_set_params
    use cam_logfile,       only: iulog

    integer :: l, m, mm, ndx, nspec, n
    character(len=32) :: spec_name
    integer, parameter :: mymodes(7) = (/ 2, 1, 3, 4, 5, 6, 7 /) ! tricky order ...

    !dust_nbin = ndst ! not this number any more -Longlei Li
    !dust_nnum = ndst

#if  ( defined MODAL_AERO_3MODE )
    dust_nbin = 2
    dust_nnum = 2
!#elif ( defined MODAL_AERO_4MODE )
!#if ( defined DUSTTYPES )
!    dust_nbin = 24 ! 8Mine * 3mode Longlei Li
!    dust_nnum = 3
!    mine_nbin = 33
!#else
!    dust_nbin = 3
!    dust_nnum = 3
!#endif
!#elif ( defined MODAL_AERO_7MODE )
!    dust_nbin = 2
!    dust_nnum = 2
!++ Longlei Li MAM10 ++ 
#elif ( defined MODAL_AERO_10MODE )
#if ( defined DUSTTYPES )
    dust_nbin = 40 ! 8Mine * 5mode Longlei Li 
    dust_nnum = 5
    mine_nbin = 55 ! 11 * 5mode Longlei Li MAM10
#endif
!-- Longlei Li MAM10 --
#endif

    !allocate( dust_names(2*ndst) )  ! demension changed -Longlei Li
    !allocate( dust_indices(2*ndst) )
    !allocate( dust_dmt_grd(ndst+1) )
    !allocate( dust_emis_sclfctr(ndst) )
    !allocate( dust_dmt_vwr(ndst) )
    !allocate( dust_stk_crc(ndst) )

    allocate( dust_names(dust_nbin+dust_nnum) )
    allocate( dust_indices(dust_nbin+dust_nnum) )
    allocate( dust_dmt_grd(dust_nnum+1) )
    allocate( dust_emis_sclfctr(dust_nbin) ) 
#if ( defined DUSTTYPES )
    allocate( dust_dmt_vwr(dust_nnum) )
    allocate( dust_stk_crc(dust_nnum) )
    allocate( brit_emis_sclfctr(mine_nbin) )
#else
    allocate( dust_dmt_vwr(dust_nbin) )
    allocate( dust_stk_crc(dust_nbin) )
#endif

    if ( ntot_amode == 3 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.011_r8,0.989_r8 /)
    elseif ( ntot_amode == 4 ) then
       dust_dmt_grd(:) = (/ 0.01e-6_r8, 0.1e-6_r8, 1.0e-6_r8, 10.0e-6_r8/) ! Aitken dust
!#if ( defined DUSTTYPES )
!#if ( defined MODAL_AERO_4MODE )
!
!       dust_emis_sclfctr(:) = &
!                              (/ 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, 0.001_r8, &
!                                 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, 0.011_r8, &
!                                 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8, 0.989_r8 /)
!
!       brit_emis_sclfctr(:) = &
!     (/ 1.000_r8, 1.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 0.000_r8,& ! aitken
!        1.000_r8, 1.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 0.000_r8,& ! accum
!        0.695_r8, 0.695_r8, 0.695_r8, 0.305_r8, 0.695_r8, 0.000_r8, 0.695_r8, 0.305_r8, 0.695_r8, 0.305_r8, 0.305_r8/) ! coarse
!!        ^illite^   ^mont^    ^kaol^   ^felds^   ^hem_c^  ^hem_s    ^qua_c    ^qua_s    ^cal_c    ^cal_s   ^ gypsu
!
!       dust_names(:) = (/ 'dst1_a2', 'dst2_a2', 'dst3_a2', 'dst4_a2', 'dst5_a2', 'dst6_a2', 'dst7_a2', 'dst8_a2', &
!                          'dst1_a1', 'dst2_a1', 'dst3_a1', 'dst4_a1', 'dst5_a1', 'dst6_a1', 'dst7_a1', 'dst8_a1', &
!                          'dst1_a3', 'dst2_a3', 'dst3_a3', 'dst4_a3', 'dst5_a3', 'dst6_a3', 'dst7_a3', 'dst8_a3', &
!                          'num_a2 ', 'num_a1 ', 'num_a3 ' /)
!#elif ( defined MODAL_AERO_9MODE )
#if ( defined DUSTTYPES )
#if ( defined MODAL_AERO_10MODE )
    elseif ( ntot_amode == 10 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 1.0e-6_r8, 5.0e-6_r8, &
                                   10.0e-6_r8, 20.0e-6_r8, 70.0e-6_r8/)
       ! Per Meng., et al. (2022, GRL)
       ! derechoFullRunMAM10Size60040804_tuningFeOx01
       dust_emis_sclfctr(:) = &
                          (/ 0.0068_r8, 0.0068_r8, 0.0068_r8, 0.0068_r8, 0.0068_r8, 0.0068_r8, 0.0068_r8, 0.0068_r8, & ! fine
                             0.1488_r8, 0.1488_r8, 0.1488_r8, 0.1488_r8, 0.1488_r8, 0.1488_r8, 0.1488_r8, 0.1488_r8, & ! coarse01
                             0.1721_r8, 0.1721_r8, 0.1721_r8, 0.1721_r8, 0.1721_r8, 0.1721_r8, 0.1721_r8, 0.1721_r8, & ! coarse02
                             0.2405_r8, 0.2405_r8, 0.2405_r8, 0.2405_r8, 0.2405_r8, 0.2405_r8, 0.2405_r8, 0.2405_r8, & ! coarse03
                             0.4318_r8, 0.4318_r8, 0.4318_r8, 0.4318_r8, 0.4318_r8, 0.4318_r8, 0.4318_r8, 0.4318_r8 /) ! giant

     ! tuned hem a bit for derechoFullRunMAM10Size60040804_tuningFeOx03
       brit_emis_sclfctr(:) = &
     (/ 1.000_r8, 1.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 1.000_r8, 0.000_r8, 0.000_r8,& ! fine
        0.970_r8, 0.970_r8, 0.970_r8, 0.030_r8, 0.570_r8, 0.430_r8, 0.970_r8, 0.030_r8, 0.970_r8, 0.030_r8, 0.030_r8,& ! coarse01 
        0.871_r8, 0.871_r8, 0.871_r8, 0.129_r8, 0.571_r8, 0.429_r8, 0.871_r8, 0.129_r8, 0.871_r8, 0.129_r8, 0.871_r8,& ! coarse02
        0.871_r8, 0.871_r8, 0.871_r8, 0.871_r8, 0.079_r8, 0.921_r8, 0.871_r8, 0.129_r8, 0.871_r8, 0.129_r8, 0.871_r8,& ! coarse03 
        0.871_r8, 0.871_r8, 0.871_r8, 0.871_r8, 0.079_r8, 0.921_r8, 0.129_r8, 0.871_r8, 0.129_r8, 0.871_r8, 0.871_r8/) ! giant 
!        ^illite^   ^mont^    ^kaol^   ^felds^   ^hem_c^  ^hem_s    ^qua_c    ^qua_s    ^cal_c    ^cal_s   ^ gypsu

       dust_names(:) = (/ 'dst1_a5', 'dst2_a5', 'dst3_a5', 'dst4_a5', 'dst5_a5', 'dst6_a5', 'dst7_a5', 'dst8_a5', &
                          'dst1_a7', 'dst2_a7', 'dst3_a7', 'dst4_a7', 'dst5_a7', 'dst6_a7', 'dst7_a7', 'dst8_a7', &
                          'dst1_a8', 'dst2_a8', 'dst3_a8', 'dst4_a8', 'dst5_a8', 'dst6_a8', 'dst7_a8', 'dst8_a8', &
                          'dst1_a9', 'dst2_a9', 'dst3_a9', 'dst4_a9', 'dst5_a9', 'dst6_a9', 'dst7_a9', 'dst8_a9', &
                          'dst1_a10', 'dst2_a10', 'dst3_a10', 'dst4_a10', 'dst5_a10', 'dst6_a10', 'dst7_a10', 'dst8_a10', &
                          'num_a5 ', 'num_a7 ', 'num_a8 ', 'num_a9 ', 'num_a10' /)

#endif

#else
       dust_emis_sclfctr(:) = (/ 1.65E-05_r8, 0.011_r8, 0.989_r8 /) ! Aitken dust
#endif

    else if( ntot_amode == 7 ) then
       dust_dmt_grd(:) = (/ 0.1e-6_r8, 2.0e-6_r8, 10.0e-6_r8/)
       dust_emis_sclfctr(:) = (/ 0.13_r8, 0.87_r8 /)
    endif

    do n = 1, dust_nbin
               !write( iulog, '(i2)' ) n

               !write( iulog, '(1x,a)' ) &
               !      'dust_names=',dust_names(n)

       call cnst_get_ind(dust_names(n), dust_indices(n))
    end do


    do n = 1, dust_nnum
       call cnst_get_ind(dust_names(dust_nbin+n), dust_indices(dust_nbin+n))
    enddo

    dust_active = any(dust_indices(:) > 0)
    if (.not.dust_active) return
   

    call  soil_erod_init( dust_emis_fact, soil_erod_file )
#if ( defined DUSTTYPES )
    call dust_set_params( dust_nnum, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )
#else
    call dust_set_params( dust_nbin, dust_dmt_grd, dust_dmt_vwr, dust_stk_crc )
#endif

  end subroutine dust_init

  !===============================================================================
  !===============================================================================
  subroutine dust_emis( ncol, lchnk, dust_flux_in, cflx, soil_erod )
    use soil_erod_mod, only : soil_erod_fact
    use soil_erod_mod, only : soil_erodibility
#if ( defined DUSTTYPES )
    ! Longlei Li - BRIFT already applied offline
    use soil_erod_mod, only :  illite_frac_accum, illite_frac_aitken, illite_frac_coarse01
    use soil_erod_mod, only :  kaolinite_frac_accum, kaolinite_frac_aitken, kaolinite_frac_coarse01
    use soil_erod_mod, only :  mont_frac_accum, mont_frac_aitken, mont_frac_coarse01
    use soil_erod_mod, only :  hematite_frac_accum, hematite_frac_aitken, hematite_frac_coarse01
    use soil_erod_mod, only :  quartz_frac_accum, quartz_frac_aitken, quartz_frac_coarse01
    use soil_erod_mod, only :  calcite_frac_accum, calcite_frac_aitken, calcite_frac_coarse01
    use soil_erod_mod, only :  feldspar_frac_accum, feldspar_frac_aitken, feldspar_frac_coarse01
    use soil_erod_mod, only :  gypsum_frac_accum, gypsum_frac_aitken, gypsum_frac_coarse01

    use soil_erod_mod, only :  illite_frac_coarse02, illite_frac_giant
    use soil_erod_mod, only :  kaolinite_frac_coarse02, kaolinite_frac_giant
    use soil_erod_mod, only :  mont_frac_coarse02, mont_frac_giant 
    use soil_erod_mod, only :  hematite_frac_coarse02, hematite_frac_giant
    use soil_erod_mod, only :  quartz_frac_coarse02, quartz_frac_giant
    use soil_erod_mod, only :  calcite_frac_coarse02, calcite_frac_giant
    use soil_erod_mod, only :  feldspar_frac_coarse02, feldspar_frac_giant
    use soil_erod_mod, only :  gypsum_frac_coarse02, gypsum_frac_giant

    ! apply brittle theory online -Longlei Li
    !use soil_erod_mod, only :  illite_clay,    kaolinite_clay,   mont_clay
    !use soil_erod_mod, only :  hematite_clay,  hematite_silt,    quartz_clay,   quartz_silt
    !use soil_erod_mod, only :  calcite_clay,   calcite_silt,     feldspar_silt, gypsum_silt
#endif
    use mo_constants,  only : dust_density
    use physconst,     only : pi
    use cam_history,   only : outfld
    use ppgrid,        only : pcols
    use cam_logfile,       only: iulog

  ! args
    integer,  intent(in)    :: ncol, lchnk
    real(r8), intent(in)    :: dust_flux_in(:,:)
    real(r8), intent(inout) :: cflx(:,:)
    real(r8), intent(out)   :: soil_erod(:)

  ! local vars
    integer :: i, m, idst, inum
    real(r8) :: x_mton
    real(r8),parameter :: soil_erod_threshold = 0.1_r8
#if ( defined DUSTTYPES )

! ++ Longlei Li MAM10 ++ 
#if ( defined MODAL_AERO_10MODE )
    real(r8) ::  frac_I_accum(pcols), frac_I_aitken(pcols), frac_I_coarse01(pcols), frac_I_coarse02(pcols), frac_I_giant(pcols)
    real(r8) ::  frac_K_accum(pcols), frac_K_aitken(pcols), frac_K_coarse01(pcols), frac_K_coarse02(pcols), frac_K_giant(pcols)
    real(r8) ::  frac_M_accum(pcols), frac_M_aitken(pcols), frac_M_coarse01(pcols), frac_M_coarse02(pcols), frac_M_giant(pcols)
    real(r8) ::  frac_H_accum(pcols), frac_H_aitken(pcols), frac_H_coarse01(pcols), frac_H_coarse02(pcols), frac_H_giant(pcols)
    real(r8) ::  frac_Q_accum(pcols), frac_Q_aitken(pcols), frac_Q_coarse01(pcols), frac_Q_coarse02(pcols), frac_Q_giant(pcols)
    real(r8) ::  frac_C_accum(pcols), frac_C_aitken(pcols), frac_C_coarse01(pcols), frac_C_coarse02(pcols), frac_C_giant(pcols)
    real(r8) ::  frac_F_accum(pcols), frac_F_aitken(pcols), frac_F_coarse01(pcols), frac_F_coarse02(pcols), frac_F_giant(pcols)
    real(r8) ::  frac_G_accum(pcols), frac_G_aitken(pcols), frac_G_coarse01(pcols), frac_G_coarse02(pcols), frac_G_giant(pcols)
    real(r8) ::  dust_flux_sum(pcols)
#endif
! -- Longlei Li -- 

#endif

#if ( defined DUSTTYPES )

#if ( defined MODAL_AERO_10MODE )
    do i = 1, ncol
       ! Longlei Li: BRIFT applied already offline
       frac_I_aitken(i)   = illite_frac_aitken(i,lchnk)
       frac_I_accum(i)    = illite_frac_accum(i,lchnk)
       frac_I_coarse01(i) = illite_frac_coarse01(i,lchnk)
       frac_I_coarse02(i) = illite_frac_coarse02(i,lchnk)
       frac_I_giant(i)    = illite_frac_giant(i,lchnk)

       frac_M_aitken(i)   = mont_frac_aitken(i,lchnk)
       frac_M_accum(i)    = mont_frac_accum(i,lchnk)
       frac_M_coarse01(i) = mont_frac_coarse01(i,lchnk)
       frac_M_coarse02(i) = mont_frac_coarse02(i,lchnk)
       frac_M_giant(i)    = mont_frac_giant(i,lchnk)

       frac_K_aitken(i)   = kaolinite_frac_aitken(i,lchnk)
       frac_K_accum(i)    = kaolinite_frac_accum(i,lchnk)
       frac_K_coarse01(i) = kaolinite_frac_coarse01(i,lchnk)
       frac_K_coarse02(i) = kaolinite_frac_coarse02(i,lchnk)
       frac_K_giant(i)    = kaolinite_frac_giant(i,lchnk)

       frac_F_aitken(i)   = feldspar_frac_aitken(i,lchnk)
       frac_F_accum(i)    = feldspar_frac_accum(i,lchnk)
       frac_F_coarse01(i) = feldspar_frac_coarse01(i,lchnk)
       frac_F_coarse02(i) = feldspar_frac_coarse02(i,lchnk)
       frac_F_giant(i)    = feldspar_frac_giant(i,lchnk)

       frac_H_aitken(i)   = hematite_frac_aitken(i,lchnk)
       frac_H_accum(i)    = hematite_frac_accum(i,lchnk) 
       frac_H_coarse01(i) = hematite_frac_coarse01(i,lchnk)
       frac_H_coarse02(i) = hematite_frac_coarse02(i,lchnk)
       frac_H_giant(i)    = hematite_frac_giant(i,lchnk)

       frac_Q_aitken(i)   = quartz_frac_aitken(i,lchnk)
       frac_Q_accum(i)    = quartz_frac_accum(i,lchnk)
       frac_Q_coarse01(i) = quartz_frac_coarse01(i,lchnk)
       frac_Q_coarse02(i) = quartz_frac_coarse02(i,lchnk)
       frac_Q_giant(i)    = quartz_frac_giant(i,lchnk)

       frac_C_aitken(i)   = calcite_frac_aitken(i,lchnk)
       frac_C_accum(i)    = calcite_frac_accum(i,lchnk)
       frac_C_coarse01(i) = calcite_frac_coarse01(i,lchnk)
       frac_C_coarse02(i) = calcite_frac_coarse02(i,lchnk)
       frac_C_giant(i)    = calcite_frac_giant(i,lchnk)

       frac_G_aitken(i)   = gypsum_frac_aitken(i,lchnk)
       frac_G_accum(i)    = gypsum_frac_accum(i,lchnk)
       frac_G_coarse01(i) = gypsum_frac_coarse01(i,lchnk)
       frac_G_coarse02(i) = gypsum_frac_coarse02(i,lchnk)
       frac_G_giant(i)    = gypsum_frac_giant(i,lchnk)
    end do

    call outfld('Illite_accum',frac_I_accum(:),pcols,lchnk)
    call outfld('Illite_aitken',frac_I_aitken(:),pcols,lchnk)
    call outfld('Illite_coarse01',frac_I_coarse01(:),pcols,lchnk)
    call outfld('Illite_coarse02',frac_I_coarse02(:),pcols,lchnk)
    call outfld('Illite_giant',frac_I_giant(:),pcols,lchnk)

    call outfld('Kaolinite_accum',frac_K_accum(:),pcols,lchnk)
    call outfld('Kaolinite_aitken',frac_K_aitken(:),pcols,lchnk)
    call outfld('Kaolinite_coarse01',frac_K_coarse01(:),pcols,lchnk)
    call outfld('Kaolinite_coarse02',frac_K_coarse02(:),pcols,lchnk)
    call outfld('Kaolinite_giant',frac_K_giant(:),pcols,lchnk)

    call outfld('Mont_accum',frac_M_accum(:),pcols,lchnk)
    call outfld('Mont_aitken',frac_M_aitken(:),pcols,lchnk)
    call outfld('Mont_coarse01',frac_M_coarse01(:),pcols,lchnk)
    call outfld('Mont_coarse02',frac_M_coarse02(:),pcols,lchnk)
    call outfld('Mont_giant',frac_M_giant(:),pcols,lchnk)

    call outfld('Hematite_accum',frac_H_accum(:),pcols,lchnk)
    call outfld('Hematite_aitken',frac_H_aitken(:),pcols,lchnk)
    call outfld('Hematite_coarse01',frac_H_coarse01(:),pcols,lchnk)
    call outfld('Hematite_coarse02',frac_H_coarse02(:),pcols,lchnk)
    call outfld('Hematite_giant',frac_H_giant(:),pcols,lchnk)

    call outfld('Quartz_accum',frac_Q_accum(:),pcols,lchnk)
    call outfld('Quartz_aitken',frac_Q_aitken(:),pcols,lchnk)
    call outfld('Quartz_coarse01',frac_Q_coarse01(:),pcols,lchnk)
    call outfld('Quartz_coarse02',frac_Q_coarse02(:),pcols,lchnk)
    call outfld('Quartz_giant',frac_Q_giant(:),pcols,lchnk)

    call outfld('Calcite_accum',frac_C_accum(:),pcols,lchnk)
    call outfld('Calcite_aitken',frac_C_aitken(:),pcols,lchnk)
    call outfld('Calcite_coarse01',frac_C_coarse01(:),pcols,lchnk)
    call outfld('Calcite_coarse02',frac_C_coarse02(:),pcols,lchnk)
    call outfld('Calcite_giant',frac_C_giant(:),pcols,lchnk)

    call outfld('Feldspar_accum',frac_F_accum(:),pcols,lchnk)
    call outfld('Feldspar_aitken',frac_F_aitken(:),pcols,lchnk)
    call outfld('Feldspar_coarse01',frac_F_coarse01(:),pcols,lchnk)
    call outfld('Feldspar_coarse02',frac_F_coarse02(:),pcols,lchnk)
    call outfld('Feldspar_giant',frac_F_giant(:),pcols,lchnk)

    call outfld('Gypsum_accum',frac_G_accum(:),pcols,lchnk)
    call outfld('Gypsum_aitken',frac_G_aitken(:),pcols,lchnk)
    call outfld('Gypsum_coarse01',frac_G_coarse01(:),pcols,lchnk)
    call outfld('Gypsum_coarse02',frac_G_coarse02(:),pcols,lchnk)
    call outfld('Gypsum_giant',frac_G_giant(:),pcols,lchnk)

#endif

#endif

    ! set dust emissions
#if ( defined DUSTTYPES )

#if ( defined MODAL_AERO_10MODE )
    col_loop: do i =1,ncol
       dust_flux_sum(i) = sum( -dust_flux_in(i,:) )

       !soil_erod(i) = soil_erodibility( i, lchnk )
       soil_erod(i) = 1._r8 ! Kok scheme

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..

       ! Longlei Li: BRIFT applied already offline

       ! Aitken mode minerals
       idst = dust_indices(1)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_I_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(2)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_K_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(3)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_M_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(4)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_H_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(5)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_Q_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(6)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_C_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(7)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_F_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(8)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_G_aitken(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       ! Accumulation mode minerals
       idst = dust_indices(9)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_I_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(10)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_K_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(11)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_M_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(12)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_H_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(13)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_Q_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(14)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_C_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(15)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_F_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(16)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_G_accum(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       ! Coarse01 mode minerals
       idst = dust_indices(17)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_I_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(18)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_K_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(19)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_M_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(20)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_H_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(21)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_Q_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(22)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_C_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(23)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_F_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(24)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_G_coarse01(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       ! Coarse02 mode minerals
       idst = dust_indices(25)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_I_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(26)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_K_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(27)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_M_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(28)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_H_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(29)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_Q_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(30)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_C_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(31)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_F_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(32)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_G_coarse02(i)*soil_erod(i)/soil_erod_fact*1.15_r8


       ! Giant mode minerals
       idst = dust_indices(33)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_I_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(34)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_K_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(35)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_M_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(36)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_H_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(37)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_Q_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(38)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_C_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(39)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_F_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8
       idst = dust_indices(40)
       cflx(i,idst) = sum( -dust_flux_in(i,:) )*frac_G_giant(i)*soil_erod(i)/soil_erod_fact*1.15_r8

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(1)**3._r8))
       inum = dust_indices(1+dust_nbin)
       cflx(i,inum) = (cflx(i,dust_indices(1)) + cflx(i,dust_indices(2)) + cflx(i,dust_indices(3)) + cflx(i,dust_indices(4)) + &
                       cflx(i,dust_indices(5)) + cflx(i,dust_indices(6)) + cflx(i,dust_indices(7)) + cflx(i,dust_indices(8))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(2)**3._r8))
       inum = dust_indices(2+dust_nbin)
       cflx(i,inum) = (cflx(i,dust_indices(9)) + cflx(i,dust_indices(10)) + cflx(i,dust_indices(11)) + cflx(i,dust_indices(12)) + &
                       cflx(i,dust_indices(13))+ cflx(i,dust_indices(14)) + cflx(i,dust_indices(15)) + cflx(i,dust_indices(16))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(3)**3._r8))
       inum = dust_indices(3+dust_nbin)
       cflx(i,inum) =(cflx(i,dust_indices(17))+ cflx(i,dust_indices(18)) + cflx(i,dust_indices(19)) + cflx(i,dust_indices(20)) + &
                      cflx(i,dust_indices(21))+ cflx(i,dust_indices(22)) + cflx(i,dust_indices(23)) + cflx(i,dust_indices(24))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(4)**3._r8))
       inum = dust_indices(4+dust_nbin)
       cflx(i,inum) =(cflx(i,dust_indices(25))+ cflx(i,dust_indices(26)) + cflx(i,dust_indices(27)) + cflx(i,dust_indices(28)) + &
                      cflx(i,dust_indices(29))+ cflx(i,dust_indices(30)) + cflx(i,dust_indices(31)) + cflx(i,dust_indices(32))   &
                       )*x_mton

       x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(5)**3._r8))
       inum = dust_indices(5+dust_nbin)
       cflx(i,inum) =(cflx(i,dust_indices(33))+ cflx(i,dust_indices(34)) + cflx(i,dust_indices(35)) + cflx(i,dust_indices(36)) + &
                      cflx(i,dust_indices(37))+ cflx(i,dust_indices(38)) + cflx(i,dust_indices(39)) + cflx(i,dust_indices(40))   &
                       )*x_mton

    end do col_loop
#endif

#else
    col_loop: do i =1,ncol

       soil_erod(i) = soil_erodibility( i, lchnk )

       if( soil_erod(i) .lt. soil_erod_threshold ) soil_erod(i) = 0._r8

       ! rebin and adjust dust emissons..
       do m = 1,dust_nbin

          idst = dust_indices(m)

          cflx(i,idst) = sum( -dust_flux_in(i,:) ) &
               * dust_emis_sclfctr(m)*soil_erod(i)/soil_erod_fact*1.15_r8

          x_mton = 6._r8 / (pi * dust_density * (dust_dmt_vwr(m)**3._r8))                

          inum = dust_indices(m+dust_nbin)

          cflx(i,inum) = cflx(i,idst)*x_mton

       enddo

    end do col_loop
#endif

#if ( defined DUSTTYPES )
  call outfld('dust_flux',dust_flux_sum(:),pcols,lchnk)
#endif

  end subroutine dust_emis

end module dust_model
