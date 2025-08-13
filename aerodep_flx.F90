!-------------------------------------------------------------------
! Manages reading and interpolation of prescribed aerosol deposition 
! fluxes.  These are the deposition fluxes sent to the surface.
!
! Created by: Francis Vitt
!-------------------------------------------------------------------
module aerodep_flx

  use shr_kind_mod,     only : r8 => shr_kind_r8
  use cam_abortutils,   only : endrun
  use spmd_utils,       only : masterproc
  use tracer_data,      only : trfld, trfile
  use cam_logfile,      only : iulog
  use ppgrid,           only : pcols, pver, begchunk, endchunk

  implicit none
  private
  save 

  type(trfld), pointer :: fields(:)
  type(trfile)         :: file

  public :: aerodep_flx_init
  public :: aerodep_flx_adv
  public :: aerodep_flx_readnl
  public :: aerodep_flx_prescribed

  logical :: has_aerodep_flx = .false.
  integer, parameter, public :: N_BULK = 14
#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
  integer, parameter, public :: N_MODAL = 174 !142 
#else
  integer, parameter, public :: N_MODAL = 78
#endif
#else
  integer, parameter, public :: N_MODAL = 22
#endif
  integer :: number_flds

  character(len=256) :: filename = 'NONE'
  character(len=256) :: filelist = ' '
  character(len=256) :: datapath = ' '
  character(len=32)  :: datatype = 'SERIAL'
  logical            :: rmv_file = .false.
  integer            :: cycle_yr = 0
  integer            :: fixed_ymd = 0
  integer            :: fixed_tod = 0
  character(len=32)  :: specifier(N_MODAL) = ' '

  ! for bulk aerosol fluxes

  character(len=12), parameter :: bulk_names(N_BULK) = (/ &
       'BCDEPWET    ', 'BCPHODRY    ', 'BCPHIDRY    ',  &
       'OCDEPWET    ', 'OCPHODRY    ', 'OCPHIDRY    ',  &
       'DSTX01DD    ', 'DSTX02DD    ', 'DSTX03DD    ', 'DSTX04DD    ', &
       'DSTX01WD    ', 'DSTX02WD    ', 'DSTX03WD    ', 'DSTX04WD    ' /)

  integer :: index_bulk_map(N_BULK)

  integer :: ibcphiwet,ibcphidry,ibcphodry
  integer :: iocphiwet,iocphidry,iocphodry

  integer :: idstdry1,idstdry2,idstdry3,idstdry4
  integer :: idstwet1,idstwet2,idstwet3,idstwet4

  ! for modal aerosol fluxes
#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
  character(len=13), parameter :: modal_names(N_MODAL) = (/ &
       'bc_a1DDF     ', 'bc_c1DDF     ', 'pom_a1DDF    ', 'pom_c1DDF    ',  &
       'soa_a1DDF    ', 'soa_c1DDF    ', 'soa_a2DDF    ', 'soa_c2DDF    ',  &
       'dst1_a5DDF   ', 'dst1_c5DDF   ', 'dst1_a7DDF   ', 'dst1_c7DDF   ',  &
       'dst2_a5DDF   ', 'dst2_c5DDF   ', 'dst2_a7DDF   ', 'dst2_c7DDF   ',  &
       'dst3_a5DDF   ', 'dst3_c5DDF   ', 'dst3_a7DDF   ', 'dst3_c7DDF   ',  &
       'dst4_a5DDF   ', 'dst4_c5DDF   ', 'dst4_a7DDF   ', 'dst4_c7DDF   ',  &
       'dst5_a5DDF   ', 'dst5_c5DDF   ', 'dst5_a7DDF   ', 'dst5_c7DDF   ',  &
       'dst6_a5DDF   ', 'dst6_c5DDF   ', 'dst6_a7DDF   ', 'dst6_c7DDF   ',  &
       'dst7_a5DDF   ', 'dst7_c5DDF   ', 'dst7_a7DDF   ', 'dst7_c7DDF   ',  &
       'dst8_a5DDF   ', 'dst8_c5DDF   ', 'dst8_a7DDF   ', 'dst8_c7DDF   ',  &

       'dst1_a8DDF   ', 'dst1_c8DDF   ', 'dst1_a9DDF   ', 'dst1_c9DDF   ',  &
       'dst2_a8DDF   ', 'dst2_c8DDF   ', 'dst2_a9DDF   ', 'dst2_c9DDF   ',  &
       'dst3_a8DDF   ', 'dst3_c8DDF   ', 'dst3_a9DDF   ', 'dst3_c9DDF   ',  &
       'dst4_a8DDF   ', 'dst4_c8DDF   ', 'dst4_a9DDF   ', 'dst4_c9DDF   ',  &
       'dst5_a8DDF   ', 'dst5_c8DDF   ', 'dst5_a9DDF   ', 'dst5_c9DDF   ',  &
       'dst6_a8DDF   ', 'dst6_c8DDF   ', 'dst6_a9DDF   ', 'dst6_c9DDF   ',  &
       'dst7_a8DDF   ', 'dst7_c8DDF   ', 'dst7_a9DDF   ', 'dst7_c9DDF   ',  &
       'dst8_a8DDF   ', 'dst8_c8DDF   ', 'dst8_a9DDF   ', 'dst8_c9DDF   ',  &

       'dst1_a10DDF  ', 'dst1_c10DDF  ',  &
       'dst2_a10DDF  ', 'dst2_c10DDF  ',  &
       'dst3_a10DDF  ', 'dst3_c10DDF  ',  &
       'dst4_a10DDF  ', 'dst4_c10DDF  ',  &
       'dst5_a10DDF  ', 'dst5_c10DDF  ',  &
       'dst6_a10DDF  ', 'dst6_c10DDF  ',  &
       'dst7_a10DDF  ', 'dst7_c10DDF  ',  &
       'dst8_a10DDF  ', 'dst8_c10DDF  ',  &

       'bc_a1SFWET   ', 'bc_c1SFWET   ',  'pom_a1SFWET  ', 'pom_c1SFWET  ', &
       'soa_a1SFWET  ', 'soa_c1SFWET  ',                                    &
       'dst1_a5SFWET ', 'dst1_c5SFWET ',  'dst1_a7SFWET ', 'dst1_c7SFWET ', &
       'dst2_a5SFWET ', 'dst2_c5SFWET ',  'dst2_a7SFWET ', 'dst2_c7SFWET ', &
       'dst3_a5SFWET ', 'dst3_c5SFWET ',  'dst3_a7SFWET ', 'dst3_c7SFWET ', &
       'dst4_a5SFWET ', 'dst4_c5SFWET ',  'dst4_a7SFWET ', 'dst4_c7SFWET ', &
       'dst5_a5SFWET ', 'dst5_c5SFWET ',  'dst5_a7SFWET ', 'dst5_c7SFWET ', &
       'dst6_a5SFWET ', 'dst6_c5SFWET ',  'dst6_a7SFWET ', 'dst6_c7SFWET ', &
       'dst7_a5SFWET ', 'dst7_c5SFWET ',  'dst7_a7SFWET ', 'dst7_c7SFWET ', &
       'dst8_a5SFWET ', 'dst8_c5SFWET ',  'dst8_a7SFWET ', 'dst8_c7SFWET ', &

       'dst1_a8SFWET ', 'dst1_c8SFWET ',  'dst1_a9SFWET ', 'dst1_c9SFWET ', &
       'dst2_a8SFWET ', 'dst2_c8SFWET ',  'dst2_a9SFWET ', 'dst2_c9SFWET ', &
       'dst3_a8SFWET ', 'dst3_c8SFWET ',  'dst3_a9SFWET ', 'dst3_c9SFWET ', &
       'dst4_a8SFWET ', 'dst4_c8SFWET ',  'dst4_a9SFWET ', 'dst4_c9SFWET ', &
       'dst5_a8SFWET ', 'dst5_c8SFWET ',  'dst5_a9SFWET ', 'dst5_c9SFWET ', &
       'dst6_a8SFWET ', 'dst6_c8SFWET ',  'dst6_a9SFWET ', 'dst6_c9SFWET ', &
       'dst7_a8SFWET ', 'dst7_c8SFWET ',  'dst7_a9SFWET ', 'dst7_c9SFWET ', &
       'dst8_a8SFWET ', 'dst8_c8SFWET ',  'dst8_a9SFWET ', 'dst8_c9SFWET ', &

       'dst1_a10SFWET', 'dst1_c10SFWET',  &
       'dst2_a10SFWET', 'dst2_c10SFWET',  &
       'dst3_a10SFWET', 'dst3_c10SFWET',  &
       'dst4_a10SFWET', 'dst4_c10SFWET',  &
       'dst5_a10SFWET', 'dst5_c10SFWET',  &
       'dst6_a10SFWET', 'dst6_c10SFWET',  &
       'dst7_a10SFWET', 'dst7_c10SFWET',  &
       'dst8_a10SFWET', 'dst8_c10SFWET'/)

#endif
#else
  character(len=12), parameter :: modal_names(N_MODAL) = (/ &
       'bc_a1DDF    ', 'bc_c1DDF    ', 'pom_a1DDF   ', 'pom_c1DDF   ',  &
       'soa_a1DDF   ', 'soa_c1DDF   ', 'soa_a2DDF   ', 'soa_c2DDF   ',  &
       'dst_a1DDF   ', 'dst_c1DDF   ', 'dst_a3DDF   ', 'dst_c3DDF   ',  &
       'bc_a1SFWET  ', 'bc_c1SFWET  ', 'pom_a1SFWET ', 'pom_c1SFWET ',  &
       'soa_a1SFWET ', 'soa_c1SFWET ', 'dst_a1SFWET ', 'dst_c1SFWET ',  &
       'dst_a3SFWET ', 'dst_c3SFWET ' /)
#endif

  integer :: index_modal_map(N_MODAL)

#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4

  integer, parameter :: idx_dst15 = 5
  integer, parameter :: idx_dst25 = 6
  integer, parameter :: idx_dst35 = 7
  integer, parameter :: idx_dst45 = 8
  integer, parameter :: idx_dst55 = 9 
  integer, parameter :: idx_dst65 = 10
  integer, parameter :: idx_dst75 = 11
  integer, parameter :: idx_dst85 = 12

  integer, parameter :: idx_dst17 = 13 
  integer, parameter :: idx_dst27 = 14
  integer, parameter :: idx_dst37 = 15
  integer, parameter :: idx_dst47 = 16
  integer, parameter :: idx_dst57 = 17
  integer, parameter :: idx_dst67 = 18 
  integer, parameter :: idx_dst77 = 19 
  integer, parameter :: idx_dst87 = 20

  integer, parameter :: idx_dst18 = 21
  integer, parameter :: idx_dst28 = 22
  integer, parameter :: idx_dst38 = 23
  integer, parameter :: idx_dst48 = 24
  integer, parameter :: idx_dst58 = 25
  integer, parameter :: idx_dst68 = 26
  integer, parameter :: idx_dst78 = 27
  integer, parameter :: idx_dst88 = 28 

  integer, parameter :: idx_dst19 = 29
  integer, parameter :: idx_dst29 = 30
  integer, parameter :: idx_dst39 = 31
  integer, parameter :: idx_dst49 = 32
  integer, parameter :: idx_dst59 = 33
  integer, parameter :: idx_dst69 = 34
  integer, parameter :: idx_dst79 = 35
  integer, parameter :: idx_dst89 = 36

  integer, parameter :: idx_dst110 = 37 
  integer, parameter :: idx_dst210 = 38
  integer, parameter :: idx_dst310 = 39
  integer, parameter :: idx_dst410 = 40
  integer, parameter :: idx_dst510 = 41
  integer, parameter :: idx_dst610 = 42
  integer, parameter :: idx_dst710 = 43
  integer, parameter :: idx_dst810 = 44

  integer, parameter :: idx_ncl1 = 45 
  integer, parameter :: idx_so41 = 46 

#else
  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4
  integer, parameter :: idx_dst11 = 5
  integer, parameter :: idx_dst21 = 6
  integer, parameter :: idx_dst31 = 7
  integer, parameter :: idx_dst41 = 8
  integer, parameter :: idx_dst51 = 9
  integer, parameter :: idx_dst61 = 10
  integer, parameter :: idx_dst71 = 11
  integer, parameter :: idx_dst81 = 12
  integer, parameter :: idx_dst13 = 13
  integer, parameter :: idx_dst23 = 14
  integer, parameter :: idx_dst33 = 15
  integer, parameter :: idx_dst43 = 16
  integer, parameter :: idx_dst53 = 17
  integer, parameter :: idx_dst63 = 18
  integer, parameter :: idx_dst73 = 19
  integer, parameter :: idx_dst83 = 20
  integer, parameter :: idx_ncl3 = 21
  integer, parameter :: idx_so43 = 22
#endif
#else
  integer, parameter :: idx_bc1 = 1
  integer, parameter :: idx_pom1 = 2
  integer, parameter :: idx_soa1 = 3
  integer, parameter :: idx_soa2 = 4
  integer, parameter :: idx_dst1 = 5
  integer, parameter :: idx_dst3 = 6
  integer, parameter :: idx_ncl3 = 7
  integer, parameter :: idx_so43 = 8
#endif

#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
  integer, parameter :: nmodal_idxs = 46 
#else
  integer, parameter :: nmodal_idxs = 22
#endif
#else
  integer, parameter :: nmodal_idxs = 8
#endif

  integer :: idx_bc1_dryis = -1
  integer :: idx_bc1_drycw = -1
  integer :: idx_pom1_dryis = -1
  integer :: idx_pom1_drycw = -1
  integer :: idx_soa1_dryis = -1
  integer :: idx_soa1_drycw = -1
  integer :: idx_soa2_dryis = -1
  integer :: idx_soa2_drycw = -1
#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
  integer :: idx_dst15_dryis = -1
  integer :: idx_dst15_drycw = -1
  integer :: idx_dst25_dryis = -1
  integer :: idx_dst25_drycw = -1
  integer :: idx_dst35_dryis = -1
  integer :: idx_dst35_drycw = -1
  integer :: idx_dst45_dryis = -1
  integer :: idx_dst45_drycw = -1
  integer :: idx_dst55_dryis = -1
  integer :: idx_dst55_drycw = -1
  integer :: idx_dst65_dryis = -1
  integer :: idx_dst65_drycw = -1
  integer :: idx_dst75_dryis = -1
  integer :: idx_dst75_drycw = -1
  integer :: idx_dst85_dryis = -1
  integer :: idx_dst85_drycw = -1

  integer :: idx_dst17_dryis = -1
  integer :: idx_dst17_drycw = -1
  integer :: idx_dst27_dryis = -1
  integer :: idx_dst27_drycw = -1
  integer :: idx_dst37_dryis = -1
  integer :: idx_dst37_drycw = -1
  integer :: idx_dst47_dryis = -1
  integer :: idx_dst47_drycw = -1
  integer :: idx_dst57_dryis = -1
  integer :: idx_dst57_drycw = -1
  integer :: idx_dst67_dryis = -1
  integer :: idx_dst67_drycw = -1
  integer :: idx_dst77_dryis = -1
  integer :: idx_dst77_drycw = -1
  integer :: idx_dst87_dryis = -1
  integer :: idx_dst87_drycw = -1

  integer :: idx_dst18_dryis = -1
  integer :: idx_dst18_drycw = -1
  integer :: idx_dst28_dryis = -1
  integer :: idx_dst28_drycw = -1
  integer :: idx_dst38_dryis = -1
  integer :: idx_dst38_drycw = -1
  integer :: idx_dst48_dryis = -1
  integer :: idx_dst48_drycw = -1
  integer :: idx_dst58_dryis = -1
  integer :: idx_dst58_drycw = -1
  integer :: idx_dst68_dryis = -1
  integer :: idx_dst68_drycw = -1
  integer :: idx_dst78_dryis = -1
  integer :: idx_dst78_drycw = -1
  integer :: idx_dst88_dryis = -1
  integer :: idx_dst88_drycw = -1

  integer :: idx_dst19_dryis = -1
  integer :: idx_dst19_drycw = -1
  integer :: idx_dst29_dryis = -1
  integer :: idx_dst29_drycw = -1
  integer :: idx_dst39_dryis = -1
  integer :: idx_dst39_drycw = -1
  integer :: idx_dst49_dryis = -1
  integer :: idx_dst49_drycw = -1
  integer :: idx_dst59_dryis = -1
  integer :: idx_dst59_drycw = -1
  integer :: idx_dst69_dryis = -1
  integer :: idx_dst69_drycw = -1
  integer :: idx_dst79_dryis = -1
  integer :: idx_dst79_drycw = -1
  integer :: idx_dst89_dryis = -1
  integer :: idx_dst89_drycw = -1

  integer :: idx_dst110_dryis = -1
  integer :: idx_dst110_drycw = -1
  integer :: idx_dst210_dryis = -1
  integer :: idx_dst210_drycw = -1
  integer :: idx_dst310_dryis = -1
  integer :: idx_dst310_drycw = -1
  integer :: idx_dst410_dryis = -1
  integer :: idx_dst410_drycw = -1
  integer :: idx_dst510_dryis = -1
  integer :: idx_dst510_drycw = -1
  integer :: idx_dst610_dryis = -1
  integer :: idx_dst610_drycw = -1
  integer :: idx_dst710_dryis = -1
  integer :: idx_dst710_drycw = -1
  integer :: idx_dst810_dryis = -1
  integer :: idx_dst810_drycw = -1
#endif
#else  
  integer :: idx_dst1_dryis = -1
  integer :: idx_dst1_drycw = -1
  integer :: idx_dst3_dryis = -1
  integer :: idx_dst3_drycw = -1
#endif

  integer :: idx_bc1_wetis = -1
  integer :: idx_bc1_wetcw = -1
  integer :: idx_pom1_wetis = -1
  integer :: idx_pom1_wetcw = -1
  integer :: idx_soa1_wetis = -1
  integer :: idx_soa1_wetcw = -1
#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
  integer :: idx_dst15_wetis = -1
  integer :: idx_dst15_wetcw = -1
  integer :: idx_dst25_wetis = -1
  integer :: idx_dst25_wetcw = -1
  integer :: idx_dst35_wetis = -1
  integer :: idx_dst35_wetcw = -1
  integer :: idx_dst45_wetis = -1
  integer :: idx_dst45_wetcw = -1
  integer :: idx_dst55_wetis = -1
  integer :: idx_dst55_wetcw = -1
  integer :: idx_dst65_wetis = -1
  integer :: idx_dst65_wetcw = -1
  integer :: idx_dst75_wetis = -1
  integer :: idx_dst75_wetcw = -1
  integer :: idx_dst85_wetis = -1
  integer :: idx_dst85_wetcw = -1

  integer :: idx_dst17_wetis = -1
  integer :: idx_dst17_wetcw = -1
  integer :: idx_dst27_wetis = -1
  integer :: idx_dst27_wetcw = -1
  integer :: idx_dst37_wetis = -1
  integer :: idx_dst37_wetcw = -1
  integer :: idx_dst47_wetis = -1
  integer :: idx_dst47_wetcw = -1
  integer :: idx_dst57_wetis = -1
  integer :: idx_dst57_wetcw = -1
  integer :: idx_dst67_wetis = -1
  integer :: idx_dst67_wetcw = -1
  integer :: idx_dst77_wetis = -1
  integer :: idx_dst77_wetcw = -1
  integer :: idx_dst87_wetis = -1
  integer :: idx_dst87_wetcw = -1

  integer :: idx_dst18_wetis = -1
  integer :: idx_dst18_wetcw = -1
  integer :: idx_dst28_wetis = -1
  integer :: idx_dst28_wetcw = -1
  integer :: idx_dst38_wetis = -1
  integer :: idx_dst38_wetcw = -1
  integer :: idx_dst48_wetis = -1
  integer :: idx_dst48_wetcw = -1
  integer :: idx_dst58_wetis = -1
  integer :: idx_dst58_wetcw = -1
  integer :: idx_dst68_wetis = -1
  integer :: idx_dst68_wetcw = -1
  integer :: idx_dst78_wetis = -1
  integer :: idx_dst78_wetcw = -1
  integer :: idx_dst88_wetis = -1
  integer :: idx_dst88_wetcw = -1

  integer :: idx_dst19_wetis = -1
  integer :: idx_dst19_wetcw = -1
  integer :: idx_dst29_wetis = -1
  integer :: idx_dst29_wetcw = -1
  integer :: idx_dst39_wetis = -1
  integer :: idx_dst39_wetcw = -1
  integer :: idx_dst49_wetis = -1
  integer :: idx_dst49_wetcw = -1
  integer :: idx_dst59_wetis = -1
  integer :: idx_dst59_wetcw = -1
  integer :: idx_dst69_wetis = -1
  integer :: idx_dst69_wetcw = -1
  integer :: idx_dst79_wetis = -1
  integer :: idx_dst79_wetcw = -1
  integer :: idx_dst89_wetis = -1
  integer :: idx_dst89_wetcw = -1

  integer :: idx_dst110_wetis = -1
  integer :: idx_dst110_wetcw = -1
  integer :: idx_dst210_wetis = -1
  integer :: idx_dst210_wetcw = -1
  integer :: idx_dst310_wetis = -1
  integer :: idx_dst310_wetcw = -1
  integer :: idx_dst410_wetis = -1
  integer :: idx_dst410_wetcw = -1
  integer :: idx_dst510_wetis = -1
  integer :: idx_dst510_wetcw = -1
  integer :: idx_dst610_wetis = -1
  integer :: idx_dst610_wetcw = -1
  integer :: idx_dst710_wetis = -1
  integer :: idx_dst710_wetcw = -1
  integer :: idx_dst810_wetis = -1
  integer :: idx_dst810_wetcw = -1
#endif
#else  
  integer :: idx_dst1_wetis = -1
  integer :: idx_dst1_wetcw = -1
  integer :: idx_dst3_wetis = -1
  integer :: idx_dst3_wetcw = -1
#endif
  logical :: modal_fluxes = .false.

contains

!-------------------------------------------------------------------
! parses the list of dep fluxes specified in aerodep_flx_specifier namelist
! variable and sets up index variables
!-------------------------------------------------------------------
  subroutine aerodep_flx_init()

    use tracer_data, only : trcdata_init
    use cam_history, only : addfld, horiz_only
    use physics_buffer, only : physics_buffer_desc
    use modal_aero_deposition, only : modal_aero_deposition_init

    implicit none

    integer :: ndx, istat, i

    if ( has_aerodep_flx ) then
       if ( masterproc ) then
          write(iulog,*) 'aero dep fluxes are prescribed in :'//trim(filename)
       endif
    else
       return
    endif

    allocate(file%in_pbuf(size(specifier)))
    file%in_pbuf(:) = .false.
    call trcdata_init( specifier, filename, filelist, datapath, fields, file, &
                       rmv_file, cycle_yr, fixed_ymd, fixed_tod, datatype)

    number_flds = 0
    if (associated(fields)) number_flds = size( fields )

    if( number_flds < 1 ) then
       has_aerodep_flx = .false.
       if (masterproc) then
          write(iulog,*) 'aerodep_flx_init: no aerosol deposition fluxes have been specified'
       endif
       return
    end if

    index_bulk_map(:) = -1
    index_modal_map(:) = -1

    do i = 1,number_flds

       ndx = get_ndx( fields(i)%fldnam, bulk_names )
       if (ndx >0) then
          index_bulk_map(ndx) = i
       else
          ndx = get_ndx( fields(i)%fldnam, modal_names )
          if (ndx >0) then
             index_modal_map(ndx) = i
          endif
       endif
       if (ndx>0) then
          call addfld(trim(fields(i)%fldnam)//'_D', horiz_only, 'A',fields(i)%units, 'prescribed aero dep' )
       else
          call endrun('aerodep_flx_init: aerosol flux name not recognized: '//trim(fields(i)%fldnam))
       endif
    enddo

    modal_fluxes = any(index_modal_map(:)>0)

    if (modal_fluxes) then

#if (defined DUSTTYPES)
#if (defined MODAL_AERO_10MODE)
       idx_bc1_dryis  = index_modal_map(1)
       idx_bc1_drycw  = index_modal_map(2)
       idx_pom1_dryis = index_modal_map(3)
       idx_pom1_drycw = index_modal_map(4)
       idx_soa1_dryis = index_modal_map(5)
       idx_soa1_drycw = index_modal_map(6)
       idx_soa2_dryis = index_modal_map(7)
       idx_soa2_drycw = index_modal_map(8)

       idx_dst15_dryis = index_modal_map(9)
       idx_dst15_drycw = index_modal_map(10)
       idx_dst25_dryis = index_modal_map(11)
       idx_dst25_drycw = index_modal_map(12)
       idx_dst35_dryis = index_modal_map(13)
       idx_dst35_drycw = index_modal_map(14)
       idx_dst45_dryis = index_modal_map(15)
       idx_dst45_drycw = index_modal_map(16)
       idx_dst55_dryis = index_modal_map(17)
       idx_dst55_drycw = index_modal_map(18)
       idx_dst65_dryis = index_modal_map(19)
       idx_dst65_drycw = index_modal_map(20)
       idx_dst75_dryis = index_modal_map(21)
       idx_dst75_drycw = index_modal_map(22)
       idx_dst85_dryis = index_modal_map(23)
       idx_dst85_drycw = index_modal_map(24)

       idx_dst17_dryis = index_modal_map(25)
       idx_dst17_drycw = index_modal_map(26)
       idx_dst27_dryis = index_modal_map(27)
       idx_dst27_drycw = index_modal_map(28)
       idx_dst37_dryis = index_modal_map(29)
       idx_dst37_drycw = index_modal_map(30)
       idx_dst47_dryis = index_modal_map(31)
       idx_dst47_drycw = index_modal_map(32)
       idx_dst57_dryis = index_modal_map(33)
       idx_dst57_drycw = index_modal_map(34)
       idx_dst67_dryis = index_modal_map(35)
       idx_dst67_drycw = index_modal_map(36)
       idx_dst77_dryis = index_modal_map(37)
       idx_dst77_drycw = index_modal_map(38)
       idx_dst87_dryis = index_modal_map(39)
       idx_dst87_drycw = index_modal_map(40)

       idx_dst18_dryis = index_modal_map(41)
       idx_dst18_drycw = index_modal_map(42)
       idx_dst28_dryis = index_modal_map(43)
       idx_dst28_drycw = index_modal_map(44)
       idx_dst38_dryis = index_modal_map(45)
       idx_dst38_drycw = index_modal_map(46)
       idx_dst48_dryis = index_modal_map(47)
       idx_dst48_drycw = index_modal_map(48)
       idx_dst58_dryis = index_modal_map(49)
       idx_dst58_drycw = index_modal_map(50)
       idx_dst68_dryis = index_modal_map(51)
       idx_dst68_drycw = index_modal_map(52)
       idx_dst78_dryis = index_modal_map(53)
       idx_dst78_drycw = index_modal_map(54)
       idx_dst88_dryis = index_modal_map(55)
       idx_dst88_drycw = index_modal_map(56)

       idx_dst19_dryis = index_modal_map(57)
       idx_dst19_drycw = index_modal_map(58)
       idx_dst29_dryis = index_modal_map(59)
       idx_dst29_drycw = index_modal_map(60)
       idx_dst39_dryis = index_modal_map(61)
       idx_dst39_drycw = index_modal_map(62)
       idx_dst49_dryis = index_modal_map(63)
       idx_dst49_drycw = index_modal_map(64)
       idx_dst59_dryis = index_modal_map(65)
       idx_dst59_drycw = index_modal_map(66)
       idx_dst69_dryis = index_modal_map(67)
       idx_dst69_drycw = index_modal_map(68)
       idx_dst79_dryis = index_modal_map(69)
       idx_dst79_drycw = index_modal_map(70)
       idx_dst89_dryis = index_modal_map(71)
       idx_dst89_drycw = index_modal_map(72)

       idx_dst110_dryis = index_modal_map(73)
       idx_dst110_drycw = index_modal_map(74)
       idx_dst210_dryis = index_modal_map(75)
       idx_dst210_drycw = index_modal_map(76)
       idx_dst310_dryis = index_modal_map(77)
       idx_dst310_drycw = index_modal_map(78)
       idx_dst410_dryis = index_modal_map(79)
       idx_dst410_drycw = index_modal_map(80)
       idx_dst510_dryis = index_modal_map(81)
       idx_dst510_drycw = index_modal_map(82)
       idx_dst610_dryis = index_modal_map(83)
       idx_dst610_drycw = index_modal_map(84)
       idx_dst710_dryis = index_modal_map(85)
       idx_dst710_drycw = index_modal_map(86)
       idx_dst810_dryis = index_modal_map(87)
       idx_dst810_drycw = index_modal_map(88)

       idx_bc1_wetis  = index_modal_map(89)
       idx_bc1_wetcw  = index_modal_map(90)
       idx_pom1_wetis = index_modal_map(91)
       idx_pom1_wetcw = index_modal_map(92)
       idx_soa1_wetis = index_modal_map(93)
       idx_soa1_wetcw = index_modal_map(94)

       idx_dst15_wetis = index_modal_map(95)
       idx_dst15_wetcw = index_modal_map(96)
       idx_dst25_wetis = index_modal_map(97)
       idx_dst25_wetcw = index_modal_map(98)
       idx_dst35_wetis = index_modal_map(99)
       idx_dst35_wetcw = index_modal_map(100)
       idx_dst45_wetis = index_modal_map(101)
       idx_dst45_wetcw = index_modal_map(102)
       idx_dst55_wetis = index_modal_map(103)
       idx_dst55_wetcw = index_modal_map(104)
       idx_dst65_wetis = index_modal_map(105)
       idx_dst65_wetcw = index_modal_map(106)
       idx_dst75_wetis = index_modal_map(107)
       idx_dst75_wetcw = index_modal_map(108)
       idx_dst85_wetis = index_modal_map(109)
       idx_dst85_wetcw = index_modal_map(110)

       idx_dst17_wetis = index_modal_map(111)
       idx_dst17_wetcw = index_modal_map(112)
       idx_dst27_wetis = index_modal_map(113)
       idx_dst27_wetcw = index_modal_map(114)
       idx_dst37_wetis = index_modal_map(115)
       idx_dst37_wetcw = index_modal_map(116)
       idx_dst47_wetis = index_modal_map(117)
       idx_dst47_wetcw = index_modal_map(118)
       idx_dst57_wetis = index_modal_map(119)
       idx_dst57_wetcw = index_modal_map(120)
       idx_dst67_wetis = index_modal_map(121)
       idx_dst67_wetcw = index_modal_map(122)
       idx_dst77_wetis = index_modal_map(123)
       idx_dst77_wetcw = index_modal_map(124)
       idx_dst87_wetis = index_modal_map(125)
       idx_dst87_wetcw = index_modal_map(126)

       idx_dst18_wetis = index_modal_map(127)
       idx_dst18_wetcw = index_modal_map(128)
       idx_dst28_wetis = index_modal_map(129)
       idx_dst28_wetcw = index_modal_map(130)
       idx_dst38_wetis = index_modal_map(131)
       idx_dst38_wetcw = index_modal_map(132)
       idx_dst48_wetis = index_modal_map(133)
       idx_dst48_wetcw = index_modal_map(134)
       idx_dst58_wetis = index_modal_map(135)
       idx_dst58_wetcw = index_modal_map(136)
       idx_dst68_wetis = index_modal_map(137)
       idx_dst68_wetcw = index_modal_map(138)
       idx_dst78_wetis = index_modal_map(139)
       idx_dst78_wetcw = index_modal_map(140)
       idx_dst88_wetis = index_modal_map(141)
       idx_dst88_wetcw = index_modal_map(142)

       idx_dst19_wetis = index_modal_map(143)
       idx_dst19_wetcw = index_modal_map(144)
       idx_dst29_wetis = index_modal_map(145)
       idx_dst29_wetcw = index_modal_map(146)
       idx_dst39_wetis = index_modal_map(147)
       idx_dst39_wetcw = index_modal_map(148)
       idx_dst49_wetis = index_modal_map(149)
       idx_dst49_wetcw = index_modal_map(150)
       idx_dst59_wetis = index_modal_map(151)
       idx_dst59_wetcw = index_modal_map(152)
       idx_dst69_wetis = index_modal_map(153)
       idx_dst69_wetcw = index_modal_map(154)
       idx_dst79_wetis = index_modal_map(155)
       idx_dst79_wetcw = index_modal_map(156)
       idx_dst89_wetis = index_modal_map(157)
       idx_dst89_wetcw = index_modal_map(158)

       idx_dst110_wetis = index_modal_map(159)
       idx_dst110_wetcw = index_modal_map(160)
       idx_dst210_wetis = index_modal_map(161)
       idx_dst210_wetcw = index_modal_map(162)
       idx_dst310_wetis = index_modal_map(163)
       idx_dst310_wetcw = index_modal_map(164)
       idx_dst410_wetis = index_modal_map(165)
       idx_dst410_wetcw = index_modal_map(166)
       idx_dst510_wetis = index_modal_map(167)
       idx_dst510_wetcw = index_modal_map(168)
       idx_dst610_wetis = index_modal_map(169)
       idx_dst610_wetcw = index_modal_map(170)
       idx_dst710_wetis = index_modal_map(171)
       idx_dst710_wetcw = index_modal_map(172)
       idx_dst810_wetis = index_modal_map(173)
       idx_dst810_wetcw = index_modal_map(174)

       call modal_aero_deposition_init( &
                 bcphi_indices=(/idx_bc1/), &
                 ocphi_indices=(/idx_pom1,idx_soa1/), &
                 ocpho_indices=(/idx_soa2/), &
                 fine_dust_indices=(/idx_dst15,idx_dst25,idx_dst35,idx_dst45,idx_dst55,idx_dst65,idx_dst75,idx_dst85/),&
                 crse_dust_indices=(/idx_dst18,idx_dst28,idx_dst38,idx_dst48,idx_dst58,idx_dst68,idx_dst78,idx_dst88/) )
                 ! MAM9 here I may need to add dst7 and dst9; double check

       ! Longlei Li MAM10


       write(iulog,*)'  Longlei debugging (aerodep_flx.F90: test01)'

#else 
       idx_bc1_dryis  = index_modal_map(1)
       idx_bc1_drycw  = index_modal_map(2)
       idx_pom1_dryis = index_modal_map(3)
       idx_pom1_drycw = index_modal_map(4)
       idx_soa1_dryis = index_modal_map(5)
       idx_soa1_drycw = index_modal_map(6)
       idx_soa2_dryis = index_modal_map(7)
       idx_soa2_drycw = index_modal_map(8)
       idx_dst1_dryis = index_modal_map(9)
       idx_dst1_drycw = index_modal_map(10)
       idx_dst3_dryis = index_modal_map(11)
       idx_dst3_drycw = index_modal_map(12)

       idx_bc1_wetis  = index_modal_map(13)
       idx_bc1_wetcw  = index_modal_map(14)
       idx_pom1_wetis = index_modal_map(15)
       idx_pom1_wetcw = index_modal_map(16)
       idx_soa1_wetis = index_modal_map(17)
       idx_soa1_wetcw = index_modal_map(18)
       idx_dst1_wetis = index_modal_map(19)
       idx_dst1_wetcw = index_modal_map(20)
       idx_dst3_wetis = index_modal_map(21)
       idx_dst3_wetcw = index_modal_map(22)

       call modal_aero_deposition_init( bcphi_indices=(/idx_bc1/), &
                                        ocphi_indices=(/idx_pom1,idx_soa1/), &
                                        ocpho_indices=(/idx_soa2/), &
                                        fine_dust_indices=(/idx_dst1/),&
                                        crse_dust_indices=(/idx_dst3/) )

#endif
#endif
    else

       ibcphiwet = index_bulk_map(1)
       ibcphodry = index_bulk_map(2)
       ibcphidry = index_bulk_map(3)
       iocphiwet = index_bulk_map(4)
       iocphodry = index_bulk_map(5)
       iocphidry = index_bulk_map(6)
       idstdry1  = index_bulk_map(7)
       idstdry2  = index_bulk_map(8)
       idstdry3  = index_bulk_map(9)
       idstdry4  = index_bulk_map(10)
       idstwet1  = index_bulk_map(11)
       idstwet2  = index_bulk_map(12)
       idstwet3  = index_bulk_map(13)
       idstwet4  = index_bulk_map(14)

    endif

  end subroutine aerodep_flx_init
!-------------------------------------------------------------------
! sets namelist options
!-------------------------------------------------------------------
subroutine aerodep_flx_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'aerodep_flx_readnl'

   character(len=32)  :: aerodep_flx_specifier(N_MODAL)
   character(len=256) :: aerodep_flx_file
   character(len=256) :: aerodep_flx_filelist
   character(len=256) :: aerodep_flx_datapath
   character(len=32)  :: aerodep_flx_type
   logical            :: aerodep_flx_rmfile
   integer            :: aerodep_flx_cycle_yr
   integer            :: aerodep_flx_fixed_ymd
   integer            :: aerodep_flx_fixed_tod

   namelist /aerodep_flx_nl/ &
      aerodep_flx_specifier, &
      aerodep_flx_file,      &
      aerodep_flx_filelist,  &
      aerodep_flx_datapath,  &
      aerodep_flx_type,      &
      aerodep_flx_rmfile,    &
      aerodep_flx_cycle_yr,  &
      aerodep_flx_fixed_ymd, &
      aerodep_flx_fixed_tod      
   !-----------------------------------------------------------------------------

   ! Initialize namelist variables from local module variables.
   aerodep_flx_specifier= specifier
   aerodep_flx_file     = filename
   aerodep_flx_filelist = filelist
   aerodep_flx_datapath = datapath
   aerodep_flx_type     = datatype
   aerodep_flx_rmfile   = rmv_file
   aerodep_flx_cycle_yr = cycle_yr
   aerodep_flx_fixed_ymd= fixed_ymd
   aerodep_flx_fixed_tod= fixed_tod

   ! Read namelist
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'aerodep_flx_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, aerodep_flx_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(aerodep_flx_specifier,len(aerodep_flx_specifier(1))*N_MODAL,     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_file,     len(aerodep_flx_file),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_filelist, len(aerodep_flx_filelist), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_datapath, len(aerodep_flx_datapath), mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_type,     len(aerodep_flx_type),     mpichar, 0, mpicom)
   call mpibcast(aerodep_flx_rmfile,   1, mpilog,  0, mpicom)
   call mpibcast(aerodep_flx_cycle_yr, 1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_ymd,1, mpiint,  0, mpicom)
   call mpibcast(aerodep_flx_fixed_tod,1, mpiint,  0, mpicom)
#endif

   ! Update module variables with user settings.
   specifier  = aerodep_flx_specifier
   filename   = aerodep_flx_file
   filelist   = aerodep_flx_filelist
   datapath   = aerodep_flx_datapath
   datatype   = aerodep_flx_type
   rmv_file   = aerodep_flx_rmfile
   cycle_yr   = aerodep_flx_cycle_yr
   fixed_ymd  = aerodep_flx_fixed_ymd
   fixed_tod  = aerodep_flx_fixed_tod

   ! Turn on prescribed volcanics if user has specified an input dataset.
   if (len_trim(filename) > 0 .and. filename.ne.'NONE' ) has_aerodep_flx = .true.

end subroutine aerodep_flx_readnl

!-------------------------------------------------------------------
! sets the aerosol deposition fluxes in the cam_out structure 
! to be sent to the surface models
!-------------------------------------------------------------------
  subroutine aerodep_flx_set( cam_out, ncol, lchnk )
    use camsrfexch,       only : cam_out_t     

    type(cam_out_t),     intent(inout) :: cam_out
    integer,             intent(in)    :: ncol, lchnk
    
    if( .not. has_aerodep_flx ) return
    
    if (modal_fluxes) then
       call set_modal_fluxes( cam_out, ncol, lchnk )
    else
       call set_bulk_fluxes( cam_out, ncol, lchnk )
    endif

  end subroutine aerodep_flx_set

!-------------------------------------------------------------------
! advances the prescribed fluxes to the current time step
!-------------------------------------------------------------------
  subroutine aerodep_flx_adv( state, pbuf2d, cam_out )

    use tracer_data,      only : advance_trcdata
    use physics_types,    only : physics_state
    use camsrfexch,       only : cam_out_t
    use physics_buffer, only : physics_buffer_desc

    implicit none

    type(physics_state), intent(in)    :: state(begchunk:endchunk)                 
    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    integer :: c, ncol
    
    if( .not. has_aerodep_flx ) return

    call advance_trcdata( fields, file, state, pbuf2d  )

!$OMP PARALLEL DO PRIVATE (C, NCOL)
    do c = begchunk, endchunk
       ncol = state(c)%ncol
       call aerodep_flx_set( cam_out(c), ncol, c )
    enddo

  end subroutine aerodep_flx_adv

!-------------------------------------------------------------------
! returns true if aerosol dep fluxes are prescribed from dataset
!-------------------------------------------------------------------
  function aerodep_flx_prescribed()
    logical :: aerodep_flx_prescribed
    aerodep_flx_prescribed = has_aerodep_flx
  endfunction aerodep_flx_prescribed

! private methods
!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_bulk_fluxes( cam_out, ncol, lchnk )
    use camsrfexch,            only : cam_out_t     

    ! Arguments
    type(cam_out_t), intent(inout) :: cam_out
    integer,         intent(in)    :: ncol, lchnk

    call set_fluxes( cam_out%bcphiwet, ibcphiwet, ncol, lchnk )
    call set_fluxes( cam_out%bcphidry, ibcphidry, ncol, lchnk )
    call set_fluxes( cam_out%bcphodry, ibcphodry, ncol, lchnk )

    call set_fluxes( cam_out%ocphiwet, iocphiwet, ncol, lchnk )
    call set_fluxes( cam_out%ocphidry, iocphidry, ncol, lchnk )
    call set_fluxes( cam_out%ocphodry, iocphodry, ncol, lchnk )

    call set_fluxes( cam_out%dstdry1, idstdry1, ncol, lchnk )
    call set_fluxes( cam_out%dstdry2, idstdry2, ncol, lchnk )
    call set_fluxes( cam_out%dstdry3, idstdry3, ncol, lchnk )
    call set_fluxes( cam_out%dstdry4, idstdry4, ncol, lchnk )

    call set_fluxes( cam_out%dstwet1, idstwet1, ncol, lchnk )
    call set_fluxes( cam_out%dstwet2, idstwet2, ncol, lchnk )
    call set_fluxes( cam_out%dstwet3, idstwet3, ncol, lchnk )
    call set_fluxes( cam_out%dstwet4, idstwet4, ncol, lchnk )

  end subroutine set_bulk_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_modal_fluxes( cam_out, ncol, lchnk )
    use camsrfexch,            only : cam_out_t     
    use modal_aero_deposition, only : set_srf_drydep, set_srf_wetdep

    ! Arguments
    type(cam_out_t), intent(inout) :: cam_out
    integer,         intent(in)    :: ncol, lchnk

    ! local vars
    integer :: i
    real(r8) :: aerdepdryis(pcols,nmodal_idxs)
    real(r8) :: aerdepdrycw(pcols,nmodal_idxs)
    real(r8) :: aerdepwetis(pcols,nmodal_idxs)
    real(r8) :: aerdepwetcw(pcols,nmodal_idxs)

    ! bin the fluxes as using modal_aero_deposition...

    aerdepdryis(:,:) = 0._r8
    aerdepdrycw(:,:) = 0._r8
    aerdepwetis(:,:) = 0._r8
    aerdepwetcw(:,:) = 0._r8

    call set_fluxes( aerdepwetis(:ncol,idx_bc1 ), idx_bc1_wetis , ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_bc1 ), idx_bc1_wetcw , ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_pom1), idx_pom1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_pom1), idx_pom1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_soa1), idx_soa1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_soa1), idx_soa1_wetcw, ncol, lchnk )
#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
    call set_fluxes( aerdepwetis(:ncol,idx_dst15), idx_dst15_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst15), idx_dst15_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst25), idx_dst25_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst25), idx_dst25_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst35), idx_dst35_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst35), idx_dst35_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst45), idx_dst45_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst45), idx_dst45_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst55), idx_dst55_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst55), idx_dst55_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst65), idx_dst65_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst65), idx_dst65_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst75), idx_dst75_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst75), idx_dst75_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst85), idx_dst85_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst85), idx_dst85_wetcw, ncol, lchnk )

    call set_fluxes( aerdepwetis(:ncol,idx_dst17), idx_dst17_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst17), idx_dst17_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst27), idx_dst27_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst27), idx_dst27_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst37), idx_dst37_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst37), idx_dst37_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst47), idx_dst47_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst47), idx_dst47_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst57), idx_dst57_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst57), idx_dst57_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst67), idx_dst67_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst67), idx_dst67_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst77), idx_dst77_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst77), idx_dst77_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst87), idx_dst87_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst87), idx_dst87_wetcw, ncol, lchnk )

    call set_fluxes( aerdepwetis(:ncol,idx_dst18), idx_dst18_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst18), idx_dst18_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst28), idx_dst28_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst28), idx_dst28_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst38), idx_dst38_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst38), idx_dst38_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst48), idx_dst48_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst48), idx_dst48_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst58), idx_dst58_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst58), idx_dst58_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst68), idx_dst68_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst68), idx_dst68_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst78), idx_dst78_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst78), idx_dst78_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst88), idx_dst88_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst88), idx_dst88_wetcw, ncol, lchnk )

    call set_fluxes( aerdepwetis(:ncol,idx_dst19), idx_dst19_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst19), idx_dst19_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst29), idx_dst29_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst29), idx_dst29_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst39), idx_dst39_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst39), idx_dst39_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst49), idx_dst49_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst49), idx_dst49_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst59), idx_dst59_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst59), idx_dst59_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst69), idx_dst69_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst69), idx_dst69_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst79), idx_dst79_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst79), idx_dst79_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst89), idx_dst89_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst89), idx_dst89_wetcw, ncol, lchnk )

    call set_fluxes( aerdepwetis(:ncol,idx_dst110), idx_dst110_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst110), idx_dst110_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst210), idx_dst210_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst210), idx_dst210_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst310), idx_dst310_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst310), idx_dst310_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst410), idx_dst410_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst410), idx_dst410_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst510), idx_dst510_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst510), idx_dst510_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst610), idx_dst610_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst610), idx_dst610_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst710), idx_dst710_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst710), idx_dst710_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst810), idx_dst810_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst810), idx_dst810_wetcw, ncol, lchnk )
#endif
#else   
    call set_fluxes( aerdepwetis(:ncol,idx_dst1), idx_dst1_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst1), idx_dst1_wetcw, ncol, lchnk )
    call set_fluxes( aerdepwetis(:ncol,idx_dst3), idx_dst3_wetis, ncol, lchnk )
    call set_fluxes( aerdepwetcw(:ncol,idx_dst3), idx_dst3_wetcw, ncol, lchnk )
#endif

    call set_fluxes( aerdepdryis(:ncol,idx_bc1 ), idx_bc1_dryis , ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_bc1 ), idx_bc1_drycw , ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_pom1), idx_pom1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_pom1), idx_pom1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_soa1), idx_soa1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_soa1), idx_soa1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_soa2), idx_soa2_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_soa2), idx_soa2_drycw, ncol, lchnk )
#if ( defined DUSTTYPES )
#if (defined MODAL_AERO_10MODE)
    call set_fluxes( aerdepdryis(:ncol,idx_dst15), idx_dst15_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst15), idx_dst15_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst25), idx_dst25_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst25), idx_dst25_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst35), idx_dst35_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst35), idx_dst35_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst45), idx_dst45_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst45), idx_dst45_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst55), idx_dst55_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst55), idx_dst55_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst65), idx_dst65_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst65), idx_dst65_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst75), idx_dst75_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst75), idx_dst75_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst85), idx_dst85_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst85), idx_dst85_drycw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_dst17), idx_dst17_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst17), idx_dst17_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst27), idx_dst27_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst27), idx_dst27_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst37), idx_dst37_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst37), idx_dst37_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst47), idx_dst47_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst47), idx_dst47_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst57), idx_dst57_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst57), idx_dst57_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst67), idx_dst67_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst67), idx_dst67_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst77), idx_dst77_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst77), idx_dst77_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst87), idx_dst87_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst87), idx_dst87_drycw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_dst18), idx_dst18_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst18), idx_dst18_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst28), idx_dst28_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst28), idx_dst28_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst38), idx_dst38_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst38), idx_dst38_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst48), idx_dst48_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst48), idx_dst48_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst58), idx_dst58_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst58), idx_dst58_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst68), idx_dst68_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst68), idx_dst68_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst78), idx_dst78_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst78), idx_dst78_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst88), idx_dst88_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst88), idx_dst88_drycw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_dst19), idx_dst19_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst19), idx_dst19_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst29), idx_dst29_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst29), idx_dst29_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst39), idx_dst39_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst39), idx_dst39_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst49), idx_dst49_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst49), idx_dst49_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst59), idx_dst59_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst59), idx_dst59_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst69), idx_dst69_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst69), idx_dst69_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst79), idx_dst79_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst79), idx_dst79_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst89), idx_dst89_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst89), idx_dst89_drycw, ncol, lchnk )

    call set_fluxes( aerdepdryis(:ncol,idx_dst110), idx_dst110_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst110), idx_dst110_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst210), idx_dst210_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst210), idx_dst210_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst310), idx_dst310_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst310), idx_dst310_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst410), idx_dst410_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst410), idx_dst410_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst510), idx_dst510_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst510), idx_dst510_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst610), idx_dst610_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst610), idx_dst610_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst710), idx_dst710_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst710), idx_dst710_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst810), idx_dst810_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst810), idx_dst810_drycw, ncol, lchnk )
#endif
#else   
    call set_fluxes( aerdepdryis(:ncol,idx_dst1), idx_dst1_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst1), idx_dst1_drycw, ncol, lchnk )
    call set_fluxes( aerdepdryis(:ncol,idx_dst3), idx_dst3_dryis, ncol, lchnk )
    call set_fluxes( aerdepdrycw(:ncol,idx_dst3), idx_dst3_drycw, ncol, lchnk )
#endif

    call set_srf_drydep(aerdepdryis, aerdepdrycw, cam_out)
    call set_srf_wetdep(aerdepwetis, aerdepwetcw, cam_out)

  end  subroutine set_modal_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  subroutine set_fluxes( fluxes, fld_indx, ncol, lchnk )
    use cam_history,  only : outfld

    real(r8), intent(inout) :: fluxes(:)
    integer,  intent(in)    :: fld_indx, ncol, lchnk

    integer :: i

    if (fld_indx<1) return

    do i = 1,ncol
       ! modal aero wet dep history fields are negative
       fluxes(i) = fields(fld_indx)%data(i,1,lchnk)
    enddo

    call outfld(trim(fields(fld_indx)%fldnam)//'_D', fluxes(:ncol), ncol, lchnk )

  endsubroutine set_fluxes

!-------------------------------------------------------------------
!-------------------------------------------------------------------
  integer function get_ndx( name, list )

    implicit none
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: list(:)

    integer :: i
    integer :: maxnum

    maxnum = size(list)

    get_ndx = -1
    do i = 1, maxnum
      if ( trim(name) == trim(list(i)) ) then
        get_ndx = i
        return
      endif
    enddo

  end function get_ndx

end module aerodep_flx
