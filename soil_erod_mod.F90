!===============================================================================
!===============================================================================
module soil_erod_mod
  use shr_kind_mod,     only: r8 => shr_kind_r8, cl => shr_kind_cl
  use cam_logfile,      only: iulog
  use spmd_utils,       only: masterproc
  use cam_abortutils,   only: endrun

  implicit none
  private

  public :: soil_erod_init
  public :: soil_erodibility
  public :: soil_erod_fact

#if ( defined DUSTTYPES )
  ! Longlei Li - BRIFT already applied offline
  public ::  illite_frac_accum
  public ::  illite_frac_aitken
  public ::  illite_frac_coarse01
  public ::  illite_frac_coarse02
  public ::  illite_frac_giant
  public ::  kaolinite_frac_accum
  public ::  kaolinite_frac_aitken
  public ::  kaolinite_frac_coarse01
  public ::  kaolinite_frac_coarse02
  public ::  kaolinite_frac_giant
  public ::  mont_frac_accum
  public ::  mont_frac_aitken
  public ::  mont_frac_coarse01
  public ::  mont_frac_coarse02
  public ::  mont_frac_giant
  public ::  hematite_frac_accum
  public ::  hematite_frac_aitken
  public ::  hematite_frac_coarse01
  public ::  hematite_frac_coarse02
  public ::  hematite_frac_giant
  public ::  quartz_frac_accum
  public ::  quartz_frac_aitken
  public ::  quartz_frac_coarse01
  public ::  quartz_frac_coarse02
  public ::  quartz_frac_giant
  public ::  calcite_frac_accum
  public ::  calcite_frac_aitken
  public ::  calcite_frac_coarse01
  public ::  calcite_frac_coarse02
  public ::  calcite_frac_giant
  public ::  feldspar_frac_accum
  public ::  feldspar_frac_aitken
  public ::  feldspar_frac_coarse01
  public ::  feldspar_frac_coarse02
  public ::  feldspar_frac_giant
  public ::  gypsum_frac_accum
  public ::  gypsum_frac_aitken
  public ::  gypsum_frac_coarse01
  public ::  gypsum_frac_coarse02
  public ::  gypsum_frac_giant

#endif
  public :: asphericity_factor

  real(r8), allocatable ::  soil_erodibility(:,:)  ! soil erodibility factor
#if ( defined DUSTTYPES )
  ! Longlei Li - BRIFT already applied offline
  real(r8), allocatable ::  illite_frac_accum(:,:)
  real(r8), allocatable ::  illite_frac_aitken(:,:)
  real(r8), allocatable ::  illite_frac_coarse01(:,:)
  real(r8), allocatable ::  illite_frac_coarse02(:,:)
  real(r8), allocatable ::  illite_frac_giant(:,:)
  real(r8), allocatable ::  kaolinite_frac_accum(:,:)
  real(r8), allocatable ::  kaolinite_frac_aitken(:,:)
  real(r8), allocatable ::  kaolinite_frac_coarse01(:,:)
  real(r8), allocatable ::  kaolinite_frac_coarse02(:,:)
  real(r8), allocatable ::  kaolinite_frac_giant(:,:)
  real(r8), allocatable ::  mont_frac_accum(:,:)
  real(r8), allocatable ::  mont_frac_aitken(:,:)
  real(r8), allocatable ::  mont_frac_coarse01(:,:)
  real(r8), allocatable ::  mont_frac_coarse02(:,:)
  real(r8), allocatable ::  mont_frac_giant(:,:)
  real(r8), allocatable ::  hematite_frac_accum(:,:)
  real(r8), allocatable ::  hematite_frac_aitken(:,:)
  real(r8), allocatable ::  hematite_frac_coarse01(:,:)
  real(r8), allocatable ::  hematite_frac_coarse02(:,:)
  real(r8), allocatable ::  hematite_frac_giant(:,:)
  real(r8), allocatable ::  quartz_frac_accum(:,:)
  real(r8), allocatable ::  quartz_frac_aitken(:,:)
  real(r8), allocatable ::  quartz_frac_coarse01(:,:)
  real(r8), allocatable ::  quartz_frac_coarse02(:,:)
  real(r8), allocatable ::  quartz_frac_giant(:,:)
  real(r8), allocatable ::  calcite_frac_accum(:,:)
  real(r8), allocatable ::  calcite_frac_aitken(:,:)
  real(r8), allocatable ::  calcite_frac_coarse01(:,:)
  real(r8), allocatable ::  calcite_frac_coarse02(:,:)
  real(r8), allocatable ::  calcite_frac_giant(:,:)
  real(r8), allocatable ::  feldspar_frac_accum(:,:)
  real(r8), allocatable ::  feldspar_frac_aitken(:,:)
  real(r8), allocatable ::  feldspar_frac_coarse01(:,:)
  real(r8), allocatable ::  feldspar_frac_coarse02(:,:)
  real(r8), allocatable ::  feldspar_frac_giant(:,:)
  real(r8), allocatable ::  gypsum_frac_accum(:,:)
  real(r8), allocatable ::  gypsum_frac_aitken(:,:)
  real(r8), allocatable ::  gypsum_frac_coarse01(:,:)
  real(r8), allocatable ::  gypsum_frac_coarse02(:,:)
  real(r8), allocatable ::  gypsum_frac_giant(:,:)
#endif
!++Longlei Li
 real(r8), allocatable ::  asphericity_factor(:)
!--Longlei Li

  real(r8) :: soil_erod_fact                       ! tuning parameter for dust emissions

contains

  !=============================================================================
  !=============================================================================
  subroutine soil_erod_init( dust_emis_fact, soil_erod_file )
    use interpolate_data, only: lininterp_init, lininterp, lininterp_finish, interp_type
    use ppgrid,           only: begchunk, endchunk, pcols
    use mo_constants,     only: pi, d2r
    use pio,              only: file_desc_t,pio_inq_dimid,pio_inq_dimlen,pio_get_var,pio_inq_varid, PIO_NOWRITE
    use phys_grid,        only: get_ncols_p, get_rlat_all_p, get_rlon_all_p
    use cam_pio_utils,    only: cam_pio_openfile
    use ioFileMod,        only: getfil

    real(r8),         intent(in) :: dust_emis_fact
    character(len=*), intent(in) :: soil_erod_file

    real(r8), allocatable ::  soil_erodibility_in(:,:)  ! temporary input array
    real(r8), allocatable :: dst_lons(:)
    real(r8), allocatable :: dst_lats(:)
#if ( defined DUSTTYPES )
  ! Longlei Li - BRIFT already applied offline
    real(r8), allocatable :: illite_frac_in_accum(:,:)
    real(r8), allocatable :: illite_frac_in_aitken(:,:)
    real(r8), allocatable :: illite_frac_in_coarse01(:,:)
    real(r8), allocatable :: illite_frac_in_coarse02(:,:)
    real(r8), allocatable :: illite_frac_in_giant(:,:)
    real(r8), allocatable :: kaolinite_frac_in_accum(:,:)
    real(r8), allocatable :: kaolinite_frac_in_aitken(:,:)
    real(r8), allocatable :: kaolinite_frac_in_coarse01(:,:)
    real(r8), allocatable :: kaolinite_frac_in_coarse02(:,:)
    real(r8), allocatable :: kaolinite_frac_in_giant(:,:)
    real(r8), allocatable :: mont_frac_in_accum(:,:)
    real(r8), allocatable :: mont_frac_in_aitken(:,:)
    real(r8), allocatable :: mont_frac_in_coarse01(:,:)
    real(r8), allocatable :: mont_frac_in_coarse02(:,:)
    real(r8), allocatable :: mont_frac_in_giant(:,:)
    real(r8), allocatable :: hematite_frac_in_accum(:,:)
    real(r8), allocatable :: hematite_frac_in_aitken(:,:)
    real(r8), allocatable :: hematite_frac_in_coarse01(:,:)
    real(r8), allocatable :: hematite_frac_in_coarse02(:,:)
    real(r8), allocatable :: hematite_frac_in_giant(:,:)
    real(r8), allocatable :: quartz_frac_in_accum(:,:)
    real(r8), allocatable :: quartz_frac_in_aitken(:,:)
    real(r8), allocatable :: quartz_frac_in_coarse01(:,:)
    real(r8), allocatable :: quartz_frac_in_coarse02(:,:)
    real(r8), allocatable :: quartz_frac_in_giant(:,:)
    real(r8), allocatable :: calcite_frac_in_accum(:,:)
    real(r8), allocatable :: calcite_frac_in_aitken(:,:)
    real(r8), allocatable :: calcite_frac_in_coarse01(:,:)
    real(r8), allocatable :: calcite_frac_in_coarse02(:,:)
    real(r8), allocatable :: calcite_frac_in_giant(:,:)
    real(r8), allocatable :: feldspar_frac_in_accum(:,:)
    real(r8), allocatable :: feldspar_frac_in_aitken(:,:)
    real(r8), allocatable :: feldspar_frac_in_coarse01(:,:)
    real(r8), allocatable :: feldspar_frac_in_coarse02(:,:)
    real(r8), allocatable :: feldspar_frac_in_giant(:,:)
    real(r8), allocatable :: gypsum_frac_in_accum(:,:)
    real(r8), allocatable :: gypsum_frac_in_aitken(:,:)
    real(r8), allocatable :: gypsum_frac_in_coarse01(:,:)
    real(r8), allocatable :: gypsum_frac_in_coarse02(:,:)
    real(r8), allocatable :: gypsum_frac_in_giant(:,:)
#endif
!--++Longlei--
    real(r8), allocatable ::  correct_factor(:,:)
!----Longlei--

    character(len=cl)     :: infile
    integer :: did, vid, nlat, nlon
    type(file_desc_t) :: ncid

    type(interp_type) :: lon_wgts, lat_wgts
    real(r8) :: to_lats(pcols), to_lons(pcols)
    integer :: c, ncols, ierr
    real(r8), parameter :: zero=0._r8, twopi=2._r8*pi

    soil_erod_fact = dust_emis_fact

    ! Summary to log file
    if (masterproc) then
       write(iulog,*) 'soil_erod_mod: soil erodibility dataset: ', trim(soil_erod_file)
       write(iulog,*) 'soil_erod_mod: soil_erod_fact = ', soil_erod_fact
    end if

    ! for soil erodibility in mobilization, apply inside CAM instead of lsm.
    ! read in soil erodibility factors, similar to Zender's boundary conditions

    ! Get file name.  
    call getfil(soil_erod_file, infile, 0)
    call cam_pio_openfile (ncid, trim(infile), PIO_NOWRITE)

    ! Get input data resolution.
    ierr = pio_inq_dimid( ncid, 'lon', did )
    ierr = pio_inq_dimlen( ncid, did, nlon )

    ierr = pio_inq_dimid( ncid, 'lat', did )
    ierr = pio_inq_dimlen( ncid, did, nlat )

    allocate(dst_lons(nlon))
    allocate(dst_lats(nlat))
    allocate(soil_erodibility_in(nlon,nlat))

#if ( defined DUSTTYPES )
  ! Longlei Li - BRIFT already applied offline
    allocate(illite_frac_in_accum(nlon,nlat))
    allocate(illite_frac_in_aitken(nlon,nlat))
    allocate(illite_frac_in_coarse01(nlon,nlat))
    allocate(illite_frac_in_coarse02(nlon,nlat))
    allocate(illite_frac_in_giant(nlon,nlat))

    allocate(kaolinite_frac_in_accum(nlon,nlat))
    allocate(kaolinite_frac_in_aitken(nlon,nlat))
    allocate(kaolinite_frac_in_coarse01(nlon,nlat))
    allocate(kaolinite_frac_in_coarse02(nlon,nlat))
    allocate(kaolinite_frac_in_giant(nlon,nlat))

    allocate(mont_frac_in_accum(nlon,nlat))
    allocate(mont_frac_in_aitken(nlon,nlat))
    allocate(mont_frac_in_coarse01(nlon,nlat))
    allocate(mont_frac_in_coarse02(nlon,nlat))
    allocate(mont_frac_in_giant(nlon,nlat))

    allocate(hematite_frac_in_accum(nlon,nlat))
    allocate(hematite_frac_in_aitken(nlon,nlat))
    allocate(hematite_frac_in_coarse01(nlon,nlat))
    allocate(hematite_frac_in_coarse02(nlon,nlat))
    allocate(hematite_frac_in_giant(nlon,nlat))

    allocate(quartz_frac_in_accum(nlon,nlat))
    allocate(quartz_frac_in_aitken(nlon,nlat))
    allocate(quartz_frac_in_coarse01(nlon,nlat))
    allocate(quartz_frac_in_coarse02(nlon,nlat))
    allocate(quartz_frac_in_giant(nlon,nlat))

    allocate(calcite_frac_in_accum(nlon,nlat))
    allocate(calcite_frac_in_aitken(nlon,nlat))
    allocate(calcite_frac_in_coarse01(nlon,nlat))
    allocate(calcite_frac_in_coarse02(nlon,nlat))
    allocate(calcite_frac_in_giant(nlon,nlat))

    allocate(feldspar_frac_in_accum(nlon,nlat))
    allocate(feldspar_frac_in_aitken(nlon,nlat))
    allocate(feldspar_frac_in_coarse01(nlon,nlat))
    allocate(feldspar_frac_in_coarse02(nlon,nlat))
    allocate(feldspar_frac_in_giant(nlon,nlat))

    allocate(gypsum_frac_in_accum(nlon,nlat))
    allocate(gypsum_frac_in_aitken(nlon,nlat))
    allocate(gypsum_frac_in_coarse01(nlon,nlat))
    allocate(gypsum_frac_in_coarse02(nlon,nlat))
    allocate(gypsum_frac_in_giant(nlon,nlat))

#endif
    allocate(correct_factor(nlon,nlat))


    ierr = pio_inq_varid( ncid, 'lon', vid )
    ierr = pio_get_var( ncid, vid, dst_lons  )

    ierr = pio_inq_varid( ncid, 'lat', vid )
    ierr = pio_get_var( ncid, vid, dst_lats  )

    ierr = pio_inq_varid( ncid, 'mbl_bsn_fct_geo', vid )
    ierr = pio_get_var( ncid, vid, soil_erodibility_in )

#if ( defined DUSTTYPES )
! mineralogy...
  ! Longlei Li - BRIFT already applied offline
    ierr = pio_inq_varid( ncid, 'FracI1', vid )            !FracI1
    ierr = pio_get_var( ncid, vid, illite_frac_in_aitken )
    ierr = pio_inq_varid( ncid, 'FracI2', vid )            !FracI2
    ierr = pio_get_var( ncid, vid, illite_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracI3', vid )            !FracI3
    ierr = pio_get_var( ncid, vid, illite_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracI4', vid )            !FracI4
    ierr = pio_get_var( ncid, vid, illite_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracI5', vid )            !FracI5
    ierr = pio_get_var( ncid, vid, illite_frac_in_giant )

    ierr = pio_inq_varid( ncid, 'FracK1', vid )            !FracK1
    ierr = pio_get_var( ncid, vid, kaolinite_frac_in_aitken )    
    ierr = pio_inq_varid( ncid, 'FracK2', vid )            !FracK2
    ierr = pio_get_var( ncid, vid, kaolinite_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracK3', vid )            !FracK3
    ierr = pio_get_var( ncid, vid, kaolinite_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracK4', vid )            !FracK4
    ierr = pio_get_var( ncid, vid, kaolinite_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracK5', vid )            !FracK5
    ierr = pio_get_var( ncid, vid, kaolinite_frac_in_giant )
        
    ierr = pio_inq_varid( ncid, 'FracS1', vid )            !FracS1
    ierr = pio_get_var( ncid, vid, mont_frac_in_aitken ) 
    ierr = pio_inq_varid( ncid, 'FracS2', vid )            !FracS2
    ierr = pio_get_var( ncid, vid, mont_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracS3', vid )            !FracS3
    ierr = pio_get_var( ncid, vid, mont_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracS4', vid )            !FracS4
    ierr = pio_get_var( ncid, vid, mont_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracS5', vid )            !FracS5
    ierr = pio_get_var( ncid, vid, mont_frac_in_giant )
        
    ierr = pio_inq_varid( ncid, 'FracH1', vid )            !FracH1
    ierr = pio_get_var( ncid, vid, hematite_frac_in_aitken )     
    ierr = pio_inq_varid( ncid, 'FracH2', vid )            !FracH2
    ierr = pio_get_var( ncid, vid, hematite_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracH3', vid )            !FracH3
    ierr = pio_get_var( ncid, vid, hematite_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracH4', vid )            !FracH4
    ierr = pio_get_var( ncid, vid, hematite_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracH5', vid )            !FracH5
    ierr = pio_get_var( ncid, vid, hematite_frac_in_giant )
        
    ierr = pio_inq_varid( ncid, 'FracQ1', vid )            !FracQ1
    ierr = pio_get_var( ncid, vid, quartz_frac_in_aitken )       
    ierr = pio_inq_varid( ncid, 'FracQ2', vid )            !FracQ2
    ierr = pio_get_var( ncid, vid, quartz_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracQ3', vid )            !FracQ3
    ierr = pio_get_var( ncid, vid, quartz_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracQ4', vid )            !FracQ4
    ierr = pio_get_var( ncid, vid, quartz_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracQ5', vid )            !FracQ5
    ierr = pio_get_var( ncid, vid, quartz_frac_in_giant )
        
    ierr = pio_inq_varid( ncid, 'FracC1', vid )            !FracC1
    ierr = pio_get_var( ncid, vid, calcite_frac_in_aitken )      
    ierr = pio_inq_varid( ncid, 'FracC2', vid )            !FracC2
    ierr = pio_get_var( ncid, vid, calcite_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracC3', vid )            !FracC3
    ierr = pio_get_var( ncid, vid, calcite_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracC4', vid )            !FracC4
    ierr = pio_get_var( ncid, vid, calcite_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracC5', vid )            !FracC5
    ierr = pio_get_var( ncid, vid, calcite_frac_in_giant )
        
    ierr = pio_inq_varid( ncid, 'FracF1', vid )            !FracF1
    ierr = pio_get_var( ncid, vid, feldspar_frac_in_aitken )     
    ierr = pio_inq_varid( ncid, 'FracF2', vid )            !FracF2
    ierr = pio_get_var( ncid, vid, feldspar_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracF3', vid )            !FracF3
    ierr = pio_get_var( ncid, vid, feldspar_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracF4', vid )            !FracF4
    ierr = pio_get_var( ncid, vid, feldspar_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracF5', vid )            !FracF5
    ierr = pio_get_var( ncid, vid, feldspar_frac_in_giant )
        
    ierr = pio_inq_varid( ncid, 'FracG1', vid )            !FracG1  
    ierr = pio_get_var( ncid, vid, gypsum_frac_in_aitken )       
    ierr = pio_inq_varid( ncid, 'FracG2', vid )            !FracG2
    ierr = pio_get_var( ncid, vid, gypsum_frac_in_accum )
    ierr = pio_inq_varid( ncid, 'FracG3', vid )            !FracG3
    ierr = pio_get_var( ncid, vid, gypsum_frac_in_coarse01 )
    ierr = pio_inq_varid( ncid, 'FracG4', vid )            !FracG4
    ierr = pio_get_var( ncid, vid, gypsum_frac_in_coarse02 )
    ierr = pio_inq_varid( ncid, 'FracG5', vid )            !FracG5
    ierr = pio_get_var( ncid, vid, gypsum_frac_in_giant )

#endif
!++Longlei Li--
    ierr = pio_inq_varid( ncid, 'corrFact', vid )            
    ierr = pio_get_var( ncid, vid, correct_factor )
!--Longlei Li--

    !-----------------------------------------------------------------------
    !     	... convert to radians and setup regridding
    !-----------------------------------------------------------------------
    dst_lats(:) = d2r * dst_lats(:)
    dst_lons(:) = d2r * dst_lons(:)

    allocate( soil_erodibility(pcols,begchunk:endchunk), stat=ierr )
#if ( defined DUSTTYPES )
  ! Longlei Li - BRIFT already applied offline
    allocate( illite_frac_accum(pcols,begchunk:endchunk) )
        allocate( illite_frac_aitken(pcols,begchunk:endchunk) )
    allocate( illite_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( illite_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( illite_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( kaolinite_frac_accum(pcols,begchunk:endchunk) )
        allocate( kaolinite_frac_aitken(pcols,begchunk:endchunk) )
    allocate( kaolinite_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( kaolinite_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( kaolinite_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( mont_frac_accum(pcols,begchunk:endchunk) )
        allocate( mont_frac_aitken(pcols,begchunk:endchunk) )
    allocate( mont_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( mont_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( mont_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( hematite_frac_accum(pcols,begchunk:endchunk) )
        allocate( hematite_frac_aitken(pcols,begchunk:endchunk) )
    allocate( hematite_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( hematite_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( hematite_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( quartz_frac_accum(pcols,begchunk:endchunk) )
        allocate( quartz_frac_aitken(pcols,begchunk:endchunk) )
    allocate( quartz_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( quartz_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( quartz_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( calcite_frac_accum(pcols,begchunk:endchunk) )
        allocate( calcite_frac_aitken(pcols,begchunk:endchunk) )
    allocate( calcite_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( calcite_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( calcite_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( feldspar_frac_accum(pcols,begchunk:endchunk) )
        allocate( feldspar_frac_aitken(pcols,begchunk:endchunk) )
    allocate( feldspar_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( feldspar_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( feldspar_frac_giant(pcols,begchunk:endchunk) )
        
    allocate( gypsum_frac_accum(pcols,begchunk:endchunk) )
        allocate( gypsum_frac_aitken(pcols,begchunk:endchunk) )
    allocate( gypsum_frac_coarse01(pcols,begchunk:endchunk) )
    allocate( gypsum_frac_coarse02(pcols,begchunk:endchunk) )
    allocate( gypsum_frac_giant(pcols,begchunk:endchunk) )
!++Longlei Li--
#endif
!++Longlei Li--
   allocate( asphericity_factor(pcols) )
!--Longlei Li--

    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to allocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to allocate soil_erodibility_in')
    end if

    !-----------------------------------------------------------------------
    !     	... regrid ..
    !-----------------------------------------------------------------------
    do c=begchunk,endchunk
       ncols = get_ncols_p(c)
       call get_rlat_all_p(c, pcols, to_lats)
       call get_rlon_all_p(c, pcols, to_lons)

       call lininterp_init(dst_lons, nlon, to_lons, ncols, 2, lon_wgts, zero, twopi)
       call lininterp_init(dst_lats, nlat, to_lats, ncols, 1, lat_wgts)

       call lininterp(soil_erodibility_in(:,:), nlon,nlat , soil_erodibility(:,c), ncols, lon_wgts,lat_wgts)

#if ( defined DUSTTYPES )
! Longlei Li - BRIFT already applied offline
       call lininterp(illite_frac_in_accum(:,:), nlon, nlat, illite_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(illite_frac_in_aitken(:,:), nlon, nlat, illite_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(illite_frac_in_coarse01(:,:), nlon, nlat, illite_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(illite_frac_in_coarse02(:,:), nlon, nlat, illite_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(illite_frac_in_giant(:,:), nlon, nlat, illite_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(kaolinite_frac_in_accum(:,:), nlon, nlat, kaolinite_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(kaolinite_frac_in_aitken(:,:), nlon, nlat, kaolinite_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(kaolinite_frac_in_coarse01(:,:), nlon, nlat, kaolinite_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(kaolinite_frac_in_coarse02(:,:), nlon, nlat, kaolinite_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(kaolinite_frac_in_giant(:,:), nlon, nlat, kaolinite_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(mont_frac_in_accum(:,:), nlon, nlat, mont_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(mont_frac_in_aitken(:,:), nlon, nlat, mont_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(mont_frac_in_coarse01(:,:), nlon, nlat, mont_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(mont_frac_in_coarse02(:,:), nlon, nlat, mont_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(mont_frac_in_giant(:,:), nlon, nlat, mont_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(hematite_frac_in_accum(:,:), nlon, nlat, hematite_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(hematite_frac_in_aitken(:,:), nlon, nlat, hematite_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(hematite_frac_in_coarse01(:,:), nlon, nlat, hematite_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(hematite_frac_in_coarse02(:,:), nlon, nlat, hematite_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(hematite_frac_in_giant(:,:), nlon, nlat, hematite_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(quartz_frac_in_accum(:,:), nlon, nlat, quartz_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(quartz_frac_in_aitken(:,:), nlon, nlat, quartz_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(quartz_frac_in_coarse01(:,:), nlon, nlat, quartz_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(quartz_frac_in_coarse02(:,:), nlon, nlat, quartz_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(quartz_frac_in_giant(:,:), nlon, nlat, quartz_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(calcite_frac_in_accum(:,:), nlon, nlat, calcite_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(calcite_frac_in_aitken(:,:), nlon, nlat, calcite_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(calcite_frac_in_coarse01(:,:), nlon, nlat, calcite_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(calcite_frac_in_coarse02(:,:), nlon, nlat, calcite_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(calcite_frac_in_giant(:,:), nlon, nlat, calcite_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(feldspar_frac_in_accum(:,:), nlon, nlat, feldspar_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(feldspar_frac_in_aitken(:,:), nlon, nlat, feldspar_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(feldspar_frac_in_coarse01(:,:), nlon, nlat, feldspar_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(feldspar_frac_in_coarse02(:,:), nlon, nlat, feldspar_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(feldspar_frac_in_giant(:,:), nlon, nlat, feldspar_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
        
       call lininterp(gypsum_frac_in_accum(:,:), nlon, nlat, gypsum_frac_accum(:,c), ncols, lon_wgts, lat_wgts)
           call lininterp(gypsum_frac_in_aitken(:,:), nlon, nlat, gypsum_frac_aitken(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(gypsum_frac_in_coarse01(:,:), nlon, nlat, gypsum_frac_coarse01(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(gypsum_frac_in_coarse02(:,:), nlon, nlat, gypsum_frac_coarse02(:,c), ncols, lon_wgts, lat_wgts)
       call lininterp(gypsum_frac_in_giant(:,:), nlon, nlat, gypsum_frac_giant(:,c), ncols, lon_wgts, lat_wgts)
!++Longlei Li--
#endif
!++Longlei Li--
       call lininterp(correct_factor(:,:),    nlon, nlat, asphericity_factor(:),    ncols, lon_wgts, lat_wgts)
!--Longlei Li--

       call lininterp_finish(lat_wgts)
       call lininterp_finish(lon_wgts)
    end do
    deallocate( soil_erodibility_in, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'soil_erod_init: failed to deallocate soil_erodibility_in, ierr = ',ierr
       call endrun('soil_erod_init: failed to deallocate soil_erodibility_in')
    end if
#if ( defined DUSTTYPES )
  ! Longlei Li - BRIFT already applied offline
    deallocate( illite_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate illite_frac_in_accum, ierr = ',ierr
       call endrun
    end if
    deallocate( illite_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate illite_frac_in_aitken, ierr = ',ierr
       call endrun
    end if      
    deallocate( illite_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate illite_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( illite_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate illite_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( illite_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate illite_frac_in_giant, ierr = ',ierr
       call endrun
    end if

    deallocate( kaolinite_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( kaolinite_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( kaolinite_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( kaolinite_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( kaolinite_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate kaolinite_frac_in_giant, ierr = ',ierr
       call endrun
    end if

    deallocate( mont_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate mont_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( mont_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate mont_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( mont_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate mont_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( mont_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate mont_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( mont_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate mont_frac_in_giant, ierr = ',ierr
       call endrun
    end if

    deallocate( hematite_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( hematite_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( hematite_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( hematite_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( hematite_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate hematite_frac_in_giant, ierr = ',ierr
       call endrun
    end if

    deallocate( quartz_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( quartz_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( quartz_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( quartz_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( quartz_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate quartz_frac_in_giant, ierr = ',ierr
       call endrun
    end if
        
    deallocate( calcite_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( calcite_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( calcite_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( calcite_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( calcite_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate calcite_frac_in_giant, ierr = ',ierr
       call endrun
    end if

    deallocate( feldspar_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( feldspar_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( feldspar_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( feldspar_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( feldspar_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate feldspar_frac_in_giant, ierr = ',ierr
       call endrun
    end if
        
    deallocate( gypsum_frac_in_accum, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate gypsum_frac_in_accum, ierr = ',ierr
       call endrun
    end if
        deallocate( gypsum_frac_in_aitken, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate gypsum_frac_in_aitken, ierr = ',ierr
       call endrun
    end if
    deallocate( gypsum_frac_in_coarse01, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate gypsum_frac_in_coarse01, ierr = ',ierr
       call endrun
    end if
    deallocate( gypsum_frac_in_coarse02, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate gypsum_frac_in_coarse02, ierr = ',ierr
       call endrun
    end if
    deallocate( gypsum_frac_in_giant, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate gypsum_frac_in_giant, ierr = ',ierr
       call endrun
    end if
!++Longlei Li--
#endif
!++Longlei Li--
    deallocate( correct_factor, stat=ierr )
    if( ierr /= 0 ) then
       write(iulog,*) 'dust_initialize: failed to deallocate correct_factor, ierr = ',ierr
       call endrun
    end if
!++Longlei Li--

    deallocate( dst_lats )
    deallocate( dst_lons )

  end  subroutine soil_erod_init

end module soil_erod_mod
