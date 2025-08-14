module ndrop

!---------------------------------------------------------------------------------
! Purpose:
!   CAM Interface for droplet activation by modal aerosols
!
! ***N.B.*** This module is currently hardcoded to recognize only the modes that
!            affect the climate calculation.  This is implemented by using list
!            index 0 in all the calls to rad_constituent interfaces.
!---------------------------------------------------------------------------------

use shr_kind_mod,     only: r8 => shr_kind_r8
use spmd_utils,       only: masterproc
use ppgrid,           only: pcols, pver, pverp
use physconst,        only: pi, rhoh2o, mwh2o, r_universal, rh2o, &
                            gravit, latvap, cpair, epsilo, rair
use constituents,     only: pcnst, cnst_get_ind, cnst_name, cnst_spec_class_gas, cnst_species_class
use physics_types,    only: physics_state, physics_ptend, physics_ptend_init
use physics_buffer,   only: physics_buffer_desc, pbuf_get_index, pbuf_get_field

use wv_saturation,    only: qsat
use phys_control,     only: phys_getopts
use ref_pres,         only: top_lev => trop_cloud_top_lev
use shr_spfn_mod,     only: erf => shr_spfn_erf
use rad_constituents, only: rad_cnst_get_info, rad_cnst_get_mode_num, rad_cnst_get_aer_mmr, &
                            rad_cnst_get_aer_props, rad_cnst_get_mode_props,                &
                            rad_cnst_get_mam_mmr_idx, rad_cnst_get_mode_num_idx
use cam_history,      only: addfld, add_default, horiz_only, fieldname_len, outfld
use cam_abortutils,   only: endrun
use cam_logfile,      only: iulog

implicit none
private
save

public ndrop_init, dropmixnuc, activate_modal, loadaer

real(r8), allocatable :: alogsig(:)     ! natl log of geometric standard dev of aerosol
real(r8), allocatable :: exp45logsig(:)
real(r8), allocatable :: f1(:)          ! abdul-razzak functions of width
real(r8), allocatable :: f2(:)          ! abdul-razzak functions of width

real(r8) :: t0            ! reference temperature
real(r8) :: aten
real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
real(r8) :: alog2, alog3, alogaten
real(r8) :: third, twothird, sixth, zero
real(r8) :: sq2, sqpi

! CCN diagnostic fields
!#if (defined IRONPROC)
integer,  parameter :: psat=13    ! number of supersaturations to calc ccn concentration
real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
     (/ 0.05_r8, 0.075_r8, 0.1_r8, 0.15_r8, 0.2_r8, 0.3_r8, &
     0.4_r8, 0.5_r8, 0.6_r8, 0.7_r8, 0.8_r8, 0.9_r8, 1.0_r8 /)
character(len=8) :: ccn_name(psat)= &
     (/'CCN005','CCN0075','CCN01','CCN015','CCN02','CCN03','CCN04', &
     'CCN05','CCN06','CCN07','CCN08','CCN09','CCN10'/)

integer, parameter :: pdiam=3    ! number of diameters to calc ccn concs/kappa 
real(r8), parameter :: diam(pdiam) = (/ 50.E-9_r8, 80.E-9_r8, 120.E-9_r8 /)
character(len=8) :: ccn_name_diam(pdiam) = (/ 'N50','N80','N120'/)
character(len=8) :: kap_name(pdiam)= (/ 'HG50','HG80','HG120' /)

! PM diagnostics
integer,  parameter :: pPM=7  ! number of PM species
character(len=8) :: PM_name(pPM)= &
     (/ 'SO4','POM','SOA','ORG','BC','SS','DU' /)
integer,  parameter :: pPM1=7 ! number of PM1 species
character(len=8) :: PM1_name(pPM1)= &
     (/ 'SO4PM1','pORGPM1','sORGPM1','ORGPM1','BCPM1','SSPM1','DUPM1' /)
integer,  parameter :: pPM25=7 ! number of PM2.5aerodynamic = PM1.6ge0metric species
character(len=8) :: PM25_name(pPM25)= &
     (/ 'SO4PM25','pORGPM25','sORGPM25','ORGPM25','BCPM25','SSPM25','DUPM25' /)
integer,  parameter :: pPM10=7 ! number of PM10 species = PM6.9ge0metric
character(len=8) :: PM10_name(pPM10)= &
     (/ 'SO4PM10','pORGPM10','sORGPM10','ORGPM10','BCPM10','SSPM10','DUPM10' /)
integer,  parameter :: pPMt=7 ! number of PMtotal species
character(len=9) :: tot_name(pPMt)= &
     (/ 'SO4PMtot','pORGPMtot','sORGPMtot','ORGPMtot','BCPMtot','SSPMtot','DUPMtot' /)
!#else
!integer,  parameter :: psat=6    ! number of supersaturations to calc ccn concentration
!real(r8), parameter :: supersat(psat)= & ! supersaturation (%) to determine ccn concentration
!                       (/ 0.02_r8, 0.05_r8, 0.1_r8, 0.2_r8, 0.5_r8, 1.0_r8 /)
!character(len=8) :: ccn_name(psat)= &
!                    (/'CCN1','CCN2','CCN3','CCN4','CCN5','CCN6'/)
!#endif

! indices in state and pbuf structures
integer :: numliq_idx = -1
integer :: kvh_idx    = -1

! description of modal aerosols
integer               :: ntot_amode     ! number of aerosol modes
integer,  allocatable :: nspec_amode(:) ! number of chemical species in each aerosol mode
real(r8), allocatable :: sigmag_amode(:)! geometric standard deviation for each aerosol mode
real(r8), allocatable :: dgnumlo_amode(:)
real(r8), allocatable :: dgnumhi_amode(:)
real(r8), allocatable :: voltonumblo_amode(:)
real(r8), allocatable :: voltonumbhi_amode(:)

logical :: history_aerosol      ! Output the MAM aerosol tendencies
character(len=fieldname_len), allocatable :: fieldname(:)    ! names for drop nuc tendency output fields
character(len=fieldname_len), allocatable :: fieldname_cw(:) ! names for drop nuc tendency output fields

! local indexing for MAM
integer, allocatable :: mam_idx(:,:) ! table for local indexing of modal aero number and mmr
integer :: ncnst_tot                  ! total number of mode number conc + mode species

! Indices for MAM species in the ptend%q array.  Needed for prognostic aerosol case.
integer, allocatable :: mam_cnst_idx(:,:)


! ptr2d_t is used to create arrays of pointers to 2D fields
type ptr2d_t
   real(r8), pointer :: fld(:,:)
end type ptr2d_t

! modal aerosols
logical :: prog_modal_aero     ! true when modal aerosols are prognostic
logical :: lq(pcnst) = .false. ! set flags true for constituents with non-zero tendencies
                               ! in the ptend object

!===============================================================================
contains
!===============================================================================

subroutine ndrop_init

   integer  :: ii, l, lptr, m, mm
   integer  :: nspec_max            ! max number of species in a mode
   character(len=32)   :: tmpname
   character(len=32)   :: tmpname_cw
   character(len=128)  :: long_name
   character(len=8)    :: unit
   logical :: history_amwg         ! output the variables used by the AMWG diag package

   !-------------------------------------------------------------------------------

   ! get indices into state%q and pbuf structures
   call cnst_get_ind('NUMLIQ', numliq_idx)

   kvh_idx      = pbuf_get_index('kvh')

   zero     = 0._r8
   third    = 1._r8/3._r8
   twothird = 2._r8*third
   sixth    = 1._r8/6._r8
   sq2      = sqrt(2._r8)
   sqpi     = sqrt(pi)

   t0       = 273._r8
   surften  = 0.076_r8
   aten     = 2._r8*mwh2o*surften/(r_universal*t0*rhoh2o)
   alogaten = log(aten)
   alog2    = log(2._r8)
   alog3    = log(3._r8)

   ! get info about the modal aerosols
   ! get ntot_amode
   call rad_cnst_get_info(0, nmodes=ntot_amode)

   allocate( &
      nspec_amode(ntot_amode),  &
      sigmag_amode(ntot_amode), &
      dgnumlo_amode(ntot_amode), &
      dgnumhi_amode(ntot_amode), &
      alogsig(ntot_amode),      &
      exp45logsig(ntot_amode),  &
      f1(ntot_amode),           &
      f2(ntot_amode),           &
      voltonumblo_amode(ntot_amode), &
      voltonumbhi_amode(ntot_amode)  )

   do m = 1, ntot_amode
      ! use only if width of size distribution is prescribed

      ! get mode info
      call rad_cnst_get_info(0, m, nspec=nspec_amode(m))

      ! get mode properties
      call rad_cnst_get_mode_props(0, m, sigmag=sigmag_amode(m),  &
         dgnumhi=dgnumhi_amode(m), dgnumlo=dgnumlo_amode(m))

      alogsig(m)     = log(sigmag_amode(m))
      exp45logsig(m) = exp(4.5_r8*alogsig(m)*alogsig(m))
      f1(m)          = 0.5_r8*exp(2.5_r8*alogsig(m)*alogsig(m))
      f2(m)          = 1._r8 + 0.25_r8*alogsig(m)

      voltonumblo_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
                             (dgnumlo_amode(m)**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
      voltonumbhi_amode(m) = 1._r8 / ( (pi/6._r8)*                          &
                             (dgnumhi_amode(m)**3._r8)*exp(4.5_r8*alogsig(m)**2._r8) )
   end do
      
   ! Init the table for local indexing of mam number conc and mmr.
   ! This table uses species index 0 for the number conc.

   ! Find max number of species in all the modes, and the total
   ! number of mode number concentrations + mode species
   nspec_max = nspec_amode(1)
   ncnst_tot = nspec_amode(1) + 1
   do m = 2, ntot_amode
      nspec_max = max(nspec_max, nspec_amode(m))
      ncnst_tot = ncnst_tot + nspec_amode(m) + 1
   end do

   allocate( &
      mam_idx(ntot_amode,0:nspec_max),      &
      mam_cnst_idx(ntot_amode,0:nspec_max), &
      fieldname(ncnst_tot),                 &
      fieldname_cw(ncnst_tot)               )

   ! Local indexing compresses the mode and number/mass indicies into one index.
   ! This indexing is used by the pointer arrays used to reference state and pbuf
   ! fields.
   ii = 0
   do m = 1, ntot_amode
      do l = 0, nspec_amode(m)
         ii = ii + 1
         mam_idx(m,l) = ii
      end do
   end do

   ! Add dropmixnuc tendencies for all modal aerosol species

   call phys_getopts(history_amwg_out = history_amwg, &
                     history_aerosol_out = history_aerosol, &
                     prog_modal_aero_out=prog_modal_aero)


   do m = 1, ntot_amode
      do l = 0, nspec_amode(m)   ! loop over number + chem constituents

         mm = mam_idx(m,l)

         unit = 'kg/m2/s'
         if (l == 0) then   ! number
            unit = '#/m2/s'
         end if

         if (l == 0) then   ! number
            call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
         else
            call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
         end if

         fieldname(mm)    = trim(tmpname) // '_mixnuc1'
         fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'

         if (prog_modal_aero) then

            ! To set tendencies in the ptend object need to get the constituent indices
            ! for the prognostic species
            if (l == 0) then   ! number
               call rad_cnst_get_mode_num_idx(m, lptr)
            else
               call rad_cnst_get_mam_mmr_idx(m, l, lptr)
            end if
            mam_cnst_idx(m,l) = lptr
            lq(lptr)          = .true.

            ! Add tendency fields to the history only when prognostic MAM is enabled.
            long_name = trim(tmpname) // ' dropmixnuc mixnuc column tendency'
            call addfld(fieldname(mm),    horiz_only, 'A', unit, long_name)

            long_name = trim(tmpname_cw) // ' dropmixnuc mixnuc column tendency'
            call addfld(fieldname_cw(mm), horiz_only, 'A', unit, long_name)

            if (history_aerosol) then
               call add_default(fieldname(mm), 1, ' ')
               call add_default(fieldname_cw(mm), 1, ' ')
            end if



         end if
            
      end do
   end do
!#if (defined IRONPROC)
   call addfld('CCN005  ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.05%')
   call addfld('CCN0075 ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.075%')
   call addfld('CCN01   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.1%')
   call addfld('CCN015  ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.15%')
   call addfld('CCN02   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.2%')
   call addfld('CCN03   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.3%')
   call addfld('CCN04   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.4%')
   call addfld('CCN05   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.5%')
   call addfld('CCN06   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.6%')
   call addfld('CCN07   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.7%')
   call addfld('CCN08   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.8%')
   call addfld('CCN09   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=0.9%')
   call addfld('CCN10   ',(/ 'lev' /), 'A','#/cm3   ','CCN concentration at S=1.0%')

   call addfld('N50     ',(/ 'lev' /), 'A','#/cm3   ','Number conc. > than 50nm' )
   call addfld('N80     ',(/ 'lev' /), 'A','#/cm3   ','Number conc. > than 80nm' )
   call addfld('N120    ',(/ 'lev' /), 'A','#/cm3   ','Number conc. > than 120nm')

   call addfld('HG50    ',(/ 'lev' /), 'A','#       ','Aerosol hygroscopicity at 50nm')
   call addfld('HG80    ',(/ 'lev' /), 'A','#       ','Aerosol hygroscopicity at 80nm')
   call addfld('HG120   ',(/ 'lev' /), 'A','#       ','Aerosol hygroscopicity at 120nm')

   call addfld('SO4PM1  ', horiz_only, 'A','ug/m3   ', 'SO4 mass for PM1 aerosols'      )
   call addfld('pORGPM1 ', horiz_only, 'A','ug/m3   ', 'POA mass for PM1 aerosols'      )
   call addfld('sORGPM1 ', horiz_only, 'A','ug/m3   ', 'SOA mass for PM1 aerosols'      )
   call addfld('ORGPM1  ', horiz_only, 'A','ug/m3   ', 'Organic mass for PM1 aerosols'  )
   call addfld('BCPM1   ', horiz_only, 'A','ug/m3   ', 'BC mass for PM1 aerosols'       )
   call addfld('SSPM1   ', horiz_only, 'A','ug/m3   ', 'Sea salt mass for PM1 aerosols' )
   call addfld('DUPM1   ', horiz_only, 'A','ug/m3   ', 'Dust mass for PM1 aerosols'     )

   call addfld('ALLPM25 ', horiz_only, 'A','ug/m3   ', 'All aerosols PM2.5 aerosols'      )
   call addfld('SO4PM25 ', horiz_only, 'A','ug/m3   ', 'SO4 mass for PM2.5 aerosols'      )
   call addfld('pORGPM25', horiz_only, 'A','ug/m3   ', 'POA mass for PM2.5 aerosols'      )
   call addfld('sORGPM25', horiz_only, 'A','ug/m3   ', 'SOA mass for PM2.5 aerosols'      )
   call addfld('ORGPM25 ', horiz_only, 'A','ug/m3   ', 'Organic mass for PM2.5 aerosols'  )
   call addfld('BCPM25  ', horiz_only, 'A','ug/m3   ', 'BC mass for PM2.5 aerosols'       )
   call addfld('SSPM25  ', horiz_only, 'A','ug/m3   ', 'Sea salt mass for PM2.5 aerosols' )
   call addfld('DUPM25  ', horiz_only, 'A','ug/m3   ', 'Dust mass for PM2.5 aerosols'     )

   call addfld('ALLPM10 ', horiz_only, 'A','ug/m3   ', 'All aerosols PM10 aerosols'      )
   call addfld('SO4PM10 ', horiz_only, 'A','ug/m3   ', 'SO4 mass for PM10 aerosols'      )
   call addfld('pORGPM10', horiz_only, 'A','ug/m3   ', 'POA mass for PM10 aerosols'      )
   call addfld('sORGPM10', horiz_only, 'A','ug/m3   ', 'SOA mass for PM10 aerosols'      )
   call addfld('ORGPM10 ', horiz_only, 'A','ug/m3   ', 'Organic mass for PM10 aerosols'  )
   call addfld('BCPM10  ', horiz_only, 'A','ug/m3   ', 'BC mass for PM10 aerosols'       )
   call addfld('SSPM10  ', horiz_only, 'A','ug/m3   ', 'Sea salt mass for PM10 aerosols' )
   call addfld('DUPM10  ', horiz_only, 'A','ug/m3   ', 'Dust mass for PM10 aerosols'     )

   call addfld('SO4PMtot ', horiz_only, 'A','ug/m3   ', 'SO4 total mass'      )
   call addfld('pORGPMtot', horiz_only, 'A','ug/m3   ', 'POA total mass'      )
   call addfld('sORGPMtot', horiz_only, 'A','ug/m3   ', 'SOA total mass'      )
   call addfld('ORGPMtot ', horiz_only, 'A','ug/m3   ', 'Organic total mass'  )
   call addfld('BCPMtot  ', horiz_only, 'A','ug/m3   ', 'BC total mass'       )
   call addfld('SSPMtot  ', horiz_only, 'A','ug/m3   ', 'Sea salt total mass' )
   call addfld('DUPMtot  ', horiz_only, 'A','ug/m3   ', 'Dust total mass'     )

   call addfld('cs       ',(/ 'lev' /), 'A','kg/m3   ', 'Air Density'                     )
   call addfld('dz       ',(/ 'lev' /), 'A','m       ', 'geometric thickness of layers (m)')
   call addfld('gph      ',(/ 'lev' /), 'A','m       ', 'Geopotential Height (above surf)')

   ! set the add_default fields  
   if (history_amwg) then
      call add_default('CCN03', 1, ' ')
   endif
!#else
!   call addfld('CCN1',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.02%')
!   call addfld('CCN2',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.05%')
!   call addfld('CCN3',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.1%')
!   call addfld('CCN4',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.2%')
!   call addfld('CCN5',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=0.5%')
!   call addfld('CCN6',(/ 'lev' /), 'A','#/cm3','CCN concentration at S=1.0%')
!
!   ! set the add_default fields  
!   if (history_amwg) then
!      call add_default('CCN3', 1, ' ')
!   endif
!#endif

   call addfld('WTKE',     (/ 'lev' /), 'A', 'm/s', 'Standard deviation of updraft velocity')
   call addfld('NDROPMIX', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number mixing')
   call addfld('NDROPSRC', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number source')
   call addfld('NDROPSNK', (/ 'lev' /), 'A', '#/kg/s', 'Droplet number loss by microphysics')
   call addfld('NDROPCOL', horiz_only,  'A', '#/m2', 'Column droplet number')

   ! set the add_default fields  
   !if (history_amwg) then
   !   call add_default('CCN3', 1, ' ')
   !endif  !-lli

   if (history_aerosol .and. prog_modal_aero) then
     do m = 1, ntot_amode
        do l = 0, nspec_amode(m)   ! loop over number + chem constituents
           mm = mam_idx(m,l)
           if (l == 0) then   ! number
              call rad_cnst_get_info(0, m, num_name=tmpname, num_name_cw=tmpname_cw)
           else
              call rad_cnst_get_info(0, m, l, spec_name=tmpname, spec_name_cw=tmpname_cw)
           end if
           fieldname(mm)    = trim(tmpname) // '_mixnuc1'
           fieldname_cw(mm) = trim(tmpname_cw) // '_mixnuc1'
        end do
     end do
   endif


end subroutine ndrop_init

!===============================================================================

subroutine dropmixnuc( &
   state, ptend, dtmicro, pbuf, wsub, &
   cldn, cldo, cldliqf, tendnd, factnum, from_spcam)

   ! vertical diffusion and nucleation of cloud droplets
   ! assume cloud presence controlled by cloud fraction
   ! doesn't distinguish between warm, cold clouds

   ! arguments
   type(physics_state), target, intent(in)    :: state
   type(physics_ptend),         intent(out)   :: ptend
   real(r8),                    intent(in)    :: dtmicro     ! time step for microphysics (s)

   type(physics_buffer_desc), pointer :: pbuf(:)

   ! arguments
   real(r8), intent(in) :: wsub(pcols,pver)    ! subgrid vertical velocity
   real(r8), intent(in) :: cldn(pcols,pver)    ! cloud fraction
   real(r8), intent(in) :: cldo(pcols,pver)    ! cloud fraction on previous time step
   real(r8), intent(in) :: cldliqf(pcols,pver) ! liquid cloud fraction (liquid / (liquid + ice))
   logical,  intent(in),optional :: from_spcam ! value insignificant - if variable present, is called from spcam

   ! output arguments
   real(r8), intent(out) :: tendnd(pcols,pver) ! change in droplet number concentration (#/kg/s)
   real(r8), intent(out) :: factnum(:,:,:)     ! activation fraction for aerosol number
   !--------------------Local storage-------------------------------------

   integer  :: lchnk               ! chunk identifier
   integer  :: ncol                ! number of columns

   real(r8), pointer :: ncldwtr(:,:) ! droplet number concentration (#/kg)
   real(r8), pointer :: temp(:,:)    ! temperature (K)
   real(r8), pointer :: omega(:,:)   ! vertical velocity (Pa/s)
   real(r8), pointer :: pmid(:,:)    ! mid-level pressure (Pa)
   real(r8), pointer :: pint(:,:)    ! pressure at layer interfaces (Pa)
   real(r8), pointer :: pdel(:,:)    ! pressure thickess of layer (Pa)
   real(r8), pointer :: rpdel(:,:)   ! inverse of pressure thickess of layer (/Pa)
   real(r8), pointer :: zm(:,:)      ! geopotential height of level (m)

   real(r8), pointer :: kvh(:,:)     ! vertical diffusivity (m2/s)

   type(ptr2d_t), allocatable :: raer(:)     ! aerosol mass, number mixing ratios
   type(ptr2d_t), allocatable :: qqcw(:)
   real(r8) :: raertend(pver)  ! tendency of aerosol mass, number mixing ratios
   real(r8) :: qqcwtend(pver)  ! tendency of cloudborne aerosol mass, number mixing ratios


   real(r8), parameter :: zkmin = 0.01_r8, zkmax = 100._r8
   real(r8), parameter :: wmixmin = 0.1_r8        ! minimum turbulence vertical velocity (m/s)
   real(r8) :: sq2pi

   integer  :: i, k, l, m, mm, n
   integer  :: km1, kp1
   integer  :: nnew, nsav, ntemp
   integer  :: lptr
   integer  :: nsubmix, nsubmix_bnd
   integer, save :: count_submix(100)
   integer  :: phase ! phase of aerosol

   real(r8) :: arg
   real(r8) :: dtinv
   real(r8) :: dtmin, tinv, dtt
   real(r8) :: lcldn(pcols,pver)
   real(r8) :: lcldo(pcols,pver)

   real(r8) :: zs(pver) ! inverse of distance between levels (m)
   real(r8) :: qcld(pver) ! cloud droplet number mixing ratio (#/kg)
   real(r8) :: qncld(pver)     ! droplet number nucleated on cloud boundaries
   real(r8) :: srcn(pver)       ! droplet source rate (/s)
   real(r8) :: cs(pcols,pver)      ! air density (kg/m3)
   real(r8) :: csbot(pver)       ! air density at bottom (interface) of layer (kg/m3)
   real(r8) :: csbot_cscen(pver) ! csbot(i)/cs(i,k)
   real(r8) :: dz(pcols,pver)      ! geometric thickness of layers (m)

   real(r8) :: wtke(pcols,pver)     ! turbulent vertical velocity at base of layer k (m/s)
   real(r8) :: wtke_cen(pcols,pver) ! turbulent vertical velocity at center of layer k (m/s)
   real(r8) :: wbar, wmix, wmin, wmax

   real(r8) :: zn(pver)   ! g/pdel (m2/g) for layer
   real(r8) :: flxconv    ! convergence of flux into lowest layer

   real(r8) :: wdiab           ! diabatic vertical velocity
   real(r8) :: ekd(pver)       ! diffusivity for droplets (m2/s)
   real(r8) :: ekk(0:pver)     ! density*diffusivity for droplets (kg/m3 m2/s)
   real(r8) :: ekkp(pver)      ! zn*zs*density*diffusivity
   real(r8) :: ekkm(pver)      ! zn*zs*density*diffusivity

   real(r8) :: dum, dumc
   real(r8) :: tmpa
   real(r8) :: dact
   real(r8) :: fluxntot         ! (#/cm2/s)
   real(r8) :: dtmix
   real(r8) :: alogarg
   real(r8) :: overlapp(pver), overlapm(pver) ! cloud overlap

   real(r8) :: nsource(pcols,pver)            ! droplet number source (#/kg/s)
   real(r8) :: ndropmix(pcols,pver)           ! droplet number mixing (#/kg/s)
   real(r8) :: ndropcol(pcols)               ! column droplet number (#/m2)
   real(r8) :: cldo_tmp, cldn_tmp
   real(r8) :: tau_cld_regenerate
   real(r8) :: zeroaer(pver)
   real(r8) :: taumix_internal_pver_inv ! 1/(internal mixing time scale for k=pver) (1/s)


   real(r8), allocatable :: nact(:,:)  ! fractional aero. number  activation rate (/s)
   real(r8), allocatable :: mact(:,:)  ! fractional aero. mass    activation rate (/s)

   real(r8), allocatable :: raercol(:,:,:)    ! single column of aerosol mass, number mixing ratios
   real(r8), allocatable :: raercol_cw(:,:,:) ! same as raercol but for cloud-borne phase


   real(r8) :: na(pcols), va(pcols), hy(pcols)
   real(r8), allocatable :: naermod(:)  ! (1/m3)
   real(r8), allocatable :: hygro(:)    ! hygroscopicity of aerosol mode
   real(r8), allocatable :: vaerosol(:) ! interstit+activated aerosol volume conc (cm3/cm3)

   real(r8) :: source(pver)

   real(r8), allocatable :: fn(:)              ! activation fraction for aerosol number
   real(r8), allocatable :: fm(:)              ! activation fraction for aerosol mass

   real(r8), allocatable :: fluxn(:)           ! number  activation fraction flux (cm/s)
   real(r8), allocatable :: fluxm(:)           ! mass    activation fraction flux (cm/s)
   real(r8)              :: flux_fullact(pver) ! 100%    activation fraction flux (cm/s)
   !     note:  activation fraction fluxes are defined as 
   !     fluxn = [flux of activated aero. number into cloud (#/cm2/s)]
   !           / [aero. number conc. in updraft, just below cloudbase (#/cm3)]


   real(r8), allocatable :: coltend(:,:)       ! column tendency for diagnostic output
   real(r8), allocatable :: coltend_cw(:,:)    ! column tendency
   real(r8) :: ccn(pcols,pver,psat)    ! number conc of aerosols activated at supersat

!#if (defined IRONPROC)
   real(r8) :: ccn_ND(pcols,pver,pdiam)  ! number conc of aerosols above given diamater
   real(r8) :: kappa(pcols,pver,pdiam)   ! kappa of aerosols above diam threshold (#/m3)
   real(r8) :: PMt(pcols,pver,pPM)       ! PM total of aerosols 
   real(r8) :: PM1(pcols,pver,pPM)       ! >PM1 of aerosols 
   real(r8) :: PM25(pcols,pver,pPM)      ! >PM25 of aerosols 
   real(r8) :: PMtout(pcols,pPM)         ! surface PM total of aerosols for output
   real(r8) :: PM1out(pcols,pPM)         ! surface PM1 of aerosols for output
   real(r8) :: PM25out(pcols,pPM)        ! surface PM25 of aerosols for output
   real(r8) :: PM25ALL(pcols)            ! surface PM25 of aerosols for output

   real(r8) :: PM10(pcols,pver,pPM)      ! >PM10 of aerosols 
   real(r8) :: PM10out(pcols,pPM)        ! surface PM10 of aerosols for output
   real(r8) :: PM10ALL(pcols)            ! surface PM10 of aerosols for output

   real(r8) :: massSO4(pcols,pver,ntot_amode)
   real(r8) :: masspORG(pcols,pver,ntot_amode)
   real(r8) :: masssORG(pcols,pver,ntot_amode)
   real(r8) :: massORG(pcols,pver,ntot_amode)
   real(r8) :: massBC(pcols,pver,ntot_amode)
   real(r8) :: massSS(pcols,pver,ntot_amode)
   real(r8) :: massDU(pcols,pver,ntot_amode)
   real(r8) :: massDU1(pcols,pver,ntot_amode)
   real(r8) :: massDU2(pcols,pver,ntot_amode)
   real(r8) :: massDU3(pcols,pver,ntot_amode)
   real(r8) :: massDU4(pcols,pver,ntot_amode)
   real(r8) :: massDU5(pcols,pver,ntot_amode)
   real(r8) :: massDU6(pcols,pver,ntot_amode)
   real(r8) :: massDU7(pcols,pver,ntot_amode)
   real(r8) :: massDU8(pcols,pver,ntot_amode)
!#endif 

   !for gas species turbulent mixing
   real(r8), pointer :: rgas(:, :, :)
   real(r8), allocatable :: rgascol(:, :, :)
   real(r8), allocatable :: coltendgas(:)
   real(r8) :: zerogas(pver)
   character*200 fieldnamegas

   logical  :: called_from_spcam
   !-------------------------------------------------------------------------------

   sq2pi = sqrt(2._r8*pi)

   lchnk = state%lchnk
   ncol  = state%ncol

   ncldwtr  => state%q(:,:,numliq_idx)
   temp     => state%t
   omega    => state%omega
   pmid     => state%pmid
   pint     => state%pint
   pdel     => state%pdel
   rpdel    => state%rpdel
   zm       => state%zm

   call pbuf_get_field(pbuf, kvh_idx, kvh)

   ! Create the liquid weighted cloud fractions that were passsed in
   ! before. This doesn't seem like the best variable, since the cloud could
   ! have liquid condensate, but the part of it that is changing could be the
   ! ice portion; however, this is what was done before.
   lcldo(:ncol,:)  = cldo(:ncol,:)  * cldliqf(:ncol,:)
   lcldn(:ncol,:) = cldn(:ncol,:) * cldliqf(:ncol,:)


   arg = 1.0_r8
   if (abs(0.8427_r8 - erf(arg))/0.8427_r8 > 0.001_r8) then
      write(iulog,*) 'erf(1.0) = ',ERF(arg)
      call endrun('dropmixnuc: Error function error')
   endif
   arg = 0.0_r8
   if (erf(arg) /= 0.0_r8) then
      write(iulog,*) 'erf(0.0) = ',erf(arg)
      write(iulog,*) 'dropmixnuc: Error function error'
      call endrun('dropmixnuc: Error function error')
   endif

   dtinv = 1._r8/dtmicro

   allocate( &
      nact(pver,ntot_amode),          &
      mact(pver,ntot_amode),          &
      raer(ncnst_tot),                &
      qqcw(ncnst_tot),                &
      raercol(pver,ncnst_tot,2),      &
      raercol_cw(pver,ncnst_tot,2),   &
      coltend(pcols,ncnst_tot),       &
      coltend_cw(pcols,ncnst_tot),    &
      naermod(ntot_amode),            &
      hygro(ntot_amode),              &
      vaerosol(ntot_amode),           &
      fn(ntot_amode),                 &
      fm(ntot_amode),                 &
      fluxn(ntot_amode),              &
      fluxm(ntot_amode)               )

   ! Init pointers to mode number and specie mass mixing ratios in 
   ! intersitial and cloud borne phases.
   do m = 1, ntot_amode
      mm = mam_idx(m, 0)
      call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, raer(mm)%fld)
      call rad_cnst_get_mode_num(0, m, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      do l = 1, nspec_amode(m)
         mm = mam_idx(m, l)
         call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, raer(mm)%fld)
         call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw(mm)%fld)  ! cloud-borne aerosol
      end do
   end do

   called_from_spcam = (present(from_spcam)) 

   if (called_from_spcam) then
      rgas  => state%q
      allocate(rgascol(pver, pcnst, 2))
      allocate(coltendgas(pcols))
   endif

   factnum = 0._r8
   wtke = 0._r8

   if (prog_modal_aero) then
      ! aerosol tendencies
      call physics_ptend_init(ptend, state%psetcols, 'ndrop', lq=lq)
   else
      ! no aerosol tendencies
      call physics_ptend_init(ptend, state%psetcols, 'ndrop')
   end if

   ! overall_main_i_loop
   do i = 1, ncol

      do k = top_lev, pver-1
         zs(k) = 1._r8/(zm(i,k) - zm(i,k+1))
      end do
      zs(pver) = zs(pver-1)

      ! load number nucleated into qcld on cloud boundaries

      do k = top_lev, pver

         qcld(k)  = ncldwtr(i,k)
         qncld(k) = 0._r8
         srcn(k)  = 0._r8
         cs(i,k)  = pmid(i,k)/(rair*temp(i,k))        ! air density (kg/m3)
         dz(i,k)  = 1._r8/(cs(i,k)*gravit*rpdel(i,k)) ! layer thickness in m

         do m = 1, ntot_amode
            nact(k,m) = 0._r8
            mact(k,m) = 0._r8
         end do

         zn(k) = gravit*rpdel(i,k)

         if (k < pver) then
            ekd(k)   = kvh(i,k+1)
            ekd(k)   = max(ekd(k), zkmin)
            ekd(k)   = min(ekd(k), zkmax)
            csbot(k) = 2.0_r8*pint(i,k+1)/(rair*(temp(i,k) + temp(i,k+1)))
            csbot_cscen(k) = csbot(k)/cs(i,k)
         else
            ekd(k)   = 0._r8
            csbot(k) = cs(i,k)
            csbot_cscen(k) = 1.0_r8
         end if

         ! rce-comment - define wtke at layer centers for new-cloud activation
         !    and at layer boundaries for old-cloud activation
         !++ag
         wtke_cen(i,k) = wsub(i,k)
         wtke(i,k)     = wsub(i,k)
         !--ag
         wtke_cen(i,k) = max(wtke_cen(i,k), wmixmin)
         wtke(i,k)     = max(wtke(i,k), wmixmin)

         nsource(i,k) = 0._r8

      end do

      nsav = 1
      nnew = 2
      do m = 1, ntot_amode
         mm = mam_idx(m,0)
         raercol_cw(:,mm,nsav) = 0.0_r8
         raercol(:,mm,nsav)    = 0.0_r8
         raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(i,top_lev:pver)
         raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(i,top_lev:pver)
         do l = 1, nspec_amode(m)
            mm = mam_idx(m,l)
            raercol_cw(top_lev:pver,mm,nsav) = qqcw(mm)%fld(i,top_lev:pver)
            raercol(top_lev:pver,mm,nsav)    = raer(mm)%fld(i,top_lev:pver)
         end do
      end do


      if (called_from_spcam) then
      !
      ! In the MMF model, turbulent mixing for tracer species are turned off.
      ! So the turbulent for gas species mixing are added here.
      ! (Previously, it had the turbulent mixing for aerosol species)
      !
         do m=1, pcnst
            if (cnst_species_class(m) == cnst_spec_class_gas) rgascol(:,m,nsav) = rgas(i,:,m)
         end do

      endif

      ! droplet nucleation/aerosol activation

      ! tau_cld_regenerate = time scale for regeneration of cloudy air 
      !    by (horizontal) exchange with clear air
      tau_cld_regenerate = 3600.0_r8 * 3.0_r8 

      if (called_from_spcam) then
      ! when this is called  in the MMF part, no cloud regeneration and decay.
      ! set the time scale be very long so that no cloud regeneration.
           tau_cld_regenerate = 3600.0_r8 * 24.0_r8 * 365.0_r8
      endif


      ! k-loop for growing/shrinking cloud calcs .............................
      ! grow_shrink_main_k_loop: &
      do k = top_lev, pver

         ! This code was designed for liquid clouds, but the cloudbourne
         ! aerosol can be either from liquid or ice clouds. For the ice clouds,
         ! we do not do regeneration, but as cloud fraction decreases the
         ! aerosols should be returned interstitial. The lack of a liquid cloud
         ! should not mean that all of the aerosol is realease. Therefor a
         ! section has been added for shrinking ice clouds and checks were added
         ! to protect ice cloudbourne aerosols from being released when no
         ! liquid cloud is present.

         ! shrinking ice cloud ......................................................
         cldo_tmp = cldo(i,k) * (1._r8 - cldliqf(i,k))
         cldn_tmp = cldn(i,k) * (1._r8 - cldliqf(i,k))

         if (cldn_tmp < cldo_tmp) then

            ! convert activated aerosol to interstitial in decaying cloud

            dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * (1._r8 - cldliqf(i,k))
            do m = 1, ntot_amode
               mm = mam_idx(m,0)
               dact   = raercol_cw(k,mm,nsav)*dumc
               raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
               raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               do l = 1, nspec_amode(m)
                  mm = mam_idx(m,l)
                  dact    = raercol_cw(k,mm,nsav)*dumc
                  raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               end do
            end do
         end if

         ! shrinking liquid cloud ......................................................
         !    treat the reduction of cloud fraction from when cldn(i,k) < cldo(i,k)
         !    and also dissipate the portion of the cloud that will be regenerated
         cldo_tmp = lcldo(i,k)
         cldn_tmp = lcldn(i,k) * exp( -dtmicro/tau_cld_regenerate )
         !    alternate formulation
         !    cldn_tmp = cldn(i,k) * max( 0.0_r8, (1.0_r8-dtmicro/tau_cld_regenerate) )

         ! fraction is also provided. 
         if (cldn_tmp < cldo_tmp) then
            !  droplet loss in decaying cloud
            !++ sungsup
            nsource(i,k) = nsource(i,k) + qcld(k)*(cldn_tmp - cldo_tmp)/cldo_tmp*cldliqf(i,k)*dtinv
            qcld(k)      = qcld(k)*(1._r8 + (cldn_tmp - cldo_tmp)/cldo_tmp)
            !-- sungsup

            ! convert activated aerosol to interstitial in decaying cloud

            dumc = (cldn_tmp - cldo_tmp)/cldo_tmp * cldliqf(i,k)
            do m = 1, ntot_amode
               mm = mam_idx(m,0)
               dact   = raercol_cw(k,mm,nsav)*dumc
               raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact   ! cloud-borne aerosol
               raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               do l = 1, nspec_amode(m)
                  mm = mam_idx(m,l)
                  dact    = raercol_cw(k,mm,nsav)*dumc
                  raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               end do
            end do
         end if

         ! growing liquid cloud ......................................................
         !    treat the increase of cloud fraction from when cldn(i,k) > cldo(i,k)
         !    and also regenerate part of the cloud 
         cldo_tmp = cldn_tmp
         cldn_tmp = lcldn(i,k)

         if (cldn_tmp-cldo_tmp > 0.01_r8) then

            ! rce-comment - use wtke at layer centers for new-cloud activation
            wbar  = wtke_cen(i,k)
            wmix  = 0._r8
            wmin  = 0._r8
            wmax  = 10._r8
            wdiab = 0._r8

            ! load aerosol properties, assuming external mixtures

            phase = 1 ! interstitial
            do m = 1, ntot_amode
!#if (defined IRONPROC)
           call loadaer( &
              state, pbuf, i, i, k, &
              m, cs, phase, na, va, &
              hy, massSO4, masspORG, masssORG, massORG, &
              massBC, massSS, massDU, massDU1,massDU2,massDU3,massDU4,massDU5,massDU6,massDU7,massDU8)
!#else
!               call loadaer( &
!                  state, pbuf, i, i, k, &
!                  m, cs, phase, na, va, &
!                  hy)
!#endif
               naermod(m)  = na(i)
               vaerosol(m) = va(i)
               hygro(m)    = hy(i)
            end do

            call activate_modal( &
               wbar, wmix, wdiab, wmin, wmax,                       &
               temp(i,k), cs(i,k), naermod, ntot_amode,             &
               vaerosol, hygro, fn, fm, fluxn,     &
               fluxm,flux_fullact(k))

            factnum(i,k,:) = fn

            dumc = (cldn_tmp - cldo_tmp)
            do m = 1, ntot_amode
               mm = mam_idx(m,0)
               dact   = dumc*fn(m)*raer(mm)%fld(i,k) ! interstitial only
               qcld(k) = qcld(k) + dact
               nsource(i,k) = nsource(i,k) + dact*dtinv
               raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
               raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               dum = dumc*fm(m)
               do l = 1, nspec_amode(m)
                  mm = mam_idx(m,l)
                  dact    = dum*raer(mm)%fld(i,k) ! interstitial only
                  raercol_cw(k,mm,nsav) = raercol_cw(k,mm,nsav) + dact  ! cloud-borne aerosol
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) - dact
               enddo
            enddo
         endif

      enddo  ! grow_shrink_main_k_loop
      ! end of k-loop for growing/shrinking cloud calcs ......................

      ! ......................................................................
      ! start of k-loop for calc of old cloud activation tendencies ..........
      !
      ! rce-comment
      !    changed this part of code to use current cloud fraction (cldn) exclusively
      !    consider case of cldo(:)=0, cldn(k)=1, cldn(k+1)=0
      !    previous code (which used cldo below here) would have no cloud-base activation
      !       into layer k.  however, activated particles in k mix out to k+1,
      !       so they are incorrectly depleted with no replacement

      ! old_cloud_main_k_loop
      do k = top_lev, pver
         kp1 = min0(k+1, pver)
         taumix_internal_pver_inv = 0.0_r8

         if (lcldn(i,k) > 0.01_r8) then

            wdiab = 0._r8
            wmix  = 0._r8                       ! single updraft
            wbar  = wtke(i,k)                   ! single updraft
            if (k == pver) wbar = wtke_cen(i,k) ! single updraft
            wmax  = 10._r8
            wmin  = 0._r8

            if (lcldn(i,k) - lcldn(i,kp1) > 0.01_r8 .or. k == pver) then

               ! cloud base

               ! ekd(k) = wtke(i,k)*dz(i,k)/sq2pi
               ! rce-comments
               !   first, should probably have 1/zs(k) here rather than dz(i,k) because
               !      the turbulent flux is proportional to ekd(k)*zs(k),
               !      while the dz(i,k) is used to get flux divergences
               !      and mixing ratio tendency/change
               !   second and more importantly, using a single updraft velocity here
               !      means having monodisperse turbulent updraft and downdrafts.
               !      The sq2pi factor assumes a normal draft spectrum.
               !      The fluxn/fluxm from activate must be consistent with the
               !      fluxes calculated in explmix.
               ekd(k) = wbar/zs(k)

               alogarg = max(1.e-20_r8, 1._r8/lcldn(i,k) - 1._r8)
               wmin    = wbar + wmix*0.25_r8*sq2pi*log(alogarg)
               phase   = 1   ! interstitial

               do m = 1, ntot_amode
                  ! rce-comment - use kp1 here as old-cloud activation involves 
                  !   aerosol from layer below
!#if (defined IRONPROC)
              call loadaer( &
                 state, pbuf, i, i, kp1, &
                 m, cs, phase, na, va, &
                 hy, massSO4, masspORG, masssORG, massORG, &
                 massBC, massSS, massDU, massDU1,massDU2,massDU3,massDU4,massDU5,massDU6,massDU7,massDU8)
!#else
!                  call loadaer( &
!                     state, pbuf, i, i, kp1,  &
!                     m, cs, phase, na, va,   &
!                     hy)
!#endif
                  naermod(m)  = na(i)
                  vaerosol(m) = va(i)
                  hygro(m)    = hy(i)
               end do

               call activate_modal( &
                  wbar, wmix, wdiab, wmin, wmax,                       &
                  temp(i,k), cs(i,k), naermod, ntot_amode,             &
                  vaerosol, hygro,  fn, fm, fluxn,     &
                  fluxm, flux_fullact(k))

               factnum(i,k,:) = fn

               if (k < pver) then
                  dumc = lcldn(i,k) - lcldn(i,kp1)
               else
                  dumc = lcldn(i,k)
               endif

               fluxntot = 0

               ! rce-comment 1
               !    flux of activated mass into layer k (in kg/m2/s)
               !       = "actmassflux" = dumc*fluxm*raercol(kp1,lmass)*csbot(k)
               !    source of activated mass (in kg/kg/s) = flux divergence
               !       = actmassflux/(cs(i,k)*dz(i,k))
               !    so need factor of csbot_cscen = csbot(k)/cs(i,k)
               !                   dum=1./(dz(i,k))
               dum=csbot_cscen(k)/(dz(i,k))

               ! rce-comment 2
               !    code for k=pver was changed to use the following conceptual model
               !    in k=pver, there can be no cloud-base activation unless one considers
               !       a scenario such as the layer being partially cloudy, 
               !       with clear air at bottom and cloudy air at top
               !    assume this scenario, and that the clear/cloudy portions mix with 
               !       a timescale taumix_internal = dz(i,pver)/wtke_cen(i,pver)
               !    in the absence of other sources/sinks, qact (the activated particle 
               !       mixratio) attains a steady state value given by
               !          qact_ss = fcloud*fact*qtot
               !       where fcloud is cloud fraction, fact is activation fraction, 
               !       qtot=qact+qint, qint is interstitial particle mixratio
               !    the activation rate (from mixing within the layer) can now be
               !       written as
               !          d(qact)/dt = (qact_ss - qact)/taumix_internal
               !                     = qtot*(fcloud*fact*wtke/dz) - qact*(wtke/dz)
               !    note that (fcloud*fact*wtke/dz) is equal to the nact/mact
               !    also, d(qact)/dt can be negative.  in the code below
               !       it is forced to be >= 0
               !
               ! steve -- 
               !    you will likely want to change this.  i did not really understand 
               !       what was previously being done in k=pver
               !    in the cam3_5_3 code, wtke(i,pver) appears to be equal to the
               !       droplet deposition velocity which is quite small
               !    in the cam3_5_37 version, wtke is done differently and is much
               !       larger in k=pver, so the activation is stronger there
               !
               if (k == pver) then
                  taumix_internal_pver_inv = flux_fullact(k)/dz(i,k)
               end if

               do m = 1, ntot_amode
                  mm = mam_idx(m,0)
                  fluxn(m) = fluxn(m)*dumc
                  fluxm(m) = fluxm(m)*dumc
                  nact(k,m) = nact(k,m) + fluxn(m)*dum
                  mact(k,m) = mact(k,m) + fluxm(m)*dum
                  if (k < pver) then
                     ! note that kp1 is used here
                     fluxntot = fluxntot &
                        + fluxn(m)*raercol(kp1,mm,nsav)*cs(i,k)
                  else
                     tmpa = raercol(kp1,mm,nsav)*fluxn(m) &
                          + raercol_cw(kp1,mm,nsav)*(fluxn(m) &
                          - taumix_internal_pver_inv*dz(i,k))
                     fluxntot = fluxntot + max(0.0_r8, tmpa)*cs(i,k)
                  end if
               end do
               srcn(k)      = srcn(k) + fluxntot/(cs(i,k)*dz(i,k))
               nsource(i,k) = nsource(i,k) + fluxntot/(cs(i,k)*dz(i,k))

            endif  ! (cldn(i,k) - cldn(i,kp1) > 0.01 .or. k == pver)

         else

            ! no liquid cloud
            nsource(i,k) = nsource(i,k) - qcld(k)*dtinv
            qcld(k)      = 0

            if (cldn(i,k) < 0.01_r8) then
               ! no ice cloud either

               ! convert activated aerosol to interstitial in decaying cloud

               do m = 1, ntot_amode
                  mm = mam_idx(m,0)
                  raercol(k,mm,nsav)    = raercol(k,mm,nsav) + raercol_cw(k,mm,nsav)  ! cloud-borne aerosol
                  raercol_cw(k,mm,nsav) = 0._r8

                  do l = 1, nspec_amode(m)
                     mm = mam_idx(m,l)
                     raercol(k,mm,nsav)    = raercol(k,mm,nsav) + raercol_cw(k,mm,nsav) ! cloud-borne aerosol
                     raercol_cw(k,mm,nsav) = 0._r8
                  end do
               end do
            end if
         end if

      end do  ! old_cloud_main_k_loop

      ! switch nsav, nnew so that nnew is the updated aerosol
      ntemp = nsav
      nsav  = nnew
      nnew  = ntemp

      ! load new droplets in layers above, below clouds

      dtmin     = dtmicro
      ekk(top_lev-1)    = 0.0_r8
      ekk(pver) = 0.0_r8
      do k = top_lev, pver-1
         ! rce-comment -- ekd(k) is eddy-diffusivity at k/k+1 interface
         !   want ekk(k) = ekd(k) * (density at k/k+1 interface)
         !   so use pint(i,k+1) as pint is 1:pverp 
         !           ekk(k)=ekd(k)*2.*pint(i,k)/(rair*(temp(i,k)+temp(i,k+1)))
         !           ekk(k)=ekd(k)*2.*pint(i,k+1)/(rair*(temp(i,k)+temp(i,k+1)))
         ekk(k) = ekd(k)*csbot(k)
      end do

      do k = top_lev, pver
         km1     = max0(k-1, top_lev)
         ekkp(k) = zn(k)*ekk(k)*zs(k)
         ekkm(k) = zn(k)*ekk(k-1)*zs(km1)
         tinv    = ekkp(k) + ekkm(k)

         ! rce-comment -- tinv is the sum of all first-order-loss-rates
         !    for the layer.  for most layers, the activation loss rate
         !    (for interstitial particles) is accounted for by the loss by
         !    turb-transfer to the layer above.
         !    k=pver is special, and the loss rate for activation within 
         !    the layer must be added to tinv.  if not, the time step
         !    can be too big, and explmix can produce negative values.
         !    the negative values are reset to zero, resulting in an 
         !    artificial source.
         if (k == pver) tinv = tinv + taumix_internal_pver_inv

         if (tinv .gt. 1.e-6_r8) then
            dtt   = 1._r8/tinv
            dtmin = min(dtmin, dtt)
         end if
      end do

      dtmix   = 0.9_r8*dtmin
      nsubmix = dtmicro/dtmix + 1
      if (nsubmix > 100) then
         nsubmix_bnd = 100
      else
         nsubmix_bnd = nsubmix
      end if
      count_submix(nsubmix_bnd) = count_submix(nsubmix_bnd) + 1
      dtmix = dtmicro/nsubmix

      do k = top_lev, pver
         kp1 = min(k+1, pver)
         km1 = max(k-1, top_lev)
         ! maximum overlap assumption
         if (cldn(i,kp1) > 1.e-10_r8) then
            overlapp(k) = min(cldn(i,k)/cldn(i,kp1), 1._r8)
         else
            overlapp(k) = 1._r8
         end if
         if (cldn(i,km1) > 1.e-10_r8) then
            overlapm(k) = min(cldn(i,k)/cldn(i,km1), 1._r8)
         else
            overlapm(k) = 1._r8
         end if
      end do


      ! rce-comment
      !    the activation source(k) = mact(k,m)*raercol(kp1,lmass)
      !       should not exceed the rate of transfer of unactivated particles
      !       from kp1 to k which = ekkp(k)*raercol(kp1,lmass)
      !    however it might if things are not "just right" in subr activate
      !    the following is a safety measure to avoid negatives in explmix
      do k = top_lev, pver-1
         do m = 1, ntot_amode
            nact(k,m) = min( nact(k,m), ekkp(k) )
            mact(k,m) = min( mact(k,m), ekkp(k) )
         end do
      end do


      ! old_cloud_nsubmix_loop
      do n = 1, nsubmix
         qncld(:) = qcld(:)
         ! switch nsav, nnew so that nsav is the updated aerosol
         ntemp   = nsav
         nsav    = nnew
         nnew    = ntemp
         srcn(:) = 0.0_r8

         do m = 1, ntot_amode
            mm = mam_idx(m,0)

            ! update droplet source
            ! rce-comment- activation source in layer k involves particles from k+1
            !	       srcn(:)=srcn(:)+nact(:,m)*(raercol(:,mm,nsav))
            srcn(top_lev:pver-1) = srcn(top_lev:pver-1) + nact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))

            ! rce-comment- new formulation for k=pver
            !              srcn(  pver  )=srcn(  pver  )+nact(  pver  ,m)*(raercol(  pver,mm,nsav))
            tmpa = raercol(pver,mm,nsav)*nact(pver,m) &
                 + raercol_cw(pver,mm,nsav)*(nact(pver,m) - taumix_internal_pver_inv)
            srcn(pver) = srcn(pver) + max(0.0_r8,tmpa)
         end do
         call explmix(  &
            qcld, srcn, ekkp, ekkm, overlapp,  &
            overlapm, qncld, zero, zero, pver, &
            dtmix, .false.)

         ! rce-comment
         !    the interstitial particle mixratio is different in clear/cloudy portions
         !    of a layer, and generally higher in the clear portion.  (we have/had
         !    a method for diagnosing the the clear/cloudy mixratios.)  the activation
         !    source terms involve clear air (from below) moving into cloudy air (above).
         !    in theory, the clear-portion mixratio should be used when calculating 
         !    source terms
         do m = 1, ntot_amode
            mm = mam_idx(m,0)
            ! rce-comment -   activation source in layer k involves particles from k+1
            !	              source(:)= nact(:,m)*(raercol(:,mm,nsav))
            source(top_lev:pver-1) = nact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))
            ! rce-comment - new formulation for k=pver
            !               source(  pver  )= nact(  pver,  m)*(raercol(  pver,mm,nsav))
            tmpa = raercol(pver,mm,nsav)*nact(pver,m) &
                 + raercol_cw(pver,mm,nsav)*(nact(pver,m) - taumix_internal_pver_inv)
            source(pver) = max(0.0_r8, tmpa)
            flxconv = 0._r8

            call explmix( &
               raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
               overlapm, raercol_cw(:,mm,nsav), zero, zero, pver,   &
               dtmix, .false.)

            call explmix( &
               raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
               overlapm, raercol(:,mm,nsav), zero, flxconv, pver, &
               dtmix, .true., raercol_cw(:,mm,nsav))

            do l = 1, nspec_amode(m)
               mm = mam_idx(m,l)
               ! rce-comment -   activation source in layer k involves particles from k+1
               !	          source(:)= mact(:,m)*(raercol(:,mm,nsav))
               source(top_lev:pver-1) = mact(top_lev:pver-1,m)*(raercol(top_lev+1:pver,mm,nsav))
               ! rce-comment- new formulation for k=pver
               !                 source(  pver  )= mact(  pver  ,m)*(raercol(  pver,mm,nsav))
               tmpa = raercol(pver,mm,nsav)*mact(pver,m) &
                    + raercol_cw(pver,mm,nsav)*(mact(pver,m) - taumix_internal_pver_inv)
               source(pver) = max(0.0_r8, tmpa)
               flxconv = 0._r8

               call explmix( &
                  raercol_cw(:,mm,nnew), source, ekkp, ekkm, overlapp, &
                  overlapm, raercol_cw(:,mm,nsav), zero, zero, pver,   &
                  dtmix, .false.)

               call explmix( &
                  raercol(:,mm,nnew), source, ekkp, ekkm, overlapp,  &
                  overlapm, raercol(:,mm,nsav), zero, flxconv, pver, &
                  dtmix, .true., raercol_cw(:,mm,nsav))

            end do
         end do

         if (called_from_spcam) then
         !
         ! turbulent mixing for gas species .
         !
               do m=1, pcnst
                  if (cnst_species_class(m) == cnst_spec_class_gas) then
                    flxconv = 0.0_r8
                    zerogas(:) = 0.0_r8
                    call explmix(rgascol(1,m,nnew),zerogas,ekkp,ekkm,overlapp,overlapm,  &
                                 rgascol(1,m,nsav),zero, flxconv, pver,dtmix,&
                                   .true., zerogas)
                  end if
               end do
         endif

      end do ! old_cloud_nsubmix_loop

      ! evaporate particles again if no cloud (either ice or liquid)

      do k = top_lev, pver
         if (cldn(i,k) == 0._r8) then
            ! no ice or liquid cloud
            qcld(k)=0._r8

            ! convert activated aerosol to interstitial in decaying cloud
            do m = 1, ntot_amode
               mm = mam_idx(m,0)
               raercol(k,mm,nnew)    = raercol(k,mm,nnew) + raercol_cw(k,mm,nnew)
               raercol_cw(k,mm,nnew) = 0._r8

               do l = 1, nspec_amode(m)
                  mm = mam_idx(m,l)
                  raercol(k,mm,nnew)    = raercol(k,mm,nnew) + raercol_cw(k,mm,nnew)
                  raercol_cw(k,mm,nnew) = 0._r8
               end do
            end do
         end if
      end do

      ! droplet number

      ndropcol(i) = 0._r8
      do k = top_lev, pver
         ndropmix(i,k) = (qcld(k) - ncldwtr(i,k))*dtinv - nsource(i,k)
         tendnd(i,k)   = (max(qcld(k), 1.e-6_r8) - ncldwtr(i,k))*dtinv
         ndropcol(i)   = ndropcol(i) + ncldwtr(i,k)*pdel(i,k)
      end do
      ndropcol(i) = ndropcol(i)/gravit

      if (prog_modal_aero) then

         raertend = 0._r8
         qqcwtend = 0._r8

         do m = 1, ntot_amode
            do l = 0, nspec_amode(m)

               mm   = mam_idx(m,l)
               lptr = mam_cnst_idx(m,l)

               raertend(top_lev:pver) = (raercol(top_lev:pver,mm,nnew) - raer(mm)%fld(i,top_lev:pver))*dtinv
               qqcwtend(top_lev:pver) = (raercol_cw(top_lev:pver,mm,nnew) - qqcw(mm)%fld(i,top_lev:pver))*dtinv

               coltend(i,mm)    = sum( pdel(i,:)*raertend )/gravit
               coltend_cw(i,mm) = sum( pdel(i,:)*qqcwtend )/gravit

               ptend%q(i,:,lptr) = 0.0_r8
               ptend%q(i,top_lev:pver,lptr) = raertend(top_lev:pver)           ! set tendencies for interstitial aerosol
               qqcw(mm)%fld(i,:) = 0.0_r8
               qqcw(mm)%fld(i,top_lev:pver) = raercol_cw(top_lev:pver,mm,nnew) ! update cloud-borne aerosol
            end do
         end do

      end if

      if (called_from_spcam) then
      !
      ! Gas tendency
      !
           do m=1, pcnst
              if (cnst_species_class(m) == cnst_spec_class_gas) then
                ptend%lq(m) = .true.
                ptend%q(i, :, m) = (rgascol(:,m,nnew)-rgas(i,:,m)) * dtinv
              end if
           end do
      endif

   end do  ! overall_main_i_loop
   ! end of main loop over i/longitude ....................................

   call outfld('NDROPCOL', ndropcol, pcols, lchnk)
   call outfld('NDROPSRC', nsource,  pcols, lchnk)
   call outfld('NDROPMIX', ndropmix, pcols, lchnk)
   call outfld('WTKE    ', wtke,     pcols, lchnk)

   if(called_from_spcam) then  
        call outfld('SPLCLOUD  ', cldn    , pcols, lchnk   )
        call outfld('SPKVH     ', kvh     , pcols, lchnk   )
   endif

!#if (defined IRONPROC)
   call ccncalc(state, pbuf, cs, ccn, ccn_ND, kappa, PMtout, PM1out, PM25out, PM25ALL, PM10out, PM10ALL)

   do l = 1, psat
      call outfld(ccn_name(l), ccn(:,:,l), pcols, lchnk)
   enddo

   do l = 1, pdiam
      call outfld(ccn_name_diam(l), ccn_ND(:,:,l), pcols, lchnk)
      call outfld(kap_name(l), kappa(:,:,l), pcols, lchnk)
   enddo

   do l = 1, pPM
      call outfld(tot_name(l), PMtout(:,l), pcols, lchnk)
      call outfld(PM1_name(l), PM1out(:,l), pcols, lchnk)
      call outfld(PM25_name(l), PM25out(:,l), pcols, lchnk)
      call outfld(PM10_name(l), PM10out(:,l), pcols, lchnk)
   enddo

   call outfld('ALLPM25 ', PM25ALL, pcols, lchnk)
   call outfld('cs      ', cs, pcols, lchnk)
   call outfld('dz      ', dz, pcols, lchnk)
   call outfld('gph     ', zm, pcols, lchnk)

!#else
!   call ccncalc(state, pbuf, cs, ccn)
!   do l = 1, psat
!      call outfld(ccn_name(l), ccn(:,:,l), pcols, lchnk)
!   enddo
!#endif

   ! do column tendencies
   if (prog_modal_aero) then
      do m = 1, ntot_amode
         do l = 0, nspec_amode(m)
            mm = mam_idx(m,l)
            call outfld(fieldname(mm),    coltend(:,mm),    pcols, lchnk)
            call outfld(fieldname_cw(mm), coltend_cw(:,mm), pcols, lchnk)
         end do
      end do
   end if

   if(called_from_spcam) then
   !
   ! output column-integrated Gas tendency (this should be zero)
   !
        do m=1, pcnst
            if (cnst_species_class(m) == cnst_spec_class_gas) then
              do i=1, ncol
                 coltendgas(i) = sum( pdel(i,:)*ptend%q(i,:,m) )/gravit
              end do
              fieldnamegas = trim(cnst_name(m)) // '_mixnuc1sp'
              call outfld( trim(fieldnamegas), coltendgas, pcols, lchnk)
            end if
        end do
        deallocate(rgascol, coltendgas)
   end if

   deallocate( &
      nact,       &
      mact,       &
      raer,       &
      qqcw,       &
      raercol,    &
      raercol_cw, &
      coltend,    &
      coltend_cw, &
      naermod,    &
      hygro,      &
      vaerosol,   &
      fn,         &
      fm,         &
      fluxn,      &
      fluxm       )

end subroutine dropmixnuc

!===============================================================================

subroutine explmix( q, src, ekkp, ekkm, overlapp, overlapm, &
   qold, surfrate, flxconv, pver, dt, is_unact, qactold )

   !  explicit integration of droplet/aerosol mixing
   !     with source due to activation/nucleation


   integer, intent(in) :: pver ! number of levels
   real(r8), intent(out) :: q(pver) ! mixing ratio to be updated
   real(r8), intent(in) :: qold(pver) ! mixing ratio from previous time step
   real(r8), intent(in) :: src(pver) ! source due to activation/nucleation (/s)
   real(r8), intent(in) :: ekkp(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
   ! below layer k  (k,k+1 interface)
   real(r8), intent(in) :: ekkm(pver) ! zn*zs*density*diffusivity (kg/m3 m2/s) at interface
   ! above layer k  (k,k+1 interface)
   real(r8), intent(in) :: overlapp(pver) ! cloud overlap below
   real(r8), intent(in) :: overlapm(pver) ! cloud overlap above
   real(r8), intent(in) :: surfrate ! surface exchange rate (/s)
   real(r8), intent(in) :: flxconv ! convergence of flux from surface
   real(r8), intent(in) :: dt ! time step (s)
   logical, intent(in) :: is_unact ! true if this is an unactivated species
   real(r8), intent(in),optional :: qactold(pver)
   ! mixing ratio of ACTIVATED species from previous step
   ! *** this should only be present
   !     if the current species is unactivated number/sfc/mass

   integer k,kp1,km1

   if ( is_unact ) then
      !     the qactold*(1-overlap) terms are resuspension of activated material
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)
         q(k) = qold(k) + dt*( - src(k) + ekkp(k)*(qold(kp1) - qold(k) +       &
            qactold(kp1)*(1.0_r8-overlapp(k)))               &
            + ekkm(k)*(qold(km1) - qold(k) +     &
            qactold(km1)*(1.0_r8-overlapm(k))) )
         !        force to non-negative
         !        if(q(k)<-1.e-30)then
         !           write(iulog,*)'q=',q(k),' in explmix'
         q(k)=max(q(k),0._r8)
         !        endif
      end do

      !     diffusion loss at base of lowest layer
      q(pver)=q(pver)-surfrate*qold(pver)*dt+flxconv*dt
      !        force to non-negative
      !        if(q(pver)<-1.e-30)then
      !           write(iulog,*)'q=',q(pver),' in explmix'
      q(pver)=max(q(pver),0._r8)
      !        endif
   else
      do k=top_lev,pver
         kp1=min(k+1,pver)
         km1=max(k-1,top_lev)
         q(k) = qold(k) + dt*(src(k) + ekkp(k)*(overlapp(k)*qold(kp1)-qold(k)) +      &
            ekkm(k)*(overlapm(k)*qold(km1)-qold(k)) )
         !        force to non-negative
         !        if(q(k)<-1.e-30)then
         !           write(iulog,*)'q=',q(k),' in explmix'
         q(k)=max(q(k),0._r8)
         !        endif
      end do
      !     diffusion loss at base of lowest layer
      q(pver)=q(pver)-surfrate*qold(pver)*dt+flxconv*dt
      !        force to non-negative
      !        if(q(pver)<-1.e-30)then
      !           write(iulog,*)'q=',q(pver),' in explmix'
      q(pver)=max(q(pver),0._r8)

   end if

end subroutine explmix

!===============================================================================

subroutine activate_modal(wbar, sigw, wdiab, wminf, wmaxf, tair, rhoair,  &
   na, nmode, volume, hygro,  &
   fn, fm, fluxn, fluxm, flux_fullact, smax_prescribed, in_cloud_in, smax_f)

   !      calculates number, surface, and mass fraction of aerosols activated as CCN
   !      calculates flux of cloud droplets, surface area, and aerosol mass into cloud
   !      assumes an internal mixture within each of up to nmode multiple aerosol modes
   !      a gaussiam spectrum of updrafts can be treated.

   !      mks units

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.


   !      input

   real(r8), intent(in) :: wbar          ! grid cell mean vertical velocity (m/s)
   real(r8), intent(in) :: sigw          ! subgrid standard deviation of vertical vel (m/s)
   real(r8), intent(in) :: wdiab         ! diabatic vertical velocity (0 if adiabatic)
   real(r8), intent(in) :: wminf         ! minimum updraft velocity for integration (m/s)
   real(r8), intent(in) :: wmaxf         ! maximum updraft velocity for integration (m/s)
   real(r8), intent(in) :: tair          ! air temperature (K)
   real(r8), intent(in) :: rhoair        ! air density (kg/m3)
   real(r8), intent(in) :: na(:)      ! aerosol number concentration (/m3)
   integer,  intent(in) :: nmode      ! number of aerosol modes
   real(r8), intent(in) :: volume(:)  ! aerosol volume concentration (m3/m3)
   real(r8), intent(in) :: hygro(:)   ! hygroscopicity of aerosol mode

   !      output

   real(r8), intent(out) :: fn(:)      ! number fraction of aerosols activated
   real(r8), intent(out) :: fm(:)      ! mass fraction of aerosols activated
   real(r8), intent(out) :: fluxn(:)   ! flux of activated aerosol number fraction into cloud (cm/s)
   real(r8), intent(out) :: fluxm(:)   ! flux of activated aerosol mass fraction into cloud (cm/s)
   real(r8), intent(out) :: flux_fullact   ! flux of activated aerosol fraction assuming 100% activation (cm/s)
   !    rce-comment
   !    used for consistency check -- this should match (ekd(k)*zs(k))
   !    also, fluxm/flux_fullact gives fraction of aerosol mass flux
   !       that is activated
  
   !      optional
   real(r8), optional, intent(in) :: smax_prescribed  ! prescribed max. supersaturation for secondary activation
   logical,  optional, intent(in) :: in_cloud_in      ! switch to modify calculations when above cloud base
   real(r8), optional, intent(in) :: smax_f           ! droplet and rain size distr factor in the smax calculation
                                                      ! used when in_cloud=.true.

   !      local

   integer, parameter:: nx=200
   integer iquasisect_option, isectional
   real(r8) integ,integf
   real(r8), parameter :: p0 = 1013.25e2_r8    ! reference pressure (Pa)
   real(r8) xmin(nmode),xmax(nmode) ! ln(r) at section interfaces
   real(r8) volmin(nmode),volmax(nmode) ! volume at interfaces
   real(r8) tmass ! total aerosol mass concentration (g/cm3)
   real(r8) sign(nmode)    ! geometric standard deviation of size distribution
   real(r8) rm ! number mode radius of aerosol at max supersat (cm)
   real(r8) pres ! pressure (Pa)
   real(r8) path ! mean free path (m)
   real(r8) diff ! diffusivity (m2/s)
   real(r8) conduct ! thermal conductivity (Joule/m/sec/deg)
   real(r8) diff0,conduct0
   real(r8) es ! saturation vapor pressure
   real(r8) qs ! water vapor saturation mixing ratio
   real(r8) dqsdt ! change in qs with temperature
   real(r8) dqsdp ! change in qs with pressure
   real(r8) g ! thermodynamic function (m2/s)
   real(r8) zeta(nmode), eta(nmode)
   real(r8) lnsmax ! ln(smax)
   real(r8) alpha
   real(r8) gamma
   real(r8) beta
   real(r8) sqrtg
   real(r8) :: amcube(nmode) ! cube of dry mode radius (m)
   real(r8) :: smcrit(nmode) ! critical supersatuation for activation
   real(r8) :: lnsm(nmode) ! ln(smcrit)
   real(r8) smc(nmode) ! critical supersaturation for number mode radius
   real(r8) sumflx_fullact
   real(r8) sumflxn(nmode)
   real(r8) sumflxm(nmode)
   real(r8) sumfn(nmode)
   real(r8) sumfm(nmode)
   real(r8) fnold(nmode)   ! number fraction activated
   real(r8) fmold(nmode)   ! mass fraction activated
   real(r8) wold,gold
   real(r8) alogam
   real(r8) rlo,rhi,xint1,xint2,xint3,xint4
   real(r8) wmin,wmax,w,dw,dwmax,dwmin,wnuc,dwnew,wb
   real(r8) dfmin,dfmax,fnew,fold,fnmin,fnbar,fsbar,fmbar
   real(r8) alw,sqrtalw
   real(r8) smax
   real(r8) x,arg
   real(r8) xmincoeff,xcut,volcut,surfcut
   real(r8) z,z1,z2,wf1,wf2,zf1,zf2,gf1,gf2,gf
   real(r8) etafactor1,etafactor2(nmode),etafactor2max
   real(r8) grow
   character(len=*), parameter :: subname='activate_modal'

   logical :: in_cloud
   integer m,n
   !      numerical integration parameters
   real(r8), parameter :: eps=0.3_r8,fmax=0.99_r8,sds=3._r8

   real(r8), parameter :: namin=1.e6_r8   ! minimum aerosol number concentration (/m3)

   integer ndist(nx)  ! accumulates frequency distribution of integration bins required
   data ndist/nx*0/
   save ndist

   if (present(in_cloud_in)) then
      if (.not. present(smax_f)) call endrun('activate_modal error: smax_f must be supplied when in_cloud is used')
      in_cloud = in_cloud_in
   else
      in_cloud = .false.
   end if

   fn(:)=0._r8
   fm(:)=0._r8
   fluxn(:)=0._r8
   fluxm(:)=0._r8
   flux_fullact=0._r8

   if(nmode.eq.1.and.na(1).lt.1.e-20_r8)return

   if(sigw.le.1.e-5_r8.and.wbar.le.0._r8)return

   if ( present( smax_prescribed ) ) then
      if (smax_prescribed <= 0.0_r8) return
   end if

   pres=rair*rhoair*tair
   diff0=0.211e-4_r8*(p0/pres)*(tair/t0)**1.94_r8
   conduct0=(5.69_r8+0.017_r8*(tair-t0))*4.186e2_r8*1.e-5_r8 ! convert to J/m/s/deg
   call qsat(tair, pres, es, qs)
   dqsdt=latvap/(rh2o*tair*tair)*qs
   alpha=gravit*(latvap/(cpair*rh2o*tair*tair)-1._r8/(rair*tair))
   gamma=(1.0_r8+latvap/cpair*dqsdt)/(rhoair*qs)
   etafactor2max=1.e10_r8/(alpha*wmaxf)**1.5_r8 ! this should make eta big if na is very small.

   grow  = 1._r8/(rhoh2o/(diff0*rhoair*qs)  &
           + latvap*rhoh2o/(conduct0*tair)*(latvap/(rh2o*tair) - 1._r8))
   sqrtg = sqrt(grow)
   beta  = 2._r8*pi*rhoh2o*grow*gamma

   do m=1,nmode

      if(volume(m).gt.1.e-39_r8.and.na(m).gt.1.e-39_r8)then
         !            number mode radius (m)
         !           write(iulog,*)'alogsig,volc,na=',alogsig(m),volc(m),na(m)
         amcube(m)=(3._r8*volume(m)/(4._r8*pi*exp45logsig(m)*na(m)))  ! only if variable size dist
         !           growth coefficent Abdul-Razzak & Ghan 1998 eqn 16
         !           should depend on mean radius of mode to account for gas kinetic effects
         !           see Fountoukis and Nenes, JGR2005 and Meskhidze et al., JGR2006
         !           for approriate size to use for effective diffusivity.
         etafactor2(m)=1._r8/(na(m)*beta*sqrtg)
         if(hygro(m).gt.1.e-10_r8)then
            smc(m)=2._r8*aten*sqrt(aten/(27._r8*hygro(m)*amcube(m))) ! only if variable size dist
         else
            smc(m)=100._r8
         endif
         !	    write(iulog,*)'sm,hygro,amcube=',smcrit(m),hygro(m),amcube(m)
      else
         smc(m)=1._r8
         etafactor2(m)=etafactor2max ! this should make eta big if na is very small.
      endif
      lnsm(m)=log(smc(m)) ! only if variable size dist
      !	 write(iulog,'(a,i4,4g12.2)')'m,na,amcube,hygro,sm,lnsm=', &
      !                   m,na(m),amcube(m),hygro(m),sm(m),lnsm(m)
   enddo

   if(sigw.gt.1.e-5_r8)then ! spectrum of updrafts

      wmax=min(wmaxf,wbar+sds*sigw)
      wmin=max(wminf,-wdiab)
      wmin=max(wmin,wbar-sds*sigw)
      w=wmin
      dwmax=eps*sigw
      dw=dwmax
      dfmax=0.2_r8
      dfmin=0.1_r8
      if (wmax <= w) return
      do m=1,nmode
         sumflxn(m)=0._r8
         sumfn(m)=0._r8
         fnold(m)=0._r8
         sumflxm(m)=0._r8
         sumfm(m)=0._r8
         fmold(m)=0._r8
      enddo
      sumflx_fullact=0._r8

      fold=0._r8
      wold=0._r8
      gold=0._r8

      dwmin = min( dwmax, 0.01_r8 )
      do n = 1, nx

100      wnuc=w+wdiab
         !           write(iulog,*)'wnuc=',wnuc
         alw=alpha*wnuc
         sqrtalw=sqrt(alw)
         etafactor1=alw*sqrtalw

         do m=1,nmode
            eta(m)=etafactor1*etafactor2(m)
            zeta(m)=twothird*sqrtalw*aten/sqrtg
         enddo

         if ( present( smax_prescribed ) ) then
            smax = smax_prescribed
         else
            call maxsat(zeta,eta,nmode,smc,smax)
         endif
         !	      write(iulog,*)'w,smax=',w,smax

         lnsmax=log(smax)

         x=twothird*(lnsm(nmode)-lnsmax)/(sq2*alogsig(nmode))
         fnew=0.5_r8*(1._r8-erf(x))


         dwnew = dw
         if(fnew-fold.gt.dfmax.and.n.gt.1)then
            !              reduce updraft increment for greater accuracy in integration
            if (dw .gt. 1.01_r8*dwmin) then
               dw=0.7_r8*dw
               dw=max(dw,dwmin)
               w=wold+dw
               go to 100
            else
               dwnew = dwmin
            endif
         endif

         if(fnew-fold.lt.dfmin)then
            !              increase updraft increment to accelerate integration
            dwnew=min(1.5_r8*dw,dwmax)
         endif
         fold=fnew

         z=(w-wbar)/(sigw*sq2)
         g=exp(-z*z)
         fnmin=1._r8
         xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3

         do m=1,nmode
            !              modal
            x=twothird*(lnsm(m)-lnsmax)/(sq2*alogsig(m))
            fn(m)=0.5_r8*(1._r8-erf(x))
            fnmin=min(fn(m),fnmin)
            !               integration is second order accurate
            !               assumes linear variation of f*g with w
            fnbar=(fn(m)*g+fnold(m)*gold)
            arg=x-1.5_r8*sq2*alogsig(m)
            fm(m)=0.5_r8*(1._r8-erf(arg))
            fmbar=(fm(m)*g+fmold(m)*gold)
            wb=(w+wold)
            if(w.gt.0._r8)then
               sumflxn(m)=sumflxn(m)+sixth*(wb*fnbar           &
                  +(fn(m)*g*w+fnold(m)*gold*wold))*dw
               sumflxm(m)=sumflxm(m)+sixth*(wb*fmbar           &
                  +(fm(m)*g*w+fmold(m)*gold*wold))*dw
            endif
            sumfn(m)=sumfn(m)+0.5_r8*fnbar*dw
            !	       write(iulog,'(a,9g10.2)')'lnsmax,lnsm(m),x,fn(m),fnold(m),g,gold,fnbar,dw=',lnsmax,lnsm(m),x,fn(m),fnold(m),g,gold,fnbar,dw
            fnold(m)=fn(m)
            sumfm(m)=sumfm(m)+0.5_r8*fmbar*dw
            fmold(m)=fm(m)
         enddo
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact &
            + sixth*(wb*(g+gold) + (g*w+gold*wold))*dw
         !            sumg=sumg+0.5_r8*(g+gold)*dw
         gold=g
         wold=w
         dw=dwnew
         if (n > 1 .and. (w > wmax .or. fnmin > fmax)) exit
         w=w+dw
         if (n == nx) then
            write(iulog,*)'do loop is too short in activate'
            write(iulog,*)'wmin=',wmin,' w=',w,' wmax=',wmax,' dw=',dw
            write(iulog,*)'wbar=',wbar,' sigw=',sigw,' wdiab=',wdiab
            write(iulog,*)'wnuc=',wnuc
            write(iulog,*)'na=',(na(m),m=1,nmode)
            write(iulog,*)'fn=',(fn(m),m=1,nmode)
            !   dump all subr parameters to allow testing with standalone code
            !   (build a driver that will read input and call activate)
            write(iulog,*)'wbar,sigw,wdiab,tair,rhoair,nmode='
            write(iulog,*) wbar,sigw,wdiab,tair,rhoair,nmode
            write(iulog,*)'na=',na
            write(iulog,*)'volume=', (volume(m),m=1,nmode)
            write(iulog,*)'hydro='
            write(iulog,*) hygro
            call endrun(subname)
         end if

      enddo

      ndist(n)=ndist(n)+1
      if(w.lt.wmaxf)then

         !            contribution from all updrafts stronger than wmax
         !            assuming constant f (close to fmax)
         wnuc=w+wdiab

         z1=(w-wbar)/(sigw*sq2)
         z2=(wmaxf-wbar)/(sigw*sq2)
         g=exp(-z1*z1)
         integ=sigw*0.5_r8*sq2*sqpi*(erf(z2)-erf(z1))
         !            consider only upward flow into cloud base when estimating flux
         wf1=max(w,zero)
         zf1=(wf1-wbar)/(sigw*sq2)
         gf1=exp(-zf1*zf1)
         wf2=max(wmaxf,zero)
         zf2=(wf2-wbar)/(sigw*sq2)
         gf2=exp(-zf2*zf2)
         gf=(gf1-gf2)
         integf=wbar*sigw*0.5_r8*sq2*sqpi*(erf(zf2)-erf(zf1))+sigw*sigw*gf

         do m=1,nmode
            sumflxn(m)=sumflxn(m)+integf*fn(m)
            sumfn(m)=sumfn(m)+fn(m)*integ
            sumflxm(m)=sumflxm(m)+integf*fm(m)
            sumfm(m)=sumfm(m)+fm(m)*integ
         enddo
         !           same form as sumflxm but replace the fm with 1.0
         sumflx_fullact = sumflx_fullact + integf
         !            sumg=sumg+integ
      endif


      do m=1,nmode
         fn(m)=sumfn(m)/(sq2*sqpi*sigw)
         !            fn(m)=sumfn(m)/(sumg)
         if(fn(m).gt.1.01_r8)then
            write(iulog,*)'fn=',fn(m),' > 1 in activate'
            write(iulog,*)'w,m,na,amcube=',w,m,na(m),amcube(m)
            write(iulog,*)'integ,sumfn,sigw=',integ,sumfn(m),sigw
            call endrun('activate')
         endif
         fluxn(m)=sumflxn(m)/(sq2*sqpi*sigw)
         fm(m)=sumfm(m)/(sq2*sqpi*sigw)
         !            fm(m)=sumfm(m)/(sumg)
         if(fm(m).gt.1.01_r8)then
            write(iulog,*)'fm=',fm(m),' > 1 in activate'
         endif
         fluxm(m)=sumflxm(m)/(sq2*sqpi*sigw)
      enddo
      !        same form as fluxm
      flux_fullact = sumflx_fullact/(sq2*sqpi*sigw)

   else

      !        single updraft
      wnuc=wbar+wdiab

      if(wnuc.gt.0._r8)then

         w=wbar
              
         if(in_cloud) then

            if (smax_f > 0._r8) then
               smax = alpha*w/(2.0_r8*pi*rhoh2o*grow*gamma*smax_f)
            else
               smax = 1.e-20_r8
            end if

         else ! at cloud base
            alw        = alpha*wnuc
            sqrtalw    = sqrt(alw)
            etafactor1 = alw*sqrtalw

            do m = 1, nmode
               eta(m)  = etafactor1*etafactor2(m)
               zeta(m) = twothird*sqrtalw*aten/sqrtg
            end do
            if ( present(smax_prescribed) ) then
               smax = smax_prescribed
            else
               call maxsat(zeta, eta, nmode, smc, smax)
            end if
         end if

         lnsmax=log(smax)
         xmincoeff=alogaten-twothird*(lnsmax-alog2)-alog3


         do m=1,nmode
            !                 modal
            x=twothird*(lnsm(m)-lnsmax)/(sq2*alogsig(m))
            fn(m)=0.5_r8*(1._r8-erf(x))
            arg=x-1.5_r8*sq2*alogsig(m)
            fm(m)=0.5_r8*(1._r8-erf(arg))
            if(wbar.gt.0._r8)then
               fluxn(m)=fn(m)*w
               fluxm(m)=fm(m)*w
            endif
         enddo
         flux_fullact = w
      endif

   endif

end subroutine activate_modal

!===============================================================================

subroutine maxsat(zeta,eta,nmode,smc,smax)

   !      calculates maximum supersaturation for multiple
   !      competing aerosol modes.

   !      Abdul-Razzak and Ghan, A parameterization of aerosol activation.
   !      2. Multiple aerosol types. J. Geophys. Res., 105, 6837-6844.

   integer,  intent(in)  :: nmode ! number of modes
   real(r8), intent(in)  :: smc(nmode) ! critical supersaturation for number mode radius
   real(r8), intent(in)  :: zeta(nmode)
   real(r8), intent(in)  :: eta(nmode)
   real(r8), intent(out) :: smax ! maximum supersaturation
   integer  :: m  ! mode index
   real(r8) :: sum, g1, g2, g1sqrt, g2sqrt

   do m=1,nmode
      if(zeta(m).gt.1.e5_r8*eta(m).or.smc(m)*smc(m).gt.1.e5_r8*eta(m))then
         !            weak forcing. essentially none activated
         smax=1.e-20_r8
      else
         !            significant activation of this mode. calc activation all modes.
         exit
      endif
          ! No significant activation in any mode.  Do nothing.
      if (m == nmode) return

   enddo

   sum=0.0_r8
   do m=1,nmode
      if(eta(m).gt.1.e-20_r8)then
         g1=zeta(m)/eta(m)
         g1sqrt=sqrt(g1)
         g1=g1sqrt*g1
         g2=smc(m)/sqrt(eta(m)+3._r8*zeta(m))
         g2sqrt=sqrt(g2)
         g2=g2sqrt*g2
         sum=sum+(f1(m)*g1+f2(m)*g2)/(smc(m)*smc(m))
      else
         sum=1.e20_r8
      endif
   enddo

   smax=1._r8/sqrt(sum)

end subroutine maxsat

!===============================================================================
!#if (defined IRONPROC)
subroutine ccncalc(state, pbuf, cs, ccn, ccn_ND, kappa, PMtout, PM1out, PM25out, PM25ALL, PM10out, PM10ALL)
!#else
!subroutine ccncalc(state, pbuf, cs, ccn)
!#endif
   ! calculates number concentration of aerosols activated as CCN at
   ! supersaturation supersat.
   ! assumes an internal mixture of a multiple externally-mixed aerosol modes
   ! cgs units

   ! Ghan et al., Atmos. Res., 1993, 198-221.

   ! arguments

   type(physics_state), target, intent(in)    :: state
   type(physics_buffer_desc),   pointer       :: pbuf(:)


   real(r8), intent(in)  :: cs(pcols,pver)       ! air density (kg/m3)
   real(r8), intent(out) :: ccn(pcols,pver,psat) ! number conc of aerosols activated at supersat (#/m3)

!#if (defined IRONPROC)
   real(r8), intent(out) :: ccn_ND(pcols,pver,pdiam) ! number conc of aerosols above diam threshold (#/m3)
   real(r8), intent(out) :: kappa(pcols,pver,pdiam) ! kappa of aerosols above diam threshold (#/m3)
   real(r8), intent(out) :: PMtout(pcols,pPM)
   real(r8), intent(out) :: PM1out(pcols,pPM)
   real(r8), intent(out) :: PM25out(pcols,pPM)
   real(r8), intent(out) :: PM25ALL(pcols)
   real(r8), intent(out) :: PM10out(pcols,pPM)
   real(r8), intent(out) :: PM10ALL(pcols)
!#endif

   ! local

   integer :: lchnk ! chunk index
   integer :: ncol  ! number of columns
   real(r8), pointer :: tair(:,:)     ! air temperature (K)

   real(r8) naerosol(pcols) ! interstit+activated aerosol number conc (/m3)
   real(r8) vaerosol(pcols) ! interstit+activated aerosol volume conc (m3/m3)

   real(r8) amcube(pcols)
   real(r8) super(psat) ! supersaturation
   real(r8), allocatable :: amcubecoef(:)
   real(r8), allocatable :: argfactor(:)
   real(r8), allocatable :: Nargfactor(:) !as argfactor but for diam calc   
   real(r8), allocatable :: lnsig(:)
   real(r8) :: surften       ! surface tension of water w/respect to air (N/m)
   real(r8) surften_coef
   real(r8) a(pcols) ! surface tension parameter
   real(r8) hygro(pcols)  ! aerosol hygroscopicity
   real(r8) sm(pcols)  ! critical supersaturation at mode radius
   real(r8) arg(pcols)

!#if (defined IRONPROC)
   real(r8) rbar(pcols)  ! mode radius
   real(r8) arghigh(pcols)
   real(r8) arglow(pcols)
   real(r8) argPM1(pcols)
   real(r8) argPM25(pcols)
   real(r8) argPM10(pcols)
   real(r8) Nvol(pcols)
   real(r8) Nvolhigh(pcols)
   real(r8) Nvollow(pcols)
   real(r8) ka(pcols,pver,pdiam)
   real(r8) kb(pcols,pver,pdiam)
   real(r8) zz(pcols)

   real(r8) massSO4(pcols,pver,ntot_amode)
   real(r8) masspORG(pcols,pver,ntot_amode)
   real(r8) masssORG(pcols,pver,ntot_amode)
   real(r8) massORG(pcols,pver,ntot_amode)
   real(r8) massBC(pcols,pver,ntot_amode)
   real(r8) massSS(pcols,pver,ntot_amode)
   real(r8) massDU1(pcols,pver,ntot_amode)
   real(r8) massDU2(pcols,pver,ntot_amode)
   real(r8) massDU3(pcols,pver,ntot_amode)
   real(r8) massDU4(pcols,pver,ntot_amode)
   real(r8) massDU5(pcols,pver,ntot_amode)
   real(r8) massDU6(pcols,pver,ntot_amode)
   real(r8) massDU7(pcols,pver,ntot_amode)
   real(r8) massDU8(pcols,pver,ntot_amode)
   real(r8) massDU(pcols,pver,ntot_amode)
   !real(r8) massFEDU(pcols,pver,ntot_amode)
   !real(r8) massFEAN(pcols,pver,ntot_amode)
   !real(r8) massFEBB(pcols,pver,ntot_amode)
   real(r8) DUSTTOT(pcols,pver,ntot_amode)
   !real(r8) FEDUTOT(pcols,pver,ntot_amode)
   !real(r8) FEDUSOL(pcols,pver,ntot_amode)
   !real(r8) FEANTOT(pcols,pver,ntot_amode)
   !real(r8) FEANSOL(pcols,pver,ntot_amode)
   !real(r8) FEBBTOT(pcols,pver,ntot_amode)
   !real(r8) FEBBSOL(pcols,pver,ntot_amode)
   !real(r8) FETOT(pcols,pver,ntot_amode)
   !real(r8) FESOL(pcols,pver,ntot_amode)
   real(r8) PM25sum(pcols,pver)
   real(r8) PM10sum(pcols,pver)
   real(r8) PMtsum(pcols,pver)
   real(r8) PMt(pcols,pver,pPM)
   real(r8) PM1(pcols,pver,pPM)
   real(r8) PM25(pcols,pver,pPM)
   real(r8) PM10(pcols,pver,pPM)
   real(r8) totdu(pcols)
   real(r8) soldu(pcols)
   real(r8) totco(pcols)
   real(r8) solco(pcols)
   real(r8) totbb(pcols)
   real(r8) solbb(pcols)
   real(r8) tot(pcols)
   real(r8) sol(pcols)

   character(len=8) :: PM_name(pPM)= &
     (/ 'SO4PM','pORGPM','sORGPM','ORGPM','BCPM','SSPM','DUPM'/)
!#endif
   !     mathematical constants
   real(r8) onethird,twothird,sq2
   integer l,m,n,i,k
   real(r8) log,cc
   real(r8) smcoefcoef,smcoef(pcols)
   integer phase ! phase of aerosol
   !-------------------------------------------------------------------------------

   lchnk = state%lchnk
   ncol  = state%ncol
   tair  => state%t

   allocate( &
      amcubecoef(ntot_amode), &
      argfactor(ntot_amode),  &
      Nargfactor(ntot_amode), &   
      lnsig(ntot_amode)       )

   super(:)=supersat(:)*0.01_r8
   sq2=sqrt(2._r8)
   onethird=1._r8/3._r8
   twothird=2._r8/3._r8
   surften=0.076_r8
   surften_coef=2._r8*mwh2o*surften/(r_universal*rhoh2o)
   smcoefcoef=2._r8/sqrt(27._r8)

   do m=1,ntot_amode
      amcubecoef(m)=3._r8/(4._r8*pi*exp45logsig(m))
      argfactor(m)=twothird/(sq2*alogsig(m))
      Nargfactor(m)=1._r8/(sq2*alogsig(m))
      lnsig(m)=(3._r8*sq2/2._r8)*alogsig(m)
   end do

   ccn = 0._r8
!#if (defined IRONPROC)
   ccn_ND  = 0._r8
   ka      = 0._r8
   kb      = 0._r8
   massSO4 = 0._r8
   massBC  = 0._r8
   massORG = 0._r8
   masspORG= 0._r8
   masssORG= 0._r8
   massSS  = 0._r8
   massDU1 = 0._r8
   massDU2 = 0._r8
   massDU3 = 0._r8
   massDU4 = 0._r8
   massDU5 = 0._r8
   massDU6 = 0._r8
   massDU7 = 0._r8
   massDU8 = 0._r8
   massDU  = 0._r8
   !massFEDU= 0._r8
   !massFEAN= 0._r8
   !massFEBB= 0._r8
   PM1     = 0._r8
   PM25    = 0._r8
   PM10    = 0._r8
   PM10sum = 0._r8
   PMt     = 0._r8
   PMtsum  = 0._r8
   PM1out  = 0._r8
   PM25out = 0._r8
   PM10out = 0._r8
   totdu   = 0._r8
   soldu   = 0._r8
   totco   = 0._r8
   solco   = 0._r8
   totbb   = 0._r8
   solbb   = 0._r8
   tot     = 0._r8
   sol     = 0._r8
!#endif

   do k=top_lev,pver

      do i=1,ncol
         a(i)=surften_coef/tair(i,k)
         smcoef(i)=smcoefcoef*a(i)*sqrt(a(i))
      end do

      do m=1,ntot_amode

         phase=3 ! interstitial+cloudborne

!#if (defined IRONPROC)
         call loadaer( &
            state, pbuf, 1, ncol, k, &
            m, cs, phase, naerosol, vaerosol, &
            hygro, massSO4, masspORG, masssORG, massORG, &
            massBC, massSS, massDU, massDU1,massDU2,massDU3,massDU4,massDU5,massDU6,massDU7,massDU8)
!#else
!         call loadaer( &
!            state, pbuf, 1, ncol, k, &
!            m, cs, phase, naerosol, vaerosol, &
!            hygro)
!#endif

!#if (defined IRONPROC)

       ! Longlei Li
       if (m==1) then
            where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) ! critical supersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=0.164E-6 !a/a value should not matter
               rbar(:ncol)=0.034E-6 !a/a value should not matter
            endwhere
         elseif (m==2) then
            where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) !critical supersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=0.034E-6 !values taken as the mid-point
               rbar(:ncol)=0.164E-6 !values taken as the mid-point
            endwhere
         elseif (m==3) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!critical supersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=2.225E-6 !in global annual mean size range
               rbar(:ncol)=0.164E-6 !in global annual mean size range
            endwhere
         elseif (m==4) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=0.164E-6 !see: X.Lie et al. 2012 MAM3 GMD
               rbar(:ncol)=2.225E-6 !see: X.Lie et al. 2012 MAM3 GMD
            endwhere
         elseif (m==5) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=0.164E-6 !see: X.Lie et al. 2012 MAM3 GMD
               rbar(:ncol)=0.50E-6 !see: L.Li et al. 2024 MAM10 JAMES; Drakaki et al. 2022 ACP 
            endwhere
         elseif (m==6) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               rbar(:ncol)=2.225E-6 !see: X.Lie et al. 2012 MAM3 GMD
            endwhere
         elseif (m==7) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=2.225E-6 !see: X.Lie et al. 2012 MAM3 GMD
               rbar(:ncol)=1.25E-6 !see: L.Li et al. 2024 MAM10 JAMES; Drakaki et al. 2022 ACP  
            endwhere
         elseif (m==8) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=5.225E-6 !see: L.Li et al. 2024 MAM10, JAMES 
               rbar(:ncol)=3.5E-6 !see: L.Li et al. 2024 MAM10 JAMES; Drakaki et al. 2022 ACP  
            endwhere
         elseif (m==9) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               !rbar(:ncol)=10.225E-6 !see: L.Li et al. 2024 MAM10, JAMES
               rbar(:ncol)=11.0E-6 !see: L.Li et al. 2024 MAM10 JAMES; Drakaki et al. 2022 ACP  
            endwhere
         elseif (m==10) then
             where(naerosol(:ncol)>1.e-3_r8)
               amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
               sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol))!criticalsupersaturation
               rbar(:ncol)=(amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol))**onethird
            elsewhere
               sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
               rbar(:ncol)=25.0E-6 !see: L.Li et al. 2024 MAM10 JAMES; Drakaki et al. 2022 ACP  
            endwhere
         endif

         do l=1,pdiam
            do i=1,ncol
               arg(i)=Nargfactor(m)*log((diam(l)*0.5_r8)/rbar(i))
               arghigh(i)=Nargfactor(m)*log(((diam(l)+5.0E-9_r8)*0.5_r8)/rbar(i))
               arglow(i)=Nargfactor(m)*log(((diam(l)-5.0E-9_r8)*0.5_r8)/rbar(i))
               Nvolhigh(i)=Nvolhigh(i)+vaerosol(i)*0.5_r8*(1._r8-erf(arghigh(i)-lnsig(m))) !vol conc at high N
               Nvollow(i)=Nvollow(i)+vaerosol(i)*0.5_r8*(1._r8-erf(arglow(i)-lnsig(m))) !vol conc at low N
               Nvol(i)=Nvolhigh(i)-Nvollow(i)   ! difference of high-low N for calc hygro @ N(diam) 
               ka(i,k,l)=ka(i,k,l)+Nvol(i)*hygro(i)
               kb(i,k,l)=kb(i,k,l)+Nvol(i)
               ccn_ND(i,k,l)=ccn_ND(i,k,l)+naerosol(i)*0.5_r8*(1._r8-erf(arg(i)))
            enddo
         enddo

         do l=1,pPM  !species to hold PM mass > cut off (NB:PM# is then totalPM-this mass)
            do i=1,ncol
               argPM1 (i) = Nargfactor(m)*log((1.e-6_r8*0.5_r8)/rbar(i))
               !argPM25(i) = Nargfactor(m)*log((2.5e-6_r8*0.5_r8)/rbar(i))
               argPM25(i) = Nargfactor(m)*log((1.6e-6_r8*0.5_r8)/rbar(i))
               argPM10(i) = Nargfactor(m)*log((6.9e-6_r8*0.5_r8)/rbar(i))
               if (PM_name(l) == 'SO4PM') then
                 PM1 (i,k,l)=PM1 (i,k,l)+massSO4(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+massSO4(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+massSO4(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
                 PMt (i,k,l)=PMt (i,k,l)+massSO4(i,k,m)
               elseif (PM_name(l) == 'pORGPM') then
                 PM1 (i,k,l)=PM1 (i,k,l)+masspORG(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+masspORG(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+masspORG(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
                 PMt (i,k,l)=PMt (i,k,l)+masspORG(i,k,m)
               elseif (PM_name(l) == 'sORGPM') then
                 PM1 (i,k,l)=PM1 (i,k,l)+masssORG(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+masssORG(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+masssORG(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
                 PMt (i,k,l)=PMt (i,k,l)+masssORG(i,k,m)
               elseif (PM_name(l) == 'ORGPM') then    ! no PMt as already counted in prim and sec org above
                 PM1 (i,k,l)=PM1 (i,k,l)+massORG(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+massORG(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+massORG(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
               elseif (PM_name(l) == 'BCPM') then
                 PM1 (i,k,l)=PM1 (i,k,l)+massBC(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+massBC(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+massBC(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
                 PMt (i,k,l)=PMt (i,k,l)+massBC(i,k,m)
               elseif (PM_name(l) == 'SSPM') then
                 PM1 (i,k,l)=PM1 (i,k,l)+massSS(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+massSS(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+massSS(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
                 PMt (i,k,l)=PMt (i,k,l)+massSS(i,k,m)
               elseif (PM_name(l) == 'DUPM') then
                 PM1 (i,k,l)=PM1 (i,k,l)+massDU(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
                 PM25(i,k,l)=PM25(i,k,l)+massDU(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
                 PM10(i,k,l)=PM10(i,k,l)+massDU(i,k,m)*0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
                 PMt (i,k,l)=PMt (i,k,l)+massDU(i,k,m)
               !elseif (PM_name(l) == 'FEDUPM') then
               !  PM1 (i,k,l)=PM1 (i,k,l)+massFEDU(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
               !  PM25(i,k,l)=PM25(i,k,l)+massFEDU(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
               !  PMt (i,k,l)=PMt (i,k,l)+massFEDU(i,k,m)
               !elseif (PM_name(l) == 'FEANPM') then
               !  PM1 (i,k,l)=PM1 (i,k,l)+massFEAN(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
               !  PM25(i,k,l)=PM25(i,k,l)+massFEAN(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
               !  PMt (i,k,l)=PMt (i,k,l)+massFEAN(i,k,m)
               !elseif (PM_name(l) == 'FEBBPM') then
               !  PM1 (i,k,l)=PM1 (i,k,l)+massFEBB(i,k,m)*0.5_r8*(1._r8-erf(argPM1(i)-lnsig(m)))
               !  PM25(i,k,l)=PM25(i,k,l)+massFEBB(i,k,m)*0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
               !  PMt (i,k,l)=PMt (i,k,l)+massFEBB(i,k,m)
               endif
            enddo
         enddo

         !do i=1,ncol
         !  PMtsum(i,k) =PMtsum(i,k)+(massSO4(i,k,m)+massORG(i,k,m)+massBC(i,k,m)+ &
         !                            massSS(i,k,m)+massDU(i,k,m)+massFEDU(i,k,m)+ &
         !                            massFEAN(i,k,m)+massFEBB(i,k,m))
         !  PM25sum(i,k)=PM25sum(i,k)+(massSO4(i,k,m)+massORG(i,k,m)+massBC(i,k,m)+ &
         !                             massSS(i,k,m)+massDU(i,k,m)+massFEDU(i,k,m)+ &
         !                             massFEAN(i,k,m)+massFEBB(i,k,m)) &
         !                             *0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))
         !enddo

         do i=1,ncol
           PMtsum(i,k) =PMtsum(i,k)+(massSO4(i,k,m)+massORG(i,k,m)+massBC(i,k,m)+ &
                                     massSS(i,k,m)+massDU(i,k,m) &
                                    )
           PM25sum(i,k)=PM25sum(i,k)+(massSO4(i,k,m)+massORG(i,k,m)+massBC(i,k,m)+ &
                                      massSS(i,k,m)+massDU(i,k,m) &
                                     ) &
                                      *0.5_r8*(1._r8-erf(argPM25(i)-lnsig(m)))

           PM10sum(i,k)=PM10sum(i,k)+(massSO4(i,k,m)+massORG(i,k,m)+massBC(i,k,m)+ &
                                      massSS(i,k,m)+massDU(i,k,m) &
                                     ) &
                                      *0.5_r8*(1._r8-erf(argPM10(i)-lnsig(m)))
         enddo
!#else
         !where(naerosol(:ncol)>1.e-3_r8)
         !   amcube(:ncol)=amcubecoef(m)*vaerosol(:ncol)/naerosol(:ncol)
         !   sm(:ncol)=smcoef(:ncol)/sqrt(hygro(:ncol)*amcube(:ncol)) ! critical supersaturation
         !elsewhere
         !   sm(:ncol)=1._r8 ! value shouldn't matter much since naerosol is small
         !endwhere
!#endif
         do l=1,psat
            do i=1,ncol
               arg(i)=argfactor(m)*log(sm(i)/super(l))
               ccn(i,k,l)=ccn(i,k,l)+naerosol(i)*0.5_r8*(1._r8-erf(arg(i)))
            enddo
         enddo
      enddo
   enddo
   ccn(:ncol,:,:)=ccn(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3

!#if (defined IRONPROC)
   ccn_ND(:ncol,:,:) = ccn_ND(:ncol,:,:)*1.e-6_r8 ! convert from #/m3 to #/cm3
   PM1out(:ncol,:)   = (PMt(:ncol,pver,:)-PM1(:ncol,pver,:))*1.e9_r8  !kg/m3 -> ug/m3
   PM25out(:ncol,:)  = (PMt(:ncol,pver,:)-PM25(:ncol,pver,:))*1.e9_r8  !kg/m3 -> ug/m3
   PM10out(:ncol,:)  = (PMt(:ncol,pver,:)-PM10(:ncol,pver,:))*1.e9_r8  !kg/m3 -> ug/m3
   PMtout(:ncol,:)   = PMt(:ncol,pver,:)*1.e9_r8  !kg/m3 -> ug/m3
   PM25ALL(:ncol)    = (PMtsum(:ncol,pver)-PM25sum(:ncol,pver))*1.e9_r8  !kg/m3 -> ug/m3
   PM10ALL(:ncol)    = (PMtsum(:ncol,pver)-PM10sum(:ncol,pver))*1.e9_r8  !kg/m3 -> ug/m3
   kappa = 0._r8
   kappa = ka/kb  ! kappa at each N
!#endif
   deallocate( &
      amcubecoef, &
      argfactor , &
      Nargfactor, &
      lnsig       )

end subroutine ccncalc

!===============================================================================
!#if (defined IRONPROC)
subroutine loadaer( &
   state, pbuf, istart, istop, k, &
   m, cs, phase, naerosol, &
   vaerosol, hygro, massSO4, masspORG, masssORG, massORG, &
   massBC, massSS, massDU, massDU1,massDU2,massDU3,massDU4,massDU5,massDU6,massDU7,massDU8)
!#else
!subroutine loadaer( &
!   state, pbuf, istart, istop, k, &
!   m, cs, phase, naerosol, &
!   vaerosol, hygro)
!#endif
   ! return aerosol number, volume concentrations, and bulk hygroscopicity

   ! input arguments
   type(physics_state), target, intent(in) :: state
   type(physics_buffer_desc),   pointer    :: pbuf(:)

   integer,  intent(in) :: istart      ! start column index (1 <= istart <= istop <= pcols)
   integer,  intent(in) :: istop       ! stop column index  
   integer,  intent(in) :: m           ! mode index
   integer,  intent(in) :: k           ! level index
   real(r8), intent(in) :: cs(:,:)     ! air density (kg/m3)
   integer,  intent(in) :: phase       ! phase of aerosol: 1 for interstitial, 2 for cloud-borne, 3 for sum

   ! output arguments
   real(r8), intent(out) :: naerosol(:)  ! number conc (1/m3)
   real(r8), intent(out) :: vaerosol(:)  ! volume conc (m3/m3)
   real(r8), intent(out) :: hygro(:)     ! bulk hygroscopicity of mode

!#if (defined IRONPROC)
! dsh -mass conc of each species- for PM calcs
   real(r8), intent(inout) :: massSO4(:,:,:)
   real(r8), intent(inout) :: masspORG(:,:,:)
   real(r8), intent(inout) :: masssORG(:,:,:)
   real(r8), intent(inout) :: massORG(:,:,:)
   real(r8), intent(inout) :: massBC(:,:,:)
   real(r8), intent(inout) :: massSS(:,:,:)
   real(r8), intent(inout) :: massDU(:,:,:)
   real(r8), intent(inout) :: massDU1(:,:,:)
   real(r8), intent(inout) :: massDU2(:,:,:)
   real(r8), intent(inout) :: massDU3(:,:,:)
   real(r8), intent(inout) :: massDU4(:,:,:)
   real(r8), intent(inout) :: massDU5(:,:,:)
   real(r8), intent(inout) :: massDU6(:,:,:)
   real(r8), intent(inout) :: massDU7(:,:,:)
   real(r8), intent(inout) :: massDU8(:,:,:)

   real(r8), pointer :: specmmr_a(:,:) !mmr of species - interstitial
   real(r8), pointer :: specmmr_c(:,:) !mmr of species -cloud bourne
   character*32      :: spectype
!#endif

   ! internal
   integer  :: lchnk               ! chunk identifier

   real(r8), pointer :: raer(:,:) ! interstitial aerosol mass, number mixing ratios
   real(r8), pointer :: qqcw(:,:) ! cloud-borne aerosol mass, number mixing ratios
   real(r8) :: specdens, spechygro

   real(r8) :: vol(pcols) ! aerosol volume mixing ratio
   integer  :: i, l
   !-------------------------------------------------------------------------------

   lchnk = state%lchnk

   do i = istart, istop
      vaerosol(i) = 0._r8
      hygro(i)    = 0._r8
!#if (defined IRONPROC)
      massSO4(i,k,m)  = 0._r8
      masspORG(i,k,m) = 0._r8
      masssORG(i,k,m) = 0._r8
      massORG(i,k,m)  = 0._r8
      massBC(i,k,m)   = 0._r8
      massSS(i,k,m)   = 0._r8
      massDU(i,k,m)   = 0._r8
      massDU1(i,k,m)  = 0._r8
      massDU2(i,k,m)  = 0._r8
      massDU3(i,k,m)  = 0._r8
      massDU4(i,k,m)  = 0._r8
      massDU5(i,k,m)  = 0._r8
      massDU6(i,k,m)  = 0._r8
      massDU7(i,k,m)  = 0._r8
      massDU8(i,k,m)  = 0._r8
!#endif
   end do

   do l = 1, nspec_amode(m)

     call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, raer)
     call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, qqcw)
!#if (defined IRONPROC)
     call rad_cnst_get_aer_props(0, m, l, density_aer=specdens, hygro_aer=spechygro, spectype=spectype)
!#else
!     call rad_cnst_get_aer_props(0, m, l, density_aer=specdens, hygro_aer=spechygro)
!#endif

     if (phase == 3) then
!#if (defined IRONPROC)

!total mmr (cloud and inter)
      call rad_cnst_get_aer_mmr(0, m, l, 'a', state, pbuf, specmmr_a)
      call rad_cnst_get_aer_mmr(0, m, l, 'c', state, pbuf, specmmr_c)
 
      if (trim(spectype) == 'sulfate') then
         do i = istart, istop
            massSO4(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
         enddo
      elseif (trim(spectype) == 'p-organic') then
         do i = istart, istop
            masspORG(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massORG(i,k,m)=massORG(i,k,m) + masspORG(i,k,m)
         enddo
      elseif (trim(spectype) == 's-organic') then
         do i = istart, istop
            masssORG(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massORG(i,k,m)=massORG(i,k,m) + masssORG(i,k,m)
         enddo
      elseif (trim(spectype) == 'black-c') then
         do i = istart, istop
            massBC(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
         enddo
      elseif (trim(spectype) == 'seasalt') then
         do i = istart, istop
            massSS(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
         enddo
      !elseif (trim(spectype) == 'dust') then
      !   do i = istart, istop
      !      massDU(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !      massDU(i,k,m)=massDU(i,k,m) + massDU(i,k,m)
      !   enddo
      elseif (trim(spectype) == 'dust-type1') then
         do i = istart, istop
            massDU1(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU1(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type2') then
         do i = istart, istop
            massDU2(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU2(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type3') then
         do i = istart, istop
            massDU3(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU3(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type4') then
         do i = istart, istop
            massDU4(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU4(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type5') then
         do i = istart, istop
            massDU5(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU5(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type6') then
         do i = istart, istop
            massDU6(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU6(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type7') then
         do i = istart, istop
            massDU7(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU7(i,k,m)
         enddo
      elseif (trim(spectype) == 'dust-type8') then
         do i = istart, istop
            massDU8(i,k,m)=max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
            massDU(i,k,m)=massDU(i,k,m) + massDU8(i,k,m)
         enddo
      !elseif (trim(spectype) == 'fess') then
      !   do i = istart, istop
      !      massFEDU(i,k,m)=massFEDU(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fesm') then
      !   do i = istart, istop
      !      massFEDU(i,k,m)=massFEDU(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fets') then
      !   do i = istart, istop
      !      massFEDU(i,k,m)=massFEDU(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fetm') then
      !   do i = istart, istop
      !      massFEDU(i,k,m)=massFEDU(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fesc') then
      !   do i = istart, istop
      !      massFEAN(i,k,m)=massFEAN(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fetc') then
      !   do i = istart, istop
      !      massFEAN(i,k,m)=massFEAN(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fesf') then
      !   do i = istart, istop
      !      massFEBB(i,k,m)=massFEBB(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      !elseif (trim(spectype) == 'fetf') then
      !   do i = istart, istop
      !      massFEBB(i,k,m)=massFEBB(i,k,m) + max(specmmr_a(i,k) + specmmr_c(i,k), 0._r8)*cs(i,k)
      !   enddo
      endif
!#endif 
         do i = istart, istop
            vol(i) = max(raer(i,k) + qqcw(i,k), 0._r8)/specdens
         end do
      else if (phase == 2) then
         do i = istart, istop
            vol(i) = max(qqcw(i,k), 0._r8)/specdens
         end do
      else if (phase == 1) then
         do i = istart, istop
            vol(i) = max(raer(i,k), 0._r8)/specdens
         end do
      else
         write(iulog,*)'phase=',phase,' in loadaer'
         call endrun('phase error in loadaer')
      end if

      do i = istart, istop
         vaerosol(i) = vaerosol(i) + vol(i)
         hygro(i)    = hygro(i) + vol(i)*spechygro
      end do

   end do

   do i = istart, istop
      if (vaerosol(i) > 1.0e-30_r8) then   ! +++xl add 8/2/2007
         hygro(i)    = hygro(i)/(vaerosol(i))
         vaerosol(i) = vaerosol(i)*cs(i,k)
      else
         hygro(i)    = 0.0_r8
         vaerosol(i) = 0.0_r8
      end if
   end do

   ! aerosol number
   call rad_cnst_get_mode_num(0, m, 'a', state, pbuf, raer)
   call rad_cnst_get_mode_num(0, m, 'c', state, pbuf, qqcw)
   if (phase == 3) then
      do i = istart, istop
         naerosol(i) = (raer(i,k) + qqcw(i,k))*cs(i,k)
      end do
   else if (phase == 2) then
      do i = istart, istop
         naerosol(i) = qqcw(i,k)*cs(i,k)
      end do
   else
      do i = istart, istop
         naerosol(i) = raer(i,k)*cs(i,k)
      end do
   end if
   ! adjust number so that dgnumlo < dgnum < dgnumhi
   do i = istart, istop
      naerosol(i) = max(naerosol(i), vaerosol(i)*voltonumbhi_amode(m))
      naerosol(i) = min(naerosol(i), vaerosol(i)*voltonumblo_amode(m))
   end do

end subroutine loadaer

!===============================================================================

end module ndrop




