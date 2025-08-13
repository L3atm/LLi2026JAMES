
      module mo_sim_dat

      private
      public :: set_sim_dat

      contains

      subroutine set_sim_dat

      use chem_mods,     only : clscnt, cls_rxt_cnt, clsmap, permute, adv_mass, fix_mass, crb_mass
      use chem_mods,     only : diag_map
      use chem_mods,     only : phtcnt, rxt_tag_cnt, rxt_tag_lst, rxt_tag_map
      use chem_mods,     only : pht_alias_lst, pht_alias_mult
      use chem_mods,     only : extfrc_lst, inv_lst, slvd_lst
      use chem_mods,     only : enthalpy_cnt, cph_enthalpy, cph_rid, num_rnts, rxntot
      use cam_abortutils,only : endrun
      use mo_tracname,   only : solsym
      use chem_mods,     only : frc_from_dataset
      use chem_mods,     only : is_scalar, is_vector
      use shr_kind_mod,  only : r8 => shr_kind_r8
      use cam_logfile,   only : iulog

      implicit none

!--------------------------------------------------------------
!      ... local variables
!--------------------------------------------------------------
      integer :: ios

      is_scalar = .true.
      is_vector = .false.

      clscnt(:) = (/      0,     0,     0,    46,     0 /)

      cls_rxt_cnt(:,4) = (/      1,     6,     0,    46 /)

      solsym(: 75) = (/ 'H2O2            ','H2SO4           ','SO2             ','DMS             ','SOAG            ', &
                        'so4_a1          ','pom_a1          ','soa_a1          ','bc_a1           ','ncl_a1          ', &
                        'num_a1          ','so4_a2          ','soa_a2          ','ncl_a2          ','num_a2          ', &
                        'pom_a3          ','bc_a3           ','num_a3          ','ncl_a4          ','so4_a4          ', &
                        'num_a4          ','dst1_a5         ','dst2_a5         ','dst3_a5         ','dst4_a5         ', &
                        'dst5_a5         ','dst6_a5         ','dst7_a5         ','dst8_a5         ','so4_a5          ', &
                        'num_a5          ','ncl_a6          ','so4_a6          ','num_a6          ','dst1_a7         ', &
                        'dst2_a7         ','dst3_a7         ','dst4_a7         ','dst5_a7         ','dst6_a7         ', &
                        'dst7_a7         ','dst8_a7         ','so4_a7          ','num_a7          ','dst1_a8         ', &
                        'dst2_a8         ','dst3_a8         ','dst4_a8         ','dst5_a8         ','dst6_a8         ', &
                        'dst7_a8         ','dst8_a8         ','so4_a8          ','num_a8          ','dst1_a9         ', &
                        'dst2_a9         ','dst3_a9         ','dst4_a9         ','dst5_a9         ','dst6_a9         ', &
                        'dst7_a9         ','dst8_a9         ','so4_a9          ','num_a9          ','dst1_a10        ', &
                        'dst2_a10        ','dst3_a10        ','dst4_a10        ','dst5_a10        ','dst6_a10        ', &
                        'dst7_a10        ','dst8_a10        ','so4_a10         ','num_a10         ','H2O             ' /)

      adv_mass(: 75) = (/    34.013600_r8,    98.078400_r8,    64.064800_r8,    62.132400_r8,    12.011000_r8, &
                            115.107340_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,    58.442468_r8, &
                              1.007400_r8,   115.107340_r8,    12.011000_r8,    58.442468_r8,     1.007400_r8, &
                             12.011000_r8,    12.011000_r8,     1.007400_r8,    58.442468_r8,   115.107340_r8, &
                              1.007400_r8,   389.340000_r8,   258.160000_r8,   549.070000_r8,   159.690000_r8, &
                             60.080000_r8,   100.090000_r8,   278.330000_r8,   172.170000_r8,   115.107340_r8, &
                              1.007400_r8,    58.442468_r8,   115.107340_r8,     1.007400_r8,   389.340000_r8, &
                            258.160000_r8,   549.070000_r8,   159.690000_r8,    60.080000_r8,   100.090000_r8, &
                            278.330000_r8,   172.170000_r8,   115.107340_r8,     1.007400_r8,   389.340000_r8, &
                            258.160000_r8,   549.070000_r8,   159.690000_r8,    60.080000_r8,   100.090000_r8, &
                            278.330000_r8,   172.170000_r8,   115.107340_r8,     1.007400_r8,   389.340000_r8, &
                            258.160000_r8,   549.070000_r8,   159.690000_r8,    60.080000_r8,   100.090000_r8, &
                            278.330000_r8,   172.170000_r8,   115.107340_r8,     1.007400_r8,   389.340000_r8, &
                            258.160000_r8,   549.070000_r8,   159.690000_r8,    60.080000_r8,   100.090000_r8, &
                            278.330000_r8,   172.170000_r8,   115.107340_r8,     1.007400_r8,    18.014200_r8 /)

      crb_mass(: 75) = (/     0.000000_r8,     0.000000_r8,     0.000000_r8,    24.022000_r8,    12.011000_r8, &
                              0.000000_r8,    12.011000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8, &
                             12.011000_r8,    12.011000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,    12.011000_r8,    12.011000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8, &
                              0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8,     0.000000_r8 /)

      fix_mass(:  7) = (/ 0.00000000_r8, 28.0134800_r8, 31.9988000_r8, 47.9982000_r8, 17.0068000_r8, &
                          62.0049400_r8, 33.0062000_r8 /)

      clsmap(: 75,4) = (/    1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                            31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                            41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                            51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                            61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                            71,  72,  73,  74,  75 /)

      permute(: 75,4) = (/   1,   2,   3,   4,   5,   6,   7,   8,   9,  10, &
                            11,  12,  13,  14,  15,  16,  17,  18,  19,  20, &
                            21,  22,  23,  24,  25,  26,  27,  28,  29,  30, &
                            31,  32,  33,  34,  35,  36,  37,  38,  39,  40, &
                            41,  42,  43,  44,  45,  46,  47,  48,  49,  50, &
                            51,  52,  53,  54,  55,  56,  57,  58,  59,  60, &
                            61,  62,  63,  64,  65,  66,  67,  68,  69,  70, &
                            71,  72,  73,  74,  75 /) 

      diag_map(: 75) = (/    1,   2,   4,   6,   7,   8,   9,  10,  11,  12, &
                            13,  14,  15,  16,  17,  18,  19,  20,  21,  22, &
                            23,  24,  25,  26,  27,  28,  29,  30,  31,  32, &
                            33,  34,  35,  36,  37,  38,  39,  40,  41,  42, &
                            43,  44,  45,  46,  47,  48,  49,  50,  51,  52, &
                            53,  54,  55,  56,  57,  58,  59,  60,  61,  62, &
                            63,  64,  65,  66,  67,  68,  69,  70,  71,  72, &
                            73,  74,  75,  76,  77/)

      extfrc_lst(:  8) = (/ 'SO2             ','so4_a1          ','so4_a2          ','pom_a3          ','bc_a3           ', &
                            'num_a1          ','num_a2          ','num_a3          '/)

      frc_from_dataset(:  8) = (/ .true., .true., .true., .true., .true., &
                                  .true., .true., .true. /)

      inv_lst(:  7) = (/ 'M               ', 'N2              ', 'O2              ', 'O3              ', 'OH              ', &
                         'NO3             ', 'HO2             ' /)

      if( allocated( rxt_tag_lst ) ) then
         deallocate( rxt_tag_lst )
      end if
      allocate( rxt_tag_lst(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_lst; error = ',ios
         call endrun
      end if
      if( allocated( rxt_tag_map ) ) then
         deallocate( rxt_tag_map )
      end if
      allocate( rxt_tag_map(rxt_tag_cnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate rxt_tag_map; error = ',ios
         call endrun
      end if
      rxt_tag_lst(     1:     4) = (/ 'jh2o2                           ', 'usr_HO2_HO2                     ', &
                                      'usr_SO2_OH                      ', 'usr_DMS_OH                      ' /)
      rxt_tag_map(:rxt_tag_cnt) = (/    1,   2,   4,   6 /)
      if( allocated( pht_alias_lst ) ) then
         deallocate( pht_alias_lst )
      end if
      allocate( pht_alias_lst(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_lst; error = ',ios
         call endrun
      end if
      if( allocated( pht_alias_mult ) ) then
         deallocate( pht_alias_mult )
      end if
      allocate( pht_alias_mult(phtcnt,2),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate pht_alias_mult; error = ',ios
         call endrun
      end if
      pht_alias_lst(:,1) = (/ '                ' /)
      pht_alias_lst(:,2) = (/ '                ' /)
      pht_alias_mult(:,1) = (/ 1._r8 /)
      pht_alias_mult(:,2) = (/ 1._r8 /)
      allocate( num_rnts(rxntot-phtcnt),stat=ios )
      if( ios /= 0 ) then
         write(iulog,*) 'set_sim_dat: failed to allocate num_rnts; error = ',ios
         call endrun
      end if
      num_rnts(:) = (/      2,     2,     2,     2,     2,     2 /)

      end subroutine set_sim_dat

      end module mo_sim_dat
