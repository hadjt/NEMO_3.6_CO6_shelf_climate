MODULE diaregmean 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Timeseries of Regional Means 
   !!======================================================================
   !! History :  3.6  !  11/2016  (J Tinker)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wrk_nemo        ! working arrays
   USE diapea          ! PEA
   USE zdfmxl          ! MLD
   USE sbc_oce
#if defined key_diaar5 
   USE diaar5
#endif


#if defined key_fabm
   USE trc
   USE par_fabm
#endif

   IMPLICIT NONE
   PRIVATE

   LOGICAL , PUBLIC ::   ln_diaregmean  ! region mean calculation
   PUBLIC   dia_regmean_init            ! routine called by nemogcm.F90
   PUBLIC   dia_regmean                 ! routine called by diawri.F90   
   PUBLIC   dia_calctmb_region_mean     ! routine called by diatmb.F90

   
   
   LOGICAL :: ln_diaregmean_ascii       ! region mean calculation ascii output
   LOGICAL :: ln_diaregmean_bin         ! region mean calculation binary output
   LOGICAL :: ln_diaregmean_nc          ! region mean calculation netcdf output
   LOGICAL :: ln_diaregmean_diaar5      ! region mean calculation including AR5 SLR terms
   LOGICAL :: ln_diaregmean_diasbc      ! region mean calculation including Surface BC
   LOGICAL :: ln_diaregmean_karamld     ! region mean calculation including kara mld terms
   LOGICAL :: ln_diaregmean_pea         ! region mean calculation including pea terms


   LOGICAL :: ln_diaregmean_bgc         ! region mean calculation including BGC terms


   
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:,:)  ::   tmp_region_mask_real   ! tempory region_mask of reals
   INTEGER,  SAVE, ALLOCATABLE,   DIMENSION(:,:,:)  ::   region_mask            ! region_mask matrix
   INTEGER                                          ::   nmasks                 ! Number of mask files in region_mask.nc file - 
   INTEGER,  SAVE, ALLOCATABLE,   DIMENSION(:)      ::   nreg_mat               ! Number of regions in each mask
   
   REAL(wp),  ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_field_mat !: temporary region_mask
   REAL(wp),  ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_field_AR5_mat !: temporary region_mask
   REAL(wp),  ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_field_SBC_mat !: temporary region_mask
   INTEGER  ::   tmp_field_cnt                                   ! tmp_field_cnt integer
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_regmean_init 
      !!---------------------------------------------------------------------------
      !!                  ***  ROUTINE dia_regmean_init  ***
      !!     
      !! ** Purpose: Initialization of region mask namelist 
      !!        
      !! ** Method : Read namelist
      !!   History
      !!   3.6  !  11-16  (J Tinker) Routine to initialize dia_regmean
      !!---------------------------------------------------------------------------
      !!
      INTEGER ::   ios                  ! Local integer output status for namelist read
      INTEGER  ::  inum                 ! temporary logical unit ! copied from DOM/domzgr.F90
      INTEGER  ::   ierr                ! error integer for IOM_get
      INTEGER  ::   idmaskvar           ! output of iom_varid
      INTEGER  ::   maskno              ! counter for number of masks
      INTEGER  ::   jj,ji               ! i and j index
      INTEGER  ::   tmpint              ! temporary integer
      REAL(wp),  ALLOCATABLE,   DIMENSION(:,:) ::   tmpregion !: temporary region_mask
      INTEGER, DIMENSION(3) ::   zdimsz   ! number of elements in each of the 3 dimensions (i.e., lon, lat, no of masks, 297,  375,  4) for an array
      INTEGER               ::   zndims   ! number of dimensions in an array (i.e. 3, )


#if defined key_fabm
      INTEGER               ::   js,jl,jn, tmp_dummy

      CHARACTER (len=120) ::    tmp_name,tmp_long_name, tmp_unit

      INTEGER               ::   BGC_nlevs,nBGC_output, bgci
      CHARACTER(len = 10), ALLOCATABLE, DIMENSION(:)       ::   BGC_stat_name(:),BGC_lev_name(:),BGC_output_var(:)
#endif

      !
      NAMELIST/nam_diaregmean/ ln_diaregmean,ln_diaregmean_ascii,ln_diaregmean_bin,ln_diaregmean_nc,&
        & ln_diaregmean_karamld, ln_diaregmean_pea,ln_diaregmean_diaar5,ln_diaregmean_diasbc,ln_diaregmean_bgc
      
      
      ! read in Namelist. 
      !!----------------------------------------------------------------------
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in referdiatmbence namelist : TMB diagnostics
      READ   ( numnam_ref, nam_diaregmean, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaregmean in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_diaregmean, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_diaregmean in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_diaregmean )

      IF(lwp) THEN                   ! Control print
          WRITE(numout,*)
          WRITE(numout,*) 'dia_regmean_init : Output regional mean Diagnostics'
          WRITE(numout,*) '~~~~~~~~~~~~'
          WRITE(numout,*) 'Namelist nam_regmean : set regmeanoutputs '
          WRITE(numout,*) 'Switch for regmean diagnostics (T) or not (F)  ln_diaregmean  = ', ln_diaregmean
          WRITE(numout,*) 'Switch for regmean ascii output (T) or not (F)  ln_diaregmean_ascii  = ', ln_diaregmean_ascii
          WRITE(numout,*) 'Switch for regmean binary output (T) or not (F)  ln_diaregmean_bin  = ', ln_diaregmean_bin
          WRITE(numout,*) 'Switch for regmean netcdf output (T) or not (F)  ln_diaregmean_nc  = ', ln_diaregmean_nc
          WRITE(numout,*) 'Switch for regmean kara mld terms (T) or not (F)  ln_diaregmean_karamld  = ', ln_diaregmean_karamld
          WRITE(numout,*) 'Switch for regmean PEA terms (T) or not (F)  ln_diaregmean_pea  = ', ln_diaregmean_pea
          WRITE(numout,*) 'Switch for regmean AR5 SLR terms (T) or not (F)  ln_diaregmean_diaar5  = ', ln_diaregmean_diaar5
          WRITE(numout,*) 'Switch for regmean Surface forcing terms (T) or not (F)  ln_diaregmean_diasbc  = ', ln_diaregmean_diasbc
          WRITE(numout,*) 'Switch for regmean BioGeoChemistry terms (T) or not (F)  ln_diaregmean_bgc  = ', ln_diaregmean_bgc
      ENDIF
      
      
      ALLOCATE( tmp_field_mat(jpi,jpj,19),  STAT= ierr ) !SS/NB/DT/ZA/VA T/S, SSH, MLD, PEA, PEAT, PEAS
          IF( ierr /= 0 )   CALL ctl_stop( 'tmp_field_mat: failed to allocate tmp_field_mat array' )
      tmp_field_mat(:,:,:) = 0.
      tmp_field_cnt = 0
      
      IF(ln_diaregmean_diaar5) THEN   
        ALLOCATE( tmp_field_AR5_mat(jpi,jpj,4),  STAT= ierr ) !SLR terms
            IF( ierr /= 0 )   CALL ctl_stop( 'tmp_field_AR5_mat: failed to allocate tmp_field_AR5_mat array' )
        tmp_field_AR5_mat(:,:,:) = 0.
      ENDIF
      
      IF(ln_diaregmean_diasbc) THEN   
        ALLOCATE( tmp_field_SBC_mat(jpi,jpj,7),  STAT= ierr ) !SBC terms
            IF( ierr /= 0 )   CALL ctl_stop( 'tmp_field_SBC_mat: failed to allocate tmp_field_SBC_mat array' )
        tmp_field_SBC_mat(:,:,:) = 0.
      ENDIF


#if defined key_fabm
      ! as there are so many BGC variables, write out the necessary iodef.xml and field_def.xml entries into ocean.output
      
      IF(ln_diaregmean_bgc) THEN   
            IF(lwp) THEN                   ! Control print

                BGC_nlevs = 5
                ALLOCATE( BGC_stat_name(6),BGC_lev_name(BGC_nlevs))
                nBGC_output = 16
                ALLOCATE( BGC_output_var(nBGC_output))

                BGC_output_var(1) = 'N1_p'
                BGC_output_var(2) = 'N3_n'
                BGC_output_var(3) = 'N4_n'
                BGC_output_var(4) = 'N5_s'
                BGC_output_var(5) = 'O2_o'
                BGC_output_var(6) = 'P1_Chl'
                BGC_output_var(7) = 'P2_Chl'
                BGC_output_var(8) = 'P3_Chl'
                BGC_output_var(9) = 'P4_Chl'
                BGC_output_var(10) = 'P1_c'
                BGC_output_var(11) = 'P2_c'
                BGC_output_var(12) = 'P3_c'
                BGC_output_var(13) = 'P4_c'
                BGC_output_var(14) = 'Z4_c'
                BGC_output_var(15) = 'Z5_c'
                BGC_output_var(16) = 'Z6_c'
                
                BGC_stat_name(1) = '_ave'
                BGC_stat_name(2) = '_tot'
                BGC_stat_name(3) = '_var'
                BGC_stat_name(4) = '_cnt'
                BGC_stat_name(5) = '_reg_id'
                BGC_stat_name(6) = '_mask_id'
                BGC_lev_name(1) = 'top'
                BGC_lev_name(2) = 'bot'
                BGC_lev_name(3) = 'dif'
                BGC_lev_name(4) = 'zav'
                BGC_lev_name(5) = 'vol'


                WRITE(numout,*) ''
                WRITE(numout,*) 'diaregmean BGC field_def.xml entries'
                WRITE(numout,*) ''
                
                
                DO jn=1,jp_fabm ! State loop
                    DO js=1,6
                        DO jl=1,BGC_nlevs
               
                            tmp_name=TRIM( TRIM("reg_")//TRIM(BGC_lev_name(jl))//TRIM("_")//TRIM(ctrcnm(jn))// TRIM(BGC_stat_name(js)) )

                            tmp_long_name  = TRIM(ctrcln(jn))
                            tmp_unit       = TRIM(ctrcun(jn))
                            
                            ! Where using volume integrated values, change units... 

                            IF ((jl .EQ. 5) .AND. (js .EQ. 2)) then
                                SELECT CASE (trim(tmp_unit))
                                    CASE ('mg C/m^3') ;     tmp_unit = 'Mg C (T C)' !'mg C/m^3'
                                    CASE ('mg/m^3') ;       tmp_unit = 'Mg (T)'     !'mg/m^3'
                                    CASE ('mmol C/m^3') ;   tmp_unit = 'Mmol C'     !'mmol C/m^3'
                                    CASE ('mmol N/m^3') ;   tmp_unit = 'Mmol N'     !'mmol N/m^3'
                                    CASE ('mmol O_2/m^3') ; tmp_unit = 'Mmol O'     !'mmol O_2/m^3'
                                    CASE ('mmol P/m^3') ;   tmp_unit = 'Mmol P'     !'mmol P/m^3'
                                    CASE ('mmol Si/m^3') ;  tmp_unit = 'Mmol S'     !'mmol Si/m^3'
                                    CASE ('umol/kg') ;      tmp_unit = 'Mmol'       !'umol/kg'   =  mmol/m^3
                                   ! CASE ('1/m') ;      cycle
                                    CASE DEFAULT
                                           tmp_unit = TRIM(TRIM(tmp_unit)//TRIM('x 1e9 m^3'))
                                END SELECT
                            ENDIF

                            WRITE(numout,*) TRIM(TRIM('<field id="')//TRIM(tmp_name)//TRIM('"         long_name="')// &
         &                        TRIM(BGC_lev_name(jl))//TRIM('_')//TRIM(tmp_long_name)//TRIM(BGC_stat_name(js))// &
         &                        TRIM('"     unit="'//TRIM(tmp_unit)  //'" />'))

                        END DO
                    END DO
                END DO
        
                WRITE(numout,*) ''
                WRITE(numout,*) 'diaregmean BGC iodef.xml entries'
                WRITE(numout,*) ''
                DO js=1,6

                    DO jn=1,jp_fabm ! State loop
                        
                        DO bgci=1,nBGC_output!
                            if (trim(ctrcnm(jn)) == TRIM(BGC_output_var(bgci))) CYCLE
                        ENDDO
                        DO jl=1,BGC_nlevs
                            ! only print out area averages for ss, nb, diff, and depth averaged, and total values for volume integrated
                            IF ((jl .EQ. 5) .AND. (js .NE. 2)) CYCLE ! cycle if vol, and not tot.
                            IF ((jl .NE. 5) .AND. (js .NE. 1)) CYCLE ! cycle if other levels, and not ave. 
               
                            tmp_name=TRIM(TRIM("reg_")//TRIM(BGC_lev_name(jl))//TRIM("_")//TRIM(ctrcnm(jn))// TRIM(BGC_stat_name(js))) 
                            tmp_long_name  = TRIM(ctrcln(jn))

                            WRITE(numout,*) TRIM(TRIM('<field field_ref="')//TRIM(tmp_name)//TRIM('"/>'))

                        END DO !level
                    END DO ! State loop
                END DO !statistic
                WRITE(numout,*) ''
                DEALLOCATE( BGC_stat_name,BGC_lev_name)

            ENDIF   ! Control print

      ENDIF !ln_diaregmean_bgc
      
#endif

      
      IF (ln_diaregmean) THEN
      
          ! Open region mask for region means, and retrieve the size of the mask (number of levels)          
          CALL iom_open ( 'region_mask.nc', inum )
          idmaskvar = iom_varid( inum, 'mask', kdimsz=zdimsz, kndims=zndims, ldstop = .FALSE.)          
          nmasks = zdimsz(3)
          
          ! read in the region mask (which contains floating point numbers) into a temporary array of reals.
          ALLOCATE( tmp_region_mask_real(jpi,jpj,nmasks),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_regmean_init: failed to allocate tmp_region_mask_real array' )
          
          ! Use jpdom_unknown to read in a n-layer mask.
          tmp_region_mask_real(:,:,:) = 0
          CALL iom_get( inum, jpdom_unknown, 'mask', tmp_region_mask_real(1:nlci,1:nlcj,1:nmasks),   &
              &          kstart = (/ mig(1),mjg(1),1 /), kcount = (/ nlci,nlcj,nmasks /) )
          
          CALL iom_close( inum )
          
          !Convert the region mask of reals into one of integers. 
          
          ALLOCATE( region_mask(jpi,jpj,nmasks),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_regmean_init: failed to allocate region_mask array' )
          region_mask(:,:,:) = 0
          region_mask = int(tmp_region_mask_real(:,:,:))
          DEALLOCATE( tmp_region_mask_real)
          
          
          ALLOCATE( nreg_mat(nmasks),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_regmean_init: failed to allocate nreg_mat array' )

          ! work out the number of regions in each mask, asssuming land is 0, and the regions are consectively numbered, 
          ! without missing any number, so the number of regions is the maximum number + 1 (for land). mpp_max across the 
          ! processors to get the global maxima
          DO maskno = 1,nmasks
              tmpint = maxval(region_mask(:,:,maskno))
              CALL mpp_max( tmpint )
              nreg_mat(maskno) = tmpint + 1
          END DO
        
          IF(lwp) THEN 
              ! if writing out as binary and text, open the files. 
              IF ( ln_diaregmean_bin ) THEN
                  ! Open binary for region means
                  CALL ctl_opn( numdct_reg_bin  ,'region_mean_timeseries.dat'  , 'NEW', 'UNFORMATTED', 'SEQUENTIAL', -1, numout,  .TRUE. )
              ENDIF
              
              IF ( ln_diaregmean_ascii ) THEN
                  ! Open text files for region means
                  CALL ctl_opn( numdct_reg_txt  ,'region_mean_timeseries.txt'  , 'NEW', 'FORMATTED', 'SEQUENTIAL', -1, numout,  .TRUE. )
              ENDIF
          ENDIF
     ENDIF

   END SUBROUTINE dia_regmean_init

   SUBROUTINE dia_calctmb_region_mean( pinfield,pouttmb )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_calctmb_region_mean  ***
      !!                   
      !! ** Purpose :    Find the Top, Bottom and Top minus Bottom fields of water Column
      !!            :    and depth average, and volume and mass intergated values. 

      !!
      !! ** Method  :   
      !!      use mbathy to find surface, mid and bottom of model levels
      !!
      !! History :
      !!   3.6  !  08-14  (E. O'Dea) Routine based on dia_wri_foam
      !!----------------------------------------------------------------------
      !! * Modules used

      ! Routine to map 3d field to top, middle, bottom
      IMPLICIT NONE


      ! Routine arguments
      REAL(wp), DIMENSION(jpi, jpj, jpk), INTENT(IN   ) :: pinfield    ! Input 3d field and mask
      REAL(wp), DIMENSION(jpi, jpj, 6  ), INTENT(  OUT) :: pouttmb     ! Output top, bottom and surface minus bed, zav, vol int, mass int

      ! Local variables
      INTEGER :: ji,jj,jk  ! Dummy loop indices

      ! Local Real
      REAL(wp)                         ::   zmdi  !  set masked values
      ! for depth int
      REAL(wp)                         ::   tmpnumer,tmpnumer_mass,tmpdenom ,z_av_val,vol_int_val

      zmdi=1.e+20 !missing data indicator for masking

      !zmdi=0 !missing data indicator for masking

      ! Calculate top
      pouttmb(:,:,1) = pinfield(:,:,1)*tmask(:,:,1)  + zmdi*(1.0-tmask(:,:,1))

     ! Calculate middle
      !DO jj = 1,jpj
      !    DO ji = 1,jpi
      !        jk              = max(1,mbathy(ji,jj)/2)
      !        pouttmb(ji,jj,2) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))
      !    END DO
      !END DO

      ! Calculate bottom, and top minus bottom
      DO jj = 1,jpj
          DO ji = 1,jpi
              IF ( tmask(ji,jj,1) .EQ. 1) THEN ! if land
               
                  jk              = max(1,mbathy(ji,jj) - 1)
                  pouttmb(ji,jj,2) = pinfield(ji,jj,jk)*tmask(ji,jj,jk)  + zmdi*(1.0-tmask(ji,jj,jk))

                  pouttmb(ji,jj,3) = (pouttmb(ji,jj,1) - pouttmb(ji,jj,2))*tmask(ji,jj,1)  + zmdi*(1.0-tmask(ji,jj,1))

                  !Depth and volume integral:
                  !---------------------------
                  !Vol int = Concentration * vol of grid box, summed over depth.
                  !Mass int = Concentration * vol of grid box * density of water, summed over depth.
                  !Depth Average = Vol int divided by * (vol of grid box summed over depth).

                  tmpnumer = 0.
                  tmpnumer_mass = 0.
                  tmpdenom = 0.
                  DO jk = 1,jpk
                     tmpnumer = tmpnumer + pinfield(ji,jj,jk)*tmask(ji,jj,jk)*e1t(ji,jj)*e2t(ji,jj)*e3t_n(ji,jj,jk)
                     tmpnumer_mass = tmpnumer_mass + pinfield(ji,jj,jk)*tmask(ji,jj,jk)*e1t(ji,jj)*e2t(ji,jj)*e3t_n(ji,jj,jk)*rhop(ji,jj,jk)
                     tmpdenom = tmpdenom +                    tmask(ji,jj,jk)*e1t(ji,jj)*e2t(ji,jj)*e3t_n(ji,jj,jk)
                  END DO
                  !z_av_val = tmpnumer/tmpdenom
                  !vol_int_val = tmpnumer
                  !mass_int_val = tmpnumer*density

                  pouttmb(ji,jj,4) = tmpnumer/tmpdenom ! depth averaged
                  pouttmb(ji,jj,5) = tmpnumer          ! Vol integrated 
                  pouttmb(ji,jj,6) = tmpnumer_mass     ! Mass integrated (for heat and salt calcs)
              ELSE
                    pouttmb(ji,jj,1) = zmdi
                    pouttmb(ji,jj,2) = zmdi
                    pouttmb(ji,jj,3) = zmdi
                    pouttmb(ji,jj,4) = zmdi
                    pouttmb(ji,jj,5) = zmdi
                    pouttmb(ji,jj,6) = zmdi
              ENDIF
          END DO
      END DO

   END SUBROUTINE dia_calctmb_region_mean


   SUBROUTINE dia_regmean( kt ) 
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE dia_regmean  ***
      !! ** Purpose :   Produce regional mean diagnostics
      !!
      !! ** Method  :   calls dia_wri_region_mean to calculate and write the regional means for a number of variables, 
      !!                (calling dia_calctmb_region_mean where necessary).
      !!                
      !!                Closes all text and binary files on last time step
      !!                
      !!      
      !!      
      !!
      !! History :
      !!   3.6  !  11-16  (J. Tinker) 
      !!         
      !!--------------------------------------------------------------------
      REAL(wp), POINTER, DIMENSION(:,:,:) :: tmp1mat   ! temporary array of 1's
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmbT    ! temporary T workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmbS    ! temporary S workspace 
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmb1    ! temporary density workspace 
      REAL(wp)                            ::   zmdi       ! set masked values
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      
      REAL(wp)                         ::   zdt  ! temporary reals
      INTEGER                          ::   i_steps, ierr         ! no of timesteps per hour, allocation error index
      INTEGER                          ::   maskno,jj,ji,jk,jm,nreg ! indices of mask, i and j, and number of regions

#if defined key_fabm
      INTEGER                          ::  jn ,tmp_dummy     ! set masked values
      REAL(wp)                         ::  tmp_val           ! tmp value, to allow min and max value clamping (not implemented)
      INTEGER                          ::  jl 
      CHARACTER (len=60) ::    tmp_name_bgc_top,tmp_name_bgc_bot,tmp_name_bgc_dif, tmp_name_bgc_zav, tmp_name_bgc_vol
      CHARACTER (len=60) ::    tmp_output_filename
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zwtmbBGC    ! temporary BGC workspace 

      LOGICAL       ::   verbose 
      verbose = .FALSE.
      tmp_val = 0
#endif
      zmdi=1.e+20 !missing data indicator for maskin

      IF (ln_diaregmean) THEN
        ! If regional mean calculations required by namelist
        ! -----------------
        ! identify hourly time steps (not used)
        zdt = rdt
        IF( nacc == 1 ) zdt = rdtmin

        IF( MOD( 3600,INT(zdt) ) == 0 ) THEN
            i_steps = 3600/INT(zdt)
        ELSE
            CALL ctl_stop('STOP', 'dia_regmean: timestep must give MOD(3600,rdt) = 0 otherwise no hourly values are possible')
        ENDIF
        
        ! Every time step, add physical, SBC, PEA, MLD terms to create hourly sums.
        ! Every hour, then hourly sums are divided by the number of timesteps in the hour to make hourly means
        ! These hourly mean values are then used to caluclate the regional means, and output with IOM. 
#if defined key_fabm
        ! BGC values are not averaged up over the hour, but are output as hourly instantaneous values. 
#endif

        
        !Extract 2d fields from 3d T and S with dia_calctmb_region_mean
        CALL wrk_alloc( jpi , jpj, 6 , zwtmbT )
        CALL wrk_alloc( jpi , jpj, 6 , zwtmbS )
        CALL wrk_alloc( jpi , jpj, 6 , zwtmb1 )
            
        CALL dia_calctmb_region_mean(  tsn(:,:,:,jp_tem),zwtmbT)
        CALL dia_calctmb_region_mean(  tsn(:,:,:,jp_sal),zwtmbS)

        ! To calc regional mean time series of int vol and mass, run region mean code on array of 1's...
        !   - then when multplying by volume, gives volume,
        !   - then when multplying by volume*density, gives mass

        CALL wrk_alloc( jpi , jpj, jpk , tmp1mat )
        DO jj = 1,jpj
            DO ji = 1,jpi
                DO jk = 1,jpk
                    tmp1mat(ji,jj,jk) = 1
                END DO
            END DO
        END DO

        CALL dia_calctmb_region_mean(  tmp1mat,zwtmb1)
        CALL wrk_dealloc( jpi , jpj, jpk , tmp1mat )
 
        ! Add 2d fields every time step to the hourly total.
            
        tmp_field_mat(:,:,1) = tmp_field_mat(:,:,1) + (zwtmbT(:,:,1)*tmask(:,:,1)) !sst
        tmp_field_mat(:,:,2) = tmp_field_mat(:,:,2) + (zwtmbT(:,:,2)*tmask(:,:,1)) !nbt
        tmp_field_mat(:,:,3) = tmp_field_mat(:,:,3) + (zwtmbT(:,:,3)*tmask(:,:,1)) !dft

        tmp_field_mat(:,:,4) = tmp_field_mat(:,:,4) + (zwtmbT(:,:,4)*tmask(:,:,1)) !zat
        tmp_field_mat(:,:,5) = tmp_field_mat(:,:,5) + (zwtmbT(:,:,5)*tmask(:,:,1)) !vat
        tmp_field_mat(:,:,6) = tmp_field_mat(:,:,6) + ((zwtmbT(:,:,6)*tmask(:,:,1)*4.2e3))! heat

        tmp_field_mat(:,:,7) = tmp_field_mat(:,:,7) + (zwtmbS(:,:,1)*tmask(:,:,1)) !sss
        tmp_field_mat(:,:,8) = tmp_field_mat(:,:,8) + (zwtmbS(:,:,2)*tmask(:,:,1)) !nbs
        tmp_field_mat(:,:,9) = tmp_field_mat(:,:,9) + (zwtmbS(:,:,3)*tmask(:,:,1)) !dfs

        tmp_field_mat(:,:,10) = tmp_field_mat(:,:,10) + (zwtmbS(:,:,4)*tmask(:,:,1)) !zas
        tmp_field_mat(:,:,11) = tmp_field_mat(:,:,11) + (zwtmbS(:,:,5)*tmask(:,:,1)) !vas
        tmp_field_mat(:,:,12) = tmp_field_mat(:,:,12) + (zwtmbS(:,:,6)*tmask(:,:,1)) !salt

        tmp_field_mat(:,:,13) = tmp_field_mat(:,:,13) + (zwtmb1(:,:,5)*tmask(:,:,1))!vol
        tmp_field_mat(:,:,14) = tmp_field_mat(:,:,14) + (zwtmb1(:,:,6)*tmask(:,:,1))!mass

        tmp_field_mat(:,:,15) = tmp_field_mat(:,:,15) + (sshn(:,:)*tmask(:,:,1)) !ssh
        
        CALL wrk_dealloc( jpi , jpj, 6 , zwtmbT )
        CALL wrk_dealloc( jpi , jpj, 6 , zwtmbS )
        CALL wrk_dealloc( jpi , jpj, 6 , zwtmb1 )
        
        IF( ln_diaregmean_karamld  ) THEN
            tmp_field_mat(:,:,16) = tmp_field_mat(:,:,16) + (hmld_kara(:,:)*tmask(:,:,1)) !mldkara
        ENDIF

        IF( ln_diaregmean_pea  ) THEN
            tmp_field_mat(:,:,17) = tmp_field_mat(:,:,17) + (pea(:,:)*tmask(:,:,1))  !pea
            tmp_field_mat(:,:,18) = tmp_field_mat(:,:,18) + (peat(:,:)*tmask(:,:,1)) !peat
            tmp_field_mat(:,:,19) = tmp_field_mat(:,:,19) + (peas(:,:)*tmask(:,:,1)) !peas
        ENDIF
          
        IF( ln_diaregmean_diaar5  ) THEN
            tmp_field_AR5_mat(:,:,1) = tmp_field_AR5_mat(:,:,1) + (sshsteric_mat(:,:)*tmask(:,:,1))
            tmp_field_AR5_mat(:,:,2) = tmp_field_AR5_mat(:,:,2) + (sshthster_mat(:,:)*tmask(:,:,1))
            tmp_field_AR5_mat(:,:,3) = tmp_field_AR5_mat(:,:,3) + (sshhlster_mat(:,:)*tmask(:,:,1))
            tmp_field_AR5_mat(:,:,4) = tmp_field_AR5_mat(:,:,4) + (zbotpres_mat(:,:)*tmask(:,:,1))
        ENDIF
        
        IF( ln_diaregmean_diasbc  ) THEN
            tmp_field_SBC_mat(:,:,1) = tmp_field_SBC_mat(:,:,1) + ((qsr  + qns)*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,2) = tmp_field_SBC_mat(:,:,2) + (qsr*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,3) = tmp_field_SBC_mat(:,:,3) + (qns*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,4) = tmp_field_SBC_mat(:,:,4) + (emp*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,5) = tmp_field_SBC_mat(:,:,5) + (wndm*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,6) = tmp_field_SBC_mat(:,:,6) + (pressnow*tmask(:,:,1))
            tmp_field_SBC_mat(:,:,7) = tmp_field_SBC_mat(:,:,7) + (rnf*tmask(:,:,1))
        ENDIF


        
        tmp_field_cnt = tmp_field_cnt + 1

        ! On the hour, calculate hourly means from the hourly total,and process the regional means. 
        
        IF ( MOD( kt, i_steps ) == 0 .and. kt .ne. nn_it000 ) THEN


            CALL dia_wri_region_mean(kt, "sst" , tmp_field_mat(:,:,1)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "nbt" , tmp_field_mat(:,:,2)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "dft" , tmp_field_mat(:,:,3)/real(tmp_field_cnt,wp))

            CALL dia_wri_region_mean(kt, "zat" , tmp_field_mat(:,:,4)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "vat" , tmp_field_mat(:,:,5)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "heat" , tmp_field_mat(:,:,6)/real(tmp_field_cnt,wp)/1e12)

            CALL dia_wri_region_mean(kt, "sss" , tmp_field_mat(:,:,7)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "nbs" , tmp_field_mat(:,:,8)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "dfs" , tmp_field_mat(:,:,9)/real(tmp_field_cnt,wp))

            CALL dia_wri_region_mean(kt, "zas" , tmp_field_mat(:,:,10)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "vas" , tmp_field_mat(:,:,11)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "salt" , tmp_field_mat(:,:,12)/real(tmp_field_cnt,wp)/1e12)

            CALL dia_wri_region_mean(kt, "vol" , tmp_field_mat(:,:,13)/real(tmp_field_cnt,wp))
            CALL dia_wri_region_mean(kt, "mass" , tmp_field_mat(:,:,14)/real(tmp_field_cnt,wp))

            CALL dia_wri_region_mean(kt, "ssh" , tmp_field_mat(:,:,15)/real(tmp_field_cnt,wp))


            IF( ln_diaregmean_karamld  ) THEN          
                CALL dia_wri_region_mean(kt, "mldkara" , tmp_field_mat(:,:,16)/real(tmp_field_cnt,wp)) ! tm
            ENDIF

            IF( ln_diaregmean_pea  ) THEN          
                CALL dia_wri_region_mean(kt, "pea"  , tmp_field_mat(:,:,17)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "peat" , tmp_field_mat(:,:,18)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "peas" , tmp_field_mat(:,:,19)/real(tmp_field_cnt,wp)) ! tmb
            ENDIF

            tmp_field_mat(:,:,:) = 0.

            IF( ln_diaregmean_diaar5  ) THEN

                CALL dia_wri_region_mean(kt, "ssh_steric" ,      tmp_field_AR5_mat(:,:,1)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "ssh_thermosteric", tmp_field_AR5_mat(:,:,2)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "ssh_halosteric" ,  tmp_field_AR5_mat(:,:,3)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "bot_pres" ,        tmp_field_AR5_mat(:,:,4)/real(tmp_field_cnt,wp))
                tmp_field_AR5_mat(:,:,:) = 0.
            ENDIF

            IF( ln_diaregmean_diasbc  ) THEN

                CALL dia_wri_region_mean(kt, "qt"   , tmp_field_SBC_mat(:,:,1)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "qsr"  , tmp_field_SBC_mat(:,:,2)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "qns"  , tmp_field_SBC_mat(:,:,3)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "emp"  , tmp_field_SBC_mat(:,:,4)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "wspd" , tmp_field_SBC_mat(:,:,5)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "mslp" , tmp_field_SBC_mat(:,:,6)/real(tmp_field_cnt,wp))
                CALL dia_wri_region_mean(kt, "rnf"  , tmp_field_SBC_mat(:,:,7)/real(tmp_field_cnt,wp))
                tmp_field_SBC_mat(:,:,:) = 0.
            ENDIF

#if defined key_fabm
            !ADD Biogeochemistry
            
            IF( ln_diaregmean_bgc  ) THEN !ln_diaregmean_bgc
                
                ! Loop through 3d BGC tracers
                DO jn=1,jp_fabm ! State loop

                    ! get variable name for different levels
                    tmp_name_bgc_top=TRIM(TRIM("top_")//TRIM(ctrcnm(jn)))
                    tmp_name_bgc_bot=TRIM(TRIM("bot_")//TRIM(ctrcnm(jn)))
                    tmp_name_bgc_dif=TRIM(TRIM("dif_")//TRIM(ctrcnm(jn)))
                    tmp_name_bgc_zav=TRIM(TRIM("zav_")//TRIM(ctrcnm(jn)))
                    tmp_name_bgc_vol=TRIM(TRIM("vol_")//TRIM(ctrcnm(jn)))

                    ! print out names if verbose
                    IF(verbose .AND. lwp) THEN 
                        WRITE(numout,*)
                        WRITE(numout,*) 'dia_regmean tmp_name_bgc_top : ',TRIM(tmp_name_bgc_top)
                        WRITE(numout,*) 'dia_regmean tmp_name_bgc_bot : ',TRIM(tmp_name_bgc_bot)
                        WRITE(numout,*) 'dia_regmean tmp_name_bgc_dif : ',TRIM(tmp_name_bgc_dif)
                        WRITE(numout,*) 'dia_regmean tmp_name_bgc_zav : ',TRIM(tmp_name_bgc_zav)
                        WRITE(numout,*) 'dia_regmean tmp_name_bgc_vol : ',TRIM(tmp_name_bgc_vol)
                        CALL FLUSH(numout) 

                    ENDIF

                    !Allocate working array, and get surface, bed etc fields. 
                    CALL wrk_alloc( jpi , jpj,  6 , zwtmbBGC )
                    CALL dia_calctmb_region_mean(  trn(:,:,:,jn),zwtmbBGC )


                    !Print out 2d fields to ascii text files to check values if verbose. (24MB per time step, per BGC variable)
                    IF (verbose) THEN 

                        WRITE (tmp_output_filename, "(A4,I3.3,A1,I6.6,A1,I3.3,A4)") "bgc_",jn,"_",kt,"_",narea,".txt"
                        WRITE (*,*) tmp_output_filename
                        OPEN(UNIT=74,FILE=TRIM(tmp_output_filename))

                        DO ji = 1,jpi
                            DO jj = 1,jpj
                                WRITE(74,FMT="(I4,I4,F3,F25.5,F25.5,F25.5,F25.5,F25.5)") nimpp+ji, njmpp+jj,tmask(ji,jj,1),&
                                      & zwtmbBGC(ji,jj,1),zwtmbBGC(ji,jj,2),zwtmbBGC(ji,jj,3),zwtmbBGC(ji,jj,4),zwtmbBGC(ji,jj,5)/1e9
                            END DO
                        END DO
                        CLOSE(74)
                    ENDIF
                    
                    ! Do region means
                    CALL dia_wri_region_mean(kt, TRIM(tmp_name_bgc_top)   , zwtmbBGC(:,:,1))
                    CALL dia_wri_region_mean(kt, TRIM(tmp_name_bgc_bot)   , zwtmbBGC(:,:,2))
                    CALL dia_wri_region_mean(kt, TRIM(tmp_name_bgc_dif)   , zwtmbBGC(:,:,3))
                    CALL dia_wri_region_mean(kt, TRIM(tmp_name_bgc_zav)   , zwtmbBGC(:,:,4))
                    CALL dia_wri_region_mean(kt, TRIM(tmp_name_bgc_vol)   , zwtmbBGC(:,:,5)/1e9)


                    !Deallocate working array
                    CALL wrk_dealloc( jpi , jpj,  6 , zwtmbBGC )
                ENDDO ! State loop
            ENDIF !ln_diaregmean_bgc

#endif
              
            tmp_field_cnt = 0
  
        ENDIF ! ( MOD( kt, i_steps ) == 0  .and. kt .ne. nn_it000 )
        
        
        ! If on the last time step, close binary and ascii files. 
        IF( kt == nitend ) THEN
            IF(lwp) THEN
                IF ( ln_diaregmean_bin ) THEN
                    !Closing binary files for regional mean time series.
                    CLOSE(numdct_reg_bin)
                ENDIF
                IF ( ln_diaregmean_ascii ) THEN
                    !Closing text files for regional mean time series.
                    CLOSE(numdct_reg_txt)
                ENDIF

                DEALLOCATE( region_mask, nreg_mat, tmp_field_mat)
                IF( ln_diaregmean_diaar5  ) DEALLOCATE( tmp_field_AR5_mat)
                IF( ln_diaregmean_diasbc  ) DEALLOCATE( tmp_field_SBC_mat)
            ENDIF
        ENDIF
          
          
      ELSE
        CALL ctl_warn('dia_regmean: regmean diagnostic is set to false you should not have seen this')
      ENDIF
      
   END SUBROUTINE dia_regmean
   
   
   SUBROUTINE dia_wri_region_mean(kt, tmp_name,         infield  )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dia_tmb  ***
      !!                   
      !! ** Purpose :   Calculate and write region mean time series for 2d arrays
      !!
      !! ** Method  :   
      !!      use 
      !!
      !! History :
      !!   ??  !  15/10/2015  (JTinker) Routine taken from old dia_wri_foam
      !!----------------------------------------------------------------------
      !! * Modules used
      !use lib_mpp
      !use lib_fortr
      IMPLICIT NONE
      
      INTEGER, INTENT(in) ::   kt
      CHARACTER (len=*) , INTENT(IN   ) ::    tmp_name
      REAL(wp), DIMENSION(jpi, jpj), INTENT(IN   ) :: infield    ! Input 3d field and mask
      
      ! Local variables
      INTEGER, DIMENSION(jpi, jpj) :: internal_region_mask    ! Input 3d field and mask 
      REAL(wp), DIMENSION(jpi, jpj) :: internal_infield    ! Internal data field
      REAL(wp), ALLOCATABLE, DIMENSION(:) ::   zrmet_ave,zrmet_tot,zrmet_var,zrmet_cnt,zrmet_mask_id,zrmet_reg_id  ,zrmet_min,zrmet_max
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   zrmet_out
      REAL(wp), ALLOCATABLE,   DIMENSION(:) ::   ave_mat,tot_mat,num_mat,var_mat,ssq_mat,cnt_mat,reg_id_mat,mask_id_mat ,min_mat,max_mat   !: region_mask
      
      REAL(wp)                         ::   zmdi, zrmet_val      ! set masked values
      INTEGER :: maskno,nreg  ! ocean time-step indexocean time step            
      INTEGER :: ji,jj,jk,ind,jm ! Dummy loop indices
      INTEGER :: reg_ind_cnt ! Dummy loop indices
      
      INTEGER  ::   ierr      
      REAL(wp)  :: tmpreal
      CHARACTER(LEN=180) :: FormatString,nreg_string,tmp_name_iom
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   dummy_zrmet
      LOGICAL       ::   verbose     
      verbose = .False.


      zmdi=1.e+20 !missing data indicator for maskin
      
      !Allocate output arrays for iomput, set to zmdi, and set a region counter = 1
      ALLOCATE( zrmet_ave(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_ave array' )
      ALLOCATE( zrmet_tot(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_tot array' )
      ALLOCATE( zrmet_var(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_var array' )
      ALLOCATE( zrmet_cnt(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_cnt array' )
      ALLOCATE( zrmet_mask_id(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_mask_id array' )
      ALLOCATE( zrmet_reg_id(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_reg_id array' )


      ALLOCATE( zrmet_min(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_min array' )
      ALLOCATE( zrmet_max(n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_max array' )
      
      ALLOCATE( zrmet_out(jpi,jpj,n_regions_output),  STAT= ierr )
        IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate zrmet_reg_id array' )

  
      
        IF(lwp .AND. verbose) THEN 
              WRITE(numout,*)
              WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//';'
              WRITE(numout,*)
        ENDIF
        
      DO ji = 1,jpi
          DO jj = 1,jpj
            internal_infield(ji,jj) = infield(ji,jj)
          END DO
      END DO
        
      ! Check for NANS # JT   03/09/2018
      DO ji = 1,jpi
          DO jj = 1,jpj
              IF ( tmask(ji,jj,1) == 1.0_wp ) THEN
                  IF ( internal_infield(ji,jj) .ne. internal_infield(ji,jj) ) THEN
                      WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//' Nan at (kt,i,j): ',kt,ji - (-jpizoom+1-nimpp+1),jj - (-jpjzoom+1-njmpp+1)
                      internal_infield(ji,jj) = 0.
                  ENDIF
              ELSE                
                  IF ( internal_infield(ji,jj) .ne. internal_infield(ji,jj) ) THEN
                      WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//' Masked Nan at (kt,i,j): ',kt,ji - (-jpizoom+1-nimpp+1),jj - (-jpjzoom+1-njmpp+1)
                      internal_infield(ji,jj) = 0.
                  ENDIF
              ENDIF
          END DO
      END DO
      
      
      zrmet_ave(:) = zmdi
      zrmet_tot(:) = zmdi
      zrmet_var(:) = zmdi
      zrmet_cnt(:) = zmdi
      zrmet_mask_id(:) = zmdi
      zrmet_reg_id(:) = zmdi
      
      zrmet_min(:) = zmdi
      zrmet_max(:) = zmdi
      reg_ind_cnt = 1
      
      
      ! loop though the masks
      DO maskno = 1,nmasks
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; begin mask loops: ',maskno
          
          
          ! For each mask, get the number of regions (nreg), and a local copy of the region. 
          nreg = nreg_mat(maskno)
          internal_region_mask = region_mask(:,:,maskno)
          
          ! allocate temporary stat arrays, and set to zero
          ALLOCATE( ave_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate ave_mat array' )
          ALLOCATE( tot_mat(nreg),      STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate tot_mat array' )
          ALLOCATE( num_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate num_mat array' )
          ALLOCATE( var_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate var_mat array' )
          ALLOCATE( ssq_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate ssq_mat array' )
          ALLOCATE( cnt_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate cnt_mat array' )

          ALLOCATE( min_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate min_mat array' )
          ALLOCATE( max_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate max_mat array' )

          ALLOCATE( reg_id_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate reg_id_mat array' )
          ALLOCATE( mask_id_mat(nreg),  STAT= ierr )
          IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate mask_id_mat array' )


          
          ave_mat(:) = 0.
          tot_mat(:) = 0.
          num_mat(:) = 0.
          var_mat(:) = 0.
          cnt_mat(:) = 0.
          ssq_mat(:) = 0.

          min_mat(:) = zmdi
          max_mat(:) = -zmdi
          reg_id_mat(:) = 0.
          mask_id_mat(:) = 0.
          
          ! loop though the array. for each sea grid box where tmask == 1), 
          ! read which region the grid box is in, add the value of the gridbox (and its square) 
          ! to the total for that region, and then increment the counter for that region.
          !CALL cpu_time(start_reg_mean_loop)
          !WRITE(numout,*) kt,start_reg_mean_loop
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; begin spatial loops: '
          DO ji = nldi,nlei
              DO jj = nldj,nlej
                    IF ( tmask(ji,jj,1) == 1.0_wp ) THEN
                        ind = internal_region_mask(ji,jj)+1
                        tot_mat(ind) = tot_mat(ind) + (internal_infield(ji,jj))
                        ssq_mat(ind) = ssq_mat(ind) + ( internal_infield(ji,jj) *  internal_infield(ji,jj))
                        cnt_mat(ind) = cnt_mat(ind) + 1.

                        min_mat(ind) = min(min_mat(ind),internal_infield(ji,jj))
                        max_mat(ind) = max(max_mat(ind),internal_infield(ji,jj))
                    ENDIF
              END DO
          END DO
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; finish spatial loops: '
          ! sum the totals, the counts, and the squares across the processors          
          CALL mpp_sum( tot_mat,nreg )
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; finished mpp_sum 1'
          CALL mpp_sum( ssq_mat,nreg )
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; finished mpp_sum 2'
          CALL mpp_sum( cnt_mat,nreg )
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; finished mpp_sum 2'

          CALL mpp_min( min_mat,nreg )
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; finished mpp_min'
          CALL mpp_max( max_mat,nreg )
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_wri_region_mean : '//tmp_name//'; finished mpp_max'
          
          
          !calculate the mean and variance from the total, sum of squares and the count. 
          
          ave_mat = tot_mat(:)/cnt_mat(:)
          var_mat = ssq_mat(:)/cnt_mat(:) - (ave_mat(:)*ave_mat(:))
          
          
          !mask array of mask and region number. 
          DO jj = 1,nreg
              reg_id_mat(jj) = real(jj-1)
              mask_id_mat(jj) = real(maskno)
          END DO
          
          
          !write text and binary, and note region statistics for current mask for later iom_put
          IF( lwp ) THEN 
          
              !Write out ascii and binary if requred
              IF ( ln_diaregmean_bin ) THEN
                  !Writing out regional mean time series to binary files
                  WRITE(numdct_reg_bin) tmp_name,kt,maskno,n_regions_output
                  WRITE(numdct_reg_bin) ave_mat
                  WRITE(numdct_reg_bin) tot_mat
                  WRITE(numdct_reg_bin) var_mat
                  WRITE(numdct_reg_bin) ssq_mat
                  WRITE(numdct_reg_bin) cnt_mat
                  WRITE(numdct_reg_bin) min_mat
                  WRITE(numdct_reg_bin) max_mat
              ENDIF
              
              IF ( ln_diaregmean_ascii  ) THEN
                  !Writing out regional mean time series to text files

                  WRITE(nreg_string, "(I5)") nreg
                  FormatString = "(A30,"//trim(nreg_string)//"F25.3)"
                  WRITE(numdct_reg_txt, FMT="(A30,I6,I6)") tmp_name,kt,maskno            
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"ave_mat:", ave_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"tot_mat:", tot_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"var_mat:", var_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"ssq_mat:", ssq_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"cnt_mat:", cnt_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"min_mat:", min_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"max_mat:", max_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"reg_mat:", reg_id_mat
                  WRITE(numdct_reg_txt, FMT=trim(FormatString)) trim(tmp_name)//" "//"msk_mat:", mask_id_mat

              ENDIF
              
              DO jm = 1,nreg
                  zrmet_ave(    reg_ind_cnt) =     ave_mat(jm)
                  zrmet_tot(    reg_ind_cnt) =     tot_mat(jm)
                  zrmet_var(    reg_ind_cnt) =     var_mat(jm)
                  zrmet_cnt(    reg_ind_cnt) =     cnt_mat(jm)
                  zrmet_min(    reg_ind_cnt) =     min_mat(jm)
                  zrmet_max(    reg_ind_cnt) =     max_mat(jm)
                  zrmet_reg_id( reg_ind_cnt) =  reg_id_mat(jm)
                  zrmet_mask_id(reg_ind_cnt) = mask_id_mat(jm)
                
                  reg_ind_cnt = reg_ind_cnt + 1 
              END DO
          
          ENDIF
        
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean about to deallocated arrays for ',kt,maskno
          DEALLOCATE(ave_mat,tot_mat,num_mat,var_mat,ssq_mat,cnt_mat,min_mat,max_mat,reg_id_mat,mask_id_mat)

          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean deallocated arrays for ',kt,maskno
          IF(lwp)CALL FLUSH(numdct_reg_txt)
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean flushed region mean text for ',kt,maskno
      END DO

      IF(lwp .AND. verbose) THEN                   ! Control print
         WRITE(numout,*) 'dia_regmean ready to start iom_put'
         CALL FLUSH(numout)
      ENDIF
      
      !With current field_def.xml and iodef.xml, these fields must be output, so set to dummy values if not required.
      
      IF ( ln_diaregmean_nc ) THEN
      
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''

          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean ready to start iom_put: ',trim(tmp_name)
          
      
          DO jm = 1,n_regions_output
            zrmet_val = zrmet_ave(jm)
!            if (zrmet_val .LT. -1e16) zrmet_val = -1e16
!            if (zrmet_val .GT. 1e16) zrmet_val = 1e16
            if (zrmet_val .NE. zrmet_val) zrmet_val = 1e20
            zrmet_out(:,:,jm) = zrmet_val
          END DO      
          tmp_name_iom =  trim(trim("reg_") // trim(tmp_name) // trim('_ave'))
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean iom_put tmp_name_iom : ',trim(tmp_name_iom)
          CALL iom_put(trim(tmp_name_iom), zrmet_out )
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''
      
          DO jm = 1,n_regions_output
            zrmet_val = zrmet_tot(jm)
!            if (zrmet_val .LT. -1e16) zrmet_val = -1e16
!            if (zrmet_val .GT. 1e16) zrmet_val = 1e16
            if (zrmet_val .NE. zrmet_val) zrmet_val = 1e20
            zrmet_out(:,:,jm) = zrmet_val
          END DO
          tmp_name_iom =  trim(trim("reg_") // trim(tmp_name) // trim('_tot'))
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean iom_put tmp_name_iom : ',trim(tmp_name_iom)          
          CALL iom_put( trim(tmp_name_iom), zrmet_out )
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''
      
          DO jm = 1,n_regions_output
            zrmet_val = zrmet_var(jm)
!            if (zrmet_val .LT. -1e16) zrmet_val = -1e16
!            if (zrmet_val .GT. 1e16) zrmet_val = 1e16
            if (zrmet_val .NE. zrmet_val) zrmet_val = 1e20
            zrmet_out(:,:,jm) = zrmet_val
          END DO
          tmp_name_iom =  trim(trim("reg_") // trim(tmp_name) // trim('_var'))
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean iom_put tmp_name_iom : ',trim(tmp_name_iom)          
          CALL iom_put( trim(tmp_name_iom), zrmet_out )
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''
      
          DO jm = 1,n_regions_output
            zrmet_val = zrmet_cnt(jm)
!            if (zrmet_val .LT. -1e16) zrmet_val = -1e16
!            if (zrmet_val .GT. 1e16) zrmet_val = 1e16
            if (zrmet_val .NE. zrmet_val) zrmet_val = 1e20
            zrmet_out(:,:,jm) = zrmet_val
          END DO
          tmp_name_iom =  trim(trim("reg_") // trim(tmp_name) // trim('_cnt'))
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean iom_put tmp_name_iom : ',trim(tmp_name_iom)          
          CALL iom_put( trim(tmp_name_iom), zrmet_out )
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''
      
          DO jm = 1,n_regions_output
            zrmet_val = zrmet_reg_id(jm)
!            if (zrmet_val .LT. -1e16) zrmet_val = -1e16
!            if (zrmet_val .GT. 1e16) zrmet_val = 1e16
            if (zrmet_val .NE. zrmet_val) zrmet_val = 1e20
            zrmet_out(:,:,jm) = zrmet_val
          END DO
          tmp_name_iom =  trim(trim("reg_") // trim(tmp_name) // trim('_reg_id'))
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean iom_put tmp_name_iom : ',trim(tmp_name_iom)          
          CALL iom_put( trim(tmp_name_iom), zrmet_out )
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''
      
          DO jm = 1,n_regions_output
            zrmet_val = zrmet_mask_id(jm)
!            if (zrmet_val .LT. -1e16) zrmet_val = -1e16
!            if (zrmet_val .GT. 1e16) zrmet_val = 1e16
            if (zrmet_val .NE. zrmet_val) zrmet_val = 1e20
            zrmet_out(:,:,jm) = zrmet_val
          END DO
          tmp_name_iom =  trim(trim("reg_") // trim(tmp_name) // trim('_mask_id'))
          IF(lwp .AND. verbose) WRITE(numout,*) 'dia_regmean iom_put tmp_name_iom : ',trim(tmp_name_iom)          
          CALL iom_put( trim(tmp_name_iom), zrmet_out )
          zrmet_out(:,:,:) = 0
          zrmet_val = 0
          tmp_name_iom = ''
      ELSE
        
          ALLOCATE( dummy_zrmet(jpi,jpj,n_regions_output),  STAT= ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'dia_wri_region_mean: failed to allocate dummy_zrmet array' )

          DO jm = 1,n_regions_output
              dummy_zrmet(:,:,jm) =     real(jm,wp)
          END DO

          DO jm = 1,9
              CALL iom_put( trim(trim("reg_") // trim(tmp_name) // trim('_ave')), dummy_zrmet )
              CALL iom_put( trim(trim("reg_") // trim(tmp_name) // trim('_tot')), dummy_zrmet )
              CALL iom_put( trim(trim("reg_") // trim(tmp_name) // trim('_var')), dummy_zrmet )
              CALL iom_put( trim(trim("reg_") // trim(tmp_name) // trim('_cnt')), dummy_zrmet )
              CALL iom_put( trim(trim("reg_") // trim(tmp_name) // trim('_reg_id')), dummy_zrmet )
              CALL iom_put( trim(trim("reg_") // trim(tmp_name) // trim('_mask_id')), dummy_zrmet )
          END DO
    
          DEALLOCATE( dummy_zrmet)
      ENDIF
      
      DEALLOCATE(zrmet_ave,zrmet_tot,zrmet_var,zrmet_cnt,zrmet_mask_id,zrmet_reg_id,zrmet_min,zrmet_max,zrmet_out)

      IF(lwp .AND. verbose) THEN                   ! Control print
         WRITE(numout,*) 
         WRITE(numout,*) 'dia_wri_region_mean finished for ', trim(tmp_name)
         WRITE(numout,*) 
         CALL FLUSH(numout)
      ENDIF
      
   END SUBROUTINE dia_wri_region_mean


   !!======================================================================
END MODULE diaregmean
