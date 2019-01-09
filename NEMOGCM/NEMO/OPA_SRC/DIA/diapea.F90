MODULE diapea
   !!======================================================================
   !!                       ***  MODULE  diapea  ***
   !! Potential Energy Anomaly 
   !!======================================================================
   !! History :  3.6  !  12/2016  (J Tinker)  Original code
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain
   USE in_out_manager  ! I/O units
   USE iom             ! I/0 library
   USE wrk_nemo        ! working arrays
   USE eosbn2          ! Equation of state - in situ and potential density
   USE phycst          ! physical constant

   IMPLICIT NONE
   PRIVATE 
   
   PUBLIC   dia_pea_init            ! routine called by nemogcm.F90
   PUBLIC   dia_pea                 ! routine called by diawri.F90
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE,   DIMENSION(:,:)  ::   pea,peat,peas   
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:,:)  ::   wgt_co_mat   ! Weighting array for proportion of grid shallower than cut off depth
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:) ::   t_zmean, s_zmean  !Depth mean temperature and salinity: 2d fields
   REAL(wp), SAVE, ALLOCATABLE,   DIMENSION(:,:,:) ::   t_zmean_mat, s_zmean_mat  !Depth mean temperature and salinity: 3d fields
   REAL(wp) ::   zcutoff
   LOGICAL , PUBLIC ::   ln_pea  ! region mean calculation
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.6 , NEMO Consortium (2014)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_pea_init
      ! Local variables
      INTEGER :: ji,jj,jk  ! Dummy loop indices
      REAL(wp) :: sumz,tmpsumz
      INTEGER  :: ierr                ! error integer for IOM_get
      INTEGER ::   ios                  ! Local integer output status for namelist read
      
      zcutoff = 200.!200m
      
      NAMELIST/nam_pea/ ln_pea
      
      
      !
      REWIND ( numnam_ref )              ! Read Namelist nam_diatmb in referdiatmbence namelist : TMB diagnostics
      READ   ( numnam_ref, nam_pea, IOSTAT=ios, ERR= 901 )
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_pea in reference namelist', lwp )

      REWIND( numnam_cfg )              ! Namelist nam_diatmb in configuration namelist  TMB diagnostics
      READ  ( numnam_cfg, nam_pea, IOSTAT = ios, ERR = 902 )
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'nam_pea in configuration namelist', lwp )
      IF(lwm) WRITE ( numond, nam_pea )

      IF(lwp) THEN                   ! Control print
          WRITE(numout,*)
          WRITE(numout,*) 'dia_pea_init : Output potential energy anomaly Diagnostics'
          WRITE(numout,*) '~~~~~~~~~~~~'
          WRITE(numout,*) 'Namelist nam_pea : set pea output '
          WRITE(numout,*) 'Switch for pea diagnostics (T) or not (F)  ln_diaregmean  = ', ln_pea
      ENDIF
      
      
   ALLOCATE( pea(jpi,jpj),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate pea array' )
   
   ALLOCATE( peat(jpi,jpj),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate peat array' )
   
   ALLOCATE( peas(jpi,jpj),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate peas array' )   
      
   ALLOCATE( t_zmean_mat(jpi,jpj,jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate t_zmean_mat array' )
   
   ALLOCATE( s_zmean_mat(jpi,jpj,jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate s_zmean_mat array' )
   
   ALLOCATE( t_zmean(jpi,jpj),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate t_zmean array' )
   
   ALLOCATE( s_zmean(jpi,jpj),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate s_zmean array' )
   
   ALLOCATE( wgt_co_mat(jpi,jpj,jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate wgt_co_mat array' )
   
   
   pea(:,:) = 0.
   peat(:,:) = 0.
   peas(:,:) = 0.
    
   if ( ln_pea ) THEN
          
      ! create wgt_co_mat mat, with the proportion of the grid (gdept_0) below cut off (200m)
      
      DO jj = 1,jpj
          DO ji = 1,jpi
              IF ( tmask(ji,jj,1) == 1.0_wp ) THEN
                sumz = 0.
                DO jk = 1,jpk
                  IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                    tmpsumz = sumz + e3t_n(ji,jj,jk)
                    IF (sumz .ge. zcutoff) THEN
                      ! Already too deep
                      wgt_co_mat(ji,jj,jk) = 0.
                    ELSE IF (tmpsumz .le. zcutoff) THEN
                      ! Too shallow
                      wgt_co_mat(ji,jj,jk) = 1.
                    ELSE
                      !proprotion of grid box above cut off depth
                      wgt_co_mat(ji,jj,jk) = (zcutoff-Sumz)/e3t_n(ji,jj,jk)
                    END IF
                    sumz = tmpsumz
                  endif
                END DO
                
              ELSE
                !if land, set to 0.
                DO jk = 1,jpk
                  wgt_co_mat(ji,jj,jk) = 0.
                END DO
                  
              ENDIF
          END DO
       END DO
   ENDIF
   

   END SUBROUTINE dia_pea_init
   

   SUBROUTINE dia_pea(kt)
   INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
    
   INTEGER :: ji,jj,jk,ierr  ! Dummy loop indices
   REAL(wp) :: tmpdenom, tmpnum, maxz
   !rau0
   
   REAL(wp), ALLOCATABLE,   DIMENSION(:,:,:,:) ::   ts_pea_mat,ts_pea_mat_TS_mean,ts_pea_mat_S_mean
   REAL(wp), ALLOCATABLE,   DIMENSION(:,:,:) ::   tmp_pea_rho,tmp_pea_TS_mean_rho,tmp_pea_S_mean_rho
   REAL(wp) ::   int_y_pea,int_y_pea_t
   
   
      
   ALLOCATE( ts_pea_mat(jpi,jpj,jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat array' )    
   ALLOCATE( ts_pea_mat_TS_mean(jpi,jpj,jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat_TS_mean array' )   
   ALLOCATE( ts_pea_mat_S_mean(jpi,jpj,jpk,2),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate ts_pea_mat_S_mean array' )
   ALLOCATE( tmp_pea_rho(jpi,jpj,jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_rho array' )
   ALLOCATE( tmp_pea_TS_mean_rho(jpi,jpj,jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_TS_mean_rho array' )
   ALLOCATE( tmp_pea_S_mean_rho(jpi,jpj,jpk),  STAT= ierr )
   IF( ierr /= 0 )   CALL ctl_stop( 'dia_pea_init: failed to allocate tmp_pea_S_mean_rho array' )

    pea(:,:)=0.0d0
    peaT(:,:)=0.0d0
    peaS(:,:)=0.0d0
    
    
    
    !calculate the depth mean temperature and salinity of the upper 200m. Save this into a 3d array. Set the value where tmask=0 to be tsn. 
    
     DO jj = 1,jpj
         DO ji = 1,jpi
            IF ( tmask(ji,jj,1) == 1.0_wp ) THEN  ! if a sea point. 
            
             !Depth mean temperature
             tmpdenom = 0.
             tmpnum = 0.
             DO jk = 1,jpk
                IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                  tmpnum = tmpnum + (wgt_co_mat(ji,jj,jk)*e3t_n(ji,jj,jk)*tsn(ji,jj,jk,jp_tem))
                  tmpdenom = tmpdenom + (wgt_co_mat(ji,jj,jk)*e3t_n(ji,jj,jk))
                endif
             END DO
             t_zmean(ji,jj) = tmpnum/tmpdenom
             
             !Depth mean salinity
             tmpdenom = 0.
             tmpnum = 0.
             DO jk = 1,jpk
                IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                  tmpnum     = tmpnum + (wgt_co_mat(ji,jj,jk)*e3t_n(ji,jj,jk)*tsn(ji,jj,jk,jp_sal))
                  tmpdenom = tmpdenom + (wgt_co_mat(ji,jj,jk)*e3t_n(ji,jj,jk))
                endif
             END DO
             s_zmean(ji,jj) = tmpnum/tmpdenom
             
             !save into a 3d grid
             DO jk = 1,jpk
               IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN
                 t_zmean_mat(ji,jj,jk) = t_zmean(ji,jj)
                 s_zmean_mat(ji,jj,jk) = s_zmean(ji,jj)
               else
                 t_zmean_mat(ji,jj,jk) = tsn(ji,jj,jk,jp_tem)
                 s_zmean_mat(ji,jj,jk) = tsn(ji,jj,jk,jp_sal)
               endif
             END DO
            else
              t_zmean(ji,jj) = tsn(ji,jj,1,jp_tem)
              s_zmean(ji,jj) = tsn(ji,jj,1,jp_sal)
              DO jk = 1,jpk
                t_zmean_mat(ji,jj,jk) = tsn(ji,jj,jk,jp_tem)
                s_zmean_mat(ji,jj,jk) = tsn(ji,jj,jk,jp_sal)
              END DO
            endif
            
         END DO
     END DO
      
   !Calculate the density from the depth varying, and depth average temperature and salinity
   !-----------------------------
   !-----------------------------
   
   ts_pea_mat(:,:,:,:) = tsn(:,:,:,:)
   
   ts_pea_mat_TS_mean(:,:,:,1) = t_zmean_mat(:,:,:)
   ts_pea_mat_TS_mean(:,:,:,2) = s_zmean_mat(:,:,:)
   
   ts_pea_mat_S_mean(:,:,:,1) = t_zmean_mat(:,:,:)
   ts_pea_mat_S_mean(:,:,:,2) = tsn(:,:,:,jp_sal)
   
   CALL eos ( ts_pea_mat,         tmp_pea_rho,         gdept_n(:,:,:) )
   CALL eos ( ts_pea_mat_TS_mean, tmp_pea_TS_mean_rho, gdept_n(:,:,:) )
   CALL eos ( ts_pea_mat_S_mean,  tmp_pea_S_mean_rho,  gdept_n(:,:,:) )
   tmp_pea_rho = (tmp_pea_rho * rau0) + rau0
   tmp_pea_TS_mean_rho = (tmp_pea_TS_mean_rho * rau0) + rau0
   tmp_pea_S_mean_rho = (tmp_pea_S_mean_rho * rau0) + rau0
   
   
   ! to test the density calculation
   !CALL iom_put( "tmp_pea_rho" , tmp_pea_rho )                 ! pea
   !CALL iom_put( "tmp_pea_TS_mean_rho" , tmp_pea_TS_mean_rho )                 ! pea
   
   
   ! Caluclation of the PEA.
    DO jj = 1,jpj
      DO ji = 1,jpi
        pea(ji,jj) = 0.
        peat(ji,jj) = 0.
        peas(ji,jj) = 0.
        maxz = 0.
        int_y_pea = 0.
        int_y_pea_t = 0.
        IF ( tmask(ji,jj,1) == 1.0_wp ) THEN ! for sea points
        
          ! the depth integrated calculation is summed up over the depths, and then divided by the depth
          DO jk = 1,jpk
            !for each level...
            
            IF ( tmask(ji,jj,jk) == 1.0_wp ) THEN ! if above the sea bed...
              int_y_pea =   -((tmp_pea_TS_mean_rho(ji,jj,jk)) - (tmp_pea_rho(ji,jj,jk)))*9.81*gdept_n(ji,jj,jk)*wgt_co_mat(ji,jj,jk)*e3t_n(ji,jj,jk)
              int_y_pea_t = -((tmp_pea_S_mean_rho(ji,jj,jk))  - (tmp_pea_rho(ji,jj,jk)))*9.81*gdept_n(ji,jj,jk)*wgt_co_mat(ji,jj,jk)*e3t_n(ji,jj,jk)
            else
              int_y_pea = 0.
              int_y_pea_t = 0.
            endif
            
            ! check that the sum is not NaN. 
            if ( int_y_pea .ne.  int_y_pea    ) int_y_pea = 0.
            if ( int_y_pea_t .ne. int_y_pea_t )  int_y_pea_t = 0.
            !if ( (int_y_pea*int_y_pea    ) .gt. 1.0e6 ) int_y_pea = 0.
            !if ( (int_y_pea_t*int_y_pea_t) .gt. 1.0e6 ) int_y_pea_t = 0.
            
            pea(ji,jj) =  pea(ji,jj) + int_y_pea
            peat(ji,jj) =  peat(ji,jj) + int_y_pea_t
            maxz = maxz + (e3t_n(ji,jj,jk)*wgt_co_mat(ji,jj,jk))
          enddo
            
            
          !divide by the depth
          pea(ji,jj) = pea(ji,jj)/maxz
          peat(ji,jj) = peat(ji,jj)/maxz
          peas(ji,jj) = pea(ji,jj) - peat(ji,jj)
            
            
          else
            pea(ji,jj) = 0.
            peat(ji,jj) = 0.
            peas(ji,jj) = 0.
          endif
       enddo
    enddo
!    
     CALL iom_put( "pea" , pea )                 ! pea
     CALL iom_put( "peat" , peat )               ! pea
     CALL iom_put( "peas" , peas )               ! pea
           

    DEALLOCATE(ts_pea_mat,ts_pea_mat_TS_mean,ts_pea_mat_S_mean,tmp_pea_rho,tmp_pea_TS_mean_rho,tmp_pea_S_mean_rho)
   
   END SUBROUTINE dia_pea

END MODULE diapea