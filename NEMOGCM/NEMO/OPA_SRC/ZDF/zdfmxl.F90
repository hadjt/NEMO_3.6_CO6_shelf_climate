MODULE zdfmxl
   !!======================================================================
   !!                       ***  MODULE  zdfmxl  ***
   !! Ocean physics: mixed layer depth 
   !!======================================================================
   !! History :  1.0  ! 2003-08  (G. Madec)  original code
   !!            3.2  ! 2009-07  (S. Masson, G. Madec)  IOM + merge of DO-loop
   !!            3.7  ! 2012-03  (G. Madec)  make public the density criteria for trdmxl 
   !!             -   ! 2014-02  (F. Roquet)  mixed layer depth calculated using N2 instead of rhop 
   !!----------------------------------------------------------------------
   !!   zdf_mxl      : Compute the turbocline and mixed layer depths.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE phycst          ! physical constants
   USE iom             ! I/O library
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! work arrays
   USE timing          ! Timing
   USE trc_oce, ONLY : lk_offline ! offline flag
   USE lbclnk 
   USE eosbn2          ! Equation of state

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_mxl       ! called by step.F90
   PUBLIC   zdf_mxl_tref  ! called by asminc.F90
   PUBLIC   zdf_mxl_alloc ! Used in zdf_tke_init

   INTEGER , PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   nmln    !: number of level in the mixed layer (used by TOP)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmld    !: mixing layer depth (turbocline)      [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlp    !: mixed layer depth  (rho=rho0+zdcrit) [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmlpt   !: mixed layer depth at t-points        [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmld_tref  !: mixed layer depth at t-points - temperature criterion [m]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::   hmld_kara  !: mixed layer depth of Kara et al   [m]

   REAL(wp), PUBLIC ::   rho_c = 0.01_wp    !: density criterion for mixed layer depth
   REAL(wp)         ::   avt_c = 5.e-4_wp   ! Kz criterion for the turbocline depth

   ! Namelist variables for  namzdf_karaml 
 
   LOGICAL   :: ln_kara            ! Logical switch for calculating kara mixed
                                     ! layer
   LOGICAL   :: ln_kara_write25h   ! Logical switch for 25 hour outputs
   INTEGER   :: jpmld_type         ! mixed layer type             
   REAL(wp)  :: ppz_ref            ! depth of initial T_ref 
   REAL(wp)  :: ppdT_crit          ! Critical temp diff  
   REAL(wp)  :: ppiso_frac         ! Fraction of ppdT_crit used 
   
   !Used for 25h mean
   LOGICAL, PRIVATE :: kara_25h_init = .TRUE.   !Logical used to initalise 25h 
                                                !outputs. Necissary, because we need to
                                                !initalise the kara_25h on the zeroth
                                                !timestep (i.e in the nemogcm_init call)
   REAL, PRIVATE, ALLOCATABLE, DIMENSION(:,:) :: hmld_kara_25h

   !! * Substitutions
#  include "domzgr_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_mxl_alloc()
      !!----------------------------------------------------------------------
      !!               ***  FUNCTION zdf_mxl_alloc  ***
      !!----------------------------------------------------------------------
      zdf_mxl_alloc = 0      ! set to zero if no array to be allocated
      IF( .NOT. ALLOCATED( nmln ) ) THEN
         ALLOCATE( nmln(jpi,jpj), hmld(jpi,jpj), hmlp(jpi,jpj), hmlpt(jpi,jpj), &
         &                           hmld_tref(jpi,jpj), STAT= zdf_mxl_alloc )
         !
         IF( lk_mpp             )   CALL mpp_sum ( zdf_mxl_alloc )
         IF( zdf_mxl_alloc /= 0 )   CALL ctl_warn('zdf_mxl_alloc: failed to allocate arrays.')
         !
      ENDIF
   END FUNCTION zdf_mxl_alloc


   SUBROUTINE zdf_mxl( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdfmxl  ***
      !!                   
      !! ** Purpose :   Compute the turbocline depth and the mixed layer depth
      !!              with density criteria.
      !!
      !! ** Method  :   The mixed layer depth is the shallowest W depth with 
      !!      the density of the corresponding T point (just bellow) bellow a
      !!      given value defined locally as rho(10m) + rho_c
      !!               The turbocline depth is the depth at which the vertical
      !!      eddy diffusivity coefficient (resulting from the vertical physics
      !!      alone, not the isopycnal part, see trazdf.F) fall below a given
      !!      value defined locally (avt_c here taken equal to 5 cm/s2 by default)
      !!
      !! ** Action  :   nmln, hmld, hmlp, hmlpt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      INTEGER  ::   iikn, iiki, ikt, imkt   ! local integer
      REAL(wp) ::   zN2_c        ! local scalar
      INTEGER, POINTER, DIMENSION(:,:) ::   imld   ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_mxl')
      !
      CALL wrk_alloc( jpi,jpj, imld )

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'zdf_mxl : mixed layer depth'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
         !                             ! allocate zdfmxl arrays
         IF( zdf_mxl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl : unable to allocate arrays' )
      ENDIF

      ! w-level of the mixing and mixed layers
      nmln(:,:)  = nlb10               ! Initialization to the number of w ocean point
      hmlp(:,:)  = 0._wp               ! here hmlp used as a dummy variable, integrating vertically N^2
      zN2_c = grav * rho_c * r1_rau0   ! convert density criteria into N^2 criteria
      DO jk = nlb10, jpkm1
         DO jj = 1, jpj                ! Mixed layer level: w-level 
            DO ji = 1, jpi
               ikt = mbkt(ji,jj)
               hmlp(ji,jj) = hmlp(ji,jj) + MAX( rn2b(ji,jj,jk) , 0._wp ) * fse3w(ji,jj,jk)
               IF( hmlp(ji,jj) < zN2_c )   nmln(ji,jj) = MIN( jk , ikt ) + 1   ! Mixed layer level
            END DO
         END DO
      END DO
      !
      ! w-level of the turbocline
      imld(:,:) = mbkt(:,:) + 1        ! Initialization to the number of w ocean point
      DO jk = jpkm1, nlb10, -1         ! from the bottom to nlb10 
         DO jj = 1, jpj
            DO ji = 1, jpi
               imkt = mikt(ji,jj)
               IF( avt (ji,jj,jk) < avt_c )   imld(ji,jj) = MAX( imkt, jk )      ! Turbocline 
            END DO
         END DO
      END DO
      ! depth of the mixing and mixed layers

      CALL zdf_mxl_kara( kt ) ! kara MLD

      DO jj = 1, jpj
         DO ji = 1, jpi
            iiki = imld(ji,jj)
            iikn = nmln(ji,jj)
            imkt = mikt(ji,jj)
            hmld (ji,jj) = ( fsdepw(ji,jj,iiki  ) - fsdepw(ji,jj,imkt ) ) * ssmask(ji,jj)    ! Turbocline depth 
            hmlp (ji,jj) = ( fsdepw(ji,jj,iikn  ) - fsdepw(ji,jj,imkt ) ) * ssmask(ji,jj)    ! Mixed layer depth
            hmlpt(ji,jj) = ( fsdept(ji,jj,iikn-1) - fsdepw(ji,jj,imkt ) ) * ssmask(ji,jj)    ! depth of the last T-point inside the mixed layer
         END DO
      END DO
      IF( .NOT.lk_offline ) THEN            ! no need to output in offline mode
         CALL iom_put( "mldr10_1", hmlp )   ! mixed layer depth
         CALL iom_put( "mldkz5"  , hmld )   ! turbocline depth
      ENDIF
      
      IF(ln_ctl)   CALL prt_ctl( tab2d_1=REAL(nmln,wp), clinfo1=' nmln : ', tab2d_2=hmlp, clinfo2=' hmlp : ', ovlap=1 )
      !
      CALL wrk_dealloc( jpi,jpj, imld )
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_mxl')
      !
   END SUBROUTINE zdf_mxl


   SUBROUTINE zdf_mxl_tref()
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_mxl_tref  ***
      !!                   
      !! ** Purpose :   Compute the mixed layer depth with temperature criteria.
      !!
      !! ** Method  :   The temperature-defined mixed layer depth is required
      !!                   when assimilating SST in a 2D analysis. 
      !!
      !! ** Action  :   hmld_tref
      !!----------------------------------------------------------------------
      !
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   t_ref               ! Reference temperature  
      REAL(wp) ::   temp_c = 0.2        ! temperature criterion for mixed layer depth  
      !!----------------------------------------------------------------------
      !
      ! Initialise array
      IF( zdf_mxl_alloc() /= 0 )   CALL ctl_stop( 'STOP', 'zdf_mxl_tref : unable to allocate arrays' )
      
      !For the AMM model assimiation uses a temperature based mixed layer depth  
      !This is defined here  
      DO jj = 1, jpj  
         DO ji = 1, jpi  
           hmld_tref(ji,jj)=fsdept(ji,jj,1  )   
           IF(ssmask(ji,jj) > 0.)THEN  
             t_ref=tsn(ji,jj,1,jp_tem) 
             DO jk=2,jpk  
               IF(ssmask(ji,jj)==0.)THEN  
                  hmld_tref(ji,jj)=fsdept(ji,jj,jk )  
                  EXIT  
               ELSEIF( ABS(tsn(ji,jj,jk,jp_tem)-t_ref) < temp_c)THEN  
                  hmld_tref(ji,jj)=fsdept(ji,jj,jk )  
               ELSE  
                  EXIT  
               ENDIF  
             ENDDO  
           ENDIF  
         ENDDO  
      ENDDO

   END SUBROUTINE zdf_mxl_tref

   SUBROUTINE zdf_mxl_kara( kt ) 
      !!---------------------------------------------------------------------------------- 
      !!                    ***  ROUTINE zdf_mxl_kara  *** 
      !                                                                        
      !   Calculate mixed layer depth according to the definition of          
      !   Kara et al, 2000, JGR, 105, 16803.  
      ! 
      !   If mld_type=1 the mixed layer depth is calculated as the depth at which the  
      !   density has increased by an amount equivalent to a temperature difference of  
      !   0.8C at the surface. 
      ! 
      !   For other values of mld_type the mixed layer is calculated as the depth at  
      !   which the temperature differs by 0.8C from the surface temperature.  
      !                                                                        
      !   Original version: David Acreman                                      
      ! 
      !!-----------------------------------------------------------------------------------
     
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index 
 
      NAMELIST/namzdf_karaml/ ln_kara,jpmld_type, ppz_ref, ppdT_crit, &
      &                       ppiso_frac, ln_kara_write25h 
 
      ! Local variables                                                        
      REAL, DIMENSION(jpi,jpk) :: ppzdep      ! depth for use in calculating d(rho) 
      REAL(wp), DIMENSION(jpi,jpj,jpts) :: ztsn_2d  !Local version of tsn 
      LOGICAL :: ll_found(jpi,jpj)              ! Is T_b to be found by interpolation ? 
      LOGICAL :: ll_belowml(jpi,jpj,jpk)        ! Flag points below mixed layer when ll_found=F 
      INTEGER :: ji, jj, jk                     ! loop counter 
      INTEGER :: ik_ref(jpi,jpj)                ! index of reference level 
      INTEGER :: ik_iso(jpi,jpj)                ! index of last uniform temp level 
      REAL    :: zT(jpi,jpj,jpk)                ! Temperature or denisty 
      REAL    :: zT_ref(jpi,jpj)                ! reference temperature 
      REAL    :: zT_b                           ! base temperature 
      REAL    :: zdTdz(jpi,jpj,jpk-2)           ! gradient of zT 
      REAL    :: zmoddT(jpi,jpj,jpk-2)          ! Absolute temperature difference 
      REAL    :: zdz                            ! depth difference 
      REAL    :: zdT                            ! temperature difference 
      REAL    :: zdelta_T(jpi,jpj)              ! difference critereon 
      REAL    :: zRHO1(jpi,jpj), zRHO2(jpi,jpj) ! Densities
      INTEGER, SAVE :: idt, i_steps
      INTEGER, SAVE :: i_cnt_25h     ! Counter for 25 hour means
      INTEGER :: ios                 ! Local integer output status for namelist read

     
 
      !!------------------------------------------------------------------------------------- 
 
      IF( kt == nit000 ) THEN 
         ! Default values 
         ln_kara      = .FALSE.
         ln_kara_write25h = .FALSE. 
         jpmld_type   = 1     
         ppz_ref      = 10.0 
         ppdT_crit    = 0.2  
         ppiso_frac   = 0.1   
         ! Read namelist
         REWIND ( numnam_ref )
         READ   ( numnam_ref, namzdf_karaml, IOSTAT=ios, ERR= 901 )
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_karaml in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist nam_diadiaop in configuration namelist  3D hourly diagnostics
         READ  ( numnam_cfg,  namzdf_karaml, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namzdf_karaml in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namzdf_karaml )
 

         WRITE(numout,*) '===== Kara mixed layer =====' 
         WRITE(numout,*) 'ln_kara = ',    ln_kara
         WRITE(numout,*) 'jpmld_type = ', jpmld_type 
         WRITE(numout,*) 'ppz_ref = ',    ppz_ref 
         WRITE(numout,*) 'ppdT_crit = ',  ppdT_crit 
         WRITE(numout,*) 'ppiso_frac = ', ppiso_frac
         WRITE(numout,*) 'ln_kara_write25h = ', ln_kara_write25h
         WRITE(numout,*) '============================' 
      
         IF ( .NOT.ln_kara ) THEN
            WRITE(numout,*) "ln_kara not set; Kara mixed layer not calculated" 
         ELSE
            IF (.NOT. ALLOCATED(hmld_kara) ) ALLOCATE( hmld_kara(jpi,jpj) )
            IF ( ln_kara_write25h .AND. kara_25h_init ) THEN
               i_cnt_25h = 0
               IF (.NOT. ALLOCATED(hmld_kara_25h) ) &
               &   ALLOCATE( hmld_kara_25h(jpi,jpj) )
               hmld_kara_25h = 0._wp
               IF( nacc == 1 ) THEN
                  idt = INT(rdtmin)
               ELSE
                  idt = INT(rdt)
               ENDIF
               IF( MOD( 3600,idt) == 0 ) THEN 
                  i_steps = 3600 / idt  
               ELSE 
                  CALL ctl_stop('STOP', &
                  & 'zdf_mxl_kara: timestep must give MOD(3600,rdt)'// &
                  & ' = 0 otherwise no hourly values are possible') 
               ENDIF  
            ENDIF
         ENDIF
      ENDIF
      
      IF ( ln_kara ) THEN
         
         !set critical ln_kara
         ztsn_2d = tsn(:,:,1,:)
         ztsn_2d(:,:,jp_tem) = ztsn_2d(:,:,jp_tem) + ppdT_crit
     
         ! Set the mixed layer depth criterion at each grid point 
         ppzdep = 0._wp
         IF( jpmld_type == 1 ) THEN                                         
            CALL eos ( tsn(:,:,1,:), &
            &          ppzdep(:,:), zRHO1(:,:) ) 
            CALL eos ( ztsn_2d(:,:,:), &
            &           ppzdep(:,:), zRHO2(:,:) ) 
            zdelta_T(:,:) = abs( zRHO1(:,:) - zRHO2(:,:) ) * rau0 
            ! RHO from eos (2d version) doesn't calculate north or east halo: 
            CALL lbc_lnk( zdelta_T, 'T', 1. ) 
            zT(:,:,:) = rhop(:,:,:) 
         ELSE 
            zdelta_T(:,:) = ppdT_crit                      
            zT(:,:,:) = tsn(:,:,:,jp_tem)                           
         ENDIF 
     
         ! Calculate the gradient of zT and absolute difference for use later 
         DO jk = 1 ,jpk-2 
            zdTdz(:,:,jk)  =    ( zT(:,:,jk+1) - zT(:,:,jk) ) / fse3w(:,:,jk+1) 
            zmoddT(:,:,jk) = abs( zT(:,:,jk+1) - zT(:,:,jk) ) 
         END DO 
     
         ! Find density/temperature at the reference level (Kara et al use 10m).          
         ! ik_ref is the index of the box centre immediately above or at the reference level 
         ! Find ppz_ref in the array of model level depths and find the ref    
         ! density/temperature by linear interpolation.                                   
         ik_ref = -1
         DO jk = jpkm1, 2, -1 
            WHERE( fsdept(:,:,jk) > ppz_ref ) 
               ik_ref(:,:) = jk - 1 
               zT_ref(:,:) = zT(:,:,jk-1) + &
               &             zdTdz(:,:,jk-1) * ( ppz_ref - fsdept(:,:,jk-1) ) 
            ENDWHERE 
         END DO
         IF ( ANY( ik_ref  < 0 ) .OR. ANY( ik_ref  > jpkm1 ) ) THEN
            CALL ctl_stop( "STOP", &
            & "zdf_mxl_kara: unable to find reference level for kara ML" ) 
         ELSE
            ! If the first grid box centre is below the reference level then use the 
            ! top model level to get zT_ref 
            WHERE( fsdept(:,:,1) > ppz_ref )  
               zT_ref = zT(:,:,1) 
               ik_ref = 1 
            ENDWHERE 
     
            ! Search for a uniform density/temperature region where adjacent levels          
            ! differ by less than ppiso_frac * deltaT.                                      
            ! ik_iso is the index of the last level in the uniform layer  
            ! ll_found indicates whether the mixed layer depth can be found by interpolation 
            ik_iso(:,:)   = ik_ref(:,:) 
            ll_found(:,:) = .false. 
            DO jj = 1, nlcj 
               DO ji = 1, nlci 
                 !CDIR NOVECTOR 
                  DO jk = ik_ref(ji,jj), mbathy(ji,jj)-1 
                     IF( zmoddT(ji,jj,jk) > ( ppiso_frac * zdelta_T(ji,jj) ) ) THEN 
                        ik_iso(ji,jj)   = jk 
                        ll_found(ji,jj) = ( zmoddT(ji,jj,jk) > zdelta_T(ji,jj) ) 
                        EXIT 
                     ENDIF 
                  END DO 
               END DO 
            END DO 
     
            ! Use linear interpolation to find depth of mixed layer base where possible 
            hmld_kara(:,:) = ppz_ref 
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  IF( ll_found(ji,jj) .and. tmask(ji,jj,1) == 1.0 ) THEN 
                     zdz =  abs( zdelta_T(ji,jj) / zdTdz(ji,jj,ik_iso(ji,jj)) ) 
                     hmld_kara(ji,jj) = fsdept(ji,jj,ik_iso(ji,jj)) + zdz 
                  ENDIF 
               END DO 
            END DO 
     
            ! If ll_found = .false. then calculate MLD using difference of zdelta_T    
            ! from the reference density/temperature 
     
            ! Prevent this section from working on land points 
            WHERE( tmask(:,:,1) /= 1.0 ) 
               ll_found = .true. 
            ENDWHERE 
     
            DO jk = 1, jpk 
               ll_belowml(:,:,jk) = abs( zT(:,:,jk) - zT_ref(:,:) ) >= &
               & zdelta_T(:,:)
            END DO 
     
            ! Set default value where interpolation cannot be used (ll_found=false)  
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  IF( .NOT. ll_found(ji,jj) )  &
                  &   hmld_kara(ji,jj) = fsdept(ji,jj,mbathy(ji,jj)) 
               END DO 
            END DO 
     
            DO jj = 1, jpj 
               DO ji = 1, jpi 
                  !CDIR NOVECTOR 
                  DO jk = ik_ref(ji,jj)+1, mbathy(ji,jj) 
                     IF( ll_found(ji,jj) ) EXIT 
                     IF( ll_belowml(ji,jj,jk) ) THEN                
                        zT_b = zT_ref(ji,jj) + zdelta_T(ji,jj) * &
                        &      SIGN(1.0, zdTdz(ji,jj,jk-1) ) 
                        zdT  = zT_b - zT(ji,jj,jk-1)                                      
                        zdz  = zdT / zdTdz(ji,jj,jk-1)                                       
                        hmld_kara(ji,jj) = fsdept(ji,jj,jk-1) + zdz 
                        EXIT                                                   
                     ENDIF 
                  END DO 
               END DO 
            END DO 
     
            hmld_kara(:,:) = hmld_kara(:,:) * tmask(:,:,1) 
 
            IF(  ln_kara_write25h  ) THEN
               !Double IF required as i_steps not defined if ln_kara_write25h =
               ! FALSE
               IF ( ( MOD( kt, i_steps ) == 0 ) .OR.  kara_25h_init ) THEN
                  hmld_kara_25h = hmld_kara_25h + hmld_kara
                  i_cnt_25h = i_cnt_25h + 1
                  IF ( kara_25h_init ) kara_25h_init = .FALSE.
               ENDIF
            ENDIF
 
#if defined key_iomput 
            IF( kt /= nit000 ) THEN 
               CALL iom_put( "mldkara"  , hmld_kara )   
               IF( ( MOD( i_cnt_25h, 25) == 0 ) .AND.  ln_kara_write25h ) &
                  CALL iom_put( "kara25h"  , ( hmld_kara_25h / 25._wp ) )
            ENDIF
#endif
 
         ENDIF
      ENDIF
       
   END SUBROUTINE zdf_mxl_kara 

   !!======================================================================
END MODULE zdfmxl
