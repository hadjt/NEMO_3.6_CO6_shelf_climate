MODULE sbcflx
   !!======================================================================
   !!                       ***  MODULE  sbcflx  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  1.0  !  2006-06  (G. Madec)  Original code
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namflx   : flux formulation namlist
   !!   sbc_flx  : flux formulation as ocean surface boundary condition (forced mode, fluxes read in NetCDF files)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE sbcdcy          ! surface boundary condition: diurnal cycle on qsr
   USE phycst          ! physical constants
   USE fldread         ! read input fields
   USE iom             ! IOM library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_flx       ! routine called by step.F90

   INTEGER , PARAMETER ::   jpfld   = 6   ! maximum number of files to read 
   INTEGER , PARAMETER ::   jp_utau = 1   ! index of wind stress (i-component) file
   INTEGER , PARAMETER ::   jp_vtau = 2   ! index of wind stress (j-component) file
   INTEGER , PARAMETER ::   jp_qtot = 3   ! index of total (non solar+solar) heat file
   INTEGER , PARAMETER ::   jp_qsr  = 4   ! index of solar heat file
   INTEGER , PARAMETER ::   jp_emp  = 5   ! index of evaporation-precipation file
   INTEGER , PARAMETER ::   jp_press = 6  ! index of pressure for UKMO shelf fluxes
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf    ! structure of input fields (file informations, fields read)
   LOGICAL , PUBLIC    ::   ln_shelf_flx = .FALSE. ! UKMO SHELF specific flux flag
   INTEGER             ::   jpfld_local   ! maximum number of files to read (locally modified depending on ln_shelf_flx) 

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO-consortium (2010) 
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_flx( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau  (N/m2)
      !!                   j-component of the stress              vtau  (N/m2)
      !!                   net downward heat flux                 qtot  (watt/m2)
      !!                   net downward radiative flux            qsr   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!                 - the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   update at each time-step
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns         non solar heat flux including heat flux due to emp
      !!              - qsr         solar heat flux
      !!              - emp         upward mass flux (evap. - precip.)
      !!              - sfx         salt flux; set to zero at nit000 but possibly non-zero
      !!                            if ice is present (computed in limsbc(_2).F90)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, jf            ! dummy indices
      INTEGER  ::   ierror                ! return error code
      INTEGER  ::   ios                   ! Local integer output status for namelist read
      REAL(wp) ::   zfact                 ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      REAL     ::   cs                    ! UKMO SHELF: Friction co-efficient at surface
      REAL     ::   totwindspd            ! UKMO SHELF: Magnitude of wind speed vector

      REAL(wp) ::   rhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   cdrag = 1.5e-3       ! drag coefficient 
      !!
      CHARACTER(len=100) ::  cn_dir                               ! Root directory for location of flx files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                    ! array of namelist information structures
      TYPE(FLD_N) ::   sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp, sn_press  !  informations about the fields to be read
      LOGICAL     ::   ln_foam_flx  = .FALSE.                     ! UKMO FOAM specific flux flag
      NAMELIST/namsbc_flx/ cn_dir, sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp,   &
      &                    ln_foam_flx, sn_press, ln_shelf_flx
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN                ! First call kt=nit000  
         ! set file information
         REWIND( numnam_ref )              ! Namelist namsbc_flx in reference namelist : Files for fluxes
         READ  ( numnam_ref, namsbc_flx, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_flx in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_flx in configuration namelist : Files for fluxes
         READ  ( numnam_cfg, namsbc_flx, IOSTAT = ios, ERR = 902 )
902      IF( ios /= 0 ) CALL ctl_nam ( ios , 'namsbc_flx in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_flx ) 
         !
         !                                         ! check: do we plan to use ln_dm2dc with non-daily forcing?
         IF( ln_dm2dc .AND. sn_qsr%nfreqh /= 24 )   &
            &   CALL ctl_stop( 'sbc_blk_core: ln_dm2dc can be activated only with daily short-wave forcing' ) 
         !
         !                                         ! store namelist information in an array
         slf_i(jp_utau) = sn_utau   ;   slf_i(jp_vtau) = sn_vtau
         slf_i(jp_qtot) = sn_qtot   ;   slf_i(jp_qsr ) = sn_qsr 
         slf_i(jp_emp ) = sn_emp
         !
            ALLOCATE( sf(jpfld), STAT=ierror )        ! set sf structure
            IF( ln_shelf_flx ) slf_i(jp_press) = sn_press
   
            ! define local jpfld depending on shelf_flx logical
            IF( ln_shelf_flx ) THEN
               jpfld_local = jpfld
            ELSE
               jpfld_local = jpfld-1
            ENDIF
            !
         IF( ierror > 0 ) THEN   
            CALL ctl_stop( 'sbc_flx: unable to allocate sf structure' )   ;   RETURN  
         ENDIF
         DO ji= 1, jpfld_local
            ALLOCATE( sf(ji)%fnow(jpi,jpj,1) )
            IF( slf_i(ji)%ln_tint ) ALLOCATE( sf(ji)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_flx', 'flux formulation for ocean surface boundary condition', 'namsbc_flx' )
         !
         sfx(:,:) = 0.0_wp                         ! salt flux due to freezing/melting (non-zero only if ice is present; set in limsbc(_2).F90)
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )                            ! input fields provided at the current time-step
     
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency

         IF( ln_dm2dc ) THEN   ;   qsr(:,:) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) )   ! modify now Qsr to include the diurnal cycle
         ELSE                  ;   qsr(:,:) =          sf(jp_qsr)%fnow(:,:,1)
         ENDIF
!CDIR COLLAPSE
            !!UKMO SHELF effect of atmospheric pressure on SSH
            ! If using ln_apr_dyn, this is done there so don't repeat here.
            IF( ln_shelf_flx .AND. .NOT. ln_apr_dyn) THEN
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     apgu(ji,jj) = (-1.0/rau0)*(sf(jp_press)%fnow(ji+1,jj,1)-sf(jp_press)%fnow(ji,jj,1))/e1u(ji,jj)
                     apgv(ji,jj) = (-1.0/rau0)*(sf(jp_press)%fnow(ji,jj+1,1)-sf(jp_press)%fnow(ji,jj,1))/e2v(ji,jj)
                  END DO
               END DO
            ENDIF ! ln_shelf_flx
      
         DO jj = 1, jpj                                           ! set the ocean fluxes from read fields
            DO ji = 1, jpi
                   IF( ln_shelf_flx ) THEN
                      !! UKMO SHELF - need atmospheric pressure to calculate Haney forcing
                      pressnow(ji,jj) = sf(jp_press)%fnow(ji,jj,1)
                      !! UKMO SHELF flux files contain wind speed not wind stress
                      totwindspd = sqrt((sf(jp_utau)%fnow(ji,jj,1))**2.0 + (sf(jp_vtau)%fnow(ji,jj,1))**2.0)
                      cs = 0.63 + (0.066 * totwindspd)
                      utau(ji,jj) = cs * (rhoa/rau0) * sf(jp_utau)%fnow(ji,jj,1) * totwindspd
                      vtau(ji,jj) = cs * (rhoa/rau0) * sf(jp_vtau)%fnow(ji,jj,1) * totwindspd
                   ELSE
                      utau(ji,jj) = sf(jp_utau)%fnow(ji,jj,1)
                      vtau(ji,jj) = sf(jp_vtau)%fnow(ji,jj,1)
                   ENDIF
                   qsr (ji,jj) = sf(jp_qsr )%fnow(ji,jj,1)
                   IF( ln_foam_flx .OR. ln_shelf_flx ) THEN
                      !! UKMO FOAM flux files contain non-solar heat flux (qns) rather than total heat flux (qtot)
                      qns (ji,jj) = sf(jp_qtot)%fnow(ji,jj,1)
                      !! UKMO FOAM flux files contain the net DOWNWARD freshwater flux P-E rather then E-P
                      emp (ji,jj) = -1. * sf(jp_emp )%fnow(ji,jj,1)
                   ELSE
                      qns (ji,jj) = sf(jp_qtot)%fnow(ji,jj,1) - sf(jp_qsr)%fnow(ji,jj,1)
                      emp (ji,jj) = sf(jp_emp )%fnow(ji,jj,1)
                   ENDIF
            END DO
         END DO
         !                                                        ! add to qns the heat due to e-p
         qns(:,:) = qns(:,:) - emp(:,:) * sst_m(:,:) * rcp        ! mass flux is at SST
         !
   
            !! UKMO FOAM wind fluxes need lbc_lnk calls owing to a bug in interp.exe
            IF( ln_foam_flx ) THEN
               CALL lbc_lnk( utau(:,:), 'U', -1. )
               CALL lbc_lnk( vtau(:,:), 'V', -1. )
            ENDIF
    
         !                                                        ! module of wind stress and wind speed at T-point
         zcoef = 1. / ( zrhoa * zcdrag )
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
!CDIR NOVERRCHK
            DO ji = fs_2, fs_jpim1   ! vect. opt.
               ztx = utau(ji-1,jj  ) + utau(ji,jj) 
               zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
               zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
               taum(ji,jj) = zmod
               wndm(ji,jj) = SQRT( zmod * zcoef )
            END DO
         END DO
         taum(:,:) = taum(:,:) * tmask(:,:,1) ; wndm(:,:) = wndm(:,:) * tmask(:,:,1)
         CALL lbc_lnk( taum(:,:), 'T', 1. )   ;   CALL lbc_lnk( wndm(:,:), 'T', 1. )

         IF( nitend-nit000 <= 100 .AND. lwp ) THEN                ! control print (if less than 100 time-step asked)
            WRITE(numout,*) 
            WRITE(numout,*) '        read daily momentum, heat and freshwater fluxes OK'
            DO jf = 1, jpfld_local
               IF( jf == jp_utau .OR. jf == jp_vtau )   zfact =     1.
               IF( jf == jp_qtot .OR. jf == jp_qsr  )   zfact =     0.1
               IF( jf == jp_emp                     )   zfact = 86400.
               WRITE(numout,*) 
               WRITE(numout,*) ' day: ', ndastp , TRIM(sf(jf)%clvar), ' * ', zfact
               CALL prihre( sf(jf)%fnow, jpi, jpj, 1, jpi, 20, 1, jpj, 10, zfact, numout )
            END DO
            CALL FLUSH(numout)
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_flx

   !!======================================================================
END MODULE sbcflx
