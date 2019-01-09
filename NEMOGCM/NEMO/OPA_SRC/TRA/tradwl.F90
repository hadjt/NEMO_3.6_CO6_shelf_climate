MODULE tradwl
   !!======================================================================
   !!                       ***  MODULE  tradwl  ***
   !! Ocean physics: solar radiation penetration in the top ocean levels
   !!======================================================================
   !! History :  POLCOMS  !  1996-10  (J. Holt)  Original code
   !!   NEMO     3.2  !  2010-03  (E. O'Dea)  Import to Nemo for use in Shelf Model
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_dwl      : trend due to the solar radiation penetration
   !!   tra_dwl_init : solar radiation penetration initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean
   USE trc_oce         ! share SMS/Ocean variables
   USE trd_oce         ! trends: ocean variables
   USE trdtra          ! ocean active tracers trends 
   USE in_out_manager  ! I/O manager
   USE phycst          ! physical constants
   USE prtctl          ! Print control
   USE iom             ! I/O manager
   USE fldread         ! read input fields

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_dwl        ! routine called by step.F90 (ln_tradwl=T)

   !                                           !!* Namelist namtra_qsr: penetrative solar radiation
   LOGICAL , PUBLIC ::   ln_tradwl  = .TRUE.    ! light absorption (dwl) flag
   LOGICAL , PUBLIC ::   ln_vary_lambda  = .TRUE.    ! vary Lambda or not flag
   
   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009) 
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_dwl( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!      penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!
      !!          Jason Holt Oct 96
      !!
      !!          Calculates change in temperature due to penetrating
      !!          radiation, with cooling at the surface layer
      !!
      !!          rad=rad0*exp(lambda*z)
      !!
      !!       Heat input into box is between z=K and z=K+1 is RAD(K)-RAD(K+1)
      !!
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!              - save the trend in ttrd ('key_trdtra')
      !!
      !!----------------------------------------------------------------------
      USE oce, ONLY :   ztrdt => ua   ! use ua as 3D workspace   
      USE oce, ONLY :   ztrds => va   ! use va as 3D workspace   
      !!
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !!
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      INTEGER  ::   irgb                 ! temporary integers
      REAL(wp) ::   zchl, zcoef, zsi0r   ! temporary scalars
      REAL(wp) ::   zc0, zc1, zc2, zc3   !    -         -
      REAL(wp), DIMENSION(jpi,jpj)     ::   zekb, zekg, zekr            ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze0, ze1 , ze2, ze3, zea    ! 3D workspace
      !!----------------------------------------------------------------------
      !! HERE GO VARIABLES USED IN POLCOMS CLEAN UP LATER

      integer i,j,k
!      real*8 dtmp(n-1)
      real*8 dtmp(jpkm1)
      real*8 z1,z2,Rad0,Rad1,Rad2,rD,SurfOut,cp
      logical first
      save first
      data first/.true./
      !!--------------------------End of POLCOMS variables Note instead of using saves
      !!--------------------------Could shift this into initial code

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_dwl : penetration of the surface solar radiation'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         CALL tra_dwl_init
         IF( .NOT.ln_tradwl )   RETURN
      ENDIF

      IF( l_trdtra ) THEN      ! Save ta and sa trends
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem)
         ztrds(:,:,:) = 0.e0
      ENDIF
!--------------------------------------------------------------------
!  Set transmissivity
!--------------------------------------------------------------------
!
!  Normal value
!
!---------------------------------------------------------------------------
!
!  Convert Heat fluxes to units used in POL subroutine dwl
!
!---------------------------------------------------------------------------
      cp=3986.0d0

    DO jj = 2, jpj
         DO ji = fs_2, fs_jpim1
           qsr(ji,jj)  = qsr(ji,jj)  * (r1_rau0_rcp)
         ENDDO       !ji
    ENDDO            !jj
!--------------------------------------------------------------------------------
      if ( first ) then
        rlambda2(:,:) = 0.0
        first=.false.
           if ( ln_vary_lambda ) then
!        do j=1,jesub    ! Original Polcoms style Loop
!          do i=1,iesub  ! Original Polcoms style Loop

        do jj=2,jpjm1
          do ji = fs_2, fs_jpim1   ! vector opt.

!            if(ipexb(i,j).ne. 0) then  (Mask, use Tmask instead)

              rlambda2(ji,jj)=-0.033*log(hbatt(ji,jj))+0.2583    ! JIAs formula
              rlambda2(ji,jj)=max(0.05,rlambda2(ji,jj))     ! limit in deep water
          enddo ! ji
        enddo ! jj
      rlambda = 0.0
       else
      rLambda=0.154
       endif ! If vary lambda
      endif ! If first

!      do j=1,jesub      ! Original Polcoms Style Loop
!        do i=1,iesub    ! Original Polcoms Style Loop
              do jk=2,jpk
        DO jj=2,jpjm1
            DO ji = fs_2, fs_jpim1   ! vector opt.

!--------------------------------------------------------------------
! Calculate change in temperature
!--------------------------------------------------------------------
!
!        rad0 = hfl_in(i,j)   ! change hfl_in to qsr I assume

               rad0 = qsr(ji,jj)
                     rD = rLambda2(ji,jj)  +rLambda      !  Transmissivity to be used here
!       if rlambda 0 then rlambda2 not zer and vica versa 

                z2=gdepw_0(ji,jj,jk-1)    ! grid box is from z=z1 to z=z2
                z1=gdepw_0(ji,jj,jk)

                Rad2=Rad0*(exp(-z2*rD)) ! radiation entering box
                Rad1=Rad0*(exp(-z1*rD)) ! radiation leaving box


                dtmp(jk)=1.0/(e3t_0(ji,jj,jk))*(Rad2-Rad1) !change in temperature
                tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + dtmp(jk)
         enddo  ! ji
         enddo  ! jj
             enddo !jk


      IF( l_trdtra ) THEN     ! qsr tracers trends saved for diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         !CEODCALL trd_mod( ztrdt, ztrds, jptra_trd_qsr, 'TRA', kt )
      ENDIF
      !                       ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' qsr  - Ta: ', mask1=tmask, clinfo3='tra-ta' )
      !
   END SUBROUTINE tra_dwl


   SUBROUTINE tra_dwl_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dwl_init  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation for Downwell routine
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (rn_si0,rn_si1) and a ratio
      !!      (rn_abs). These parameters are read in the namtra_qsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov' 
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize rn_si0, rn_si1 and rn_abs
      !!
      !! Reference : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!----------------------------------------------------------------------
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      INTEGER  ::   ios                   ! Local integer output status for namelist read
      INTEGER  ::   irgb, ierror          ! temporary integer
      INTEGER  ::   ioptio, nqsr          ! temporary integer
      REAL(wp) ::   zc0  , zc1            ! temporary scalars
      REAL(wp) ::   zc2  , zc3  , zchl    !    -         -
      REAL(wp) ::   zsi0r, zsi1r, zcoef   !    -         -
      REAL(wp), DIMENSION(jpi,jpj)     ::   zekb, zekg, zekr              ! 2D workspace
      REAL(wp), DIMENSION(jpi,jpj,jpk) ::   ze0 , ze1 , ze2 , ze3 , zea   ! 3D workspace
      !!
      CHARACTER(len=100) ::   cn_dir   ! Root directory for location of ssr files
      TYPE(FLD_N)        ::   sn_chl   ! informations about the chlorofyl field to be read
      NAMELIST/namtra_dwl/  ln_tradwl, ln_vary_lambda
      !!----------------------------------------------------------------------

      REWIND( numnam_ref )            ! Read Namelist namtra_dwl in reference namelist :
      READ  ( numnam_ref, namtra_dwl, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_qsr in reference namelist', lwp )
      
      REWIND( numnam_cfg )            ! Read Namelist namtra_dwl in configuration namelist :
      READ  ( numnam_cfg, namtra_dwl, IOSTAT = ios, ERR = 902)
902   IF( ios /= 0 ) CALL ctl_nam ( ios , 'namtra_qsr in configuration namelist', lwp )
      !
      IF(lwp) THEN                ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_dwl_init : '
         WRITE(numout,*) '~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_dwl : set the parameter of penetration'
         WRITE(numout,*) '      Light penetration (T) or not (F)         ln_tradwl  = ', ln_tradwl
         WRITE(numout,*) '      Vary Lambda  (T) or not (F))             ln_vary_lambda  = ', ln_vary_lambda
      ENDIF

   END SUBROUTINE tra_dwl_init

   !!======================================================================
END MODULE tradwl
