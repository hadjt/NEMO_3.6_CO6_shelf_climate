MODULE trcwri_my_trc
   !!======================================================================
   !!                       *** MODULE trcwri ***
   !!    my_trc :   Output of my_trc tracers
   !!======================================================================
   !! History :   1.0  !  2009-05 (C. Ethe)  Original code
   !!----------------------------------------------------------------------
#if defined key_top && defined key_my_trc && defined key_iomput
   !!----------------------------------------------------------------------
   !!   'key_my_trc'                                           my_trc model
   !!----------------------------------------------------------------------
   !! trc_wri_my_trc   :  outputs of concentration fields
   !!----------------------------------------------------------------------
   USE trc         ! passive tracers common variables 
   USE iom         ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC trc_wri_my_trc 
#if defined key_tracer_budget
   REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:), SAVE :: trb_temp ! slwa
#endif


#  include "top_substitute.h90"
CONTAINS

#if defined key_tracer_budget
   SUBROUTINE trc_wri_my_trc (kt, fl) ! slwa
#else
   SUBROUTINE trc_wri_my_trc
#endif
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE trc_wri_trc  ***
      !!
      !! ** Purpose :   output passive tracers fields 
      !!---------------------------------------------------------------------
#if defined key_tracer_budget
      INTEGER, INTENT( in ), OPTIONAL     :: fl 
      INTEGER, INTENT( in )               :: kt
      REAL(wp), DIMENSION(jpi,jpj,jpk)    :: trpool !tracer pool temporary output
#endif
      CHARACTER (len=20)   :: cltra
      INTEGER              :: jn,jk
      !!---------------------------------------------------------------------
 
      ! write the tracer concentrations in the file
      ! ---------------------------------------


#if defined key_tracer_budget
      IF( PRESENT(fl)) THEN
! depth integrated
! for strict budgetting write this out at end of timestep as an average between 'now' and 'after' at kt
         DO jn = jp_myt0, jp_myt1 
            trpool(:,:,:) = 0.5 * ( trn(:,:,:,jn) * fse3t_a(:,:,:) +  &
                                        trb_temp(:,:,:,jn) * fse3t(:,:,:) )
!
            cltra = TRIM( ctrcnm(jn) )                  ! output of tracer density 
            CALL iom_put( cltra, trpool(:,:,:) / (0.5* (fse3t(:,:,:) + fse3t_a(:,:,:) ) ) )
!
            cltra = TRIM( ctrcnm(jn) )//"_pool"     ! volume integrated output
            DO jk = 1, jpk
               trpool(:,:,jk) = trpool(:,:,jk) * e1t(:,:) * e2t(:,:)
            END DO
            CALL iom_put( cltra, trpool)

!           cltra = TRIM( ctrcnm(jn) )//"_pool"     ! volume integrated output
!           DO jk = 1, jpk
!              trpool(:,:,jk) = 0.5 * ( trn(:,:,jk,jn) * fse3t_a(:,:,jk) +  & 
!                                       trb_temp(:,:,jk,jn) * fse3t(:,:,jk) ) * & 
!                                       e1t(:,:) * e2t(:,:)
!           END DO
!           CALL iom_put( cltra, trpool)
!           cltra = TRIM( ctrcnm(jn) )                  ! output of tracer density 
!           CALL iom_put( cltra, trpool(:,:,:) / (0.5* (fse3t(:,:,:) + fse3t_a(:,:,:) ) ) )
         END DO
         CALL iom_put( "DEPTH" , 0.5* (fse3t(:,:,:) + fse3t_a(:,:,:) ) )  !  equivalent 'depth' at same time as tracer pool output
      ELSE

         IF( kt == nittrc000 ) THEN
           ALLOCATE(trb_temp(jpi,jpj,jpk,jptra))  ! slwa
         ENDIF
         trb_temp(:,:,:,:)=trn(:,:,:,:) ! slwa save for tracer budget (unfiltered trn)

!        DO jn = jp_myt0, jp_myt1
!           cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
!           CALL iom_put( cltra, trn(:,:,:,jn) ) 
!        END DO
! write out depths and areas in double precision for tracer budget calculations
         CALL iom_put( "AREA" , e1t(:,:) * e2t(:,:))
!        CALL iom_put( "DEPTH" , fse3t(:,:,:) )  ! need depth at same time as tracer output

      END IF
#else
      DO jn = jp_myt0, jp_myt1
         cltra = TRIM( ctrcnm(jn) )                  ! short title for tracer
         CALL iom_put( cltra, trn(:,:,:,jn) )
      END DO
#endif
      !
   END SUBROUTINE trc_wri_my_trc

#else
   !!----------------------------------------------------------------------
   !!  Dummy module :                                     No passive tracer
   !!----------------------------------------------------------------------
   PUBLIC trc_wri_my_trc
CONTAINS
   SUBROUTINE trc_wri_my_trc                     ! Empty routine  
   END SUBROUTINE trc_wri_my_trc
#endif

   !!----------------------------------------------------------------------
   !! NEMO/TOP 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE trcwri_my_trc
