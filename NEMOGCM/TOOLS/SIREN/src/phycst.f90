!----------------------------------------------------------------------
! NEMO system team, System and Interface for oceanic RElocable Nesting
!----------------------------------------------------------------------
!
! MODULE: phycst
!
! DESCRIPTION:
!> @brief This module defines physical constant.
!
!> @author
!> J.paul
! REVISION HISTORY:
!> @date November, 2013 - Initial Version
!
!> @note Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
!----------------------------------------------------------------------
MODULE phycst
   USE kind                         ! F90 kind parameter

   IMPLICIT NONE
   ! NOTE_avoid_public_variables_if_possible

   PUBLIC :: dp_pi      !< pi
   PUBLIC :: dp_eps     !< epsilon value
   PUBLIC :: dp_rearth  !< earth radius (km)
   PUBLIC :: dp_deg2rad !< degree to radian ratio 
   PUBLIC :: dp_rad2deg !< radian to degree ratio 
   PUBLIC :: dp_delta   !<  

   REAL(dp), PARAMETER :: dp_pi = 3.14159274101257_dp
   REAL(dp), PARAMETER :: dp_eps = EPSILON(1._dp)
   REAL(dp), PARAMETER :: dp_rearth = 6871._dp
   REAL(dp), PARAMETER :: dp_deg2rad = dp_pi/180.0
   REAL(dp), PARAMETER :: dp_rad2deg = 180.0/dp_pi

   REAL(dp), PARAMETER :: dp_delta=1.e-6
END MODULE phycst

