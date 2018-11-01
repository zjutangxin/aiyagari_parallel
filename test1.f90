! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Test Functions
!  
! Author:
!     Commented by: Xin Tang @ IMF, Summer 2018
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      08/01/2018:                 Original Code
!
! =============================================================================
PROGRAM TEST

IMPLICIT NONE

! -----------------------------------------------------------------------------
!                             DATA DICTIONARY
! -----------------------------------------------------------------------------
integer, parameter :: dp = kind(1.0d0)
integer, parameter :: ni = 14000, nproc = 3
integer :: nii

nii = int(real(ni-1, dp)/real(nproc,dp))+1

write (*,*) 'nii =', nii, 'raw =', int(real(ni, dp)/real(nproc, dp))

END PROGRAM TEST