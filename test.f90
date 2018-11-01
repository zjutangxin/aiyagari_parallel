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

! Model Parameters
    real(dp), parameter :: beta = 0.96_dp
    real(dp), parameter :: alpha = 0.36_dp
    real(dp), parameter :: gamma = 2.0_dp
    real(dp), parameter :: delta = 0.1_dp
    
! Computational Parameters
    integer, parameter :: ne = 7
    integer, parameter :: na = 500
    real(dp), parameter :: eps = 1.0d-9        
    real(dp), parameter :: a_inf = eps
    real(dp), parameter :: a_sup = 50.0_dp
    real(dp), parameter :: grid_curv = 2.0_dp
    real(dp), dimension(na) :: agrid, agrid_temp
    real(dp) :: K, N
    
! Partial Equilibrium Prices
    real(dp), parameter :: r = 1.0_dp/beta-1
    
! Variables
    integer, parameter :: ni = na*ne
    real(dp), dimension(ne) :: lgrid, invp
    real(dp), dimension(ne,ne) :: pi_l, pi_inv    
    real(dp), dimension(ni) :: dist
    integer, dimension(1) :: adummy
    
! Auxilliary Variables
    integer, dimension(ni) :: iafun, iefun
    integer, dimension(ne, na) :: iifun
    integer :: inda, inde, indi, indtr
    
! -----------------------------------------------------------------------------
!                             PREPARATION
! -----------------------------------------------------------------------------

! Earnings Shock
lgrid = (/ -1.20014, -0.75736, -0.36941, 0.00000, & 
    0.36941, 0.75736, 1.20014 /)
pi_l(1,:) = (/ 0.2019378, 0.5034201, 0.2568107, 0.0363963, 0.0014232, &
    0.0000119, 0.0000000/)
pi_l(2,:) = (/ 0.0411231, 0.3233636, 0.4513294, 0.1667877, 0.0170059, &
    0.0003894,0.0000010 /)
pi_l(3,:) = (/ 0.0057942, 0.1246568, 0.4202567, 0.3596454, 0.0849178, &
    0.0046970,0.0000321 /)
pi_l(4,:) = (/ 0.0005483, 0.0307571, 0.2401232, 0.4571429, 0.2401232, &
    0.0307571,0.0005483 /)
pi_l(5,:) = (/ 0.0000321, 0.0046970, 0.0849178, 0.3596454, 0.4202567, &
    0.1246568, 0.0057942 /)
pi_l(6,:) = (/ 0.0000010, 0.0003894, 0.0170059, 0.1667877, 0.4513294, &
    0.3233636, 0.0411231 /)
pi_l(7,:) = (/ 0.0000000, 0.0000119, 0.0014232, 0.0363963, 0.2568107, &
    0.5034201, 0.2019378 /)

! Stationary Distribution
pi_inv = pi_l
do indtr = 1,500
    pi_inv = MATMUL(pi_inv, pi_l)
enddo
invp = pi_inv(1,:)
! do indtr = 1,ne
!     write (*,*) pi_inv(indtr,:)
! enddo

write (*,*) 'Stationary Distribution:'
write (*,*) invp

! Asset Grid
agrid(1) = a_inf
agrid(na) = a_sup
do inda = 2,na-1
    agrid(inda) = a_inf + (a_sup-a_inf)*&
        (real(inda-1,dp)/real(na-1,dp))**grid_curv
enddo

! Conversion Function between i and (e,a)
indi = 0
do inde = 1,ne
    do inda = 1,na
        indi = indi + 1
        iifun(inde,inda) = indi
        iefun(indi) = inde
        iafun(indi) = inda
    enddo
enddo

! Steady State Capital
N = sum(exp(lgrid)*invp)
K = ((r+delta)/alpha)**(1.0_dp/(alpha-1.0_dp))*N
write (*,*) 'K=', K, 'N=', N

agrid_temp = abs(K-agrid)
adummy = minloc(agrid_temp)
write (*,*) 'grid_loc=', adummy(1)

dist(adummy(1):ni:na) = invp

K = 0.0_dp
do inde = 1,ne
    K = K + sum(agrid*dist((inde-1)*na+1:inde*na))
enddo
write (*,*) 'summed capital = ', K

open(1,file='dist.txt',form='formatted')
do indi = 1,ni
    write(1,101) dist(indi)
101 format(f15.9)
enddo        
close(1)

! Write Intermediate Results to Hard Drive
! open(1,file='agrid.txt',form='formatted')
! do inda = 1,na
!     write(1,101) agrid(inda)
! 101 format(f15.9)
! enddo        
! close(1)

! open(1,file='shock_pi.txt',form='formatted')
! do inde = 1,ne
!     write(1,102) pi_l(inde,:)
! 102 format(7f12.6)
! enddo        
! close(1)

! open(1,file='inv_dist.txt',form='formatted')
!     write(1,102) lgrid     
!     write(1,102) invp
! close(1)

END PROGRAM TEST