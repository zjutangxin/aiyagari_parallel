! -----------------------------------------------------------------------------
!                             PROGRAM DESCRIPTION
! -----------------------------------------------------------------------------
!   
! Purpose:
!     - Main function replicating Aiyagari (1994)
!     - Parallel Computing
!  
! Author:
!     Xin Tang @ IMF, Summer 2018
!  
! Record of Revisions:
!         Date:                 Description of Changes
!     ===========        =================================
!      08/02/2018:                 Original Code
!
! Compiler Used:
!   Intel(R) Visual Fortran Compiler XE 14.0.0.103 [IA-64]
!   Commercial Release: Intel Parallel Studio XE 2013
!   Integrated with Microsoft Visual Studio 2012
!
! Library Used:
!   Intel MPI
! =============================================================================

PROGRAM Aiyagari_GE_Pa

IMPLICIT NONE

include 'mpif.h'

! -----------------------------------------------------------------------------
!                             DATA DICTIONARY
! -----------------------------------------------------------------------------
integer, parameter :: dp = kind(1.0d0)
! Model Parameters
    real(dp), parameter :: beta = 0.96_dp
    real(dp), parameter :: alpha = 0.36_dp
    real(dp), parameter :: gamma = 2.0_dp
    real(dp), parameter :: delta = 0.10_dp
    
! Computational Parameters
    integer, parameter :: ne = 7
    integer, parameter :: na = 1000
    real(dp), parameter :: eps = 1.0d-9    
    real(dp), parameter :: a_inf = eps
    real(dp), parameter :: a_sup = 50.0_dp
    real(dp), parameter :: grid_curv = 2.0_dp
    real(dp), dimension(na) :: agrid
    
! General Equilibrium Prices and Aggregates
    real(dp) :: r, w, ky0, ky1, kl, cap, lab
    
! Variables
    integer, parameter :: ni = na*ne
    real(dp), dimension(ni) :: val0, val1, polc, pola, dist0, dist1
    real(dp), dimension(ne, na) :: vale
    integer, dimension(ni) :: astar
    real(dp), dimension(ne) :: lgrid, invp
    real(dp), dimension(ne,ne) :: pi_l, pi_inv    
    
! Auxilliary Variables
    integer, dimension(ni) :: iafun, iefun
    integer, dimension(ne, na) :: iifun
    integer :: inda, inde, inde1, indi, iter_p, iter_v, iter_d
    integer :: indtr
    integer, dimension(1) :: adummy, ddummy
    integer, parameter :: maxiter_p = 200
    integer, parameter :: maxiter_v = 2000
    integer, parameter :: maxiter_d = 5000
    real(dp), parameter :: maxerr_p = 1.0d-3
    real(dp), parameter :: maxerr_v = 1.0d-5
    real(dp), parameter :: maxerr_d = 1.0d-9
    real(dp), parameter :: update_p = 0.05_dp
    real(dp) :: tempc, etime, btime, err_p, err_v, err_d
    real(dp), dimension(na) :: cvec, vvec
    logical, dimension(na) :: cfea

! Parallel Computation Variables    
    integer, parameter :: root = 0
    integer :: ierr, myrank, nproc, ni_i, ibegin, iend
    real(dp), dimension(:), allocatable :: val_ind, val_agg
    real(dp), dimension(:), allocatable :: astar_ind, astar_agg    
    real(dp), dimension(ni) :: dist_temp
    
! -----------------------------------------------------------------------------
!                             PREPARATION
! -----------------------------------------------------------------------------
! Initialize MPI environment
call MPI_INIT(ierr)
call MPI_COMM_SIZE(mpi_comm_world, nproc, ierr)
call MPI_COMM_RANK(mpi_comm_world, myrank, ierr)

call cpu_time(btime)

! Distribution of Tasks
ni_i = int(real(ni-1, dp)/real(nproc, dp))+1
allocate(val_ind(ni_i), val_agg(ni_i*nproc))
allocate(astar_ind(ni_i), astar_agg(ni_i*nproc))
ibegin = myrank*ni_i +1
iend = min((myrank+1)*ni_i,ni)

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

! Asset Grid
agrid(1) = a_inf
agrid(na) = a_sup
do inda = 2,na-1
    agrid(inda) = a_inf + (a_sup-a_inf)* &
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

! -----------------------------------------------------------------------------
!                             LOOP OVER PRICES
! -----------------------------------------------------------------------------
! Initialization
ky0 = 2.7_dp
kl = ky0**(1.0_dp/(1.0_dp - alpha))
r = alpha*kl**(alpha-1.0_dp) - delta
w = (1.0_dp-alpha)*kl**alpha

do indi = 1,ni
    tempc = max(exp(lgrid(iefun(indi)))*w + agrid(iafun(indi))*(1+r) &
        - a_inf, eps)
    val0(indi) = tempc**(1.0_dp - gamma)/(1.0_dp - gamma)
enddo

dist0 = 0.0_dp
! Initialize dist0, everyone with representative ss. capital
! Labor distribution as invariant distribution
lab = sum(exp(lgrid)*invp)
cap = (((1.0_dp/beta-1)+delta)/alpha)**(1.0_dp/(alpha-1.0_dp))*lab
ddummy = minloc(abs(cap-agrid))
dist0(ddummy(1):ni:na) = invp

err_p = 1000.0_dp
iter_p = 1

do
    ! ===================== Value Function Iteration ==========================
    err_v = 1000.0_dp
    iter_v = 1
    
    do
        ! Construct Continuation Value
        vale = 0.0_dp
        do inde = 1,ne
            do inda = 1,na
                do inde1 = 1,ne
                    vale(inde,inda) = vale(inde,inda) + &
                        pi_l(inde,inde1)*val0(iifun(inde1,inda))
                enddo
            enddo
        enddo
    
        ! Iterate over Current State
        ! Parallel across Asset Space
        do indi = ibegin, iend
            cvec = exp(lgrid(iefun(indi)))*w + agrid(iafun(indi))*(1+r) &
                - agrid
            cfea = (cvec > eps)
            where (cfea .eqv. .TRUE.)
                vvec = cvec**(1.0_dp-gamma)/(1.0_dp-gamma) + &
                 beta*vale(iefun(indi),:)
            end where
            adummy = maxloc(vvec, cfea)
            astar_ind(indi-ibegin+1) = adummy(1)
            val_ind(indi-ibegin+1) = vvec(adummy(1))
        enddo
    
        ! Gather val\_ind from all processors
        call mpi_allgather(val_ind(1),ni_i,mpi_double_precision,&
            val_agg(1), ni_i, mpi_double_precision, mpi_comm_world, ierr)
        call mpi_allgather(astar_ind(1),ni_i,mpi_double_precision,&
            astar_agg(1), ni_i, mpi_double_precision, mpi_comm_world, ierr)
        val1(1:ni) = val_agg(1:ni)
        astar(1:ni) = astar_agg(1:ni)
        
        err_v = maxval(abs(val0-val1))
        val0 = val1
        iter_v = iter_v + 1
    
        if ((mod(iter_v,100) == 0) .and. (myrank == root)) then
            write (*,'(A7,I4,A10,F12.6)') 'Iter_v=', iter_v, 'err_v=', err_v    
        endif
    
        if (err_v < maxerr_v) exit
        if (iter_v > maxiter_v) then 
            call mpi_finalize(ierr)
            stop
        endif
        
    enddo

    ! ===================== Stationary Distribution =========================
    err_d = 1000.0_dp
    iter_d = 1
    
    do
        dist1 = 0.0_dp
        
        !!! Sequential Version
        !do indi = 1,ni
        !    do inde1 = 1,ne
        !        dist1(iifun(inde1,astar(indi))) = & 
        !            dist1(iifun(inde1,astar(indi))) &
        !            + pi_l(iefun(indi),inde1)*dist0(indi)
        !    enddo
        !enddo
        
        ! Parallel Version        
        do indi = ibegin, iend
            do inde1 = 1,ne
                dist1(iifun(inde1,astar(indi))) = & 
                    dist1(iifun(inde1,astar(indi))) &
                    + pi_l(iefun(indi),inde1)*dist0(indi)
            enddo
        enddo
        call mpi_allreduce(dist1(1),dist_temp(1),ni,mpi_double_precision, &
            mpi_sum, mpi_comm_world, ierr)
        dist1 = dist_temp
        
        err_d = maxval(abs(dist0-dist1))
        dist0 = dist1
        iter_d = iter_d + 1
        
        if ((mod(iter_d,100) == 0) .and. (myrank == root)) then
            write (*,'(A7,I4,A10,F14.10)') 'Iter_d=', iter_d, 'err_d=', err_d
        endif
    
        if (err_d < maxerr_d) exit
        if (iter_d > maxiter_d) stop        
    enddo
    ! ===================== Aggregating Results =============================
    cap = 0.0_dp
    do inde = 1,ne
        cap = cap + sum(agrid*dist0((inde-1)*na+1:inde*na))
    enddo
    ky1 = (cap/lab)**(1.0_dp - alpha)
    
iter_p = iter_p + 1    
err_p = abs(ky0-ky1)
if (err_p < maxerr_p) exit
if (iter_p > maxiter_p) stop

if (myrank == root) then
    write (*,*) '===================================================='
    write (*,'(A7,I4,A10,F12.6)') 'Iter_p=', iter_p, 'err_p =', err_p
endif

! Update Guess
ky0 = (1-update_p)*ky0 + update_p*ky1
kl = ky0**(1.0_dp/(1.0_dp - alpha))
r = alpha*kl**(alpha-1.0_dp) - delta
w = (1.0_dp-alpha)*kl**alpha

enddo

if (myrank == root) then
    write (*,'(A10, F7.4, A12, F7.4, A8, F7.4)') 'capital =', cap, &
            'interest =', 1+r, 'wage =', w
    write (*,'(A10, F7.4)') 'kyratio= ', ky0
endif

call cpu_time(etime)
if (myrank == root) then
    write (*,'(A10, F6.2)') 'Time =', etime - btime
endif

! -----------------------------------------------------------------------------
!                             EXPORT RESULTS
! -----------------------------------------------------------------------------
if (myrank == root) then
    open(1,file='agrid.txt',form='formatted')
    do inda = 1,na
        write(1,101) agrid(inda)
    101 format(ES14.6)
    enddo        
    close(1)

    open(1,file='shock_pi.txt',form='formatted')
    do inde = 1,ne
        write(1,102) pi_l(inde,:)
    102 format(7f12.6)
    enddo        
    close(1)

    open(1,file='shock_state.txt',form='formatted')
        write(1,102) lgrid     
    close(1)

    open(1,file='policies.txt',form='formatted')
    103 format(ES14.6, I2, 2ES14.6, ES14.6, ES14.6)
    do indi = 1,ni
        pola(indi) = agrid(astar(indi)) 
        polc(indi) = exp(lgrid(iefun(indi)))*w + agrid(iafun(indi))*(1+r) &
                - pola(indi)
        write(1,103) agrid(iafun(indi)), iefun(indi), polc(indi), pola(indi), & 
            dist0(indi), val0(indi)
    enddo
    close(1)

    open(1,file='distribution.txt',form='formatted')
    105 format(ES14.6)
    do indi = 1,ni
        write(1,105) dist0(indi)
    enddo
    close(1)    

    open(1,file='value_function.txt',form='formatted')
    do indi = 1,ni
        write(1,105) val0(indi)
    enddo
    close(1)    

    open(1,file='ss_compu.txt',form='formatted')
        write(1,104) na, ne
    104 format(I5, I2)    
    close(1)    

    open(1,file='ss_vec.txt',form='formatted')
    write(1,106) cap, ky0, delta, lab
    106 format(4ES14.6)    
    close(1) 

endif    

call mpi_finalize(ierr)

END PROGRAM Aiyagari_GE_Pa