module mod_laminarBL
    use mod_types
    use mod_utils, only: linear_interp
    implicit none
    private
    
    ! Make all public interfaces explicit
<<<<<<< HEAD
    public :: initBL
=======
    public :: initBL, velocityBL, solve_blasius_ls_rk3
    
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
    ! Define constants as parameters
    real(rp), parameter :: INITIAL_GUESS = 0.332057336215195_rp
    real(rp), parameter :: A2 = -5.0_rp/9.0_rp
    real(rp), parameter :: A3 = -153.0_rp/128.0_rp
    real(rp), parameter :: B1 = 1.0_rp/3.0_rp
    real(rp), parameter :: B2 = 15.0_rp/16.0_rp
    real(rp), parameter :: B3 = 8.0_rp/15.0_rp

<<<<<<< HEAD
    contains

    subroutine initBL(n_points,eta_max, thick0, u, z, visc, ubulk,n)
        integer, intent(in) :: n_points,n(3)
        real(rp), intent(in) :: eta_max, thick0, visc, ubulk
        real(rp), dimension(:,0:,0:):: u
        real(rp), intent(in) :: z(:)
=======
contains
    subroutine initBL(n_points, eta_max, thick0, u, z, visc, ubulk)
        integer, intent(in) :: n_points
        real(rp), intent(in) :: eta_max, thick0, visc, ubulk
        real(rp), intent(inout) :: u(:,:,0:), z(:)  ! Fixed array bounds
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
        
        ! Local variables
        integer :: i
        real(rp) :: deta
        real(rp), allocatable :: eta(:), f(:), fp(:), fpp(:)
        
        ! Input validation
        if (n_points <= 1) error stop "n_points must be greater than 1"
        if (eta_max <= 0) error stop "eta_max must be positive"
        
        allocate(eta(n_points), f(n_points), fp(n_points), fpp(n_points), stat=i)
        if (i /= 0) error stop "Allocation failed in initBL"
        
        deta = eta_max / (n_points - 1)
        
        do i = 1, n_points
            eta(i) = (i-1) * deta
        end do
        
        ! Solve Blasius equation
        call solve_blasius_ls_rk3(eta, f, fp, fpp, n_points, deta)
        
        ! Output results
<<<<<<< HEAD
        call velocityBL(eta, f, fp, n_points,thick0, u, z, visc, ubulk)
=======
        call velocityBL(eta, f, fp, n_points, thick0, u, z, visc, ubulk)
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
        
        deallocate(eta, f, fp, fpp)
    end subroutine initBL
    
    subroutine velocityBL(eta, f, fp, n_points, thick0, u, z, visc, ubulk)
        real(rp), dimension(:), intent(in) :: eta, f, fp, z
        real(rp), intent(in) :: thick0, visc, ubulk
<<<<<<< HEAD
        real(rp), intent(inout), dimension(:,0:,0:) :: u  ! Fixed array bounds
=======
        real(rp), intent(inout), dimension(:,:,0:) :: u  ! Fixed array bounds
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
        integer, intent(in) :: n_points
        
        ! Local variables
        real(rp), allocatable :: etaEval(:), fEval(:), fpEval(:)
        integer :: n, i, alloc_stat
        real(rp) :: x
        
        ! Input validation
        if (thick0 <= 0) error stop "thick0 must be positive"
        if (visc <= 0) error stop "visc must be positive"
        
<<<<<<< HEAD
        n = size(z, 1)-1
=======
        n = size(z, 1)
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
        allocate(etaEval(n), fpEval(n), fEval(n), stat=alloc_stat)
        if (alloc_stat /= 0) error stop "Allocation failed in velocityBL"
        
        etaEval = z/thick0
        call linear_interp(eta, f, n_points, etaEval, fEval, n)
        call linear_interp(eta, fp, n_points, etaEval, fpEval, n)
        
        x = thick0**2 * ubulk/visc
        
        ! Boundary check should be added here
        if (x <= 0) error stop "Invalid x calculation in velocityBL"
        
        do i = 1, n
            u(1,:,i) = ubulk * fpEval(i)
            u(2,:,:) = 0.0_rp  ! This assignment looks suspicious
            u(3,:,i) = 0.5_rp * sqrt(visc*ubulk/x) * (etaEval(i)*fpEval(i) - fEval(i))
        end do
        
        deallocate(etaEval, fpEval, fEval)
    end subroutine velocityBL
    
    subroutine solve_blasius_ls_rk3(eta, f, fp, fpp, n, deta)
        integer, intent(in) :: n
        real(rp), intent(in) :: eta(n), deta
        real(rp), intent(out) :: f(n), fp(n), fpp(n)
        
        ! Local variables
        real(rp) :: y(3), dy(3)
        real(rp) :: q(3)  ! Single storage register for RK stages
        integer :: i
        
        ! Input validation
        if (deta <= 0) error stop "deta must be positive"
        
        ! Initial conditions
        f(1) = 0.0_rp    ! f(0) = 0
        fp(1) = 0.0_rp   ! f'(0) = 0
        fpp(1) = INITIAL_GUESS  ! f''(0) = guess
        
        ! Low storage RK3 integration
        do i = 1, n-1
            y = [f(i), fp(i), fpp(i)]
            q = 0.0_rp  ! Initialize storage register
            
            ! Stage 1
            call derivatives(y, dy)
            q = deta * dy
            y = y + B1 * q
            
            ! Stage 2
            call derivatives(y, dy)
            q = A2 * q + deta * dy
            y = y + B2 * q
            
            ! Stage 3
            call derivatives(y, dy)
            q = A3 * q + deta * dy
            y = y + B3 * q
            
            ! Store results
            f(i+1) = y(1)
            fp(i+1) = y(2)
            fpp(i+1) = y(3)
        end do
    end subroutine solve_blasius_ls_rk3
    
<<<<<<< HEAD
    subroutine derivatives(y, dydt)
        real(rp), intent(in) :: y(3)
        real(rp), intent(out) :: dydt(3)
        
=======
    pure subroutine derivatives(y, dydt)
        real(rp), intent(in) :: y(3)
        real(rp), intent(out) :: dydt(3)
        
        ! Blasius equation in first-order form
        ! y(1) = f
        ! y(2) = f'
        ! y(3) = f''
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
        dydt(1) = y(2)        ! df/dη = f'
        dydt(2) = y(3)        ! df'/dη = f''
        dydt(3) = -(y(1)*y(3))*0.5  ! df''/dη = -f*f''
    end subroutine derivatives
<<<<<<< HEAD


=======
>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
end module mod_laminarBL