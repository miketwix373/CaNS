module mod_laminarBL
use mod_types
 implicit none
 private
 public 
contains
    subroutine initBL(npoints,eta_max)
        integer, intent(in) :: n_points
        real(rp), intent(in) :: eta_max
        integer :: i
        real(rp) :: initial_guess,deta,a2,a3,b1,b2,b3
        real(rp), allocatable :: eta, f, fp, fpp
        
        allocate(eta(n_points),f(n_points),fp(n_points),fpp(n_points))
        initial_guess = 0.332057336215195_rp  
        a2 = -5.0_rp/9.0_rp
        a3 = -153.0_rp/128.0_rp
        b1 = 1.0_rp/3.0_rp
        b2 = 15.0_rp/16.0_rp
        b3 = 8.0_rp/15.0_rp  

        deta = eta_max / (n_points - 1)
        do i = 1, n_points
            eta(i) = (i-1) * deta
        end do
        
        ! Solve Blasius equation
        call solve_blasius_ls_rk3(eta, f, fp, fpp, n_points, deta)
        
        ! Output results
        !call output_results(eta, f, fp, fpp, n_points)

    end subroutine(initBL)


    

    subroutine solve_blasius_ls_rk3(eta, f, fp, fpp, n, deta)
        integer, intent(in) :: n
        real(rp), intent(in) :: eta(n), deta
        real(rp), intent(out) :: f(n), fp(n), fpp(n)
        
        ! Local variables
        real(rp) :: y(3), dy(3)
        real(rp) :: q(3)  ! Single storage register for RK stages
        integer :: i
        
        ! Initial conditions
        f(1) = 0.0_rp          ! f(0) = 0
        fp(1) = 0.0_rp         ! f'(0) = 0
        fpp(1) = initial_guess ! f''(0) = guess
        
        ! Low storage RK3 integration
        do i = 1, n-1
            y = [f(i), fp(i), fpp(i)]
            q = 0.0_rp  ! Initialize storage register
            
            ! Stage 1
            call derivatives(y, dy)
            q = deta * dy
            y = y + b1 * q
            
            ! Stage 2
            call derivatives(y, dy)
            q = a2 * q + deta * dy
            y = y + b2 * q
            
            ! Stage 3
            call derivatives(y, dy)
            q = a3 * q + deta * dy
            y = y + b3 * q
            
            ! Store results
            f(i+1) = y(1)
            fp(i+1) = y(2)
            fpp(i+1) = y(3)
        end do
    end subroutine solve_blasius_ls_rk3
    
    subroutine derivatives(y, dydt)
        real(rp), intent(in) :: y(3)
        real(rp), intent(out) :: dydt(3)
        
        ! Blasius equation in first-order form
        ! y(1) = f
        ! y(2) = f'
        ! y(3) = f''
        dydt(1) = y(2)           ! df/dη = f'
        dydt(2) = y(3)           ! df'/dη = f''
        dydt(3) = -y(1)*y(3)     ! df''/dη = -f*f''
    end subroutine derivatives
end module mod_laminarBL