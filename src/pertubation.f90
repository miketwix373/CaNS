module mod_perturbation
    use mod_types
    use mpi
    use mod_common_mpi, only: ierr, myid
    implicit none
    public pert_force
    private gen_wave,pert_param

contains
    subroutine pert_param(uinf,thick0,at,lx,lz,x0,ts,zs)
        real(rp), intent(inout) :: at, lx ,lz , x0, ts, zs
        real(rp), intent(in) :: thick0,uinf

        x0 = 10*thick0
        lx = 4*thick0
        lz = thick0
        ts = 4*thick0/uinf
        zs = 1.7* thick0
        at = 1.0
    end subroutine pert_param
    subroutine pert_force(t,pf,n,dl,zcg,lo,l,past_cycle,init_flag,h0,h1,thick0,uinf)

        real(rp):: at, lx ,lz , x0, ts, zs

        integer, intent(in):: lo(3),n(3)
        integer, intent(inout):: past_cycle
        real(rp), intent(in):: zcg(0:),t, l(3),thick0,uinf,dl(3)
        real(rp), intent(inout) :: pf(0:,0:,0:),h0(:),h1(:)
        logical, intent(in):: init_flag
        real(rp):: p,b,x,z
        real(rp), allocatable:: g(:)
        integer :: cycles,i,j,k

        call pert_param(uinf,thick0,at,lx,lz,x0,ts,zs)

        if (init_flag) then 
            call gen_wave(h0,n,dl(2),lo,l,zs)
            call gen_wave(h1,n,dl(2),lo,l,zs)
        end if


        ! ----- Calculations
        cycles = floor(t/ts)

        if (cycles>past_cycle) then
            h0 = h1
            call gen_wave(h1,n,dl(2),lo,l,zs)
        end if

        p = t/ts - real(cycles)
        b = 3*p**2-2*p**3

        allocate(g(0:n(2)+1))

        g = at*((1-b)*h0 + b*h1)

        pf = 0.0

        do i=1,n(1)
            do j=1,n(2)
                do k=1,n(3)
                    x = (lo(1)-1+i)*dl(1)
                    z = zcg(k)
                    pf(i,j,k) = exp(-(((x-x0)/lx)**2+(z/lz)**2))*g(j)
                end do
            end do
        end do

        past_cycle = cycles
    end subroutine pert_force


    subroutine gen_wave(h,n,dl,lo,l,zs)
        real(rp), parameter :: PI = 3.14159265358979323846
 
        ! ------ Input defintion --------------
        integer, intent(in):: lo(3),n(3)
        real(rp), intent(in):: l(3),zs,dl
        real(rp), intent(inout) :: h(:)
        ! ------ Allocation of local variables 
        real(rp):: amp,phase,kz,y(0:n(2)+1),g(0:n(2)+1)
        integer :: nModes,i
        
        nModes = int(floor(l(2)/zs))
        h = 0.0
        call random_seed()

        do i=0,n(2)
            y(i) = dl*(lo(2)-1+i)
        end do
        
        do i=1,nModes     
            if (myid == 0) then       
                call random_number(amp)
                call random_number(phase)
            end if

            call MPI_BCAST(amp,1,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)
            call MPI_BCAST(phase,1,MPI_REAL_RP,0,MPI_COMM_WORLD,ierr)

            phase   = 2*PI*phase
            kz      = 2.0 * PI * i/ l(2) 
            if (kz < 2.0 * PI/zs) then 
                h = h + amp * cos(kz * y + phase)
            end if
        end do

        h = h / (real(nModes,dp))
    end subroutine gen_wave   

end module mod_perturbation

    