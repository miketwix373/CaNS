module mod_adv
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_types
  use mod_utils       ,only:store_snap, store_field
  use mod_param       ,only:datadir
  contains
  subroutine adv(dt,dli,u,u_old,uBC,uMean)
  real(rp), intent(in), dimension(3) :: dli
  real(rp), intent(in),dimension(0:,0:, 0:)::  u
  real(rp), intent(in),dimension(:,0:,0:):: u_old
  real(rp), intent(in),dimension(0:,0:):: uMean
  real(rp) :: uBC(0:,0:)
  integer :: i, j, k, n1, n2, nx
  real(rp), intent(in) :: dt
  real(rp):: c
  ! Initialize upast with fixed values
   ! Print the matrices
  n1 = size(uBC,1)
  n2 = size(uBC,2)
  nx = size(u,1)
  do i =0,n1-1
    do j=0,n2-1
            c             = uMean(i,j)*dt/dli(1)
            uBC(i,j)      = u_old(2,i,j)*(1-c) + c*u_old(1,i,j)
    end do
  end do  
  end subroutine adv

subroutine get_Umean(ng, lo, hi, l, dl, n, p, uMean)
    implicit none
    integer, dimension(3) :: coord, startIter, endIter, maxblock
    integer, allocatable :: coord3d(:,:)
    real(rp), intent(in), dimension(0:,0:,0:) :: p ! Import the field
    integer, intent(in), dimension(3) :: ng, lo, hi, n ! Import the local and global number of cells and the limits of each stencil
    real(rp), intent(in), dimension(3) :: dl, l ! Import the dimensions of the domain
    integer :: i, j, k, dim, np, rank ! Create the iterators
    real(rp), allocatable, dimension(:,:) :: uMeanLoc, uMeanGlob
    real(rp), allocatable, dimension(:,:,:) :: uMeanLoc3d
    real(rp), intent(inout) :: uMean(0:,0:)
    integer :: ierr, dimn

    ! Identify MPI communicators
    call MPI_Comm_rank(MPI_COMM_WORLD, rank, ierr)
    call MPI_Comm_size(MPI_COMM_WORLD, np, ierr)

    allocate(uMeanLoc(0:n(2)+1,0:n(3)+1))
    allocate(uMeanLoc3d(0:n(2)+1,0:n(3)+1,np))
    allocate(coord3d(3,np))
    allocate(uMeanGlob(0:ng(2)+1,0:ng(3)+1))


    ! initialize the mean local value
    uMeanLoc = 0.0_rp
    uMeanLoc3d = 0.0_rp
    uMeanGlob = 0.0_rp
    uMean =0.0_rp 

    ! Compute the mean flow in all the stencil
    do j=0,n(2)+1
        do k = 0, n(3)+1
            do i=1,n(1)
                uMeanLoc(j,k) = uMeanLoc(j,k) + p(i,j,k)
            end do
        end do
    end do
    uMeanLoc = uMeanLoc * dl(1) / l(1)

    coord = hi / n
    
    dim = (n(2)+2)*(n(3)+2)

    call MPI_Allgather(uMeanLoc,dim,MPI_REAL_RP, uMeanLoc3d, dim*np,MPI_REAL_RP, MPI_COMM_WORLD,ierr)
    call MPI_Allgather(coord, 3, MPI_INTEGER, coord3d, 3, MPI_INTEGER, MPI_COMM_WORLD, ierr)

    maxblock = maxval(coord3d, dim=2)

    print *, n

    do i=1,np
        startIter = [0,1,1]
        endIter = [0,n(2),n(3)]
        do dim = 2,3
            if (coord3d(dim,i) == 1) then
                startIter(dim) = 0
            end if
            if (coord3d(dim,i) == maxblock(dim)) then
                endIter(dim) = n(dim) + 1
            end if
        end do
        print *, 'startIter'
        print *, startIter
        print *, 'endIter'
        print *, endIter
        print *, 'coord'
        print *, coord3d

        do k = startIter(3), endIter(3)
            do j = startIter(2), endIter(2)
                uMeanGlob((coord3d(2,i)-1)*n(2)+j,(coord3d(3,i)-1)*n(3)+k) = &
                    uMeanGlob((coord3d(2,i)-1)*n(2)+j,(coord3d(3,i)-1)*n(3)+k) + &
                    uMeanLoc3d(j,k,i)        
            end do
        end do
    end do


    uMean = uMeanGlob(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)


end subroutine get_Umean

end module mod_adv
