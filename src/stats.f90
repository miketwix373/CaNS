module mod_stats
  use mod_types
  use mpi
  use mod_common_mpi, only: ierr, myid
  implicit none
  private
  public mean2D

contains

  subroutine mean2D(ng, n, lo, hi, dl, l, p, pMean)
    integer, intent(in) :: ng(3), n(3), lo(3), hi(3)
    integer :: coord(3), i, j, k, nblocks(3), rank, np, startIter(2), endIter(2), dim
    integer, dimension(:,:), allocatable :: coord3, n3d
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:) :: pMean
    real(rp), intent(in), dimension(3) :: dl, l
    real(rp), dimension(:,:), allocatable :: pMeanGlob, pMeanLoc
    real(rp), dimension(:,:,:), allocatable :: pMeanLoc3d
    real(rp) :: localsum

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

    ! Allocate variables
    allocate(pMeanLoc(0:n(2)+1, 0:n(3)+1))
    allocate(pMeanGlob(0:ng(2)+1, 0:ng(3)+1))
    allocate(pMeanLoc3d(0:ng(2)+1, 0:ng(3)+1, np))
    allocate(coord3(3, np))
    allocate(n3d(3, np))

    ! Initialize variables
    pMeanGlob = 0.0_rp
    pMeanLoc = 0.0_rp
    pMeanLoc3d = 0.0_rp
    coord = 0
    coord3 = 0

    ! Compute the block ID positions of each of the stencils
    coord = hi / n

    ! Compute the value of pMean in the stencil face
    do k = 0, n(3)+1
      do j = 0, n(2)+1
        localsum = 0.0_rp
        do i = 1, n(1)
          localsum = localsum + p(i,j,k) * dl(1) / l(1)
        end do
        pMeanLoc(j,k) = localsum
      end do
    end do

    call MPI_Allgather(pMeanLoc, (n(2)+2)*(n(3)+2), MPI_REAL_RP, &
                       pMeanLoc3d, (n(2)+2)*(n(3)+2), MPI_REAL_RP, &
                       MPI_COMM_WORLD, ierr)
    call MPI_Allgather(coord, 3, MPI_INTEGER, coord3, 3, MPI_INTEGER, &
                       MPI_COMM_WORLD, ierr)
    call MPI_Allgather(n, 3, MPI_INTEGER, n3d, 3, MPI_INTEGER, &
                       MPI_COMM_WORLD, ierr)

    nblocks = maxval(coord3, dim=2)

    do i = 1, np
      startIter = 1
      endIter = n3d(2:3, i)
      do dim = 1, 2
        if (coord3(dim+1, i) == 1) then
          startIter(dim) = 0
        end if
        if (coord3(dim+1, i) == nblocks(dim+1)) then
          endIter(dim) = endIter(dim) + 1
        end if
      end do
      do k = startIter(2), endIter(2)
        do j = startIter(1), endIter(1)
          pMeanGlob((coord3(2,i)-1)*n3d(2,i)+j, (coord3(3,i)-1)*n3d(3,i)+k) = &
            pMeanGlob((coord3(2,i)-1)*n3d(2,i)+j, (coord3(3,i)-1)*n3d(3,i)+k) + &
            pMeanLoc3d(j,k,i)
        end do
      end do
    end do

    pMean = pMeanGlob(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)

    ! Deallocate variables
    deallocate(pMeanLoc, pMeanGlob, pMeanLoc3d, coord3, n3d)

  end subroutine mean2D

end module mod_stats