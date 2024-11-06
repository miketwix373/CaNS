module mod_stats
  use mod_types
  use mpi
  use mod_common_mpi, only: ierr, myid
  use mod_utils, only: trapezoidal_integral, MeanFlow2D
  implicit none
  private
  public mean2D, fluctuations

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

  subroutine fluctuations (p,pFluct,n,ng,lo,hi,dl,l)
    real(rp), dimension(0:,0:,0:,:), intent(in) :: p
    real(rp), dimension(0:,0:,:), intent(inout) :: pFluct
    real(rp), allocatable :: pMeanInstant(:,:,:), flowInstant(:,:,:)
    real(rp), allocatable :: pMean (:,:), pMeanTemp (:,:)
    integer:: nSnaps,i
    integer, intent(in):: n(3), ng(3), lo(3), hi(3)
    real(rp), intent(in):: l(3), dl(3)

    nSnaps = size(p, dim=4)

    allocate(pMeanInstant(0:n(3)+1,0:n(1)+1,nSnaps))
    allocate(flowInstant(0:n(1)+1,0:n(2)+1,0:n(3)+1))
    allocaTE(pMean(0:n(3)+1,0:n(1)+1))
    allocate(pMeanTemp(0:n(3)+1,0:n(1)+1))
    do i=1,nSnaps
      flowInstant = p(:,:,:,i)
      call MeanFlow2D(ng,n,lo,hi,dl,l,flowInstant,pMeanTemp,2,.true.)
      pMeanInstant(:,:,i) = pMeanTemp
    end do

    pMean = sum(pMeanInstant, dim=3)/nSnaps

    do i=1,nSnaps
      pFluct(:,:,i) = pMeanInstant(:,:,i)-pMean
    end do

  end subroutine fluctuations

  subroutine displ_thickness(p,n,ng,lo,hi,dl,l,delta99,zcg)
    real(rp), dimension(0:,0:,0:), intent(in):: p
    real(rp), intent(in) :: dl(3),l(3),zcg(:),delta99(:)
    integer,intent(in) :: n(3),ng(3),lo(3),hi(3)
    real(rp), allocatable :: uinf(:),pMean(:,:),sliceVel(:)
    integer:: i

    allocate(pMean(0:ng(3)+1,0:ng(1)+1))
    allocate(uinf(0:ng(1)+1))
    allocate(sliceVel(0:ng(3)+1))
    call MeanFlow2D(ng,n,lo,hi,dl,l,p,pMean,2,.false.)

    uinf = pMean(ng(3),:)

    do i=1,ng(1)
      sliceVel = pMean(:,i)
      call linear_interp(sliceVel(1:ng(3)),zcg,ng(3),0.99*uinf(i),delta99(i),1)
    end do
  end subroutine displ_thickness

  subroutine get_wss(p,n,ng,lo,hi,dl,l,wss,uTau,zcg,visc)
    real(rp), dimension(0:,0:,0:), intent(in):: p
    real(rp), intent(in) :: dl(3),l(3),zcg(:),visc
    real(rp), intent(inout) ::wss(:),uTau(:)
    integer,intent(in) :: n(3),ng(3),lo(3),hi(3)
    real(rp), allocatable :: uinf(:) ,pMean(:,:),sliceVel(:)
    integer:: i

    allocate(pMean(0:ng(3)+1,0:ng(1)+1))
    allocate(uinf(0:ng(1)+1))
    allocate(sliceVel(0:ng(3)+1))
    call MeanFlow2D(ng,n,lo,hi,dl,l,p,pMean,2,.false.)

    do i=1,ng(1)
      wss(i) = (pMean(2,i)-pMean(1,i))/(zcg(2)-zcg(1))
      uTau(i) = sqrt(visc*wss(i))
    end do
  end subroutine get_wss

end module mod_stats