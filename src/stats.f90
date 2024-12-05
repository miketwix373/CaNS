module mod_stats
<<<<<<< HEAD

use mod_types
use mpi
use mod_common_mpi, only: ierr, myid
use mod_utils, only: trapezoidal_integral, MeanFlow2D, sort_couple, linear_interp

implicit none
private
public mean2D, fluctuations, bl_stats, displ_thickness

character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
character(len=*), parameter :: fmt_rp = fmt_dp
#else
character(len=*), parameter :: fmt_rp = fmt_sp
#endif

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

  allocate(pMeanLoc(0:n(2)+1, 0:n(3)+1))
  allocate(pMeanGlob(0:ng(2)+1, 0:ng(3)+1))
  allocate(pMeanLoc3d(0:ng(2)+1, 0:ng(3)+1, np))
  allocate(coord3(3, np))
  allocate(n3d(3, np))

  pMeanGlob = 0.0_rp
  pMeanLoc = 0.0_rp
  pMeanLoc3d = 0.0_rp
  coord = 0
  coord3 = 0

  coord = hi / n

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
        pMeanLoc3d, (n(2)+2)*(n(3)+2), MPI_REAL_RP, MPI_COMM_WORLD, ierr)
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
         pMeanGlob((coord3(2,i)-1)*n3d(2,i)+j, (coord3(3,i)-1)*n3d(3,i)+k) &
         + pMeanLoc3d(j,k,i)
      end do
    end do
  end do

  pMean = pMeanGlob(lo(2)-1:hi(2)+1,lo(3)-1:hi(3)+1)
  deallocate(pMeanLoc, pMeanGlob, pMeanLoc3d, coord3, n3d)
end subroutine mean2D

subroutine fluctuations(p, pFluct, n, ng, lo, hi, dl, l)
  real(rp), dimension(0:,0:,0:,:), intent(in) :: p
  real(rp), dimension(0:,0:,:), intent(inout) :: pFluct
  real(rp), allocatable :: pMeanInstant(:,:,:), flowInstant(:,:,:)
  real(rp), allocatable :: pMean(:,:), pMeanTemp(:,:)
  integer :: nSnaps, i
  integer, intent(in) :: n(3), ng(3), lo(3), hi(3)
  real(rp), intent(in) :: l(3), dl(3)
  nSnaps = size(p, dim=4)
  allocate(pMeanInstant(0:n(3)+1,0:n(1)+1,nSnaps))
  allocate(flowInstant(0:n(1)+1,0:n(2)+1,0:n(3)+1))
  allocate(pMean(0:n(3)+1,0:n(1)+1))
  allocate(pMeanTemp(0:n(3)+1,0:n(1)+1))

  pMeanInstant = 0.0_rp
  flowInstant = 0.0_rp
  pMean = 0.0_rp
  pMeanTemp = 0.0_rp

  do i=1,nSnaps
    flowInstant = p(:,:,:,i)
    call MeanFlow2D(ng,n,lo,hi,dl,l,flowInstant,pMeanTemp,2,.true.)
    pMeanInstant(:,:,i) = pMeanTemp
  end do
  pMean = sum(pMeanInstant, dim=3)/nSnaps
  do i=1,nSnaps
    pFluct(:,:,i) = pMeanInstant(:,:,i)-pMean
  end do

  deallocate(pMeanInstant, flowInstant, pMean, pMeanTemp)

end subroutine fluctuations

subroutine displ_thickness(p, n, ng, lo, hi, dl, l, delta99, zcg, uinf)
  real(rp), dimension(0:,0:,0:), intent(in) :: p
  real(rp), intent(in) :: dl(3), l(3), zcg(:)
  integer, intent(in) :: n(3), ng(3), lo(3), hi(3)
  real(rp), intent(out) :: delta99(:), uinf(:)
  real(rp), allocatable :: pMean(:,:), sliceVel(:), velTemp(:), thickTemp(:)
  integer :: i
  allocate(pMean(0:ng(3)+1, 0:ng(1)+1))
  allocate(sliceVel(0:ng(3)+1))
  allocate(velTemp(1), thickTemp(1))

  pMean = 0.0_rp
  sliceVel = 0.0_rp
  velTemp = 0.0_rp
  thickTemp = 0.0_rp

  print *, ng(1)
  call MeanFlow2D(ng, n, lo, hi, dl, l, p, pMean, 2, .false.)
  do i = 1, ng(1)
    sliceVel = pMean(:,i)
    velTemp(1) = 0.99_rp * pMean(ng(3),i)
    call linear_interp(sliceVel(1:ng(3)), zcg, ng(3), velTemp, thickTemp, 1)
    delta99(i) = thickTemp(1)
  end do

  deallocate(pMean, sliceVel, velTemp, thickTemp)
end subroutine displ_thickness

subroutine get_wss(p, n, ng, lo, hi, dl, l, wss, uTau, zcg, visc)
  real(rp), dimension(0:,0:,0:), intent(in) :: p
  real(rp), intent(in) :: dl(3), l(3), zcg(:), visc
  real(rp), intent(inout) :: wss(:), uTau(:)
  integer, intent(in) :: n(3), ng(3), lo(3), hi(3)
  real(rp), allocatable :: uinf(:), pMean(:,:), sliceVel(:)
  integer :: i
  allocate(pMean(0:ng(3)+1,0:ng(1)+1))
  allocate(uinf(0:ng(1)+1))
  allocate(sliceVel(0:ng(3)+1))

  pMean = 0.0_rp
  uinf = 0.0_rp
  sliceVel = 0.0_rp

  call MeanFlow2D(ng,n,lo,hi,dl,l,p,pMean,2,.false.)
  do i=1,ng(1)
    wss(i) = (pMean(2,i)-pMean(1,i))/(zcg(2)-zcg(1))
    uTau(i) = sqrt(abs(visc*wss(i))) * sign(1.0_rp, wss(i))
  end do
  deallocate(pMean, uinf, sliceVel)
end subroutine get_wss

subroutine bl_stats(fname1d, datadir, n, ng, lo, hi, dl, l, u, v, w, zcg, visc, myid, fldnum)
  character(len=*) :: fname1d, datadir, fldnum
  character(len=100) :: fname2d
  integer, intent(in) :: n(3), ng(3), lo(3), hi(3), myid
  real(rp), intent(in) :: visc, dl(:), l(3)
  real(rp), intent(in) :: u(0:,0:,0:), v(0:,0:,0:), w(0:,0:,0:), zcg(:)
  integer :: iunit, rank, node
  real(rp), allocatable :: delta99(:), wss(:), uTau(:), lscale(:), delta99plus(:)
  real(rp), allocatable :: yPlus(:,:), uPlus(:,:), uMean(:,:), uInf(:)
  real(rp), allocatable :: yPlusTemp(:), uPlusTemp(:), integrand(:)
  real(rp), allocatable :: reDelta(:), reTheta(:), velTemp(:), thickTemp(:)
  integer :: i, npoints, idLimit, j
  allocate(delta99(ng(1)), wss(ng(1)), uTau(ng(1)), lscale(ng(1)))
  allocate(yPlus(ng(3),ng(1)), uPlus(ng(3),ng(1)))
  allocate(uMean(0:ng(3)+1,0:ng(1)+1))
  allocate(uInf(ng(1)))
  allocate(yPlusTemp(ng(3)), uPlusTemp(ng(3)))
  allocate(reDelta(ng(1)), reTheta(ng(1)))
  allocate(delta99plus(ng(1)))
  allocate(velTemp(1), thickTemp(1))

  yPlusTemp = 0.0_rp
  uPlusTemp = 0.0_rp
  delta99 = 0.0_rp
  wss = 0.0_rp
  uTau = 0.0_rp
  uMean = 0.0_rp

  call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD, node, ierr)
  call Meanflow2D(ng, n, lo, hi, dl, l, u, uMean, 2, .false.)
  call displ_thickness(u, n, ng, lo, hi, dl, l, delta99, zcg, uInf)
  call get_wss(u, n, ng, lo, hi, dl, l, wss, uTau, zcg, visc)

  where (uTau /= 0.0_rp)
    lscale = visc/abs(uTau)
    delta99plus = delta99/lscale
  elsewhere
    lscale = 0.0_rp
    delta99plus = 0.0_rp
  end where

  do i = 1, ng(1)
    if (uTau(i) /= 0.0_rp) then
      yPlus(:,i) = zcg(1:ng(3))/lscale(i)
      uPlus(:,i) = uMean(1:ng(3),i)/uTau(i)
      yPlusTemp = yPlus(:,i)
      uPlusTemp = uPlus(:,i)
      idLimit = count(yPlusTemp <= delta99(i))
        if (idLimit > 0) then
          allocate(integrand(idLimit))
          integrand = 1.0_rp - (uPlusTemp(1:idLimit)/(uInf(i)/uTau(i)))
          call trapezoidal_integral(yPlusTemp(1:idLimit), integrand, idLimit, reDelta(i))
          integrand = (1.0_rp - (uPlusTemp(1:idLimit)/(uInf(i)/uTau(i)))) * (uPlusTemp(1:idLimit)/(uInf(i)/uTau(i)))
          call trapezoidal_integral(yPlusTemp(1:idLimit), integrand, idLimit, reTheta(i))
          deallocate(integrand)
        else
          reDelta(i) = 0.0_rp
          reTheta(i) = 0.0_rp
        end if
    else
      yPlus(:,i) = 0.0_rp
      uPlus(:,i) = 0.0_rp
      reDelta(i) = 0.0_rp
      reTheta(i) = 0.0_rp
    end if
  end do

  if (myid == 0) then
    open(newunit=iunit,file=fname1d)
    do i=1,ng(1)
      write(iunit,fmt_rp) wss(i), reDelta(i), reTheta(i)
    end do
  close(iunit)

  fname2d = trim(datadir)//'yPlus'//fldnum//'.out'

  open(newunit=iunit,file=fname2d)
  do i=1,ng(3)
    write(iunit,fmt_rp) yPlus(i,:)
  end do
  close(iunit)

  fname2d = trim(datadir)//'uPlus'//fldnum//'.out'

  open(newunit=iunit,file=fname2d)
  do i=1,ng(3)
    write(iunit,fmt_rp) uPlus(i,:)
  end do
  close(iunit)
  end if
  deallocate(delta99, wss, uTau, lscale, delta99plus)
  deallocate(yPlus, uPlus, uMean, uInf)
  deallocate(yPlusTemp, uPlusTemp)
  deallocate(reDelta, reTheta)
  deallocate(velTemp, thickTemp)
=======
  use mod_types
  use mpi
  use mod_common_mpi, only: ierr, myid
  use mod_utils, only: trapezoidal_integral, MeanFlow2D,sort_couple, linear_interp
  implicit none
  private
  public mean2D, fluctuations, bl_stats, displ_thickness
  character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                 fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
  character(len=*), parameter :: fmt_rp = fmt_dp
#else
  character(len=*), parameter :: fmt_rp = fmt_sp
#endif
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

subroutine displ_thickness(p, n, ng, lo, hi, dl, l, delta99, zcg,uinf)
    real(rp), dimension(0:,0:,0:), intent(in) :: p
    real(rp), intent(in) :: dl(3), l(3), zcg(:)
    integer, intent(in) :: n(3), ng(3), lo(3), hi(3)
    real(rp), intent(out) :: delta99(:),uinf(:)
    
    ! Local variables
    real(rp), allocatable :: pMean(:,:), sliceVel(:), velTemp(:), thickTemp(:)
    integer :: i
    
    ! Allocate arrays
    allocate(pMean(0:ng(3)+1, 0:ng(1)+1))
    allocate(sliceVel(0:ng(3)+1))
    allocate(velTemp(1), thickTemp(1))
    
    ! Calculate mean flow field
    call MeanFlow2D(ng, n, lo, hi, dl, l, p, pMean, 2, .false.)
    
    ! Get free-stream velocity
    uinf = pMean(ng(3),1:n(1))
    
    ! Calculate boundary layer thickness for each streamwise location
    do i = 1, ng(1)
        ! Extract velocity profile at current x-location
        sliceVel = pMean(:,i)
        
        ! Target velocity for δ99 (99% of free-stream)
        velTemp(1) = 0.99_rp * uinf(i)
        
        ! Interpolate to find δ99
        call linear_interp(sliceVel(1:ng(3)), zcg, ng(3), velTemp, thickTemp, 1)
        
        delta99(i) = thickTemp(1)
    end do
    
    ! Deallocate arrays
    deallocate(pMean, sliceVel, velTemp, thickTemp)
    
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

  subroutine bl_stats(fname1d,datadir, n, ng, lo, hi, dl, l, u, v, w, zcg, visc,myid,fldnum)
    character(len=*) :: fname1d,datadir,fldnum
    character(len=100) :: fname2d
    integer, intent(in) :: n(3), ng(3), lo(3), hi(3), myid
    real(rp), intent(in) :: visc, dl(:), l(3)
    real(rp), intent(in) :: u(0:,0:,0:), v(0:,0:,0:), w(0:,0:,0:), zcg(:)
    integer:: iunit, rank, node
    ! Local variables
    real(rp), allocatable :: delta99(:), wss(:), uTau(:), lscale(:), delta99plus(:)
    real(rp), allocatable :: yPlus(:,:), uPlus(:,:), uMean(:,:), uInf(:)
    real(rp), allocatable :: yPlusTemp(:), uPlusTemp(:), integrand(:)
    real(rp), allocatable :: reDelta(:), reTheta(:),velTemp(:), thickTemp(:)
    integer :: i, npoints, idLimit, j
    real(rp) :: delta99u


    ! Fix 2: Correct allocation sizes
    allocate(delta99(ng(1)), wss(ng(1)), uTau(ng(1)), lscale(ng(1)), delta99plus(ng(1)))
    allocate(yPlus(ng(3),ng(1)))
    allocate(uPlus(ng(3),ng(1)), uMean(0:ng(3)+1,0:ng(1)+1), uInf(ng(1)))
    allocate(yPlusTemp(ng(3)+1), uPlusTemp(ng(3)+1))
    allocate(reDelta(ng(1)), reTheta(ng(1)))
    allocate(velTemp(1), thickTemp(1))

    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, node, ierr)
    ! Fix 3: Initialize arrays
    yPlusTemp = 0.0_rp
    uPlusTemp = 0.0_rp
    delta99 = 0.0_rp
    wss = 0.0_rp
    uTau = 0.0_rp

    ! Calculate mean flow and boundary layer parameters
    call Meanflow2D(ng, n, lo, hi, dl, l, u, uMean, 2, .false.)
    call displ_thickness(u, n, ng, lo, hi, dl, l, delta99, zcg,uInf)
    call get_wss(u, n, ng, lo, hi, dl, l, wss, uTau, zcg, visc)
    lscale = visc/uTau
    delta99plus = delta99/lscale

    do i = 1, ng(1)
        ! Fix 5: Clear temporary arrays for each iteration
        uPlusTemp = 0.0_rp
        yPlusTemp = 0.0_rp
        
        yPlus(:,i) = zcg(1:n(3))/lscale(i)
        uPlus(:,i) = uMean(1:ng(3),i)/uTau(i)  ! Fix 6: Correct array indexing
        
        yPlusTemp(1:ng(3)) = yPlus(:,i)
        uPlusTemp(1:ng(3)) = uPlus(:,i)
        
        thickTemp(1) = delta99(i)
        
        call linear_interp(yPlusTemp, uPlusTemp, ng(3)+1, thickTemp, velTemp, 1)
        
        yPlusTemp(ng(3)+1) = delta99(i)
        uPlusTemp(ng(3)+1) = velTemp(1)

        ! Fix 8: Add error checking for sorting
        call sort_couple(yPlusTemp, uPlusTemp)

        ! Fix 9: Reset idLimit for each iteration
        idLimit = 0
        do j = 1, ng(3)+1
            if (yPlusTemp(j) <= delta99(i)) then
                idLimit = idLimit + 1
            end if
        end do
        allocate(integrand(idLimit))

        ! Calculate displacement and momentum thickness
        integrand = 1.0_rp - (uPlusTemp(1:idLimit)/(uInf(i)/uTau(i)))
        call trapezoidal_integral(yPlusTemp(1:idLimit), integrand, idLimit, reDelta(i))
        
        integrand = (1.0_rp - (uPlusTemp(1:idLimit)/(uInf(i)/uTau(i)))) * &
                   (uPlusTemp(1:idLimit)/(uInf(i)/uTau(i)))
        call trapezoidal_integral(yPlusTemp(1:idLimit), integrand, idLimit, reTheta(i))
        deallocate(integrand)
    end do
    if(myid == 0) then
        open(newunit=iunit,file=fname1d)
          do i=1,ng(1)
            write(iunit,fmt_rp) wss(i), reDelta(i), reTheta(i)
          end do
        close(iunit)

        fname2d = trim(datadir)//'yPlus'//fldnum//'.out'
        open(newunit=iunit,file=fname2d)
          do i=1,ng(3)
            write(iunit,fmt_rp) yPlus(i,:)
          end do
        close(iunit)
        fname2d = trim(datadir)//'uPlus'//fldnum//'.out'
        open(newunit=iunit,file=fname2d)
          do i=1,ng(3)
            write(iunit,fmt_rp) uPlus(i,:)
          end do
        close(iunit)
    end if
    ! Fix 11: Deallocate arrays before exiting
    !deallocate(delta99, wss, uTau, lscale, delta99plus)
    !deallocate(yPlus, uPlus, uMean, uInf)
    !deallocate(yPlusTemp, uPlusTemp)
    !deallocate(reDelta, reTheta)

>>>>>>> abb693ff53f5b8359f88929d9fd23b8557292714
end subroutine bl_stats
end module mod_stats