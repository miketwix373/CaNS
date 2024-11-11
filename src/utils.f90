! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_utils
  use mod_types
  use mpi
  use mod_common_mpi, only: ierr, myid
  implicit none
  private 
  public bulk_mean,f_sizeof,swap,advection,trapezoidal_integral,&
   identify_fringe,fringeForce,linear_interp,MeanFlow2D,sort_couple,check_init_profile


  !@acc public device_memory_footprint
contains
  subroutine check_init_profile (profile,n,zcg,fname,dir)
    character(len=*), intent(in) :: fname
    real(rp),intent(in)::profile(:,:,:), zcg(:)
    integer, intent(in):: n(3),dir
    integer:: i,iunit
    character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                 fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
  character(len=*), parameter :: fmt_rp = fmt_dp
#else
  character(len=*), parameter :: fmt_rp = fmt_sp
#endif

    if(myid == 0) then
        open(newunit=iunit,file=fname)
          do i=1,n(3)
            write(iunit,fmt_rp) zcg(i),profile(dir,1,i)
          end do
        close(iunit)
    end if
  end subroutine check_init_profile 

  subroutine sort_couple(vec_a, vec_b)
      real(rp), dimension(:), intent(inout) :: vec_a  ! Vector to sort by
      real(rp), dimension(:), intent(inout) :: vec_b  ! Related vector
      integer :: i, j, n
      real(rp) :: temp_a, temp_b
      
      ! Check if vectors have same size
      n = size(vec_a)
      
      ! Bubble sort implementation
      do i = 1, n-1
          do j = 1, n-i
              if (vec_a(j) > vec_a(j+1)) then
                  ! Swap elements in vec_a
                  temp_a = vec_a(j)
                  vec_a(j) = vec_a(j+1)
                  vec_a(j+1) = temp_a
                  
                  ! Swap corresponding elements in vec_b
                  temp_b = vec_b(j)
                  vec_b(j) = vec_b(j+1)
                  vec_b(j+1) = temp_b
              end if
          end do
      end do
  end subroutine sort_couple

  subroutine linear_interp(x, y, n, x_new, y_new, n_new)
      integer, intent(in) :: n, n_new
      real(dp), intent(in) :: x(n), y(n), x_new(n_new)
      real(dp), intent(out) :: y_new(n_new)
      
      ! Local variables
      integer :: i, j
      real(dp) :: t
      logical :: found
      
      
      ! Perform interpolation
      do i = 1, n_new
          ! Check bounds
          if (x_new(i) < x(1) .or. x_new(i) > x(n)) then
              return
          endif
          
          ! Find interval
          found = .false.
          do j = 1, n-1
              if (x_new(i) >= x(j) .and. x_new(i) <= x(j+1)) then
                  ! Linear interpolation formula
                  t = (x_new(i) - x(j)) / (x(j+1) - x(j))
                  y_new(i) = y(j) + t * (y(j+1) - y(j))
                  found = .true.
                  exit
              endif
          end do
          
          ! Handle exact match with last point
          if (.not. found .and. abs(x_new(i) - x(n)) < tiny(1.0_dp)) then
              y_new(i) = y(n)
          endif
      end do
      
  end subroutine linear_interp

  subroutine cosine_blend_weight(x, N, weight)      
      ! Input/Output variables
      real(rp), intent(in)  :: x  ! Input value
      real(rp), intent(in)  :: N    ! Lower bound (0 < N < 1)
      real(rp), intent(out) :: weight ! Output weight
      real(rp) :: pi, x_normalized
      
      pi = ACOS(-1.0)
      ! Normalize x from [N,1] to [0,1]
      x_normalized = (x - N)/(1.0_rp - N)
      
      ! Compute cosine blend weight: (1 - cos(π*x))/2
      weight = (1.0_rp - cos(x_normalized * pi))/2.0_rp
      
  end subroutine cosine_blend_weight

  subroutine fringeForce (bforce,isFringe,dt,u,utarget,lo,fringeLim,L,dir)
    real(rp), intent(inout), dimension(0:,0:,0:):: bforce,u
    logical, intent(in), dimension(0:,0:,0:):: isFringe
    real(rp), intent(in) :: utarget(:,0:,0:), dt
    integer :: i,j,k, lo(3),fringeLim,n,L,dir
    real(rp):: weight,x,fringeStart
    n = size(isFringe,1)
    do i = 0, (n-1)
      if (isFringe(i,1,1)) then
        x = (lo(1)-1+ i)/L
        fringeStart = fringeLim/L
        call cosine_blend_weight(x,fringeStart,weight)
        bforce(i,:,:) = -weight*(u(i,:,:)-utarget(dir,:,:))/dt
      end if
    end do
  end subroutine fringeForce

  subroutine  identify_fringe(isFringe,loLimFringe,lo)
    logical, dimension(:,:,:), intent(inout) :: isFringe
    integer, intent(in) :: loLimFringe, lo(3)
    integer :: n,i

    n = size(isFringe,1)
    do i = 0, n-1   
      if (lo(1)-1+i > loLimFringe) then
        isFringe(i,:,:) = .true.
      end if
    end do 
  end subroutine identify_fringe

  subroutine trapezoidal_integral(x, fx, n, integral) 
      integer, intent(in) :: n
      real(rp), dimension(:), intent(in) :: x, fx
      real(rp) :: integral
      integer :: i
      
      ! Initialize integral with first and last points (half contribution)
      integral = 0.5_rp * (fx(1) + fx(n))
      
      ! Add contribution from middle points
      do i = 2, n-1
          integral = integral + fx(i)
      end do
      
      ! Multiply by dx (assumes uniform spacing)
      integral = integral * (x(n) - x(1)) / (n - 1)
      
  end subroutine trapezoidal_integral

subroutine MeanFlow2D(ng, n, lo, hi, dl, l, p, pMean, idir, localFlag)
    integer, intent(in) :: ng(3), n(3), lo(3), hi(3), idir
    integer :: coord(3), i, j, k, rank, np, startIter(2), endIter(2)
    integer :: nrows, ncols, ngrows, ngcols, nmean, ngmean
    integer :: blockSelect(2), row, nblock(2), nblocksCol, nblocksRow
    integer, dimension(:,:), allocatable :: coord3, n3d
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), intent(inout), dimension(0:,0:) :: pMean
    real(rp), intent(in), dimension(3) :: dl, l
    real(rp), dimension(:,:), allocatable :: pMeanGlob, pMeanLoc
    real(rp), dimension(:,:,:), allocatable :: pMeanLoc3d
    real(rp) :: localsum
    logical, intent(in) :: localFlag
    integer, dimension(:), allocatable :: maxn
    
    ! Get MPI info
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, np, ierr)

    ! Set dimensions based on direction
    select case (idir) 
    case(1)
        ncols = n(2)
        nrows = n(3)
        ngcols = ng(2)
        ngrows = ng(3)
        nmean = n(1)
        ngmean = ng(1)
        blockSelect = (/3,2/)
    case(2)
        ncols = n(1)
        nrows = n(3)
        ngcols = ng(1)
        ngrows = ng(3)
        nmean = n(2)
        ngmean = ng(2)
        blockSelect = (/3,1/)
    case(3)
        ncols = n(2)
        nrows = n(1)
        ngcols = ng(2)
        ngrows = ng(1)
        nmean = n(3)
        ngmean = ng(3)
        blockSelect = (/1,2/)
    end select

    ! Allocate arrays
    allocate(pMeanLoc(0:nrows+1, 0:ncols+1))
    allocate(pMeanGlob(0:ngrows+1, 0:ngcols+1))
    allocate(pMeanLoc3d(0:nrows+1, 0:ncols+1, np))
    allocate(coord3(3, np))
    allocate(n3d(3, np))
    allocate(maxn(np))

    ! Initialize arrays
    pMeanGlob = 0.0_rp
    pMeanLoc = 0.0_rp
    pMeanLoc3d = 0.0_rp
    coord = 0
    coord3 = 0

    ! Compute block ID positions
    coord = hi / n

    ! Compute mean in stencil face
    do i = 0, nrows+1
        do j = 0, ncols+1
            localsum = 0.0_rp
            do k = 1, nmean
                select case (idir)
                case (1)
                    localsum = localsum + p(k,j,i) * dl(1) / l(1)
                case (2)
                    localsum = localsum + p(j,k,i) * dl(2) / l(2)
                case (3)
                    localsum = localsum + p(i,j,k) * dl(3) / l(3)
                end select
            end do
            pMeanLoc(i,j) = localsum
        end do
    end do

    ! Gather data across processes
    call MPI_Allgather(pMeanLoc, (nrows+2)*(ncols+2), MPI_REAL_RP, &
                       pMeanLoc3d, (nrows+2)*(ncols+2), MPI_REAL_RP, &
                       MPI_COMM_WORLD, ierr)
    call MPI_Allgather(coord, 3, MPI_INTEGER, coord3, 3, MPI_INTEGER, &
                       MPI_COMM_WORLD, ierr)
    call MPI_Allgather(n, 3, MPI_INTEGER, n3d, 3, MPI_INTEGER, &
                       MPI_COMM_WORLD, ierr)

    ! Find max block dimensions
    maxn = coord3(blockSelect(1),:)
    nblocksRow = maxval(maxn)
    maxn = coord3(blockSelect(2),:)
    nblocksCol = maxval(maxn)
    nblock = (/nblocksRow, nblocksCol/)

    ! Assemble global mean
    do k = 1, np
        startIter = 1
        endIter = n3d(blockSelect, k)

        ! Adjust iteration bounds for edge blocks
        do row = 1, 2
            if (coord3(blockSelect(row), k) == 1) then
                startIter(row) = 0
            end if
            if (coord3(blockSelect(row), k) == nblock(row)) then
                endIter(row) = endIter(row) + 1
            end if
        end do

        ! Assemble global array
        do i = startIter(1), endIter(1)
            do j = startIter(2), endIter(2)
                select case (idir)
                case(1)
                    pMeanGlob((coord3(3,k)-1)*n3d(3,k)+i, (coord3(2,k)-1)*n3d(2,k)+j) = &
                        pMeanGlob((coord3(3,k)-1)*n3d(3,k)+i, (coord3(2,k)-1)*n3d(2,k)+j) + &
                        pMeanLoc3d(i,j,k)
                case(2)
                    pMeanGlob((coord3(3,k)-1)*n3d(3,k)+i, (coord3(1,k)-1)*n3d(1,k)+j) = &
                        pMeanGlob((coord3(3,k)-1)*n3d(3,k)+i, (coord3(1,k)-1)*n3d(1,k)+j) + &
                        pMeanLoc3d(i,j,k)
                case(3)
                    pMeanGlob((coord3(1,k)-1)*n3d(1,k)+i, (coord3(2,k)-1)*n3d(2,k)+j) = &
                        pMeanGlob((coord3(1,k)-1)*n3d(1,k)+i, (coord3(2,k)-1)*n3d(2,k)+j) + &
                        pMeanLoc3d(i,j,k)
                end select
            end do
        end do
    end do

    ! Set output based on localFlag
    if (localFlag) then
        pMean = pMeanGlob(lo(blockSelect(1))-1:hi(blockSelect(1))+1, &
                         lo(blockSelect(2))-1:hi(blockSelect(2))+1)
    else
        pMean = pMeanGlob
    end if

    ! Deallocate arrays
    deallocate(pMeanLoc, pMeanGlob, pMeanLoc3d, coord3, n3d, maxn)

end subroutine MeanFlow2D

  subroutine advection (n,dt,dl,upast,uMean,u_adv)
    real(rp), intent(in), dimension(:,0:,0:) :: upast
    real(rp), intent(inout), dimension(0:,0:) :: u_adv
    real(rp), intent(in), dimension(0:,0:) :: uMean
    real(rp), intent(in) :: dt, dl
    integer,intent(in) :: n(3)
    real(rp):: c
    integer:: j,k

    do j=0,n(2)+1
      do k=0,n(3)+1
        !c = uMean(j,k) *dt/dl
        c = dt/dl
        u_adv(j,k) = upast(2,j,k)*(1-c)+c*upast(1,j,k)
      end do
    end do


  end subroutine advection

  subroutine bulk_mean(n,grid_vol_ratio,p,mean)
    !
    ! compute the mean value of an observable over the entire domain
    !
    use mpi
    use mod_types
    implicit none
    integer , intent(in), dimension(3) :: n
    real(rp), intent(in), dimension(0:) :: grid_vol_ratio
    real(rp), intent(in), dimension(0:,0:,0:) :: p
    real(rp), intent(out) :: mean
    integer :: i,j,k
    integer :: ierr
    mean = 0.
    !$acc data copy(mean) async(1)
    !$acc parallel loop collapse(3) default(present) reduction(+:mean) async(1)
    !$OMP PARALLEL DO   COLLAPSE(3) DEFAULT(shared)  REDUCTION(+:mean)
    do k=1,n(3)
      do j=1,n(2)
        do i=1,n(1)
          mean = mean + p(i,j,k)*grid_vol_ratio(k)
        end do
      end do
    end do
    !$acc end data
    !$acc wait(1)
    call MPI_ALLREDUCE(MPI_IN_PLACE,mean,1,MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine bulk_mean
  pure integer function f_sizeof(val) result(isize)
    !
    ! returns storage size of the scalar argument val in bytes
    !
    implicit none
    class(*), intent(in) :: val
    isize = storage_size(val)/8
  end function f_sizeof
  subroutine swap(arr1,arr2)
    use mod_types, only: rp
    implicit none
    real(rp), intent(inout), pointer, contiguous, dimension(:,:,:) :: arr1,arr2
    real(rp),                pointer, contiguous, dimension(:,:,:) :: tmp
    tmp  => arr1
    arr1 => arr2
    arr2 => tmp
  end subroutine swap
#if defined(_OPENACC)
  function device_memory_footprint(n,n_z) result(itotal)
    !
    ! estimate GPU memory footprint, assuming one MPI task <-> one GPU
    !
    use mod_types, only: i8,rp
    integer, intent(in), dimension(3) :: n,n_z
    integer :: nh(3)
    integer(i8) :: itotal,itemp,rp_size
    rp_size = f_sizeof(1._rp)
    itotal = 0
    !
    ! 1. 'main' arrays: u,v,w,p,pp
    !
    nh(:) = 1
    itotal = itotal + product(n(:)+2*nh(:))*rp_size*5
    !
    ! 2. grids arrays: zc,zf,dzc,dzf,dzci,dzfi,grid_vol_ratio_c,grid_vol_ratio_f (tiny footprint)
    !
    nh(:) = 1
    itotal = itotal + (n(3)+2*nh(3))*rp_size*8
    !
    ! 3. solver eigenvalues and Gauss elimination coefficient arrays (small footprint)
    !    rhs?%[x,y,z] arrays, lambdaxy? arrays, and a?,b?,c? arrays
    !
    block
      integer(i8) :: itemp1,itemp1_(3),itemp2,itemp3
      itemp1_(:) = [n_z(2)*n_z(3)*2,n_z(1)*n_z(3)*2,n_z(1)*n_z(2)*2]
      itemp1 = sum(itemp1_(:))   ! rhs
      itemp2 = product(n_z(1:2)) ! lambdaxy
      itemp3 = n_z(3)*3          ! a,b,c
#if   !defined(_IMPDIFF)
      !
      ! rhsbp, lambdaxyp, ap,bp,cp
      !
      itotal = itotal + itemp1*rp_size                + itemp2*rp_size   + itemp3*rp_size
#elif  defined(_IMPDIFF_1D)
      !
      ! rhsbp,rhsb[u,v,w,buf]%z, lambdaxyp, a?,b?,c? [p,u,v,w,buf]
      !
      itotal = itotal + (itemp1+itemp1_(3)*4)*rp_size + itemp2*rp_size   + itemp3*rp_size*5
#else
      !
      ! rhsbp,rhsb[u,v,w,buf]%[x,y,z], lambdaxy[p,u,v,w], (a?,b?,c?)[p,u,v,w,buf]
      !
      itotal = itotal + itemp1*rp_size*(1+4)          + itemp2*rp_size*5 + itemp3*rp_size*5
#endif
    end block
    !
    ! 4. prediction velocity arrays arrays d[u,v,w]dtrk_t, d[u,v,w]dtrko_t
    !
    itemp  = product(n(:))*rp_size
    itotal = itotal + itemp*6
#if defined(_IMPDIFF)
    itotal = itotal + itemp*3
#endif
    !
    ! 5. transpose & FFT buffer arrays, halo buffer arrays, and solver arrays
    !    taken directly from `mod_common_cudecomp`
    !
    block
      use mod_common_cudecomp, only: work,work_halo,solver_buf_0,solver_buf_1,pz_aux_1,pz_aux_2
      itemp = storage_size(work        ,i8)*size(work        ) + storage_size(work_halo   ,i8)*size(work_halo   ) + &
              storage_size(solver_buf_0,i8)*size(solver_buf_0) + storage_size(solver_buf_1,i8)*size(solver_buf_1) + &
              storage_size(pz_aux_1    ,i8)*size(pz_aux_1    ) + storage_size(pz_aux_2    ,i8)*size(pz_aux_2    )
      itotal = itotal + itemp/8
    end block
  end function device_memory_footprint
#endif
end module mod_utils
