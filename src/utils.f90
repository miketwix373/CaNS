! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_utils
  use mod_types
  implicit none
  private cosine_blend_weight
  public bulk_mean,f_sizeof,swap,advection,trapezoidal_integral, identify_fringe,fringeForce


  !@acc public device_memory_footprint
contains

contains
  subroutine cosine_blend_weight(x, N, weight)      
      ! Input/Output variables
      real(rp), intent(in)  :: x    ! Input value
      real(rp), intent(in)  :: N    ! Lower bound (0 < N < 1)
      real(rp), intent(out) :: weight ! Output weight
      
      ! Normalize x from [N,1] to [0,1]
      x_normalized = (x - N)/(1.0_rp - N)
      
      ! Compute cosine blend weight: (1 - cos(π*x))/2
      weight = (1.0_rp - cos(x_normalized * pi))/2.0_rp
      
  end subroutine cosine_blend_weight

  subroutine fringeForce (bforce,isFringe,dt,u,utarget,lo,fringeLim,L)
  real(rp), intent(inout), dimension(0:,0:,0:):: bforce,u
  logical, intent(in), dimension(0:,0:,0:):: isFringe
  real(rp), intent(in) :: utarget(0:,0:,1), dt
  integer :: i,j,k, lo(3),fringeLim,n(3),x,L
  real(rp):: weight
  n = size(bforce)
  do i = 0, n(1)
    if (isFringe(i,1,1)) then
      x = (lo(1)-1+ i)
      call cosine_blend_weight(x/L,fringeLim/L,weight)
      bforce(i,:,:) = weight*(u-utarget)/dt
    end if
  end do
  end subroutine fringeForce

  subroutine  identify_fringe(isFringe,loLimFringe,lo)
    logical, dimension(:,:,:), intent(inout) :: isFringe
    integer, intent(in) :: loLimFringe, lo(3)
    integer :: n(3),i

    n = size(isFringe,3)-2
    do i = 0, n(1)
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
