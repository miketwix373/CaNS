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

subroutine get_Umean(ng, lo, hi, l, dl, n, p, uMeanOut)
  implicit none
  ! Input/Output arguments
  integer, intent(in), dimension(3) :: n, ng, lo, hi ! Grid dimensions and bounds
  real(rp), intent(in), dimension(3) :: l, dl ! Domain length and grid spacing
  real(rp), intent(in), dimension(0:,0:,0:) :: p ! Input field
  real(rp), allocatable, dimension(:,:) :: uMean ! Output mean (2D array)
  real(rp), intent(out), dimension(0:,0:) :: uMeanOut ! Output mean (2D array)

  ! Local variables
  integer :: np, rank, i, j, k, ierr
  real(rp), allocatable, dimension(:,:,:) :: uMean3D
  real(rp) :: dx, local_sum
  integer :: comm
  ! MPI initialization
  comm = MPI_COMM_WORLD
  call MPI_Comm_rank(comm, rank, ierr)
  call MPI_Comm_size(comm, np, ierr)
  ! Allocate arrays
  allocate(uMean3D(0:ng(2)+1,0:ng(3)+1,0:np-1))
  allocate(uMean(0:ng(2)+1,0:ng(3)+1))
  ! Calculate x-direction element size
  dx = dl(1)
  ! Compute sum over local domain along x-direction
  !$acc parallel loop collapse(2) private(local_sum) default(present)
  do k = 0, n(3) +1
    do j = 1, n(2) +1
      local_sum = 0._rp
      !$acc loop reduction(+:local_sum)
      do i = 0, n(1)+1
        local_sum = local_sum + p(i,j,k)
      end do
      uMean3D(j+lo(2)-1, k+lo(3)-1, rank) = local_sum * dx
    end do
  end do
  ! MPI reduction to get global sum
  call MPI_ALLREDUCE(MPI_IN_PLACE, uMean3D, size(uMean3D), MPI_REAL_RP, &
                    MPI_SUM, comm, ierr)
  ! Compute final average
    uMean = sum(uMean3D, dim=3) / l(1)
  deallocate(uMean3D)
  do k = 0, n(3)+1
    do j = 0, n(2)+1
      uMeanOut(j,k) = uMean(j + lo(2)-1,k+ lo(3)-1)
    end do
  end do
end subroutine get_Umean

end module mod_adv
