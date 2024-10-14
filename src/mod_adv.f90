module mod_adv
  use mpi
  use decomp_2d_io
  use mod_common_mpi, only:ierr,myid
  use mod_types
  implicit none 
    character(len=*), parameter :: fmt_dp = '(*(es24.16e3,1x))', &
                                 fmt_sp = '(*(es15.8e2,1x))'
#if !defined(_SINGLE_PRECISION)
  character(len=*), parameter :: fmt_rp = fmt_dp
#else
  character(len=*), parameter :: fmt_rp = fmt_sp
#endif

  contains
  subroutine adv(dt,dli,u,u_old,uBC,uMean)
  real(rp), intent(in), dimension(3) :: dli
  real(rp), intent(in),dimension(0:,0:, 0:)::  u,u_old
  real(rp), intent(in):: uMean
  real(rp) :: uBC(:,:)
  integer :: i, j, k, n1, n2, nx
  real(rp), intent(in) :: dt
  real(rp):: c
  ! Initialize upast with fixed values
   ! Print the matrices
  n1 = size(uBC,1)
  n2 = size(uBC,2)
  nx = size(u,1)
  do i =1,n1
    do j=1,n2
            c             = uMean*dt/dli(1)
            uBC(i,j)        = u_old(1,i,j)/c + (c-1)/c*u(nx-2,i,j)
    end do
  end do
  end subroutine adv


subroutine get_Umean(ng, lo, hi, l, dl, p, uMean)
    implicit none
    ! Input/Output arguments
    integer, intent(in), dimension(3) :: ng, lo, hi      ! Grid dimensions and bounds
    real(rp), intent(in), dimension(3) :: l, dl          ! Domain length and grid spacing
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p  ! Input field
    real(rp), intent(out) :: uMean                       ! Output mean (scalar)
    
    ! Local variables
    integer :: i, j, k
    integer :: ierr
    real(rp) :: local_sum, global_sum
    real(rp) :: volume_element, total_volume
    
    ! Calculate volume element and total volume
    volume_element = dl(1) * dl(2) * dl(3)
    total_volume = l(1) * l(2) * l(3)
    
    ! Initialize local sum
    local_sum = 0._rp
    
    ! Compute sum over local domain
    !$acc parallel loop collapse(3) reduction(+:local_sum) default(present)
    do k = lo(3), hi(3)
        do j = lo(2), hi(2)
            do i = lo(1), hi(1)
                local_sum = local_sum + p(i,j,k) * volume_element
            end do
        end do
    end do
    
    ! MPI reduction to get global sum
    call MPI_ALLREDUCE(local_sum, global_sum, 1, MPI_REAL_RP, &
                       MPI_SUM, MPI_COMM_WORLD, ierr)
    
    ! Compute final average
    uMean = global_sum / total_volume
    
end subroutine get_Umean
end module mod_adv
