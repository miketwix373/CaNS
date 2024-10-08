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
  subroutine adv(n,dt,dli,u,u_old,uBC,uMean)
  integer, intent(in), dimension(3) :: n
  real, intent(in),dimension(0:,0:, 0:)::  u,u_old,dli
  real, intent(in), dimension(0:,0:):: uMean
  real, allocatable :: uBC(:,:)
  integer :: i, j, k
  real, intent(in) :: dt
  real:: c
  ! Initialize upast with fixed values
   ! Print the matrices
  allocate(uBC(n(2),n(3)))
     do i =1,n(2)
        do j=1,n(3)
                c             = uMean(i,j)*dt/dli(n(1),i,j)
                uBC(i,j)        = u_old(1,i,j)/c + (c-1)/c*u(n(1)-1,i,j)
        end do
        
     end do
  end subroutine adv


  subroutine get_Umean(ng,lo,hi,l,dl,p,uMean)
    implicit none
    integer , intent(in), dimension(3) :: ng,lo,hi
    integer , intent(in) :: idir
    real(rp), intent(in), dimension(3) :: l,dl
    real(rp), intent(in), dimension(lo(1)-1:,lo(2)-1:,lo(3)-1:) :: p
    real(rp), allocatable, dimension(:,:) :: p1d
    character(len=*), intent(in) :: fileName 
    integer :: i,j,k
    integer :: iunit
    real(rp) :: grid_length_ratio,p1d_s
    allocate(p1d(ng(2),ng(3)))
    !$acc enter data create(p1d)
    !$acc kernels default(present)
    p1d(:,:) = 0._rp
    !$acc end kernels
    grid_length_ratio = dl(1)/(l(1))
    !$acc parallel loop gang default(present) private(p1d_s)
    do k=lo(3),hi(3)
      do j=lo(2),hi(2)
      p1d_s = 0._rp
      !$acc loop collapse(2) reduction(+:p1d_s)
        do i=lo(1),hi(1)
          p1d_s = p1d_s + p(i,j,k)*grid_length_ratio
        end do
        p1d(j,k) = p1d_s
      end do
    end do
    !$acc exit data copyout(p1d)
    call MPI_ALLREDUCE(MPI_IN_PLACE,p1d(1,1),ng(2)*ng(3),MPI_REAL_RP,MPI_SUM,MPI_COMM_WORLD,ierr)
  end subroutine get_Umean
end module mod_adv
