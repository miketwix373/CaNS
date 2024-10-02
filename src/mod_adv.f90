module mod_adv
  implicit none 
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
                uBC(i,j)        = u_old(n(1)-1,i,j)/c + (c-1)/c*u(n(1)-1,i,j)
        end do
        
     end do
  end subroutine adv


  subroutine get_Umean(u,n,uMean)
    real, intent(in),dimension(0:,0:, 0:)::  u
    integer, intent(in), dimension(3):: n
    real, allocatable, intent(out):: uMean(:,:)
    integer :: i,j

    allocate(uMean(n(2),n(3)))
    do i=1,n(2)
      do j=1,n(3)
        uMean(i,j)      = sum(u(:,i,j), dim =1)/n(1)  
      end do
    end do
    
  end subroutine get_Umean
end module mod_adv
