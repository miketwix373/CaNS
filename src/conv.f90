module mod_conv
  use mod_types
  use mod_utils
  implicit none
  private 
  public outflow

  contains
 subroutine outflow(n,ng,lo,hi,dt,l,dl,u,outf)
    real(rp), intent(in):: dt,dl(3),l(3)
    integer, intent(in):: n(3),ng(3),lo(3),hi(3)
    real(rp), intent(in) :: u(0:,0:,0:)
    real(rp), intent(inout) :: outf(0:,0:)
    real(rp), allocatable, dimension(:,:):: uMean 
    integer:: i,j
    real(rp):: c

    allocate(Umean(0:n(2)+1,0:n(3)+1))
    uMean = 0.0_rp

    if (hi(1)==ng(1)) then
        if (n(1)<= 20) then
            do i=1,n(1)
                uMean(:,:) = uMean(:,:) + u(i,:,:) 
            end do
            uMean = uMean/n(1)
        else
            do i=1,20
                uMean(:,:) = uMean(:,:) + u(i,:,:) 
            end do
            uMean = uMean/n(1)
        end if
    end if

    do i=0,n(2)+1
        do j=0,n(3)+1
            c = uMean(i,j)*dt/dl(1)
            outf(i,j)= u(n(1),i,j)*(1-c) + u(n(1)-1,i,j)*c
        end do
    end do






 end subroutine outflow

end module mod_conv