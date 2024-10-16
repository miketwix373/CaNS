! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_utils
  implicit none
  private
  public bulk_mean,f_sizeof,swap,allocate_bc_vel,store_field, store_snap
  !@acc public device_memory_footprint
contains
  subroutine allocate_bc_vel(mat,bcvel)
    use mod_types
    implicit none
    real(rp), intent(inout), dimension (0:,0:) :: mat
    integer n1,n2,i,j
    real(rp) :: bcvel

    n1 = size(mat,1)
    n2 = size(mat,2)
    do i = 0, n1-1
        do j = 0, n2-1
            mat(i,j) = bcvel
        end do
    end do

  end subroutine allocate_bc_vel

 subroutine store_field(fname,nskip,p)
  use mod_common_mpi, only: ipencil => ipencil_axis
  use mpi
  use mod_types
  implicit none
  character(len=*), intent(in) :: fname
  integer, intent(in), dimension(3) :: nskip
  real(rp), intent(in), dimension(:,:,:) :: p
  integer :: fh, ierr, i, j, k
  integer(kind=MPI_OFFSET_KIND) :: disp
  character(len=30) :: fmt_str
  character(len=20) :: num_str
  character(len=:), allocatable :: line_buffer

  ! Open the file
  call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
       MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)

  ! Set the initial displacement to 0
  disp = 0_MPI_OFFSET_KIND

  ! Allocate buffer for a line of output
  allocate(character(len=20*size(p,1)) :: line_buffer)

  ! Write data in decimal format
  do k = 1, size(p,3), nskip(3)
    do j = 1, size(p,2), nskip(2)
      line_buffer = ''
      do i = 1, size(p,1), nskip(1)
        write(num_str, '(F20.6)') p(i,j,k)  ! Convert to decimal with 6 decimal places
        line_buffer = line_buffer // trim(adjustl(num_str)) // ' '
      end do
      line_buffer = trim(line_buffer) // new_line('a')
      
      ! Write the line
      call MPI_FILE_WRITE_AT_ALL(fh, disp, line_buffer, len_trim(line_buffer), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)
      
      ! Update displacement
      disp = disp + len_trim(line_buffer)
    end do
  end do

  ! Deallocate buffer
  deallocate(line_buffer)

  ! Close the file
  call MPI_FILE_CLOSE(fh, ierr)
end subroutine store_field

subroutine store_snap(fname, nskip, p)
    use mod_common_mpi, only: ipencil => ipencil_axis
    use mpi
    use mod_types
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(in), dimension(2) :: nskip
    real(rp), intent(in), dimension(:,:) :: p
    integer :: fh, ierr, i, j
    integer(kind=MPI_OFFSET_KIND) :: disp
    character(len=30) :: fmt_str
    character(len=20) :: num_str
    character(len=:), allocatable :: line_buffer

    ! Open the file
    call MPI_FILE_OPEN(MPI_COMM_WORLD, fname, &
        MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, fh, ierr)

    ! Set the initial displacement to 0
    disp = 0_MPI_OFFSET_KIND

    ! Allocate buffer for a line of output
    allocate(character(len=20*size(p,1)) :: line_buffer)

    ! Write data in decimal format
    do j = 1, size(p,2), nskip(2)
        line_buffer = ''
        do i = 1, size(p,1), nskip(1)
            write(num_str, '(F20.6)') p(i,j) ! Convert to decimal with 6 decimal places
            line_buffer = line_buffer // trim(adjustl(num_str)) // ' '
        end do
        line_buffer = trim(line_buffer) // new_line('a')

        ! Write the line
        call MPI_FILE_WRITE_AT_ALL(fh, disp, line_buffer, len_trim(line_buffer), MPI_CHARACTER, MPI_STATUS_IGNORE, ierr)

        ! Update displacement
        disp = disp + len_trim(line_buffer)
    end do

    ! Deallocate buffer
    deallocate(line_buffer)

    ! Close the file
    call MPI_FILE_CLOSE(fh, ierr)
end subroutine store_snap

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
