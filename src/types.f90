! -
!
! SPDX-FileCopyrightText: Copyright (c) 2017-2022 Pedro Costa and the CaNS contributors. All rights reserved.
! SPDX-License-Identifier: MIT
!
! -
module mod_types
  use mpi, only: MPI_REAL,MPI_DOUBLE_PRECISION
  integer, parameter, public :: sp = selected_real_kind(6 , 37), &
                                dp = selected_real_kind(15,307), &
                                i8 = selected_int_kind(18)
#if defined(_SINGLE_PRECISION)
  integer, parameter, public :: rp = sp
  integer, parameter, public :: MPI_REAL_RP = MPI_REAL
#else
  integer, parameter, public :: rp = dp
  integer, parameter, public :: MPI_REAL_RP = MPI_DOUBLE_PRECISION
#endif
  type  flow_data
    real(rp), allocatable :: inf(:,:)
    real(rp), allocatable :: outf(:,:)
  end type flow_data

  type xyz_case
    type(flow_data) :: x, y, z
  end type xyz_case

  type  bc_direct
    type(xyz_case) :: u, v, w
  end type bc_direct
end module mod_types
