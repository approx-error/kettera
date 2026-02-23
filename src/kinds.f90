! This file is part of kettera - a numerical solver and visualizer for 
! the time-dependent Schrödinger equation
! Copyright (C) 2026 Juuso Kaarela
! 
! Kettera is free software: you can redistribute it and/or modify it under the terms of
! the GNU General Public License as published by the Free Software Foundation, either
! version 3 of the License, or (at your option) any later version.
!
! Kettera is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
! without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with kettera.
! If not, see <https://www.gnu.org/licenses/>.

module kinds
  implicit none

  private
  ! Kind parameters for integers used as exit codes and parameter labels
  integer, parameter, public :: excode = selected_int_kind(2) ! Corresponds to an 8-bit int
  integer, parameter, public :: label = selected_int_kind(4) ! Corresponds to a 16-bit int
  ! Kind parameters for 32, 64 and 128 bit integers
  integer, parameter, public :: i32 = selected_int_kind(9)
  integer, parameter, public :: i64 = selected_int_kind(18)
  integer, parameter, public :: i128 = selected_int_kind(38)
  ! Kind parameters for single precision floats (32 bits), double precision floats (64 bits),
  ! extended precision floats (80 bits) and quadruple precision floats (128 bits)
  integer, parameter, public :: r32 = selected_real_kind(P=6, R=37)
  integer, parameter, public :: r64 = selected_real_kind(P=15, R=307)
  integer, parameter, public :: r80 = selected_real_kind(P=18, R=4931)
  integer, parameter, public :: r128 = selected_real_kind(P=33, R=4931)
end module kinds
