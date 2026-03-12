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

module exit_codes
  use kinds

  implicit none

  private

  integer(excode), parameter, public :: SUCCESS = 0_excode
  integer(excode), parameter, public :: FILE_ERROR = 10_excode
  integer(excode), parameter, public :: USER_INPUT_ERROR = 20_excode
  integer(excode), parameter, public :: USER_INPUT_WARNING = 21_excode
  integer(excode), parameter, public :: GIVEN_PARAMETER_ERROR = 30_excode
  integer(excode), parameter, public :: DERIVED_PARAMETER_ERROR = 31_excode
  integer(excode), parameter, public :: ALLOCATION_ERROR = 40_excode
  integer(excode), parameter, public :: DEALLOCATION_WARNING = 41_excode
  integer(excode), parameter, public :: NOT_IMPLEMENTED_ERROR = 50_excode

end module exit_codes
