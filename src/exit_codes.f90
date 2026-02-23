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

  integer(excode), parameter, public :: SUCCESS = 0
  integer(excode), parameter, public :: FILE_ERROR = 10
  integer(excode), parameter, public :: GIVEN_PARAMETER_ERROR = 20
  integer(excode), parameter, public :: DERIVED_PARAMETER_ERROR = 21
  integer(excode), parameter, public :: ALLOCATION_ERROR = 30
  integer(excode), parameter, public :: DEALLOCATION_WARNING = 31

end module exit_codes
