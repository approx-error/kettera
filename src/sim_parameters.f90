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

module sim_parameters
  ! Module contents: Parameters used in the simulations and derived types used for storing
  ! the parameters and arrays of data
  use kinds
  use io_parameters, only: MAX_STR_LEN

  implicit none

  ! "Enumerations" for the different solving methods, types of potential function and
  ! different types of wavefunction

  ! NOTE: Enumeration is in quotes since Fortran does not have an actual enumeration data type
  ! that has unmutable members. The closest we have are derived types with default value fields
  ! but this is not the same as a proper enumeration type for two reasons:
  !   1. Proper enumerations are immutable but derived type fields cannot be immutable
  !      (ie. they cannot be fortran parameters)
  !   2. Proper enumerations are their own distinct type meaning that if we had an enumeration
  !      type 'ExitCode', a function could require that the argument it is given is of type
  !      'ExitCode' even if the exit codes defined in 'ExitCode' are just simple integers.
  !      Fortran doesn't support distinguishing between an integer stored in a derived type
  !      and a normal integer so a function that requires a certain type of integer would still
  !      accept any integer that is of the same integer kind as the one in the derived type field
  !      definition. 

  type IterationMethod
    integer(label) :: INVALID = -1_label
    integer(label) :: CRANK_NICOLSON = 0_label
    integer(label) :: SPLIT_STEP = 1_label
  end type IterationMethod

  type WavefunctionType
    integer(label) :: INVALID = -1_label
    integer(label) :: GAUSSIAN = 0_label
    integer(label) :: SINC = 1_label
  end type WavefunctionType

  type PotentialType   
    integer(label) :: INVALID = -1_label
    integer(label) :: ZERO = 0_label
    integer(label) :: BOX = 1_label
    integer(label) :: BARRIER = 2_label
    integer(label) :: WELL = 3_label
    integer(label) :: DOUBLE_WELL = 4_label
    integer(label) :: HARMONIC = 5_label
    integer(label) :: LINEAR = 6_label
    integer(label) :: HALF_LINEAR = 7_label
    integer(label) :: ABSOLUTE_VALUE = 8_label
    integer(label) :: LOGARITHMIC = 9_label
    integer(label) :: HYPERBOLIC_COSINE = 10_label
  end type PotentialType

  ! Defining instances of the PotentialType and WavefunctionType derived types.
  ! These will be the only instances ever used and will simply hold the different
  ! potential types and wavefcuntion types for use by other pieces of code
  type(IterationMethod), public :: Method
  type(WavefunctionType), public :: WaveType
  type(PotentialType), public :: PotType

  ! Derived type for storing the parameters fed into the simulation
  type, public :: SimulationParams
    character(len=MAX_STR_LEN) :: param_file
    character(len=MAX_STR_LEN) :: input_file
    character(len=MAX_STR_LEN) :: output_file
    character(len=MAX_STR_LEN) :: log_file
    integer(label) :: iter_method
    logical :: imag_time
    logical :: unit_bounds
    integer(i32) :: step_count
    integer(i32) :: write_interval
    real(r64) :: delta_t
    integer(i32) :: point_count
    real(r64) :: x_max
    real(r64) :: delta_x
    real(r64) :: delta_p
    integer(label) :: wave_type
    real(r64) :: mass
    real(r64) :: charge
    real(r64) :: momentum
    real(r64) :: wave_offset
    real(r64) :: wave_width
    integer(label) :: pot_type
    real(r64) :: pot_offset
    real(r64) :: pot_width
    real(r64) :: pot_strength
  end type SimulationParams

  type(SimulationParams), public :: SimParams

  ! Derived type for storing all the relevant arrays used in Crank - Nicolson method iteration
  type, public :: CrankNicolsonArrays
    real(r64), allocatable :: x_space(:)
    complex(r64), allocatable :: wavefunction(:)
    real(r64), allocatable :: potential(:)
    real(r64), allocatable :: hamiltonian(:,:)
    complex(r64), allocatable :: forward_op(:,:)
    complex(r64), allocatable :: backward_op(:,:)
  end type CrankNicolsonArrays

  type(CrankNicolsonArrays), public :: CNArrays

  ! Derived type for storing all the relevant arrays used in split step method iteration
  type, public :: SplitStepArrays
    real(r64), allocatable :: x_space(:)
    real(r64), allocatable :: p_space(:)
    complex(r64), allocatable :: wavefunction(:)
    real(r64), allocatable :: potential(:)
    real(r64), allocatable :: kinetic(:)
    real(r64), allocatable :: hamiltonian(:,:)
    complex(r64), allocatable :: x_space_op(:)
    complex(r64), allocatable :: p_space_op(:)
  end type SplitStepArrays

  type(SplitStepArrays), public :: SSArrays

  ! Derived type for storing all the arrays that contain the output data produced by the simulation
  type, public :: OutputArrays
    real(r64), allocatable :: real_part(:)
    real(r64), allocatable :: imag_part(:)
    real(r64), allocatable :: amplitude(:)
    real(r64), allocatable :: density(:)
  end type OutputArrays

  type(OutputArrays), public :: OutArrays

  ! Derived type for storing all the arrays that contain the logging data produced by the simulation
  type, public :: LoggingArrays
    integer(i32), allocatable :: timestamps(:)
    real(r64), allocatable :: energies(:)
    real(r64), allocatable :: squared_norms(:)
  end type LoggingArrays

  type(LoggingArrays), public :: LogArrays

  ! Numerical parameters
  real(r64), parameter, public :: PI_CONST = 4.0_r64 * atan2(1.0_r64, 1.0_r64)
  real(r64), parameter, public :: HBAR_CONST = 1.0545718176461565e-34_r64 ! Reduced Planck constant h-bar
  complex(r64), parameter, public :: IMAG_UNIT = (0.0_r64, 1.0_r64) ! imaginary unit i 
  real(r64), parameter, public :: EPS = 1e-5_r64 ! Used to check if a real variable is close to zero

  ! Limits for parameters
  real(r64), parameter, public :: DELTA_X_MAX = 1.0_r64
  !TODO: figure out a good DELTA_X_MAX: real(r64), parameter, public :: DELTA_X_MAX = 0.02_r64
  real(r64), parameter, public :: X_BOUND_MIN = 0.1_r64
  real(r64), parameter, public :: MASS_MIN = 0.001_r64
  real(r64), parameter, public :: TIME_RANGE_MIN = 10.0_r64

end module sim_parameters
