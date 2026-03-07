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

program main

  use kinds
  use exit_codes
  use io_parameters
  use sim_parameters, only: SimParams, CNArrays, SSArrays, OutArrays, LogArrays 
  use user_io

  use alloc_dealloc ! Maybe remove later
  use initialize ! Maybe remove later
  use read_write ! Add later
  use iterate ! Maybe remove later 

  implicit none

  integer(label) :: output_mode
  integer(excode) :: exit_code

  integer(i32) :: N, frames

  call init_params(SimParams)

  call parse_user_input(SimParams, output_mode, exit_code)

  if (exit_code /= SUCCESS .or. exit_code /= USER_INPUT_WARNING) then
    stop int(exit_code, i32)
  end if

  stop 0

  call read_param_file(SimParams, exit_code)
  if (exit_code /= SUCCESS) then
    stop int(exit_code, i32)
  end if
  
  call alloc_cn_arrays(N, CNArrays, exit_code)
  call alloc_output_arrays(N, OutArrays, exit_code)
  call alloc_log_arrays(frames, LogArrays, exit_code)
  if (exit_code /= SUCCESS) then
    call dealloc_cn_arrays(CNArrays, exit_code)
    call dealloc_output_arrays(OutArrays, exit_code)
    call dealloc_log_arrays(LogArrays, exit_code)
    stop int(exit_code, i32)
  end if

  call read_input_file(SimParams, CNArrays%orthogonal, exit_code)
  if (exit_code /= SUCCESS) then
    call dealloc_cn_arrays(CNArrays, exit_code)
    call dealloc_output_arrays(OutArrays, exit_code)
    call dealloc_log_arrays(LogArrays, exit_code)
    stop int(exit_code, i32)
  end if


  call init_cn_iteration(SimParams, CNArrays, exit_code)
  if (exit_code /= SUCCESS) then
    call dealloc_cn_arrays(CNArrays, exit_code)
    call dealloc_output_arrays(OutArrays, exit_code)
    call dealloc_log_arrays(LogArrays, exit_code)
    print *, 'exit_code', exit_code
    stop int(exit_code, i32)
  end if

  call iter_cn_method(SimParams, CNArrays, OutArrays, LogArrays, exit_code)
  if (exit_code /= SUCCESS) then
    call dealloc_cn_arrays(CNArrays, exit_code)
    call dealloc_output_arrays(OutArrays, exit_code)
    call dealloc_log_arrays(LogArrays, exit_code)
    stop int(exit_code, i32)
  end if

  call dealloc_cn_arrays(CNArrays, exit_code)
  call dealloc_output_arrays(OutArrays, exit_code)
  call dealloc_log_arrays(LogArrays, exit_code)
  stop int(exit_code, i32)

end program main
