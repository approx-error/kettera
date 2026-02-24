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
  use sim_parameters, only: Method, PotType, WaveType, SimParams, CNArrays, SSArrays, &
    OutArrays, LogArrays 
  !use io_parameters

  use alloc_dealloc ! Maybe remove later
  use initialize ! Maybe remove later
  !use read_write ! Add later
  use iterate ! Maybe remove later 

  implicit none

  ! TESTING

  
  integer(excode) :: exit_code

  integer(i32) :: N, frames

  SimParams%param_file = 'not used'
  SimParams%input_file = 'not used'
  SimParams%output_file = 'wave.out'
  SimParams%log_file = 'ket.out'

  SimParams%iter_method = Method%CRANK_NICOLSON
  SimParams%imag_time = .false.!.true.!.false.
  SimParams%unit_bounds = .true. 

  SimParams%step_count = 20000_i32!10000_i32!10000_i32
  SimParams%write_interval = 100_i32
  SimParams%delta_t = 1e-3_r64

  SimParams%point_count = 1024_i32!512_i32
  SimParams%x_max = 10.0_r64

  SimParams%wave_type = WaveType%GAUSSIAN
  SimParams%mass = 1.0_r64
  SimParams%charge = 1.0_r64
  SimParams%wave_offset = -1.0_r64
  SimParams%wave_width = 1.0_r64
  SimParams%wave_momentum = 0.0_r64

  SimParams%pot_type = PotType%HARMONIC
  SimParams%pot_offset = 0.0_r64
  SimParams%pot_width = 1.0_r64
  SimParams%pot_strength = 1.0_r64

  N = SimParams%point_count
  frames = SimParams%step_count / SimParams%write_interval + 1 
  
  call alloc_cn_arrays(N, CNArrays, exit_code)
  call alloc_output_arrays(N, OutArrays, exit_code)
  call alloc_log_arrays(frames, LogArrays, exit_code)

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
