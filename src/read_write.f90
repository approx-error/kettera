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

module read_write
  ! Subroutines and functions for reading data from files and writing data to files or stdout
  use kinds
  use exit_codes
  use io_parameters
  use sim_parameters

  implicit none

  public 

  contains

    subroutine init_output_file(params, exit_code)
      implicit none

      type(SimulationParams), intent(in) :: params
      integer(excode), intent(out) :: exit_code

      character(len=MAX_STR_LEN) :: iomessage
      integer(i32) :: iostatus
      character(len=MAX_STR_LEN) :: filename
      integer(i32) :: framecount
      integer(label) :: method, pot, imag
      real(r64) :: x_bound


      filename = trim(params%output_file)
      framecount = params%step_count / params%write_interval + 1
      method = params%iter_method
      x_bound = params%x_max
      pot = params%pot_type
      if (params%imag_time) then
        imag = 1
      else
        imag = 0
      end if

      open(unit=OUTPUT_FILE_UNIT, file=filename, status='replace', action='write', iostat=iostatus, iomsg=iomessage)

      if (iostatus /= 0) then
        print '(a,a,a,/,a,a)', 'init_output_file: ERROR: Could not open file "', filename, '"', 'Message: ', trim(iomessage) 
        exit_code = FILE_ERROR
        return
      end if 

      write(OUTPUT_FILE_UNIT, fmt='(A,/,A,/,A)') trim(OUTPUT_TITLE_1), trim(OUTPUT_TITLE_2), &
        trim(OUTPUT_TITLE_3)
      write(OUTPUT_FILE_UNIT, fmt=trim(OUTPUT_HEADER_FMT)) framecount, method, x_bound, pot, imag
      write(OUTPUT_FILE_UNIT, fmt='(a)') trim(OUTPUT_DATA_FIELDS)

      close(OUTPUT_FILE_UNIT)

      exit_code = SUCCESS
      return
    end subroutine init_output_file

    subroutine write_output_file(params, frame, out_arrays, x_values, potential, exit_code)
      implicit none

      type(SimulationParams), intent(in) :: params
      integer(i32), intent(in) :: frame
      type(OutputArrays), intent(in) :: out_arrays
      real(r64), intent(in) :: x_values(:)
      real(r64), intent(in) :: potential(:)
      integer(excode), intent(out) :: exit_code

      integer(i32) :: N

      character(len=MAX_STR_LEN) :: iomessage
      integer(i32) :: i, iostatus
      character(len=MAX_STR_LEN) :: filename

      logical, save :: first_write = .true.

      N = params%point_count
      
      ! On the first call to 'write_output_file' we open the file in the correct mode
      ! and leave it open for future calls
      if (first_write) then
        filename = trim(params%output_file)

        open(unit=OUTPUT_FILE_UNIT, file=filename, status='old', action='write', position='append', iostat=iostatus, iomsg=iomessage)

        if (iostatus /= 0) then
          print '(a,a,a,/,a,a)', 'write_output_file: ERROR: Could not open file "', filename, '"', 'Message: ', trim(iomessage) 
          exit_code = FILE_ERROR
          return
        end if 
        first_write = .false.
      end if

      ! On all calls to 'write_output_file' we write the frame number and all relevant
      ! data about the system
      write(OUTPUT_FILE_UNIT, fmt=OUTPUT_FRAME_FMT) frame
      do i = 1,N
        write(OUTPUT_FILE_UNIT, fmt=OUTPUT_DATA_FMT) x_values(i), out_arrays%amplitude(i), &
          out_arrays%real_part(i), out_arrays%imag_part(i), out_arrays%density(i), potential(i) 
      end do

      exit_code = SUCCESS
      return
      
    end subroutine write_output_file

    subroutine close_output_file()
       ! Simple subroutine to close the file associated with the output file unit
       implicit none
       close(OUTPUT_FILE_UNIT)
    end subroutine close_output_file

    subroutine write_log_file(params, log_arrays, exit_code)
      ! Opens the log file, writes the relevant data to it and closes the file.
      ! This function is only called once at the end of every simulation
      implicit none

      type(SimulationParams), intent(in) :: params
      type(LoggingArrays), intent(in) :: log_arrays
      integer(excode), intent(out) :: exit_code

      character(len=MAX_STR_LEN) :: iomessage
      integer(i32) :: iostatus
      character(len=MAX_STR_LEN) :: filename
      integer(i32) :: framecount, i
      integer(label) :: imag

      filename = trim(params%log_file)
      framecount = params%step_count / params%write_interval + 1
      if (params%imag_time) then
        imag = 1
      else
        imag = 0
      end if

      open(unit=LOG_FILE_UNIT, file=filename, status='replace', action='write', iostat=iostatus, iomsg=iomessage)

      if (iostatus /= 0) then
        print '(A,A,A,/,A,A)', 'write_log_file: ERROR: Could not open file "', filename, '"', 'Message: ', trim(iomessage) 
        exit_code = FILE_ERROR
        return
      end if 

      write(LOG_FILE_UNIT, fmt='(A,/,A,/,A)') trim(LOG_TITLE_1), trim(LOG_TITLE_2), trim(LOG_TITLE_3)
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_FILES_TITLE)
      write(LOG_FILE_UNIT, fmt='(3(A,1X,A,/),A,1X,A)') &
        '# parameter file:', trim(params%param_file), &
        '# input file:    ', trim(params%input_file), &
        '# output file:   ', trim(params%output_file), &
        '# log file:      ', trim(params%log_file)
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_PARAMS_TITLE)
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')  '# method:        ', params%iter_method
      write(LOG_FILE_UNIT, fmt='(A,1X,L)')   '# imaginary time:', params%imag_time
      write(LOG_FILE_UNIT, fmt='(A,1X,L,/)') '# unit bounds:   ', params%unit_bounds
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')   '# time steps:', params%step_count
      write(LOG_FILE_UNIT, fmt='(A,1X,F8.6)') '# delta t:   ', params%delta_t
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')   '# interval:  ', params%write_interval
      write(LOG_FILE_UNIT, fmt='(A,1X,I0,/)') '# frames:    ', framecount
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')       '# points:    ', params%point_count
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')    '# x bound:   ', params%x_max
      write(LOG_FILE_UNIT, fmt='(A,1X,F8.6)')     '# delta x:   ', params%delta_x
      if (params%iter_method == Method%SPLIT_STEP) then
        write(LOG_FILE_UNIT, fmt='(A,1X,f8.6,/)') '# delta p:   ', params%delta_p
      end if
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')      '# wavefunction type:    ', params%wave_type
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction mass:    ', params%mass
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction charge:  ', params%charge
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction offset:  ', params%wave_offset
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction width:   ', params%wave_width
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6,/)') '# wavefunction momentum:', params%wave_momentum
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')    '# potential type:    ', params%pot_type
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)') '# potential offset:  ', params%pot_offset
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)') '# potential width:   ', params%pot_width
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)') '# potential strength:', params%pot_strength
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_RESULTS_TITLE)
      write(LOG_FILE_UNIT, fmt='(A)') '# TBD' ! TODO: cpu time etc
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_QUANTITIES_TITLE)
      write(LOG_FILE_UNIT, fmt=LOG_HEADER_FMT) framecount, params%iter_method, params%pot_type, imag
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_DATA_FIELDS)
      do i = 1, framecount
        write(LOG_FILE_UNIT, fmt=LOG_DATA_FMT) log_arrays%timestamps(i), log_arrays%energies(i), log_arrays%squared_norms(i)
      end do 

      close(LOG_FILE_UNIT)

      exit_code = SUCCESS
      return
    end subroutine write_log_file
    

end module read_write
