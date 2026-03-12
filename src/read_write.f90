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
  use iso_fortran_env, only: OUTPUT_UNIT
  use kinds
  use exit_codes
  use io_parameters
  use sim_parameters
  use calculate, only: normalize

  implicit none

  private

  public :: read_param_file, read_input_file, init_output_file, write_output_file, &
    close_output_file, write_log_file

  contains

    function strip(string) result(stripped)
      ! Simple function to remove leading and trailing whitespace for a string
      implicit none
      
      character(*), intent(in) :: string
      character(len=MAX_STR_LEN) :: stripped

      stripped = trim(adjustl(string))
      return
    end function strip

    function parse_method(iter_method) result(method_label)
      implicit none
      
      character(*), intent(in) :: iter_method
      integer(label) :: method_label

      select case (iter_method)
        case (PARAM_CN)
          method_label = Method%CRANK_NICOLSON 
        case (PARAM_SS)
          method_label = Method%SPLIT_STEP
        case default
          method_label = Method%INVALID
      end select

      return
    end function parse_method

    function parse_wave_type(wave_type) result(wave_label)
      implicit none
      
      character(*), intent(in) :: wave_type
      integer(label) :: wave_label

      select case (wave_type)
        case (PARAM_GAUSS)
          wave_label = WaveType%GAUSSIAN 
        case (PARAM_SINC)
          wave_label = WaveType%SINC
        case default
          wave_label = WaveType%INVALID
      end select

      return
    end function parse_wave_type

    function parse_pot_type(pot_type) result(pot_label)
      implicit none
      
      character(*), intent(in) :: pot_type
      integer(label) :: pot_label

      select case (pot_type)
        case (PARAM_ZERO)
          pot_label = PotType%ZERO  
        case (PARAM_BOX)
          pot_label = PotType%BOX
        case (PARAM_BARRIER)
          pot_label = PotType%BARRIER
        case (PARAM_WELL)
          pot_label = PotType%WELL
        case (PARAM_DWELL)
          pot_label = PotType%DOUBLE_WELL
        case (PARAM_HARMONIC)
          pot_label = PotType%HARMONIC
        case (PARAM_LINEAR)
          pot_label = PotType%LINEAR
        case (PARAM_HLINEAR)
          pot_label = PotType%HALF_LINEAR
        case (PARAM_ABS)
          pot_label = PotType%ABSOLUTE_VALUE
        case (PARAM_LOG)
          pot_label = PotType%LOGARITHMIC
        case (PARAM_COSH)
          pot_label = PotType%HYPERBOLIC_COSINE
        case (PARAM_MORSE)
          pot_label = PotType%MORSE
        case default
          pot_label = PotType%INVALID
      end select

      return
    end function parse_pot_type

    subroutine read_param_file(params, output_mode, exit_code)
      implicit none
      
      type(SimulationParams), intent(inout) :: params
      integer(label), intent(in) :: output_mode
      integer(excode), intent(out) :: exit_code

      character(len=MAX_STR_LEN) :: filename

      integer(i32) :: iostatus, line_num
      character(len=MAX_STR_LEN) :: iomessage
      logical :: break, invalid_type, invalid_value

      character(len=MAX_STR_LEN) :: line
      character(len=MAX_STR_LEN) :: param_type
      character(len=MAX_STR_LEN) :: param_name
      character(len=MAX_STR_LEN) :: param_value

      integer(label) :: given_type, type_label
      integer(i32) :: int_holder, separator_id, equals_id, comment_id
      real(r64) :: float_holder
      logical :: bool_holder
      character(len=PARAM_STR_LEN) :: str_holder

      filename = trim(params%param_file)

      open(unit=PARAM_FILE_UNIT, file=filename, status='old', action='read', iostat=iostatus, iomsg=iomessage)
      if (iostatus /= 0) then
        print '(A,/,A)', 'kettera: ERROR:', trim(iomessage) 
        exit_code = FILE_ERROR
        return
      end if 

      line_num = 0
      invalid_type = .false.
      invalid_value = .false.
      break = .false.

      print '(A)', '----- BEGIN READING PARAMETER FILE -----'

      do while (.not. break)
        
        line_num = line_num + 1
        given_type = NO_TYPE

        read(unit=PARAM_FILE_UNIT, fmt='(A)', iostat=iostatus, iomsg=iomessage) line

        if (is_iostat_end(iostatus) .or. is_iostat_eor(iostatus)) then
          if (output_mode == VERBOSE_OUTPUT) then
            print '(A,1X,I0,A,1X,I0)', 'kettera: INFO: EOF or EOR reached on line', line_num, &
              ', status code:', iostatus
          end if
          break = .true.
          cycle
        end if

        ! Removing leading and trailing whitespace from line
        line = strip(line)

        ! Skipping commented lines and empty lines
        if ((line(1:1) == PARAM_COMMENT) .or. (len_trim(line) == 0)) then
          cycle
        end if

        ! Determining the positions of the various data fields
        separator_id = index(line, PARAM_NAME_SEPARATOR, kind=i32)
        equals_id = index(line, PARAM_VALUE_SEPARATOR, kind=i32)
        comment_id = index(line, PARAM_COMMENT, kind=i32)

        if (separator_id == 0) then
          if (output_mode /= QUIET_OUTPUT) then
            print '(A,1X,I0,A,1X,A,/,A)', &
              'kettera: WARNING: Invalid input parameter specification on line', &
              line_num, ': Missing separator', PARAM_NAME_SEPARATOR, '(IGNORING THE LINE)'
          end if
            cycle
        end if

        if (equals_id == 0) then
          if (output_mode /= QUIET_OUTPUT) then
            print '(A,1X,I0,A,1X,A,/,A)', &
              'kettera: WARNING: Invalid input parameter specification on line', &
              line_num, ': Missing separator', PARAM_VALUE_SEPARATOR, '(IGNORING THE LINE)'
          end if
            cycle
        end if  

        ! Truncate the line to ignore any trailing comments if applicable
        if (comment_id /= 0) then
          line = line(:comment_id-1)
        end if
  
        param_type = strip(line(:separator_id-1))
        param_name = strip(line(separator_id+2:equals_id-1))
        param_value = strip(line(equals_id+2:))

        ! Check to see if param_value is 'DEF' in which case cycle to use the default value
        ! TODO: Make it so that the code recognizes whether the variable name is valid
        ! so that the program doesn't say it read in a variable that doesn't exist
        if (param_value == PARAM_DEF) then
          if (output_mode == VERBOSE_OUTPUT) then
            print '(A,1X,A)', &
              'kettera: INFO: Using default value for parameter', trim(param_name)
          end if
          cycle
        end if

        select case (param_type)
          case (PARAM_INT)
            read(param_value, fmt=*, iostat=iostatus, iomsg=iomessage) int_holder 
            given_type = INT_TYPE
          case (PARAM_FLOAT)
            read(param_value, fmt=*, iostat=iostatus, iomsg=iomessage) float_holder 
            given_type = FLOAT_TYPE
          case (PARAM_BOOL)
            read(param_value, fmt='(L)', iostat=iostatus, iomsg=iomessage) bool_holder
            given_type = BOOL_TYPE
          case (PARAM_STR)
            read(param_value, fmt='(A)', iostat=iostatus, iomsg=iomessage) str_holder
            given_type = STR_TYPE
          case default
            if (output_mode /= QUIET_OUTPUT) then
              print '(A,1X,I0,A,A,A,/,A)', &
                'kettera: WARNING: Invalid input parameter specification on line', &
                line_num, ': Unknown type "', trim(param_type), '"', '(IGNORING THE LINE)'
            end if
            cycle
        end select

        if (iostatus /= 0) then
          if (output_mode /= QUIET_OUTPUT) then
            print '(A,1X,I0,A,A,A,/,A,1X,A,/,A)', &
              'kettera: WARNING: Invalid input parameter specification on line', &
              line_num, ': Unreadable value "', trim(param_value), '"', 'Message:', trim(iomessage), &
              '(IGNORING THE LINE)'
          end if
          cycle
        end if

        name_select: select case (param_name)
          case (PARAM_METHOD)
            type_label = parse_method(str_holder)
            if (given_type /= STR_TYPE) then
              invalid_type = .true.
              exit name_select
            else if (type_label == Method%INVALID) then
              invalid_value = .true.
              exit name_select
            else
              params%iter_method = type_label
            end if
          case (PARAM_IMAG_TIME)
            if (given_type /= BOOL_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%imag_time = bool_holder
            end if
          case (PARAM_NORMAL)
            if (given_type /= BOOL_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%normal = bool_holder
            end if
          case (PARAM_ORTHO)
            if (given_type /= BOOL_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%ortho = bool_holder
            end if
          case (PARAM_UNIT_BOUNDS)
            if (given_type /= BOOL_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%unit_bounds = bool_holder
            end if
          case (PARAM_STEP_COUNT)
            if (given_type /= INT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%step_count = int_holder
            end if
          case (PARAM_WRITE_INTERVAL)
            if (given_type /= INT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%write_interval = int_holder
            end if
          case (PARAM_DELTA_T)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%delta_t = float_holder
            end if
          case (PARAM_POINT_COUNT)
            if (given_type /= INT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%point_count = int_holder
            end if
          case (PARAM_X_MAX)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%x_max = float_holder
            end if
          case (PARAM_WAVE_TYPE)
            type_label = parse_wave_type(str_holder)
            if (given_type /= STR_TYPE) then
              invalid_type = .true.
              exit name_select
            else if (type_label == WaveType%INVALID) then
              invalid_value = .true.
              exit name_select
            else
              params%wave_type = type_label
            end if
          case (PARAM_MASS)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%mass = float_holder
            end if
          case (PARAM_CHARGE)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%charge = float_holder
            end if
          case (PARAM_MOMENTUM)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%momentum = float_holder
            end if
          case (PARAM_WAVE_OFFSET)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%wave_offset = float_holder
            end if
          case (PARAM_WAVE_WIDTH)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%wave_width = float_holder
            end if
          case (PARAM_POT_TYPE)
            type_label = parse_pot_type(str_holder)
            if (given_type /= STR_TYPE) then
              invalid_type = .true.
              exit name_select
            else if (type_label == PotType%INVALID) then
              invalid_value = .true.
              exit name_select
            else
              params%pot_type = type_label
            end if
          case (PARAM_POT_OFFSET)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%pot_offset = float_holder
            end if
          case (PARAM_POT_WIDTH)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%pot_width = float_holder
            end if
          case (PARAM_POT_STRENGTH)
            if (given_type /= FLOAT_TYPE) then
              invalid_type = .true.
              exit name_select
            else
              params%pot_strength = float_holder
            end if
          case default
            if (output_mode /= QUIET_OUTPUT) then
              print '(A,A,A,1x,I0,1X,A)', &
                'kettera: WARNING: Invalid variable name "', &
                trim(param_name), '" on line', line_num, '(IGNORING THE LINE)'
            end if
            cycle
        end select name_select
        
        if (invalid_type) then
          if (output_mode /= QUIET_OUTPUT) then
            print '(A,1X,A,1X,A,1X,I0,1X,A)', &
              'kettera: WARNING: Invalid type specification for variable', &
              trim(param_name), 'on line', line_num, '(IGNORING THE LINE)'
          end if

          invalid_type = .false.
          cycle
        end if

        if (invalid_value) then
          if (output_mode /= QUIET_OUTPUT) then
            print '(A,1X,A,1X,A,1X,I0,1X,A)', &
              'kettera: WARNING: Invalid value for variable', &
              trim(param_name), 'on line', line_num, '(IGNORING THE LINE)'
          end if

          invalid_value = .false.
          cycle
        end if
        
        if (output_mode == VERBOSE_OUTPUT) then
          write(OUTPUT_UNIT, fmt='(A,1X,A,1X,A,1X)', advance='no') &
            'kettera: INFO: Read in parameter', trim(param_name), ':='
          select case (given_type)
            case (INT_TYPE)
              print '(I0)', int_holder
            case (FLOAT_TYPE)
              print '(F10.4)', float_holder
            case (BOOL_TYPE)
              print '(L)', bool_holder
            case (STR_TYPE)
              print '(A,1x,A,I0,A)', trim(str_holder), '(=', type_label, ')'
            case default
              print '(A)', 'kettera: INFO: Type not recognized'
          end select
        end if

      end do

      close(PARAM_FILE_UNIT)

      print '(A)', '----- END READING PARAMETER FILE -----'

      exit_code = SUCCESS
      return
    end subroutine read_param_file

    subroutine read_input_file(params, wavefunction, output_mode, exit_code)
      implicit none
      
      type(SimulationParams), intent(inout) :: params
      complex(r64), intent(inout) :: wavefunction(params%point_count)
      integer(label), intent(in) :: output_mode
      integer(excode), intent(out) :: exit_code

      character(len=MAX_STR_LEN) :: filename

      integer(i32) :: iostatus, line_num, array_index
      character(len=MAX_STR_LEN) :: iomessage
      logical :: break, invalid_value

      character(len=MAX_STR_LEN) :: line
      character(len=MAX_STR_LEN) :: point_count_str, x_max_str
      integer(i32) :: point_count_int
      real(r64) :: x_max_float
      integer(i32) :: equals1_id, equals2_id, semicolon_id
      logical :: params_read
      character(len=MAX_STR_LEN) :: real_str, imag_str
      real(r64) :: real_holder, imag_holder

      filename = trim(params%input_file)

      if (trim(filename) == NO_INPUT_FILE) then
        exit_code = SUCCESS
        return
      end if

      open(unit=INPUT_FILE_UNIT, file=filename, status='old', action='read', iostat=iostatus, iomsg=iomessage)
      if (iostatus /= 0) then
        print '(A,/,A)', 'kettera: ERROR:', trim(iomessage) 
        exit_code = FILE_ERROR
        return
      end if 

      line_num = 0
      array_index = 1
      invalid_value = .false.
      params_read = .false.
      break = .false.

      print '(A)', '----- BEGIN READING INPUT FILE -----'

      do while (.not. break)
        
        line_num = line_num + 1

        read(unit=INPUT_FILE_UNIT, fmt='(A)', iostat=iostatus, iomsg=iomessage) line

        if (is_iostat_end(iostatus) .or. is_iostat_eor(iostatus)) then
          if (output_mode == VERBOSE_OUTPUT) then
            print '(A,1X,I0,A,1X,I0)', 'kettera: INFO: EOF or EOR reached on line', line_num, &
              ', status code:', iostatus
          end if
          break = .true.
          cycle
        end if

        ! Removing leading and trailing whitespace from line
        line = strip(line)

        ! Skipping commented lines, header line and empty lines
        if ((line(1:1) == INPUT_COMMENT) .or. (line == INPUT_HEADER) .or. (len_trim(line) == 0)) then
          cycle
        end if

        ! Checking that the parameters in the input file match the given parameters
        ! so that the wave function can fit properly into the space
        if (.not. params_read) then
          equals1_id = index(line, INPUT_VALUE_SEPARATOR, kind=i32)
          equals2_id = index(line, INPUT_VALUE_SEPARATOR, kind=i32, back=.true.)
          semicolon_id = index(line, INPUT_FIELD_SEPARATOR, kind=i32)

          if ((equals1_id == 0) .or. (equals2_id == 0) .or. (semicolon_id == 0)) then
            print '(A,1X,I0,1X,A,/,A)', &
              'kettera: ERROR: Invalid parameter specification line', line_num, &
              'in input file:', line
            exit_code = GIVEN_PARAMETER_ERROR
            return
          end if

          point_count_str = strip(line(equals1_id+1:semicolon_id-1)) 
          x_max_str = strip(line(equals2_id+1:))

          read(point_count_str, fmt=*, iostat=iostatus, iomsg=iomessage) point_count_int
          if (iostatus /= 0) then
            print '(A,1X,I0,A,A,A,/,A,1X,A)', &
              'kettera: ERROR: Invalid parameter specification line', &
              line_num, ': Unreadable value "', trim(point_count_str), '"', 'Message:', trim(iomessage)
            exit_code = GIVEN_PARAMETER_ERROR
            return
          end if
          read(x_max_str, fmt=*, iostat=iostatus, iomsg=iomessage) x_max_float
          if (iostatus /= 0) then
            print '(A,1X,I0,A,A,A,/,A,1X,A)', &
              'kettera: ERROR: Invalid parameter specification line', &
              line_num, ': Unreadable value "', trim(x_max_str), '"', 'Message:', trim(iomessage)
            exit_code = GIVEN_PARAMETER_ERROR
            return
          end if
  
          if (point_count_int /= params%point_count) then
            print '(A,1X,I0,1X,A,1X,I0)', 'kettera: ERROR: Input file PointCount', &
              point_count_int, 'is different from parameter PointCount', params%point_count
            exit_code = GIVEN_PARAMETER_ERROR
            return
          end if

          if (abs(x_max_float - params%x_max) >= EPS) then
            print '(A,1X,F6.2,1X,A,1X,F6.2)', 'kettera: ERROR: Input file XMax', &
              x_max_float, 'is different from parameter XMax', params%x_max
            exit_code = GIVEN_PARAMETER_ERROR
            return
          end if

          params_read = .true.
          cycle
        end if
        
        semicolon_id = index(line, INPUT_FIELD_SEPARATOR, kind=i32)

        if (semicolon_id == 0) then
          print '(A,1X,I0,1X,A,/,A)', &
            'kettera: ERROR: Invalid data line', line_num, &
            'in input file:', line
          exit_code = GIVEN_PARAMETER_ERROR
          return
        end if

        real_str = strip(line(:semicolon_id-1))
        imag_str = strip(line(semicolon_id+1:))


        read(real_str, fmt=*, iostat=iostatus, iomsg=iomessage) real_holder
        if (iostatus /= 0) then
          print '(A,1X,I0,A,A,A,/,A,1X,A)', &
            'kettera: ERROR: Invalid data line', &
            line_num, ': Unreadable value "', trim(real_str), '"', 'Message:', trim(iomessage)
          exit_code = GIVEN_PARAMETER_ERROR
          return
        end if
        read(imag_str, fmt=*, iostat=iostatus, iomsg=iomessage) imag_holder
        if (iostatus /= 0) then
          print '(A,1X,I0,A,A,A,/,A,1X,A)', &
            'kettera: ERROR: Invalid parameter specification line', &
            line_num, ': Unreadable value "', trim(imag_str), '"', 'Message:', trim(iomessage)
          exit_code = GIVEN_PARAMETER_ERROR
          return
        end if

        if (array_index <= params%point_count) then
          wavefunction(array_index) = cmplx(real_holder, imag_holder, kind=r64)
          array_index = array_index + 1
        else
          print '(A)', 'kettera: ERROR: Input file contains more data that can fit inside allocated space!'
          exit_code = GIVEN_PARAMETER_ERROR
          return
        end if

      end do

      close(INPUT_FILE_UNIT)

      print '(A)', '----- END READING INPUT FILE -----'

      if (params%normal) then
        call normalize(wavefunction)
      end if
      
      exit_code = SUCCESS
      return

    end subroutine read_input_file

    subroutine init_output_file(params, exit_code)
      implicit none

      type(SimulationParams), intent(in) :: params
      integer(excode), intent(out) :: exit_code

      character(len=MAX_STR_LEN) :: iomessage
      integer(i32) :: iostatus
      character(len=MAX_STR_LEN) :: filename
      integer(i32) :: framecount, pointcount
      integer(label) :: method, pot, imag
      real(r64) :: x_bound


      filename = trim(params%output_file)
      framecount = params%frame_count
      pointcount = params%point_count
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
        print '(a,a,a,/,a,a)', 'kettera: ERROR: Could not open file "', filename, '"', 'Message: ', trim(iomessage) 
        exit_code = FILE_ERROR
        return
      end if 

      write(OUTPUT_FILE_UNIT, fmt='(A,/,A,/,A)') trim(OUTPUT_TITLE_1), trim(OUTPUT_TITLE_2), &
        trim(OUTPUT_TITLE_3)
      write(OUTPUT_FILE_UNIT, fmt=trim(OUTPUT_HEADER_FMT)) framecount, method, pointcount, x_bound, pot, imag
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
          print '(a,a,a,/,a,a)', 'kettera: ERROR: Could not open file "', filename, '"', 'Message: ', trim(iomessage) 
          exit_code = FILE_ERROR
          return
        end if 
        first_write = .false.
      end if

      ! On all calls to 'write_output_file' we write the frame number and all relevant
      ! data about the system
      write(OUTPUT_FILE_UNIT, fmt=OUTPUT_FRAME_FMT) frame
      do i = 1,N
        write(OUTPUT_FILE_UNIT, fmt=OUTPUT_DATA_FMT) x_values(i), out_arrays%real_part(i), &
          out_arrays%imag_part(i), out_arrays%amplitude(i), out_arrays%density(i), potential(i) 
      end do

      exit_code = SUCCESS
      return
      
    end subroutine write_output_file

    subroutine close_output_file()
       ! Simple subroutine to close the file associated with the output file unit
       implicit none
       close(OUTPUT_FILE_UNIT)
    end subroutine close_output_file

    subroutine write_log_file(params, log_arrays, fin_results, exit_code)
      ! Opens the log file, writes the relevant data to it and closes the file.
      ! This function is only called once at the end of every simulation
      implicit none

      type(SimulationParams), intent(in) :: params
      type(LoggingArrays), intent(in) :: log_arrays
      type(FinalResults), intent(in) :: fin_results
      integer(excode), intent(out) :: exit_code

      character(len=MAX_STR_LEN) :: iomessage
      integer(i32) :: iostatus
      character(len=MAX_STR_LEN) :: filename
      integer(i32) :: i
      integer(label) :: imag

      filename = trim(params%log_file)
      if (params%imag_time) then
        imag = 1
      else
        imag = 0
      end if


      open(unit=LOG_FILE_UNIT, file=filename, status='replace', action='write', iostat=iostatus, iomsg=iomessage)

      if (iostatus /= 0) then
        print '(A,A,A,/,A,A)', 'kettera: ERROR: Could not open file "', filename, '"', 'Message: ', trim(iomessage) 
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
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')  '# method:           ', params%iter_method
      write(LOG_FILE_UNIT, fmt='(A,1X,L)')   '# imaginary time:   ', params%imag_time
      write(LOG_FILE_UNIT, fmt='(A,1X,L)')   '# normalization:    ', params%normal
      write(LOG_FILE_UNIT, fmt='(A,1X,L)')   '# orthogonalization:', params%ortho
      write(LOG_FILE_UNIT, fmt='(A,1X,L,/)') '# unit bounds:   ', params%unit_bounds
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')   '# time steps:', params%step_count
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')   '# interval:  ', params%write_interval
      write(LOG_FILE_UNIT, fmt='(A,1X,I0,/)') '# frames:    ', params%frame_count
      write(LOG_FILE_UNIT, fmt='(A,1X,F8.6)') '# delta t:   ', params%delta_t
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')       '# points:    ', params%point_count
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')    '# x max:     ', params%x_max
      write(LOG_FILE_UNIT, fmt='(A,1X,F8.6)')     '# delta x:   ', params%delta_x
      if (params%iter_method == Method%SPLIT_STEP) then
        write(LOG_FILE_UNIT, fmt='(A,1X,f8.6,/)') '# delta p:   ', params%delta_p
      else
        write(LOG_FILE_UNIT, fmt='()')
      end if
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')      '# wavefunction type:    ', params%wave_type
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction mass:    ', params%mass
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction charge:  ', params%charge
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction momentum:', params%momentum
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)')   '# wavefunction offset:  ', params%wave_offset
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6,/)') '# wavefunction width:   ', params%wave_width
      write(LOG_FILE_UNIT, fmt='(A,1X,I0)')    '# potential type:    ', params%pot_type
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)') '# potential offset:  ', params%pot_offset
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)') '# potential width:   ', params%pot_width
      write(LOG_FILE_UNIT, fmt='(A,1X,F10.6)') '# potential strength:', params%pot_strength
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_RESULTS_TITLE)
      if (params%imag_time) then
        write(LOG_FILE_UNIT, fmt='(A)') &
          '# NOTE: imaginary time was used so change in energy is due to exponential decay'
        write(LOG_FILE_UNIT, fmt='(A)') &
          '# and change in norm squared is negligible due to constant normalization during the iteration'
      else
        write(LOG_FILE_UNIT, fmt='(A)') '# NOTE: change in energy/norm squared is due to numerical error'
      end if
      write(LOG_FILE_UNIT, fmt='(A,1X,F11.6,1X,A)') '# loop time elapsed:     ', &
        fin_results%loop_time_elapsed, '(seconds)'
      write(LOG_FILE_UNIT, fmt='(A,1X,F20.15,1X,A)')    '# initial energy:        ', &
        fin_results%initial_energy, '(hbar=1 units)'
      write(LOG_FILE_UNIT, fmt='(A,1X,F20.15,1X,A)')    '# final energy:          ', &
        fin_results%final_energy, '(hbar=1 units)'
      write(LOG_FILE_UNIT, fmt='(A,1X,F20.15,1X,A)')    '# change in energy:      ', &
        fin_results%delta_energy, '(hbar=1 units)'
      write(LOG_FILE_UNIT, fmt='(A,1X,F20.15,1X,A)')    '# initial norm squared:  ', &
        fin_results%initial_norm_squared, '(unitless)'
      write(LOG_FILE_UNIT, fmt='(A,1X,F20.15,1X,A)')    '# final norm squared:    ', &
        fin_results%final_norm_squared, '(unitless)'
      write(LOG_FILE_UNIT, fmt='(A,1X,F20.15,1X,A)')    '# change in norm squared:', &
        fin_results%delta_norm_squared, '(unitless)'
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_QUANTITIES_TITLE)
      write(LOG_FILE_UNIT, fmt=LOG_HEADER_FMT) params%frame_count, params%iter_method, params%pot_type, imag
      write(LOG_FILE_UNIT, fmt='(A)') trim(LOG_DATA_FIELDS)
      do i = 1, params%frame_count
        write(LOG_FILE_UNIT, fmt=LOG_DATA_FMT) log_arrays%timestamps(i), log_arrays%energies(i), log_arrays%squared_norms(i)
      end do 

      close(LOG_FILE_UNIT)

      exit_code = SUCCESS
      return
    end subroutine write_log_file

    

end module read_write
