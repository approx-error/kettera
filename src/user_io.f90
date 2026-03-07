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

module user_io
  ! Subroutines and functions for parsing user input and showing output to the user
  use kinds
  use exit_codes
  use io_parameters
  use sim_parameters

  implicit none

  private

  public :: parse_user_input

  contains

    subroutine show_usage
      implicit none
      
      print '(A)', 'Usage: kettera -<s|c>[vq] [-p <param-file>] [-i <input-file>] [-o <output-file>] [-l <log-file>]'
      print '(A)', '   or: kettera [--help | --version]'
      print '(A)', "Run 'kettera --help' for more information"

    end subroutine show_usage

    subroutine show_help
      implicit none
      
      print '(A)', 'Usage: kettera -<s|c>[vq] [-p <param-file>] [-i [input-file]] [-o <output-file>] [-l <log-file>]'
      print '(A)', '   or: kettera [--help | --version]'
      print '(A)', 'Numerically solve the time-dependent Schrödinger equation using the'
      print '(A)', 'provided parameters and output the results into an output file and a log file'
      print '(/,A)', 'Options surrounded by square brackets [] are optional while options surrounded'
      print '(A)',   'by angle brackets <> are mandatory. x | y means either x or y must be chosen'
      print '(/,A)', 'If a file is left unspecified, the default filenames will be searched for in the current directory'
      print '(/,A)', 'Options:'
      print '(A)',   '  -s, --simulate         use the input file(s) to run a simulation and output'
      print '(A)',   '                         the results into the output files'
      print '(A)',   '  -c, --check            check the validity of the input file(s) and exit'
      print '(A)',   '  -v, --verbose          output more information during execution'
      print '(A)',   '  -q, --quiet            show output only when processes start/finish or if'
      print '(A)',   '                         an error occurs'
      print '(A)',   '  -p, --params <file>    choose <file> as the input parameter file' 
      print '(A)',   "                         default parameter file name is 'params.in'"
      print '(A)',   '  -i, --input [file]     choose [file] as the input wavefunction file'
      print '(A)',   "                         default input wavefunction file name is 'wave.in'"
      print '(A)',   '                         but default behaviour is no input file at all so'
      print '(A)',   '                         to use the default input wavefunction the flag -i'
      print '(A)',   '                         or --input must be provided'
      print '(A)',   '  -o, --output <file>    choose <file> as the output wavefunction file'
      print '(A)',   "                         default output wavefunction file name is 'wave.out'"
      print '(A)',   '  -l, --log <file>       choose <file> as the output log file'
      print '(A)',   "                         default output log file name is 'wave.log'"
      print '(A)',   '      --help             show this help and exit'
      print '(A)',   '      --version          show version information and exit'
     
    end subroutine show_help

    subroutine show_version
      implicit none
      
      print '(A)', 'kettera v0.5.0'
      print '(A)', 'Copyright (C) 2026 Juuso Kaarela'
      print '(/,A)', 'Kettera is a numerical solver and visualizer for'
      print '(A)',   'the time-dependent Schrödinger equation'
      print '(A)', 'License GPLv3+: GNU GPL version 3 or later <https://gnu.org/licenses/gpl.html>'
      print '(A)', 'This is free software: you are free to change and redistribute it'
      print '(A)', 'There is NO WARRANTY, to the extent permitted by the law'

    end subroutine show_version
  
    subroutine parse_user_input(params, output_mode, exit_code)
      implicit none
      
      type(SimulationParams), intent(out) :: params
      integer(label), intent(out) :: output_mode
      integer(excode), intent(out) :: exit_code

      integer(i32) :: argc, i
      integer(i32) :: arg_len
      integer(i32) :: arg_stat
      character(MAX_STR_LEN) :: arg, prev_arg, next_arg

      logical :: quiet_flag, verbose_flag, params_flag, input_flag, output_flag, log_flag
      
      logical :: argument

      argc = command_argument_count()

      if (argc == 0) then
        call show_usage
        exit_code = SUCCESS
        return
      end if

      exit_code = SUCCESS

      verbose_flag = .false.
      quiet_flag = .false.
      params_flag = .false.
      input_flag = .false.
      output_flag = .false.
      log_flag = .false.

      argument = .false.
      prev_arg = ''

      do i = 1, argc
        call get_command_argument(i, arg, arg_len, arg_stat)

        if (arg_stat > 0) then
          print '(A)', 'parse_user_input: ERROR: Failed to retrieve command arguments'
          exit_code = USER_INPUT_ERROR
          return
        end if

        if (arg_stat == -1) then
          print '(A,I0,A,1X,I0)', 'parse_user_input: ERROR: Maximum length for arguments (', &
            MAX_STR_LEN, ' characters) exceeded in argument', i
          exit_code = USER_INPUT_ERROR
          return
        end if

        if (argument) then

          if (arg(:1) == '-') then
            print '(A,A,A,/,A)', "parse_user_input: ERROR: Expected file name but got flag '", &
              trim(arg), "' instead", "Run 'kettera --help' to see correct usage"
            exit_code = USER_INPUT_ERROR
            return
          end if

          if (trim(prev_arg) == PARAMS_SHORT .or. trim(prev_arg) == PARAMS_LONG) then
            print '(A,1X,A)', 'DEBUG: param file:', trim(arg)
            params%param_file = arg
          else if (trim(prev_arg) == INPUT_SHORT .or. trim(prev_arg) == INPUT_LONG) then
            print '(A,1X,A)', 'DEBUG: input file:', trim(arg)
            params%input_file = arg
          else if (trim(prev_arg) == OUTPUT_SHORT .or. trim(prev_arg) == OUTPUT_LONG) then
            print '(A,1X,A)', 'DEBUG: output file:', trim(arg)
            params%output_file = arg
          else if (trim(prev_arg) == LOG_SHORT .or. trim(prev_arg) == LOG_LONG) then
            print '(A,1X,A)', 'DEBUG: log file:', trim(arg)
            params%log_file = arg
          else
            print '(A,1X,A,1X,A)', 'parse_user_input: ERROR: Option', trim(prev_arg), &
              'does not take any arguments'
            exit_code = USER_INPUT_ERROR
            return
          end if
          argument = .false.
        end if

        if (arg(:1) == '-')
          if (trim(arg) == VERBOSE_SHORT .or. trim(arg) == VERBOSE_LONG) then
            print '(A)', 'DEBUG: verbose flag'
            if (verbose_flag) then
              print '(A)', 'parse_user_input: WARNING: Verbose flag set multiple times'
              exit_code = USER_INPUT_WARNING
            end if

            if (quiet_flag) then
              print '(A)', 'parse_user_input: WARNING: Setting verbose flag overrides quiet flag'
              exit_code = USER_INPUT_WARNING
            end if

            verbose_flag = .true.

          else if (trim(arg) == QUIET_SHORT .or. trim(arg) == QUIET_LONG) then
            print '(A)', 'DEBUG: quiet flag'
            if (quiet_flag) then
              print '(A)', 'parse_user_input: WARNING: Quiet flag set multiple times'
              exit_code = USER_INPUT_WARNING
            end if

            if (verbose_flag) then
              print '(A)', 'parse_user_input: WARNING: Quiet flag will be ignored since verbose flag was already set'
              exit_code = USER_INPUT_WARNING
            end if
            
            quiet_flag = .true.
            
          else if (trim(arg) == PARAMS_SHORT .or. trim(arg) == PARAMS_LONG) then
            print '(A)', 'DEBUG: params flag'
            if (params_flag) then
              print '(A)', 'parse_user_input: ERROR: Cannot set parameter file multiple times'
              exit_code = USER_INPUT_ERROR
              return
            end if

            params_flag = .true.
            argument = .true.
            prev_arg = arg
          else if (trim(arg) == INPUT_SHORT .or. trim(arg) == INPUT_LONG) then
            print '(A)', 'DEBUG: input flag'
            if (input_flag) then
              print '(A)', 'parse_user_input: ERROR: Cannot set input file multiple times'
              exit_code = USER_INPUT_ERROR
              return
            end if

            input_flag = .true.
            argument = .false.
            params%input_file = DEFAULT_INPUT_FILE
            
            if (i < argc) then
              call get_command_argument(i+1, next_arg, arg_len, arg_stat)

              if (arg_stat > 0) then
                print '(A)', 'parse_user_input: ERROR: Failed to retrieve command arguments'
                exit_code = USER_INPUT_ERROR
                return
              end if

              if (arg_stat == -1) then
                print '(A,I0,A,1X,I0)', 'parse_user_input: ERROR: Maximum length for arguments (', &
                  MAX_STR_LEN, ' characters) exceeded in argument', i+1
                exit_code = USER_INPUT_ERROR
                return
              end if

              if (next_arg(:1) /= '-') then
                argument = .true.
                prev_arg = arg
              end if
            end if

          else if (trim(arg) == OUTPUT_SHORT .or. trim(arg) == OUTPUT_LONG) then
            print '(A)', 'DEBUG: output flag'
            if (output_flag) then
              print '(A)', 'parse_user_input: ERROR: Cannot set output file multiple times'
              exit_code = USER_INPUT_ERROR
              return
            end if

            output_flag = .true.
            argument = .true.
            prev_arg = arg
          else if (trim(arg) == LOG_SHORT .or. trim(arg) == LOG_LONG) then
            print '(A)', 'DEBUG: log flag'
            if (log_flag) then
              print '(A)', 'parse_user_input: ERROR: Cannot set log file multiple times'
              exit_code = USER_INPUT_ERROR
              return
            end if

            log_flag = .true.
            argument = .true.
            prev_arg = arg
          else
            print '(A,A,A,/,A)', "parse_user_input: ERROR: Unrecognized flag '", trim(arg), "'", &
              "Run 'kettera --help' to see a list of valid flags"
            exit_code = USER_INPUT_ERROR
            return
          end if 
        else
          print '(A,A,A,/,A)', "parse_user_input: ERROR: Unrecognized argument '", trim(arg), "'", &
            "Run 'kettera --help' to see a list of valid arguments"
          exit_code = USER_INPUT_ERROR
          return
        end if

        if (i == argc .and. argument) then
          print '(A,A,A,/,A)', "parse_user_input: ERROR: Expected argument for flag '", trim(arg), &
            "' but reached end of arguments", "Run 'kettera --help' to see correct usage"
          exit_code = USER_INPUT_ERROR
        end if
      end do

      if (verbose_flag) then
        output_mode = VERBOSE_OUTPUT
      else if (quiet_flag) then
        output_mode = QUIET_OUTPUT
      else
        output_mode = NORMAL_OUTPUT
      end if

      return
    end subroutine parse_user_input
end module user_io
