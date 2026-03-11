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
  
    subroutine parse_user_input(params, output_mode, execution_mode, exit_code)
      implicit none
      
      type(SimulationParams), intent(out) :: params
      integer(label), intent(out) :: output_mode
      integer(label), intent(out) :: execution_mode
      integer(excode), intent(out) :: exit_code

      integer(i32) :: argc, i, j
      integer(i32) :: arg_len
      integer(i32) :: arg_stat
      character(MAX_STR_LEN) :: arg, prev_arg, next_arg

      logical :: quiet_flag, verbose_flag, simulate_flag, check_flag, &
      params_flag, input_flag, output_flag, log_flag
      
      logical :: argument, long_flag, short_flag, multiple_flags

      argc = command_argument_count()
      print '(A,1X,I0)', 'DEBUG: argc = ', argc

      if (argc == 0) then
        call show_usage
        exit_code = SUCCESS
        return
      end if

      exit_code = SUCCESS

      simulate_flag = .false.
      check_flag = .false.
      verbose_flag = .false.
      quiet_flag = .false.
      params_flag = .false.
      input_flag = .false.
      output_flag = .false.
      log_flag = .false.

      argument = .false.
      prev_arg = ''

      do i = 1, argc
        print '(A)', 'DEBUG: new arg'
        short_flag = .false.
        long_flag = .false.
        multiple_flags = .false.

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

          select case (trim(prev_arg))
            case (PARAMS_SHORT)
              print '(A,1X,A)', 'DEBUG: param file:', trim(arg)
              params%param_file = arg
            case (INPUT_SHORT) 
              print '(A,1X,A)', 'DEBUG: input file:', trim(arg)
              params%input_file = arg
            case (OUTPUT_SHORT)
              print '(A,1X,A)', 'DEBUG: output file:', trim(arg)
              params%output_file = arg
            case (LOG_SHORT)
              print '(A,1X,A)', 'DEBUG: log file:', trim(arg)
              params%log_file = arg
            case default
              print '(A,1X,A,1X,A)', 'parse_user_input: ERROR: Option', trim(prev_arg), &
                'does not take any arguments'
              exit_code = USER_INPUT_ERROR
              return
          end select
          argument = .false.
          cycle
        end if

        if (arg(1:2) == '--') then
          arg = arg(3:)
          select case (trim(arg))
            case (HELP_FLAG)
              print '(A)', 'DEBUG: help'
              call show_help
              exit_code = SUCCESS
              return
            case (VERSION_FLAG)
              print '(A)', 'DEBUG: version'
              call show_version
              exit_code = SUCCESS
              return
            case (SIMULATE_LONG)
              print '(A)', 'DEBUG: simulate'
              arg = SIMULATE_SHORT
            case (CHECK_LONG)
              print '(A)', 'DEBUG: check'
              arg = CHECK_SHORT
            case (VERBOSE_LONG)
              print '(A)', 'DEBUG: verbose'
              arg = VERBOSE_SHORT
            case (QUIET_SHORT)
              print '(A)', 'DEBUG: quiet'
              arg = QUIET_SHORT
            case (PARAMS_LONG)
              print '(A)', 'DEBUG: params'
              arg = PARAMS_SHORT
            case (INPUT_LONG)
              print '(A)', 'DEBUG: input'
              arg = INPUT_SHORT
            case (OUTPUT_LONG)
              print '(A)', 'DEBUG: output'
              arg = OUTPUT_SHORT
            case (LOG_LONG)
              print '(A)', 'DEBUG: log'
              arg = LOG_SHORT
            case default
              print '(A,A,A,/,A)', "parse_user_input: ERROR: Unrecognized flag '", trim(arg), "'", &
                "Run 'kettera --help' to see a list of valid flags"
              exit_code = USER_INPUT_ERROR
              return
          end select
          
          long_flag = .true.
        end if

        if (arg(:1) == '-') then
          arg = arg(2:)
          if (len_trim(arg) == 1) then
            short_flag = .true.
          else
            multiple_flags = .true.
          end if
        end if

        if (long_flag .or. short_flag) then
          select case (trim(arg))
            case (SIMULATE_SHORT)
              print '(A)', 'DEBUG: s'
              if (simulate_flag) then
                print '(A)', "parse_user_input: WARNING: Flag 'simulate' was set multiple times"
                exit_code = USER_INPUT_WARNING
              end if

              if (check_flag) then
                print '(A)', "parse_user_input: ERROR: Cannot set flag 'simulate' as flag 'check' was already set"
                exit_code = USER_INPUT_ERROR
                return
              end if

              simulate_flag = .true.

            case (CHECK_SHORT)
              print '(A)', 'DEBUG: c'
              if (check_flag) then
                print '(A)', "parse_user_input: WARNING: Flag 'check' was set multiple times"
                exit_code = USER_INPUT_WARNING
              end if

              if (simulate_flag) then
                print '(A)', "parse_user_input: ERROR: Cannot set flag 'check' as flag 'simulate' was already set"
                exit_code = USER_INPUT_ERROR
                return
              end if

              check_flag = .true.

            case (VERBOSE_SHORT)
              print '(A)', 'DEBUG: v'
              if (verbose_flag) then
                print '(A)', "parse_user_input: WARNING: Flag 'verbose' was set multiple times"
                exit_code = USER_INPUT_WARNING
              end if

              if (quiet_flag) then
                print '(A)', "parse_user_input: WARNING: Setting flag 'verbose' flag overrides flag 'quiet'"
                exit_code = USER_INPUT_WARNING
              end if

              verbose_flag = .true.

            case (QUIET_SHORT)
              print '(A)', 'DEBUG: q'
              if (quiet_flag) then
                print '(A)', "parse_user_input: WARNING: Flag 'quiet' was set multiple times"
                exit_code = USER_INPUT_WARNING
              end if

              if (verbose_flag) then
                print '(A)', "parse_user_input: WARNING: Flag 'quiet' will be ignored since flag 'verbose' was already set"
                exit_code = USER_INPUT_WARNING
              end if

              quiet_flag = .true.

            case (PARAMS_SHORT)
              print '(A)', 'DEBUG: p'
              if (params_flag) then
                print '(A)', 'parse_user_input: ERROR: Cannot set parameter file multiple times'
                exit_code = USER_INPUT_ERROR
                return
              end if

              params_flag = .true.
              argument = .true.
              prev_arg = arg

            case (INPUT_SHORT)
              print '(A)', 'DEBUG: i'
              if (input_flag) then
                print '(A)', 'parse_user_input: ERROR: Cannot set input file multiple times'
                exit_code = USER_INPUT_ERROR
                return
              end if

              input_flag = .true.
              argument = .false.
              params%input_file = DEFAULT_INPUT_FILE
              
              ! If there are arguments after the -i/--input flag check whether they are 
              ! flags in which case we use the default input file or whether they are
              ! not flags in which case the next argument should be a file name
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

            case (OUTPUT_SHORT)
              print '(A)', 'DEBUG: o'
              if (output_flag) then
                print '(A)', 'parse_user_input: ERROR: Cannot set output file multiple times'
                exit_code = USER_INPUT_ERROR
                return
              end if

              output_flag = .true.
              argument = .true.
              prev_arg = arg
            
            case (LOG_SHORT)
              print '(A)', 'DEBUG: l'
              if (log_flag) then
                print '(A)', 'parse_user_input: ERROR: Cannot set log file multiple times'
                exit_code = USER_INPUT_ERROR
                return
              end if

              log_flag = .true.
              argument = .true.
              prev_arg = arg

            case default
              print '(A,A,A,/,A)', "parse_user_input: ERROR: Unrecognized flag '", trim(arg), "'", &
                "Run 'kettera --help' to see a list of valid flags"
              exit_code = USER_INPUT_ERROR
              return
          end select 

        end if

        if (multiple_flags) then
          
          do j = 1, len_trim(arg)
            select case (arg(j:j))
              case (SIMULATE_SHORT)
                print '(A)', 'DEBUG: s'
                if (simulate_flag) then
                  print '(A)', "parse_user_input: WARNING: Flag 'simulate' was set multiple times"
                  exit_code = USER_INPUT_WARNING
                end if

                if (check_flag) then
                  print '(A)', "parse_user_input: ERROR: Cannot set flag 'simulate' as flag 'check' was already set"
                  exit_code = USER_INPUT_ERROR
                  return
                end if

                simulate_flag = .true.
              case (CHECK_SHORT)
                print '(A)', 'DEBUG: c'
                if (check_flag) then
                  print '(A)', "parse_user_input: WARNING: Flag 'check' was set multiple times"
                  exit_code = USER_INPUT_WARNING
                end if

                if (simulate_flag) then
                  print '(A)', "parse_user_input: ERROR: Cannot set flag 'check' as flag 'simulate' was already set"
                  exit_code = USER_INPUT_ERROR
                  return
                end if

                check_flag = .true.
              case (VERBOSE_SHORT)
                print '(A)', 'DEBUG: v'
                if (verbose_flag) then
                  print '(A)', "parse_user_input: WARNING: Flag 'verbose' was set multiple times"
                  exit_code = USER_INPUT_WARNING
                end if

                if (quiet_flag) then
                  print '(A)', "parse_user_input: WARNING: Setting flag 'verbose' flag overrides flag 'quiet'"
                  exit_code = USER_INPUT_WARNING
                end if

                verbose_flag = .true.
              case (QUIET_SHORT)
                print '(A)', 'DEBUG: q'
                if (quiet_flag) then
                  print '(A)', "parse_user_input: WARNING: Flag 'quiet' was set multiple times"
                  exit_code = USER_INPUT_WARNING
                end if

                if (verbose_flag) then
                  print '(A)', "parse_user_input: WARNING: Flag 'quiet' will be ignored since flag 'verbose' was already set"
                  exit_code = USER_INPUT_WARNING
                end if

                quiet_flag = .true.

              case (INPUT_SHORT)
                print '(A)', 'DEBUG: i'
                if (input_flag) then
                  print '(A)', 'parse_user_input: ERROR: Cannot set input file multiple times'
                  exit_code = USER_INPUT_ERROR
                  return
                end if
                
                input_flag = .true.
                argument = .false.
                params%input_file = DEFAULT_INPUT_FILE

                ! If the -i flag is provided in a list with multiple flags we still allow the user
                ! to set a custom input file by checking the next argument as before
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
                    ! We set the previous argument to INPUT_SHORT = i
                    ! instead of the whole previous arg
                    prev_arg = INPUT_SHORT 
                  end if
                end if
                
              case default
                print '(A,A,A,/,A)', "parse_user_input: ERROR: Flag '", trim(arg(j:j)), & 
                  "'is either not allowed in flag lists or is unrecognized", &
                  "Run 'kettera --help' to see a list of valid flags"
                exit_code = USER_INPUT_ERROR
                return
            end select
          end do 

        end if

        if ((.not. short_flag) .and. (.not. long_flag) .and. (.not. multiple_flags)) then
          print '(A,A,A,/,A)', "parse_user_input: ERROR: Unrecognized argument '", trim(arg), "'", &
            "Run 'kettera --help' to see a list of valid arguments"
          exit_code = USER_INPUT_ERROR
          return
        end if

        if (i == argc .and. argument) then
          print '(A,A,A,/,A)', "parse_user_input: ERROR: Expected argument for flag '", trim(arg), &
            "' but reached end of arguments", "Run 'kettera --help' to see correct usage"
          exit_code = USER_INPUT_ERROR
          return
        end if

      end do

      if (verbose_flag) then
        output_mode = VERBOSE_OUTPUT
      else if (quiet_flag) then
        output_mode = QUIET_OUTPUT
      else
        output_mode = NORMAL_OUTPUT
      end if

      if (simulate_flag) then
        execution_mode = SIMULATE_MODE
      else if (check_flag) then
        execution_mode = CHECK_MODE
      else
        print '(A,/,A)', "parse_user_input: ERROR: Execution mode ('simulate' or 'check') was not provided.", &
          "Run 'kettera --help' to see correct usage"
        exit_code = USER_INPUT_ERROR
        return
      end if

      return
    end subroutine parse_user_input
end module user_io
