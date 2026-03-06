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

  implicit none

  private

  contains

    subroutine show_usage
      implicit none
      
      print '(A)', 'Usage: kettera [-p <param-file>] [-i <input-file>] [-o <output-file>] [-l <log-file>]'
      print '(A)', "Run 'kettera --help' for more information"

    end subroutine show_usage

    subroutine show_help
      implicit none
      
      print '(A)', 'Usage: kettera [-v | -q] [-p <param-file>] [-i [input-file]] [-o <output-file>] [-l <log-file>]'
      print '(A)', '   or: kettera [--help | --version]'
      print '(A)', 'Simulate a quantum mechanical system using provided parameters and output the results into files'
      print '(/,A)', 'Options surrounded by square brackets [] are optional while options surrounded'
      print '(A)',   'by angle brackets <> are mandatory'
      print '(/,A)', 'If a file is left unspecified, the default filenames will be searched for in the current directory'
      print '(/,A)', 'Options:'
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
  
end module user_io
