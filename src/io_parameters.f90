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

module io_parameters
  use kinds

  implicit none

  ! ----- General numerical i/o parameters -----
  integer(label), parameter, public :: MAX_STR_LEN = 200_label
  integer(label), parameter, public :: PARAM_STR_LEN = 20_label
  integer(label), parameter, public :: PARAM_FILE_UNIT = 10_label
  integer(label), parameter, public :: INPUT_FILE_UNIT = 20_label
  integer(label), parameter, public :: OUTPUT_FILE_UNIT = 30_label
  integer(label), parameter, public :: LOG_FILE_UNIT = 40_label

  ! ----- General format strings -----
  character(*), parameter, public :: SINGLE_DIGIT_DOUBLE_FMT = 'F17.15'
  character(*), parameter, public :: SINGLE_FMT = 'ES13.6E2'
  character(*), parameter, public :: DOUBLE_FMT = 'ES23.15E3'
  character(*), parameter, public :: EXTENDED_FMT = 'ES27.18E4'
  character(*), parameter, public :: QUADRUPLE_FMT = 'ES42.33E4'
  character(*), parameter, public :: COMPLEX_FMT = '('//DOUBLE_FMT//',SP,'//DOUBLE_FMT//',"*i")'

  ! ----- Default file names -----
  character(*), parameter, public :: DEFAULT_PARAM_FILE = 'params.in'
  character(*), parameter, public :: DEFAULT_INPUT_FILE = 'wave.in'
  character(*), parameter, public :: DEFAULT_OUTPUT_FILE = 'wave.out'
  character(*), parameter, public :: DEFAULT_LOG_FILE = 'wave.log'

  ! ----- Parameter file definitions -----

  ! Integer labels for parameter types
  integer(label), parameter, public :: NO_TYPE = -1_label
  integer(label), parameter, public :: INT_TYPE = 100_label
  integer(label), parameter, public :: FLOAT_TYPE = 110_label
  integer(label), parameter, public :: BOOL_TYPE = 120_label
  integer(label), parameter, public :: STR_TYPE = 130_label

  ! User given parameter names
  character(*), parameter, public :: PARAM_METHOD = 'Method'
  character(*), parameter, public :: PARAM_IMAG_TIME = 'ImagTime'
  character(*), parameter, public :: PARAM_NORMAL = 'Normal'
  character(*), parameter, public :: PARAM_ORTHO = 'Ortho'
  character(*), parameter, public :: PARAM_UNIT_BOUNDS  = 'UnitBounds'
  character(*), parameter, public :: PARAM_STEP_COUNT = 'StepCount'
  character(*), parameter, public :: PARAM_WRITE_INTERVAL = 'WriteInterval'
  character(*), parameter, public :: PARAM_DELTA_T = 'DeltaT'
  character(*), parameter, public :: PARAM_POINT_COUNT = 'PointCount'
  character(*), parameter, public :: PARAM_X_MAX = 'XMax'
  character(*), parameter, public :: PARAM_WAVE_TYPE = 'WaveType'
  character(*), parameter, public :: PARAM_MASS = 'Mass'
  character(*), parameter, public :: PARAM_CHARGE = 'Charge'
  character(*), parameter, public :: PARAM_MOMENTUM = 'Momentum'
  character(*), parameter, public :: PARAM_WAVE_OFFSET = 'WaveOffset'
  character(*), parameter, public :: PARAM_WAVE_WIDTH = 'WaveWidth'
  character(*), parameter, public :: PARAM_POT_TYPE = 'PotType'
  character(*), parameter, public :: PARAM_POT_OFFSET = 'PotOffset'
  character(*), parameter, public :: PARAM_POT_WIDTH = 'PotWidth'
  character(*), parameter, public :: PARAM_POT_STRENGTH = 'PotStrength'

  ! Parameter file syntax keywords
  character(*), parameter, public :: PARAM_COMMENT = '#'
  character(*), parameter, public :: PARAM_NAME_SEPARATOR = '::'
  character(*), parameter, public :: PARAM_VALUE_SEPARATOR = ':='

  character(*), parameter, public :: PARAM_INT = 'int'
  character(*), parameter, public :: PARAM_FLOAT = 'float'
  character(*), parameter, public :: PARAM_BOOL = 'bool'
  character(*), parameter, public :: PARAM_STR = 'str'

  character(*), parameter, public :: PARAM_DEF = 'DEF'

  ! Accepted values for parameters that expect string values 
  ! Iteration methods
  character(*), parameter, public :: PARAM_CN = 'crank-nicolson'
  character(*), parameter, public :: PARAM_SS = 'split-step'
  ! Wavefunction types
  character(*), parameter, public :: PARAM_GAUSS = 'gaussian'
  character(*), parameter, public :: PARAM_SINC = 'sinc'
  ! Potential types
  character(*), parameter, public :: PARAM_ZERO = 'zero'
  character(*), parameter, public :: PARAM_BOX = 'box'
  character(*), parameter, public :: PARAM_BARRIER = 'barrier'
  character(*), parameter, public :: PARAM_WELL = 'well'
  character(*), parameter, public :: PARAM_DWELL = 'double-well'
  character(*), parameter, public :: PARAM_HARMONIC = 'harmonic'
  character(*), parameter, public :: PARAM_LINEAR = 'linear'
  character(*), parameter, public :: PARAM_HLINEAR = 'half-linear'
  character(*), parameter, public :: PARAM_ABS = 'abs'
  character(*), parameter, public :: PARAM_LOG = 'log'
  character(*), parameter, public :: PARAM_COSH = 'cosh'

  ! ----- Input file definitions -----
  !character(*), parameter, public :: INPUT_POINT_COUNT = 'PointCount'
  !character(*), parameter, public :: INPUT_X_MAX = 'XMax'
  character(*), parameter, public :: INPUT_HEADER = 'Real;Imag'

  character(*), parameter, public :: INPUT_COMMENT = '#'
  character(*), parameter, public :: INPUT_FIELD_SEPARATOR = ';'
  character(*), parameter, public :: INPUT_VALUE_SEPARATOR = '='

  ! ----- Output file strings and format strings -----
  character(*), parameter, public :: OUTPUT_TITLE_1 = &
    '# This file was generated by kettera and contains the output data from a simulation'
  character(*), parameter, public :: OUTPUT_TITLE_2 = &
    '# as frames that represent the state of the system at a given time. The frames can be'
  character(*), parameter, public :: OUTPUT_TITLE_3 = &
    "# animated or a single frame plotted using the program 'ketanim' supplied with kettera"
  character(*), parameter, public :: OUTPUT_HEADER_FMT = '("FrameCount=",I0,";Method=",I0,";PointCount=",I0,";XMax=",F0.2,";PotType=",I0,";ImagTime=",I0)'
  character(*), parameter, public :: OUTPUT_DATA_FIELDS = 'XCoord;Real;Imag;Amplitude;Probability;Potential'
  character(*), parameter, public :: OUTPUT_FRAME_FMT = '("Frame=",I0)'
  character(*), parameter, public :: OUTPUT_DATA_FMT = '(5('//DOUBLE_FMT//',";"),'//DOUBLE_FMT//')'

  ! ----- Log file strings and format strings -----
  character(*), parameter, public :: LOG_TITLE_1 = &
    '# This file was generated by kettera and contains the parameters used to run a simulation'
  character(*), parameter, public :: LOG_TITLE_2 = &
    '# as well as the calculated physical quantities for each time step. The quantities can be'
  character(*), parameter, public :: LOG_TITLE_3 = &
    "# plotted using the program 'ketplot' supplied with kettera"
  character(*), parameter, public :: LOG_FILES_TITLE = '# FILES USED/GENERATED:'
  character(*), parameter, public :: LOG_PARAMS_TITLE = '# SIMULATION PARAMETERS USED/GENERATED:'
  character(*), parameter, public :: LOG_RESULTS_TITLE = '# SIMULATION RESULTS:'
  character(*), parameter, public :: LOG_QUANTITIES_TITLE = '# CALCULATED PHYSICAL QUANTITIES:'
  character(*), parameter, public :: LOG_HEADER_FMT = '("FrameCount=",I0,";Method=",I0,";PotType=",I0,";ImagTime=",I0)'
  character(*), parameter, public :: LOG_DATA_FIELDS = 'TimeStep;Energy;NormSquared'
  character(*), parameter, public :: LOG_DATA_FMT = '(I8,";",'//DOUBLE_FMT//',";",'//DOUBLE_FMT//')'


end module io_parameters
