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

module initialize
  use iso_fortran_env, only: OUTPUT_UNIT
  use kinds
  use exit_codes
  use sim_parameters
  use calculate

  implicit none

  private

  integer(label), parameter :: CN_FORWARD = +1_label
  integer(label), parameter :: CN_BACKWARD = -1_label
  integer(label), parameter :: SS_POSITION = 1_label
  integer(label), parameter :: SS_MOMENTUM = 2_label

  public init_params, init_cn_iteration, init_ss_iteration

  ! All private module elements have a leading underscore eg. init_delta_x
  ! whereas public module elements don't eg. init_cn_iteration

  contains

    subroutine init_params(params)
      ! Subroutine to initialize the simulation parameters with default values
      implicit none

      type(SimulationParams), intent(inout) :: params
      
      params%iter_method = Method%CRANK_NICOLSON
      params%imag_time = .false.
      params%normalize = .true.
      params%unit_bounds = .true.

      params%step_count = 10000_i32
      params%write_interval = 100_i32
      params%delta_t = 1e-3_r64
      params%point_count = 1024_i32
      params%x_max = 10.0_r64

      params%wave_type = WaveType%GAUSSIAN
      params%mass = 1.0_r64
      params%charge = 1.0_r64
      params%momentum = 0.0_r64
      params%wave_offset = -3.0_r64
      params%wave_width = 1.0_r64

      params%pot_type = PotType%HARMONIC
      params%pot_offset = 0.0_r64
      params%pot_width = 1.0_r64
      params%pot_strength = 0.5_r64
      
      return
    end subroutine init_params

    subroutine init_delta_x(params, exit_code)
      ! Determine discretization of x-space
      implicit none

      type(SimulationParams), intent(inout) :: params
      integer(excode), intent(out) :: exit_code

      real(r64) :: dx
      real(r64) :: one_over_point_count

      one_over_point_count = 1.0_r64 / params%point_count
      dx = 2.0_r64 * params%x_max * one_over_point_count
      if (dx > DELTA_X_MAX) then
        print '(a,1x,f10.6,1x,a,1x,f10.6)', 'init_delta_x: ERROR: Calculated value for dx was', &
          dx, 'which is too large! dx must be at most', DELTA_X_MAX
        exit_code = DERIVED_PARAMETER_ERROR
        return
      end if

      params%delta_x = dx

      exit_code = SUCCESS
      return
    end subroutine init_delta_x

    subroutine init_delta_p(params, exit_code)
      ! Determine discretization of p-space
      implicit none

      type(SimulationParams), intent(inout) :: params
      integer(excode), intent(out) :: exit_code

      real(r64) :: dx, dp

      dx = params%delta_x

      if (abs(dx) <= EPS) then
        print '(a)', 'init_delta_p: ERROR: dx has not yet been set meaning that dp cannot be set!'
        exit_code = DERIVED_PARAMETER_ERROR
        return
      end if

      dp = PI_CONST / dx
      params%delta_p = dp

      exit_code = SUCCESS
      return
    end subroutine init_delta_p

    function init_x_space(params, exit_code) result(x_space)
      ! Construct the real space the wavefunction lives in
      implicit none

      type(SimulationParams), intent(in) :: params
      integer(excode), intent(out) :: exit_code
      real(r64) :: x_space(params%point_count)
      real(r64) :: bound, dx
      integer(i32) :: N, i

      bound = params%x_max
      dx = params%delta_x

      if (abs(dx) <= EPS) then
        print '(a)', 'init_x_space: ERROR: dx has not yet been set meaning that x-space cannot be defined!'
        exit_code = DERIVED_PARAMETER_ERROR
        return
      end if

      if (bound < X_BOUND_MIN) then
        print '(a,1x,f10.6,1x,a,1x,f3.1)', 'init_x_space: ERROR: Given value for x bound', &
          bound, 'is too small! X bound must be at least', X_BOUND_MIN
        exit_code = GIVEN_PARAMETER_ERROR
        x_space = 0.0_r64
        return
      end if

      N = params%point_count

      do i = 1, N
        x_space(i) = -bound + 0.5_r64*dx + (i-1)*dx
      end do

      exit_code = SUCCESS
      
      return
    end function init_x_space

    function init_p_space(params, exit_code) result(p_space)
      ! Construct the momentum space the wavefunction lives in
      implicit none

      type(SimulationParams), intent(in) :: params
      integer(excode), intent(out) :: exit_code
      real(r64) :: p_space(params%point_count)
      real(r64) :: dp
      integer(i32) :: N, i

      dp = params%delta_p

      if (abs(dp) <= EPS) then
        print '(a)', 'init_p_space: ERROR: dp has not yet been set meaning that p-space cannot be defined!'
        exit_code = DERIVED_PARAMETER_ERROR
        return
      end if

      N = params%point_count

      if (mod(N, 2) == 0) then
        do i = 0, N/2-1
          p_space(i+1) = i
        end do
        do i = N/2, N-1
          p_space(i+1) = -(N - i)
        end do
      else
        do i = 0, (N-1)/2
          p_space(i+1) = i
        end do
        do i = (N+1)/2, N-1
          p_space(i+1) = -(N - i)
        end do
      end if
      
      p_space = dp * p_space
      return
    end function init_p_space

    function init_wavefunction(params, x_space, exit_code) result(wavefunction)
      implicit none

      type(SimulationParams), intent(in) :: params
      real(r64), intent(in) :: x_space(:)
      integer(excode), intent(out) :: exit_code
      complex(r64) :: wavefunction(params%point_count)

      real(r64) :: p, w, off
      complex(r64) :: momentum_term(params%point_count)
      real(r64) :: factor
      integer(label) :: wave

      p = params%momentum
      w = params%wave_width
      off = params%wave_offset
      wave = params%wave_type
  
      momentum_term = exp(IMAG_UNIT * p * x_space)

      if (wave == WaveType%GAUSSIAN) then

        factor = 1.0_r64 / sqrt(2.0_r64 * PI_CONST * w**2)
        wavefunction = factor * momentum_term * exp(-(x_space - off)**2 / (2 * w**2))

      else if (wave == WaveType%SINC) then

        factor = PI_CONST / w
        wavefunction = momentum_term * sinc(factor * (x_space - off))

      else

        print '(a,1x,I0,1x,a)', 'init_wavefunction: ERROR: Given value for wavefunction type', &
          params%wave_type, 'is not recognized!'
        wavefunction = 0.0_r64 * x_space
        exit_code = GIVEN_PARAMETER_ERROR
        return

      end if

      if (params%normalize) then
        call normalize(wavefunction)
      end if

      exit_code = SUCCESS
      return
    end function init_wavefunction


    function init_potential(params, x_space, exit_code) result(potential)
      implicit none

      type(SimulationParams), intent(in) :: params
      real(r64), intent(in) :: x_space(:)
      integer(excode), intent(out) :: exit_code
      real(r64) :: potential(params%point_count)

      real(r64) :: m, q, off, w, s, half_w
      integer(label) :: pot

      m = params%mass
      q = params%charge
      off = params%pot_offset
      w = params%pot_width
      s = params%pot_strength
      pot = params%pot_type

      half_w = w / 2.0_r64

      if (pot == PotType%ZERO) then

        potential = 0.0_r64 * x_space

      else if (pot == PotType%BOX) then

        where ((x_space <= -off) .or. (x_space >= off))
          potential = s
        elsewhere
          potential = 0.0_r64
        end where

      else if (pot == PotType%BARRIER) then

        where ((x_space >= off) .and. (x_space <= off + w))
          potential = s
        elsewhere
          potential = 0.0_r64
        end where

      else if (pot == PotType%WELL) then

        where ((x_space >= off - half_w) .and. (x_space <= off + half_w))
          potential = -abs(s) ! Making sure that the well has negative potential
        elsewhere
          potential = 0.0_r64
        end where

      else if (pot == PotType%DOUBLE_WELL) then

        where (((x_space >= -off - half_w) .and. (x_space <= -off + half_w)) &
            .or. ((x_space >= off - half_w) .and. (x_space <= off + half_w)))
          potential = -abs(s) ! Making sure that the wells have negative potential
        elsewhere
          potential = 0.0_r64
        end where

      else if (pot == PotType%HARMONIC) then

        potential = 0.5_r64 * m * s**2 * (x_space - off)**2 / w

      else if (pot == PotType%LINEAR) then

        potential = m * s * (x_space - off)

      else if (pot == PotType%HALF_LINEAR) then

        where (x_space < off)
          potential = 0.0_r64
        elsewhere
          potential = m * s * (x_space - off)
        end where

      else if (pot == PotType%ABSOLUTE_VALUE) then

        potential = m * s * abs(x_space - off) / w

      else if (pot == PotType%LOGARITHMIC) then

        potential = q * s * log(x_space - off) / w

      else if (pot == PotType%HYPERBOLIC_COSINE) then

        potential = q * s * cosh(x_space - off) / w

      else

        print '(a,1x,I0,1x,a)', 'init_potential: ERROR: Given value for potential type', &
          params%pot_type, 'is not recognized!'
        potential = 0.0_r64 * x_space
        exit_code = GIVEN_PARAMETER_ERROR
        return

      end if

      exit_code = SUCCESS
      return
    end function init_potential

    function init_ss_kinetic_operator(params, p_space) result(ss_kinetic)
      implicit none
      
      type(SimulationParams), intent(in) :: params
      real(r64), intent(in) :: p_space(:)
      real(r64) :: ss_kinetic(params%point_count)

      real(r64) :: m
      
      m = params%mass

      ss_kinetic = (1.0_r64 / (2.0_r64 * m)) * p_space**2
      
      return
    end function init_ss_kinetic_operator

    function init_hamiltonian(params, potential) result(hamiltonian)
      implicit none
      
      type(SimulationParams), intent(in) :: params
      real(r64), intent(in) :: potential(:)
      real(r64) :: hamiltonian(params%point_count,params%point_count)
      
      real(r64) :: dx, m, kinetic
      integer(i32) :: N
      logical :: unitb

      integer :: i

      dx = params%delta_x
      m = params%mass
      kinetic = 1.0_r64 / (2.0_r64 * m * dx**2)
      N = params%point_count
      unitb = params%unit_bounds

      hamiltonian = 0.0_r64

      do i = 2, N-1
        hamiltonian(i,i-1) = -kinetic
        hamiltonian(i,i+1) = -kinetic
        hamiltonian(i,i) = 2.0_r64 * kinetic + potential(i)
      end do

      ! The bounds are determined separately
      if (unitb) then
        !hamiltonian(1,1) = 0.0_r64
        !hamiltonian(N,N) = 0.0_r64
        hamiltonian(1,1) = 1.0_r64
        hamiltonian(N,N) = 1.0_r64
      else
        hamiltonian(1,1) = 2.0_r64 * kinetic + potential(1)
        hamiltonian(N,N) = 2.0_r64 * kinetic + potential(N)
      end if

      return
    end function init_hamiltonian

    function init_cn_evolution_operator(params, hamiltonian, direction, exit_code) result(cn_operator) 
      implicit none
      
      type(SimulationParams), intent(in) :: params
      real(r64), intent(in) :: hamiltonian(:,:)
      integer(label), intent(in) :: direction
      integer(excode), intent(out) :: exit_code

      complex(r64) :: cn_operator(params%point_count,params%point_count)
      real(r64) :: identity(params%point_count,params%point_count)

      real(r64) :: dt, halfstep
      integer(i32) :: N
      logical :: imag_t

      integer(i32) :: i,j

      N = params%point_count

      ! Creating an NxN identity operator
      do j = 1, N
        do i =  1, N
          if (i == j) then
            identity(i,j) = 1.0_r64
          else
            identity(i,j) = 0.0_r64
          end if
        end do
      end do

      dt = params%delta_t
      halfstep = dt / 2.0_r64
      N = params%point_count
      imag_t = params%imag_time

      if (direction == CN_FORWARD) then
        if (imag_t) then
          cn_operator = identity - (halfstep * hamiltonian)
        else
          cn_operator = identity - (IMAG_UNIT * halfstep * hamiltonian)
        end if
      else if (direction == CN_BACKWARD) then
        if (imag_t) then
          cn_operator = identity + (halfstep * hamiltonian)
        else
          cn_operator = identity + (IMAG_UNIT * halfstep * hamiltonian)
        end if
      else
        print '(a,1x,I0,1x,a)', 'init_cn_evolution_operator: ERROR: Given value for direction', &
          direction, 'is not recognized!'
        cn_operator = 0
        exit_code = DERIVED_PARAMETER_ERROR
        return
      end if

      if (N <= 12) then
        print '(A)', 'init_cn_evolution_operator: DEBUG: cn_operator:'
        do j = 1, N
          do i = 1, N
            write(OUTPUT_UNIT, fmt='(F8.5,SP,F8.5,"*i",2X)', advance='no') cn_operator(j,i)
          end do
          print *
        end do
      end if

      exit_code = SUCCESS

      return
    end function init_cn_evolution_operator

    function init_ss_evolution_operator(params, values, which_type, exit_code) result(ss_operator)
      implicit none
      
      type(SimulationParams), intent(in) :: params
      real(r64), intent(in) :: values(:) ! values is either the kinetic operator or the potential operator
      integer(label), intent(in) :: which_type
      integer(excode), intent(out) :: exit_code

      complex(r64) :: ss_operator(params%point_count)

      real(r64) :: dt, halfstep
      integer(i32) :: N
      logical :: imag_t

      dt = params%delta_t
      halfstep = dt / 2.0_r64
      N = params%point_count
      imag_t = params%imag_time

      if (which_type == SS_POSITION) then
        if (imag_t) then
          ss_operator = exp(-halfstep * values) ! position space case -> halfstep
        else
          ss_operator = exp(-IMAG_UNIT * halfstep * values) ! position space case -> halfstep
        end if
      else if (which_type == SS_MOMENTUM) then
        if (imag_t) then
          ss_operator = exp(-dt * values) ! momentum space case -> dt
        else
          ss_operator = exp(-IMAG_UNIT * dt * values) ! momentum space case -> dt
        end if
      else
        print '(a,1x,I0,1x,a)', 'init_ss_evolution_operator: ERROR: Given value for which_type', &
          which_type, 'is not recognized!'
        ss_operator = 0
        exit_code = DERIVED_PARAMETER_ERROR
        return
      end if

      exit_code = SUCCESS
      return
    end function init_ss_evolution_operator

    subroutine init_cn_iteration(params, arrays, extern_input, exit_code)
      ! This subroutine is used to initialize the variables and
      ! parameters required to perform Crank - Nicolson method iteration
      implicit none
      
      type(SimulationParams), intent(inout) :: params
      type(CrankNicolsonArrays), intent(inout) :: arrays
      logical, intent(in) :: extern_input
      integer(excode), intent(out) :: exit_code 

      call init_delta_x(params, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%x_space = init_x_space(params, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      if (.not. extern_input) then
        arrays%wavefunction = init_wavefunction(params, arrays%x_space, exit_code)
        if (exit_code /= SUCCESS) then
          return
        end if
      end if
      arrays%potential =  init_potential(params, arrays%x_space, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%hamiltonian = init_hamiltonian(params, arrays%potential)
      arrays%forward_op = init_cn_evolution_operator(params, arrays%hamiltonian, CN_FORWARD, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%backward_op = init_cn_evolution_operator(params, arrays%hamiltonian, CN_BACKWARD, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if

      exit_code = SUCCESS
      return
    end subroutine init_cn_iteration

    subroutine init_ss_iteration(params, arrays, extern_input, exit_code)
      ! This subroutine is used to initialize the variables and
      ! parameters required to perform split step method iteration
      implicit none
      
      type(SimulationParams), intent(inout) :: params
      type(SplitStepArrays), intent(inout) :: arrays
      logical, intent(in) :: extern_input
      integer(excode), intent(out) :: exit_code 

      call init_delta_x(params, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%x_space = init_x_space(params, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      call init_delta_p(params, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%p_space = init_p_space(params, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      if (.not. extern_input) then
        arrays%wavefunction = init_wavefunction(params, arrays%x_space, exit_code)
        if (exit_code /= SUCCESS) then
          return
        end if
      end if
      arrays%potential =  init_potential(params, arrays%x_space, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%kinetic = init_ss_kinetic_operator(params, arrays%p_space)
      arrays%hamiltonian = init_hamiltonian(params, arrays%potential)
      arrays%x_space_op = init_ss_evolution_operator(params, arrays%potential, SS_POSITION, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if
      arrays%p_space_op = init_ss_evolution_operator(params, arrays%kinetic, SS_MOMENTUM, exit_code)
      if (exit_code /= SUCCESS) then
        return
      end if

      exit_code = SUCCESS
      return
      
    end subroutine init_ss_iteration

end module initialize

