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

module iterate

  use kinds
  use exit_codes
  use io_parameters, only: QUIET_OUTPUT
  use sim_parameters, only: SimulationParams, CrankNicolsonArrays, SplitStepArrays, &
    OutputArrays, LoggingArrays 
  use read_write
  use calculate
  use fftw3

  implicit none
  
  public

  contains
  
    ! Crank - Nicolson method iteration

    subroutine iter_cn_method(params, cn_arrays, out_arrays, log_arrays, fin_results, &
        output_mode, exit_code)
      implicit none
      
      type(SimulationParams), intent(in) :: params
      type(CrankNicolsonArrays), intent(inout) :: cn_arrays
      type(OutputArrays), intent(inout) :: out_arrays
      type(LoggingArrays), intent(inout) :: log_arrays
      type(FinalResults), intent(inout) :: fin_results
      integer(label), intent(in) :: output_mode
      integer(excode), intent(out) :: exit_code

      complex(r64), dimension(params%point_count) :: sub, main, sup

      integer(i32) :: N, steps, interval 
      logical :: imag, orto

      integer(excode) :: return_status

      integer(i32) :: t, frame

      real(r64) :: cputime_start, cputime_end

      N = params%point_count
      steps = params%step_count
      interval = params%write_interval
      imag = params%imag_time
      orto = params%ortho

      call init_output_file(params, return_status)
      if (return_status /= SUCCESS) then
         exit_code = return_status
         return
      end if

      call vectorize_tridiagonal_matrix(params, cn_arrays%backward_op, sub, main, sup)

      print '(A)', '----- BEGIN CRANK - NICOLSON METHOD ITERATION -----'
      if (output_mode /= QUIET_OUTPUT) then
        print '(A8,1X,A20,1X,A20)', 'TimeStep', 'Energy', 'NormSquared'
      end if

      frame = 0
      call cpu_time(cputime_start)
      do t = 0, steps
        
        if (modulo(t, interval) == 0) then
          call calculate_output_quantities(cn_arrays%wavefunction, out_arrays)

          call calculate_logging_quantities(frame+1, t, cn_arrays%wavefunction, cn_arrays%hamiltonian, log_arrays)

          if (output_mode /= QUIET_OUTPUT) then
            print '(I8,1X,F20.15,1X,F20.15)', t, log_arrays%energies(frame+1), log_arrays%squared_norms(frame+1)
          end if
          call write_output_file(params, frame, out_arrays, cn_arrays%x_space, cn_arrays%potential, return_status) 
          if (return_status /= SUCCESS) then
            exit_code = return_status
            return
          end if
          
          frame = frame + 1
        end if

        cn_arrays%wavefunction = matmul(cn_arrays%forward_op, cn_arrays%wavefunction)

        cn_arrays%wavefunction = solve_tridiagonal_system(params, sub, main, sup, cn_arrays%wavefunction)

        if (imag) then
          call normalize(cn_arrays%wavefunction)
        end if

        if (orto) then
          call orthogonalize(cn_arrays%wavefunction, cn_arrays%orthogonal)
          ! We only normalize if not iterating in imaginary time since the imaginary
          ! time state was normalized just before orthogonalization and so we don't want
          ! to waste time doing that again
          if (.not. imag) then
            call normalize(cn_arrays%wavefunction)
          end if
        end if

      end do
      call cpu_time(cputime_end)

      print '(a)', '----- END CRANK - NICOLSON METHOD ITERATION -----'

      call close_output_file()

      call calculate_result_quantities(params%frame_count, log_arrays, cputime_start, cputime_end, fin_results)

      call write_log_file(params, log_arrays, fin_results, return_status)
      if (return_status /= SUCCESS) then
        exit_code = return_status
        return
      end if

      exit_code = SUCCESS
      return
    end subroutine iter_cn_method

    ! Split operator method iteration

    subroutine iter_ss_method(params, arrays, output_mode, exit_code)
      implicit none
      
      type(SimulationParams), intent(in) :: params
      type(SplitStepArrays), intent(inout) :: arrays
      integer(label), intent(in) :: output_mode
      integer(excode), intent(out) :: exit_code

      complex(C_DOUBLE_COMPLEX), pointer, dimension(:) :: ft_input, ft_output
      type(C_PTR) :: input_holder, output_holder
      type(C_PTR) :: forward_plan, backward_plan
      real(r64) :: energy, norm_squared

      integer(i32) :: N, steps, interval 
      logical :: imag, orto

      integer(excode) :: return_status

      integer(i32) :: t, frame

      N = params%point_count
      steps = params%step_count
      interval = params%write_interval
      imag = params%imag_time
      orto = params%ortho

      call init_output_file(params, return_status)
      if (return_status /= SUCCESS) then
         exit_code = return_status
         return
      end if

      call associate_holders(int(N, C_SIZE_T), ft_input, ft_output, input_holder, output_holder)

      call get_dft_plans(int(N, C_INT), ft_input, ft_output, forward_plan, backward_plan)

      print '(a)', '----- BEGIN SPLIT STEP METHOD ITERATION -----'
      if (output_mode /= QUIET_OUTPUT) then
        print '(a)', 'TimeStep Energy               NormSquared'
      end if

      frame = 0
      do t = 0, steps
        
        if (modulo(t, interval) == 0) then

          if (output_mode /= QUIET_OUTPUT) then
            print '(I8,1X,F20.15,1X,F17.15)', t, energy, norm_squared
          end if
          if (return_status /= SUCCESS) then
            call destroy_dft_plans(forward_plan, backward_plan)
            call free_holders(input_holder, output_holder)      

            exit_code = return_status
            return
          end if
          
          frame = frame + 1
        end if

        ft_input = arrays%x_space_op * arrays%wavefunction

        call fftw_execute_dft(forward_plan, ft_input, ft_output)

        ft_input = arrays%p_space_op * ft_output

        call fftw_execute_dft(backward_plan, ft_input, ft_output)

        arrays%wavefunction = arrays%x_space_op * arrays%wavefunction

        if (imag) then
          call normalize(arrays%wavefunction)
        end if

        if (orto) then
          call orthogonalize(arrays%wavefunction, arrays%orthogonal)
          call normalize(arrays%wavefunction)
        end if

      end do

      call close_output_file()

      !TODO: Add a call to 'write_log_file(params, log_data)' once that's implemented

      print '(a)', '----- END SPLIT STEP METHOD ITERATION -----'

      call destroy_dft_plans(forward_plan, backward_plan)
      call free_holders(input_holder, output_holder)      

      exit_code = SUCCESS
      return
    end subroutine iter_ss_method

end module iterate
