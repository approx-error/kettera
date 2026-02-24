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

module calculate

  use kinds
  use sim_parameters
  use fftw3

  implicit none

  public

  contains

    ! ----- General calculation subroutines and functions -----

    elemental function sinc(x) result(res)
      ! Elemental implementation of the sinc function
      !              __
      !             /
      !             |  sin(x) / x,  when x != 0
      ! sinc(x) := <   1         ,  when x == 0
      !             |
      !             \__
      implicit none
      
      real(r64), intent(in) :: x
      real(r64) :: res

      if (abs(x) <= EPS) then
        res = 1.0_r64
      else
        res = sin(x) / x
      end if

      return
    end function sinc

    subroutine normalize(wavefunction)
      ! Normalizes the given wavefunction so that |wavefunction|^2 = 1
      implicit none

      complex(r64), intent(inout) :: wavefunction(:)
      real(r64) :: norm_squared
      complex(r64) :: test ! TODO: remove this later if confirmed that norm_squared has zero complex component

      test = dot_product(wavefunction, wavefunction)
      if (abs(aimag(test)) >= EPS) then
        print '(A,F18.15)', 'normalize: INFO: norm squared has nonzero imaginary component: ', aimag(test)
      endif 
      norm_squared = real(test)
      wavefunction = wavefunction / sqrt(norm_squared)
    end subroutine normalize

    subroutine calculate_output_quantities(wavefunction, out_arrays)
      ! Calculates the amplitude, real part, imaginary part and probability density
      ! of the given wavefunction
      implicit none
      
      complex(r64), intent(in) :: wavefunction(:)
      type(OutputArrays), intent(inout) :: out_arrays

      out_arrays%density = abs(wavefunction ** 2)
      out_arrays%amplitude = sqrt(out_arrays%density)
      out_arrays%real_part = real(wavefunction)
      out_arrays%imag_part = aimag(wavefunction)

      return
    end subroutine calculate_output_quantities

    subroutine calculate_logging_quantities(array_index, time, wavefunction, hamiltonian, log_arrays)
      ! Calculates the energy and norm squared of the given wavefunction and associates
      ! a time stamp with the calculations
      implicit none

      integer(i32), intent(in) :: array_index
      integer(i32), intent(in) :: time
      complex(r64), intent(in) :: wavefunction(:)
      real(r64), intent(in) :: hamiltonian(:,:)
      type(LoggingArrays), intent(inout) :: log_arrays

      real(r64) :: energy, norm_squared
      complex(r64) :: test ! TODO: remove this later if confirmed that energy has zero complex component

      test = dot_product(wavefunction, matmul(hamiltonian, wavefunction)) 
      if (abs(aimag(test)) >= EPS) then
        print '(A,F18.15)', 'calculate_logging_quantities: INFO: energy has nonzero imaginary component: ', aimag(test)
      endif 
      energy = real(test)

      test = dot_product(wavefunction, wavefunction)
      if (abs(aimag(test)) >= EPS) then
        print '(A,F18.15)', 'calculate_logging_quantities: INFO: norm squared has nonzero imaginary component: ', aimag(test)
      endif 
      norm_squared = real(test)

      log_arrays%timestamps(array_index) = time
      log_arrays%energies(array_index) = energy
      log_arrays%squared_norms(array_index) = norm_squared

      return 
    end subroutine calculate_logging_quantities
  
    ! ----- Crank - Nicolson method subroutines and functions -----

    subroutine vectorize_tridiagonal_matrix(params, matrix, sub, main, sup)
      ! Expresses a tridiagonal matrix as the three diagonals that make it up
      implicit none

      type(SimulationParams), intent(in) :: params
      complex(r64), intent(in) :: matrix(:,:)
      complex(r64), dimension(params%point_count), intent(out) :: sub, main, sup

      integer(i32) :: N

      integer :: i

      N = params%point_count

      ! bounds
      sub(1) = 0.0_r64
      sub(N) = matrix(N,N-1)
      main(1) = matrix(1,1)
      main(N) = matrix(N,N)
      sup(1) = matrix(1,2)
      sup(N) = 0.0_r64

      ! other values
      do i = 2, N-1
        sub(i) = matrix(i,i-1)
        main(i) = matrix(i,i)
        sup(i) = matrix(i,i+1)
      end do
    end subroutine vectorize_tridiagonal_matrix

    function solve_tridiagonal_system(params, sub, main, sup, vec) result(solution)
      ! Solves x from a matrix equation of the form Tx = y where T is a tridiagonal matrix
      ! using the Thomas algorithm. The matrix is given as three arrays sub, main and sup
      ! which correspond to the subdiagonal, main diagonal and supradiagonal respectively.
      ! The vector y is given as an array called vec. This implementations keeps the
      ! inputs the same as they will be needed every time this function is called
      ! TODO: Check to see if this can be optimized since every term on the sub- and supradiagonals
      ! is identical
      implicit none
      
      type(SimulationParams), intent(in) :: params
      complex(r64), dimension(:), intent(in) :: sub, main, sup, vec
      complex(r64) :: solution(params%point_count) 

      complex(r64) :: sup_new(params%point_count)
      complex(r64) :: vec_new(params%point_count)
      complex(r64) :: numerator, denominator
     
      integer(i32) :: N

      integer :: i

      N = params%point_count

      ! Initial values
      sup_new(1) = sup(1) / main(1)
      vec_new(1) = vec(1) / main(1)

      ! Solving the rest through recurrence
      do i = 2, N
        numerator = vec(i) - sub(i)*vec_new(i-1)
        denominator = main(i) - sub(i)*sup_new(i-1)    

        sup_new(i) = sup(i) / denominator
        vec_new(i) = numerator / denominator
      end do
      ! The cost of doing the calculations in one for loop
      ! is that new_sup(N) is calculated even though it is not
      ! used when solving

      ! Backsubstitution to get the final solution
      solution(N) = vec_new(N)
      do i = N-1, 1, -1
        solution(i) = vec_new(i) - sup_new(i)*solution(i+1)
      end do
       
      return
    end function solve_tridiagonal_system


    ! ----- Split operator method subroutines and functions -----

    subroutine associate_holders(N, input_array, output_array, input_holder, output_holder)
      ! This subroutine associates to the fortran arrays 'input_array' and 'output_array'
      ! C pointers 'input_holder' and 'output_holder' by using a custom allocator from
      ! FFTW so that memory alignment is guaranteed and thus performance is improved
      implicit none
      
      integer(C_SIZE_T), intent(in) :: N
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:), intent(in) :: input_array, output_array
      type(C_PTR), intent(out) :: input_holder, output_holder

      input_holder = fftw_alloc_complex(N)
      output_holder = fftw_alloc_complex(N)

      call c_f_pointer(input_holder, input_array, [N])
      call c_f_pointer(output_holder, output_array, [N])
      
    end subroutine associate_holders

    subroutine free_holders(input_holder, output_holder)
      ! This subroutine frees the C pointers allocated by 'associate_c_to_f' after
      ! we are done using them
      implicit none
      
      type(C_PTR), intent(in) :: input_holder, output_holder

      call fftw_free(input_holder)
      call fftw_free(output_holder)
    end subroutine free_holders

    subroutine get_dft_plans(N, input_array, output_array, forward_plan, backward_plan)
      ! This subroutine produces FFTW plans for the forward and backward fourier transforms
      implicit none
      
      integer(C_INT), intent(in) :: N
      complex(C_DOUBLE_COMPLEX), pointer, dimension(:), intent(in) :: input_array, output_array
      type(C_PTR), intent(out) :: forward_plan, backward_plan

      forward_plan = fftw_plan_dft_1d(N, input_array, output_array, FFTW_FORWARD, FFTW_ESTIMATE)
      backward_plan = fftw_plan_dft_1d(N, input_array, output_array, FFTW_BACKWARD, FFTW_ESTIMATE)
      
    end subroutine get_dft_plans

    subroutine destroy_dft_plans(forward_plan, backward_plan)
      ! This subroutine destroys the FFTW plans after we are done using them
      implicit none
      
      type(C_PTR), intent(in) :: forward_plan, backward_plan
      
      call fftw_destroy_plan(forward_plan)
      call fftw_destroy_plan(backward_plan)

    end subroutine destroy_dft_plans

    
end module calculate
