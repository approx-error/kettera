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

module alloc_dealloc
  use kinds
  use exit_codes
  use sim_parameters
  use io_parameters

  implicit none
  
  public

  contains

    subroutine alloc_cn_arrays(N, arrays, exit_code)

      implicit none
      
      integer(i32), intent(in) :: N
      type(CrankNicolsonArrays), intent(out) :: arrays 
      integer(excode), intent(out) :: exit_code

      integer(i32) :: error_status
      character(len=MAX_STR_LEN) :: error_message

      allocate(arrays%x_space(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate x-space. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%wavefunction(N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate wavefunction. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%orthogonal(N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate orthogonal wavefunction. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%potential(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate potential. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%hamiltonian(N,N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate hamiltonian. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%forward_op(N,N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate forward evolution operator. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%backward_op(N,N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_cn_arrays: ERROR: Could not allocate backward evolution operator. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if

      exit_code = SUCCESS
      return
    end subroutine alloc_cn_arrays

    subroutine dealloc_cn_arrays(arrays, exit_code)
      implicit none

      type(CrankNicolsonArrays), intent(inout) :: arrays
      integer(excode), intent(out) :: exit_code

      exit_code = SUCCESS

      if (allocated(arrays%x_space)) then
        deallocate(arrays%x_space)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: x-space was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%wavefunction)) then
        deallocate(arrays%wavefunction)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: wavefunction was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%orthogonal)) then
        deallocate(arrays%orthogonal)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: orthogonal wavefunction was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if ((allocated(arrays%potential))) then
        deallocate(arrays%potential)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: potential was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%hamiltonian)) then
        deallocate(arrays%hamiltonian)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: hamiltonian was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%forward_op)) then
        deallocate(arrays%forward_op)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: forward evolution operator was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%backward_op)) then
        deallocate(arrays%backward_op)
      else
        print '(a)', 'dealloc_cn_arrays: WARNING: backward evolution operator was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
  
      return
    end subroutine dealloc_cn_arrays


    subroutine alloc_ss_arrays(N, arrays, exit_code)

      implicit none
      
      integer(i32), intent(in) :: N
      type(SplitStepArrays), intent(out) :: arrays
      integer(excode), intent(out) :: exit_code

      integer(i32) :: error_status
      character(len=MAX_STR_LEN) :: error_message

      allocate(arrays%x_space(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate x-space. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%p_space(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate p-space. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%wavefunction(N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate wavefunction. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%orthogonal(N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate orthogonalwavefunction. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%potential(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate potential. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%kinetic(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate kinetic energy operator. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%hamiltonian(N,N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate hamiltonian. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%x_space_op(N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate x-space evolution operator. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%p_space_op(N), source=(0.0_r64,0.0_r64), stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_ss_arrays: ERROR: Could not allocate p-space evolution operator. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if

      exit_code = SUCCESS
      return
    end subroutine alloc_ss_arrays

    subroutine dealloc_ss_arrays(arrays, exit_code)
      implicit none

      type(SplitStepArrays), intent(inout) :: arrays
      integer(excode), intent(out) :: exit_code

      exit_code = SUCCESS

      if (allocated(arrays%x_space)) then
        deallocate(arrays%x_space)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: x-space was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%p_space)) then
        deallocate(arrays%p_space)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: p-space was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%wavefunction)) then
        deallocate(arrays%wavefunction)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: wavefunction was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%orthogonal)) then
        deallocate(arrays%orthogonal)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: orthogonal wavefunction was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%potential)) then
        deallocate(arrays%potential)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: potential was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%kinetic)) then
        deallocate(arrays%kinetic)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: kinetic energy operator was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%hamiltonian)) then
        deallocate(arrays%hamiltonian)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: hamiltonian was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%x_space_op)) then
        deallocate(arrays%x_space_op)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: x-space evolution operator was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%p_space_op)) then
        deallocate(arrays%p_space_op)
      else
        print '(a)', 'dealloc_ss_arrays: WARNING: p-space evolution operator was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
  
      return
    end subroutine dealloc_ss_arrays

    subroutine alloc_output_arrays(N, arrays, exit_code)
      implicit none

      integer(i32), intent(in) :: N
      type(OutputArrays), intent(out) :: arrays
      integer(excode), intent(out) :: exit_code

      integer(i32) :: error_status
      character(len=MAX_STR_LEN) :: error_message

      allocate(arrays%real_part(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_output_arrays: ERROR: Could not allocate real part array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%amplitude(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_output_arrays: ERROR: Could not allocate amplitude array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%imag_part(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_output_arrays: ERROR: Could not allocate imaginary part array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%density(N), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_output_arrays: ERROR: Could not allocate density array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      
      exit_code = SUCCESS
      return
    end subroutine alloc_output_arrays

    subroutine dealloc_output_arrays(arrays, exit_code)
      implicit none
      
      type(OutputArrays), intent(inout) :: arrays
      integer(excode), intent(out) :: exit_code

      exit_code = SUCCESS

      if (allocated(arrays%real_part)) then
        deallocate(arrays%real_part)
      else
        print '(a)', 'dealloc_output_arrays: WARNING: real part array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%imag_part)) then
        deallocate(arrays%imag_part)
      else
        print '(a)', 'dealloc_output_arrays: WARNING: imaginary part array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%amplitude)) then
        deallocate(arrays%amplitude)
      else
        print '(a)', 'dealloc_output_arrays: WARNING: amplitude array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%density)) then
        deallocate(arrays%density)
      else
        print '(a)', 'dealloc_output_arrays: WARNING: density array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if

      return
    end subroutine dealloc_output_arrays

    subroutine alloc_log_arrays(frames, arrays, exit_code)
      implicit none
      
      integer(i32), intent(in) :: frames
      type(LoggingArrays), intent(out) :: arrays
      integer(excode), intent(out) :: exit_code

      integer(i32) :: error_status
      character(len=MAX_STR_LEN) :: error_message

      allocate(arrays%timestamps(frames), source=0_i32, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_log_arrays: ERROR: Could not allocate timestamp array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%energies(frames), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_log_arrays: ERROR: Could not allocate energy array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      allocate(arrays%squared_norms(frames), source=0.0_r64, stat=error_status, errmsg=error_message)
      if (error_status /= 0) then
        print '(a,1/,a)', 'alloc_log_arrays: ERROR: Could not allocate norm squared array. Message: ', trim(error_message)
        exit_code = ALLOCATION_ERROR
        return
      end if
      
      exit_code = SUCCESS
      return
    end subroutine alloc_log_arrays

    subroutine dealloc_log_arrays(arrays, exit_code)
      implicit none
      
      type(LoggingArrays), intent(inout) :: arrays
      integer(excode), intent(out) :: exit_code

      exit_code = SUCCESS

      if (allocated(arrays%timestamps)) then
        deallocate(arrays%timestamps)
      else
        print '(a)', 'dealloc_log_arrays: WARNING: timestamp array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%energies)) then
        deallocate(arrays%energies)
      else
        print '(a)', 'dealloc_log_arrays: WARNING: energy array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if
      if (allocated(arrays%squared_norms)) then
        deallocate(arrays%squared_norms)
      else
        print '(a)', 'dealloc_log_arrays: WARNING: norm squared array was already deallocated!'
        exit_code = DEALLOCATION_WARNING
      end if

      return
    end subroutine dealloc_log_arrays
  
end module alloc_dealloc
  
