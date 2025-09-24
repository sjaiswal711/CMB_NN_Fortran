module constants
    implicit none
    integer, parameter :: long_kind = selected_int_kind(18)

    ! Define the parameters
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), parameter :: theta1 = 7.5 * pi / 180.0d0
    real(kind=8), parameter :: theta2 = 85.0 * pi / 180.0d0
    real(kind=8), parameter :: w1 = 2.0d0 * pi  ! rad/min
    real(kind=8), parameter :: w3 = 0.000011954d0  ! rad/min
   !  real(kind=8), parameter :: w2 = 2.0d0 * w3  ! rad/min
    real(kind=8), parameter :: w2 = 2.0d0 * w1  ! rad/min
 
    integer, parameter :: nside = 1024
    integer(long_kind), parameter :: npix = 12 * nside * nside
 
    real(kind=8), parameter :: scan_time = sqrt(4.0d0 * pi / real(npix,kind=8)) / w1
 
    real(kind=8), parameter :: start_time = 0.0d0
    real(kind=8), parameter :: duration = 1.0d0 * 365.0d0 * 60.0d0 * 24.0d0
 
    integer(long_kind), parameter :: steps = int(duration / scan_time, kind=long_kind) - 1_long_kind
 
    real(kind=8), parameter :: dt = duration / real(steps - 1, kind=8)
 
 end module constants
 
 module subroutines
    use math_library
    use healpix_modules
    use constants
    implicit none
 
 contains
 
    subroutine get_vector(node_id, t, pixel, result_r)

         real(kind=8), intent(in) :: t
         integer, intent(in) :: node_id
         real(kind=8), intent(out) :: result_r(3,1)
         integer, intent(out) :: pixel
   
         real(kind=8) :: cos_theta1, sin_theta1, cos_theta2, sin_theta2
         real(kind=8) :: cos_w1t, sin_w1t, cos_w2t, sin_w2t, cos_w3t, sin_w3t
         real(kind=8), dimension(3, 3) :: A, B, C, result1
         real(kind=8), dimension(3,1) :: D_R, D_S
   
         cos_theta1 = cos(theta1)
         sin_theta1 = sin(theta1)
   
         cos_theta2 = cos(theta2)
         sin_theta2 = sin(theta2)
   
         ! Calculate trigonometric terms
         cos_w1t = cos(w1 * t)
         sin_w1t = sin(w1 * t)
   
         cos_w2t = cos(w2 * t)
         sin_w2t = sin(w2 * t)
   
         cos_w3t = cos(w3 * t)
         sin_w3t = sin(w3 * t)
   
         ! Define the matrices
         A = reshape([ &
            cos_w3t, sin_w3t, 0.0d0, &
            -sin_w3t, cos_w3t, 0.0d0, &
            0.0d0, 0.0d0, 1.0d0 &
            ], [3, 3], order=[2,1])

         B = reshape([ &
            1.0d0, 0.0d0, 0.0d0, &
            0.0d0, cos_w2t, sin_w2t, &
            0.0d0, -sin_w2t, cos_w2t &
            ], [3, 3], order=[2,1])
   
         C = reshape([ &
            cos_theta1, 0.0d0, sin_theta1, &
            0.0d0, 1.0d0, 0.0d0, &
            -sin_theta1, 0.0d0, cos_theta1 &
            ], [3, 3], order=[2,1])
   
         ! Define vectors
         D_R = reshape([cos_theta2, sin_theta2 * cos(w1 * t), sin_theta2 * sin(w1 * t)], [3, 1], order=[2,1])
        
         ! Perform matrix multiplications (optimized)
         result1 = matmul(matmul(A,B),C)
         result_r = matmul(result1, D_R)
   
         call vec2pix_ring(nside, reshape(result_r, [3]), pixel)   
    end subroutine get_vector
 
 end module subroutines
 
 program main
   use constants
   use healpix_modules
   use subroutines
   implicit none

   integer(long_kind) :: i, idx
   integer :: pixel
   real(kind=8) :: t
   real(kind=8) :: s_time, e_time, elapsed_time
   real(kind=8), allocatable :: results(:,:)  ! columns: pixel_id, hit_count
   real(kind=8) :: result_r(3)
   integer :: unit
   integer :: progress_done, last_progress

   last_progress = 0

   ! allocate results: column1 store pixel id, column2 store hit count
   allocate(results(npix, 2))
   results(:, :) = 0.0d0

   ! Optionally initialize pixel ids in column 1:
   do idx = 1, npix
      results(idx, 1) = idx - 1          ! store 0-based pixel id for clarity
   end do

   call cpu_time(s_time)

do i = 1, steps
   t = start_time + real(i - 1, kind=8) * dt
   call get_vector(0, t, pixel, result_r)

   ! Debug first few steps
   if (i <= 10) then
       print *, 'Step:', i, 'pixel=', pixel, 'result_r=', result_r
   end if

   ! Correct pixel range check
   if (pixel >= 0 .and. pixel < npix) then
       results(pixel+1, 2) = results(pixel+1, 2) + 1.0d0
   else
       write(*,*) 'WARNING: pixel out of range:', pixel
   end if
   progress_done = int(1000.0d0 * real(i, kind=8) / real(steps, kind=8))
   if (mod(progress_done, 10) == 0 .and. progress_done > last_progress) then
      print *, progress_done, "/1000 parts completed"
      last_progress = progress_done
  end if
end do

   call cpu_time(e_time)
   elapsed_time = e_time - s_time
   print*, 'Total execution time (seconds):', elapsed_time

   ! Write results as formatted text file
open(unit=11, file="scan_frequency.dat", form="formatted", status="replace")
! Optionally, write a header
write(11,'(A)') '# pixel_id   hit_count'

! Write each row of results
do idx = 1, npix
   write(11,'(I8,2X,I12)') int(results(idx,1)), int(results(idx,2))
end do


   deallocate(results)

end program main

 
