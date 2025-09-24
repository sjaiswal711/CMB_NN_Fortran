program main
   use constants
   use healpix_modules
   use subroutines
   use mpi
   !use omp_lib
   implicit none

   !|- -| - -| - -| - -| Declare Variables |- -| - -| - -| - -|!

   ! PARAMETERS
   real(kind=8) :: t
   integer :: pixel
   real(kind=8) :: result_r(3,1), result_s(3,1)
   integer  :: j, k
   integer(long_kind) :: i

   integer :: progress_interval, last_progress, progress_done

   !FILENAMES
   character(len=20) :: map_file
   character(len=20) :: grid_file
   character(len=20) :: result_file
   character(len=20) :: filename

   !To load_map
   integer :: nmaps
   real(kind=8), dimension(:,:), allocatable :: map
   integer :: total_pixels

   !To load data_grid
   real(kind=8), allocatable :: data_grid(:,:)
   integer, parameter :: nx = 6001, ny = 6001

   !for parallel processing
   integer, dimension(3) :: pixel_ranges_2
   integer, dimension(2) :: pixel_ranges_1
   integer, dimension(49) :: pixel_ranges
   integer, dimension(25) :: pixel_ranges_25

   integer :: pixel_start, pixel_end
   integer :: node_id, ierr
   integer :: num_nodes, rank
   character(len=10) :: rank_str

   real(kind=8) :: results_t(length+3)
   real(kind=8), allocatable :: results(:,:)
   integer :: idx
   integer :: unit
   integer :: local_pixel_count

   real(kind=8) :: cnvl_temp

   real(kind=8) :: s_time, e_time, elapsed_time
   real(kind=8), dimension(length+3) :: row
   integer(long_kind) :: nrows, ncols
   ! real(8), allocatable :: row(:)

   last_progress = 0
   total_pixels = 1024*1024*12
   ! Initialize variables
   !allocate(results(npix, length+3))
   !print*,results
   !FILENAMES
   map_file = "nested.fits"
   grid_file = "grid.txt"
   result_file = "result_year.dat"

   ! LOAD MAP
   !nmaps = 1
   !allocate(map(0:npix-1, 1:nmaps))
   !call input_map(map_file, map, npix, nmaps)

   ! LOAD GRID
   call read_grid(data_grid, nx, ny)
   !print *, "Printing a small section of the grid (3001:3010, 3001:3010):"
   !print *, "-----------------------------------------------"
   !do i = 3001, 3010
   !    do j = 3001, 3010
   !       print *, "data_grid(", i, ",", j, ") = ", data_grid(i, j)
   !	end do
   !    write(*,*)
   !end do



   call MPI_INIT(ierr)
   call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
   call MPI_COMM_SIZE(MPI_COMM_WORLD, num_nodes, ierr)
   pixel_ranges = reshape((/ 0, 49588, &
      117965, 178159, 267139, 384304, 527354, 695404, 887670, &
      1103107, 1340977, 1599972, 1879012, 2176998, 2492572, 2824460, &
      3171299, 3531365, 3903423, 4285623, 4676305, 5074091, 5476887, &
      5883452, 6291499, 6699587, 7106141, 7508964, 7906735, 8297424, &
      8679581, 9051651, 9411734, 9758509, 10090394, 10405996, 10704008, &
      10983044, 11242033, 11479857, 11695305, 11887537, 12055622, &
      12198680, 12315886, 12404854, 12465037, 12533405, 12582912 /),(/49/))
   ! do i = 1, 49
   !    pixel_ranges(i) = int( real(total_pixels) * real(i-1) / real(48) )
   ! end do

   t = 1.0d0  ! Time input, change as necessary

   ! Determine the pixel range for this node
   pixel_start = pixel_ranges(rank + 1)
   pixel_end = pixel_ranges(rank + 2)
   local_pixel_count = pixel_end - pixel_start + 1
   allocate(results(local_pixel_count, length + 3))
   results = 0.0d0  ! Initialize to zero

   print*, "Node ID: ", rank
   print*, "Processing pixel range: ", pixel_start, " to ", pixel_end
   print *, "Local pixel count:", local_pixel_count

   ! Convert rank to string
   write(rank_str, '(I0)') rank

   call cpu_time(s_time)

   ! process the map for each time step in `time_periods`
   do i = 1, steps
      t = start_time + (i-1) * dt
      call get_vectors(rank, t, pixel, result_r, result_s)
      !print*, "node_id,pixel,t,pixel_start,pixel_end"
      !print*, node_id,pixel,t,pixel_start,pixel_end
      if (pixel >= pixel_start .and. pixel < pixel_end) then
         ! Process the time step

         !print*,"t Pixel Node_id",t,pixel,node_id
         !call process_time_step(rank, data_grid, map, t, result_r, result_s,  pixel, cnvl_temp)
         call process_time_step_new(rank, data_grid, result_r, result_s, pixel, results_t)

         idx = int(results_t(2))+1 - pixel_start
         results(idx,2) = results_t(2)
         results(idx,3:) = results(idx,3:) + results_t(3:)
         !print*,results_t(3:)
         !print*,map

         !grid, input_map, time_step, R, S,  pix_ring, cnvl_temp
      end if
      ! Check if we have reached the next progress interval
      progress_done = int(1000.0d0 * i / steps)

      if (progress_done > last_progress) then
         print *, progress_done, "/ 1000 parts completed",rank
         last_progress = progress_done
      end if
   end do

   print*, "Parallel loop done"

   ! Record end time
   call cpu_time(e_time)

   ! Calculate and print elapsed time
   elapsed_time = e_time - s_time
   print*, 'Total execution time (seconds):', elapsed_time


   write(filename, '("beam_response_mat_", I0, ".dat")') rank
   ! Open the file for writing
   unit = 10 + rank  ! Each node uses a unique unit number
   open(unit=unit, file=filename, status='replace', form='formatted')

   do i = 1, pixel_end - pixel_start
      do j = 1, length+3
         if (j < 4) then
            write(unit, '(I10)', advance='no') int(results(i, j))
         else
            write(unit, '(F10.8)', advance='no') results(i, j)/results(i, 3)
         end if

         ! Add a comma separator, except after the last element
         if (j < length+3) then
            write(unit, '(",")', advance='no')
         end if
      end do
      write(unit, *)  ! Newline after each row
   end do

   ! Close the file
   close(unit)
   

   call MPI_FINALIZE(ierr)
   deallocate(data_grid)
   ! Finalize MPI


   !call get_vectors(t, pixel, result_r, result_s)


end program main
