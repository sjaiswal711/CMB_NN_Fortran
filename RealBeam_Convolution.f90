module constants
    implicit none

    ! Define the parameters
    real(kind=8), parameter :: pi = 3.141592653589793d0
    real(kind=8), parameter :: theta1 = 7.5 * pi / 180.0d0
    real(kind=8), parameter :: theta2 = 85.0 * pi / 180.0d0
    real(kind=8), parameter :: w1 = 2.0d0 * pi  ! rad/min
    real(kind=8), parameter :: w2 = 2.0d0 * w1  ! rad/min
    real(kind=8), parameter :: w3 = 0.000011954d0  ! rad/min

    integer, parameter :: nside = 1024
    integer, parameter :: npix = 12 * nside * nside

    real(kind=8), parameter :: scan_time = sqrt(4 * pi / npix) / w1
    real(kind=8), parameter :: grid_size = 2.0d0 * pi / (180.0d0 * 3600.0d0)
    real(kind=8), parameter :: Radius = (100.0d0 / 60.0d0) * (pi / 180.0d0)
    integer, parameter :: length = 2630
    integer, parameter :: centre(2) = [3001, 3004]
    
    real(kind=8), parameter :: start_time = 0.0d0
    real(kind=8), parameter :: duration = 1.0d0 * 365 * 60 *24 
    
    !integer, parameter :: steps = int(duration / scan_time) - 1
    integer, parameter :: long_kind = selected_int_kind(18)
	integer(long_kind), parameter :: steps = int(duration / scan_time, long_kind) - 1

    real(kind=8), parameter :: dt = duration / real(steps - 1)
    
end module constants

module subroutines
	use math_library
	use healpix_modules
    use constants
    implicit none
    
    contains

    subroutine read_grid(grid, nx, ny)

        integer, intent(in) :: nx, ny
        real*8, allocatable, dimension(:,:) :: grid
        integer :: i, j, ios

        ! Allocate memory for the grid
        allocate(grid(nx, ny))

        ! Open the file and read the data column by column
        open(unit=10, file='grid.txt', status='old', access='sequential', &  
            form='formatted')
        do i = 1, nx
            read(10, *, iostat=ios) (grid(i, j), j = 1, ny)
            if (ios /= 0) then
                write(*,*) 'Error reading column', i
                stop
            end if
        end do
        close(10)
    end subroutine read_grid
    
    subroutine anglev(vec1, vec2, angle)
    use math_library
    implicit none
    real(kind=8), intent(in) :: vec1(3), vec2(3)
    real(kind=8), intent(out) :: angle
    real(kind=8) :: dot_prod, clipped_dp
    integer :: i=0

    dot_prod = dot_product(vec1, vec2)

    ! Clip dot_product to the valid range for arccos to avoid NaNs
    clipped_dp = max(-1.0d0, min(1.0d0, dot_prod))
    !print*,clipped_dp

    ! Calculate the angle
    angle = acos(clipped_dp)
end subroutine anglev
    
    subroutine get_vectors(node_id, t, pixel, result_r, result_s)

        real(kind=8), intent(in) :: t
        integer, intent(in) :: node_id
        real(kind=8), intent(out) :: result_r(3,1), result_s(3,1)
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
            cos_w1t, sin_w1t, 0.0d0, &
            -sin_w1t, cos_w1t, 0.0d0, &
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
        D_R = reshape([cos_theta2, sin_theta2 * cos(w3 * t), sin_theta2 * sin(w3 * t)], [3, 1], order=[2,1])
        D_S = reshape([1.0d0, 0.0d0, 0.0d0], [3, 1], order=[2,1])

        ! Perform matrix multiplications (optimized)
        result1 = matmul(matmul(A,B),C)
        result_r = matmul(result1, D_R)
        result_s = matmul(result1, D_S)
  
        call vec2pix_ring(nside, reshape(result_r, [3]), pixel)

        !print*,"node_id, t, pixel",node_id, t, pixel
        !print*,"cos(w1 * t), sin(w1 * t), cos(w2 * t), sin(w2 * t), cos(w3 * t), sin(w3 * t)"
        !print*,cos(w1 * t), sin(w1 * t), cos(w2 * t), sin(w2 * t), cos(w3 * t), sin(w3 * t)
    end subroutine get_vectors
    

subroutine process_time_step_new(node_id, grid, R, S, pix_ring, results_t)
    use math_library 
    use constants
    use healpix_modules
    implicit none
    !real(8), dimension(:,:), intent(in) :: input_map
    integer, intent(in) :: node_id
    real*8, dimension(6001,6001), intent(in) :: grid
    real(kind=8), intent(in) :: R(3), S(3)
    integer, intent(in) :: pix_ring
    real(kind=8), intent(out) :: results_t(2*length+3)
    
    integer, allocatable :: listpix(:)
    
    integer :: nlist, i, index_x, index_y

    real(8) :: x(length)
    real(8) :: y(length)
    real(8) :: weight(length)
    real(8) :: R_i(3), ABC(3,3)
    real(8) :: Z_t(3), I_t(3), N_t(3),nor_A_i(3), nor_N_t(3)
    real(8) :: Rc(3)
    real(8) :: theta_i, alpha_i, A_i(3)
    
    real(8), dimension(length) :: neighbor_temperatures

    allocate(listpix(length))

    results_t = 0
    results_t(1) = node_id
	results_t(2)=pix_ring
	results_t(3)=1
    ! 3. Calculate Z_t, I_t, and N_t
    call cross_product(R, S, Z_t)
    call cross_product(R, Z_t, I_t)
    N_t = I_t

    ! 5. Get the vector centered at pix_ring
    call pix2vec_ring(nside, pix_ring, Rc)

    ! 6. Find neighboring pixels in RING format
    call query_disc(nside, Rc, Radius, listpix, nlist, nest=0)

    do i = 1, nlist
        ! 7. Get the vector corresponding to the neighboring pixel
        call pix2vec_ring(nside, listpix(i), R_i)
        
        ! 8. Angular separation between central pixel and the neighboring pixels
        call anglev(Rc, R_i, theta_i)

        ! 9. A_i = line joining central pixel and neighbor pixel
        A_i = Rc - R_i

        ! 10. Angle between N_t and A_i
        call normalize_vector(A_i, nor_A_i)
        call normalize_vector(N_t, nor_N_t)
        call anglev(nor_A_i, nor_N_t, alpha_i)
        
        ! 11. Calculate the weights of the neighboring pixel
        x(i) = theta_i * cos(alpha_i)
        y(i) = theta_i * sin(alpha_i)
		
        index_x = int(centre(1) + nint(x(i) / grid_size))
        index_y = int(centre(2) + nint(y(i) / grid_size))
        weight(i) = grid(index_x, index_y)

		results_t(2*i + 2)=listpix(i)
		results_t(2*i + 3)=weight(i)
		
    end do
    
end subroutine process_time_step_new

end module subroutines

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
    
    real(kind=8) :: results_t(2*length+3)
    real(kind=8), allocatable :: results(:,:)
    integer :: idx
    integer :: unit
	
	real(kind=8) :: cnvl_temp
    
    real(kind=8) :: s_time, e_time, elapsed_time
    
    last_progress = 0
    
    ! Initialize variables
    
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
                
    pixel_ranges_25 = reshape((/ 0,  &
                117965,  267139, 527354, 887670, &
                1340977, 1879012, 2492572,  &
                3171299,  3903423,  4676305,  5476887, &
                 6291499, 7106141,  7906735,  &
                8679581, 9411734, 10090394, 10704008, &
                 11242033,  11695305,  12055622, &
                 12315886, 12465037,  12582912 /),(/25/))
                
	pixel_ranges_1 = reshape( (/ 0, 49588 /), (/ 2 /)) 
	pixel_ranges_2 = reshape( (/ 0, nside*nside*6, npix /), (/ 3 /))  
    
    t = 1.0d0  ! Time input, change as necessary
    
    ! Determine the pixel range for this node
    pixel_start = pixel_ranges(rank + 1)
    pixel_end = pixel_ranges(rank + 2) 
	allocate(results((pixel_end - pixel_start + 1), 2*length+3))

    print*, "Node ID: ", rank
    print*, "Processing pixel range: ", pixel_start, " to ", pixel_end
    
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
        		progress_done = int(10000.0d0 * i / steps)

        		if (progress_done > last_progress) then
            		print *, progress_done, "/ 10000 parts completed"
            		last_progress = progress_done
        		end if
    end do

    print*, "Parallel loop done"

    ! Record end time
    call cpu_time(e_time)
    
    ! Calculate and print elapsed time
    elapsed_time = e_time - s_time
    print*, 'Total execution time (seconds):', elapsed_time
    
    
    write(filename, '("results_", I0, ".dat")') rank
        ! Open the file for writing
    unit = 10 + rank  ! Each node uses a unique unit number
    open(unit=unit, file=filename, status='replace', form='formatted')

do i = pixel_start+1, pixel_end+1
    do j = 1, 2*length+3
        if (j < 4) then
            write(unit, '(I10)', advance='no') int(results(i, j))
        else if (mod(j, 2) /= 0) then
            write(unit, '(F10.5)', advance='no') results(i, j)/results(i, 3)
        else
            write(unit, '(I10)', advance='no') int(results(i, j)/results(i, 3))
        end if
        
        ! Add a comma separator, except after the last element
        if (j < 2*length+3) then
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
