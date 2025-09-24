subroutine get_neighbour(node_id, pix_ring, results_t)
    use math_library
    use constants
    use healpix_modules
    implicit none

    integer, intent(in) :: node_id
    integer, intent(in) :: pix_ring
    real(kind=8), intent(out) :: results_t(length+3)

    real(kind=8) :: Rc(3)
    integer :: nlist, i
    integer, allocatable :: listpix(:)
    allocate(listpix(length))

    results_t = 0.0d0
    results_t(1) = node_id
    results_t(2) = pix_ring
    results_t(3) = 1.0d0

    call pix2vec_ring(nside, pix_ring, Rc)
    call query_disc(nside, Rc, Radius, listpix, nlist, nest=0)

    do i = 1, nlist
        results_t(i + 3) = listpix(i)
    end do

    deallocate(listpix)
end subroutine get_neighbour


program main
    use constants
    use healpix_modules
    use mpi
    implicit none

    real(kind=8) :: results_t(length+3)
    real(kind=8), allocatable :: results(:,:)

    integer :: progress_done, last_progress
    integer :: idx, i, j, unit
    integer :: local_pixel_count
    integer, dimension(49) :: pixel_ranges
    real(kind=8) :: s_time, e_time, elapsed_time
    integer :: pixel_start, pixel_end
    integer :: ierr, num_nodes, rank
    character(len=20) :: filename

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, num_nodes, ierr)

    pixel_ranges = reshape((/ &
        0, 49588, 117965, 178159, 267139, 384304, 527354, 695404, 887670, &
        1103107, 1340977, 1599972, 1879012, 2176998, 2492572, 2824460, &
        3171299, 3531365, 3903423, 4285623, 4676305, 5074091, 5476887, &
        5883452, 6291499, 6699587, 7106141, 7508964, 7906735, 8297424, &
        8679581, 9051651, 9411734, 9758509, 10090394, 10405996, 10704008, &
        10983044, 11242033, 11479857, 11695305, 11887537, 12055622, &
        12198680, 12315886, 12404854, 12465037, 12533405, 12582911 /), (/49/))

    pixel_start = pixel_ranges(rank+1)
    pixel_end = pixel_ranges(rank+2)
    local_pixel_count = pixel_end - pixel_start + 1

    allocate(results(local_pixel_count, length+3))
    results = 0.0d0

    print*, "Node ID:", rank
    print*, "Processing pixel range:", pixel_start, "to", pixel_end
    print*, "Local pixel count:", local_pixel_count

    last_progress = -1
    call cpu_time(s_time)

    do i = pixel_start, pixel_end
        call get_neighbour(rank, i, results_t)
        idx = i - pixel_start + 1
        results(idx,2) = results_t(2)
        results(idx,3:) = results(idx,3:) + results_t(3:)

        progress_done = int(100.0d0 * idx / local_pixel_count)
        if (progress_done > last_progress) then
            print*, progress_done, "/ 100 parts completed"
            last_progress = progress_done
        end if
    end do

    print*, "Parallel loop done"

    call cpu_time(e_time)
    elapsed_time = e_time - s_time
    print*, 'Total execution time (seconds):', elapsed_time

    write(filename, '("neighbors_mat_", I0, ".dat")') rank
    unit = 10 + rank
    open(unit=unit, file=filename, status='replace', form='formatted')

    do i = 1, local_pixel_count
        do j = 1, length+3
            write(unit, '(I10)', advance='no') int(results(i,j))
            if (j < length+3) write(unit, '(",",A)', advance='no') ','
        end do
        write(unit, *)
    end do

    close(unit)
    deallocate(results)

    call MPI_FINALIZE(ierr)
end program main
