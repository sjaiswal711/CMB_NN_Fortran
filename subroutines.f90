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



   subroutine process_time_step(node_id, grid, input_map, time_step, R, S,  pix_ring, cnvl_temp)
      use math_library
      use constants
      use healpix_modules
      implicit none
      !real(kind=8), dimension(:,:), intent(in) :: input_map
      real(kind=8), dimension(6001,6001), intent(in) :: grid
      real(kind=8), dimension(:,:), intent(in) :: input_map
      integer, intent(in) :: node_id
      real(kind=8), intent(in) :: time_step
      real(kind=8), intent(in) :: R(3), S(3)
      integer, intent(in) :: pix_ring
      real(kind=8), intent(out) :: cnvl_temp

      integer, allocatable :: listpix(:)
      integer :: nlist, i, index_x, index_y

      real(kind=8) :: x(length)
      real(kind=8) :: y(length)
      real(kind=8) :: weight(length)
      real(kind=8) :: t
      real(kind=8) :: R_i(3), ABC(3,3)
      real(kind=8) :: Z_t(3), I_t(3), N_t(3),nor_A_i(3), nor_N_t(3)
      real(kind=8) :: Rc(3)
      real(kind=8) :: theta_i, alpha_i, A_i(3)

      real(kind=8), dimension(length) :: neighbor_temperatures

      allocate(listpix(length))
      !print*,input_map
      t=time_step

      ! 1. Calculate R(t) and S(t) vectors
      !call get_vectors_ABC(t, ABC)
      !call get_vector_R(t, ABC, R)
      !call get_vector_S(t, ABC, S)

      ! 2. Calculate pixel number along R(t) vector (ring format)
      !call vec2pix_ring(nside, R, pix_ring)

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

         ! 12. Retrieve the temperature of the neighboring pixel
         neighbor_temperatures(i) = input_map(listpix(i)+1, 1)
         !print*,"neighbor_temperatures(i) : ",neighbor_temperatures(i)
         !print*,"neighbor_temperatures(i)=",listpix(i),input_map(listpix(i)+1, 1)
      end do

      ! 13. Apply convolution
      cnvl_temp = 0.0d0
      do i = 1, nlist
         cnvl_temp = cnvl_temp + (weight(i)*neighbor_temperatures(i))
         !print*,pix_ring,listpix(i)+1,neighbor_temperatures(i),weight(i)
      end do
      !print*,"node_id, time_step, pix_ring, cnvl_temp, pix_ring_temp"
      !print*,node_id, time_step, pix_ring, cnvl_temp, input_map(pix_ring+1, 1)

   end subroutine process_time_step

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
      real(kind=8), intent(out) :: results_t(length+3)

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

         !results_t(2*i + 2)=listpix(i)
         results_t(i + 3)=weight(i)

      end do

   end subroutine process_time_step_new

end module subroutines
