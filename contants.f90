module constants
   implicit none

   ! Define the parameters
   real(kind=8), parameter :: pi = 3.141592653589793d0
   real(kind=8), parameter :: theta1 = 7.5 * pi / 180.0d0
   real(kind=8), parameter :: theta2 = 85.0 * pi / 180.0d0
   real(kind=8), parameter :: w1 = 2.0d0 * pi  ! rad/min
   real(kind=8), parameter :: w3 = 0.000011954d0  ! rad/min
   real(kind=8), parameter :: w2 = 2.0d0 * w3  ! rad/min

   integer, parameter :: nside = 1024
   integer, parameter :: npix = 12 * nside * nside

   real(kind=8), parameter :: scan_time = sqrt(4 * pi / npix) / w1
   real(kind=8), parameter :: grid_size = 2.0d0 * pi / (180.0d0 * 3600.0d0)
   real(kind=8), parameter :: Radius = (99.0d0 / 60.0d0) * (pi / 180.0d0)
   integer, parameter :: length = 2800
   integer, parameter :: centre(2) = [3001, 3004]

   real(kind=8), parameter :: start_time = 0.0d0
   real(kind=8), parameter :: duration = 1.0d0 * 365 * 60 *24

   !integer, parameter :: steps = int(duration / scan_time) - 1
   integer, parameter :: long_kind = selected_int_kind(18)
   integer(long_kind), parameter :: steps = int(duration / scan_time, long_kind) - 1

   real(kind=8), parameter :: dt = duration / real(steps - 1)

end module constants
