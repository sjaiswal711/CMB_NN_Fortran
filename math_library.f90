module math_library
    implicit none
contains

    ! Subroutine for calculating the cross product of two vectors
    subroutine cross_product(vec1, vec2, result)
        implicit none
        real(8), intent(in) :: vec1(3), vec2(3)
        real(8), intent(out) :: result(3)

        ! Compute the cross product
        result(1) = vec1(2) * vec2(3) - vec1(3) * vec2(2)
        result(2) = vec1(3) * vec2(1) - vec1(1) * vec2(3)
        result(3) = vec1(1) * vec2(2) - vec1(2) * vec2(1)
    end subroutine cross_product


    ! Function for calculating the dot product of two vectors
    real(8) function dot_product(vec1, vec2)
        implicit none
        real(8), intent(in) :: vec1(3), vec2(3)
        integer :: i

        dot_product = 0.0d0
        do i = 1, 3
            dot_product = dot_product + vec1(i) * vec2(i)
        end do
    end function dot_product


    ! Subroutine for normalizing a vector
    subroutine normalize_vector(vector, normalized_vector)
        implicit none
        real(8), intent(in) :: vector(3)
        real(8), intent(out) :: normalized_vector(3)
        real(8) :: norm
        integer :: i

        ! Calculate the norm of the vector
        norm = sqrt(sum(vector**2))

        ! Normalize the vector
        if (norm /= 0.0d0) then
            do i = 1, 3
                normalized_vector(i) = vector(i) / norm
            end do
        else
            ! If the norm is zero, return the zero vector
            normalized_vector = 0.0d0
        end if
    end subroutine normalize_vector


    ! Function for calculating the magnitude of a vector
    real(8) function magnitude(vector)
        implicit none
        real(8), intent(in) :: vector(3)

        magnitude = sqrt(sum(vector**2))
    end function magnitude


    ! Function for calculating the angle between two vectors (in radians)
    real(8) function angle_between_vectors(vec1, vec2)
        implicit none
        real(8), intent(in) :: vec1(3), vec2(3)
        real(8) :: dp, mag1, mag2

        dp = dot_product(vec1, vec2)
        mag1 = magnitude(vec1)
        mag2 = magnitude(vec2)

        ! Calculate the angle using arccos of the normalized dot product
        if (mag1 /= 0.0d0 .and. mag2 /= 0.0d0) then
            angle_between_vectors = acos(max(-1.0d0, min(1.0d0, dp / (mag1 * mag2))))
        else
            angle_between_vectors = 0.0d0
        end if
    end function angle_between_vectors
    
    subroutine linspace(x, start_time, end_time, steps)
        implicit none
        real(8), intent(in) ::  end_time, start_time
        integer, intent(in) :: steps
        real(8), allocatable, intent(out) :: x(:)
        real(8) :: dx,duration
        integer :: i
        duration = end_time - start_time
        ! Calculate the interval between time steps
        dx = duration / real(steps-1)

        ! Allocate the time_periods array
        allocate(x(steps))

        ! Calculate the time periods using a loop
        do i = 1, steps
            x(i) = start_time + (i-1) * dx
        end do

    end subroutine linspace
end module math_library
