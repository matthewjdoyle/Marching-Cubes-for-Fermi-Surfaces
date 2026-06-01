module dispersion_interface_mod
    implicit none

    ! Abstract interface defining the dispersion function signature
    abstract interface
        real(8) function dispersion_fn(kx, ky, kz)
            real(8), intent(in) :: kx, ky, kz
        end function dispersion_fn

        real(8) function grad_mag_fn(kx, ky, kz)
            real(8), intent(in) :: kx, ky, kz
        end function grad_mag_fn
    end interface

end module dispersion_interface_mod
