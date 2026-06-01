!----------------------------------------program------------------------------------
!****************** evaluate band structure for multiple lattices
program tb

    use brillouin_zone_mod
    use fermi_bisection_mod
    use dispersion_interface_mod

    implicit none

    real(8) :: kx, ky, kz, t_start, t_end
    integer :: info, k, l, npts, band, i, j, NE
    complex(8) :: ef
    real(8), allocatable :: recbv(:,:)
    real(8) :: volume
    real(8), dimension(3,3) :: bravais
    real(8) :: rotmat(3,3,3)
    integer :: nkx, nky, nkz, nsymbz, kpoibz, maxmesh
    integer, allocatable :: NOFKS(:), IBK(:,:,:,:), eigsave3(:,:,:)
    real(8), allocatable :: BZKP(:,:,:), VOLCUB(:,:)
    real(8) :: tauvbz

    character(len=5) :: green_ansi, reset_ansi
    character(len=1) :: esc

    esc = achar(27)
    green_ansi = esc // '[32m'
    reset_ansi = esc // '[0m'

    ! Print initialization banner
    write(6,*) " "
    write(6,*) green_ansi // "========================================================" // reset_ansi
    write(6,*) green_ansi // "   STARTING MULTI-LATTICE HIGH-PRECISION DOS SWEEP      " // reset_ansi
    write(6,*) green_ansi // "========================================================" // reset_ansi
    write(6,*) " "

    ! Setup uniform grid dimensions for the BZ mesh (50^3 grid provides fast exact integration)
    npts = 50
    nkx = npts
    nky = npts
    nkz = npts
    maxmesh = 1
    nsymbz = 1

    ! Setup standard cubic Bravais lattice
    do i = 1, 3
        do j = 1, 3
            if (i .eq. j) then
                bravais(i,j) = 1.0d0
            else
                bravais(i,j) = 0.0d0
            end if
        end do
    end do

    ! Standard identity reciprocal lattice vectors for uniform BZ sampling in fractional coordinates
    allocate(recbv(3,3))
    recbv = 0.0d0
    recbv(1,1) = 1.0d0
    recbv(2,2) = 1.0d0
    recbv(3,3) = 1.0d0

    allocate(IBK(0:nkx, 0:nky, 0:nkz, 2))

    call cpu_time(t_start)

    ! =========================================================================
    ! 1. SIMPLE CUBIC (SC) SWEEP
    ! =========================================================================
    write(6,*) green_ansi // ">>> Running Simple Cubic (SC) Sweep..." // reset_ansi
    IBK = 0
    call kp_gen(bravais, recbv, nsymbz, rotmat, nkx, nky, nkz, kpoibz, maxmesh, NOFKS, BZKP, VOLCUB, tauvbz, .false., IBK)
    
    NE = 80  ! Number of energy points for smooth curves
    open(10, file='data/DOS_sc.txt', status='replace')
    do i = 1, NE
        ! Sweep SC band range: -5.9 eV to +5.9 eV
        ef = cmplx(-5.9d0 + (11.8d0 / (NE - 1)) * (i - 1), 0.0d0, 8)
        
        allocate(eigsave3(NOFKS(1),1,1))
        do l = 1, NOFKS(1)
            kx = BZKP(1,l,1)
            ky = BZKP(2,l,1)
            kz = BZKP(3,l,1)
            if (sc_dispersion(kx, ky, kz) .lt. dreal(ef)) then
                eigsave3(l,1,1) = 1
            else
                eigsave3(l,1,1) = 0
            end if
        end do
        
        call BISECTION1(1, 1, 1, 1, nkx, nky, nkz, nsymbz, ef, IBK, recbv, eigsave3(:,1,1), rotmat, sc_dispersion, sc_grad_mag)
        deallocate(eigsave3)
    end do
    close(10)
    deallocate(NOFKS, BZKP, VOLCUB)

    ! =========================================================================
    ! 2. BODY-CENTRED CUBIC (BCC) SWEEP
    ! =========================================================================
    write(6,*) green_ansi // ">>> Running Body-Centred Cubic (BCC) Sweep..." // reset_ansi
    IBK = 0
    call kp_gen(bravais, recbv, nsymbz, rotmat, nkx, nky, nkz, kpoibz, maxmesh, NOFKS, BZKP, VOLCUB, tauvbz, .false., IBK)
    
    open(10, file='data/DOS_bcc.txt', status='replace')
    do i = 1, NE
        ! Sweep BCC band range: -7.9 eV to +7.9 eV
        ef = cmplx(-7.9d0 + (15.8d0 / (NE - 1)) * (i - 1), 0.0d0, 8)
        
        allocate(eigsave3(NOFKS(1),1,1))
        do l = 1, NOFKS(1)
            kx = BZKP(1,l,1)
            ky = BZKP(2,l,1)
            kz = BZKP(3,l,1)
            if (bcc_dispersion(kx, ky, kz) .lt. dreal(ef)) then
                eigsave3(l,1,1) = 1
            else
                eigsave3(l,1,1) = 0
            end if
        end do
        
        call BISECTION1(1, 1, 1, 1, nkx, nky, nkz, nsymbz, ef, IBK, recbv, eigsave3(:,1,1), rotmat, bcc_dispersion, bcc_grad_mag)
        deallocate(eigsave3)
    end do
    close(10)
    deallocate(NOFKS, BZKP, VOLCUB)

    ! =========================================================================
    ! 3. FACE-CENTRED CUBIC (FCC) SWEEP
    ! =========================================================================
    write(6,*) green_ansi // ">>> Running Face-Centred Cubic (FCC) Sweep..." // reset_ansi
    IBK = 0
    call kp_gen(bravais, recbv, nsymbz, rotmat, nkx, nky, nkz, kpoibz, maxmesh, NOFKS, BZKP, VOLCUB, tauvbz, .false., IBK)
    
    open(10, file='data/DOS_fcc.txt', status='replace')
    do i = 1, NE
        ! Sweep FCC band range: -11.9 eV to +3.9 eV
        ef = cmplx(-11.9d0 + (15.8d0 / (NE - 1)) * (i - 1), 0.0d0, 8)
        
        allocate(eigsave3(NOFKS(1),1,1))
        do l = 1, NOFKS(1)
            kx = BZKP(1,l,1)
            ky = BZKP(2,l,1)
            kz = BZKP(3,l,1)
            if (fcc_dispersion(kx, ky, kz) .lt. dreal(ef)) then
                eigsave3(l,1,1) = 1
            else
                eigsave3(l,1,1) = 0
            end if
        end do
        
        call BISECTION1(1, 1, 1, 1, nkx, nky, nkz, nsymbz, ef, IBK, recbv, eigsave3(:,1,1), rotmat, fcc_dispersion, fcc_grad_mag)
        deallocate(eigsave3)
    end do
    close(10)
    deallocate(NOFKS, BZKP, VOLCUB)

    ! =========================================================================
    ! 4. SIMPLE HEXAGONAL (HEX) SWEEP
    ! =========================================================================
    write(6,*) green_ansi // ">>> Running Simple Hexagonal (HEX) Sweep..." // reset_ansi
    IBK = 0
    call kp_gen(bravais, recbv, nsymbz, rotmat, nkx, nky, nkz, kpoibz, maxmesh, NOFKS, BZKP, VOLCUB, tauvbz, .false., IBK)
    
    open(10, file='data/DOS_hex.txt', status='replace')
    do i = 1, NE
        ! Sweep HEX band range: -6.9 eV to +3.9 eV
        ef = cmplx(-6.9d0 + (10.8d0 / (NE - 1)) * (i - 1), 0.0d0, 8)
        
        allocate(eigsave3(NOFKS(1),1,1))
        do l = 1, NOFKS(1)
            kx = BZKP(1,l,1)
            ky = BZKP(2,l,1)
            kz = BZKP(3,l,1)
            if (hex_dispersion(kx, ky, kz) .lt. dreal(ef)) then
                eigsave3(l,1,1) = 1
            else
                eigsave3(l,1,1) = 0
            end if
        end do
        
        call BISECTION1(1, 1, 1, 1, nkx, nky, nkz, nsymbz, ef, IBK, recbv, eigsave3(:,1,1), rotmat, hex_dispersion, hex_grad_mag)
        deallocate(eigsave3)
    end do
    close(10)
    deallocate(NOFKS, BZKP, VOLCUB)

    call cpu_time(t_end)
    write(6,*) " "
    write(6,*) green_ansi // "========================================================" // reset_ansi
    write(6,*) green_ansi // "   COMPLETED ALL SWEEPS SUCCESSFULY!                    " // reset_ansi
    write(6,1020) t_end - t_start
1020 FORMAT("    Total execution time: ", F8.2, " seconds")
    write(6,*) green_ansi // "========================================================" // reset_ansi
    write(6,*) " "

    deallocate(recbv, IBK)

contains

    ! =========================================================================
    ! LATTICE DISPERSION AND GRADIENT MAGNITUDE FUNCTIONS
    ! =========================================================================

    ! 1. Simple Cubic (SC)
    real(8) function sc_dispersion(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi
        pi = 3.1415926535897932d0
        sc_dispersion = -2.d0 * (cos((kx-0.5d0)*2.d0*pi) + cos((ky-0.5d0)*2.d0*pi) + cos((kz-0.5d0)*2.d0*pi))
    end function sc_dispersion

    real(8) function sc_grad_mag(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi
        pi = 3.1415926535897932d0
        sc_grad_mag = 4.d0 * pi * sqrt(sin((kx-0.5d0)*2.d0*pi)**2 + &
                                       sin((ky-0.5d0)*2.d0*pi)**2 + &
                                       sin((kz-0.5d0)*2.d0*pi)**2)
        if (sc_grad_mag < 1.d-10) sc_grad_mag = 1.d-10
    end function sc_grad_mag

    ! 2. Body-Centred Cubic (BCC)
    real(8) function bcc_dispersion(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi
        pi = 3.1415926535897932d0
        bcc_dispersion = -8.d0 * cos((kx-0.5d0)*2.d0*pi) * cos((ky-0.5d0)*2.d0*pi) * cos((kz-0.5d0)*2.d0*pi)
    end function bcc_dispersion

    real(8) function bcc_grad_mag(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi, cx, cy, cz, sx, sy, sz
        pi = 3.1415926535897932d0
        cx = cos((kx-0.5d0)*2.d0*pi); sx = sin((kx-0.5d0)*2.d0*pi)
        cy = cos((ky-0.5d0)*2.d0*pi); sy = sin((ky-0.5d0)*2.d0*pi)
        cz = cos((kz-0.5d0)*2.d0*pi); sz = sin((kz-0.5d0)*2.d0*pi)
        bcc_grad_mag = 16.d0 * pi * sqrt((sx*cy*cz)**2 + (cx*sy*cz)**2 + (cx*cy*sz)**2)
        if (bcc_grad_mag < 1.d-10) bcc_grad_mag = 1.d-10
    end function bcc_grad_mag

    ! 3. Face-Centred Cubic (FCC)
    real(8) function fcc_dispersion(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi, cx, cy, cz
        pi = 3.1415926535897932d0
        cx = cos((kx-0.5d0)*2.d0*pi)
        cy = cos((ky-0.5d0)*2.d0*pi)
        cz = cos((kz-0.5d0)*2.d0*pi)
        fcc_dispersion = -4.d0 * (cx*cy + cy*cz + cz*cx)
    end function fcc_dispersion

    real(8) function fcc_grad_mag(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi, cx, cy, cz, sx, sy, sz
        pi = 3.1415926535897932d0
        cx = cos((kx-0.5d0)*2.d0*pi); sx = sin((kx-0.5d0)*2.d0*pi)
        cy = cos((ky-0.5d0)*2.d0*pi); sy = sin((ky-0.5d0)*2.d0*pi)
        cz = cos((kz-0.5d0)*2.d0*pi); sz = sin((kz-0.5d0)*2.d0*pi)
        fcc_grad_mag = 8.d0 * pi * sqrt((sx*(cy + cz))**2 + (sy*(cx + cz))**2 + (sz*(cx + cy))**2)
        if (fcc_grad_mag < 1.d-10) fcc_grad_mag = 1.d-10
    end function fcc_grad_mag

    ! 4. Simple Hexagonal (HEX)
    real(8) function hex_dispersion(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi
        pi = 3.1415926535897932d0
        hex_dispersion = -2.d0 * (cos(4.d0*pi*(kx-0.5d0)) + &
                         2.d0 * cos(2.d0*pi*(kx-0.5d0)) * cos(2.d0*pi*sqrt(3.d0)*(ky-0.5d0))) &
                         - 2.d0 * 0.5d0 * cos(2.d0*pi*(kz-0.5d0))
    end function hex_dispersion

    real(8) function hex_grad_mag(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi, kx_c, ky_c, kz_c, gx, gy, gz
        pi = 3.1415926535897932d0
        kx_c = kx - 0.5d0
        ky_c = ky - 0.5d0
        kz_c = kz - 0.5d0
        gx = 8.d0*pi*sin(4.d0*pi*kx_c) + 8.d0*pi*sin(2.d0*pi*kx_c)*cos(2.d0*pi*sqrt(3.d0)*ky_c)
        gy = 8.d0*sqrt(3.d0)*pi*cos(2.d0*pi*kx_c)*sin(2.d0*pi*sqrt(3.d0)*ky_c)
        gz = 2.d0*pi*sin(2.d0*pi*kz_c)
        hex_grad_mag = sqrt(gx**2 + gy**2 + gz**2)
        if (hex_grad_mag < 1.d-10) hex_grad_mag = 1.d-10
    end function hex_grad_mag

end program tb
