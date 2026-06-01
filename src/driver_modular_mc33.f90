!----------------------------------------program------------------------------------
!****************** evaluate band strcuture
program tb

!use hamiltonian

use brillouin_zone_mod
use fermi_bisection_mod
use dispersion_interface_mod

Implicit None

real(8)::kx,ky,kz,buffer,kweight,T
integer::info,k,l,npts,band,i,ii,j,jj,m,n
complex(8)::ef
complex(8), allocatable :: vf(:)
real(8)::evalns
real(8),allocatable::recbv(:,:)
real(8) :: volume
real(8), allocatable :: kmesh(:,:)  !QP or BZKP in gen routines
real(8), dimension(100) :: kxs, kys, kzs, Es

real(8) :: pi, t_start, t_end

real(8), dimension(3,3) :: bravais
real(8) :: rotmat(3,3,3)  ! rotation matrices of real lattice
logical LINTERFACE
integer(8):: NYSMBZ
integer, allocatable :: NOFKS(:),IBK(:,:,:,:),eigsave2(:,:,:)
real(8), allocatable :: BZKP(:,:,:), VOLCUB(:,:)
real(8) :: tauvbz
integer :: nkx, nky, nkz, nsymbz, kpoibz, maxmesh, NE
!----------------------------------------------------end of preamble------------------


volume=0.0d+00
npts=100                ! number of k points in each direction
NE = 1                ! number of fermi energies to loop over
pi = 3.141592653589793115997963468544185161590576171875d+00

allocate(recbv(1:3,1:3))
recbv(1,1)=1.00000000d0
recbv(2,1)=0.00000000d0
recbv(3,1)=0.00000000d+00
recbv(1,2)=0.00000000d0
recbv(2,2)=1.00000000d0
recbv(3,2)=0.000000000d0
recbv(1,3)=0.00000000d0
recbv(2,3)=0.00000000d0
recbv(3,3)=1.00000000d0
recbv=recbv
nkx=npts
nky=npts
nkz=npts
allocate(IBK(0:nkX,0:nky,0:nkz,2))
         IBK=0
MAXMESH=1
NSYMBZ=1

! FORM SIMPLE CUBIC BRAVAIS LATTICE.
DO i=1,3
    DO j=1,3
        IF (i.eq.j) THEN
            bravais(i,j) = 1d0
        ELSEIF (i.ne.j) THEN
            bravais(i,j) = 0d0
        ENDIF
    END DO
END DO


    call kp_gen(bravais, recbv, NSYMBZ, ROTMAT, NKX, NKY, NKZ, KPOIBZ, MAXMESH, NOFKS, BZKP, VOLCUB, TAUVBZ, .false., IBK)

    allocate(eigsave2(NOFKS(1),1,1))


OPEN(10, FILE='data/DOS.txt', STATUS = 'UNKNOWN')


call cpu_time(t_start)
DO i=1,NE
!ef=(-5.99d+00,0d+00)+12.d+00/NE*i
ef = cmplx(-4.0d0 + dble(i-1)*2.0d0, 0.0d0, 8)  ! E = -4, -2, 0, +2, +4 eV
    DO l=1,NOFKS(1)
	    kx = BZKP(1,l,1)
        ky = BZKP(2,l,1)
        kz = BZKP(3,l,1)

        evalns = my_dispersion(kx, ky, kz)
!	    write(6,*) evalns
        if(evalns.lt.DREAL(ef)) then
            EIGSAVE2(l,1,1)=1
        else
            EIGSAVE2(l,1,1)=0
        end if
    END DO

    call BISECTION1(1,1,1,1,NkX,NkY,NkZ,NSYMBZ,ef,IBK,recbv,eigsave2(:,1,1),rotmat, my_dispersion, sc_grad_mag)


END DO
call cpu_time(t_end)

contains

    ! Implementation of the abstract interface for the Simple Cubic lattice
    real(8) function my_dispersion(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi
        pi = 3.141592653589793115997963468544185161590576171875d+00
        my_dispersion = -(2*cos((kx-0.5)*2*pi) + 2*cos((ky-0.5)*2*pi) + 2*cos((kz-0.5)*2*pi))
    end function my_dispersion

    real(8) function sc_grad_mag(kx, ky, kz)
        real(8), intent(in) :: kx, ky, kz
        real(8) :: pi
        pi = 3.1415926535897932d0
        sc_grad_mag = 4.d0 * pi * sqrt(sin((kx-0.5d0)*2.d0*pi)**2 + &
                                       sin((ky-0.5d0)*2.d0*pi)**2 + &
                                       sin((kz-0.5d0)*2.d0*pi)**2)
        if (sc_grad_mag < 1.d-10) sc_grad_mag = 1.d-10
    end function sc_grad_mag

end program tb
