module brillouin_zone_mod
    implicit none

contains

    subroutine bzirr3d(ix,iy,iz,nsym,rotmat,recbv,bravais,volbz,&
           lsurf, nkp,QP,wtkp,ibk1)
        implicit none
        ! arguments
        integer ix,iy,iz,nsym,nkp
        INTEGER,allocatable :: ibk1(:,:,:,:)
        real(8) :: volbz
        real(8) :: rotmat(:,:,:), recbv(:,:),bravais(:,:),&
                               QP(:,:),wtkp(:)
        logical lsurf
        ! local
        integer :: i,jx,jy,jz,iws,iwt,is,ja,jb,jc,NK !n,k
        real(8) :: v1,x,y,z,a,b,c,ox,oy,oz
        INTEGER,allocatable :: nkxyz1(:)
        real(8),allocatable :: rq(:,:),gq(:,:)
        logical, allocatable :: ibk(:,:,:)
        !
        allocate(ibk(0:ix,0:iy,0:iz),nkxyz1(3),rq(3,3),gq(3,3))

        nkxyz1(1)=ix
        nkxyz1(2)=iy
        nkxyz1(3)=iz

        NK=ix*iy*iz

        do i = 1, 3
           gq(1:3,i)=  recbv(1:3,i)/nkxyz1(i)
           rq(1:3,i)=bravais(1:3,i)*nkxyz1(i)
        enddo

        ibk=.true.
        nkp=0
        iws=0

        do jx = 0, ix-1
           do jy = 0, iy-1
              do jz = 0, iz-1
                 if(ibk(jx,jy,jz)) then
                    nkp=nkp+1

                    iwt=0
                    x=jx*gq(1,1)+jy*gq(1,2)+jz*gq(1,3)
                    y=jx*gq(2,1)+jy*gq(2,2)+jz*gq(2,3)
                    z=jx*gq(3,1)+jy*gq(3,2)+jz*gq(3,3)

                    nsym=1
                    rotmat=0.0d+00
                    rotmat(1,1,:)=1.0d+00
                    rotmat(2,2,:)=1.0d+00
                    rotmat(3,3,:)=1.0d+00
                    do is=1,nsym 

                       ox=rotmat(1,1,is)*x+rotmat(1,2,is)*y+rotmat(1,3,is)*z
                       oy=rotmat(2,1,is)*x+rotmat(2,2,is)*y+rotmat(2,3,is)*z
                       oz=rotmat(3,1,is)*x+rotmat(3,2,is)*y+rotmat(3,3,is)*z

                       a=ox*rq(1,1)+oy*rq(2,1)+oz*rq(3,1)
                       b=ox*rq(1,2)+oy*rq(2,2)+oz*rq(3,2)
                       c=ox*rq(1,3)+oy*rq(2,3)+oz*rq(3,3)

                      ja=nint(a)
                      jb=nint(b)
                      jc=nint(c)
                      if(abs(ja-a)+abs(jb-b)+abs(jc-c).gt.1.e-3) then
                         write(6,'(''ERROR in bzirr3d!'')')
                         write(6,'(3(i3))')ja,jb,jc
                         write(6,'(3(f12.10))')a,b,c
                         write(6,*) is,z,oz
                         write(6,*) x,y,z
                         write(6,*) ox,oy,oz
                      endif

                       ja=mod(ja,ix)
                       if(ja.lt.0) ja=ja+ix
                       jb=mod(jb,iy)
                       if(jb.lt.0) jb=jb+iy
                       jc=mod(jc,iz)
                       if(jc.lt.0) jc=jc+iz
                       if(ibk(ja,jb,jc)) then
                          ibk(ja,jb,jc)=.false.
                          iwt=iwt+1
                          IBK1(ja,jb,jc,1)=NKP
                          IBK1(ja,jb,jc,2)=is
                       endif
                    enddo
                    IF(.not.allocated(IBK1)) THEN
                       do is=1,3
                          if (sqrt(x*x+y*y+z*z).gt.sqrt((x+recbv(1,is))**2+(y+recbv(2,is))**2+(z+recbv(3,is))**2)) then
                             x=x+recbv(1,is)
                             y=y+recbv(2,is)
                             z=z+recbv(3,is)
                          endif
                          if (sqrt(x*x+y*y+z*z).gt.sqrt((x-recbv(1,is))**2+(y-recbv(2,is))**2+(z-recbv(3,is))**2)) then
                             x=x-recbv(1,is)
                             y=y-recbv(2,is)
                             z=z-recbv(3,is)
                          endif
                       enddo
                    END IF

                    QP(1,nkp)=x
                    QP(2,nkp)=y
                    QP(3,nkp)=z
                    wtkp(nkp)=dble(iwt)/NK

                    iws=iws+iwt

                 endif
              enddo
           enddo
        enddo
        v1 = 0.d0
        do i=1,nkp
           wtkp(i) = wtkp(i)*volbz/dfloat(nsym)
           v1 = v1 + wtkp(i)*dfloat(nsym) ! check volume
        end do

        deallocate(ibk,nkxyz1,rq,gq)
    end subroutine bzirr3d

    SUBROUTINE kp_gen(BRAVAIS,RECBV,NSYMBZ,ROTMAT,NBX,NBY,NBZ,&
       KPOIBZ,MAXMESH,NOFKS,BZKP,VOLCUB,TAUVBZ,LINTERFACE,IBK)
    implicit none
    ! Arguments
    integer :: NSYMBZ,NBX,NBY,NBZ,KPOIBZ,MAXMESH
    real(8) :: tauvbz
    logical LINTERFACE
    real(8) :: BRAVAIS(:,:),RECBV(:,:),ROTMAT(:,:,:)
    integer,allocatable :: NOFKS(:),IBK(:,:,:,:)
    real(8),allocatable :: BZKP(:,:,:),VOLCUB(:,:)
    ! Local
    integer :: ifile,LMESH,I
    real(8) :: volbz
    real(8),allocatable :: BZKPw(:,:,:),VOLCUBw(:,:)

    volbz = dabs(recbv(1,1)*recbv(2,2)-recbv(1,2)*recbv(2,1))
    TAUVBZ=1.d0/VOLBZ

    ifile=0
    kpoibz=NBX*NBY*NBZ
    allocate (BZKPw(3,KPOIBZ,MAXMESH),&
              VOLCUBw(KPOIBZ,MAXMESH),&
              NOFKS(MAXMESH))
    BZKPw=0.d0
    VOLCUBw=0.d0

    DO  LMESH=1,MAXMESH
       IF (LMESH.GT.1) THEN
          NBX = (NBX)/1.4
          NBY = (NBY)/1.4
          NBZ = (NBZ)/1.4
          IF (NBX.eq.0) NBX  = 1
          IF (NBY.eq.0) NBY  = 1
          IF (NBZ.eq.0) NBZ  = 1
       END IF
       call BZIRR3D(NBX,NBY,NBZ,NSYMBZ,ROTMAT,RECBV,BRAVAIS,&
            volbz,LINTERFACE,&
            NOFKS(LMESH),BZKPw(:,:,LMESH),VOLCUBw(:,LMESH),ibk)
    END DO

    kpoibz=maxval(NOFKS)
    allocate (BZKP(3,KPOIBZ,MAXMESH),VOLCUB(KPOIBZ,MAXMESH))
    BZKP=BZKPw(:,1:KPOIBZ,:)
    VOLCUB=VOLCUBw(1:KPOIBZ,:)

    deallocate(BZKPw,VOLCUBw)
    end subroutine kp_gen

end module brillouin_zone_mod
