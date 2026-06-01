module fermi_bisection_mod
    use mc33_core_mod
    use dispersion_interface_mod
    implicit none
contains

    SUBROUTINE BISECTION1(NAEZ,KVREL,ISPIN,Epoint,NBX,NBY,NBZ,NSYMBZ,E2,IBK,recbv,eigiso,rotmat, disp_fn, grad_fn)

!    use type,only:doublepr
    !use const,only:CONE,CI,PI,TPI,CZERO
    !use math_vec,only: scal,CROSS          !write own functions for these
 !   use fourier
  !  use inout,only:io


    implicit none
    procedure(dispersion_fn) :: disp_fn
    procedure(grad_mag_fn) :: grad_fn

    ! ************************************************************************
    !
    ! ---> calculate section points of an isoenergetic surface (e=EFERMI)
    !      with given lines in the reciprocal space
    !
    ! ------------------------------------------------------------------------
!***************************************
    INTEGER :: valog,IELAST
    REAL(8) :: RFCTOR,evalns
    !    arrays
    INTEGER,allocatable :: KAPOFQ(:),M2YOFQ(:)
    REAL(8),allocatable :: RBASIS(:,:)
    COMPLEX(8),allocatable :: RADJZ(:,:,:),RADJPLUS(:,:,:),RADJMINUS(:,:,:), &
         NORMINT(:,:,:,:),RADFZINT(:,:,:),RADPZINT(:,:,:), &
         RADPLUS(:,:,:),RADMINUS(:,:,:),RADZ(:,:,:), RMAT(:,:,:,:),ALPHAMAT(:,:,:,:), &
         NORMINTE(:,:,:,:),RMATE(:,:,:,:,:)
    COMPLEX(8), allocatable :: PZSQ(:)
    COMPLEX(8),allocatable::TINVBERRY(:,:,:,:,:), &
					TMATBERRY(:,:,:,:,:),CL(:),vel(:,:)
    COMPLEX(8),allocatable::berry_k_save(:,:,:),berry_r_save(:,:,:),berry_v_save(:,:,:)
    COMPLEX(8),allocatable::ginp_in(:,:,:,:)
    !	locals
    logical::berry,mass,wfct

!**************************************
    !    arguments
    INTEGER :: NAEZ,KVREL,ISPIN,Epoint,NBX,NBY,NBZ,NSYMBZ
    COMPLEX(8) :: E2
    COMPLEX(8), allocatable :: vf(:)
    REAL(8), allocatable :: vfl(:)
    !
    !     .arrays
    !
    INTEGER,allocatable :: LMOFQ(:,:),REFPOT(:),LOFLM(:),CLS(:),&
         NACLS(:),EZOA(:,:),ATOM(:,:),IBK(:,:,:,:)
    INTEGER:: EIGISO(:)
    REAL(8),allocatable :: CLE_JL(:,:),BZKP(:,:,:),&
         RECBV(:,:)
    REAL(8)::ROTMAT(3,3,3)

    COMPLEX(8),allocatable :: TINVLL(:,:,:,:),DELTALL(:,:,:),&
         CITPIRR(:,:),CITPIRCLS(:,:,:),GINP(:,:,:)
    COMPLEX(8),allocatable ::DELTALLQ(:,:),IMATRIX(:,:), &
         INVIMATRIX(:,:)
    !
    !     locals
    INTEGER :: Nlines,Nkpoints,Npoints,IA,I,LMMAX,LDIM,&
         NA,NE,NUM,IL,NL,test,IE,NDIMSPIN,M_DS,MOD_SPIR,&
         jx,jy,jz,jxm1,jym1,jzm1,jxp1,jyp1,jzp1,IN,&
         dmax,dmin,dmax1,dmin1,summin,summax,npp,nperm,j,pperm,i1,i2
    REAL(8)::TIMEI,TIMES,DCLOCK,FACA,FACE,FACM,length,minl,vol,kpscal,vol1, DOS, vel1, vel2, vel3
    COMPLEX(8)::masstr,massm(3,3),plot,berrytr
    !     array
    REAL(8) :: KP(3),KPA(3),KPE(3),KP1(3),KPC(3)
    REAL(8),allocatable :: EQ(:),FACrest(:)
    COMPLEX(8),allocatable :: GLLKE(:,:) !,massm(:,:,:)
    INTEGER,allocatable :: LINES(:,:),LINESW(:,:),IKP(:),perms(:,:),&
         ppoint(:),line2ksurf(:,:),nkp(:)
    INTEGER :: edges(2,12),linplane(12)
    REAL(8),allocatable :: LINESFAC(:,:),linek(:,:,:),&
         linekw(:,:,:),gq(:,:),poiplane(:,:),ksurf(:,:,:),kweight(:,:)
    CHARACTER*12 FILEN(5,2)
    LOGICAL :: found
    DATA FILEN / 'isodown.-   ','isodown ','isodown.+ ',&
         'isodown.++','isodown.--',&
         'isoup.-     ','isoup   ','isoup.+   ',&
         'isoup.++  ','isoup.--  ' /
    DATA edges / 1,2,2,3,3,4,4,1,1,5,2,6,3,7,4,8,5,6,6,7,7,8,8,5 /

    complex(8) :: CONE, CI, CZERO
    real(8) :: PI, TPI
    character(len=100) :: estring
    integer :: vert_edge(12)
    character(len=1) :: esc
    character(len=5) :: cyan_ansi, reset_ansi
    character(len=200) :: line_buf


    LDIM=SIZE(TINVLL,1)
    LMMAX=SIZE(GINP,2)
    M_DS=1
    NDIMSPIN=1
    MOD_SPIR=1

    CONE = (1,0)
    CI = (0,1)
    CZERO = (0,0)
    PI = 3.141592653589793115997963468544185161590576171875d+00
    TPI = 2*PI
    esc = achar(27)
    cyan_ansi = esc // '[36m'
    reset_ansi = esc // '[0m'
    write(estring, *) DREAL(E2)  ! Write Fermi energy to a string so file names can be changed.

!!$    OPEN(52, FILE='lines',STATUS='OLD')
    OPEN(12,FILE='data/'//trim(estring(1:8))//FILEN(EPOINT,ISPIN),STATUS='UNKNOWN')

    write(6,*) " "
    write(6,*) cyan_ansi // "+--------------------------------------------------------+" // reset_ansi
    write(6,*) cyan_ansi // "|  MC33 SOLVER RUN: Fermi Surface Reconstruction         |" // reset_ansi
    write(6,*) cyan_ansi // "+--------------------------------------------------------+" // reset_ansi
    write(line_buf, 1010) ISPIN, FILEN(EPOINT,ISPIN)
1010 FORMAT("  Spin Channel: ", I1, "            | Output File: ", A12)
    write(6,*) cyan_ansi // "|" // reset_ansi // line_buf(1:56) // cyan_ansi // "|" // reset_ansi
    WRITE(12,*) '!! HEAD' !(HEAD(I),I=1,76)

    !TIMES = DCLOCK()

    IF(ISPIN.eq.1) THEN
       OPEN(98,FILE='data/'//trim(estring(1:8))//"isosurf.dn",STATUS='UNKNOWN')
    ELSE
       OPEN(98,FILE='data/'//trim(estring(1:8))//"isosurf.up",STATUS='UNKNOWN')
    END IF

!!$    REWIND(52)
!!$    READ (52,*) NLINES,Nkpoints
    npoints=size(Eigiso,1)
!!$    IF (Nkpoints.NE.NPOINTS) THEN
!!$       write(io,*) 'N(',Nkpoints,').NE.NPOINTS(',NPOINTS,')'
!!$    END IF



    !###########################################################################
    !
    !    Find the lines of the mesh
    !
    !###########################################################################
    allocate(linesw(2,6*npoints),IKP(0:6),&
         linekw(3,2,6*npoints),gq(3,3))
    NLINES=0
    dmin=10000
    dmax=0
    vol1=0.0d+00

    gq(:,1)=  recbv(:,1)/nbx
    gq(:,2)=  recbv(:,2)/nby
    gq(:,3)=  recbv(:,3)/nbz

    ! the end points are equal to the start points by translation
    forall(jx=0:nbx, jy=0:nby, IBK(jx,jy,nbz,1).eq.0) IBK(jx,jy,nbz,1:2)=IBK(jx,jy,0,1:2)
    forall(jy=0:nby, jz=0:nbz, IBK(nbx,jy,jz,1).eq.0) IBK(nbx,jy,jz,1:2)=IBK(0,jy,jz,1:2)
    forall(jx=0:nbx, jz=0:nbz, IBK(jx,nby,jz,1).eq.0) IBK(jx,nby,jz,1:2)=IBK(jx,0,jz,1:2)

    !
    ! ---> find all lines
    !
    do jx = 0, nbx-1
       do jy = 0, nby-1
          do jz = 0, nbz-1

             IF (IBK(jx,jy,jz,1).eq.0) then
                write(6,*) "IBK eq 0 !!!!!!!!!!!!!!!!!!!!!"
                stop
             END IF

             IKP=0
             IKP(0)=IBK(jx,jy,jz,1)

             ! 6 Neighbors
             IKP(1)=IBK(jx+1,jy,jz,1)
             if(jx.gt.0) IKP(2)=IBK(jx-1,jy,jz,1)
             IKP(3)=IBK(jx,jy+1,jz,1)
             if(jy.gt.0) IKP(4)=IBK(jx,jy-1,jz,1)
             IKP(5)=IBK(jx,jy,jz+1,1)
             if(jz.gt.0) IKP(6)=IBK(jx,jy,jz-1,1)

             do IN=1,6
                IF ((ikp(in)).eq.0) cycle
                ! Take only the necessary lines
                IF ((EIGISO(IKP(0))-EIGISO(IKP(IN))).eq.0) cycle
                found=.false.
                do il=1,nlines
                   IF ( ( (IKP(0).eq.linesw(1,il)).and.(ikp(in).eq.linesw(2,il)) ).or. &
                        ( (IKP(0).eq.linesw(2,il)).and.(ikp(in).eq.linesw(1,il)) ) ) then
                      found=.true.
                      exit
                   END IF
                end do
                if (.not.found) then
                   nlines=nlines+1
                   ! The starting point has the smaller Eigiso
                   IF ((EIGISO(IKP(0))-EIGISO(IKP(IN))).gt.0) THEN
                      i1=2
                      i2=1
                   ELSE
                      i1=1
                      i2=2
                   END IF
                   linesw(i1,nlines)=IKP(0)
                   linesw(i2,nlines)=IKP(IN)
                   linekw(:,i1,il)=jx*gq(:,1)+jy*gq(:,2)+jz*gq(:,3)
                   select case( in )
                   case( 1 )
                      linekw(:,i2,il)=(jx+1)*gq(:,1)+jy*gq(:,2)+jz*gq(:,3)
                   case( 2 )
                      linekw(:,i2,il)=(jx-1)*gq(:,1)+jy*gq(:,2)+jz*gq(:,3)
                   case( 3 )
                      linekw(:,i2,il)=jx*gq(:,1)+(jy+1)*gq(:,2)+jz*gq(:,3)
                   case( 4 )
                      linekw(:,i2,il)=jx*gq(:,1)+(jy-1)*gq(:,2)+jz*gq(:,3)
                   case( 5 )
                      linekw(:,i2,il)=jx*gq(:,1)+jy*gq(:,2)+(jz+1)*gq(:,3)
                   case( 6 )
                      linekw(:,i2,il)=jx*gq(:,1)+jy*gq(:,2)+(jz-1)*gq(:,3)
                   end select
!!$                   IF (ABS(EIGISO(IKP(0))-EIGISO(IKP(IN))).gt.dimlin) dimlin=ABS(EIGISO(IKP(0))-EIGISO(IKP(IN)))
                   dmin=min(dmin,EIGISO(IKP(0)),EIGISO(IKP(IN)))
                   dmax=max(dmax,EIGISO(IKP(0)),EIGISO(IKP(IN)))
!!$                   IF (dimlin.gt.100) then
                   IF ((dmax-dmin).gt.100) then
                      write(6,*) "jx,jy,jz",jx,jy,jz
                      write(6,*) "Neighbor",in
                      write(6,*) nlines,dmax-dmin,dmin,dmax
                      write(6,*) IKP
                      STOP
                   endif
                end if
             end do
          end do
       end do
!       write(io,*) "END",jx,NLINES
    end do
!    write(6,*) 'martin'
    ! skip the big arrays and take smaller ones
    allocate(lines(2,nlines),linek(3,2,nlines))
    lines=linesw(:,:nlines)
    linek=linekw(:,:,:nlines)
    deallocate(linesw,linekw)


    WRITE(12,1001) NLINES,NPOINTS,DREAL(E2)
    write(line_buf, 1011) NLINES, NPOINTS, dmax-dmin
1011 FORMAT("  Mesh lines: ", I8, " | Points: ", I8, " | Dim: ", I5)
    write(6,*) cyan_ansi // "|" // reset_ansi // line_buf(1:56) // cyan_ansi // "|" // reset_ansi


    !###########################################################################
    !
    !    Calculate the bisection of the lines
    !
    !###########################################################################
!!$    write(io,*) NLINES,NPOINTS,dimlin
    ! write(6,*) NLINES,NPOINTS,dmax-dmin
!!$    if ((dmax-dmin).ne.1) then
!!$       write(io,*) "BISECTION1 DOES NOT WORK WITH multiple cuttings of the iso surface"
!!$       write(io,*) "HAS TO BE CHANGED !!!"
!!$       STOP
!!$    end if
!!$    allocate(linesFAC(NLINES,(dmax-dmin)),ksurf(3,nlines*(dmax-dmin)),&
!!$         kweight(nlines*(dmax-dmin)),line2ksurf(nlines,(dmax-dmin)))
    allocate(linesFAC(NLINES,dmin:(dmax-1)),ksurf(3,nlines,dmin:(dmax-1)),&
         kweight(nlines,dmin:(dmax-1)),line2ksurf(nlines,dmin:(dmax-1)),NKP(dmin:(dmax-1)))
    linesFAC=0.0d+00
    ksurf=0.0d+00
    kweight=0.0d+00
    line2ksurf=0
    nkp=0
    !
    ! --->    loop over lines
    !
    DO  IL = 1,NLINES
       !
!!$       READ (52,*) NL,IA,IE
       NL=IL
       IA=lines(1,il)
       IE=lines(2,il)

       !
       ! --->      number of eigenvalues at begin and end of the line
       !
!       NA = EIGISO(IA)

!       NE = EIGISO(IE)




       IF (EIGISO(IA).NE.EIGISO(IE)) THEN

	    IF(EIGISO(IA).lt.EIGISO(IE)) THEN
		NA=EIGISO(IA)
		NE=EIGISO(IE)
		 KPA(1:3)=linek(:,1,il)
		 KPE(1:3)=linek(:,2,il)
	    ELSE
		NA=EIGISO(IE)
		NE=EIGISO(IA)
		KPA(1:3)=linek(:,2,il)
	        KPE(1:3)=linek(:,1,il)
	    END IF

          !
          ! --->        line is crossed by the isoenergetic plane
          !
!!$          KPA(1:3)=BZKP(1:3,IA,1)
!!$          KPE(1:3)=BZKP(1:3,IE,1)
!          KPA(1:3)=linek(:,1,il)
!          KPE(1:3)=linek(:,2,il)

          Allocate(FACrest(0:NE-NA),&
               EQ(ABS(NA-NE)))

          EQ=0.0d+00
          FACrest=0.0d+00

          DO test=0,NE-NA-1
             FACA=FACrest(test)
             FACE=1.0d+00
             !changed from 10.d-8 to 10.d-10
             DO WHILE (ABS(FACA-FACE).gt.10.0d-15)

                FACM=(FACA+FACE)/float(2)

                KP(:)=KPA(:)+(KPE(:)-KPA(:))*FACM

                NUM=0

		!HERE you need to find the eigenvalues of the Hamiltonian

	            evalns = disp_fn(kp(1), kp(2), kp(3))
!                write(6,*) evalns

	if(evalns.lt.Dreal(e2)) then
	    num=1
	else
	    num=0
	end if

!                CALL EIGENFIND(GLLKE,NUM,.FALSE.)


!**********************************************
		IF(NUM.gt.NE) THEN
		    FACE=FACM
		ELSE IF(NUM.lt.NA) THEN
		    FACA=FACM
		ELSE IF(NUM-NA.gt.test) THEN
                   FACE=FACM
                ELSE
                   FACA=FACM
                END IF


                IF(NUM-NA.eq.test+1.and.FACrest(test+1).eq.0.0d+00)&
                     FACrest(test+1)=FACE
             END DO

	      IF(EIGISO(IA).lt.EIGISO(IE)) THEN
                EQ(test+1)=FACM
                linesfac(il,NA+test)=FACM
	      ELSE
		EQ(test+1)=1.0d+00-FACM
		linesfac(il,NA+test)=1.0d+00-FACM
	      END IF


             KP=(KPA+(KPE-KPA)*FACM)

rotmat=0.0d+00
rotmat(1,1,:)=1.0d+00
rotmat(2,2,:)=1.0d+00
rotmat(3,3,:)=1.0d+00

             DO I=1,1 !NSYMBZ
                kp1(1)=rotmat(1,1,i)*kp(1)+rotmat(1,2,i)*kp(2)+rotmat(1,3,i)*kp(3)
                kp1(2)=rotmat(2,1,i)*kp(1)+rotmat(2,2,i)*kp(2)+rotmat(2,3,i)*kp(3)
                kp1(3)=rotmat(3,1,i)*kp(1)+rotmat(3,2,i)*kp(2)+rotmat(3,3,i)*kp(3)
             END DO

             nkp(NA+test)=nkp(NA+test)+1
             ksurf(:,nkp(NA+test),NA+test)=KP
             line2ksurf(il,NA+test)=nkp(NA+test)

          END DO

          DEALLOCATE(EQ,FACrest)

       END IF                  ! (NA.NE.NE)

    END DO                    ! IL=1,NLINES



    deallocate(IKP,linek)

    !###########################################################################
    !
    !    Determine the lines which are cut by the isosurface in each box
    !
    !###########################################################################
    allocate(IKP(8),poiplane(3,12))
    !
    ! loop over boxes
    !
    do jx = 0, nbx-1
       do jy = 0, nby-1
          do jz = 0, nbz-1
             !
             ! 8 corners
             !
             jxp1=jx+1
             jyp1=jy+1
             jzp1=jz+1

             IKP(1)=IBK(jx,jy,jz,1)
             IKP(2)=IBK(jxp1,jy,jz,1)
             IKP(3)=IBK(jxp1,jyp1,jz,1)
             IKP(4)=IBK(jx,jyp1,jz,1)
             IKP(5)=IBK(jx,jy,jzp1,1)
             IKP(6)=IBK(jxp1,jy,jzp1,1)
             IKP(7)=IBK(jxp1,jyp1,jzp1,1)
             IKP(8)=IBK(jx,jyp1,jzp1,1)

             vert_edge(:) = 0 ! reset for each box

             DMIN1=10000
             DMAX1=0
             DO I=1,8
                IF (EIGISO(IKP(I)).lt.DMIN1) DMIN1=EIGISO(IKP(I))
                IF (EIGISO(IKP(I)).gt.DMAX1) DMAX1=EIGISO(IKP(I))
             END DO
             IF ((DMAX1-DMIN1).gt.0) THEN
                !
                ! The iso surface cuts the box
                !
                !12 lines
                DO I1=DMIN,DMAX-1,2
                npp=0
                DO I=1,12
                   IF ((ABS(EIGISO(IKP(edges(1,i)))-EIGISO(IKP(edges(2,i)))).gt.0))  then
                   IF ( (min(EIGISO(IKP(edges(1,i))),EIGISO(IKP(edges(2,i)))).le.I1).and.&
                        (max(EIGISO(IKP(edges(1,i))),EIGISO(IKP(edges(2,i)))).gt.I1) ) THEN
                      npp=npp+1

                      Do il=1,2
                         select case( edges(il,i) )
                         case( 1 )
                            kp=jx*gq(:,1)+jy*gq(:,2)+jz*gq(:,3)
                         case( 2 )
                            kp=(jx+1)*gq(:,1)+jy*gq(:,2)+jz*gq(:,3)
                         case( 3 )
                            kp=(jx+1)*gq(:,1)+(jy+1)*gq(:,2)+jz*gq(:,3)
                         case( 4 )
                            kp=jx*gq(:,1)+(jy+1)*gq(:,2)+jz*gq(:,3)
                         case( 5 )
                            kp=jx*gq(:,1)+jy*gq(:,2)+(jz+1)*gq(:,3)
                         case( 6 )
                            kp=(jx+1)*gq(:,1)+jy*gq(:,2)+(jz+1)*gq(:,3)
                         case( 7 )
                            kp=(jx+1)*gq(:,1)+(jy+1)*gq(:,2)+(jz+1)*gq(:,3)
                         case( 8 )
                            kp=jx*gq(:,1)+(jy+1)*gq(:,2)+(jz+1)*gq(:,3)
                         end select
                         IF (IL.eq.1) KPA=KP
                         KPE=KP
                      End Do

                      do il=1,nlines
                         IF ( ( (IKP(edges(1,i)).eq.lines(1,il)).and.&
                                (IKP(edges(2,i)).eq.lines(2,il)) ).or.&
                              ( (IKP(edges(1,i)).eq.lines(2,il)).and.&
                                (IKP(edges(2,i)).eq.lines(1,il)) ) ) then
                            found=.true.
                            exit
                         END IF
                      end do
                      linplane(npp)=il
                      if (.not.found) write(6,*) 'ERROR'



                      IF ( (IKP(edges(1,i)).eq.lines(1,il)).and.&
                           (IKP(edges(2,i)).eq.lines(2,il)) ) THEN
                           poiplane(:,npp)=KPA(:)+(KPE(:)-KPA(:))*linesfac(il,i1)
                            vert_edge(npp) = edge_index(edges(1,i),edges(2,i))
                      END IF

                      IF ( (IKP(edges(1,i)).eq.lines(2,il)).and.&
                           (IKP(edges(2,i)).eq.lines(1,il)) ) THEN
                           poiplane(:,npp)=KPE(:)+(KPA(:)-KPE(:))*linesfac(il,i1)
                            vert_edge(npp) = edge_index(edges(1,i),edges(2,i))
                      END IF

                    END IF
                    END IF
                END DO

                IF(npp.eq.0) cycle ! Not this isoplane cuts the box

                IF ((npp.lt.3).or.(npp.gt.9)) THEN
                   write(6,*) "something went terribly wrong !!!"
                   write(6,*) "npp=",npp
                   write(6,'(8I5)') (EIGISO(IKP(I)),I=1,8)
                   STOP
                END IF

                !
                ! write out triangles
                !
!		write(io,*) 'here1'
                IF (npp.eq.3) THEN
		    plot = CZERO
                    vel1 = grad_fn(poiplane(1,1), poiplane(2,1), poiplane(3,1))
                    vel2 = grad_fn(poiplane(1,2), poiplane(2,2), poiplane(3,2))
                    vel3 = grad_fn(poiplane(1,3), poiplane(2,3), poiplane(3,3))
                    write(98,'(4D15.6)') poiplane(:,1), vel1     ! GNUplot FORMATTING
		    write(98,'(4D15.6)') poiplane(:,1), vel1
		    write(98,*)
		    write(98,'(4D15.6)') poiplane(:,2), vel2
		    write(98,'(4D15.6)') poiplane(:,3), vel3
		    write(98,*)
		    write(98,*)
                   ! Calculate Area
                   KPA=poiplane(:,2)-poiplane(:,1)
                   KPE=poiplane(:,3)-poiplane(:,1)
                   KP=CROSS(KPA,KPE)
                   vol=sqrt(dot_product(KP,KP))/2.d0
                   vol1=vol1+vol
                   Do i=1,3
                      kweight(line2ksurf(linplane(i),i1),i1)=&
                           kweight(line2ksurf(linplane(i),i1),i1)+vol/3.d0
                   END Do
                ELSE
                    call marching_cube(poiplane,jx,jy,jz,gq,e2,npp,vert_edge,grad_fn)
                   ! plane has more than 3 corners
                   ! cut into triangles
                   ! write(6,*) npp
                   ! find point of centre
!===============================================================!
!======== THIS IS WHERE MARCHING CUBES NEEDS TO HELP ===========!
                   KP=0.d0
                   do i=1,npp
                      KP=KP+poiplane(:,i)
                   end do
                   KP=KP/npp
                   ! find shortest track
                   nperm = factorial(npp)                                      ! GAMMA(x+1)=x!
                   allocate(PERMS(1:npp,1:nperm),ppoint(1:npp))

		    forall (i=1:npp) ppoint(i)=i
                   call permutation(ppoint,PERMS)
                   minl=1.d+10
                   DO j=1,nperm
                      ! length of this permutation
                      length=0.d0
                      kpa=poiplane(:,perms(1,j))
                      DO i=2,npp

                         kpe(1:3)=poiplane(1:3,perms(i,j))

			             KPC=kpa-kpe
			             kpscal=dot_product(kpc,kpc)
                 		 length=length + dsqrt(kpscal)
                         kpa=kpe
                      END DO
                      kpe=poiplane(:,perms(1,j))
                      length=length + dsqrt(dot_product(kpa-kpe,kpa-kpe))

                      if (length.lt.minl) then
                         minl=length
                         pperm=j
                      end if
                   END DO

                   ! pperm is the shortest one
!====================================================================!
                   DO i=1,npp-1
                                         ! Calculate Area
                      KPA=kp-poiplane(:,perms(i,pperm))
                      KPE=kp-poiplane(:,perms(i+1,pperm))
                      KP1=CROSS(KPA,KPE)
                      vol=sqrt(dot_product(KP1,KP1))/2.d0
                      vol1=vol1+vol
                      kweight(line2ksurf(linplane(perms(i,pperm)),i1),i1)=&
                           kweight(line2ksurf(linplane(perms(i,pperm)),i1),i1)+vol/2.d0
                      kweight(line2ksurf(linplane(perms(i+1,pperm)),i1),i1)=&
                           kweight(line2ksurf(linplane(perms(i+1,pperm)),i1),i1)+vol/2.d0
                   END DO
            ! WRITE FINAL TRIANGLE
                   KPA=kp-poiplane(:,perms(npp,pperm))
                   KPE=kp-poiplane(:,perms(1,pperm))
                   KP1=CROSS(KPA,KPE)
                   vol=sqrt(dot_product(KP1,KP1))/2.d0
                   vol1=vol1+vol
                   kweight(line2ksurf(linplane(perms(npp,pperm)),i1),i1)=&
                        kweight(line2ksurf(linplane(perms(npp,pperm)),i1),i1)+vol/2.d0
                   kweight(line2ksurf(linplane(perms(1,pperm)),i1),i1)=&
                        kweight(line2ksurf(linplane(perms(1,pperm)),i1),i1)+vol/2.d0
                   deallocate(PERMS,ppoint)
             END IF
             END DO ! I1=DMIN,DMAX-1
             END IF ! ((DMAX-DMIN).gt.0)


         end do ! jz
!                   write(6,*) perms

       end do   ! jy
    end do      ! jx
                    line_buf = "  Calculated bisections & topological edge-crossings.   "
                    write(6,*) cyan_ansi // "|" // reset_ansi // line_buf(1:56) // cyan_ansi // "|" // reset_ansi
    CLOSE(12)
    CLOSE(52)
    open(12,FILE='data/'//trim(estring(1:8))//'kmesh.isosurf',STATUS='UNKNOWN')
    vol=0.0d+00
    dos=0.0d+00
    DO I1=DMIN,DMAX-1
	DO i=1,nkp(i1)
    	    write(12,'(4F16.12)') ksurf(:,i,i1),kweight(i,i1)
    	    vol=vol+kweight(i,i1)
    	    dos=dos+kweight(i,i1)/grad_fn(ksurf(1,i,i1), ksurf(2,i,i1), ksurf(3,i,i1))
	END DO
    END DO
    close(12)
    write(line_buf, 1012) vol
1012 FORMAT("  Surface Area:  ", F12.6, " (2pi/a)^2                 ")
    write(6,*) cyan_ansi // "|" // reset_ansi // line_buf(1:56) // cyan_ansi // "|" // reset_ansi
   ! open(12,FILE='kmesh.3d',STATUS='OLD')
!    DO i=1,npoints
!       write(12,'(3F16.12)') BZKP(:,I,1)
!    END DO
    !close(12)



    write(line_buf, 1013) real(e2), DOS
1013 FORMAT("  E_F: ", F8.3, " eV | DOS: ", ES12.4, " states/eV       ")
    write(6,*) cyan_ansi // "|" // reset_ansi // line_buf(1:56) // cyan_ansi // "|" // reset_ansi
    write(6,*) cyan_ansi // "+--------------------------------------------------------+" // reset_ansi
    write(6,*) " "
    write(10, *) real(e2), DOS


    deallocate(lines,IKP,linesFAC,poiplane,gq,ksurf,kweight,line2ksurf)


1001 FORMAT(2I8,1p,d20.10,'   NLINES NPOINTS EFERMI')
1002 FORMAT(I8,I4,/,(I4,f16.12))
1003 FORMAT(' time in ISPIN loop :',f12.2,' sec.')


  END SUBROUTINE BISECTION1

recursive subroutine permutation(a,b)
    implicit none
    INTEGER :: A(:),B(:,:)
    INTEGER,allocatable :: C(:)
    integer i,n,k,j

    n=size(A)
    IF(n.eq.1) THEN
       B(1,1)=A(1)
    ELSE
       allocate(C(n-1))
       k=factorial(n-1)
       DO i=1,n
     forall (j=(i-1)*k+1:i*k) B(1,j)=A(i)
     IF(i.eq.1) THEN

               C(i:n-1)=A(i+1:n)
          ELSE IF(i.eq.n) THEN
                    C(1:i-1)=A(1:i-1)
      ELSE
           C(1:i-1)=A(1:i-1)
              C(i:n-1)=A(i+1:n)
      END IF
    !          C(1:i-1)=A(1:i-1)
    !          C(i:n-1)=A(i+1:n)
          call permutation(C,B(2:n,(i-1)*k+1:i*k))
       END DO
       deallocate(C)
    END IF
end subroutine permutation

recursive function factorial(n) result(fac)
    implicit none
    integer n,fac
    if (n.eq.0) then
       fac=1
    else
       fac=factorial(n-1)*n
    endif
end function factorial

FUNCTION CROSS(a, b)


    REAL(8), DIMENSION(3) :: cross
    REAL(8), DIMENSION(3), INTENT(IN) :: a, b

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
END FUNCTION CROSS

end module fermi_bisection_mod
