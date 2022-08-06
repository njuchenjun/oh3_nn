#define pesdir "/Users/chenjun/pes/OH3-2013-CXZ/pes/"

! ====================================================================
c input ntot: number of geometries
c input r(6,ntot): distances of o h h h, bohr
c output v(ntot): potential energy in hartree
c     vv0=76.84440690d0 ! uccsd(t)-f12a/avtz  oh+h2
c     vv0=76.86995525d0 ! uccsd(t)-f12a/avtz  h+h2o

      subroutine spotdriva(ntot,r,v)
      implicit none
      integer,intent(in) :: ntot
      real*4,intent(in)  :: r(6,ntot)
      real*4,intent(out) :: v(ntot)
      real*8 rx(6),vx,dv(6)
      integer i,j,idv
      idv=0
      v=0.d0
      do i=1,ntot
        rx(1:6)=r(1:6,i)
        call nninter2(rx,vx,idv,dv)  ! v (hartree), v(oh+h2)=0.0
        v(i)=vx
      enddo
      return
      end subroutine spotdriva

      subroutine dpotdriva(ntot,r,v)
      implicit none
      integer,intent(in) :: ntot
      real*8,intent(in)  :: r(6,ntot)
      real*8,intent(out) :: v(ntot)
      real*8 rx(6),vx,dv(6)
      integer i,j,idv
      idv=0
      v=0.d0
      do i=1,ntot
        rx(1:6)=r(1:6,i)
        call nninter2(rx,vx,idv,dv)  ! v (hartree), v(oh+h2)=0.0
        v(i)=vx
      enddo
      return
      end subroutine dpotdriva

      subroutine dpotdrivax(ntot,r,v) ! v(h+h2o)=0.0
      implicit none
      integer,intent(in) :: ntot
      real*8,intent(in)  :: r(6,ntot)
      real*8,intent(out) :: v(ntot)
      call dpotdriva(ntot,r,v)
      v=v-76.84440690D0+76.86995525D0
      return
      end subroutine dpotdrivax

! --------------------------------------------------------------------
! ---- simulate a neural network with two hidden layers
! ----      n0-n1-n2-1
! ---- blas routines used in this subroutine: dgemv, ddot
      subroutine nsim(r0,v,idv,dv,n0,n1,n2,rg,w1,b1,w2,b2,w3,b3,vg)
      implicit none
      integer,intent(in) :: n0,n1,n2,idv
      real*8,intent(in)  :: r0(n0),rg(2,n0),vg(2)
      real*8,intent(in)  :: w1(n0,n1),b1(n1)
      real*8,intent(in)  :: w2(n1,n2),b2(n2)
      real*8,intent(in)  :: w3(n2),b3
      real*8,intent(out) :: v,dv(n0)
      integer :: i,j,k
      real*8  :: r(n0),rgg(n0),vgg,ax(n1),bx(n2)
      real*8  :: dvtm,rtmp,rt1(n1),rt2(n2)
      real*8,external :: ddot
      v=0.d0
      r=r0
      vgg=vg(2)-vg(1)
      ! mapminmax [-1,1]
      do i=1,n0
        rgg(i)=rg(2,i)-rg(1,i)
        r(i)=2.d0*(r(i)-rg(1,i))/rgg(i)-1.d0
      end do
      ! 1st layer
      rt1=b1
      call dgemv('t',n0,n1,1.d0,w1,n0,r,1,1.d0,rt1,1)
      ax=dtanh(rt1)
      ! 2nd layer
      rt2=b2
      call dgemv('t',n1,n2,1.d0,w2,n1,ax,1,1.d0,rt2,1)
      bx=dtanh(rt2)
      ! output layer
      v=b3+ddot(n2,w3,1,bx,1)
      !reverse map
      v=vgg*(v+1.d0)/2+vg(1)
      if(idv.ne.1) return
      ! calculate first derivatives, dv(i)=dv/dr(i)
      dv=0.d0
      do i=1,n0
        do k=1,n2
          dvtm=0.d0
          do j=1,n1
            dvtm=dvtm+w2(j,k)*w1(i,j)*(1-ax(j)**2)
          enddo
          dv(i)=dv(i)+w3(k)*dvtm*(1-bx(k)**2)
        enddo
        dv(i)=dv(i)*vgg/rgg(i)
      enddo
      return
      end subroutine nsim
! --------------------------------------------------------------------

      module nnparam2

      implicit none
      save

      real*8,parameter :: pi=3.14159265358979323846d0
      integer,parameter :: ndim=6
      integer,parameter :: bondidx(4,4)=[-1,1,2,3,1,-1,4,5,2,4,-1,6,3,5,6,-1]
!     character*99 :: pesdir='/home/chenjun/oh3pes/7.seg/a-NN1/'
      character*99 :: n22f='n22-r.txt' ! 22
      character*99 :: n40f='n40-2w.txt' ! 40
      character*99 :: n13f='n13-l.txt' ! 13
      integer,parameter :: n22s1=25,n22s2=25
      integer,parameter :: n40s1=50,n40s2=50
      integer,parameter :: n13s1=20,n13s2=20

      real*8 :: n22w1(ndim,n22s1),n22b1(n22s1)
      real*8 :: n22w2(n22s1,n22s2),n22b2(n22s2)
      real*8 :: n22w3(n22s2),n22b3
      real*8 :: n22in(2,ndim),n22out(2)
      character*8 :: n22fcn(2)

      real*8 :: n40w1(ndim,n40s1),n40b1(n40s1)
      real*8 :: n40w2(n40s1,n40s2),n40b2(n40s2)
      real*8 :: n40w3(n40s2),n40b3
      real*8 :: n40in(2,ndim),n40out(2)
      character*8 :: n40fcn(2)

      real*8 :: n13w1(ndim,n13s1),n13b1(n13s1)
      real*8 :: n13w2(n13s1,n13s2),n13b2(n13s2)
      real*8 :: n13w3(n13s2),n13b3
      real*8 :: n13in(2,ndim),n13out(2)
      character*8 :: n13fcn(2)

      contains
      subroutine nninit
      implicit none
      integer :: i,j,k

c     print*,trim(pesdir)//trim(n22f)
      open(100,file=trim(pesdir)//trim(n22f),status='old',action='read')
      read(100,*)n22w1,n22b1,n22w2,n22b2,n22w3,n22b3,n22in,n22out,n22fcn
      close(100)

c     print*,trim(pesdir)//trim(n40f)
      open(100,file=trim(pesdir)//trim(n40f),status='old',action='read')
      read(100,*)n40w1,n40b1,n40w2,n40b2,n40w3,n40b3,n40in,n40out,n40fcn
      close(100)

c     print*,trim(pesdir)//trim(n13f)
      open(100,file=trim(pesdir)//trim(n13f),status='old',action='read')
      read(100,*)n13w1,n13b1,n13w2,n13b2,n13w3,n13b3,n13in,n13out,n13fcn
      close(100)

      n22out=n22out/27.2116d0
      n40out=n40out/27.2116d0
      n13out=n13out/27.2116d0

      end subroutine nninit

      end module nnparam2

c===========================================================
      subroutine nninter2(rr0,vx,dx,dv)
c  dx=0 calculate energy only
c  dx=1 calculate 1st derivative
      use nnparam2
      implicit none
      real*8 :: rr0(ndim),vx,dv(ndim)
      integer :: i,j,k,itmp,idx,dx
      integer :: iloc(1),sortidx(4)
      real*8 :: rr(ndim), rtmp,vtmp,df
      real*8 :: v13,v22,v40
      real*8 :: dv13(ndim),dv22(ndim),dv40(ndim)
      real*8 :: f13,f22,f40,fss
      real*8 :: z(6)
      real*8 :: logsig,x
      integer :: flag0,flag23,flag12
      integer,save :: init=0

      real*8 :: weight_h2o,r_h2o(3),v_h2o
      real*8,external :: potoh,pothh

      logsig(x)=1.d0/(1.d0+dexp(-x))
      if (init.eq.0) then
        init=1
        call nninit()
      end if
      flag0=0;flag23=0;flag12=0

c bond length permutation
c input: rr0(1:6)=r(oh1),r(oh2),r(oh3),r(h1h2),r(h1h3),r(h2h3)
c permute h1,h2,h3 to oh1.le.oh2.le.oh3
      sortidx(1)=1
      iloc=minloc(rr0(1:3))
      sortidx(2)=iloc(1) +1
      j=3
      do i=2,4
        if(i == sortidx(2)) cycle
        sortidx(j)=i
        j=j+1
      end do
      if(rr0(sortidx(3)-1) > rr0(sortidx(4)-1)) then
        itmp=sortidx(3)
        sortidx(3)=sortidx(4)
        sortidx(4)=itmp
      end if
      k=0
      do i=1,3
        do j=i+1,4
          k=k+1
          rr(k)=rr0( bondidx( sortidx(i),sortidx(j) ) )
        end do
      end do
c end bond length permutation

990   continue
      
      f22=logsig(50.d0*(rr(2)-3.5d0))
      f13=logsig(50.d0*(rr(3)-6.4d0))*(1.d0-f22)
      f40=1.d0-f22-f13
 
      v13=0.d0
      v22=0.d0
      v40=0.d0
      
      fss=1.0e-8
      
      if(f22.gt.fss) call nsim(rr,v22,dx,dv22,ndim,n22s1,n22s2,n22in,n22w1,n22b1,n22w2,n22b2,n22w3,n22b3,n22out)
      if(f13.gt.fss) call nsim(rr,v13,dx,dv13,ndim,n13s1,n13s2,n13in,n13w1,n13b1,n13w2,n13b2,n13w3,n13b3,n13out)
      if(f40.gt.fss) call nsim(rr,v40,dx,dv40,ndim,n40s1,n40s2,n40in,n40w1,n40b1,n40w2,n40b2,n40w3,n40b3,n40out)
      
      vx=f22*v22+f13*v13+f40*v40
!!    return ! no permutation symmetry

      if(flag0.eq.1) goto 991

      if (rr(3)-rr(2).lt.0.1.and.rr(3).ge.rr(2).and.flag23.eq.0) then
         rtmp=rr(3);rr(3)=rr(2);rr(2)=rtmp;
         rtmp=rr(4);rr(4)=rr(5);rr(5)=rtmp;
         vtmp=vx
         flag23=1
         goto 990
      endif

      if (flag23.eq.1) then
         df=0.5*(1+sin((rr(2)-rr(3))*5*pi))
         vx=vx*(1-df)+vtmp*df
         rtmp=rr(3);rr(3)=rr(2);rr(2)=rtmp;
         rtmp=rr(4);rr(4)=rr(5);rr(5)=rtmp;
      endif
991   flag0=1

      if (rr(2)-rr(1).lt.0.1.and.rr(2).ge.rr(1).and.flag12.eq.0) then
         rtmp=rr(1);rr(1)=rr(2);rr(2)=rtmp;
         rtmp=rr(6);rr(6)=rr(5);rr(5)=rtmp;
         vtmp=vx
         flag12=1
         goto 990
      endif

      if (flag12.eq.1) then
         df=0.5*(1+sin((rr(1)-rr(2))*5*pi))
         vx=vx*(1-df)+vtmp*df
         rtmp=rr(1);rr(1)=rr(2);rr(2)=rtmp;
         rtmp=rr(6);rr(6)=rr(5);rr(5)=rtmp;
      endif

!! for H+H2O asymptotic, switch to H2O's PES, 2014-8-13
      if(rr(1).gt.rr(2)) rr=rr([2,1,3,4,6,5])
      if(rr(1).gt.rr(3)) rr=rr([3,2,1,6,5,4])
      if(rr(2).gt.rr(3)) rr=rr([1,3,2,5,4,6])

      weight_h2o=logsig(25.d0*(5.d0-rr(2)))*logsig(25.d0*(rr(3)-12.d0))
      if(weight_h2o.gt.1.d-6) then !
        r_h2o=rr([1,2,4])
        ! h2o pes, a.u.
        call nnpes_h2o(r_h2o,v_h2o)
        v_h2o=v_h2o+76.84440690D0-76.86995525D0
!       if(abs(v_h2o-vx).gt.1.d-3) write(*,'(a4,10f10.5)') "Vh2o",rr
!       write(9,'(a4,10f10.5)') "Vh2o",rr,[vx,v_h2o,vx-v_h2o]*27.2114d0
      else
        v_h2o=vx
      endif
      vx=vx*(1.d0-weight_h2o)+v_h2o*weight_h2o

!! for H2+OH asymptotic, not-the-final-solution
      if(rr(2).gt.30.d0) then
        vx=potoh(rr(1))+pothh(rr(6))
        vx=vx/27.2116d0
      endif

      return
      end

! --------------------------------------------------------------------
!       O a1     a2 H3
!     r1 \---------/ r2  d
!         H2  z   H4
      subroutine jacob22_dist(r1,r2,z,a1,a2,d,r)
      implicit none
      real*8 r1,r2,z,a1,a2,d,r(6),c(3,4),m(4)
      real*8 tmp(3)
      real*8,parameter :: torad=dacos(-1.d0)/180.d0
      integer i,j,k,id
      data m/29160.D0,1837.15D0,1837.15D0,1837.15D0/
      r=0.d0
      c=0.d0
!     step 1: linear
!     H2-O----H3-H4 -->x, (COM of OH)=(000)
!       | <-z-> |
      c(1,1)=r1*m(2)/(m(1)+m(2))
      c(1,2)=-r1*m(1)/(m(1)+m(2))
      c(1,3)=z-r2*m(4)/(m(3)+m(4))
      c(1,4)=z+r2*m(3)/(m(3)+m(4))
!     step 2: rotate O-H2 around z unticlockwisely
!       O a1            cos   -sin  0
!     r1 \-----H3-H4    sin   cos   0
!         H2 z-> |      0     0     1
      do i=1,2
        tmp(1)=c(1,i)*dcos(a1*torad)-c(2,i)*dsin(a1*torad)
        tmp(2)=c(1,i)*dsin(a1*torad)+c(2,i)*dcos(a1*torad)
        c(1,i)=tmp(1)
        c(2,i)=tmp(2)
      enddo
!     step 3: translate -z in x axis
      do i=1,4
        c(1,i)=c(1,i)-z
      enddo
!     step 4: rotate H3-H4 around z axis clockwisely
!       O a1     a2 H3        cos   sin   0
!     r1 \---------/ r2       -sin  cos   0
!         H2  z   H4          0     0     1
      do i=3,4
        tmp(1)=c(1,i)*dcos(a2*torad)+c(2,i)*dsin(a2*torad)
        tmp(2)=-c(1,i)*dsin(a2*torad)+c(2,i)*dcos(a2*torad)
        c(1,i)=tmp(1) 
        c(2,i)=tmp(2)
      enddo
!     step 5: rotate H3-H4 around x axis clockwisely
!     1     0     0
!     0     cos   sin
!     0     -sin  cos
      do i=3,4
        tmp(2)=c(2,i)*dcos(d*torad)+c(3,i)*dsin(d*torad)
        tmp(3)=-c(2,i)*dsin(d*torad)+c(3,i)*dcos(d*torad)
        c(2,i)=tmp(2)
        c(3,i)=tmp(3)
      enddo

!     calculate distance
      id=0
      do i=1,3
      do j=i+1,4
        id=id+1
        r(id)=0.d0
        do k=1,3
        r(id)=r(id)+(c(k,i)-c(k,j))**2
        enddo
        r(id)=dsqrt(r(id))
      enddo
      enddo

      return
      end

! --------------------------------------------------------------------
C POTOH, O-H UCCSD(T)-F12A/AVTZ ENERGY IN EV
      FUNCTION POTOH(R)
      REAL*8 R,POTOH
      REAL*8 X,W1(5),B1(5),W2(5),B2,RA,RB,VA,VB
      INTEGER I
      DATA W1,B1,W2,B2,RA,RB,VA,VB/
     &   6.4673602067509854E+000,
     &  -1.3623917123812161E+000,
     &   1.9692891975248517E+000,
     &  -3.0842317137227271E+000,
     &   2.8753386118410593E+000,
     &  -6.2995066766419541E+000,
     &   1.1338119520445564E-001,
     &   7.3524071060910179E-001,
     &  -1.9242209086024227E+000,
     &   4.6341234161417368E+000,
     &   4.4487144796737021E-003,
     &  -7.5921384799742464E-001,
     &   4.4369462093719220E-001,
     &  -8.7516096362305562E-002,
     &  -4.6326649659746060E+001,
     &   4.5791134903018026E+001,
     &   1.2500000000000000E+000,
     &   4.1000000000000000E+000,
     &   2.5300000000000002E-004,
     &   5.2943569999999998E+000/
      SAVE W1,B1,W2,B2,RA,RB,VA,VB

      POTOH=0.D0
      X=2*(R-RA)/(RB-RA)-1
      DO I=1,5
        POTOH=POTOH+W2(I)*TANH(B1(I)+W1(I)*X)
      ENDDO
      POTOH=POTOH+B2
      POTOH=(POTOH+1)*(VB-VA)/2+VA
      RETURN
      END FUNCTION POTOH
C===========================================================
C POTHH, H-H UCCSD(T)-F12A/AVTZ ENERGY IN EV
      FUNCTION POTHH(R)
      REAL*8 R,POTHH
      REAL*8 X,W1(5),B1(5),W2(5),B2,RA,RB,VA,VB
      INTEGER I
      DATA W1,B1,W2,B2,RA,RB,VA,VB/
     &   1.6264981845926245E+000,
     &  -3.6036352665856217E+000,
     &  -1.4262429547486353E+000,
     &  -4.5389257076092262E+000,
     &  -1.0370365421774828E+001,
     &  -5.9327868896216229E-001,
     &   6.9006605801506959E-001,
     &  -9.0839635069711533E-001,
     &  -5.9790603462261309E+000,
     &  -1.3007124917554380E+001,
     &  -6.3219690048134394E-002,
     &   6.3788100794878316E-003,
     &  -1.0432656434523606E+000,
     &   2.0101618557956517E+001,
     &   4.0421850398756270E+001,
     &   5.9403385445541353E+001,
     &   6.0000000000000000E-001,
     &   6.1000000000000000E+000,
     &   1.8000000000000000E-005,
     &   1.1043358000000000E+001/
      SAVE W1,B1,W2,B2,RA,RB,VA,VB

      POTHH=0.D0
      X=2*(R-RA)/(RB-RA)-1
      DO I=1,5
        POTHH=POTHH+W2(I)*TANH(B1(I)+W1(I)*X)
      ENDDO
      POTHH=POTHH+B2
      POTHH=(POTHH+1)*(VB-VA)/2+VA
      RETURN
      END FUNCTION POTHH
! ====================================================================
! Training set : 5344 points, exchange r1 and r2 to 10688 points
! /home/chenjun/h2o/6-pes-f12a-avtz/nnfit/fit2/W06.txt 3-30-30-1
!  Total rmse_06=   0.0577 meV; train=   0.0180; validation=   0.1742

! provide bond in bohr and angle in degree, return energy in Hartree
! r1=r(O-H1), r2=r(O-H2), theta=angle(H1-O-H2)
! V(1.81176168,1.81176168,104.45471213)=0.0
! --------------------------------------------------------------------
!     program test
!     implicit none
!     real*8 r(3),v
!     r=[1.81176168d0,1.81176168d0,0.6d0]
!     do while(r(3).lt.3.6d0)
!     call nnpes_h2o(r,v)
!     print*,r(3),v*27.2114d0
!     r(3)=r(3)+0.01d0
!     enddo
!     end
! --------------------------------------------------------------------
      module h2omod
      implicit none
      save
!     character*99 :: n01fa='/home/chenjun/oh3pes/7.seg/a-NN1/nh2o_2_6.txt'
      character*99 :: n01fa='nh2o_2_6.txt'
      integer,parameter :: ndim=3
      integer,parameter :: n01s1=30,n01s2=30

!     01 double hidden layers
      real*8 :: n01w1a(ndim,n01s1),n01b1a(n01s1),
     &n01w2a(n01s1,n01s2),n01b2a(n01s2),n01w3a(n01s2),
     &n01b3a,n01ina(2,ndim),n01outa(2)
      contains
      subroutine nninit
      open(1,file=trim(pesdir)//trim(n01fa),status='old');
      read(1,*)n01w1a,n01b1a,n01w2a,n01b2a,n01w3a,n01b3a,n01ina,n01outa;
      close(1)
      return
      end subroutine nninit
      end module h2omod
! --------------------------------------------------------------------
! provide bond distances of O H1 H2 in bohr, return energy in Hartree
! r: O-H1, O-H2, H1-H2
      subroutine nnpes_h2o(r_in,v_out)
      implicit none
      integer,parameter :: kind_use=8
      real(kind=kind_use),intent(in) :: r_in(3)
      real(kind=kind_use),intent(out) :: v_out
      real*8 :: r1(3),r2(3),v1,v2

      r1([1,2,3])=r_in([1,2,3])
      r2([1,2,3])=r_in([2,1,3])
      call nnfit_h2o(r1,v1)
      call nnfit_h2o(r2,v2)
      v_out=(v1+v2)/2.D0/27.21138505D0

      return
      contains

! --------------------------------------------------------------------
      subroutine nsim(rin,vx,ndim,neu1,neu2,nxin,nxw1,nxb1,
     &                nxw2,nxb2,nxw3,nxb3,nxout)
!     for two-hidden-layer feedforward neural network, I-J-K-1
!     two blas routines used: dgemv, ddot
      implicit none
      integer,intent(in) :: ndim,neu1,neu2
      real*8,intent(in) :: rin(ndim),nxin(2,ndim),
     &                     nxw1(ndim,neu1),nxb1(neu1),
     &                     nxw2(neu1,neu2),nxb2(neu2),
     &                     nxw3(neu2),nxb3,nxout(2)
      real*8,intent(out) :: vx
      integer :: i,j,k
      real*8 :: r(ndim),ax(neu1),bx(neu2),vrange,rrange(ndim),
     &          rtmp,rt1(neu1),rt2(neu2)
      real*8,external :: ddot

      vx=0.d0
      r=rin
      vrange=nxout(2)-nxout(1)
      ! mapminmax [-1,1]
      do i=1,ndim
        rrange(i)=nxin(2,i)-nxin(1,i)
        r(i)=2.d0*(r(i)-nxin(1,i))/rrange(i)-1.d0
      end do

      ! 1st layer
      rt1=nxb1
      call dgemv('t',ndim,neu1,1.d0,nxw1,ndim,r,1,1.d0,rt1,1)
      do j=1,neu1
        ax(j)=dtanh(rt1(j))
      end do

      ! 2nd layer
      rt2=nxb2
      call dgemv('t',neu1,neu2,1.d0,nxw2,neu1,ax,1,1.d0,rt2,1)
      do k=1,neu2
        bx(k)=dtanh(rt2(k))
      end do

      ! output layer
      vx=nxb3+ddot(neu2,nxw3,1,bx,1)

      !reverse map
      vx=vrange*(vx+1.d0)/2.d0+nxout(1)

      return
      end subroutine nsim

! --------------------------------------------------------------------
      subroutine nnfit_h2o(r0,vx)
      use h2omod
      implicit none
      real*8,intent(in) :: r0(ndim)
      real*8,intent(out) :: vx
      real*8 :: r(ndim),tmp,va,vb,vc,vd,ve,vf
      integer,save :: init=0
      integer i,j,k
      if (init.eq.0) then
        init=1
        call nninit
      endif

      r(1:ndim)=r0(1:ndim)
      call nsim(r,va,ndim,n01s1,n01s2,n01ina,n01w1a,n01b1a,
     &n01w2a,n01b2a,n01w3a,n01b3a,n01outa)

      vx=va
      return
      end subroutine nnfit_h2o

! --------------------------------------------------------------------
      end subroutine nnpes_h2o

! --------------------------------------------------------------------
