!! calculate the cartesian of V and VI in JCP, 115, 174 [2001]

program main
  implicit none
  real*8,parameter :: ang=0.5291772083d0, ev=27.21138505d0, pi=dacos(-1.d0)
  real*8 :: c(3,4)
  real*8,external :: angle3,angle
  real*8 :: t,aa

  !! structure VI
  c=0.d0
  c(1,2)= 1.0160d0*cos(106.7d0/2.d0/180.d0*pi)
  c(2,2)= 1.0160d0*sin(106.7d0/2.d0/180.d0*pi)
  c(1,3)= 1.0160d0*cos(106.7d0/2.d0/180.d0*pi)
  c(2,3)=-1.0160d0*sin(106.7d0/2.d0/180.d0*pi)

  print*,angle(c(:,2)-c(:,1),c(:,3)-c(:,1))

  aa=120.d0
  t=61.2230440303916d0
  do while(aa.gt.106.7d0)
  c(1,4)=-1.0160d0*cos(t/180.d0*pi)
  c(2,4)=0.d0
  c(3,4)= 1.0160d0*sin(t/180.d0*pi)
  aa=angle(c(:,2)-c(:,1),c(:,4)-c(:,1))
  print*,t,aa
  t=t+1.d-13
  enddo

  write(*,'(3f20.12)') c

  !! structure V
  c=0.d0
  c(1,2)= 0.9864d0*cos(105.7d0/2.d0/180.d0*pi)
  c(2,2)= 0.9864d0*sin(105.7d0/2.d0/180.d0*pi)
  c(1,3)= 0.9864d0*cos(105.7d0/2.d0/180.d0*pi)
  c(2,3)=-0.9864d0*sin(105.7d0/2.d0/180.d0*pi)

  print*,angle(c(:,2)-c(:,1),c(:,3)-c(:,1))

  aa=120.d0
  t=62.664285444677d0
  do while(aa.gt.106.1d0)
  c(1,4)=-1.1986d0*cos(t/180.d0*pi)
  c(2,4)=0.d0
  c(3,4)=1.1986d0*sin(t/180.d0*pi)
  aa=angle(c(:,2)-c(:,1),c(:,4)-c(:,1))
  print*,t,aa
  t=t+1.d-12
  enddo

  write(*,'(3f20.12)') c

end program

function angle3(r1,r2,r3)
  ! triangle r1 r2 r3, angle3 face to r3
  implicit none
  real*8,intent(in) :: r1,r2,r3
  real*8 :: angle3
  real*8 :: tmp
  real*8,parameter :: pi=3.14159265358979323846d0
  tmp=(r1**2+r2**2-r3**2)/(2.d0*r1*r2)
  tmp=max(tmp,-1.d0)
  tmp=min(tmp,1.d0)
  angle3=dacos(tmp)*180.d0/pi
  return
end function

function angle(v1,v2)
  implicit none
  real*8,intent(in) :: v1(3),v2(3)
  real*8 :: angle
  real*8 :: tmp
  real*8,parameter :: pi=3.14159265358979323846d0
  real*8 :: vtmp(3),r1,r2,r3
  vtmp=v1-v2
  r1=dot_product(v1,v1)
  r2=dot_product(v2,v2)
  r3=dot_product(vtmp,vtmp)
  tmp=(r1+r2-r3)/(2.d0*dsqrt(r1*r2))
  tmp=max(tmp,-1.d0)
  tmp=min(tmp,1.d0)
  angle=dacos(tmp)*180.d0/pi
  return
end function

