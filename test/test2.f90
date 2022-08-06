!Stationary points from JCP, 115, 174 [2001]
!    II:   -0.0176
!   III:   -0.0303
!    IV:    0.2335
!     V:    0.1925
!    VI:    0.0942
! H+H2O:   -0.6952

program test
  implicit none
  real*8,parameter :: ang=0.5291772083d0, ev=27.21138505d0, pi=dacos(-1.d0)
  real*8 :: c(3,4),r(6,1),epes(1)

  write(*,'(a44)') "Stationary points from JCP, 115, 174 [2001]"

  c=reshape([ 0.d0,0.d0,0.d0, &
             -0.9703d0,0.d0,0.d0, &
             2.7476d0,0.d0,0.d0, &
             2.7476d0+0.7433d0,0.d0,0.d0,0.d0 ],[3,4])
  call cart2dist(c,r); r=r/ang
  call dpotdriva(1,r,epes)
! call OH3PES('NN1',1,r,epes)
  write(*,'(a8,f10.4)') "II:",epes*ev

  c=reshape([ 2.2504d0+0.9704d0,0.d0,0.d0, &
              2.2504d0,0.d0,0.d0, &
              0.d0,0.7438d0/2.d0,0.d0, &
              0.d0,-0.7438d0/2.d0,0.d0 ],[3,4])
  call cart2dist(c,r); r=r/ang
  call dpotdriva(1,r,epes)
! call OH3PES('NN1',1,r,epes)
  write(*,'(a8,f10.4)') "III:",epes*ev

  c=reshape([ 0.d0,0.d0,0.d0, &
              1.356d0,0.d0,0.d0, &
              0.9700d0*cos(96.5d0/180.d0*pi),0.9700d0*sin(96.5d0/180.d0*pi),0.d0, &
              1.356d0-0.8192d0*cos(163.6d0/180.d0*pi),0.8192d0*sin(163.6d0/180.d0*pi),0.d0 ],[3,4])
  call cart2dist(c,r); r=r/ang
  call dpotdriva(1,r,epes)
! call OH3PES('NN1',1,r,epes)
  write(*,'(a8,f10.4)') "IV:",epes*ev

  c=reshape([ 0.d0,0.d0,0.d0, &
              0.9864d0*cos(105.7d0/2.d0/180.d0*pi), 0.9864d0*sin(105.7d0/2.d0/180.d0*pi),0.d0, &
              0.9864d0*cos(105.7d0/2.d0/180.d0*pi),-0.9864d0*sin(105.7d0/2.d0/180.d0*pi),0.d0, &
             -1.1986d0*cos(62.664285444677d0/180.d0*pi),0.d0,1.1986*sin(62.664285444677d0/180.d0*pi) ], [3,4])
  call cart2dist(c,r); r=r/ang
  call dpotdriva(1,r,epes)
! call OH3PES('NN1',1,r,epes)
  write(*,'(a8,f10.4)') "V:",epes*ev

  c=reshape([ 0.d0,0.d0,0.d0, &
              1.0160d0*cos(106.7d0/2.d0/180.d0*pi), 1.0160d0*sin(106.7d0/2.d0/180.d0*pi),0.d0, &
              1.0160d0*cos(106.7d0/2.d0/180.d0*pi),-1.0160d0*sin(106.7d0/2.d0/180.d0*pi),0.d0, &
             -1.0160d0*cos(61.2230440303916d0/180.d0*pi),0.d0,1.0160*sin(61.2230440303916d0/180d0*pi) ], [3,4])
  call cart2dist(c,r); r=r/ang
  call dpotdriva(1,r,epes)
! call OH3PES('NN1',1,r,epes)
  write(*,'(a8,f10.4)') "VI:",epes*ev

  c=reshape([ 0.d0,0.d0,0.d0, &
              0.9589d0*cos(104.2d0/2.d0/180.d0*pi), 0.9589d0*sin(104.2d0/2.d0/180.d0*pi),0.d0, &
              0.9589d0*cos(104.2d0/2.d0/180.d0*pi),-0.9589d0*sin(104.2d0/2.d0/180.d0*pi),0.d0, &
              5.d0,5.d0,5.d0 ],[3,4])
  call cart2dist(c,r); r=r/ang
  call dpotdriva(1,r,epes)
! call OH3PES('NN1',1,r,epes)
  write(*,'(a8,f10.4)') "H+H2O:",epes*ev

  stop
end program

subroutine cart2dist(c,r)
  implicit none
  integer,parameter :: natom=4
  real*8,intent(in) :: c(3,natom)
  real*8,intent(out) :: r(natom*(natom-1)/2)
  integer i,j,k,id
  id=0
  do i=1,natom-1
  do j=i+1,natom
    id=id+1
    r(id)=0.d0
    do k=1,3
    r(id)=r(id)+(c(k,i)-c(k,j))**2
    enddo
    r(id)=sqrt(r(id))
  enddo
  enddo
  return
end subroutine
