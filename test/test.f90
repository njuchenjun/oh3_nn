program test
  implicit none
  real*4 :: c(3,4),irc,e0,r(6,1),epes(1)

  open(101,file='OH3-AB-IRC',status='old',action='read')
  do while(.true.)
    read(101,'(2f12.6,5x,3f14.8,3x,3f14.8,3x,3f14.8,3x,3f14.8)',end=901) irc,e0,c(1:3,1),c(1:3,2),c(1:3,3),c(1:3,4)
    call cart2dist(c,r)
    call spotdriva(1,r,epes)
    write(801,*) irc,e0,epes*27.2114d0
  enddo
  901 close(101)

  open(102,file='OH3-EX-IRC',status='old',action='read')
  do while(.true.)
    read(102,'(2f12.6,5x,3f14.8,3x,3f14.8,3x,3f14.8,3x,3f14.8)',end=902) irc,e0,c(1:3,1),c(1:3,2),c(1:3,3),c(1:3,4)
    call cart2dist(c,r)
    call spotdriva(1,r,epes)
    write(802,*) irc,e0,epes*27.2114d0
  enddo
  902 close(102)

  stop
end program

subroutine cart2dist(c,r)
  implicit none
  integer,parameter :: natom=4
  real*4,intent(in) :: c(3,natom)
  real*4,intent(out) :: r(natom*(natom-1)/2)
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
