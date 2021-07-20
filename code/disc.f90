module disc1
!double precision :: a
!double precision, private :: b

!!demo vectorization
! double precision, allocatable, intent(out) :: vals(:)
! allocate(vals(10))
! vals = 1

contains 

subroutine doublelayer(t,s,val)
implicit double precision (a-h,o-z)
! Compute double layer potential K(t,s)
! Input:
!  t,s - parameters over [0,2*pi]
! Output:
!  val - K(t,s)
data pi /3.141592653589793238462643383279502884d0/
eps = epsilon(0.0d0)*100
call curve(t,xt,yt,dxt,dyt,dxxt,dyyt)

if (abs(t-s) <= eps) then
   val = (-dyt*dxxt+dxt*dyyt)/(4*pi*(dxt**2+dyt**2))
else 
   call curve(s,xs,ys,dxs,dys,dxxs,dyys)
   val = ((xs-xt)*dys-(ys-yt)*dxs)/(2*pi*((xt-xs)**2+(yt-ys)**2))
end if
end subroutine

subroutine curve(t,x,y,dx,dy,dx2,dy2)
implicit double precision (a-h,o-z)
! Curve parametrization
! Input: 
!  t - parameter over [0,2*pi]
! Output:
!  x,y,dx,dy,dx2,dy2

x   = 2.0d0*cos(t)
y   = sin(t)
dx  = -2.0d0*sin(t)
dy  = cos(t)
dx2 = -2.0d0*cos(t)
dy2 = -sin(t)
end subroutine 

end module

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! program disc

! use disc1  

! implicit double precision (a-h,o-z)
! double precision, allocatable :: xs(:), whts(:)

! !!checking
! ! t = 0.5d0
! ! s = 0.0d0
! ! call doublelayer(s,t,d1)
! ! call doublelayer(s,s,d2)
! ! print *,d1,d2

! do k = 1,6
!    n = 2**k !power of two
!    allocate(xs(n), whts(n))
!    data pi /3.141592653589793238462643383279502884d0/
!    call legendre_quad(n,xs,whts)
!    x = 0.0d0
!    approx = 0.0d0

!    !!approximate \int K(x,y) dy over [0,2*pi]
!    do i = 1,n
!       y = pi*(xs(i)+1)
!       call doublelayer(x,y,val)
!       approx = approx + val*pi*whts(i)
!    end do
!    print *,n,approx,abs(approx-0.5d0)

!    deallocate(xs,whts)
! end do

! end program
