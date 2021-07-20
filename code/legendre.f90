! program test_legendre
! implicit double precision (a-h, o-z)
! double precision, allocatable :: xs(:), whts(:)

! !!!!test legendre_poly!!!
! !n = 3
! !x = 0.5d0
! !call legendre_poly(n,x,val,der)
! !!!!Reference values!!!
! !val0 = -0.4375d0
! !der0 = 0.375d
! !print *,val-val0
! !print *,der-der0 

! n = 10
! allocate(xs(n), whts(n))
! call legendre_quad(n, xs, whts)
! !Validate degree of exactness (2n-1) of xs and whts
! do k = 0,2*n-1
!    approx = 0 
!    do i = 1,n
!       approx = approx + whts(i) * xs(i)**k
!    end do
!    print*,"degree",k,"over [-1,1], approx. val.:",approx
! end do

! end program


subroutine legendre_quad(n, xs, whts)
implicit double precision (a-h, o-z)
double precision :: xs(n), whts(n)

! Return the n-point Gauss-Legendre quadrature on (-1,1)
! 
! Input parameter:
! n - length of the desired quadrature
!
! Output parameter: 
! xs - an array containing quadrature nodes
! whts - an array containing quadrature weights

!tolerance for Netwon method
eps0 = 2.0d0*epsilon(0.0d0)
pi = acos(-1.0d0)

!Find roots of P_n(x) via Newton method
!Chebyshev node as initial guess
do i=1,n
   t = (n-i+1.0d0)/(n+1.0d0)
   x0 = cos(pi * t)
   x = x0

   !Newton method - abs(delta/x) < eps0
   delta = 1.0d0
   do while (abs(delta) > abs(x)*eps0)
      call legendre_poly(n,x,val,der)
      delta = val/der
      x = x - delta
      !print *,"                  ", delta
   end do

   xs(i) = x
   !print *,i,x0,x
end do

!Compute weights
do i=1,n
   !w = 2/((1-xs^2) * P'_n(xs)^2)
   call legendre_poly(n,xs(i),val,der)
   whts(i) = 2/((1-xs(i)**2) * der**2)
   !print *,i,whts(i) 
end do
end subroutine

subroutine legendre_poly(n,x,val,der)
implicit double precision (a-h, o-z)
! Evaluate legendre polynomials and its derivative 
! Input parameters:
!  n - degggree of poly
!  x - the point at which to evaluate
!
! Output
!  val - the value of P_n(x)
!  der - the value of P_n' (x)

double precision :: pols(0:n)

!corner cases n=0 & n =1
if (n==0) then
   val = 1.0d0
   der = 0.0d0
   return 
endif

if (n==1) then
   val = x
   der = 1.0d0
   return
endif

pols(0) = 1.0d0
pols(1) = x

!Three terms recurrence relation
!(n+1)P_{n+1}(z) - (2n+1)z * P_n(z) + n P_{n-1}(z) = 0
!P_{n+1}(z) = - (2n+1) / (n+1) z * P_n(z) - n/(n+1)  P_{n-1}(z) 
do i=1,n-1
   pols(i+1) = (2*i+1.0d0) / (i+1.0d0) * x * pols(i) - i/(i+1.0d0) * pols(i-1)
end do  
val = pols(n)

!P_n'(x) = 1/(z^2-1) * n* (z P_n(z) - P_{n-1}(z))
der = 1/(x**2-1.0d0) * n * (x* pols(n) - pols(n-1))

end subroutine
