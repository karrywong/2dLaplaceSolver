program matrix
use disc1

implicit double precision (a-h,o-z)
double precision, allocatable :: amatr(:,:), x(:), y(:)
integer, allocatable :: ipiv(:) !tests II & III
double precision, allocatable :: fval(:), xq(:), whtq(:) !test III
data pi /3.141592653589793238462643383279502884d0/

!!!test I - Verify A*1 = 1, with number of quadrature nodes doubled in each run
! print *,"Compute \|A*1 - 1\| "
! do k=1,7
!    ! number of quadrature
!    n = 2**k
!    allocate(amatr(n,n),x(n),y(n))
!    x = 1.0d0 !one-vector

!    call intop(n,amatr)
!    !print*,n,amatr
!    y = matmul(amatr,x)
!    !print *,y

!    print *,"n=",n,", 2-norm error:",norm2(y-x)
!    deallocate(amatr,x,y)
! end do

!!!test II - Verify inv(A)*1 = 1, with number of quadrature nodes doubled in each run
! print *, NEW_LINE('a'), "Compute \|inv(A)*1 - 1\| "
! do k = 1,7 
!    n=2**k

!    allocate(amatr(n,n),x(n),y(n))
!    y = 1.0d0
!    call intop(n,amatr)

!    nrhs=1 !number of RHS in b
!    lda=n
!    ldb=n
!    allocate(ipiv(n))
!    ! Solve Ax = y
!    call dgesv(n, nrhs, amatr, lda, ipiv, y, ldb, info)
!    !print *,"soln = inv(A)*y: ", y

!    x = 1.0d0
!    print *,"n=",n,", 2-norm error:", norm2(x-y)
!    deallocate(amatr,x,y)
!    deallocate(ipiv)
! end do

!!! test III - solve \sigma for RHS = log|p-p0|, then verify through integral formula
n=500
allocate(xq(n), whtq(n), fval(n))
call legendre_quad(n,xq,whtq)
!Chosen fixed point outside \Omega
x2=4.0d0
y2=-4.0d0
do i=1,n
   t = pi*(xq(i)+1.0d0)
   call curve(t,xt,yt,dxt,dyt,dx2t,dy2t)
   call func(xt,yt,x2,y2,val)
   fval(i) = val !RHS vector
end do

allocate(amatr(n,n),x(n),y(n),ipiv(n))
call intop(n,amatr)
nrhs=1 !number of RHS in b
lda=n
ldb=n
! Solve Ax = y
y = fval
call dgesv(n, nrhs, amatr, lda, ipiv, y, ldb, info)
!print*,y
!print*,info

!!!compute u using integral formula and compare it with fval0
u=0
xu=0
yu=0
do i=1,n
   s = pi*(xq(i)+1.0d0)
   call doublelayer2(xu,yu,s,val)
   u  = u + val*y(i)*pi*whtq(i)
end do
!u = 0.5*u/pi

print *,"u:",u
call func(xu,yu,x2,y2,fval0)
print *,"fval0:",fval0
end program


subroutine intop(n,amatr)
! Return discretized integral operator for BVP on \Omega
! Input: n - size of square matrix amatr
! Output: amatr - 0.5*I plus matrix of double layer potential [K(xi, yj)]

use disc1
implicit double precision (a-h,o-z)
double precision :: amatr(n,n)
data pi /3.141592653589793238462643383279502884d0/
double precision, allocatable :: xs(:), whts(:)
allocate(xs(n), whts(n))
call legendre_quad(n,xs,whts)

do i=1,n
   xi = pi*(xs(i)+1.0d0)
   do j=1,n
      yj = pi*(xs(j)+1.0d0)
      call doublelayer(xi,yj,val)
      amatr(i,j) = val*pi*whts(j)
   end do
end do

do i=1,n
   amatr(i,i) = amatr(i,i)+0.5d0
end do
end subroutine

subroutine func(x1,y1,x2,y2,val)
! Return value of RHS function f=log|p1-p2| on boundary \partial\Omega
! Input: x1,y1,x2,y2 - p1=(x1,y1), p2=(x2,y2)
! Output: val - values of f
implicit double precision (a-h,o-z)
val = log(sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)))
end subroutine

subroutine doublelayer2(x1,y1,s,val)
use disc1
implicit double precision (a-h,o-z)
! Compute double layer potential K(x1,y1,s)
data pi /3.141592653589793238462643383279502884d0/
call curve(s,xs,ys,dxs,dys,dxxs,dyys)
val = ((xs-x1)*dys-(ys-y1)*dxs)/(2*pi*((x1-xs)**2+(y1-ys)**2))
end subroutine
