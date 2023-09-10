!*******************************************************************************
!> author: Oscar Mojica
!
!  Modern Fortran version of the driver for ZHANKS function, whose original 
!  version is listed in the appendix of the paper:
!
!  Walter L. Anderson, (1979), "Numerical integration of related Hankel transforms 
!  of orders 0 and 1 by adaptive digital filtering," GEOPHYSICS 44: 1287-1305.
!  https://doi.org/10.1190/1.1441007
!
!  Note that the original COMPLEX FUNCTION ZHANKS(args) was converted to 
!  function eval_hankel(args) result(zhanks)
!  The syntax for complex kernels functions C1-C5 also was changed, i.e.
!  COMPLEX FUNCTION C1(G) became function cmplx_fun1(g) result(c1)
!
!### Original header function:
!
!--DRIVER PROGRAM TO TEST SUBPROGRAM ZHANKS.
!
!--BY W.L.ANDERSON, U.S. GEOLOCICAL SURVEY, DENVER, COLORADO.

! REQUIRES ONLY A PRINTER FILE (UNIT 06).

! TEST DATA ARE PRESTORED IN DATA STATEMENTS FOR EVALUATING
! HANKEL TRANSFORMS OF ORDER N=0 OR 1 OF THE FORM
! INTEGRAL FROM 0 TO INFINITY OF FUN(G)*JN(G*B)*DG, B.GT.0.

! THE FOLLOWING (REAL) TEST INTEGRALS ARE EVALUATED USING ZHANKS

! INTEGRAL N KERNEL FUNCTION  RELATED KERNEL EXACT RESULT OF INTEGRAL
! ======== = ===============  ============== =========================
!    1     0 C1=G*EXP(-G*G)   C1 (USE NEW=1) 0.5*EXP(-.25*B*B)
!    2     1 C2=G*G*EXP(-G*G) C2=G*C1        0.25*EXP(-.25*B*B)
!    3     1 C3=EXP(-G)       C3 (USE NEW=1) (R-1)/(B*R),R=SQRT(1+B*B)
!    4     1 C4=G*EXP(-2.*G)  C4=G*C3**2     B/SQRT(4.+B*B)**3
!    5     0 C5=EXP(-2.*G)    C5=C4/G        1./SQRT(4.+B*B)

!--SUBPROGRAMS CALLED ARE ZHANKS,C1,C2,C3,C4,C5, AND SAVER.

!--NOTE USE OF DIFFERENT COMPUTERS AND WORD LENGTHS WILL FRODUCE
!  SLIGHTLY DIFFERENT ROUND-OFF AND POSSIBLY DIFFERENT NUMBER OF
!  DIRECT FUNCTION CALLS (NF) USED IN ZHANKS. HOWEVER, THE FILTERED
!  RESULTS SHOULD AGREE REASONABLY WELL WITH EXACT RESULTS WITH
!  RESPECT TO THE TOLERANCE (TOL) USED IN THESE EXAMPLES.

!--THE USER SHOULD INSERT A CALL TO SUPPRESS EXPONENT UNDERFLOW MESSAGES
!  FOR THE MACHINE SYSTEM BEING USED (SEE NOTE(1) IN ZHANKS COMMENTS).

!=======================================================================
program driver

  use tem_kinds
  use hankel_module
  use iso_fortran_env, only: stdout => output_unit

  implicit none
  
  complex(dp) :: z
  integer :: i, j, nf
  real(dp) :: droot4, filt, exact, abserr, relerr, bb
  character(len=44), parameter :: fmt1 = "(a17,12x,a3,10x,a5,11x,a8,8x,a9,7x,a9,6x,7a)"
  character(len=36), parameter :: fmt2 = "(2(5x,i1),2e13.5,4e16.8,2x,i3,3x,i1)"

  real(dp) :: b(9)=[.0001_dp,.001_dp,.005_dp,.01_dp,.05_dp,.1_dp,.5_dp,1._dp,2._dp]
  real(dp) :: tol(9)=[0.1e-5_dp,0.1e-5_dp,0.1e-6_dp,0.1e-7_dp,0.1e-7_dp, &
                      0.1e-7_dp,0.1e-7_dp,0.1e-8_dp,0.1e-8_dp]
  integer :: n(5)=[0,1,1,1,0]
  integer :: newj(5)=[1,0,1,0,0]

  write(stdout,'(a30)')'TEST results for zhanks driver'
  write(stdout,'(a)')''
  write(stdout,fmt1)' INTEGRAL  N    B','TOL','EXACT','FILTERED','ABS.ERROR','REL.ERROR','NF  NEW'
  write(stdout,'(a)')''

  do i=1,9
    bb=b(i)*b(i)
	! get integral 1 (use new=1)
    z=eval_hankel(0,b(i),cmplx_fun1,tol(i),nf,1)
    filt=z%re
    exact=.5_dp*dexp(-.25_dp*bb)
    abserr=dabs(filt-exact)
    relerr=abserr/exact
    j=1
    write(stdout,fmt2)j,n(j),b(i),tol(i),exact,filt,abserr,relerr,nf,newj(j)
    call saver(1,1)
  ! get integral 2 (new=0 in eval_hankel)
    z=eval_hankel(1,b(i),cmplx_fun2,tol(i),nf,0)
    filt=z%re
    exact=.25_dp*b(i)*dexp(-.25_dp*bb)
    abserr=dabs(filt-exact)
    relerr=abserr/exact
    j=2
    write(stdout,fmt2)j,n(j),b(i),tol(i),exact,filt,abserr,relerr,nf,newj(j)
	! get integral 3 (use new=1)
    z=eval_hankel(1,b(i),cmplx_fun3,0.1_dp*tol(i),nf,1)
    filt=z%re
    exact=(1.0_dp-1.0_dp/dsqrt(1.0_dp+bb))/b(i)
    abserr=dabs(filt-exact)
    relerr=abserr/exact
    j=3
    write(stdout,fmt2)j,n(j),b(i),0.1_dp*tol(i),exact,filt,abserr,relerr,nf,newj(j)
    call saver(1,2)
	! get integral 4 (new=0 in eval_hankel)
    z=eval_hankel(1,b(i),cmplx_fun4,0.1*tol(i),nf,0)
    filt=z%re
    droot4=dsqrt(4.0_dp+bb)
    exact=b(i)/droot4**3
    abserr=dabs(filt-exact)
    relerr=abserr/exact
    j=4
    write(stdout,fmt2)j,n(j),b(i),0.1_dp*tol(i),exact,filt,abserr,relerr,nf,newj(j)
    call saver(-1,1)
	! get integral 5 (new=0 in eval_hankel)
    z=eval_hankel(0,b(i),cmplx_fun5,0.1_dp*tol(i),nf,0)
    filt=z%re
    exact=1.0_dp/droot4
    abserr=dabs(filt-exact)
    relerr=abserr/exact
    j=5
    write(stdout,fmt2)j,n(j),b(i),0.1_dp*tol(i),exact,filt,abserr,relerr,nf,newj(j)
    write(stdout,'(a)')''
  enddo
  
  contains
  
!  This subroutine modifies the fsave array as follows
!  fsave(k)=gsave(k)**i * fsave(k)**j
!  for k=1,nsave
!
!  input parameters (i,j) may be negative, zero, or positive integers
!
!  As function above, the subroutine and the kernel functions after that,
!  are provided in the appendix of the paper:
! 
!  Walter L. Anderson, (1979), "Numerical integration of related Hankel transforms 
!  of orders 0 and 1 by adaptive digital filtering," GEOPHYSICS 44: 1287-1305.
!  https://doi.org/10.1190/1.1441007

  subroutine saver(i,j)
  
  implicit none
  integer :: i,j,k
  
  do k=1,nsave
    fsave(k)=cmplx(gsave(k)**i,0._dp, kind=dp)*(fsave(k)**j)
  enddo
	
  end subroutine saver 

  ! kernel function for integral 1

  function cmplx_fun1(g) result(c1)

  implicit none
  real(dp), intent(in) :: g
  complex(dp) :: c1

  c1%re=g*dexp(-g*g)
  c1%im=0._dp
	
  end function cmplx_fun1

  ! kernel function for integral 2

  function cmplx_fun2(g) result(c2)

  implicit none
  real(dp), intent(in) :: g
  real(dp) :: g2
  complex(dp) :: c2

  g2=g*g
  c2%re=g2*dexp(-g2)
  c2%im=0._dp

  end function cmplx_fun2

  ! kernel function for integral 3

  function cmplx_fun3(g) result(c3)

  implicit none
  real(dp), intent(in) :: g
  complex(dp) :: c3
    
  c3%re=dexp(-g)
  c3%im=0._dp
    
  end function cmplx_fun3

  ! kernel function for integral 4

  function cmplx_fun4(g) result(c4)

  implicit none
  real(dp), intent(in) :: g
  complex(dp) :: c4

  c4%re=g*dexp(-2._dp*g)
  c4%im=0._dp
    
  end function cmplx_fun4

  ! kernel function for integral 5

  function cmplx_fun5(g) result(c5)

  implicit none
  real(dp), intent(in) :: g
  complex(dp) :: c5

  c5%re=dexp(-2._dp*g)
  c5%im=0._dp
    
  end function cmplx_fun5  
  
end program driver