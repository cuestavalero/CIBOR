!
!----------------------------------------------------------------------
! Module containing all the tools to perform a linear regression
! analysis.
!
! This module needs:
! - 01b_mod_kinds.f03
! - 01c_mod_logger.f03
!
! Francisco Jose Cuesta Valero
! 2021-01-02 (Leipzig)
!----------------------------------------------------------------------

module regression_module
  use kinds_module
  use logger_module
  implicit none
  private

  type, public :: regression_type
    real(rprec), allocatable :: table(:,:) ! Input data: Y, X1, X2,..., Xn
    real(rprec), allocatable :: coefficients(:) ! 1 = intercept ; 2 = slope
    real(rprec), allocatable :: sgm_coefficients(:) ! 2sigma values
    real(rprec) :: r2 ! r2=1 means perfect fit
    real(rprec) :: pvalue ! pvalue .lt. 0.05 means significant fit
    integer(iprec) :: n ! Number of records used to get the fit
  contains
    procedure :: new_1D_abscissa, new_2D_abscissa
    generic :: new => new_1D_abscissa, new_2D_abscissa
    procedure :: delete
    procedure :: linreg
  end type

  ! For logging purposes
  type(logger_type) :: book
  integer(iprec), parameter :: log_lvl = warning

contains
  ! ------------------------
  ! Logging.
  ! ------------------------
  subroutine init_logger()
    ! This routine is called from the constructor to initialize the logger
    logical, save :: loggerInitialized = .false.

    if ( .not. loggerInitialized ) then
      call book%new(log_lvl)
      loggerInitialized = .true.
    end if
  end subroutine

  ! ------------------------
  ! Standard ADT Methods. Construction, Destruction, Copying, and Assignment.
  ! ------------------------

  subroutine new_1D_abscissa(self, y, x)
    ! Subroutine to create a new regression object. This is done by saving
    !            the input data to be fitted.
    ! - self :: the object to be created.
    ! - y :: ordinate coordinate (real 1D array)
    ! - x :: abscissa coordinates (real 1D array)
    class(regression_type), intent(out) :: self
    real(rprec), intent(in) :: y(:), x(:)

    allocate(self%table( size(y), 2 ))
    allocate(self%coefficients( size(x)+1 ))
    allocate(self%sgm_coefficients( size(x)+1 ))

    self%table(:,1) = y
    self%table(:,2) = x
    self%n = size(y)

    return
  end subroutine

  subroutine new_2D_abscissa(self, y, x)
    ! Subroutine to create a new regression object. This is done by saving
    !            the input data to be fitted.
    ! - self :: the object to be created.
    ! - y :: ordinate coordinate (real 1D array)
    ! - x :: abscissa coordinates (real 2D array)
    class(regression_type), intent(out) :: self
    real(rprec), intent(in) :: y(:)
    real(rprec), intent(in) :: x(:,:)

    allocate(self%table( size(y), (size(x,2)+1) ))
    allocate(self%coefficients( size(x,2)+1 ))
    allocate(self%sgm_coefficients( size(x,2)+1 ))

    self%table(:,1) = y
    self%table(:,2:(size(x,2)+1)) = x(:,:)
    self%n = size(y)

    return
  end subroutine

  subroutine delete(self)
    ! Subroutine to delete a regression object.
    ! - self :: object to be deleted.
    class(regression_type), intent(inout) :: self

    deallocate(self%table)
    deallocate(self%coefficients)
    deallocate(self%sgm_coefficients)

    return
  end subroutine


  ! ------------------------
  ! Other methods.
  ! ------------------------

  subroutine linreg(self)
    ! See "Statistical Analysis in Climate Research" (von Storch, Zwiers)
    ! Subroutine to perform linear regression analyses.
    ! - self = object to be fitted (record). The structure here is
    !                 (y,x1,x2,...,xn)
    class(regression_type), intent(inout) :: self

    integer(iprec) :: n, nt
    real(rprec), dimension(:), allocatable :: a ! Coefficients
    real(rprec), dimension(:), allocatable :: ca ! 2sigma values for each a
    real(rprec) :: r2, pvalue
    real(rprec), dimension(:,:), allocatable :: yy
    real(rprec), dimension(:,:), allocatable :: xi, a4
    real(rprec), dimension(:,:), allocatable :: a1, a2
    real(rprec), dimension(:,:), allocatable :: a3
    real(rprec), dimension(:,:), allocatable :: a5, u, yi
    real(rprec), dimension(:,:), allocatable :: a6
    real(rprec), dimension(:,:), allocatable :: e
    real(rprec), dimension(:,:), allocatable :: sst, ssr, sse, fp
    real(rprec) :: dfe, dfr
    real(rprec) :: f, sigma, alpha, p, t
    integer(iprec) :: i,j

    character(len=*), parameter :: pname = "mod_regression < linreg"

    call init_logger

    call book%record(pname,trace,'Starting regression analysis')

    n = size(self%table,2) - 1_iprec
    nt = size(self%table,1)

    ! Check n and nt
    call book%record(pname,debug,"Number variables - ",n)
    call book%record(pname,debug,"Number records - ",nt)

    allocate(a(n+1))
    allocate(ca(n+1))
    allocate(yy(nt,1))
    allocate(xi(nt,n+1))
    allocate(a4(nt,n+1))
    allocate(a1(n+1,n+1))
    allocate(a2(n+1,n+1))
    allocate(a3(n+1,nt))
    allocate(a5(nt,nt))
    allocate(u(nt,nt))
    allocate(yi(nt,nt))
    allocate(a6(1,nt))
    allocate(e(n+1,1))
    allocate(sst(1,1))
    allocate(ssr(1,1))
    allocate(sse(1,1))
    allocate(fp(1,1))

    if(nt.lt.3_iprec) then
      call book%record(pname,error,&
                           'ERROR - less than 3 records')
      stop
    end if

    dfr=dble(n)
    dfe=dble(nt-n-1)

    do i=1,nt
      xi(i,1) = 1.0_rprec
      do j=2,n+1
        xi(i,j) = self%table(i,j)
      end do
    end do
    do i=1,nt
      yy(i,1) = self%table(i,1)
    end do

    do i=1,nt
      do j=1,nt
        u(i,j)=1.0_rprec/dble(nt)
      end do
    end do
    yi=identity(nt)

    a1=matmul(transpose(xi),xi)
    a2=invert(n+1,a1)
    a3=matmul(a2,transpose(xi))

    a=matmul(a3,self%table(:,1))

    a4=matmul(xi,a2)
    a5=matmul(a4,transpose(xi))
    a6=matmul(transpose(yy),(a5-u))
    ssr=matmul(a6,yy)

    sst=matmul(matmul(transpose(yy),yi-u),yy)

    sse=matmul(transpose(yy),matmul(yi-a5,yy))

    r2=ssr(1,1)/sst(1,1)

    f=(ssr(1,1)/dfr)/(sse(1,1)/dfe)
    if(f.lt.0.0_rprec) then
      call book%record(pname,warning,&
                                 "WARNING - Almost Perfect Fit!!!",&
                                 [ssr(1,1),sse(1,1),f])
      f=abs(f)
      sse(1,1)=abs(sse(1,1))
      ssr(1,1)=abs(ssr(1,1))
    end if

    pvalue=p_ffisher(f,dfr,dfe)

    ! Find 95% confidence intervals for a coefficients
    sigma=sqrt(sse(1,1)/dfe)
    alpha=0.95
    ! Find the t value for the interval
    p=(1.0_rprec+alpha)/2.0_rprec
    t=q_tstudent(dfe,p)

    do j=1,n+1 ! sigma values for each coefficient (ca vector)
      e=0.0_rprec
      do i=1,n+1
        if(i.eq.j) e(i,1) = 1.0_rprec
      end do ! i
      fp=matmul(matmul(transpose(e),a2),e)
      ca(j)= t * sigma * sqrt(fp(1,1))
    end do ! j

    self%coefficients = a
    self%sgm_coefficients = ca
    self%r2 = r2
    self%pvalue = pvalue

    call book%record(pname,trace,'Endding regression analysis')

    return
  end subroutine linreg

  function p_ffisher(f,v1,v2) result(p_f)
    ! F-Fisher Probability Function
    ! See "Numerical Recipes" (Press et al.)
    ! Function to obtain the probability (quantile?) of a f-Fisher
    ! distribution given an f value and two degrees of freedom
    implicit none

    real(rprec), intent(in) :: f, v1, v2
    real(rprec) :: p_f

    real(rprec) :: x

    x=v2/(v2+v1*f)

    p_f=betai(v2/2.0_rprec,v1/2.0_rprec,x)
    return
  end function p_ffisher

  function p_tstudent(df,t) result(p_t)
    ! NOTE that Numerical Recipes is wrong by a 0.5 factor!
    ! Requires kinds module and betai function (nrmodule)
    ! Function for obtaining the quantile? of a t-Student distribution
    ! given a t value and the digrees of freedom.
    ! See "Numerical Recipes" (Press et al.)
    implicit none

    real(rprec), intent(in) :: df,t ! Degrees of freedom and t value
    real(rprec) :: p_t

    real(rprec) :: a,b,x ! Dummy constants

    a=df/2.0_rprec
    b=0.5_rprec
    x=df/(df+t**2)
    p_t=1.0_rprec-0.5_rprec*betai(a,b,x)

    if(t.lt.0.0_rprec) p_t=1.0_rprec-p_t ! Not sure of this feature

    return
  end function p_tstudent

  function q_tstudent(df,p) result(q_t)
    ! Funtion for finding the t value given a probability and a degree of
    ! freedom (t-Student distribution)
    implicit none

    real(rprec), intent(in) :: df ! Degrees of freedom
    real(rprec), intent(in) :: p ! Probability
    real(rprec) :: q_t

    real(rprec), parameter :: tol=0.00001_rprec ! Tolerance to find the t value
    real(rprec), parameter :: t0=0.0_rprec ! Initial value for searching
    real(rprec), parameter :: tf=10000.0_rprec ! Final value for searching
    real(rprec), parameter :: dz=0.00001_rprec ! Step for searching

    real(rprec) :: tsol, t ! Dummy constants
    integer(iprec) :: nn ! Dummy constants
    integer(iprec) :: i ! Indices

    tsol=-1.0_rprec
    nn=int((tf-t0)/dz)
    do i=0,nn
      t=t0+dble(i*dz)
      if(p-tol.le.p_tstudent(df,t).and.p+tol.ge.p_tstudent(df,t)) then
        tsol=t
        exit
      end if
    end do ! i
    if(tsol.eq.-1.0_rprec) stop 'Could not find an appropriate t value in q_tstudent!'
    q_t=tsol
    return
  end function q_tstudent

  function betai(a,b,x) result(bb)
    ! See Numerical Recipes
    implicit none

    real(rprec) :: bb
    real(rprec), intent(in) :: a,b,x
!   USES betacf,gammln
    real(rprec) :: bt!,betacf,gammln
    if(x.lt.0.0_rprec.or.x.gt.1.0_rprec) stop 'bad argument x in betai'
    if(x.eq.0.0_rprec.or.x.eq.1.0_rprec) then
      bt=0.0_rprec
    else
      bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0_rprec-x))
    endif
    if(x.lt.(a+1.0_rprec)/(a+b+2.0_rprec))then
      bb=bt*betacf(a,b,x)/a
      return
    else
      bb=1.0_rprec-bt*betacf(b,a,1.0_rprec-x)/b
      return
    endif
  end function betai

  function betacf(a,b,x) result(bb)
    ! See Numerical Recipes
    implicit none

    real(rprec) :: bb
    real(rprec), intent(in) :: a,b,x
    integer(iprec), parameter :: MAXIT=100
    real(rprec), parameter :: EPS=3.e-7_rprec,FPMIN=1.e-30_rprec
    integer(iprec) :: m,m2
    real(rprec) :: aa,c,d,del,h,qab,qam,qap
    qab=a+b
    qap=a+1.0_rprec
    qam=a-1.0_rprec
    c=1.0_rprec
    d=1.0_rprec-qab*x/qap
    if(abs(d).lt.FPMIN)d=FPMIN
    d=1.0_rprec/d
    h=d
    do 11 m=1,MAXIT
      m2=2*m
      aa=m*(b-m)*x/((qam+m2)*(a+m2))
      d=1.0_rprec+aa*d
      if(abs(d).lt.FPMIN)d=FPMIN
      c=1.0_rprec+aa/c
      if(abs(c).lt.FPMIN)c=FPMIN
      d=1.0_rprec/d
      h=h*d*c
      aa=-(a+m)*(qab+m)*x/((a+m2)*(qap+m2))
      d=1.0_rprec+aa*d
      if(abs(d).lt.FPMIN)d=FPMIN
      c=1.0_rprec+aa/c
      if(abs(c).lt.FPMIN)c=FPMIN
      d=1.0_rprec/d
      del=d*c
      h=h*del
      if(abs(del-1.0_rprec).lt.EPS)goto 1
11  continue
    stop 'a or b too big, or MAXIT too small in betacf'
1   bb=h
    return
  end function betacf

  function gammln(xx) result(bb)
    ! See Numerical Recipes
    implicit none

    real(rprec) :: bb
    real(rprec), intent(in) :: xx
    integer(iprec) :: j
    real(rprec) :: ser,stp,tmp,x,y,cof(6)
    save :: cof,stp
    data cof,stp/76.18009172947146_rprec,-86.50532032941677_rprec,&
    24.01409824083091_rprec,-1.231739572450155_rprec,.1208650973866179e-2_rprec,&
    -0.5395239384953e-5_rprec,2.5066282746310005_rprec/
    x=xx
    y=x
    tmp=x+5.5_rprec
    tmp=(x+0.5_rprec)*log(tmp)-tmp
    ser=1.000000000190015_rprec
    do 11 j=1,6
      y=y+1.0_rprec
      ser=ser+cof(j)/y
11  continue
    bb=tmp+log(stp*ser/x)
    return
  end function gammln

  function identity(n) result(x)
    ! Function to generate identity matrices
    ! - n = dimension of the matrix (integer)
    ! - x = identity matrix with dimension n
    implicit none

    integer(iprec), intent(in) :: n
    real(rprec), dimension(n,n) :: x

    integer(iprec) :: i, j

    do i=1,n
      do j=1,n
      if(i.eq.j) then
        x(i,j) = 1.0_rprec
      else
        x(i,j) = 0.0_rprec
      end if
      end do
    end do
    return
  end function identity

  function invert(n,x) result(m)
    ! See "Numerical Recipes" (Press et al.)
    ! Function to invert matrices
    ! - n = dimension of the matrix (integer) [legacy]
    ! - x = matrix to be inverted
    ! - m = inverted matrix
    implicit none

    integer(iprec), intent(in) :: n
    real(rprec), dimension(n,n), intent(in) :: x
    real(rprec), dimension(n,n) :: xx
    real(rprec), dimension(n,n) :: m

    real(rprec), dimension(n,n) :: y
    real(rprec), dimension(n) :: b
    integer(iprec), dimension(n) :: indx
    real(rprec) :: d
    integer(iprec) :: i,j

    y = identity(n)

    do i=1,n
      do j=1,n
        xx(i,j)=x(i,j)
      end do
    end do

    call ludcmp(xx,n,n,indx,d)

    do j=1,n
      do i=1,n
        b(i) = y(i,j)
      end do
      call lubksb(xx,n,n,indx,b)
      do  i=1,n
        m(i,j) = b(i)
      end do
    end do
    return
  end function invert

  subroutine ludcmp(a,n,np,indx,d)
    ! See Numerical Recipes
    implicit none

    integer(iprec), intent(in) :: n,np
    integer(iprec), dimension(n), intent(out) :: indx
    real(rprec), intent(out) :: d
    real(rprec), dimension(np,np), intent(inout) :: a
    integer(iprec), parameter :: NMAX=500
    real(rprec), parameter :: TINY=1.0e-20
    integer(iprec) :: i,imax,j,k
    real(rprec) :: aamax,dum,sum
    real(rprec), dimension(NMAX) :: vv(NMAX)
    d=1.0_rprec
    do 12 i=1,n
      aamax=0.0_rprec
      do 11 j=1,n
        if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11    continue
      if (aamax.eq.0.0_rprec) stop 'singular matrix in ludcmp'
      vv(i)=1.0_rprec/aamax
12  continue
    do 19 j=1,n
      do 14 i=1,j-1
        sum=a(i,j)
        do 13 k=1,i-1
          sum=sum-a(i,k)*a(k,j)
13      continue
        a(i,j)=sum
14    continue
      aamax=0.0_rprec
      do 16 i=j,n
        sum=a(i,j)
        do 15 k=1,j-1
          sum=sum-a(i,k)*a(k,j)
15      continue
        a(i,j)=sum
        dum=vv(i)*abs(sum)
        if (dum.ge.aamax) then
          imax=i
          aamax=dum
        endif
16    continue
      if (j.ne.imax)then
        do 17 k=1,n
          dum=a(imax,k)
          a(imax,k)=a(j,k)
          a(j,k)=dum
17      continue
        d=-d
        vv(imax)=vv(j)
      endif
      indx(j)=imax
      if(a(j,j).eq.0.)a(j,j)=TINY
      if(j.ne.n)then
        dum=1.0_rprec/a(j,j)
        do 18 i=j+1,n
          a(i,j)=a(i,j)*dum
18      continue
      endif
19  continue
    return
  end subroutine ludcmp

  subroutine lubksb(a,n,np,indx,b)
    ! See Numerical Recipes
    implicit none

    integer(iprec), intent(in) :: n,np
    integer(iprec), dimension(n), intent(in) :: indx
    real(rprec), dimension(np,np) :: a(np,np)
    real(rprec), dimension(n), intent(inout) :: b(n)
    integer(iprec) :: i,ii,j,ll
    real(rprec) :: sum
    ii=0
    do 12 i=1,n
      ll=indx(i)
      sum=b(ll)
      b(ll)=b(i)
      if (ii.ne.0)then
        do 11 j=ii,i-1
          sum=sum-a(i,j)*b(j)
11      continue
      else if (sum.ne.0.0_rprec) then
        ii=i
      endif
      b(i)=sum
12  continue
    do 14 i=n,1,-1
      sum=b(i)
      do 13 j=i+1,n
        sum=sum-a(i,j)*b(j)
13    continue
      b(i)=sum/a(i,i)
14  continue
    return
  end subroutine lubksb

end module regression_module

