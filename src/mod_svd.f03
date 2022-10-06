!
!----------------------------------------------------------------------
! Module to propagate surface signals through the ground using a purely
!  conductive forward model.
!
! This module needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_profile.f03 (it requires mod_regression.f03)
! - mod_fm (it requires mod_profile.f03)
!
! Francisco Jose Cuesta-Valero
! 2021-08-19 (Leipzig)
!----------------------------------------------------------------------
module svd_module
  use kinds_module
  use logger_module
  use profile_module
  use fm_module
  implicit none

  private

  type, public :: svd_type
    ! Thermal diffusivity [m2 s-1]
    real(rprec) :: alpha
    ! Number of eigenvalues to be retained [1]
    integer(iprec) :: nw
    ! Logging year [C.E.]
    real(rprec) :: logy
    ! Model for surface signal [time]
    real(rprec), allocatable :: time_series(:)
    ! Profile to be inverted [profile]
    type(profile_type) :: prof
    ! Original inversions [K]
    real(rprec), allocatable :: oinv(:,:)
    ! Inversions [K or W m-2]
    real(rprec), allocatable :: inv(:,:)
    ! Forward models [K or W m-2]
    real(rprec), allocatable :: fm(:,:)
  contains
    procedure :: new_svd, new_copy_svd
    generic :: new => new_svd, new_copy_svd
    procedure :: delete
    procedure :: inversion => inversion_svd
    procedure :: put => write_svd
  end type svd_type

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
  subroutine new_svd(self,alpha,nw,logy,time_series,prof)
    ! Subrotine to create SVD objects
    ! - self :: object to be created (svd)
    ! - alpha :: thermal diffusivity (rprec)
    ! - nw :: number of eigenvalues (iprec)
    ! - time_series :: timesteps for the surface signal (rprec(:))
    ! - prof :: profile to be inverted (profile)
    class(svd_type), intent(inout) :: self
    real(rprec), intent(in) :: alpha
    integer(iprec), intent(in) :: nw
    real(rprec), intent(in) :: logy
    real(rprec), allocatable, intent(in) :: time_series(:)
    type(profile_type), intent(in) :: prof

    self%alpha = alpha
    self%nw = nw
    self%logy = logy
    self%time_series = time_series
    self%prof = prof

    return
  end subroutine new_svd

  subroutine new_copy_svd(self,other)
    ! Subrotine to copy SVD objects
    ! - self :: object to be created (svd)
    ! - other :: object to be copied (svd)
    class(svd_type), intent(out) :: self
    class(svd_type), intent(in) :: other

    self%alpha = other%alpha
    self%nw = other%nw
    self%logy = other%logy
    self%time_series = other%time_series
    self%prof = other%prof
    if(allocated(other%inv)) self%inv = other%inv
    if(allocated(other%oinv)) self%oinv = other%oinv
    if(allocated(other%fm)) self%fm = other%fm

    return
  end subroutine new_copy_svd

  subroutine delete(self)
    ! Subroutine to delete a SVD object
    ! - self :: object to be deleted.
    class(svd_type) :: self

    if(allocated(self%time_series)) then
      deallocate(self%time_series)
    end if
    if(allocated(self%inv)) then
      deallocate(self%inv)
    end if
    if(allocated(self%oinv)) then
      deallocate(self%oinv)
    end if
    if(allocated(self%fm)) then
      deallocate(self%fm)
    end if
    if(allocated(self%prof%depth)) then
      call self%prof%delete()
    end if

    return
  end subroutine delete

  ! ------------------------
  ! Principal subroutines
  ! ------------------------
  subroutine inversion_svd(self)
    ! Subroutine to invert a profile
    ! - self :: svd object (svd)
    class(svd_type), intent(inout) :: self

    ! Forward model of the inversions
    type(fm_type) :: forward

    ! Filling value
    real(rprec), parameter :: fils = -999.9_rprec

    ! Dummy array for saving profile and ensure a good number of depths
    real(rprec), allocatable :: profi(:,:), y(:,:)
    real(rprec), allocatable :: inv(:)

    ! Dummy integers
    integer(iprec) :: ndat, nu

    ! Dummy indices
    integer(iprec) :: i

    character(len=*), parameter :: pname = "mod_svd < inversion_svd"

    call init_logger

    ndat = size(self%prof%depth)
    nu = size(self%time_series)

    ! Dimension of array to save final inversion
    call book%record(pname,trace,'Allocate solution arrays')
    allocate(self%inv(int(maxval(self%time_series),iprec),4))
    allocate(self%oinv(nu,3))
    allocate(self%fm(ndat,4))
    self%inv=fils
    self%fm=fils

    ! The algorithm needs a sufficient number of temp. records
    if(ndat.lt.nu) then
      call book%record(pname,warning,'WARNING - ndat is smaller than nu!!')

      ! Add null records to complete the anomaly profile
      allocate(profi(nu+1,2))
      profi(:,2) = 0.0_rprec
      do i = 1, nu+1
        if(i.le.ndat) then
          profi(i,1) = self%prof%depth(i)
        else
          profi(i,1) = profi(i-1,1) + &
                             ( self%prof%depth(ndat) - &
                               self%prof%depth(ndat-1) )
        end if
      end do
    else
      allocate(profi(ndat,2))
      profi(:,1) = self%prof%depth
      profi(:,2) = 0.0_rprec
    end if


    ! I have to call svd subroutine for each allocated anomaly
    if(allocated(self%prof%atemp_minus)) then
      call book%record(pname,trace,'Invert minus profile')
      do i = 1, ndat
        profi(i,2) = self%prof%atemp_minus(i)
      end do

      if(allocated(inv)) deallocate(inv)
      call svd(self,profi,inv)
      self%oinv(:,1) = inv
      call inversion2year(self%logy,inv,self%time_series,y)
      self%inv(:,1) = y(:,1)
      self%inv(:,2) = y(:,2)
      deallocate(y)

      call book%record(pname,trace,'Forward model minus inversion')
      call forward%new(self%alpha,self%time_series,inv,self%prof%depth)
      call forward%forward_model()
      self%fm(:,1) = forward%prof%depth
      self%fm(:,2) = forward%prof%atemp_original
      call forward%delete()
    else
      allocate(self%prof%atemp_minus(size(self%prof%depth)))
      self%prof%atemp_minus = fils
    end if

    if(allocated(self%prof%atemp_original)) then
      call book%record(pname,trace,'Invert original profile')
      do i = 1, ndat
        profi(i,2) = self%prof%atemp_original(i)
      end do

      if(allocated(inv)) deallocate(inv)
      call svd(self,profi,inv)
      self%oinv(:,2) = inv
      call inversion2year(self%logy,inv,self%time_series,y)
      self%inv(:,1) = y(:,1)
      self%inv(:,3) = y(:,2)
      deallocate(y)

      call book%record(pname,trace,'Forward model original inversion')
      call forward%new(self%alpha,self%time_series,inv,self%prof%depth)
      call forward%forward_model()
      self%fm(:,1) = forward%prof%depth
      self%fm(:,3) = forward%prof%atemp_original
      call forward%delete()
    else
      call book%record(pname,error,"ERROR - There are no anomalies to invert!")
      stop
    end if

    if(allocated(self%prof%atemp_plus)) then
      call book%record(pname,trace,'Invert plus profile')
      do i = 1, ndat
        profi(i,2) = self%prof%atemp_plus(i)
      end do

      if(allocated(inv)) deallocate(inv)
      call svd(self,profi,inv)
      self%oinv(:,3) = inv
      call inversion2year(self%logy,inv,self%time_series,y)
      self%inv(:,1) = y(:,1)
      self%inv(:,4) = y(:,2)
      deallocate(y)

      call book%record(pname,trace,'Forward model plus inversion')
      call forward%new(self%alpha,self%time_series,inv,self%prof%depth)
      call forward%forward_model()
      self%fm(:,1) = forward%prof%depth
      self%fm(:,4) = forward%prof%atemp_original
      call forward%delete()
    else
      allocate(self%prof%atemp_plus(size(self%prof%depth)))
      self%prof%atemp_plus = fils
    end if

    return
  end subroutine inversion_svd

  subroutine svd(self,prof,inv)
    ! Subroutine to apply the svd method.
    ! - self :: inversion object with parameters (svd).
    ! - prof :: anomly profile to be inverted (real(:,:))
    ! - inv :: array with inversion (real(:))
    class(svd_type), intent(in) :: self
    real(rprec), intent(in) :: prof(:,:)
    real(rprec), allocatable, intent(out) :: inv(:)

    real(rprec), allocatable :: tinv(:) ! time inverse
    real(rprec), allocatable :: w(:) ! eigenvalues [1]
    real(rprec), allocatable :: v(:,:) ! eigenvectors [1]
    real(rprec), allocatable :: x(:) ! solution of the inversion [K]

    real(rprec), allocatable :: a(:,:) ! Matrix for the system of equations

    integer(iprec) :: nu, ndat ! Dummy integers
    integer(iprec) :: i, j ! Dummy indices

    character(len=*), parameter :: pname = "mod_svd < svd"

    call init_logger

    ndat = size(prof,1)
    nu = size(self%time_series)

        ! Prepare matrix for system of equations
    allocate(a(ndat,nu))

    ! Prepare matrix for eigenvectors
    allocate(v(nu,nu))

    ! Prepare array for eigenvalues
    allocate(w(nu))

    ! Prepare the solution of the inversion
    allocate(x(nu))

    ! Prepare time inverse
    allocate(tinv(nu))

    ! Create the inversion model
    do i = 1, nu
      tinv(i) = dsqrt(self%alpha * self%time_series(i)) ! Required for the equations' matrix
    end do ! i

    ! Create the matrix for the system of equations
    a = 0.0_rprec
    do i = 1, ndat
      do j = 1, nu
        if(prof(i,1).eq.0.0_rprec) then
          a(i,j) = 1.0_rprec
          if(j.ne.1) then
            a(i,j) = 0.0_rprec
          end if
        else
          a(i,j) = erfcc( 0.5_rprec * prof(i,1) / tinv(j) )
          if(j.ne.1) then
            a(i,j) =  a(i,j) - erfcc( 0.5_rprec * prof(i,1) / &
                      tinv(j-1) )
          end if
        end if
      end do ! j
    end do ! i

    ! Call the singular value decomposition subroutine
    call svdcmp(a,ndat,nu,ndat,nu,w,v)

    ! Range of eigenvalues to be taken into account
    ! Select adequate number of eigenvalues
    do i = self%nw+1, nu
      w(i) = 0.0_rprec
    end do

    ! Perform the reconstruction
    call svbksb(a,w,v,ndat,nu,ndat,nu,prof(:,2),x)

    ! save the results
    allocate(inv(nu))
    inv = x

    return
  end subroutine svd

  subroutine write_svd(self,ofile,u2)
    ! Subroutine to write the results of svd objects.
    ! - self :: svd object to be written (svd)
    ! - ofile :: name of the file to write the results (char)
    ! - u2 :: unit for writing (integer, optional)
    class(svd_type), intent(in) :: self
    character(len=*), intent(in) :: ofile
    integer(iprec), optional, intent(in) :: u2

    real(rprec), parameter :: fils = -999.9_rprec
    integer(iprec), parameter :: ifils = -999_iprec

    character(90) :: f ! Format for writing

    integer(iprec) :: u ! Dummy integer
    integer(iprec) :: i ! Dummy indices

    character(len=*), parameter :: pname = "mod_svd < write_svd"

    call init_logger


    if(present(u2)) then
      u = u2
    else
      u = 90_iprec
    end if

    open(u,file=trim(ofile),action="write")

    ! First, relevant parameters of the model
    f = "(a,es12.3e3)"
    if(test_real(self%prof%zmin,fils).eq.1) then
      write(u,trim(f)) '#z1', self%prof%zmin
    end if
    if(test_real(self%prof%zmax,fils).eq.1) then
      write(u,trim(f)) '#z2', self%prof%zmax
    end if

    write(u,*) ""
    write(u,*) ""

    ! The equilibrium profile parameters
    f = "(a,es12.3e3,es12.3e3)"
    if(test_real(self%prof%eq_profile%inter,fils).eq.1.and.&
       test_real(self%prof%eq_profile%s_inter,fils).eq.1) then
      write(u,trim(f)) '#t0', self%prof%eq_profile%inter&
                                  ,self%prof%eq_profile%s_inter
    end if
    if(test_real(self%prof%eq_profile%slope,fils).eq.1.and.&
       test_real(self%prof%eq_profile%s_slope,fils).eq.1) then
      write(u,trim(f)) '#gamma', self%prof%eq_profile%slope,&
                                  self%prof%eq_profile%s_slope
    end if

    f = "(a,i6)"
    if(test_integer(self%prof%eq_profile%n,ifils).eq.1) then
      write(u,trim(f)) '#n', self%prof%eq_profile%n
    end if

    f = "(a,es12.3e3)"
    if(test_real(self%prof%eq_profile%r2,fils).eq.1) then
      write(u,trim(f)) '#r2', self%prof%eq_profile%r2
    end if
    if(test_real(self%prof%eq_profile%pvalue,fils).eq.1) then
      write(u,trim(f)) '#pvalue', self%prof%eq_profile%pvalue
    end if

    write(u,*) ""
    write(u,*) ""

    ! Write the original profile
    f = "(a,es12.3e3,es12.3e3)"
    if(test_array(self%prof%temp,fils).eq.1) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#o', self%prof%depth(i), self%prof%temp(i)
      end do

      write(u,*) ""
      write(u,*) ""

    else
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#o', self%prof%depth(i), fils
      end do

      write(u,*) ""
      write(u,*) ""
    end if

    ! Write the anomaly profiles
    f = "(a,es12.3e3,es12.3e3,es12.3e3,es12.3e3)"
    if(test_array(self%prof%atemp_minus,fils).eq.1.and.&
       test_array(self%prof%atemp_plus,fils).eq.1) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              self%prof%atemp_minus(i),&
                              self%prof%atemp_original(i),&
                              self%prof%atemp_plus(i)
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(test_array(self%prof%atemp_minus,fils).eq.1.and.&
       test_array(self%prof%atemp_plus,fils).eq.0) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              self%prof%atemp_minus(i),&
                              self%prof%atemp_original(i),&
                              fils
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(test_array(self%prof%atemp_minus,fils).eq.0.and.&
       test_array(self%prof%atemp_plus,fils).eq.1) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              fils,&
                              self%prof%atemp_original(i),&
                              self%prof%atemp_plus(i)
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(test_array(self%prof%atemp_minus,fils).eq.0.and.&
       test_array(self%prof%atemp_plus,fils).eq.0) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              fils,&
                              self%prof%atemp_original(i),&
                              fils
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(test_array(self%prof%atemp_original,fils).eq.0) then
      call book%record(pname,error,"ERROR - There is no anomaly to write")
      stop
    end if


    ! Write general solution
    f = "(a,4es12.3e3)"
    if(test_array(self%inv(:,2)).eq.0.or.test_array(self%inv(:,4)).eq.0) then
      do i = 1, size(self%inv,1)
        write(u,trim(f)) '#g', self%inv(i,1), fils,&
                       self%inv(i,3), fils
      end do

      write(u,*)
      write(u,*)
    elseif(test_array(self%inv(:,3)).eq.0) then
      call book%record(pname,error,"ERROR - There is no inversion to write")
      stop
    else
      do i = 1, size(self%inv,1)
        write(u,trim(f)) '#g', self%inv(i,1), min(self%inv(i,2),self%inv(i,4)),&
                       self%inv(i,3), max(self%inv(i,2),self%inv(i,4))
      end do

      write(u,*)
      write(u,*)
    end if


    ! Write general forward model
    f = "(a,4es12.3e3)"
    if(test_array(self%fm(:,2)).eq.0.or.test_array(self%fm(:,4)).eq.0) then
      do i = 1, size(self%fm,1)
        write(u,trim(f)) '#f', self%fm(i,1), fils,&
                       self%fm(i,3), fils
      end do
    elseif(test_array(self%fm(:,3)).eq.0) then
      call book%record(pname,error,"ERROR - There is no forward model do write")
      stop
    else
      do i = 1, size(self%fm,1)
        write(u,trim(f)) '#f', self%fm(i,1), min(self%fm(i,2),self%fm(i,4)),&
                       self%fm(i,3), max(self%fm(i,2),self%fm(i,4))
      end do
    end if

    return
  end subroutine write_svd

  ! ------------------------
  ! Auxiliar subroutines
  ! ------------------------
  function erfcc(x) result(r)
    ! IMPORTANT: this function is vital for the forward model subroutine,
    !            so I prefer to keep both the subroutine and the function
    !            in the same module
    real(rprec), intent(in) :: x

    real(rprec) :: z,t, r

    z  =  dabs(x)
    t  =  1.0_rprec/(1.0_rprec+0.5_rprec*z)
    r  =  t*dexp( - z*z - 1.26551223_rprec &
      + t*(1.00002368_rprec + t*(0.37409196_rprec &
      + t*(0.09678418_rprec + t*( - 0.18628806_rprec &
      + t*(0.27886807_rprec + t*( - 1.13520398_rprec &
      + t*(1.48851587_rprec + t*( - 0.82215223_rprec &
      + t*0.17087277_rprec)))))))))
    if (x.lt.0.0_rprec) r  =  2.0_rprec - r

    return
  end function erfcc

  subroutine inversion2year(logy,x,xtime,y)
    real(rprec), intent(in) :: logy
    real(rprec), intent(in) :: x(:), xtime(:)

    real(rprec), allocatable, intent(out) :: y(:,:)

    integer(iprec) :: n, nt
    integer(iprec) :: i, j

    character(len=*), parameter :: pname="mod_svd < inversion2year"

    n = size(xtime)
    nt = int(maxval(xtime),iprec)

    allocate(y(nt,2))

    ! Create temporal series
    do i = 1, nt
      y(i,1) = logy - real( (i-1_iprec),rprec)
    end do

    ! Trasform the inversion into annual data
    j = 1_iprec
    do i = 1, nt
      if(y(i,1).ge.(logy-xtime(j))) then
        y(i,2) = x(j)
      elseif(y(i,1).lt.(logy-xtime(j)).and.j+1.lt.n) then
        j = j + 1_iprec
        y(i,2) = x(j)
      elseif(j+1.eq.n.and.y(i,1).le.(logy-xtime(j))) then
        y(i,2) = x(j+1)
      elseif(j+1.gt.n) then
        exit
      end if
    end do

    return
  end subroutine inversion2year

  function test_array(x,fils) result(t)
    ! Function to check that 1D arrays does not include strange results.
    ! - inversion => the array to be checked (record)
    ! - fils :: filling value (real, optional)
    ! - t => 1 if all test are successful, 0 otherwise
    real(rprec), intent(in) :: x(:)
    real(rprec), optional, intent(in) :: fils
    integer(iprec) :: t

    character(len=7) :: c0,c1,c2,c3 ! Dummy characters for logger

    integer(iprec) :: i, ind ! Dummy indice

    character(len=*), parameter :: pname = "mod_svd < test_array"

    call init_logger

    call book%record(pname,trace,"Starting function")

    t = 1_iprec ! Initialize test

    ! Test if there is any filling value
    if(present(fils)) then
      ind = 0_iprec
      do i = 1, size(x)
        if(x(i).eq.fils) ind = ind + 1_iprec
      end do ! i
      if(ind.gt.0_iprec) then
        t = 0_iprec
        c0 = "z=0"
      else
        c0 = "z=1"
      end if
      call book%record(pname,debug,"Fils test -> "//trim(c0), [ind,size(x)])
    end if

    ! Test if all values are zero
    ind = 0_iprec
    do i = 1, size(x)
      if(x(i).eq.0.0_rprec) ind = ind + 1_iprec
    end do ! i
    if(ind.eq.size(x)) then
      t = 0_iprec
      c1 = "z=0"
    else
      c1 = "z=1"
    end if
    call book%record(pname,debug,"Zero test -> "//trim(c1), [ind,size(x)])

    ! Check if there are NANs in the inversion
    ind = 0_iprec
    do i = 1, size(x)
      if(x(i).ne.x(i)) then
        ind = ind + 1_iprec
        call book%record(pname,debug,"Index of NA - ", i)
      end if
    end do ! i
    if(ind.ne.0_iprec) then
      t = 0_iprec
      c2 = "NA=0"
    else
      c2 = "NA=1"
    end if
    call book%record(pname,debug,"NANs test -> "//trim(c2), ind)

    ! Check there are no infinite values
    ind = 0_iprec
    do i = 1, size(x)
      if(x(i).ge.huge(rprec)) then
        ind = ind + 1_iprec
        call book%record(pname,debug,"Index of Inf - ", i)
      end if
    end do !i
    if(ind.ne.0_iprec) then
      t = 0_iprec
      c3 = "infty=0"
    else
      c3 = "infty=1"
    end if
    call book%record(pname,debug,"Infty test -> "//trim(c3), ind)

    call book%record(pname,trace,"Ending function")

    return
  end function test_array

  function test_real(x,fils) result(t)
    ! Function to check that reals do not include strange results.
    ! - x => the real to be checked (real)
    ! - fils :: filling value (real, optional)
    ! - t => 1 if all test are successful, 0 otherwise
    real(rprec), intent(in) :: x
    real(rprec), optional, intent(in) :: fils
    integer(iprec) :: t

    character(len=7) :: c1,c2,c3 ! Dummy characters for logger

    character(len=*), parameter :: pname = "mod_svd < test_real"

    call init_logger

    call book%record(pname,trace,"Starting function")

    t = 1_iprec ! Initialize test

    ! Test if value is fils
    if(present(fils)) then
      if(x.eq.fils) then
        t = 0_iprec
        c1 = "z=0"
      else
        c1 = "z=1"
      end if
      call book%record(pname,debug,"Fils test -> "//trim(c1))
    end if

    ! Check if there are NANs in the inversion
    if(x.ne.x) then
      t = 0_iprec
      c2 = "NA=0"
    else
      c2 = "NA=1"
    end if
    call book%record(pname,debug,"NANs test -> "//trim(c2))

    ! Check there are no infinite values
    if(x.ge.huge(rprec)) then
      t = 0_iprec
      c3 = "infty=0"
    else
      c3 = "infty=1"
    end if
    call book%record(pname,debug,"Infty test -> "//trim(c3))

    call book%record(pname,trace,"Ending function")

    return
  end function test_real

  function test_integer(x,fils) result(t)
    ! Function to check that integers do not include strange results.
    ! - x => the real to be checked (real)
    ! - fils :: filling value (integer, optional)
    ! - t => 1 if all test are successful, 0 otherwise
    integer(iprec), intent(in) :: x
    integer(iprec), optional, intent(in) :: fils
    integer(iprec) :: t

    character(len=7) :: c1,c2,c3 ! Dummy characters for logger

    character(len=*), parameter :: pname = "mod_svd < test_integer"

    call init_logger

    call book%record(pname,trace,"Starting function")

    t = 1_iprec ! Initialize test

    ! Test if value is fils
    if(present(fils)) then
      if(x.eq.fils) then
        t = 0_iprec
        c1 = "z=0"
      else
        c1 = "z=1"
      end if
      call book%record(pname,debug,"Fils test -> "//trim(c1))
    end if

    ! Check if there are NANs in the inversion
    if(x.ne.x) then
      t = 0_iprec
      c2 = "NA=0"
    else
      c2 = "NA=1"
    end if
    call book%record(pname,debug,"NANs test -> "//trim(c2))

    ! Check there are no infinite values
    if(x.ge.huge(iprec)) then
      t = 0_iprec
      c3 = "infty=0"
    else
      c3 = "infty=1"
    end if
    call book%record(pname,debug,"Infty test -> "//trim(c3))

    call book%record(pname,trace,"Ending function")

    return
  end function test_integer

  subroutine svdcmp(a,m,n,mp,np,w,v)
    ! Modified to be fortran 90 friendly. From Hugo's inversion program
    implicit none

    integer(iprec), intent(in) :: m,n,mp,np
    real(rprec), dimension(mp,np) :: a
    real(rprec), dimension(np) :: w
    real(rprec), dimension(np,np) :: v

    integer(iprec), parameter :: nmax = 200_iprec
    real(rprec), dimension(nmax) :: rv1

    real(rprec) :: g, scales, anorm, f, h, s, y, z, c, x
    integer(iprec) :: i, k, j, l, nm, its

    character(len=*), parameter :: pname = "mod_svd < svdcmp"

    call init_logger

    if(n.gt.nmax) then
      call book%record(pname,error,'nu cannot be larger than nmax = 200')
      call book%record(pname,error,'Yoy cannot have more than 200 time steps')
      stop
    end if

    g = 0.0_rprec
    scales = 0.0_rprec
    anorm = 0.0_rprec
    do i = 1, n
      l = i + 1_iprec
      rv1(i) = scales * g
      g = 0.0_rprec
      s = 0.0_rprec
      scales = 0.0_rprec
      if(i.le.m) then
        do k = i, m
          scales = scales + abs( a(k,i) )
        end do ! k
        if(scales.ne.0.0_rprec) then
          do k = i, m
            a(k,i) = a(k,i) / scales
            s = s + a(k,i) * a(k,i)
          end do
          f = a(i,i)
          g = -dsign( dsqrt(s), f )
          h = f * g - s
          a(i,i) = f - g
          if(i.ne.n) then
            do j = l, n
              s = 0.0_rprec
              do k = i, m
                s = s + a(k,i) * a(k,j)
              end do ! k
              f = s / h
              DO k = i, m
                a(k,j) = a(k,j) + f * a(k,i)
              end do ! k
            end do ! j
          end if
          do k = i, m
            a(k,i) = scales * a(k,i)
          end do ! k
        end if
      end if
      w(i) = scales * g
      g = 0.0_rprec
      s = 0.0_rprec
      scales = 0.0_rprec
      if(i.le.m.and.i.ne.n) then
        do k = l, n
          scales = scales + abs( a(i,k) )
        end do ! k
        if(scales.ne.0.0_rprec) then
          do k = l, n
            a(i,k) = a(i,k) / scales
            s = s + a(i,k) * a(i,k)
          end do ! k
          f = a(i,l)
          g = -dsign( dsqrt(s), f )
          h = f * g - s
          a(i,l) = f - g
          do k = l, n
              rv1(k) = a(i,k) / h
          end do ! k
          if(i.ne.m) then
            do j = l, m
              s = 0.0_rprec
              do k = l, n
                s = s + a(j,k) * a(i,k)
              end do ! k
              do k = l, n
                a(j,k) = a(j,k) + s * rv1(k)
              end do ! k
            end do ! j
          end if
          do k = l, n
            a(i,k) = scales * a(i,k)
          end do ! k
        end if
      end if
      anorm = max( anorm, ( abs( w(i) ) + abs( rv1(i) ) ) )
    end do ! i

    do i = n, 1, -1
      if(i.lt.n) then
        if(g.ne.0.0_rprec) then
          do j = l, n
            v(j,i) = ( a(i,j) / a(i,l) ) / g
          end do ! j
          do j = l, n
            s = 0.0_rprec
            do k = l, n
              s = s + a(i,k) * v(k,j)
            end do ! k
            do k = l, n
              v(k,j) = v(k,j) + s * v(k,i)
            end do ! k
          end do ! j
        end if
        do j = l, n
          v(i,j) = 0.0_rprec
          v(j,i) = 0.0_rprec
        end do ! j
      end if
      v(i,i) = 1.0_rprec
      g = rv1(i)
      l = i
    end do ! i

    do i = n, 1, -1
      l = i + 1_iprec
      g = w(i)
      if(i.lt.n) then
        do j = l, n
          a(i,j) = 0.0_rprec
        end do ! j
      end if
      if(g.ne.0.0_rprec) then
        g = 1.0_rprec / g
        if(i.ne.n) then
          do j = l, n
            s = 0.0_rprec
            do k = l, m
              s = s + a(k,i) * a(k,j)
            end do ! k
            f = ( s / a(i,i) ) * g
            do k = i, m
              a(k,j) = a(k,j) + f * a(k,i)
            end do ! k
          end do ! j
        end if
        do j = i, m
          a(j,i) = a(j,i) * g
        end do ! j
      else
        do j = i, m
          a(j,i) = 0.0_rprec
        end do ! j
      end if
      a(i,i) = a(i,i) + 1.0_rprec
    end do ! i

    do k = n, 1, -1
      do its = 1, 30
        do l = k, 1, -1
          nm = l - 1_iprec
          if( (abs( rv1(l) ) + anorm).eq.anorm) then
            exit
          end if
          if( (abs( w(nm) ) + anorm).eq.anorm) then
            c = 0.0_rprec
            s = 1.0_rprec
            do i = l, k
              f = s * rv1(i)
              if( (abs(f) + anorm).ne.anorm) then
                g = w(i)
                h = dsqrt(f * f + g * g)
                w(i) = h
                h = 1.0_rprec / h
                c = (g * h)
                s = -(f * h)
                do j = 1, m
                  y = a(j,nm)
                  z = a(j,i)
                  a(j,nm) = (y * c) + (z * s)
                  a(j,i) = -(y * s) + (z * c)
                end do ! j
              end if
            end do ! i
            exit
          end if
        end do ! l
        z = w(k)
        if(l.eq.k) then
          if(z.lt.0.0_rprec) then
            w(k) = -z
            do j = 1, n
              v(j,k) = -v(j,k)
            end do ! j
          end if
          exit
        end if
        if(its.eq.30_iprec) stop 'No convergence in 30 iterations'
        x = w(l)
        nm = k - 1_iprec
        y = w(nm)
        g = rv1(nm)
        h = rv1(k)
        f = ( (y - z) * (y + z) + (g - h) * (g + h) ) / (2.0_rprec * h * y)
        g = dsqrt(f * f + 1.0_rprec)
        f = ( (x - z) * (x + z) + h * ( (y / (f + dsign(g, f)) ) - h ) ) / x
        c = 1.0_rprec
        s = 1.0_rprec
        do j = l, nm
          i = j + 1_iprec
          g = rv1(i)
          y = w(i)
          h = s * g
          g = c * g
          z = dsqrt(f * f + h * h)
          rv1(j) = z
          c = f / z
          s = h / z
          f = (x * c) + (g * s)
          g = -(x * s) + (g * c)
          h = y * s
          y = y * c
          do nm = 1, n
            x = v(nm,j)
            z = v(nm,i)
            v(nm,j) = (x * c) + (z * s)
            v(nm,i) = -(x * s) + (z * c)
          end do ! nm
          z = dsqrt(f * f + h * h)
          w(j) = z
          if(z.ne.0.0_rprec) then
            z = 1.0_rprec / z
            c = f * z
            s = h * z
          end if
          f = (c * g) + (s * y)
          x = -(s * g) + (c * y)
          do nm = 1, m
            y = a(nm,j)
            z = a(nm,i)
            a(nm,j) = (y * c) + (z * s)
            a(nm,i) = -(y * s) + (z * c)
          end do ! nm
        end do ! j
        rv1(l) = 0.0_rprec
        rv1(k) = f
        w(k) = x
      end do ! its
    end do ! k

    return
  end subroutine svdcmp

  subroutine svbksb(u,w,v,m,n,mp,np,b,x)
    ! Modified to be fortran 90 friendly. From Hugo's inversion program
    implicit none

    integer(iprec), intent(in) :: mp, np, n, m
    real(rprec), dimension(mp,np) :: u
    real(rprec), dimension(np) :: w, x
    real(rprec), dimension(np,np) :: v
    real(rprec), dimension(mp) :: b

    integer(iprec), parameter :: nmax = 200_iprec
    real(rprec), dimension(nmax) :: tmp

    real(rprec) :: s
    integer(iprec) :: i, j, jj

    character(len=*), parameter :: pname = "mod_svd < svdksb"

    call init_logger

    if(n.gt.nmax) then
      call book%record(pname,error,'nu cannot be larger than nmax = 200')
      call book%record(pname,error,'Yoy cannot have more than 200 time steps')
      stop
    end if

    do j = 1, n
      s = 0.0_rprec
      if(w(j).ne.0.0_rprec) then
        do i = 1, m
          s = s + u(i,j) * b(i)
        end do ! i
        s = s / w(j)
      end if
      tmp(j) = s
    end do ! j

    do j = 1, n
      s = 0.0
      do jj = 1, n
        s = s + v(j,jj) * tmp(jj)
      end do ! jj
      x(j) = s
    end do ! j

    return
  end subroutine svbksb

end module svd_module
