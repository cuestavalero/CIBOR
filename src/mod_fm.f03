!
!----------------------------------------------------------------------
! Module to propagate surface signals through the ground using a purely
!  conductive forward model.
!
! This module needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_profile.f03 (it requires mod_regression.f03)
!
! Francisco Jose Cuesta-Valero
! 2021-08-19 (Leipzig)
!----------------------------------------------------------------------
module fm_module
  use kinds_module
  use logger_module
  use profile_module
  implicit none

  private

  type, public :: fm_type
    ! Thermal diffusivity [m2 s-1]
    real(rprec) :: alpha
    ! Time in time steps before present
    real(rprec), allocatable :: timesteps(:)
    ! Surface signal anomaly [K or W m-2]
    real(rprec), allocatable :: signal(:)
    ! Vector with depths to obtain the prfile
    real(rprec), allocatable :: depths(:)
    ! Final profile
    type(profile_type) :: prof
  contains
    procedure :: new_fm, new_copy_fm
    generic :: new => new_fm, new_copy_fm
    procedure :: delete
    procedure :: forward_model
    procedure :: write_fm
    generic :: put => write_fm
  end type

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
  subroutine new_fm(self,alpha,timesteps,signal,depths)
    ! Subroutine to initialize the FM object.
    ! - self :: object to be initialized (fm).
    ! - alpha :: thermal diffusivity (real)
    ! - timesteps :: time in timesteps before present (real(:))
    ! - signal :: surface signal anomaly (real(:))
    ! - depths :: depths to estimate the profile (real(:))
    class(fm_type), intent(out) :: self
    real(rprec), intent(in) :: alpha
    real(rprec), allocatable, intent(in) :: timesteps(:)
    real(rprec), allocatable, intent(in) :: signal(:)
    real(rprec), allocatable, intent(in) :: depths(:)

    character(len=*), parameter :: pname = "mod_fm < new_fm"

    call init_logger

    self%alpha = alpha
    if(timesteps(2).le.timesteps(1)) then
      call book%record(pname,warning,"WARNING - The order of the timesteps is wrong")
      self%timesteps = timesteps(size(timesteps):1:-1)
      self%signal = signal(size(signal):1:-1)
    else
      self%timesteps = timesteps
      self%signal = signal
    end if
    self%depths = depths

    return
  end subroutine new_fm

  subroutine new_copy_fm(self,other)
    ! Subroutine to create a new FM object from a previously existing one
    ! - self :: object to be created (fm)
    ! - other :: object to be copied (fm)
    class(fm_type), intent(out) :: self
    class(fm_type), intent(in) :: other

    self%alpha = other%alpha
    self%timesteps = other%timesteps
    self%signal = other%signal
    self%depths = other%depths

    return
  end subroutine new_copy_fm

  subroutine delete(self)
    ! Subroutine to delete FM objects
    ! - self :: object to be deleted (fm)
    class(fm_type), intent(inout) :: self

    if(allocated(self%timesteps)) deallocate(self%timesteps)
    if(allocated(self%signal)) deallocate(self%signal)
    if(allocated(self%depths)) deallocate(self%depths)
    if(allocated(self%prof%depth)) call self%prof%delete()

    return
  end subroutine delete

  ! ------------------------
  ! Fundamental calculations
  ! ------------------------
  subroutine forward_model(self)
    ! Subroutine to estimate the forward model of a surface anomaly
    ! - self = forward model object (fm)
    class(fm_type), intent(inout) :: self

    real(rprec),allocatable :: z(:), ff(:) ! Depth and anomaly profile

    integer(iprec) :: nd ! Number of depth points
    integer(iprec) :: i, k, j, nt, maxtim, maxz ! Dummy integers
    real(rprec), dimension(:), allocatable :: u, tm, time ! Dummy real array
    real(rprec) :: x, y, o, t, f, alpha ! Dummy reals

    character(len=*), parameter :: pname = 'mod_fm < forward_model'

    call init_logger

    nt = size(self%timesteps)
    nd = size(self%depths)
    alpha = self%alpha

    ! Allocate dummy arrays
    allocate(u(nt))
    allocate(time(nt))
    allocate(tm(nt))
    allocate(z(nd))
    allocate(ff(nd))

    tm = self%signal
    time = self%timesteps

    if(nt.lt.2_iprec) then
      call book%record(pname,error,'ERROR - temporal series must be larger than 1!')
      stop
    end if

    ! Check that time goes from present to past
    if(time(2).lt.time(1)) then ! check time series. From present to past, in "yr before present" units
      call book%record(pname,warning,"WARNING: Check time series.&
               & From present to past, in 'timesteps before present' units")
      time(:) = time(nt:1:-1)
      tm(:) = tm(nt:1:-1)
    end if

    maxtim = nt

    ! Prepare the profile
    do i=1,nd
      z(i) = self%depths(i)
    end do
    maxz = nd

    ! Forward the anomaly
    do k = 1,maxtim
      u(k) = 0.5_rprec/dsqrt(alpha*time(k))
    enddo
    do j = 1,maxz
      f = 0.0_rprec
      t = tm(1)*erfcc(z(j)*u(1))
      do k=2,maxtim
        x = z(j)*u(k)
        y = z(j)*u(k-1)
        o = tm(k)*(erfcc(x) - erfcc(y))
        f = f + o
      end do
      ff(j) = f + t
    end do

    ! Save the results
    call self%prof%new(z)
    call self%prof%set_anoma(ff)

    return
  end subroutine forward_model

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

  ! ------------------------
  ! Writing functions
  ! ------------------------
  subroutine write_fm(self,ofile,u)
    ! Subroutine to save FM objects
    ! - self :: object to be written (fm)
    ! - ofile :: name of file to save the object (char)
    ! - u :: unit for open command (iprec)
    class(fm_type), intent(in) :: self
    character(len=*), intent(in) :: ofile
    integer(iprec), optional :: u

    integer(iprec) :: i

    character(len=180) :: f ! Format

    if(.not.present(u)) then
      u = 90_iprec
    end if

    open(u,file=trim(ofile),action="write")
    f = "(a,es12.3)"
    write(u,trim(f)) '#alpha', self%alpha

    write(u,*) ""
    write(u,*) ""

    f = "(a,2es12.3)"
    do i = 1,size(self%timesteps)
      write(u,trim(f)) '#sig', self%timesteps(i), self%signal(i)
    end do

    write(u,*) ""
    write(u,*) ""

    f = "(a,2es12.3)"
    do i = 1,size(self%prof%depth)
      write(u,trim(f)) '#fm', self%prof%depth(i), self%prof%atemp_original(i)
    end do
    close(u)

    return
  end subroutine write_fm
end module fm_module










