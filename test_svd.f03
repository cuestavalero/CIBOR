!
!----------------------------------------------------------------------
! Program to test the svd module
!
! This program needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_profile.f03 (it requires mod_regression.f03)
! - mod_svd.f03 (it requires mod_fm.f03)
!
! Francisco Jose Cuesta Valero
! 2021-09-01 (Alhama de Murcia)
!----------------------------------------------------------------------
program test_svd
  use kinds_module
  use logger_module
  use profile_module
  use svd_module
  implicit none

  ! Name of final file
  character(len = *), parameter :: ofile = 'inv_svd.dat'

  ! Prepare logger
  type(logger_type) :: book
  integer(iprec), parameter :: log_lvl = warning
  character(len=*), parameter :: pname="test_svd"

  ! Type of inversion
  type(svd_type) :: svd

  ! Profile type
  type(profile_type) :: prof

  ! Number of eigenvalues retained in the solutions
  integer(iprec), parameter :: nw = 2_iprec

  ! Thermal diffusivity (m2 s-1)
  real(rprec), parameter :: alpha = 1.0e-6_rprec

  ! Number of time steps in the surface history
  integer(iprec), parameter :: nu = 21_iprec

  ! Length of each time step (years)
  real(rprec), parameter :: dt = 30.0_rprec

  ! Minimum and maximum depths to estimate the quasi-equilibrium temperature
  !    profile (m)
  real(rprec), parameter :: az1 = 200.0_rprec
  real(rprec), parameter :: az2 = 300.0_rprec

  ! Logging year of the profile
  real(rprec), parameter :: logy = 2021.0_rprec ! C.E.

  ! Constants
  real(rprec), parameter :: cu = 3600.0_rprec * 24.0_rprec * 365.25_rprec

  real(rprec), allocatable :: signal(:), depths(:), time(:), tx(:), ty(:)

  integer(iprec) :: i


  call book%new(log_lvl)

  call book%record(pname,trace,'Start')

  call book%record(pname,trace,'Create time series')
  allocate(time(nu))
  do i = 1, nu
    time(i) = real(i,rprec) * dt
  end do

  call book%record(pname,trace,'Read temperature profile')
  call read_profile(10_iprec,'profile.txt',depths,signal)


  call book%record(pname,trace,'Create profile')
  call prof%new(depths,signal)


  call book%record(pname,trace,'Estimate anomaly')
  call prof%set_zmin(az1)
  call prof%set_zmax(az2)
  call prof%anomaly()


  call book%record(pname,trace,'Create the SVD object')
  call svd%new(alpha*cu,nw,logy,time,prof)


  call book%record(pname,trace,'Invert the anomaly')
  call svd%inversion()


  call book%record(pname,trace,'Write the results')
  call svd%put(ofile,11_iprec)
  call svd%delete()


  call book%record(pname,trace,'End')
  call book%record(pname,error,'Status = 1')

contains
  function count_lines(u,ifile) result(nlines)
    ! Funtion to count the lines of .dat files.
    ! - u :: unit to open the file to be analyzed (integer)
    ! - ifile :: name of the file of interest (character)
    ! - nlines :: number of lines of the file (integer)
    integer(iprec), intent(in) :: u
    character(len=*), intent(in) :: ifile
    integer(iprec) :: nlines

    integer :: io ! Dummy integers

    open(u,file=ifile, iostat=io, action="read")
    if (io/=0) stop 'Cannot open file! '

    nlines = 0
    do
      read(u,*,iostat=io)
      if (io/=0) exit
      nlines = nlines + 1
    end do
    close(u)

    return
  end function count_lines

  subroutine read_profile(u,ifile,x,y)
    ! Subroutine to read files containing 2D data
    ! - u :: unit to open the file to be read (integer)
    ! - ifile :: name of the file (character)
    ! - x :: array with the data in the first dimension (real(:))
    ! - y :: array with the data in the second dimension (real(:))
    character(len=*), intent(in) :: ifile
    integer(iprec), intent(in) :: u
    real(rprec), allocatable, intent(out) :: x(:), y(:)

    integer(iprec) :: n ! Dummy integers
    integer(iprec) :: i ! Dummy index

    n = count_lines(u,trim(ifile))

    allocate(x(n))
    allocate(y(n))

    open(u,file=ifile,action="read")
    do i = 1, n
      read(u,*) x(i), y(i)
    end do
    close(u)

    return
  end subroutine read_profile


end program test_svd
