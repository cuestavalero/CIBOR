!
!----------------------------------------------------------------------
! Program to test the ppi module
!
! This program needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_profile.f03 (it requires mod_regression.f03)
! - mod_ppi.f03 (it requires mod_svd.f03 and mod_fm.f03)
!
! Francisco Jose Cuesta Valero
! 2021-09-03 (Alhama de Murcia)
!----------------------------------------------------------------------
program test_ppi
  use kinds_module
  use logger_module
  use profile_module
  use ppi_module
  implicit none

  ! Name of file for saving the results
  character(len = *), parameter :: ofile = "inv_ppi.dat"

  ! Logger Information
  type(logger_type) :: book
  integer(iprec), parameter :: log_lvl = warning
  character(len = *), parameter :: pname="test_ppi"

  ! Inversion type
  type(ppi_type) :: inv

  ! Profile type
  type(profile_type) :: prof

  ! Inversion Parameters
  ! Logging year of the profile
  real(rprec), parameter :: logy = 2021.0_rprec

  ! Minimum and maximum diffusivity (m2 s-1)
  real(rprec), parameter :: alpha1 = 0.5e-6_rprec
  real(rprec), parameter :: alpha2 = 1.5e-6_rprec

  ! Number of diffusivities and conductivities
  integer(iprec), parameter :: nalpha = 100_iprec

  ! Number of eigenvalues retained in the inversions
  integer(iprec), parameter :: neigen = 2_iprec

  ! Length of the time steps for surface histories (years)
  real(rprec), parameter :: dtime = 30.0_rprec

  ! Number of time steps for surface histories
  integer(iprec), parameter :: ntime = 21_iprec

  ! Maximum change of temperature between two consecutive time steps (C)
  real(rprec), parameter :: tol_temp = 100.0_rprec

  ! Measurement error (approaximated) (C)
  real(rprec), parameter :: er = 0.05_rprec

  ! Minimum and maximum depth to estimate the quasi-equilibrium
  ! temperature profile (m)
  real(rprec), parameter :: az1 = 200.0_rprec
  real(rprec), parameter :: az2 = 300.0_rprec

  ! Constants
  real(rprec), parameter :: cu = 3600.0_rprec * 24.0_rprec * 365.25_rprec

  integer(iprec) :: w(2)
  real(rprec), allocatable :: alpha(:), time_series(:)
  real(rprec), allocatable :: x(:), y(:)

  integer(iprec) :: i

  call book%new(log_lvl)

  call book%record(pname,trace,'Start')

  call book%record(pname,trace, 'Create time series for surface history')
  allocate(alpha(nalpha))
  do i = 1, nalpha
    alpha(i) = alpha1 + real(i-1,rprec) * (alpha2-alpha1)/real(nalpha,rprec)
  end do


  call book%record(pname,trace, 'Create time series for surface history')
  allocate(time_series(ntime))
  do i = 1, ntime
    time_series(i) = real(i,rprec) * dtime
  end do


  call book%record(pname,trace, 'Create array with eigenvalue')
  w = [neigen, neigen]


  call book%record(pname,trace,'Read profile')
  call read_profile(10_iprec,'profile.txt',x,y)


  call book%record(pname,trace,'Estimate anomaly')
  call prof%new(x,y)
  call prof%set_zmin(az1)
  call prof%set_zmax(az2)
  call prof%anomaly()


  call book%record(pname,trace,'Create the PPI object')
  call inv%new(prof,logy,time_series,w,alpha*cu,tol_temp,er)


  call book%record(pname,trace,'Invert the anomaly')
  call inv%inversion()


  call book%record(pname,trace,'Write the results')
  call inv%put(ofile,11_iprec)


  call inv%delete()

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
end program test_ppi
