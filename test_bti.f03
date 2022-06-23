!
!----------------------------------------------------------------------
! Program to test the scripts to perform bootstrap inversions
!
! Francisco Jose Cuesta Valero
! 2022-06-22 (Leipzig)
!----------------------------------------------------------------------
program test_bti
  use kinds_module
  use logger_module
  use bti_module
  implicit none

  ! Name of file for saving the results
  character(len = *), parameter :: ofile = "inv_bti.dat"

  ! Logger Information
  type(logger_type) :: book
  integer(iprec), parameter :: log_lvl = warning
  character(len = *), parameter :: pname="test_bti"

  ! Inversion types
  ! Bootstrap Inversion
  type(bti_type) :: inv

  ! Inversion Parameters
  ! Minimum and maximum diffusivity (m2 s-1)
  real(rprec), parameter :: alpha1 = 0.5e-6_rprec
  real(rprec), parameter :: alpha2 = 1.5e-6_rprec

  ! Minimum and maximum conductivity ( W m-1 K-1)
  real(rprec), parameter :: lambda1 = 2.5_rprec
  real(rprec), parameter :: lambda2 = 3.5_rprec

  ! Number of diffusivities and conductivities
  integer(iprec), parameter :: nalpha = 1000_iprec
  integer(iprec), parameter :: nlambda = 1000_iprec

  ! Number of eigenvalues retained in the inversions
  integer(iprec), parameter :: neigen = 2_iprec

  ! Length of the time steps for surface histories (years)
  real(rprec), parameter :: dtime = 30.0_rprec

  ! Number of time steps for surface histories
  integer(iprec), parameter :: ntime = 21_iprec

  ! Minimum and maximum depths for artificial profile (m)
  real(rprec), parameter :: afm1 = 0.1_rprec
  real(rprec), parameter :: afm2 = 300_rprec

  ! Number of depths for artificial profile (m)
  integer(iprec), parameter :: ndepth = 30_iprec

  ! Minimum and maximum depth to estimate the quasi-equilibrium
  ! temperature profile (m)
  real(rprec), parameter :: az1 = 200.0_rprec
  real(rprec), parameter :: az2 = 300.0_rprec

  ! Number of populations in Bootstrap inversions
  integer(iprec), parameter :: n_b = 1e3_iprec

  ! Number of threads for Bootstrap inversions
  integer(iprec), parameter :: nthreads = 2_iprec

  ! Name of the profiles
  character(len=180), allocatable :: logs(:) ! List of boreholes

  ! Logging year of the profiles
  real(rprec), allocatable :: logy(:)

  ! Constants
  real(rprec), parameter :: cu = 3600.0_rprec * 24.0_rprec * 365.25_rprec

  real(rprec), allocatable :: alpha(:), time_series(:), lambda(:), depth(:)
  integer(iprec), allocatable :: eigen(:)

  integer(iprec) :: i, nbore


  call book%new(log_lvl)

  call book%record(pname,trace,'Start')

  call book%record(pname,trace, 'Name and logging year of the profile')
  nbore=1
  allocate(logs(nbore))
  allocate(logy(nbore))
  logs(1) = "profile"
  logy(1) = 2021.0_rprec


  call book%record(pname,trace, 'Create time series for surface history')
  allocate(time_series(ntime))
  do i = 1, ntime
    time_series(i) = real(i,rprec) * dtime
  end do

  call book%record(pname,trace, 'Create array with different values for &
      &diffusivity')
  allocate(alpha(nalpha))
  do i = 1, nalpha
    alpha(i) = uniform_number_real(alpha1,alpha2)
  end do
  alpha = alpha * cu ! so alpha is in the adequate units

  call book%record(pname,trace, 'Create array with different values for &
      &conductivity')
  allocate(lambda(nlambda))
  do i = 1, nlambda
    lambda(i) = uniform_number_real(lambda1,lambda2)
  end do
  lambda  = lambda * cu ! so lambda is in the adequate units

  call book%record(pname,trace, 'Create array with depths for forward &
      &model')
  allocate(depth(ndepth))
  do i = 1, ndepth
    depth(i) = real(i,rprec) * (afm2-afm1)/real(ndepth,rprec)
  end do ! i

  call book%record(pname,trace, 'Create array with eigenvalue')
  allocate(eigen(1))
  eigen(1) = 2_iprec


  call book%record(pname,trace,'Create the BTI object')
  call inv%new(az1,az2,logs,logy,n_b,nthreads,eigen,alpha,time_series,&
      depth=depth,lambda=lambda)

  call book%record(pname,trace,'Invert the anomaly')
  call inv%inversion()

  call book%record(pname,trace,'Write the results')
  call inv%put(ofile,11_iprec)

  call inv%delete()

  call book%record(pname,trace,'End')
  call book%record(pname,error,'Status = 1')

contains
    subroutine init_logger()
    ! This routine is called from the constructor to initialize the logger
    logical, save :: loggerInitialized = .false.

    if ( .not. loggerInitialized ) then
      call book%new(log_lvl)
      loggerInitialized = .true.
    end if
  end subroutine

  function gaussian_number_real(mu,sigma)
    ! Subroutine to generate random numbers normally disttributed
    real(rprec), intent(in) :: mu
    real(rprec), intent(in) :: sigma
    real(rprec) :: gaussian_number_real

    real(rprec), parameter :: pi = 2.0_rprec * asin(1.0_rprec)
    real(rprec) :: r1, r2

    call random_number(r1)
    call random_number(r2)

    gaussian_number_real = mu + sigma * sqrt( - 2.0_rprec * log(r1)) * &
                           cos(2.0_rprec * pi * r2)

    return
  end function gaussian_number_real

  function uniform_number_real(n1,n2)
    ! Function to generate a random number between n1 and n2 (uniform distribution)
    ! - n1 :: lower limit (real)
    ! - n2 :: upper limit (real)
    real(rprec), intent(in) :: n1
    real(rprec), intent(in) :: n2
    real(rprec) :: uniform_number_real

    real(rprec) :: idum

    character(len=*), parameter :: pname = 'mod_bni < uniform_number'

    call init_logger

    if(n1.gt.n2) then
      call book%record(pname,error,"ERROR - n1 larger than n2")
      stop
    elseif(n1.eq.n2) then
      call book%record(pname,warning,"WARNING - n1 equal to n2")
      uniform_number_real = n1
    else
      call random_number(idum)
      uniform_number_real = (n2-n1) * idum + n1
    end if

    return
  end function uniform_number_real
end program test_bti
