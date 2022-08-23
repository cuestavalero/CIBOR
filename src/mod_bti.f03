!
!----------------------------------------------------------------------
! Module to invert profiles using the boostrap technique.
!
! This module needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_reader.f03
! - mod_profile.f03 (it requires mod_regression.f03)
! - mod_fm.f03 (it requires mod_profile.f03)
! - mod_svd.f03 (it requires mod_profile and mod_fm.f03)
! - mod_ghf.f03
! - mod_ensembles.f03
!
! Francisco Jose Cuesta-Valero
! 2021-09-04 (Alhama de Murcia)
!----------------------------------------------------------------------
module bti_module
  use omp_lib ! Parallelization
  use kinds_module
  use logger_module
  use profile_module
  use fm_module
  use svd_module
  use ghf_module
  use ensembles_module
  implicit none

  private

  type, public :: bti_type
    ! Incoming data
    ! Total number of resamplings [1]
    integer(iprec) :: n_total
    ! Number of threads to parallelize [1]
    integer(iprec) :: n_threads
    ! Estimate anomaly using depths from az1 to az2 [m]
    real(rprec) :: az1,az2
    ! Write population of means
    logical :: means
    ! Names of the logs [1]
    character(180), allocatable :: logs(:)
    ! Logging yeras [CE]
    real(rprec), allocatable :: logy(:)
    ! Eigenvalues [1]
    integer(iprec), allocatable :: eigen(:)
    ! Thermal diffusivities [m2 s-1]
    real(rprec), allocatable :: alpha(:)
    ! Time series for surface signal [Years before present]
    real(rprec), allocatable :: time_series(:)
    ! Populations to be printed
    integer(iprec), allocatable :: population(:)

    ! Optional inputs
    ! Thermal conductivity [W m−1 K−1]
    real(rprec), allocatable :: lambda(:)
    ! Subsurface profile for forward model
    real(rprec), allocatable :: depth(:)

    ! Internals
    ! Missing value
    real(rprec) :: rfils

    ! Outputs
    ! Minimum year [Years CE]
    real(rprec) :: min_year
    ! Maximum year [Years CE]
    real(rprec) :: max_year
    ! General time [Years CE]
    real(rprec), allocatable :: gtime(:)
    ! Histogram of log usage [1]
    integer(iprec), allocatable :: hist(:)
    ! Number of different logs each year [1]
    integer(iprec), allocatable :: logyear(:)
    ! General inversion (temp) [K]
    real(rprec), allocatable :: general_temp(:,:)
    ! Number of inversions (temp) [K]
    real(rprec), allocatable :: num_temp(:)
    ! General inversion (flux) [W m-2]
    real(rprec), allocatable :: general_flux(:,:)
    ! Number of inversion (flux) [W m-2]
    real(rprec), allocatable :: num_flux(:)
    ! General forward model (temp) [K]
    real(rprec), allocatable :: general_temp_fm(:,:)
    ! General forward model (flux) [W m-2]
    real(rprec), allocatable :: general_flux_fm(:,:)


  contains
    procedure, private :: new_bti, new_copy_bti
    procedure, private :: log_temp
    generic :: new => new_bti, new_copy_bti
    procedure :: delete => delete_bti
    procedure :: inversion => inversion_bti
    procedure :: put => write_bti
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
  subroutine new_bti(self,az1,az2,logs,logy,n_total,n_threads,eigen,&
      alpha,time_series,depth,lambda,population,means)
    ! Subroutine to initialize the BTI object
    ! - self :: object to be initialized (bti)
    ! - logs :: names for the logs to be considered (char(:))
    ! - logs_dir :: address for the logs' files (char)
    ! - n_total :: number of total resamplings (integer)
    ! - n_samp :: size of individual resamblings (integer)
    ! - n_threads :: number of threads to parallelize (integer)
    ! - eigen :: eigenvalues to be considered (real(:))
    ! - alpha :: thermal diffusivities to be considered (real(:))
    ! - time_series :: time steps for surface signal (real(:))
    ! - depth :: synthetic profile for forward models (real(:))
    ! - lambda :: thermal conductivities to be considered (real(:))
    ! - land :: land surface considered (real)
    class(bti_type), intent(out) :: self
    integer(iprec), intent(in) :: n_total
    integer(iprec), intent(in) :: n_threads
    real(rprec), intent(in) :: az1,az2
    character(180), allocatable, intent(in) :: logs(:)
    integer(iprec), allocatable, intent(in) :: eigen(:)
    real(rprec), allocatable, intent(in) :: logy(:)
    real(rprec), allocatable, intent(in) :: alpha(:)
    real(rprec), allocatable, intent(in) :: time_series(:)
    real(rprec), allocatable, optional, intent(in) :: depth(:)
    real(rprec), allocatable, optional, intent(in) :: lambda(:)
    integer(iprec), allocatable, optional, intent(in) :: population(:)
    logical, optional, intent(in) :: means

    self%rfils = -999.9_rprec

    self%az1 = az1
    self%az2 = az2
    self%n_total = n_total
    self%n_threads = n_threads
    self%logs = logs
    self%logy = logy
    self%eigen = eigen
    self%alpha = alpha
    self%time_series = time_series
    if(present(depth)) self%depth = depth
    if(present(lambda)) self%lambda = lambda
    if(present(population)) self%population = population
    if(present(means)) self%means = means

    return
  end subroutine new_bti

  subroutine new_copy_bti(self,other)
    ! Subroutine to create a new BNI object from a previously existing one
    ! - self :: object to be created (bti)
    ! - other :: object to be copied (bti)
    class(bti_type), intent(out) :: self
    class(bti_type), intent(in) :: other

    self%rfils = other%rfils
    self%az1 = other%az1
    self%az2 = other%az2
    self%logs = other%logs
    self%logy = other%logy
    self%n_total = other%n_total
    self%n_threads = other%n_threads
    self%eigen = other%eigen
    self%alpha = other%alpha
    self%time_series = other%time_series
    self%means = other%means
    if(allocated(other%depth)) self%depth = other%depth
    if(allocated(other%lambda)) self%lambda = other%lambda
    if(allocated(other%population)) self%population = other%population

    return
  end subroutine new_copy_bti

  subroutine delete_bti(self)
    ! Subroutine to delete BNI objects
    ! - self :: object to be deleted (bti)
    class(bti_type), intent(inout) :: self

    if(allocated(self%logs)) deallocate(self%logs)
    if(allocated(self%logy)) deallocate(self%logy)
    if(allocated(self%eigen)) deallocate(self%eigen)
    if(allocated(self%alpha)) deallocate(self%alpha)
    if(allocated(self%time_series)) deallocate(self%time_series)
    if(allocated(self%depth)) deallocate(self%depth)
    if(allocated(self%lambda)) deallocate(self%lambda)
    if(allocated(self%hist)) deallocate(self%hist)
    if(allocated(self%logyear)) deallocate(self%logyear)
    if(allocated(self%general_temp)) deallocate(self%general_temp)
    if(allocated(self%general_flux)) deallocate(self%general_flux)
    if(allocated(self%general_temp_fm)) deallocate(self%general_temp_fm)
    if(allocated(self%general_flux_fm)) deallocate(self%general_flux_fm)
    if(allocated(self%population)) deallocate(self%population)

    return
  end subroutine delete_bti

  ! ------------------------
  ! Fundamental calculations
  ! ------------------------
  subroutine inversion_bti(self)
    ! Subroutine to estimate the inversions
    ! - self :: BNI object to produce estimates (bti)
    class(bti_type), intent(inout) :: self

    ! To transform the profiles
    type(profile_type) :: prof

    ! To obtain SVD inversions
    type(svd_type) :: svd

    ! To obtain forward models
    type(fm_type) :: fm

    ! To obtain flux estimates
    type(ghf_type) :: ghf

    ! Ensemble to save subsampling
    type(ensemble_type) :: ensemble, ensemble_fm
    type(ensemble_type) :: fensemble, fensemble_fm

    ! Ensemble to save the subsampling average
    real(rprec), allocatable :: principal(:,:), principal_fm(:,:)
    real(rprec), allocatable :: fprincipal(:,:), fprincipal_fm(:,:)

    ! Arrays for saving results
    real(rprec), allocatable :: gst(:), x(:), y(:), row(:)
    real(rprec), allocatable :: gtime(:), ay(:), tseries(:), signal(:)
    real(rprec), allocatable :: tarray(:), numinv(:)

    ! Change units
    real(rprec), parameter :: cu = 3600.0_rprec * 24.0_rprec * 365.25_rprec

    ! Critical indices
    integer(iprec) :: ii,ialpha,ieigen,ilog,ilambda

    ! Dummy character
    character(10) :: uu, pp

    ! Dummy integers
    integer(iprec) :: u, i


    character(len=*), parameter :: pname = "mod_bti < inversion_bti"

    call init_logger

    call book%record(pname,trace,'Start subroutine')

    ! Temperature inversion
    if(allocated(self%alpha)) then
      call book%record(pname,trace,"Temp Inversions")

      ! Create general time axis
      self%max_year = maxval(self%logy)
      self%min_year = minval(self%logy) - maxval(self%time_series)
      allocate(gtime(int(self%max_year,iprec)-int(self%min_year,iprec)+1_iprec))
      do i = 1, (int(self%max_year,iprec)-int(self%min_year,iprec)+1_iprec)
        gtime(i) = self%max_year - real(i,rprec) + 1.0_rprec
      end do ! i


      ! Perform the large sampling
      allocate(principal(self%n_total,size(gtime)))
      if(allocated(self%depth)) then
        ! Perform the large sampling
        allocate(principal_fm(self%n_total,size(self%depth)))
      end if

      if(allocated(self%lambda)) then
        ! Perform the large sampling
        allocate(fprincipal(self%n_total,size(gtime)))

        if(allocated(self%depth)) then
          ! Perform the large sampling
          allocate(fprincipal_fm(self%n_total,size(self%depth)))
        end if
      end if

      !$omp parallel num_threads(self%n_threads) private(u,uu,&
      !$omp ensemble,ensemble_fm,ialpha,ieigen,ilog,x,y,&
      !$omp prof,row,gst,svd,fm,tseries,signal,&
      !$omp fensemble,ilambda,fensemble_fm,ghf,ay,pp)

      !u = 90
      u = omp_get_thread_num() + 10_iprec
      call book%record(pname,debug,'Unit - ',u)
      write(uu,'(i2)') u

      !$omp do
      do ii = 1, self%n_total

        call ensemble%new(self%rfils)
        if(allocated(self%depth)) then
          call ensemble_fm%new(self%rfils)
        end if


        if(allocated(self%lambda)) then
          call fensemble%new(self%rfils)

          if(allocated(self%depth)) then
            call fensemble_fm%new(self%rfils)
          end if
        end if

        ! Perform the small sampling
        do ilog = 1, size(self%logs)
          call book%record(pname,debug,'Diffusivity')
          ialpha = uniform_number_integer(1,size(self%alpha))

          call book%record(pname,debug,'Eigenvalue')
          ieigen = uniform_number_integer(1,size(self%eigen))

          call book%record(pname,debug,'Log')

          call book%record(pname,debug,'Read the profile '//trim(self%logs(ilog)))
          call book%record(pname,debug,'Opening '//trim(self%logs(ilog))&
                           //'_'//trim(uu)//'.txt')
          call read_profile(u,trim(self%logs(ilog))&
                          //'_'//trim(uu)//'.txt',x,y)

          ! Test the data
          if(test_array(x).eq.0_iprec) then
            call book%record(pname,error,"ERROR - inconsistencies in"&
                                           //trim(self%logs(ilog))//" -> depth")
            stop
          end if
          if(test_array(y).eq.0_iprec) then
            call book%record(pname,error,"ERROR - inconsistencies in"&
                                           //trim(self%logs(ilog))//" -> temp")
            stop
          end if

          call book%record(pname,debug,'Initialize temp profile '&
              //trim(self%logs(ilog)))


          call prof%new(x,y)
          call book%record(pname,debug,'Estimate temp anomaly profiles')
          call prof%set_zmin(self%az1)
          call prof%set_zmax(self%az2)
          call prof%anomaly()


          call book%record(pname,debug,'Estimate anomaly')
          call constrain_profile(prof,2.0_rprec,ay)


          call book%record(pname,debug,'Perform inversion')
          call prof%set_anoma(ay) ! Set the anomaly to be inverted


          ! Invert the anomaly
          call svd%new(self%alpha(ialpha),self%eigen(ieigen),self%logy(ilog),&
                       self%time_series,prof)
          call svd%inversion()


          ! Sort the years
          call general_time(self%rfils,gtime,svd%inv(:,1),svd%inv(:,3),row)
          call ensemble%add(row)
          deallocate(row)


          if(allocated(self%depth)) then
            tseries = self%logy(ilog) - svd%inv(:,1)
            signal = svd%inv(:,3)
            call fm%new(self%alpha(ialpha),tseries,signal,self%depth)
            call fm%forward_model()
            call ensemble_fm%add(fm%prof%atemp_original)
            call fm%delete()
            deallocate(tseries)
            deallocate(signal)
          end if

          call prof%delete()

          if(allocated(self%lambda)) then
            call book%record(pname,debug,'Perform flux inversion')

            call book%record(pname,debug,'Diffusivity')
            ilambda = uniform_number_integer(1,size(self%lambda))

            call book%record(pname,trace,'GHF from Wang and Brass')
            call ghf%new(self%alpha(ialpha),self%lambda(ilambda))
            call ghf%flux(svd%oinv(size(svd%oinv,1):1:-1,2),&
                svd%logy-svd%time_series(size(svd%time_series,1):1:-1),cu,&
                svd%logy)

            if(test_array(ghf%gflux(:,2)).eq.0) then
              call book%record(pname,error,'ERROR - something wrong in flux estimate')
              call book%record(pname,error,'Profile: '//trim(self%logs(ilog)))
              call book%record(pname,error,'Alpha: ',self%alpha(ialpha))
              call book%record(pname,error,'Lambda: ',self%lambda(ilambda))
              call book%record(pname,error,'Times: ',ghf%gflux(:,1))
              call book%record(pname,error,'Flux: ',ghf%gflux(:,2))
              stop
            end if

            call general_time(self%rfils,gtime,ghf%gflux(:,1),ghf%gflux(:,2),row)
            call fensemble%add(row)
            deallocate(row)



            if(allocated(self%depth)) then
              tseries = self%logy(ilog) - ghf%gflux(:,1)
              tseries = tseries(size(tseries,1):1:-1)
              signal = ghf%gflux(size(ghf%gflux,1):1:-1,2)
              call fm%new(self%alpha(ialpha),tseries,signal,self%depth)
              call fm%forward_model()
              call fensemble_fm%add(fm%prof%atemp_original)
              call fm%delete()
              deallocate(tseries)
              deallocate(signal)
            end if
          end if ! lambda

          call svd%delete()
          call ghf%delete()

        end do ! ilog


        call book%record(pname,debug,'Save mean to principal ensemble')
        call ensemble%mean_cols(gst)
        principal(ii,:) = gst
        deallocate(gst)

        ! Save some populations to check quality of inversions
        if(allocated(self%population)) then
          if(count(self%population .eq. ii).gt.0) then
            write(pp,'(i8)') ii
            pp = adjustl(pp)
            if(size(self%alpha,1).eq.1) then
              call write_matrix(90,gtime(:),ensemble%ensemble(:,:),&
                  'temp_matrix_a0_'//trim(pp)//'.dat')
            else
              call write_matrix(90,gtime(:),ensemble%ensemble(:,:),&
                  'temp_matrix_an_'//trim(pp)//'.dat')
            end if
          end if
        end if

        call ensemble%delete()


        if(allocated(self%depth)) then
          call ensemble_fm%mean_cols(gst)
          principal_fm(ii,:) = gst
          deallocate(gst)

          ! Save some populations to check quality of inversions
          if(allocated(self%population)) then
            if(count(self%population .eq. ii).gt.0) then
              write(pp,'(i8)') ii
              pp = adjustl(pp)
              if(size(self%alpha,1).eq.1) then
                call write_matrix(91,self%depth(:),&
                    ensemble_fm%ensemble(:,:),&
                    'temp_fm_matrix_a0_'//trim(pp)//'.dat')
              else
                call write_matrix(91,self%depth(:),&
                    ensemble_fm%ensemble(:,:),&
                    'temp_fm_matrix_an_'//trim(pp)//'.dat')
              end if
            end if
          end if

          call ensemble_fm%delete()


        end if


        if(allocated(self%lambda)) then
          call fensemble%mean_cols(gst)
          fprincipal(ii,:) = gst

          ! Save some populations to check quality of inversions
          if(allocated(self%population)) then
            if(count(self%population .eq. ii).gt.0) then
              write(pp,'(i8)') ii
              pp = adjustl(pp)
              if(size(self%alpha,1).eq.1) then
                call write_matrix(92,gtime(:),fensemble%ensemble(:,:),&
                    'flux_matrix_a0_'//trim(pp)//'.dat')
              else
                call write_matrix(92,gtime(:),fensemble%ensemble(:,:),&
                    'flux_matrix_an_'//trim(pp)//'.dat')
              end if
            end if
          end if
          deallocate(gst)

          call fensemble%delete()


          if(allocated(self%depth)) then
            call fensemble_fm%mean_cols(gst)
            fprincipal_fm(ii,:) = gst
            deallocate(gst)

            ! Save some populations to check quality of inversions
            if(allocated(self%population)) then
              if(count(self%population .eq. ii).gt.0) then
                write(pp,'(i8)') ii
                pp = adjustl(pp)
                if(size(self%alpha,1).eq.1) then
                  call write_matrix(93,self%depth(:),&
                      fensemble_fm%ensemble(:,:),&
                      'flux_fm_matrix_a0_'//trim(pp)//'.dat')
                else
                  call write_matrix(93,self%depth(:),&
                      fensemble_fm%ensemble(:,:),&
                      'flux_fm_matrix_an_'//trim(pp)//'.dat')
                end if
              end if
            end if

            call fensemble_fm%delete()


          end if


        end if ! lambda
      end do ! ii
      !$omp end do

      !$omp end parallel

      call book%record(pname,debug,'Save inversion results')
      call manage_results(self%rfils,gtime,principal,self%general_temp,&
          self%gtime,self%num_temp)

      ! Save the population of means to check quality of inversions
      if(self%means) then
        if(size(self%logs,1).eq.1) then
          if(size(self%alpha,1).eq.1) then
            call write_matrix(94,gtime(:),principal(:,:),&
                'temp_matrix_a0_'//trim(self%logs(1))//'.dat')
          else
            call write_matrix(94,gtime(:),principal(:,:),&
                'temp_matrix_an_'//trim(self%logs(1))//'.dat')
          end if
        else
          if(size(self%alpha,1).eq.1) then
            call write_matrix(94,gtime(:),principal(:,:),&
                'temp_matrix_a0.dat')
          else
            call write_matrix(94,gtime(:),principal(:,:),&
                'temp_matrix_an.dat')
          end if
        end if
      end if


      if(allocated(self%depth)) then
        call book%record(pname,debug,'Save forward model results')
        call manage_results(self%rfils,self%depth,principal_fm,&
            self%general_temp_fm,tarray,numinv)

        ! Save the population of means to check quality of inversions
        if(self%means) then
          if(size(self%logs,1).eq.1) then
            if(size(self%alpha,1).eq.1) then
              call write_matrix(95,self%depth(:),principal_fm(:,:),&
                  'temp_fm_matrix_a0_'//trim(self%logs(1))//'.dat')
            else
              call write_matrix(95,self%depth(:),principal_fm(:,:),&
                  'temp_fm_matrix_an_'//trim(self%logs(1))//'.dat')
            end if
          else
            if(size(self%alpha,1).eq.1) then
              call write_matrix(95,self%depth(:),principal_fm(:,:),&
                  'temp_fm_matrix_a0.dat')
            else
              call write_matrix(95,self%depth(:),principal_fm(:,:),&
                  'temp_fm_matrix_an.dat')
            end if
          end if
        end if

      end if


      if(allocated(self%lambda)) then
        call book%record(pname,debug,'Save flux inversion results')
        call manage_results(self%rfils,gtime,fprincipal,self%general_flux,&
            tarray,self%num_flux)

        ! Save the population of means to check quality of inversions
        if(self%means) then
          if(size(self%logs,1).eq.1) then
            if(size(self%alpha,1).eq.1) then
              call write_matrix(96,gtime(:),fprincipal(:,:),&
                  'flux_matrix_a0_'//trim(self%logs(1))//'.dat')
            else
              call write_matrix(96,gtime(:),fprincipal(:,:),&
                  'flux_matrix_an_'//trim(self%logs(1))//'.dat')
            end if
          else
            if(size(self%alpha,1).eq.1) then
              call write_matrix(96,gtime(:),fprincipal(:,:),&
                  'flux_matrix_a0.dat')
            else
              call write_matrix(96,gtime(:),fprincipal(:,:),&
                  'flux_matrix_an.dat')
            end if
          end if
        end if


        if(allocated(self%depth)) then
          call book%record(pname,debug,'Save forward model flux results')
          call manage_results(self%rfils,self%depth,fprincipal_fm,&
              self%general_flux_fm,tarray,numinv)

          ! Save the population of means to check quality of inversions
          if(self%means) then
            if(size(self%logs,1).eq.1) then
              if(size(self%alpha,1).eq.1) then
                call write_matrix(97,self%depth(:),fprincipal_fm(:,:),&
                    'flux_fm_matrix_a0_'//trim(self%logs(1))//'.dat')
              else
                call write_matrix(97,self%depth(:),fprincipal_fm(:,:),&
                    'flux_fm_matrix_an_'//trim(self%logs(1))//'.dat')
              end if
            else
              if(size(self%alpha,1).eq.1) then
                call write_matrix(97,self%depth(:),fprincipal_fm(:,:),&
                    'flux_fm_matrix_a0.dat')
              else
                call write_matrix(97,self%depth(:),fprincipal_fm(:,:),&
                    'flux_fm_matrix_an.dat')
              end if
            end if
          end if

        end if


      end if ! lambda


      ! Estimate number of times each log has been used
      call self%log_temp()

    end if ! alpha

    call book%record(pname,trace,'End subroutine')

    return
  end subroutine inversion_bti

  subroutine manage_results(rfils,garray,imatrix,omatrix,tarray,numinv)
    ! Subroutine to organize the results of temp, flux, ftemp, fflux
    ! - rfils :: filling value of matrix
    ! - garray :: array of general years/depths
    ! - imatrix :: matrix of inversions/forward models
    ! - omatrix :: matrix with the percentiles of inversions/forward models
    ! - tarray :: array with years/depths covering inversions/ forward models
    ! - numinv :: number of inversions/forward models considered
    real(rprec), intent(in) :: rfils
    real(rprec), allocatable, intent(in) :: garray(:)
    real(rprec), intent(in) :: imatrix(:,:)
    real(rprec), allocatable, intent(out) :: omatrix(:,:)
    real(rprec), allocatable, intent(out) :: tarray(:)
    real(rprec), allocatable, intent(out) :: numinv(:)

    real(rprec), allocatable :: x1(:), x2(:), x3(:), num(:)
    integer(iprec) :: i1, i2

    integer(iprec) :: i, j

    character(len=*), parameter :: pname = 'mod_bti <- manage_results'


    call init_logger

    ! Percentiles from inversions
    allocate(x1(size(imatrix,2)))
    allocate(x2(size(imatrix,2)))
    allocate(x3(size(imatrix,2)))
    do i = 1, size(imatrix,2)
      x1(i) = quantile(rfils,imatrix(:,i),0.025_rprec)
      x2(i) = quantile(rfils,imatrix(:,i),0.5_rprec)
      x3(i) = quantile(rfils,imatrix(:,i),0.975_rprec)
    end do ! i

    ! Number of inversions considered per year
    allocate(num(size(imatrix,2)))
    num = 0.0_rprec
    do i = 1, size(imatrix,2)
      do j = 1, size(imatrix,1)
        if(imatrix(j,i).ne.rfils) then
          num(i) = num(i) + 1.0_rprec
        end if
      end do
      if(num(i).eq.0.0_rprec) then
        num(i) = rfils
      end if
    end do

    ! Clear the missing values
    do i = 1, size(x2)
      if(x2(i).ne.rfils) then
        i1 = i
        exit
      end if
    end do

    do i = size(x2), 1, -1
      if(x2(i).ne.rfils) then
        i2 = i
        exit
      end if
    end do

    if(i2.lt.i1) then
      i = i1
      i1 = i2
      i2 = i
    else if(i1.eq.i2) then
      call book%record(pname,error,"ERROR - i1 is equal to i2")
      stop
    end if

    ! Save the results
    allocate(tarray(i2 - i1 + 1_iprec))
    allocate(numinv(i2 - i1 + 1_iprec))
    allocate(omatrix(3,i2 - i1 + 1_iprec))
    tarray(:) = garray(i1:i2)
    omatrix(1,:) = x1(i1:i2)
    omatrix(2,:) = x2(i1:i2)
    omatrix(3,:) = x3(i1:i2)
    numinv(:) = num(i1:i2)
    deallocate(x1)
    deallocate(x2)
    deallocate(x3)
    deallocate(num)

    return
  end subroutine manage_results

  subroutine write_bti(self,ofile,u2)
    ! Subroutine to write the results of the BTI inversion
    ! - u2 :: unit to write the file (optional, integer)
    ! - self :: bootstrap inversion (bti)
    ! - ofile :: name of the output file (.dat)
    ! - inversion ::
    class(bti_type), intent(in) :: self
    integer(iprec), optional, intent(in) :: u2
    character(len=*), intent(in) :: ofile

    character(90) :: f ! Dummy characters

    integer(iprec) :: u ! Dummy integer
    integer(iprec) :: i ! Dummy indices

    if(present(u2)) then
      u = u2
    else
      u = 90_iprec
    end if

    open(u,file=trim(ofile),action="write")

    ! First, relevant parameters of the bootstrap
    f = "(a,i6)"
    write(u,trim(f)) '#n_total', self%n_total
    write(u,trim(f)) '#n_threads', self%n_threads
    write(u,trim(f)) '#n_original_logs', size(self%logs)

    ! Temporal parameters
    f = "(a,es12.3e3)"
    write(u,trim(f)) '#max_year', self%max_year
    write(u,trim(f)) '#min_year', self%min_year

    ! Parameters for selecting the profile
    f = "(a,es12.3e3)"
    write(u,trim(f)) '#az1', self%az1
    write(u,trim(f)) '#az2', self%az2


    write(u,*) ''
    write(u,*) ''

    ! Time series
    f = "(a,es12.3e3)"
    do i = 1, size(self%time_series)
        write(u,trim(f)) '#tseries', self%time_series(i)
    end do


    ! General temp
    if(allocated(self%general_temp)) then
      write(u,*) ''
      write(u,*) ''

      f = "(a,4es12.3e3,2i8)"
      do i = 1, size(self%general_temp,2)
          write(u,trim(f)) '#gtemp', self%gtime(i),&
                                     self%general_temp(1,i),&
                                     self%general_temp(2,i),&
                                     self%general_temp(3,i),&
                                     int(self%num_temp(i),iprec),&
                                     self%logyear(i)
      end do
    end if

    ! General flux
    if(allocated(self%general_flux)) then
      write(u,*) ''
      write(u,*) ''

      f = "(a,4es12.3e3,2i8,)"
      do i = 1, size(self%general_flux,2)
          write(u,trim(f)) '#gflux', self%gtime(i),&
                                     self%general_flux(1,i),&
                                     self%general_flux(2,i),&
                                     self%general_flux(3,i),&
                                     int(self%num_flux(i),iprec),&
                                     self%logyear(i)
      end do
    end if

    ! General temp forward model
    if(allocated(self%general_temp_fm)) then
      write(u,*) ''
      write(u,*) ''

      f = "(a,4es12.3e3)"
      do i = 1, size(self%general_temp_fm,2)
          write(u,trim(f)) '#ftemp', self%depth(i),&
                                     self%general_temp_fm(1,i),&
                                     self%general_temp_fm(2,i),&
                                     self%general_temp_fm(3,i)
      end do
    end if

    ! General flux forawd model
    if(allocated(self%general_flux_fm)) then
      write(u,*) ''
      write(u,*) ''

      f = "(a,4es12.3e3)"
      do i = 1, size(self%general_flux_fm,2)
          write(u,trim(f)) '#fflux', self%depth(i),&
                                     self%general_flux_fm(1,i),&
                                     self%general_flux_fm(2,i),&
                                     self%general_flux_fm(3,i)
      end do
    end if

    close(u)

    return
  end subroutine write_bti

  subroutine write_matrix(u,array,matrix,fname)
    ! Subroutine to print matrixces.
    ! - u :: unit for printing table. It should be different to the used in    !     the parallel area (integer).
    ! - array :: reference array. It indicates the year or depth of each
    !     column in the matrix.
    ! - matrix :: matrix with the inversions to be printed.
    ! - fname :: name for the file containing matrix and array data.
    real(rprec), intent(in) :: array(:)
    real(rprec), intent(in) :: matrix(:,:)
    integer(iprec), optional, intent(in) :: u
    character(len=*), intent(in) :: fname

    character(len=4) :: tj
    integer(iprec) :: i, j, nj, uu

    character(len=*), parameter :: pname = 'mod_bti < write_matrix'


    call init_logger

    if(present(u)) then
      uu = u
    else
      uu = 90
    end if

    if(size(array,1).ne.size(matrix,2)) then
      call book%record(pname,error,"ERROR - array dimension different"&
          &" than matrix axes")
      call book%record(pname,error,"N array = ",size(array,1))
      call book%record(pname,error,"Nj matrix = ",size(matrix,2))
      call book%record(pname,error,"Ni matrix = ",size(matrix,1))
      stop 'ERROR - array dimension different than matrix axes'
    end if

    nj = size(matrix,2)
    write(tj,'(i4)') nj
    tj = adjustl(tj)

    open(90,file=fname)
    write(90,'('//trim(tj)//'es12.3e3)') (array(j), j=1,nj)
    do i = 1, size(matrix,1)
      write(90,'('//trim(tj)//'es12.3e3)') (matrix(i,j), j=1,nj)
    end do
    close(90)

    return
  end subroutine

  ! ------------------------
  ! Auxiliar calculations
  ! ------------------------
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

  function uniform_number_real(n1,n2)
    ! Function to generate a random number between n1 and n2 (uniform distribution)
    ! - n1 :: lower limit (real)
    ! - n2 :: upper limit (real)
    real(rprec), intent(in) :: n1
    real(rprec), intent(in) :: n2
    real(rprec) :: uniform_number_real

    real(rprec) :: idum

    character(len=*), parameter :: pname = 'mod_bti < uniform_number_real'

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

  function uniform_number_integer(n1,n2)
    ! Function to generate a random number between n1 and n2 (uniform distribution)
    ! - n1 :: lower limit (integer)
    ! - n2 :: upper limit (integer)
    integer(iprec), intent(in) :: n1
    integer(iprec), intent(in) :: n2
    integer(iprec) :: uniform_number_integer

    real(rprec) :: idum

    character(len=*), parameter :: pname = 'mod_bti < uniform_number_integer'

    call init_logger

    if(n1.gt.n2) then
      call book%record(pname,error,"ERROR - n1 larger than n2")
      stop
    elseif(n1.eq.n2) then
      call book%record(pname,warning,"WARNING - n1 equal to n2")
      uniform_number_integer = n1
    else
      call random_number(idum)
      uniform_number_integer = floor(real((n2-n1+1_iprec),rprec) * idum) + n1
    end if

    return
  end function uniform_number_integer

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

  subroutine general_time(rfils,gtime,tseries,temp,y)
    ! Subroutine to sort the time of each log inversion with relation to the genral time axis
    ! - rfils :: filling value (real)
    ! - gtime :: general time axis (real(:))
    ! - tseries :: years of inversion (real(:))
    ! - temp :: annual temperature series (real(:))
    real(rprec), intent(in) :: rfils
    real(rprec), intent(in) :: gtime(:), tseries(:), temp(:)
    real(rprec), allocatable, intent(out) :: y(:)

    ! Dummy integers
    integer(iprec) :: nt, n

    ! Dummy indices
    integer(iprec) :: i, j

    nt = size(gtime)
    n = size(tseries)

    allocate(y(nt))
    y = rfils

    do i = 1, nt ! Look for the corresponding general year
      do j = 1, n ! Search through each particular year
        if(gtime(i).eq.tseries(j)) then
          y(i) = temp(j)
          exit
        end if
      end do ! j
    end do ! i

    return
  end subroutine general_time

  subroutine constrain_profile(prof,sig,ty)
    ! Subroutine to obtain a profile to a certain spread sig
    ! - prof :: profile type with the anomalies and the regression coefficientes
    ! - sig :: value, in sigma units, to limit the spread of the coefficients
    ! - ty :: array for saving the anomaly
    class(profile_type), intent(in) :: prof
    real(rprec), intent(in) :: sig
    real(rprec), allocatable, intent(out) :: ty(:)

    ! Dummy reals
    real(rprec) :: gam, t0


    if(sig.eq.0.0_rprec) then
      ty = prof%atemp_original
    else

      call select_slope_inter(prof,sig,t0,gam)

      ty = prof%temp(:) - ( t0 + gam * prof%depth(:) )

    end if

    return
  end subroutine constrain_profile

  subroutine select_slope_inter(prof,sig,t0,gam)
    ! Subroutine to select gaussian slopes and intercepts in agreement with
    !     the two extremal anomaly profiles
    ! - prof :: profile object containing all information
    ! - sig :: parameter for limiting the spread in slope and intercept
    ! - t0 :: value of the intercept result
    ! - gam :: value of the slope result
    class(profile_type), intent(in) :: prof
    real(rprec), intent(in) :: sig
    real(rprec), intent(out) :: t0
    real(rprec), intent(out) :: gam

    ! Look for candidates
    integer(iprec), parameter :: n_candi = 1000_iprec

    real(rprec) :: ct0, cgam ! Candidates
    real(rprec) :: th, zh ! H point
    !real(rprec) :: t_minus, t_plus

    integer(iprec) :: k

    character(len=*), parameter :: pname = "mod_bti < select_slope_inter"

    call init_logger

    !t_minus = prof%eq_profile%inter - sig * prof%eq_profile%s_inter
    !t_plus = prof%eq_profile%inter + sig * prof%eq_profile%s_inter


    ! Look for candidates to t0 and gam
    do k = 1, n_candi
      ct0 = gaussian_number_real(prof%eq_profile%inter, &
          prof%eq_profile%s_inter)
      !ct0 = uniform_number_real(prof%eq_profile%inter - sig * &
      !    prof%eq_profile%s_inter, prof%eq_profile%inter + sig * &
      !    prof%eq_profile%s_inter)

      !if(ct0.ge.t_minus.and.ct0.le.t_plus) then
      !  t0 = ct0
      !  exit
      !end if

      t0 = ct0
      exit

      if(k.eq.n_candi) then
        call book%record(pname,error,"ERROR: no good value for inter")
        stop "No good value for inter"
      end if
    end do

    call prof%point_h(sig,th,zh)
    !call prof%point_b(th,zh)

    cgam = ( th - ct0 ) / zh


    gam = cgam

    return
  end subroutine select_slope_inter

  !subroutine temp2flux(lambda,x,y,fx,fy)
  !  ! Subrotine to estimate the flux profile
  !  ! lambda :: thermal conductivity (real)
  !  ! x :: depths (real)
  !  ! y :: temperatures (real)
  !  ! fx :: depths for flux profile (real)
  !  ! fy :: flux profile (real)
  !  real(rprec), intent(in) :: x(:), y(:)
  !  real(rprec), intent(in) :: lambda
  !  real(rprec), allocatable, intent(out) :: fx(:), fy(:)

  !  integer(iprec) :: n ! Dummy integer
  !  integer(iprec) :: i, j ! Dummy index

  !  character(len=*), parameter :: pname = "mod_bti < temp2flux"

  !  call init_logger

  !  n = size(y)
  !  allocate(fx(n-1))
  !  allocate(fy(n-1))

  !  j = 0_iprec
  !  do i = 1, n-1
  !    if(x(i+1)-x(i).ne.0.0_rprec) then
  !      j = j + 1_iprec
  !      fx(j) = x(i) + ( x(i+1) - x(i) )/2.0_rprec
  !      fy(j) = - lambda * ( y(i+1) - y(i) ) / ( x(i+1) - x(i) )
  !    else
  !      call book%record(pname,warning,"WARNING - There is an infty at index ",i)
  !    end if
  !    !print *, x(i+1)-x(i),fy(i),i
  !  end do

  !  return
  !end subroutine temp2flux

  real(rprec) function quantile(rfils,xx,alpha)
    implicit none

    real(rprec), intent(in) :: alpha, rfils
    real(rprec), intent(in) :: xx(:)

    real(rprec), allocatable :: x(:), p(:)

    integer(iprec) :: n

    real(rprec) :: gam
    integer(iprec) :: i, j, k0

    character(len=*), parameter :: pname = 'mod_bti < quantile'

    call init_logger

    n = 0_iprec
    do i = 1, size(xx)
      if(xx(i).ne.rfils) then
        n = n + 1_iprec
      end if
    end do ! i

    if(n.eq.0_iprec) then
      quantile = rfils
      return
    end if

    allocate(x(n))
    j = 0_iprec
    do i = 1, size(xx)
      if(xx(i).ne.rfils) then
        j = j + 1_iprec
        x(j) = xx(i)
      end if
    end do ! i

    allocate(p(n))

    ! Sort the data
    call sort(n,x)

    ! Vector of percentiles
    do i = 1, n
      p(i) = ( real(i,rprec) - 1.0_rprec/3.0_rprec ) / ( real(n,rprec) + 1.0_rprec/3.0_rprec )
    end do ! i

    ! Find k0
    k0 = 0_iprec
    do i = 1, (n-1)
      if(p(i).le.alpha.and.p(i+1).gt.alpha) then
        k0 = i
        exit
      end if
    end do ! i

    if(k0.eq.0_iprec) then
      call book%record(pname,warning,'WARNING - cannot find a k0 index.&
                                & q defined as closest percentile')
      if(alpha.le.p(1)) then
        k0 = 1_iprec
        call book%record(pname,warning,'alpha, p(1) = ',[alpha,p(1)])
      end if
      if(alpha.ge.p(n)) then
        k0 = n - 1_iprec ! To find gam (see the k0+1 bit below)
        call book%record(pname,warning,'alpha, p(n) = ',[alpha,p(n)])
      end if
    end if

    ! Estimate gamma
    gam = ( alpha - p(k0) ) / ( p(k0+1) - p(k0) )

    ! Quantile
    quantile = ( 1.0_rprec - gam ) * x(k0) + gam * x(k0)

    return
  end function quantile

  ! ------------------------
  ! Quality tests
  ! ------------------------
  function test_array(x) result(t)
    ! Function to check that the individual inversions are performed
    !          without strange results.
    ! - inversion => the inversion to be checked (record)
    ! - tol => maxmum temperature difference allowed between time
    !          steps (real)
    ! - t => 1 if all test are successful, 0 otherwise
    real(rprec), intent(in) :: x(:)
    integer(iprec) :: t

    character(len=7) :: c1,c2,c3 ! Dummy characters for logger

    integer(iprec) :: i, ind ! Dummy indice

    character(len=*), parameter :: pname = "mod_bti < test_array"

    call init_logger

    call book%record(pname,trace,"Starting function")

    t = 1_iprec ! Initialize test

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

  subroutine log_temp(self)
    ! Subroutine to estimate the number of different logs contributing each year
    ! - self :: inversion object to be analyzed and to save results.
    class(bti_type), intent(inout) :: self

    integer(iprec) :: n, nt
    integer(iprec) :: l, t

    n = size(self%logs)
    nt = size(self%gtime)

    allocate(self%logyear(nt))
    self%logyear = 0_iprec

    do l = 1, n
      do t = 1, nt
        if(self%logy(l).ge.self%gtime(t).and.( self%logy(l) - &
           maxval(self%time_series) ).le.self%gtime(t)) then
          self%logyear(t) = self%logyear(t) + 1_iprec
        end if
      end do ! t
    end do ! l

    return
  end subroutine log_temp


  ! ------------------------
  ! Numerical Recipes
  ! ------------------------
  subroutine sort(n,arr)
      integer(iprec), intent(in) :: n
      real(rprec), dimension(n), intent(inout) :: arr

      integer(iprec), parameter :: M=7, NSTACK=50
      integer(iprec) :: i,ir,j,jstack,k,l
      integer(iprec), dimension(NSTACK) :: istack
      real(rprec) :: a,temp
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 12 j=l+1,ir
          a=arr(j)
          do 11 i=j-1,l,-1
            if(arr(i).le.a)goto 2
            arr(i+1)=arr(i)
11        continue
          i=l-1
2         arr(i+1)=a
12      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        temp=arr(k)
        arr(k)=arr(l+1)
        arr(l+1)=temp
        if(arr(l).gt.arr(ir))then
          temp=arr(l)
          arr(l)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l+1).gt.arr(ir))then
          temp=arr(l+1)
          arr(l+1)=arr(ir)
          arr(ir)=temp
        endif
        if(arr(l).gt.arr(l+1))then
          temp=arr(l)
          arr(l)=arr(l+1)
          arr(l+1)=temp
        endif
        i=l+1
        j=ir
        a=arr(l+1)
3       continue
          i=i+1
        if(arr(i).lt.a)goto 3
4       continue
          j=j-1
        if(arr(j).gt.a)goto 4
        if(j.lt.i)goto 5
        temp=arr(i)
        arr(i)=arr(j)
        arr(j)=temp
        goto 3
5       arr(l+1)=arr(j)
        arr(j)=a
        jstack=jstack+2
        if(jstack.gt.NSTACK) stop 'NSTACK too small in sort'
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
  end subroutine sort

end module bti_module
