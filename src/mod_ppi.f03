!
!----------------------------------------------------------------------
! Module to invert profiles.
!
! This module needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_profile.f03 (it requires mod_regression.f03)
! - mod_fm.f03 (it requires mod_profile.f03)
! - mod_svd.f03 (it requres mod_fm.f03)
!
! Francisco Jose Cuesta-Valero
! 2021-09-02 (Alhama de Murcia)
!----------------------------------------------------------------------

module ppi_module
  use kinds_module
  use logger_module
  use profile_module
  use fm_module
  use svd_module
  implicit none

  private

  type, public :: model_type
    real(rprec) :: logy ! Logging year [C.E.]
    real(rprec), allocatable :: time_series(:) ! in years before present
    integer(iprec) :: w(2) ! minimum (1) and maximum (2) eigenvalues
    real(rprec), allocatable :: alpha(:) ! series of diffusivities [m2/s]
    real(rprec) :: tol ! maximum temp difference between time steps for a inversion [K]
    real(rprec) :: err ! typical borehole measurement error [K]
  end type

  type, extends(model_type), public :: ppi_type
    type(profile_type) :: prof ! Profile to be inverted
    real(rprec), allocatable :: years(:)
    real(rprec), allocatable :: metadata(:,:) ! Save info for inversions
    real(rprec), allocatable :: svd_solutions(:,:) ! PPI inversions
    real(rprec), allocatable :: fm_solutions(:,:) ! Forward Model results
    real(rprec), allocatable :: general_solution(:,:) ! General solution
    real(rprec), allocatable :: general_fm(:,:) ! General solution
  contains
    procedure :: new_ppi, new_copy_ppi
    generic :: new => new_ppi, new_copy_ppi
    procedure :: delete => delete_ppi
    procedure :: inversion => inversion_ppi
    procedure :: add_inversion
    procedure :: general_inversion
    procedure :: general_forward_model
    procedure :: put => write_ppi
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

  subroutine new_ppi(self,prof,logy,signal_model,w,alpha,tol,err)
    ! SUbroutine to initialize an inversion object
    ! - self :: object to be initialized (ppi)
    ! - prof :: profile to be inverted (profile)
    ! - logy :: logging year (real)
    ! - signal_model :: model to reconstruct the surface signal (real)
    ! - w :: minimum and maximum eigenvalues for the inversion (integer)
    ! - alpha :: difussivities for the inversion (real)
    ! - tol :: maximum temp difference between inversion time steps (real)
    ! - err :: typical error measirements in boreholes (real)
    class(ppi_type), intent(out) :: self
    class(profile_type), intent(in) :: prof
    real(rprec), intent(in) :: logy
    real(rprec), intent(in) :: signal_model(:)
    real(rprec), intent(in) :: alpha(:)
    integer(iprec), intent(in) :: w(2)
    real(rprec), intent(in) :: tol
    real(rprec), intent(in) :: err

    self%model_type%logy = logy
    self%model_type%time_series = signal_model
    self%model_type%alpha = alpha
    self%model_type%w = w
    self%model_type%tol = tol
    self%model_type%err = err
    self%prof = prof

    return
  end subroutine new_ppi

  subroutine new_copy_ppi(self,other)
    ! Subroutine to create a new ppi object from a previously existing one
    ! - self :: object to be created (ppi)
    ! - other :: object to be copied (ppi)
    class(ppi_type), intent(out) :: self
    class(ppi_type), intent(in) :: other

    self%model_type%logy = other%model_type%logy
    self%model_type%time_series = other%model_type%time_series
    self%model_type%alpha = other%model_type%alpha
    self%model_type%w = other%model_type%w
    self%model_type%err = other%model_type%err
    self%model_type%tol = other%model_type%tol
    self%prof = other%prof

    if(allocated(other%metadata)) then
      self%metadata = other%metadata
    end if

    if(allocated(other%svd_solutions)) then
      self%svd_solutions = other%svd_solutions
    end if

    if(allocated(other%fm_solutions)) then
      self%fm_solutions = other%fm_solutions
    end if

    if(allocated(other%years)) then
      self%years = other%years
    end if
    if(allocated(other%general_solution)) then
      self%general_solution = other%general_solution
    end if
    if(allocated(other%general_fm)) then
      self%general_fm = other%general_fm
    end if

    return
  end subroutine new_copy_ppi

  subroutine delete_ppi(self)
    ! Subroutine to delete inversion objects.
    ! - self :: object to be deleted (ppi)
    class(ppi_type), intent(inout) :: self

    if(allocated(self%model_type%time_series)) then
      deallocate(self%model_type%time_series)
    end if

    if(allocated(self%model_type%alpha)) then
      deallocate(self%model_type%alpha)
    end if

    if(allocated(self%metadata)) then
      deallocate(self%metadata)
    end if

    if(allocated(self%svd_solutions)) then
      deallocate(self%svd_solutions)
    end if

    if(allocated(self%fm_solutions)) then
      deallocate(self%fm_solutions)
    end if

    if(allocated(self%years)) then
      deallocate(self%years)
    end if
    if(allocated(self%general_solution)) then
      deallocate(self%general_solution)
    end if
    if(allocated(self%general_fm)) then
      deallocate(self%general_fm)
    end if
    if(allocated(self%prof%depth)) then
      call self%prof%delete()
    end if

    return
  end subroutine delete_ppi

  subroutine delete_model(self)
    ! Subroutine to delete inversion models
    ! - self :: object to be deleted (model)
    class(model_type), intent(inout) :: self

    if(allocated(self%time_series)) deallocate(self%time_series)
    if(allocated(self%alpha)) deallocate(self%alpha)

    return
  end subroutine delete_model

  subroutine add_inversion(self,ptype,alpha,w,wtd,inversion,fm_profile)
    ! Subroutine to add new rows to store more inversions
    ! - self :: object to store the inversion (ppi)
    ! - ptype :: type of anomaly profile (real)
    ! - alpha :: thermal diffusivity (real)
    ! - w :: number of eigenvalues (integer)
    ! - wtd :: associated weight with the inversion (real)
    ! - inversion :: inversion to be included (real)
    ! - fm_profile :: forward model of the corresponding inversion (real)
    class(ppi_type), intent(inout) :: self
    real(rprec), intent(in) :: inversion(:)
    real(rprec), intent(in) :: fm_profile(:)
    real(rprec), intent(in) :: alpha
    real(rprec), intent(in) :: ptype
    integer(iprec), intent(in) :: w
    real(rprec), intent(in) :: wtd

    type(ppi_type) :: other ! Dummy ppi object
    integer(iprec) :: nrow, ncol ! Dummy integer

    character(len=*), parameter :: pname = "mod_ppi < add_inversion"

    call other%new(self)

    ! Metadata
    if(.not.allocated(self%metadata)) then
      allocate(self%metadata(1,4))
      nrow = 0_iprec
    else
      nrow = size(self%metadata,1)
      deallocate(self%metadata)
      allocate(self%metadata(nrow+1,4))
      self%metadata(1:nrow,:) = other%metadata
    end if
    self%metadata(nrow+1,1) = ptype
    self%metadata(nrow+1,2) = alpha
    self%metadata(nrow+1,3) = w
    self%metadata(nrow+1,4) = wtd

    ! SVD Solutions
    if(.not.allocated(self%svd_solutions)) then
      allocate(self%svd_solutions(1,size(inversion)))
      nrow = 0_iprec
      ncol = size(self%svd_solutions,1)
    else
      nrow = size(self%svd_solutions,1)
      ncol = size(self%svd_solutions,2)
      deallocate(self%svd_solutions)
      allocate(self%svd_solutions(nrow+1,ncol))
      self%svd_solutions(1:nrow,:) = other%svd_solutions
    end if
    self%svd_solutions(nrow+1,:) = inversion

    ! Forward Model Solutions
    if(.not.allocated(self%fm_solutions)) then
      allocate(self%fm_solutions(1,size(fm_profile)))
      nrow = 0_iprec
      ncol = size(self%fm_solutions,2)
    else
      nrow = size(self%fm_solutions,1)
      ncol = size(self%fm_solutions,2)
      deallocate(self%fm_solutions)
      allocate(self%fm_solutions(nrow+1,ncol))
      self%fm_solutions(1:nrow,:) = other%fm_solutions
    end if
    self%fm_solutions(nrow+1,:) = fm_profile

    call other%delete

    return
  end subroutine add_inversion

  subroutine inversion_ppi(self,weight)
    ! Subroutine to invert a profile using the ppi technique
    ! - self :: inversion object to save all the results (ppi)
    ! - weight :: weighting method: 0 for equally weighted or 1 for
    !             gaussian weight (optional)
    class(ppi_type), intent(inout) :: self
    integer(iprec), optional, intent(in) :: weight

    type(svd_type) :: svd ! SVD object for inversions

    real(rprec) :: prof_type, rmse, wtd ! Dummy reals

    real(rprec), parameter :: fils = -999.9_rprec ! Filling value

    integer(iprec), dimension(:), allocatable :: w_series ! Dummy allocatable
    integer(iprec) :: n_alpha, n_w, n_sol, n_depth, n_time, t ! Dummy integers
    integer(iprec) :: i, j, k, ind ! Dummy indices
    character(20) :: ni ! Dummy character

    character(len=*), parameter :: pname = "mod_ppi < inversion_ppi"

    call init_logger

    call book%record(pname,trace,'Starting subroutine')

    n_time = size(self%model_type%time_series) ! Number of time steps
    n_alpha = size(self%model_type%alpha) ! Number of diffusivities
    n_w = self%model_type%w(2) - self%model_type%w(1) + 1_iprec ! Number of different eigenvalues
    n_sol = n_alpha * n_w * 3_iprec ! Number of solutions for each profile

    n_depth = size(self%prof%depth) ! Number of depths

    ! Build a sequence of w for inversions
    allocate(w_series(n_w))
    do i = 1, n_w
      w_series(i) = (i-1_iprec) + self%model_type%w(1)
    end do ! w

    ! Then, perform the inversions using a SVD algorithm
    ! First, the original anomaly profile
    ind = 0_iprec

    do i = 1, n_alpha
      do j = 1, n_w
        ind = ind + 1_iprec
        write(ni,'(i6)') ind
        call book%record(pname,trace,"Inversion "//ni//" of ",n_alpha*n_w*3_iprec)

        ! Perform the SVD inversion
        call svd%new(self%model_type%alpha(i),w_series(j),&
                     self%model_type%logy,self%model_type%time_series,&
                     self%prof)
        call svd%inversion()

        if(i.eq.1_iprec.and.j.eq.1_iprec) then
          self%years = svd%inv(:,1)
        end if


        ! Classify the results according to the three anomalies
        do k = 1, 3
          if(k.eq.1_iprec) then
            prof_type = -1.0_rprec
            call book%record(pname,debug,"Prof_type = ",prof_type)

            ! Check all inversions, to avoid data problems
            t = test_inversion(svd%inv(:,k+1),self%model_type%tol,fils)

            if(t.eq.0_iprec) then
              call book%record(pname,warning,"WARNING - Model rejected")
              call book%record(pname,debug,"WARNING - Model rejected"&
               ,[self%model_type%alpha(i),real(w_series(j),rprec),prof_type])
            elseif(t.eq.1_iprec) then
              ! Compare the result from the forward model with the original anomaly
              rmse = rmse_prof(svd%prof%atemp_minus,svd%fm(:,k+1))
              call book%record(pname,debug,"RMSE = ",rmse)

              ! Select the weighting scheme
              if(.not.present(weight)) then
                wtd = 1.0_rprec
              elseif(present(weight).and.weight.eq.0_iprec) then
                wtd = 1.0_rprec
              elseif(present(weight).and.weight.eq.1_iprec) then
                wtd = dexp( (-rmse**2) / (self%model_type%err**2) )
              end if
              call book%record(pname,debug,"Weight = ",wtd)

              ! Save the inversions
              call self%add_inversion(prof_type,&
                              self%model_type%alpha(i),&
                              w_series(j),wtd,svd%inv(:,k+1),svd%fm(:,k+1))
            end if
          elseif(k.eq.2_iprec) then
            prof_type = 0.0_rprec
            call book%record(pname,debug,"Prof_type = ",prof_type)

            ! Check all inversions, to avoid data problems
            t = test_inversion(svd%inv(:,k+1),self%model_type%tol,fils)

            if(t.eq.0_iprec) then
              call book%record(pname,warning,"WARNING - Model rejected")
              call book%record(pname,debug,"WARNING - Model rejected"&
               ,[self%model_type%alpha(i),real(w_series(j),rprec),prof_type])
            elseif(t.eq.1_iprec) then
              ! Compare the result from the forward model with the original anomaly
              rmse = rmse_prof(svd%prof%atemp_original,svd%fm(:,k+1))
              call book%record(pname,debug,"RMSE = ",rmse)

              ! Select the weighting scheme
              if(.not.present(weight)) then
                wtd = 1.0_rprec
              elseif(present(weight).and.weight.eq.0_iprec) then
                wtd = 1.0_rprec
              elseif(present(weight).and.weight.eq.1_iprec) then
                wtd = dexp( (-rmse**2) / (self%model_type%err**2) )
              end if
              call book%record(pname,debug,"Weight = ",wtd)

              ! Save the inversions
              call self%add_inversion(prof_type,&
                              self%model_type%alpha(i),&
                              w_series(j),wtd,svd%inv(:,k+1),svd%fm(:,k+1))
            end if
          elseif(k.eq.3_iprec) then
            prof_type = 1.0_rprec
            call book%record(pname,debug,"Prof_type = ",prof_type)

            ! Check all inversions, to avoid data problems
            t = test_inversion(svd%inv(:,k+1),self%model_type%tol,fils)

            if(t.eq.0_iprec) then
              call book%record(pname,warning,"WARNING - Model rejected")
              call book%record(pname,debug,"WARNING - Model rejected"&
               ,[self%model_type%alpha(i),real(w_series(j),rprec),prof_type])
            elseif(t.eq.1_iprec) then
              ! Compare the result from the forward model with the original anomaly
              rmse = rmse_prof(svd%prof%atemp_minus,svd%fm(:,k+1))
              call book%record(pname,debug,"RMSE = ",rmse)

              ! Select the weighting scheme
              if(.not.present(weight)) then
                wtd = 1.0_rprec
              elseif(present(weight).and.weight.eq.0_iprec) then
                wtd = 1.0_rprec
              elseif(present(weight).and.weight.eq.1_iprec) then
                wtd = dexp( (-rmse**2) / (self%model_type%err**2) )
              end if
              call book%record(pname,debug,"Weight = ",wtd)

              ! Save the inversions
              call self%add_inversion(prof_type,&
                              self%model_type%alpha(i),&
                              w_series(j),wtd,svd%inv(:,k+1),svd%fm(:,k+1))
            end if
          end if
        end do

        ! Remove SVD object
        call svd%delete()

      end do ! j
    end do ! i

    ! Create general solution
    call self%general_inversion()

    ! Create general forward model
    call self%general_forward_model()

    call book%record(pname,trace,'Ending subroutine')

    return
  end subroutine inversion_ppi

  subroutine general_inversion(self)
    ! Subroutine to estimate the 5th,50, and 95th weighted percentiles
    !            of all inversions
    ! - self :: all the inversions to be considered (inversion)
    class(ppi_type), intent(inout) :: self

    real(rprec) :: q ! Dummy real
    real(rprec), allocatable :: w(:) ! Dummy real

    integer(iprec) :: ncol, nrow ! Dummy integer
    integer(iprec) :: i ! Dummy index

    character(len=*), parameter :: pname = "mod_ppi < general_inversion"

    call init_logger

    call book%record(pname,trace,'Starting subroutine')

    ! Dimensions of the matrix of solutions
    nrow = size(self%svd_solutions,1)
    ncol = size(self%svd_solutions,2)

    if(nrow.lt.2) then
      call book%record(pname,error,"ERROR - No solutions for this profile")
      stop
    end if

    ! Allocate weights
    allocate(w(nrow))

    ! Start the solutions
    allocate(self%general_solution(ncol,3))

    ! Estimate the 5th, 50th, and 95th weighted quantiles
    do i = 1, ncol
      ! Weights
      w = self%metadata(:,4)

      if(nrow.eq.3_iprec) then
        q = minval(self%svd_solutions(:,i))
        self%general_solution(i,1) = q

        q = sum(self%svd_solutions(:,i)) / size(self%svd_solutions(:,i))
        self%general_solution(i,2) = q

        call wtd_quantile(w,self%svd_solutions(:,i),0.975_rprec,q)
        q = maxval(self%svd_solutions(:,i))
        self%general_solution(i,3) = q
      else
        ! 5th quantile
        call wtd_quantile(w,self%svd_solutions(:,i),0.025_rprec,q)
        self%general_solution(i,1) = q

        ! 50th quantile
        call wtd_quantile(w,self%svd_solutions(:,i),0.5_rprec,q)
        self%general_solution(i,2) = q

        ! 95th quantile
        call wtd_quantile(w,self%svd_solutions(:,i),0.975_rprec,q)
        self%general_solution(i,3) = q
      end if
    end do ! i

    call book%record(pname,trace,'Endding subroutine')

    return
  end subroutine general_inversion

  subroutine general_forward_model(self)
    ! Subroutine to estimate the FMs from the general solutions
    ! - self :: object with the general solutions (ppi)
    class(ppi_type), intent(inout) :: self

    type(fm_type) :: fm ! Forward model object

    real(rprec), allocatable :: time(:), inv(:)

    real(rprec) :: malpha ! Dummy reals
    integer(iprec) :: i ! Dummy integer

    allocate(self%general_fm(size(self%prof%depth), 3))

    ! Estimate mean diffusivity
    malpha = 0.0_rprec
    do i = 1, size(self%model_type%alpha)
      malpha = malpha + self%model_type%alpha(i)
    end do
    malpha = malpha / real( size(self%model_type%alpha), rprec )

    ! Minus profile
    time = self%model_type%logy-self%years
    inv = self%general_solution(:,1)
    call fm%new(malpha,time,inv,self%prof%depth)
    call fm%forward_model()
    self%general_fm(:,1) = fm%prof%atemp_original
    call fm%delete()

    ! Original profile
    inv = self%general_solution(:,2)
    call fm%new(malpha,time,inv,self%prof%depth)
    call fm%forward_model()
    self%general_fm(:,2) = fm%prof%atemp_original
    call fm%delete()

    ! Plus profile
    inv = self%general_solution(:,3)
    call fm%new(malpha,time,inv,self%prof%depth)
    call fm%forward_model()
    self%general_fm(:,3) = fm%prof%atemp_original
    call fm%delete()

    return
  end subroutine general_forward_model

  subroutine write_ppi(self,ofile,u2)
    ! Subroutine to write the results of the PPI ensemble.
    ! - self :: object to be saved (ppi)
    ! - ofile :: name of the file to save the PPI results (char)
    ! - u2 :: number indicating the unit to write the file (integer, optional)
    class(ppi_type), intent(in) :: self
    character(len=*), intent(in) :: ofile
    integer(iprec), optional, intent(in) :: u2

    ! Filling values
    real(rprec), parameter :: fils = -999.9_rprec
    integer(iprec), parameter :: ifils = -999_iprec

    ! Dummy characters
    character(90) :: f, n ! Dummy characters

    integer(iprec) :: u ! Dummy integer
    integer(iprec) :: i, j ! Dummy indices

    character(len=*), parameter :: pname = "mod_ppi < write_ppi"

    call init_logger

    if(present(u2)) then
      u = u2
    else
      u = 90_iprec
    end if

    open(u,file=trim(ofile),action="write")

    ! First, relevant parameters of the model
    f = "(a,es12.3e3)"
    if(test_real(self%model_type%tol,fils).eq.1) then
      write(u,trim(f)) '#tol', self%model_type%tol
    end if
    if(test_real(self%model_type%err,fils).eq.1) then
      write(u,trim(f)) '#err', self%model_type%err
    end if
    f = "(a,i4)"
    if(test_integer(self%model_type%w(1),ifils).eq.1) then
      write(u,trim(f)) '#wmin', self%model_type%w(1)
    end if
    if(test_integer(self%model_type%w(2),ifils).eq.1) then
      write(u,trim(f)) '#wmax', self%model_type%w(2)
    end if
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
    if(allocated(self%prof%atemp_minus).and.&
       test_array(self%prof%atemp_minus,fils).eq.1.and.&
       allocated(self%prof%atemp_plus).and.&
       test_array(self%prof%atemp_plus,fils).eq.1) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              self%prof%atemp_minus(i),&
                              self%prof%atemp_original(i),&
                              self%prof%atemp_plus(i)
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(allocated(self%prof%atemp_minus).and.&
       test_array(self%prof%atemp_minus,fils).eq.1.and.&
       .not.allocated(self%prof%atemp_plus)) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              self%prof%atemp_minus(i),&
                              self%prof%atemp_original(i),&
                              fils
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(.not.allocated(self%prof%atemp_minus).and.&
       allocated(self%prof%atemp_plus).and.&
       test_array(self%prof%atemp_plus,fils).eq.1) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              fils,&
                              self%prof%atemp_original(i),&
                              self%prof%atemp_plus(i)
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(.not.allocated(self%prof%atemp_minus).and.&
       .not.allocated(self%prof%atemp_plus)) then
      do i = 1, size(self%prof%depth)
        write(u,trim(f)) '#a', self%prof%depth(i),&
                              fils,&
                              self%prof%atemp_original(i),&
                              fils
      end do

      write(u,*) ""
      write(u,*) ""
    elseif(.not.allocated(self%prof%atemp_original).or.&
      test_array(self%prof%atemp_original,fils).eq.0) then
      call book%record(pname,error,"ERROR - There is no anomaly to write")
      stop
    end if


    ! Write the key for the solutions
    f = "(a,es12.3e3,es12.3e3,es12.3e3,3x,a5,a5)"
    do i = 1, size(self%metadata, 1)
      write(n,'(i4)') i
      n = trim(adjustl(n))
      if(self%metadata(i,1).eq.-1.0_rprec) then
        write(u,trim(f)) '#k',self%metadata(i,2),&
                               self%metadata(i,3),&
                               self%metadata(i,4),n,'m'
      elseif(self%metadata(i,1).eq.0.0_rprec) then
        write(u,trim(f)) '#k',self%metadata(i,2),&
                               self%metadata(i,3),&
                               self%metadata(i,4),n,'o'
      elseif(self%metadata(i,1).eq.1.0_rprec) then
        write(u,trim(f)) '#k',self%metadata(i,2),&
                               self%metadata(i,3),&
                               self%metadata(i,4),n,'p'
      end if
    end do

    write(u,*) ""
    write(u,*) ""

   ! ! Write all the solutions
   ! tl = size(inversion%svd_solutions, 1)
   ! write(n,'(i6)') tl + 1_iprec
   ! f = "(a,"//trim(n)//"es12.3e3)"
   ! call inversion2step(inversion%model_type%time_series,&
   !                     inversion%svd_solutions,step)
   ! do i = 1, 2*size(inversion%svd_solutions, 2)
   !   write(u,trim(f)) '#s', (step(i,j),j=1,tl+1)
   ! end do

   ! write(u,*) ""
   ! write(u,*) ""

   ! ! Write all the forward model profiles
   ! tl = size(inversion%fm_solutions, 1)
   ! write(n,'(i6)') tl + 1_iprec
   ! f = "(a,"//trim(n)//"es12.3e3)"
   ! do i = 1, size(inversion%fm_solutions, 2)
   !   write(u,trim(f)) '#p', prof%depth(i),&
   !                           (inversion%fm_solutions(j,i),j=1,tl)
   ! end do

   ! write(u,*) ""
   ! write(u,*) ""

    ! Write general solution
    f = "(a,4es12.3e3)"
    do i = 1, size(self%years)
      write(u,trim(f)) '#g', self%years(i),self%general_solution(i,1),&
                             self%general_solution(i,2),&
                             self%general_solution(i,3)
    end do

    write(u,*)
    write(u,*)

    ! Write general forward model
    f = "(a,4es12.3e3)"
    do i = 1, size(self%prof%depth)
      write(u,trim(f)) '#f',self%prof%depth(i),(self%general_fm(i,j),j=1,3)
    end do

    close(u)

    return
  end subroutine write_ppi

  !-------------
  ! Other methods
  !-------------

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

    character(len=*), parameter :: pname = "mod_ppi < test_array"

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

    character(len=*), parameter :: pname = "mod_ppi < test_real"

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

    character(len=*), parameter :: pname = "mod_ppi < test_integer"

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

  function test_inversion(inversion,tol,fils) result(t)
    ! Function to check that the individual inversions are performed
    !          without strange results.
    ! - inversion => the inversion to be checked (real(:))
    ! - tol => maxmum temperature difference allowed between time
    !          steps (real)
    ! - t => 1 if all test are successful, 0 otherwise
    real(rprec), intent(in) :: inversion(:)
    real(rprec), intent(in) :: tol, fils
    integer(iprec) :: t

    character(len=7) :: c1,c2,c3,c4,c5! Dummy characters for logger

    integer(iprec) :: i, ind ! Dummy indice

    character(len=*), parameter :: pname = "mod_ppi < test_inversion"

    call init_logger

    call book%record(pname,trace,"Starting function")

    t = 1_iprec ! Initialize test

    ! Test if all values are zero
    ind = 0_iprec
    do i = 1, size(inversion)
      if(inversion(i).eq.0.0_rprec) ind = ind + 1_iprec
    end do ! i
    if(ind.eq.size(inversion)) then
      t = 0_iprec
      c1 = "z=0"
    else
      c1 = "z=1"
    end if
    call book%record(pname,trace,"Zero test -> "//trim(c1), ind)

    ! Check if there are NANs in the inversion
    ind = 0_iprec
    do i = 1, size(inversion)
      if(inversion(i).ne.inversion(i)) ind = ind + 1_iprec
    end do ! i
    if(ind.ne.0_iprec) then
      t = 0_iprec
      c2 = "NA=0"
    else
      c2 = "NA=1"
    end if
    call book%record(pname,trace,"NANs test -> "//trim(c2), ind)

    ! Check there are no infinite values
    ind = 0_iprec
    do i = 1, size(inversion)
      if(inversion(i).ge.huge(rprec)) ind = ind + 1_iprec
    end do !i
    if(ind.ne.0_iprec) then
      t = 0_iprec
      c3 = "infty=0"
    else
      c3 = "infty=1"
    end if
    call book%record(pname,trace,"Infty test -> "//trim(c3), ind)

    ! Check all differences are within the tolerance level
    do i = 1, size(inversion)-1
      if(abs(inversion(i)-inversion(i+1)).gt.tol) then
        t = 0_iprec
        c4 = "tol=0"
        exit
      else
        if(i.eq.size(inversion)-1) c4 = "tol=1"
      end if
    end do
    call book%record(pname,trace,"Tol test -> "//trim(c4))

    ! Test if there is some filling value
    ind = 0_iprec
    do i = 1, size(inversion)
      if(inversion(i).eq.fils) ind = ind + 1_iprec
    end do ! i
    if(ind.ne.0_iprec) then
      t = 0_iprec
      c5 = "z=0"
    else
      c5 = "z=1"
    end if
    call book%record(pname,trace,"Filling value test -> "//trim(c5), ind)


    call book%record(pname,trace,"Ending function")

    return
  end function test_inversion

  function rmse_prof(prof,fm_prof) result(r)
    ! Function to estimate the RMSE between two profiles
    ! - prof = reference profile (real). Typically a measured anomaly
    !        profile.
    ! - fm_prof = another profile (real). Typically from forward model.
    real(rprec), intent(in) :: prof(:), fm_prof(:)
    real(rprec) :: r

    real(rprec) :: s ! Dummy real

    integer(iprec) :: n ! Dummy integer
    integer(iprec) :: i ! Dummy index

    character(len=*), parameter :: pname = "mod_ppi < rmse_prof"

    call init_logger

    n = size(prof)

    if(size(prof).ne.size(fm_prof)) then
      call book%record(pname,error,'ERROR - Different sizes in profiles')
      stop
    end if

    s = 0.0_rprec
    do i = 1, n
      s = s + ( prof(i) - fm_prof(i) )**2
    end do ! i
    r = dsqrt(s/dble(n))

    return
  end function rmse_prof

  subroutine wtd_quantile(w,x,alpha,q,v)
    ! Subroutine to estimate weighted quantiles
    ! - w = the weights associated to the elements of x (real)
    ! - x = elements to be processed (real)
    ! - alpha = probability for the corresponding quantiles (real)
    ! - q = quantile (result)
    ! - v = type of warning (optional, integer)
    real(rprec), intent(in) :: alpha
    real(rprec), dimension(:), intent(in) :: x
    real(rprec), dimension(:), intent(inout) :: w

    real(rprec), intent(out) :: q ! Quantile
    integer(iprec), optional, intent(out) :: v ! Verbose error

    real(rprec), dimension(:), allocatable :: p

    real(rprec) :: gam, s
    integer(iprec), dimension(1) :: ei
    integer(iprec) :: n, i, k0, ind, AllocateStatus

    character(len=*), parameter :: pname = "mod_ppi < wtd_quantile"

    call init_logger

    if(present(v)) v = 0_iprec ! SO far, no errors/warnings

    ! Data length
    n = size( x )

    if(n.le.1_iprec) then
      call book%record(pname,error,"ERROR - you need more than 1 datum")
      stop
    end if

    ! Check the weights
    ind = 0_iprec
    do i = 1, size(w)
      if(w(i).eq.0.0_rprec) ind = ind + 1_iprec
    end do
    if(ind .eq. size(w)) then
      call book%record(pname,warning,"WARNING - all weights are zero!")
      w = 1.0_iprec
    end if

    ! Allocate vector of prbabilities
    allocate(p(n), stat = AllocateStatus)
    if (AllocateStatus /= 0) stop "*** Not enough memory ***"

    ! Sort the data
    call sort2(n,x,w)

    ! Vector of percentiles
    s = 0.0_rprec
    do i = 1, n
      s = s + w(i)
      p(i) = ( ( s * real(n,rprec) / sum(w) ) - 1.0_rprec/3.0_rprec ) /&
             ( real(n,rprec) + 1.0_rprec/3.0_rprec )
    end do ! i

    ! Find k0
    k0 = 0_iprec
    do i = 1, (n-1)
      if(p(i).le.alpha.and.p(i+1).gt.alpha) then
        k0 = i
        exit
      end if
    end do !  i
    if(k0.eq.0_iprec) then
      call book%record(pname,warning,'WARNING - cannot find a k0 index.&
                                & q defined as closest percentile')
      if(present(v)) v = 1_iprec ! Inform the main program that there was a warming
      if(alpha.le.p(1)) k0 = 1_iprec
      if(alpha.ge.p(n)) k0 = n - 1_iprec ! To find gam (see the k0+1 bit below)
    end if

    ! Find the x
    if(p(k0+1)-p(k0).ne.0.0_rprec) then
      ! Estimate gamma
      gam = ( alpha - p(k0) ) / ( p(k0+1) - p(k0) )

      ! Quantile
      q = ( 1.0_rprec - gam ) * x(k0) + gam * x(k0)
    else
      call book%record(pname,warning,'WARNING - there is a dominant&
         & percentile in the ensemble. q defined as x with the maximum w')
      if(present(v)) v = 2_iprec ! There is only one valid value for percentile!
      ei = maxloc(w)
      q = x(ei(1))
    end if

    if(test_real(q).eq.0_iprec) then
      call book%record(pname,error,"ERROR - something bad in quantile")
      call book%record(pname,error,"q = ",q)
      call book%record(pname,error,"alpha = ",alpha)
      call book%record(pname,error,"k0 = ",k0)
      call book%record(pname,error,"gam = ",gam)
      call book%record(pname,error,"p = ",p)
      call book%record(pname,error,"x = ",x)
      call book%record(pname,error,"w = ",w)
      stop
    end if

    return
  end subroutine wtd_quantile

  subroutine sort2(n,arr,brr)
    ! Subroutine to sort a set of elements from minor to major, and
    !            change the same elements in a companion set.
    ! - n = dimension of the sets (integer)
    ! - arr = array to be sorted (real)
    ! - brr = array that will suffer the same element alteration
    !            as arr (real)
    integer(iprec) :: n
    real(rprec) :: arr(n),brr(n)
    integer(iprec), parameter :: M=7,NSTACK=50
    integer(iprec) i,ir,j,jstack,k,l
    integer(iprec), dimension(NSTACK) :: istack
    real(rprec) :: a,b,temp
    jstack=0
    l=1
    ir=n
1   if(ir-l.lt.M)then
      do 12 j=l+1,ir
        a=arr(j)
        b=brr(j)
        do 11 i=j-1,l,-1
          if(arr(i).le.a)goto 2
          arr(i+1)=arr(i)
          brr(i+1)=brr(i)
11      continue
        i=l-1
2       arr(i+1)=a
        brr(i+1)=b
12    continue
      if(jstack.eq.0)return
      ir=istack(jstack)
      l=istack(jstack-1)
      jstack=jstack-2
    else
      k=(l+ir)/2
      temp=arr(k)
      arr(k)=arr(l+1)
      arr(l+1)=temp
      temp=brr(k)
      brr(k)=brr(l+1)
      brr(l+1)=temp
      if(arr(l).gt.arr(ir))then
        temp=arr(l)
        arr(l)=arr(ir)
        arr(ir)=temp
     temp=brr(l)
     brr(l)=brr(ir)
     brr(ir)=temp
    endif
    if(arr(l+1).gt.arr(ir))then
      temp=arr(l+1)
      arr(l+1)=arr(ir)
      arr(ir)=temp
      temp=brr(l+1)
      brr(l+1)=brr(ir)
      brr(ir)=temp
    endif
    if(arr(l).gt.arr(l+1))then
      temp=arr(l)
      arr(l)=arr(l+1)
      arr(l+1)=temp
      temp=brr(l)
      brr(l)=brr(l+1)
      brr(l+1)=temp
    endif
    i=l+1
    j=ir
    a=arr(l+1)
    b=brr(l+1)
3   continue
      i=i+1
    if(arr(i).lt.a)goto 3
4   continue
      j=j-1
    if(arr(j).gt.a)goto 4
    if(j.lt.i)goto 5
    temp=arr(i)
    arr(i)=arr(j)
    arr(j)=temp
    temp=brr(i)
    brr(i)=brr(j)
    brr(j)=temp
    goto 3
5   arr(l+1)=arr(j)
    arr(j)=a
    brr(l+1)=brr(j)
    brr(j)=b
    jstack=jstack+2
    if(jstack.gt.NSTACK) stop 'NSTACK too small in sort2'
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
  end subroutine sort2


end module ppi_module




