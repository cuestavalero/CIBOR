!
!----------------------------------------------------------------------
! Module to work with ensembles of data, particularly with ensembles of
! PPI inversions.
!
! This module needs:
! kinds_module
! logger_module
!
! Francisco Jose Cuesta-Valero
! 2021-02-15 (Leipzig)
!----------------------------------------------------------------------

module ensembles_module
  use kinds_module
  use logger_module
  implicit none

  private

  type, public :: ensemble_type
    real(rprec), public, allocatable :: ensemble(:,:) ! The ensemble
    real(rprec), public :: fils ! Missing value
  contains
    procedure :: new_fils_ensemble, new_copy_ensemble
    generic, public :: new => new_fils_ensemble, new_copy_ensemble
    procedure, public :: delete => delete_ensemble
    procedure, public :: add => add_member
    procedure, public :: mean_rows, mean_cols
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

!  subroutine new_ensemble(self,n)
!    ! Subroutine to create a new ensemble object
!    ! - self :: ensemble to be initialized (ensemble)
!    ! - n :: length of each ensemble member (integer)
!    class(ensemble_type), intent(out) :: self
!    integer(iprec), intent(in) :: n
!
!    allocate(self%ensemble(1,n))
!
!    return
!  end subroutine new_ensemble

  subroutine new_fils_ensemble(self,fils)
    ! Subroutine to create a new ensemble object
    ! - self :: ensemble to be initialized (ensemble)
    ! - fils :: missing value (real)
    class(ensemble_type), intent(out) :: self
    real(rprec), intent(in) :: fils

    self%fils = fils

    return
  end subroutine new_fils_ensemble

  subroutine new_copy_ensemble(self,other)
    ! Subroutine to create an ensemble by copying another ensemble
    ! - self :: new ensemble
    ! - other :: ensemble to be copied
    class(ensemble_type), intent(out) :: self
    class(ensemble_type), intent(in) :: other

    integer(iprec) :: ncol, nrow ! Dummy integers
    integer(iprec) :: i ! Dummy indices

    character(len=*), parameter :: pname="mod_ensembles < copy_ensemble"

    nrow = size(other%ensemble,1)
    ncol = size(other%ensemble,2)

    self%fils = other%fils
    allocate(self%ensemble(nrow,ncol))
    self%ensemble = other%ensemble

    return
  end subroutine new_copy_ensemble

  subroutine delete_ensemble(self)
    ! subroutine to delete an ensemble
    ! - self :: ensemble to be deleted (ensemble)
    class(ensemble_type), intent(inout) :: self

    if(allocated(self%ensemble)) deallocate(self%ensemble)

    return
  end subroutine delete_ensemble

  subroutine add_member(self,member)
    ! Subroutine to add a new member to the ensemble
    ! - self :: ensemble (ensemble)
    ! - member :: new member to be added to the ensemble (real array)
    class(ensemble_type), intent(inout) :: self
    real(rprec), intent(in) :: member(:)

    type(ensemble_type) :: other
    integer(iprec) :: ncol, nrow, n ! Dummy integers

    character(len=*), parameter :: pname="mod_ensembles < add_member"

    call init_logger

    if(allocated(self%ensemble)) then

      if(size(member).ne.size(self%ensemble,2)) then
        call book%record(pname,error,"ERROR - size of the member different&
                                     & to the length of the ensemble")
        stop
      end if

      nrow = size(self%ensemble,1)
      ncol = size(self%ensemble,2)

      ! Save the ensemble so far
      call other%new(self)
      deallocate(self%ensemble)
      allocate(self%ensemble(nrow+1,ncol))
      self%ensemble(1:nrow,:) = other%ensemble
    else
      n = size(member)
      allocate(self%ensemble(1,n))
      nrow = 0_iprec
    end if
    self%ensemble(nrow+1,:) = member

  end subroutine add_member

  subroutine mean_rows(self,means,nmem,coeff)
    ! Subroutine to estimate the average of each row of an ensemble
    ! - self :: ensemble (ensemble)
    ! - means :: array with the averages of each row
    ! - nmem :: employed in analysis (optional, real)
    ! - coeff :: coefficients to weight the mean (optional, real)
    class(ensemble_type), intent(in) :: self
    real(rprec), allocatable, intent(out) :: means(:)
    real(rprec), optional, allocatable, intent(out) :: nmem(:)
    real(rprec), optional, allocatable, intent(in) :: coeff(:)

    real(rprec) :: s, c ! Dummy real
    integer(iprec) :: nrow, ncol, n ! Dummy integers
    integer(iprec) :: i,j ! Dummy indices

    character(len=*), parameter :: pname="mod_ensembles < mean_rows"

    call init_logger

    nrow = size(self%ensemble,1)
    ncol = size(self%ensemble,2)

    if(present(nmem)) then
      if(.not.allocated(nmem)) then
        allocate(nmem(nrow))
      end if
    end if

    if(allocated(means)) then
      if(size(means).ne.nrow) then
        call book%record(pname,error, "ERROR - size of ensemble is not&
                                       & equal to size of means")
        stop
      end if
    else
      allocate(means(nrow))
    end if

    if(present(coeff)) then
      if(size(coeff).ne.ncol) then
        call book%record(pname,error, "ERROR - size of coeff is not equal to&
                                       & the number of cols")
        stop
      end if

      do i = 1, nrow
        c = 0.0_rprec
        s = 0.0_rprec
        n = 0_iprec
        do j = 1, ncol
          if(self%ensemble(i,j).ne.self%fils) then
            s = s + coeff(j) * self%ensemble(i,j)
            c = c + coeff(j)
            n = n + 1_iprec
          end if
        end do
        if(n.eq.0_iprec) then
          means(i) = self%fils
          if(present(nmem)) nmem(i) = 0.0_rprec
        else
          means(i) = s / c
          if(present(nmem)) nmem(i) = real(n,rprec)
        end if
      end do

    else
      do i = 1, nrow
        s = 0.0_rprec
        n = 0_iprec
        do j = 1, ncol
          if(self%ensemble(i,j).ne.self%fils) then
            s = s + self%ensemble(i,j)
            n = n + 1_iprec
          end if
        end do
        if(n.eq.0_iprec) then
          means(i) = self%fils
          if(present(nmem)) nmem(i) = 0.0_rprec
        else
          means(i) = s / real(n,rprec)
          if(present(nmem)) nmem(i) = real(n,rprec)
        end if
      end do
    end if

    return
  end subroutine mean_rows

  subroutine mean_cols(self,means,nmem,coeff)
    ! Subroutine to estimate the average of each column of an ensemble
    ! - self :: ensemble (ensemble)
    ! - means :: array with the averages of each column
    ! - nmem :: employed in analysis (optional, real)
    ! - coeff :: coefficients to weight the mean (optional, real)
    class(ensemble_type), intent(in) :: self
    real(rprec), allocatable, intent(out) :: means(:)
    real(rprec), optional, allocatable, intent(out) :: nmem(:)
    real(rprec), optional, allocatable, intent(in) :: coeff(:)

    real(rprec) :: s, c ! Dummy real
    integer(iprec) :: nrow, ncol, n ! Dummy integers
    integer(iprec) :: i,j ! Dummy indices

    character(len=*), parameter :: pname="mod_ensembles < mean_cols"

    call init_logger

    nrow = size(self%ensemble,1)
    ncol = size(self%ensemble,2)

    if(present(nmem)) then
      if(.not.allocated(nmem)) then
        allocate(nmem(ncol))
      end if
    end if

    if(allocated(means)) then
      if(size(means).ne.ncol) then
        call book%record(pname,error, "ERROR - size of ensemble is not&
                                       & equal to size of means")
        stop
      end if
    else
      allocate(means(ncol))
    end if

    if(present(coeff)) then
      if(size(coeff).ne.nrow) then
        call book%record(pname,error, "ERROR - size of coeff is not equal to&
                                       & the number of row")
        stop
      end if

      do j = 1, ncol
        s = 0.0_rprec
        c = 0.0_rprec
        n = 0_iprec
        do i = 1, nrow
          if(self%ensemble(i,j).ne.self%fils) then
            s = s + coeff(i) * self%ensemble(i,j)
            c = c + coeff(i)
            n = n + 1_iprec
          end if
        end do
        if(n.eq.0_iprec) then
          means(j) = self%fils
          if(present(nmem)) nmem(j) = 0.0_rprec
        else
          means(j) = s / c
          if(present(nmem)) nmem(j) = real(n,rprec)
        end if
      end do

    else
      do j = 1, ncol
        s = 0.0_rprec
        n = 0_iprec
        do i = 1, nrow
          if(self%ensemble(i,j).ne.self%fils) then
            s = s + self%ensemble(i,j)
            n = n + 1_iprec
          end if
        end do
        if(n.eq.0_iprec) then
          means(j) = self%fils
          if(present(nmem)) nmem(j) = 0.0_rprec
        else
          means(j) = s / real(n,rprec)
          if(present(nmem)) nmem(j) = real(n,rprec)
        end if
      end do
    end if

    return
  end subroutine mean_cols

end module ensembles_module



