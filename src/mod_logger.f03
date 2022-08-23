!
!----------------------------------------------------------------------
! Program to record a log of operations and errors.
!
! This program needs:
! - 01b_mod_kinds.f03
!
! Francisco Jose Cuesta-Valero
! 2021-01-01 (Leipzig)
!----------------------------------------------------------------------

module logger_module
  use kinds_module
  implicit none

  private

  save

  type, public :: logger_type
    integer(iprec) :: level ! Level of the messages to be recorded
    integer(iprec) :: unt
    character(128) :: logbook
  contains
    procedure, private :: new_private, new_copy_private
    generic :: new => new_private, new_copy_private
    procedure :: delete => delete_private
    !procedure :: assignment(=) => assign_private
    procedure, private :: rec_real, rec_integer, rec_real_array, rec_integer_array
    procedure, private :: rec_text_only, rec_string
    generic :: record => rec_real, rec_integer, rec_real_array,&
                         rec_integer_array, rec_text_only, rec_string
  end type

  !interface assignment(=)
  !  module procedure assign_private
  !end interface

  ! Level Definitions
  integer(iprec), public, parameter :: debug = 1_iprec ! Debug logging level
  integer(iprec), public, parameter :: trace = 2_iprec ! Trace logging level
  integer(iprec), public, parameter :: warning = 3_iprec ! Warning logging level
  integer(iprec), public, parameter :: error = 4_iprec ! Error logging level

contains

  ! ------------------------
  ! Standard ADT Methods. Construction, Destruction, Copying, and Assignment.
  ! ------------------------
  subroutine new_private(self,level)
    ! Subroutine to initialize a logger.
    ! - self :: the object to be initialized (logger_type)
    ! - level :: logging level (integer)
    class (logger_type), intent(out) :: self
    integer(iprec), intent(in) :: level

    self%level = level

    return
  end subroutine new_private

  subroutine new_copy_private(self, other)
    ! Subroutine to create a logger from another one
    ! - self :: object to be created (logger_type)
    ! - other :: object to be copied (logger_type)
    class (logger_type), intent(out) :: self
    class (logger_type), intent(in) :: other

    self%level = other%level

    return
  end subroutine new_copy_private

  subroutine delete_private(self)
    ! Subroutine to delete a logger
    ! - self :: logger to be deleted (logger_type)
    class (logger_type), intent(inout)  :: self

    close(self%unt)

    return
  end subroutine delete_private

 ! subroutine assign_private(self, other)
 !   ! Subtoutine to allow = sign to work with loggers
 !   ! - self :: logger to receive the assignement (logger_type)
 !   ! - other :: logger to be assigned (logger_type)
 !   class (logger_type), intent(inout) :: self
 !   class (logger_type), intent(in) :: other
 !
 !   self%level = other%level
 !
 !   return
 ! end subroutine assign_private

  ! ------------------------
  ! Accessors.
  ! ------------------------
  subroutine set_level(self, level)
    ! Subroutine to set the logging level in a logger
    ! - self :: logger which level will be set (logger_type)
    ! - level :: level to be set (integer)
    class (logger_type), intent(inout) :: self
    integer(iprec), intent(in) :: level

    self%level = level

    return
  end subroutine set_level

  function get_level(self)
    ! Function to check the level of a logger
    ! - self :: logger to retrive the level (logger_type)
    ! - get_level :: level to be retrieved (integer)
    class (logger_type), intent(in) :: self
    integer(iprec) :: get_level

    get_level = self%level

    return
  end function get_level


  ! ------------------------
  ! Other methods.
  ! ------------------------
  subroutine rec_text_only(self, pname, level, text)
    ! Subroutine to record text only
    ! - self :: logger that is going to record the message (logger_type)
    ! - pname :: name of the program (char)
    ! - level :: logging level of the message (integer)
    ! - text :: message to be recorded (char)
    class (logger_type), intent(in) :: self
    character(len=*), intent(in) :: pname
    integer(iprec), intent(in) :: level
    character(len=*), intent(in) :: text

    if( level >= self%level ) then
      print *, trim(pname)//': '//trim(text)
    end if

    return
  end subroutine rec_text_only

  subroutine rec_string(self, pname, level, key, string)
    ! Subroutine to record text only
    ! - self :: logger that is going to record the message (logger_type)
    ! - pname :: name of the program (char)
    ! - level :: logging level of the message (integer)
    ! - key :: message to be recorded (char)
    ! - string :: string value to be recorded (char)
    class (logger_type), intent(in) :: self
    character(len=*), intent(in) :: pname
    integer(iprec), intent(in) :: level
    character(len=*), intent(in) :: key
    character(len=*), intent(in) :: string

    if( level >= self%level ) then
      print *, trim(pname)//': '//trim(key), string
    end if

    return
  end subroutine rec_string

  subroutine rec_integer(self, pname, level, key, val)
    ! Subroutine to record integers
    ! - self :: logger that is going to record the message (logger_type)
    ! - pname :: name of the program (char)
    ! - level :: logging level of the message (integer)
    ! - key :: message to be recorded (char)
    ! - val :: integer value to be recorded (integer)
    class (logger_type), intent(in) :: self
    integer(iprec), intent(in) :: level
    character(len=*), intent(in) :: key
    integer(iprec), intent(in) :: val
    character(len=*), intent(in) :: pname

    if( level >= self%level ) then
      print *, trim(pname)//': '//trim(key), val
    end if

    return
  end subroutine rec_integer

  subroutine rec_integer_array(self, pname, level, key, vals)
    ! Subroutine to record integer arrays
    ! - self :: logger that is going to record the message (logger_type)
    ! - pname :: name of the program (char)
    ! - level :: logging level of the message (integer)
    ! - key :: message to be recorded (char)
    ! - vals :: integer array to be recorded (integer)
    class (logger_type), intent(in) :: self
    integer(iprec), intent(in) :: level
    character(len=*), intent(in) :: key
    integer(iprec), dimension(:), intent(in) :: vals
    character(len=*), intent(in) :: pname

    if ( level >= self%level ) then
      print *, trim(pname)//': '//trim(key)
      print *, vals
    endif

    return
  end subroutine rec_integer_array

  subroutine rec_real(self, pname, level, key, val)
    ! Subroutine to record reals
    ! - self :: logger that is going to record the message (logger_type)
    ! - pname :: name of the program (char)
    ! - level :: logging level of the message (integer)
    ! - key :: message to be recorded (char)
    ! - val :: real to be recorded (real)
    class (logger_type), intent(in) :: self
    integer(iprec), intent(in) :: level
    character(len=*), intent(in) :: key
    real(rprec), intent(in) :: val
    character(len=*), intent(in) :: pname

    if ( level >= self%level ) then
      print *, trim(pname)//': '//trim(key), val
    end if

    return
  end subroutine rec_real

  subroutine rec_real_array(self, pname, level, key, vals)
    ! Subroutine to record real arrays
    ! - self :: logger that is going to record the message (logger_type)
    ! - pname :: name of the program (char)
    ! - level :: logging level of the message (integer)
    ! - key :: message to be recorded (char)
    ! - vals :: real array to be recorded (real)
    class (logger_type), intent(in) :: self
    integer(iprec), intent(in) :: level
    character(len=*), intent(in) :: key
    real(rprec), dimension(:), intent(in) :: vals
    character(len=*), intent(in) :: pname

    if ( level >= self%level ) then
      print *, trim(pname)//': '//trim(key)
      print *, vals
    endif

    return
  end subroutine rec_real_array

end module logger_module
