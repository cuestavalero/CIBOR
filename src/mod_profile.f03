!
!----------------------------------------------------------------------
! Module containing the object profile
!
! This module needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_regression.f03
!
! Francisco Jose Cuesta-Valero
! 2021-01-03 (Leipzig)
!----------------------------------------------------------------------

module profile_module
  use kinds_module
  use logger_module
  use regression_module
  implicit none

  private

  type :: eq_profile_type
    real(rprec) :: inter ! [C]
    real(rprec) :: s_inter ! [C]
    real(rprec) :: slope ! [C/m]
    real(rprec) :: s_slope ! [C/m]
    real(rprec) :: pvalue ! [1]
    real(rprec) :: r2 ! [1]
    integer(iprec) :: n! [1]
  end type

  type, public :: profile_type
    real(rprec) :: zmin ! Minumim depth for estimating anomalies [m]
    real(rprec) :: zmax ! Maximum depth for estimating anomalies [m]
    real(rprec), allocatable :: depth(:) ! Depths of the profile [m]
    real(rprec), allocatable :: temp(:) ! Temperature profile [C]
    real(rprec), allocatable :: atemp_plus(:) ! Plus tempe anomaly [C]
    real(rprec), allocatable :: atemp_original(:) ! Original temp anoma [C]
    real(rprec), allocatable :: atemp_minus(:) ! Minus temp anomaly [C]
    type(eq_profile_type) :: eq_profile ! The equilibrium profile
  contains
    procedure :: new_profile, new_copy_profile, new_profile_depth
    generic :: new => new_profile, new_copy_profile, new_profile_depth
    procedure :: delete
    procedure :: set_zmin, get_zmin
    procedure :: set_zmax, get_zmax
    procedure :: set_depth, set_temp, set_three_anoma, set_one_anoma
    generic :: set_anoma => set_three_anoma, set_one_anoma
    procedure :: set_slope, set_s_slope, set_inter, set_s_inter
    procedure :: set_pvalue, set_r2, set_n
    procedure :: anomaly
    procedure :: point_h, point_b
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

  subroutine new_profile(self,depth,temp)
    ! Subroutine to create a new temperature profile
    ! - self :: profile to be created (profile_type)
    ! - depth :: depths of the profile (real)
    ! - temp :: temperatures of the profile (real)
    class(profile_type), intent(out) :: self
    real(rprec), intent(in) :: depth(:)
    real(rprec), intent(in) :: temp(:)

    integer(iprec) :: n ! Dummy integer

    n = size(temp)

    allocate(self%temp(n))
    allocate(self%depth(n))
    allocate(self%atemp_plus(n))
    allocate(self%atemp_original(n))
    allocate(self%atemp_minus(n))

    self%depth = depth
    self%temp = temp

    return
  end subroutine new_profile

  subroutine new_profile_depth(self,depth)
    ! Subroutine to create a new temperature profile
    ! - self :: profile to be created (profile_type)
    ! - depth :: depths of the profile (real)
    class(profile_type), intent(out) :: self
    real(rprec), intent(in) :: depth(:)

    integer(iprec) :: n ! Dummy integer

    n = size(depth)

    allocate(self%temp(n))
    allocate(self%depth(n))
    allocate(self%atemp_plus(n))
    allocate(self%atemp_original(n))
    allocate(self%atemp_minus(n))

    self%depth = depth

    return
  end subroutine new_profile_depth

  subroutine new_copy_profile(self,other)
    ! Subroutine to create a new temperature profile from another one
    ! - self :: profile to be created (profile_type)
    ! - other :: profile to be copied (profile_type)
    class(profile_type), intent(out) :: self
    class(profile_type), intent(in) :: other

    integer(iprec) :: n ! Dummy integer

    n = size(other%temp)

    allocate(self%temp(n))
    allocate(self%depth(n))
    allocate(self%atemp_plus(n))
    allocate(self%atemp_original(n))
    allocate(self%atemp_minus(n))

    self%depth = other%depth
    self%temp = other%temp

    return
  end subroutine new_copy_profile

  subroutine delete (self)
    ! Subroutine to delete profiles
    ! - self :: the profile to be deleted (profile_type)
    class(profile_type), intent(inout) :: self

    if(allocated(self%depth)) deallocate(self%depth)
    if(allocated(self%temp)) deallocate(self%temp)
    if(allocated(self%atemp_plus)) deallocate(self%atemp_plus)
    if(allocated(self%atemp_original)) deallocate(self%atemp_original)
    if(allocated(self%atemp_minus)) deallocate(self%atemp_minus)

    return
  end subroutine delete

  ! ------------------------
  ! Accessors.
  ! ------------------------

  subroutine set_zmin(self,zmin)
    ! Subroutine to set minimum depth to estimate anomalies
    ! - self :: profile to be given a lambda value (profile)
    ! - zmin :: minimum depth to estimate anomalies (real)
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: zmin

    self%zmin = zmin

    return
  end subroutine set_zmin

  function get_zmin(self)
    ! Function to get minimim depth to estimate anomalies
    ! - self :: profile from which to fetch minimum depth (profile)
    class(profile_type), intent(in) :: self
    real(rprec) :: get_zmin

    get_zmin = self%zmin

    return
  end function get_zmin

  subroutine set_zmax(self,zmax)
    ! Subroutine to set maximum depth to estimate anomalies
    ! - self :: profile to be given a lambda value (profile)
    ! - zmax :: maximum depth to estimate anomalies (real)
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: zmax

    self%zmax = zmax

    return
  end subroutine set_zmax

  function get_zmax(self)
    ! Function to get maximum depth to estimate anomalies
    ! - self :: profile from which to fetch minimum depth (profile)
    class(profile_type), intent(in) :: self
    real(rprec) :: get_zmax

    get_zmax = self%zmax

    return
  end function get_zmax

  subroutine set_depth(self,depth)
    ! Function to set the depth of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: depth(:)
    character(len=*), parameter :: pname = "mod_profile < set_depth"

    ! Initialize logger
    call init_logger()

    if(size(depth).eq.size(self%depth)) then
      self%depth(:) = depth(:)
    else
      call book%record(pname,error,'ERROR - Different number of depths')
      stop
    end if

    return
  end subroutine set_depth

  subroutine set_temp(self,temp)
    ! Function to set the temp of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: temp(:)
    character(len=*), parameter :: pname = "mod_profile < set_temp"

    ! Initialize logger
    call init_logger()

    if(size(temp).eq.size(self%temp)) then
      self%temp(:) = temp(:)
    else
      call book%record(pname,error,'ERROR - Different number of temps')
      stop
    end if

    return
  end subroutine set_temp

  subroutine set_three_anoma(self,anoma)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: anoma(:,:)
    character(len=*), parameter :: pname = "mod_profile < set_three_anoma"

    ! Initialize logger
    call init_logger()

    if(size(anoma,1).eq.size(self%atemp_original)) then
      self%atemp_minus(:) = anoma(:,1)
      self%atemp_original(:) = anoma(:,2)
      self%atemp_plus(:) = anoma(:,3)
    else
      call book%record(pname,error,'ERROR - Different number of anomas')
      stop
    end if

    return
  end subroutine set_three_anoma

  subroutine set_one_anoma(self,anoma)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: anoma(:)

    character(len=*), parameter :: pname = "mod_profile < set_one_anoma"

    ! Initialize logger
    call init_logger()

    if(size(anoma).eq.size(self%atemp_original)) then
      self%atemp_original(:) = anoma(:)
      if(allocated(self%atemp_minus)) deallocate(self%atemp_minus)
      if(allocated(self%atemp_plus)) deallocate(self%atemp_plus)
    else
      call book%record(pname,error,'ERROR - Different number of anomas')
      stop
    end if

    return
  end subroutine set_one_anoma

  subroutine set_inter(self,inter)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: inter
    character(len=*), parameter :: pname = "mod_profile < set_inter"

    self%eq_profile%inter = inter

    return
  end subroutine set_inter

  subroutine set_s_inter(self,s_inter)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: s_inter
    character(len=*), parameter :: pname = "mod_profile < set_s_inter"

    self%eq_profile%s_inter = s_inter

    return
  end subroutine set_s_inter

  subroutine set_slope(self,slope)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: slope
    character(len=*), parameter :: pname = "mod_profile < set_slope"

    self%eq_profile%slope = slope

    return
  end subroutine set_slope

  subroutine set_s_slope(self,s_slope)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: s_slope
    character(len=*), parameter :: pname = "mod_profile < set_s_slope"

    self%eq_profile%s_slope = s_slope

    return
  end subroutine set_s_slope

  subroutine set_pvalue(self,pvalue)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: pvalue
    character(len=*), parameter :: pname = "mod_profile < set_pvalue"

    self%eq_profile%pvalue = pvalue

    return
  end subroutine set_pvalue

  subroutine set_r2(self,r2)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    real(rprec), intent(in) :: r2
    character(len=*), parameter :: pname = "mod_profile < set_r2"

    self%eq_profile%r2 = r2

    return
  end subroutine set_r2

  subroutine set_n(self,n)
    ! Function to set the anomalies of the profile
    class(profile_type), intent(inout) :: self
    integer(iprec), intent(in) :: n
    character(len=*), parameter :: pname = "mod_profile < set_n"

    self%eq_profile%n = n

    return
  end subroutine set_n


  ! ------------------------
  ! Other methods.
  ! ------------------------

  subroutine cut_profile(z1,z2,x,y,ox,oy)
    ! Subroutine to cut a profile to a determined depth range
    ! - x :: depths of the profile (real)
    ! - y :: flux or temp of the profile (real)
    ! - z1 :: minimum depth for cutting the profile (real)
    ! - z2 :: maximum depth for cutting the profile (real)
    ! - ox :: truncated depths (real)
    ! - oy :: truncated flux or temp (real)
    real(rprec), intent(in) :: x(:), y(:)
    real(rprec), intent(in) :: z1, z2
    real(rprec), allocatable, intent(out) :: ox(:), oy(:)

    integer(iprec) :: n, i1, i2 ! Dummy integer
    integer(iprec) :: i ! Index

    character(len=*), parameter :: pname = "mod_profile < cut_profile"

    ! Initialize logger
    call init_logger()

    call book%record(pname,trace,'Starting subroutine')

    ! Check the number of depths in the range
    n = 0_iprec
    i1 = 0_iprec
    i2 = 0_iprec
    do i = 1, size(x)
      if(x(i).ge.z1.and.x(i).le.z2) then
        n = n + 1_iprec
        if(i1.eq.0_iprec) i1 = i
      end if
    end do ! i
    do i = size(x), 1, -1
      if(x(i).ge.z1.and.x(i).le.z2) then
        if(i2.eq.0_iprec) then
          i2 = i
          exit
        end if
      end if
    end do
    call book%record(pname,debug,"zmin and zmax = ",[z1,z2])
    call book%record(pname,debug,"i1 and i2 = ",[i1,i2])
    call book%record(pname,debug,"z(i1) and z(i2) = ",[x(i1),x(i2)])

    ! Select suitable records
    if(i1.lt.i2) then
      ox = x(i1:i2)
      oy = y(i1:i2)
    elseif(i2.lt.i1) then
      ox = x(i2:i1)
      oy = y(i2:i1)
    elseif(i1.eq.i2) then
      call book%record(pname,error,'ERROR - zmin is equal to zmax!')
      stop
    end if

    ! Check the resuls
    if(size(ox).ne.n) then
      call book%record(pname,error,'ERROR - Something went wrong when cutting profile')
      call book%record(pname,error,'Size cutted profile:',size(ox))
      call book%record(pname,error,'Number depths in rnage:',n)
      call book%record(pname,debug,"zmin and zmax = ",[z1,z2])
      call book%record(pname,debug,"i1 and i2 = ",[i1,i2])
      call book%record(pname,debug,"z(i1) and z(i2) = ",[x(i1),x(i2)])
      stop
    end if

    call book%record(pname,trace,'Ending subroutine')

    return
  end subroutine cut_profile

  subroutine anomaly(self)
    ! Subroutine to estmate the three anomaly profiles
    ! - self = measured profile from which to estimate the
    !                    anomaly profiles (profile)
    class(profile_type), intent(inout) :: self

    type(regression_type) :: linfit ! Dummy regression
    real(rprec), allocatable :: mx(:), my(:) ! Dummy reals
    real(rprec) :: slope, interc, s_slope, s_interc ! Dummy reals

    character(len=*), parameter :: pname = "mod_profile < anomaly"

    call init_logger

    call book%record(pname,trace,'Starting subroutine')

    call cut_profile(self%zmin,self%zmax,self%depth,self%temp,mx,my) ! Desired depth window

    ! The fit of the depth window (temp = y, depth = x)
    call linfit%new(my,mx)
    call linfit%linreg()

    ! We have some parts of the anomaly_profiles object
    slope = linfit%coefficients(2)
    interc = linfit%coefficients(1)
    s_slope = linfit%sgm_coefficients(2)
    s_interc = linfit%sgm_coefficients(1)

    self%eq_profile%inter = interc ! The regresion part
    self%eq_profile%s_inter = s_interc ! The regresion part
    self%eq_profile%slope = slope ! The regresion part
    self%eq_profile%s_slope = s_slope ! The regresion part
    self%eq_profile%pvalue = linfit%pvalue ! The regresion part
    self%eq_profile%r2 = linfit%r2 ! The regresion part
    self%eq_profile%n = linfit%n ! The regresion part

    ! The first anomaly, the one from the original profile
    call define_anomaly(self,0_iprec,2.0_rprec,self%atemp_original)

    ! Two cases now: negative slope or otherwise
    call define_anomaly(self,-1_iprec,2.0_rprec,self%atemp_minus)
    call define_anomaly(self,1_iprec,2.0_rprec,self%atemp_plus)

    call book%record(pname,trace,'Exiting subroutine')

    return
  end subroutine anomaly

  subroutine define_anomaly(self,ind,sig,ty)
    ! Subroutine to obatin anomaly profiles from a profile object
    ! - self :: profile object with the data
    ! - sig :: number of sigma
    ! - ind :: indicator to select method to create anomaly profile
    ! - ty :: array to save anomaly profile
    class(profile_type), intent(in) :: self
    integer(iprec), intent(in) :: ind
    real(rprec), intent(in) :: sig
    real(rprec), allocatable, intent(out) :: ty(:)

    real(rprec) :: inter_minus, inter_plus, slope_minus, slope_plus

    if(ind.eq.0_iprec) then
      ty = self%temp(:) - (self%eq_profile%inter + self%eq_profile%slope *&
          self%depth(:) )
    elseif(ind.eq.-1_iprec) then

      if(self%eq_profile%slope.lt.0.0_rprec) then
        inter_minus = self%eq_profile%inter - sig * self%eq_profile%s_inter
        slope_minus = self%eq_profile%slope - sig * self%eq_profile%s_slope

        ty = self%temp(:) - ( inter_minus + slope_minus * self%depth(:) )
      else
        inter_minus = self%eq_profile%inter + sig * self%eq_profile%s_inter
        slope_minus = self%eq_profile%slope - sig * self%eq_profile%s_slope

        ty = self%temp(:) - ( inter_minus + slope_minus * self%depth(:) )
      end if

    elseif(ind.eq.1_iprec) then

      if(self%eq_profile%slope.lt.0.0_rprec) then
        inter_plus = self%eq_profile%inter + sig * self%eq_profile%s_inter
        slope_plus = self%eq_profile%slope + sig * self%eq_profile%s_slope

        ty = self%temp(:) - ( inter_plus + slope_plus * self%depth(:) )
      else
        inter_plus = self%eq_profile%inter - sig * self%eq_profile%s_inter
        slope_plus = self%eq_profile%slope + sig * self%eq_profile%s_slope

        ty = self%temp(:) - ( inter_plus + slope_plus * self%depth(:) )
      end if

    end if

    return
  end subroutine define_anomaly

  subroutine point_h(self,sig,th,zh)
    ! Subroutine to estimate the corssing point for extremal profiles
    !    (H point, Hugo's point)
    ! - self :: profile object with all data.
    ! - th :: temperature of H point
    ! - zh :: depth of the H point
    class(profile_type), intent(in) :: self
    real(rprec), intent(in) :: sig
    real(rprec), intent(out) :: th
    real(rprec), intent(out) :: zh

    real(rprec) :: inter_minus, inter_plus, slope_minus, slope_plus
    real(rprec) :: t_plus, t_minus

    character(len=*), parameter :: pname = "mod_profile < point_h"

    call init_logger


    if(self%eq_profile%slope.lt.0.0_rprec) then
      inter_minus = self%eq_profile%inter - sig * self%eq_profile%s_inter
      slope_minus = self%eq_profile%slope - sig * self%eq_profile%s_slope

    else
      inter_minus = self%eq_profile%inter + sig * self%eq_profile%s_inter
      slope_minus = self%eq_profile%slope - sig * self%eq_profile%s_slope

    end if

    if(self%eq_profile%slope.lt.0.0_rprec) then
      inter_plus = self%eq_profile%inter + sig * self%eq_profile%s_inter
      slope_plus = self%eq_profile%slope + sig * self%eq_profile%s_slope

    else
      inter_plus = self%eq_profile%inter - sig * self%eq_profile%s_inter
      slope_plus = self%eq_profile%slope + sig * self%eq_profile%s_slope

    end if

    zh = (inter_minus - inter_plus) / (slope_plus - slope_minus)

    t_plus = inter_plus + slope_plus * zh
    t_minus = inter_plus + slope_plus * zh

    if(abs(t_plus-t_minus).gt.0.001_rprec) then
      call book%record(pname,error,"ERROR: cannot estimate th")
      stop "Cannot estimate th"
    end if

    th = t_minus

    return
  end subroutine point_h

  subroutine point_b(self,tb,zb)
    ! Subroutine to estimate the lowermost point of the regrssion line
    ! - self :: profile object with all data.
    ! - tb :: temperature of B point
    ! - zb :: depth of the B point
    class(profile_type), intent(in) :: self
    real(rprec), intent(out) :: tb
    real(rprec), intent(out) :: zb

    integer(iprec) :: nz

    character(len=*), parameter :: pname = "mod_profile < point_b"

    call init_logger


    nz = size(self%depth,1)
    tb = self%eq_profile%inter + self%eq_profile%slope * self%depth(nz)
    zb = self%depth(nz)

    return
  end subroutine point_b

end module profile_module




