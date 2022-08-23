!
!----------------------------------------------------------------------
! Module to estimate ground heat flux from ground temperature series.
!
! This module needs:
! - mod_kinds.f03
! - mod_logger.f03
! - mod_profile.f03 (it requires mod_regression.f03)
! - mod_fm (it requires mod_profile.f03)
!
! Francisco Jose Cuesta-Valero
! 2021-12-14 (Ritterstrasse, Leipzig)
!----------------------------------------------------------------------
module ghf_module
  use kinds_module
  use logger_module
  use fm_module
  implicit none

  private

  type, public :: ghf_type
    ! Thermal diffusivity m-2 s-1 or appropriate units
    real(rprec) :: alpha
    ! Thermal conductivity W m-1 K-1 or appropriate units
    real(rprec) :: lambda
    ! Thermal inertia J m-2 K-1 s-1/2 or appropriate units
    real(rprec) :: ti
    ! Ground heat flux
    real(rprec), allocatable :: gflux(:,:)
    ! Forward model for Surface Ground heat flux
    real(rprec), allocatable :: fflux(:,:)
  contains
    procedure :: new_ghf, new_copy_ghf
    generic :: new => new_ghf, new_copy_ghf
    procedure :: delete
    procedure :: flux_1d, flux_2d, flux_nd, flux_ndp
    generic :: flux => flux_1d, flux_2d, flux_nd, flux_ndp
    procedure :: forward
    procedure :: put => write_ghf
  end type ghf_type

  type(logger_type) :: book
  integer(iprec), parameter :: log_lvl = warning

contains
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
  subroutine new_ghf(self,alpha,lambda)
    ! Subrouine to create ghf objects
    ! - self :: object to be created (ghf)
    ! - alpha :: thermal diffusivity of the medioum (real)
    ! - lambda :: thermal conductivity of the medioum (real)
    class(ghf_type), intent(inout) :: self
    real(rprec), intent(in) :: alpha
    real(rprec), intent(in) :: lambda

    self%alpha = alpha
    self%lambda = lambda
    self%ti = lambda / sqrt(alpha)

    return
  end subroutine new_ghf

  subroutine new_copy_ghf(self,other)
    ! Subroutine to copy ghf objects
    ! - self :: object to be created (ghf)
    ! - other :: object to be copied (ghf)
    class(ghf_type), intent(out) :: self
    class(ghf_type), intent(in) :: other

    self%alpha = other%alpha
    self%lambda = other%lambda
    self%ti = other%ti
    if(allocated(other%gflux)) self%gflux = other%gflux
    if(allocated(other%fflux)) self%fflux = other%fflux

    return
  end subroutine new_copy_ghf

  subroutine delete(self)
    ! Subroutine to delete ghf objects
    ! - self :: object to be deleted
    class(ghf_type), intent(inout) :: self

    if(allocated(self%gflux)) then
      deallocate(self%gflux)
    end if

    if(allocated(self%fflux)) then
      deallocate(self%fflux)
    end if

    return
  end subroutine delete

  ! ------------------------
  ! Principal subroutines
  ! ------------------------
  function cumweitemp(temp,time) result(cwt)
    ! Subroutine to estimate the cumulative weighted temperature
    ! - temp :: temperature time series (real(:))
    ! - time :: temporal data (real(:))
    ! - cwt :: cumulative weighted temperature (real(:,:))
    real(rprec), intent(in) :: temp(:)
    real(rprec), intent(in) :: time(:)

    real(rprec), allocatable :: cwt(:,:)

    real(rprec), parameter :: fils = -999.9_rprec
    integer(iprec), parameter :: ifils = -999_iprec

    real(rprec) :: dtemp, dt1, dt2, dt3

    integer(iprec) :: t, tt, nt ! Dummy integers

    character(len=*), parameter :: pname = "mod_ghf < cumweitemp"

    call init_logger

    nt = size(time)
    allocate(cwt(nt-1,2))

    cwt(:,1) = time(2:nt)

    if(time(1).gt.time(nt)) then
      call book%record(pname,error,"ERROR - Time must go from past to present!")
      stop
    end if

    if(test_array(temp,fils).eq.0_iprec) then
      cwt(:,2) = fils
    else
      do t = 2, nt
        cwt(t-1,2) = 0.0_rprec
        do tt = 1, t-1
          dtemp = temp(tt+1) - temp(tt)
          dt1 = time(tt+1) - time(tt)
          dt2 = time(t) - time(tt)
          dt3 = time(t) - time(tt+1)
          cwt(t-1,2) = cwt(t-1,2) + ( dtemp / dt1 ) * (sqrt(dt2) - sqrt(dt3))
        end do ! tt
      end do ! t
    end if
    return
  end function cumweitemp

  subroutine flux_1d(self,temp,cu,logy)
    ! Subroutine to estimate flux from a temperature series
    ! - self :: flux object to save results (ghf)
    ! - temp :: temperature series (real(:))
    ! - cu :: factor to obtain SI units (optional, real)
    ! - logy :: logging year to obtain annual data (optional, real)
    class(ghf_type), intent(inout) :: self
    real(rprec), intent(in) :: temp(:)
    real(rprec), optional, intent(in) :: cu
    real(rprec), optional, intent(in) :: logy

    real(rprec), allocatable :: ctemp(:,:)
    real(rprec), allocatable :: ctime(:), y(:,:)

    real(rprec), parameter :: pi = 2.0_rprec * asin(1.0_rprec)

    integer(iprec) :: nt
    integer(iprec) :: t

    nt = size(temp,1)
    allocate(self%gflux(nt-1,2))
    allocate(ctime(nt))

    ! Create time
    do t = 1, nt
      ctime = real(t,rprec)
    end do ! t

    ! Cumulative weighted temperature
    ctemp = cumweitemp(temp,ctime)

    self%gflux(:,1) = ctemp(:,1)
    if(present(cu)) then
      self%gflux(:,2) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2) / cu
    else
      self%gflux(:,2) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2)
    end if

    if(present(logy)) then
      call inversion2year(logy,self%gflux((nt-1):1:-1,2),logy-self%gflux((nt-1):1:-1,1),y)
      deallocate(self%gflux)
      allocate(self%gflux(size(y,1),2))
      self%gflux(:,1) = y(:,1)
      self%gflux(:,2) = y(:,2)
      deallocate(y)
    end if

    return
  end subroutine flux_1d

  subroutine flux_2d(self,temp,time,cu,logy)
    ! Subroutine to estimate flux from a temperature and a time series
    ! - self :: flux object to save results (ghf)
    ! - temp :: temperature series (real(:))
    ! - time :: time series (real(:))
    ! - cu :: factor to obtain SI units (optional, real)
    ! - logy :: logging year to obtain annual data (optional, real)
    class(ghf_type), intent(inout) :: self
    real(rprec), intent(in) :: temp(:)
    real(rprec), intent(in) :: time(:)
    real(rprec), optional, intent(in) :: cu
    real(rprec), optional, intent(in) :: logy

    real(rprec), allocatable :: ctemp(:,:), y(:,:)

    real(rprec), parameter :: pi = 2.0_rprec * asin(1.0_rprec)

    integer(iprec) :: nt

    nt = size(temp,1)
    allocate(self%gflux(nt-1,2))

    ! Cumulative weighted temperature
    ctemp = cumweitemp(temp,time)

    self%gflux(:,1) = ctemp(:,1)
    if(present(cu)) then
      self%gflux(:,2) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2) / cu
    else
      self%gflux(:,2) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2)
    end if

    if(present(logy)) then
      call inversion2year(logy,self%gflux((nt-1):1:-1,2),logy-self%gflux((nt-1):1:-1,1),y)
      deallocate(self%gflux)
      allocate(self%gflux(size(y,1),2))
      self%gflux(:,1) = y(:,1)
      self%gflux(:,2) = y(:,2)
      deallocate(y)
    end if


    return
  end subroutine flux_2d

  subroutine flux_nd(self,temp,cu,logy)
    ! Subroutine to estimate flux from a temperature matrix where the first column is time
    ! - self :: flux object to save results (ghf)
    ! - temp :: temperature series (real(:,:))
    ! - cu :: factor to obtain SI units (optional, real)
    ! - logy :: logging year to obtain annual data (optional, real)
    class(ghf_type), intent(inout) :: self
    real(rprec), intent(in) :: temp(:,:)
    real(rprec), optional, intent(in) :: cu
    real(rprec), optional, intent(in) :: logy

    real(rprec), allocatable :: ctemp(:,:), y(:,:), yy(:,:)

    real(rprec), parameter :: pi = 2.0_rprec * asin(1.0_rprec)

    integer(iprec) :: nt, ntemp, i

    nt = size(temp,1)
    ntemp = size(temp,2)
    allocate(self%gflux(nt-1,ntemp))

    do i = 2, ntemp
      ! Cumulative weighted temperature
      ctemp = cumweitemp(temp(:,i),temp(:,1))

      if(i.eq.2) then
        self%gflux(:,1) = ctemp(:,1)
      end if
      if(present(cu)) then
        self%gflux(:,i) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2) / cu
      else
        self%gflux(:,i) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2)
      end if

      deallocate(ctemp)
    end do ! i

    if(present(logy)) then
      yy = self%gflux
      deallocate(self%gflux)
      call inversion2year(logy,yy((nt-1):1:-1,2),logy-yy((nt-1):1:-1,1),y)
      allocate(self%gflux(size(y,1),ntemp))
      self%gflux(:,1) = y(:,1)
      deallocate(y)
      do i = 2, ntemp
        call inversion2year(logy,yy((nt-1):1:-1,i),logy-yy((nt-1):1:-1,1),y)
        self%gflux(:,i) = y(:,2)
        deallocate(y)
      end do
    end if

    return
  end subroutine flux_nd

  subroutine flux_ndp(self,temp,time,cu,logy)
    ! Subroutine to estimate flux from a temperature matrix with an additional vector containing time
    ! - self :: flux object to save results (ghf)
    ! - temp :: temperature series (real(:,:))
    ! - time :: time series (real(:))
    ! - cu :: factor to obtain SI units (optional, real)
    ! - logy :: logging year to obtain annual data (optional, real)
    class(ghf_type), intent(inout) :: self
    real(rprec), intent(in) :: temp(:,:)
    real(rprec), intent(in) :: time(:)
    real(rprec), optional, intent(in) :: cu
    real(rprec), optional, intent(in) :: logy

    real(rprec), allocatable :: ctemp(:,:), y(:,:), yy(:,:)

    real(rprec), parameter :: pi = 2.0_rprec * asin(1.0_rprec)

    integer(iprec) :: nt, ntemp, i

    nt = size(time,1)
    ntemp = size(temp,2)
    allocate(self%gflux(nt-1,ntemp+1))

    do i = 1, ntemp
      ! Cumulative weighted temperature
      ctemp = cumweitemp(temp(:,i),time)

      if(i.eq.1) then
        self%gflux(:,1) = ctemp(:,1)
      end if
      if(present(cu)) then
        self%gflux(:,i+1) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2) / cu
      else
        self%gflux(:,i+1) = (2.0_rprec * self%ti / sqrt(pi)) * ctemp(:,2)
      end if

      deallocate(ctemp)
    end do ! i

    if(present(logy)) then
      yy = self%gflux
      deallocate(self%gflux)
      call inversion2year(logy,yy((nt-1):1:-1,2),logy-yy((nt-1):1:-1,1),y)
      allocate(self%gflux(size(y,1),ntemp+1))
      self%gflux(:,1) = y(:,1)
      deallocate(y)
      do i = 1, ntemp
        call inversion2year(logy,yy((nt-1):1:-1,i+1),logy-yy((nt-1):1:-1,1),y)
        self%gflux(:,i+1) = y(:,2)
        deallocate(y)
      end do
    end if

    return
  end subroutine flux_ndp

  subroutine forward(self,depth)
    ! Subroutine to obtain the forward model of the ghf surface signal
    ! - self :: object containing the surface flux signal (ghf)
    ! - depth :: depth profile for the forward model (real(:))
    class(ghf_type), intent(inout) :: self
    real(rprec), intent(in) :: depth(:)

    type(fm_type) :: fm

    real(rprec), allocatable :: tseries(:), signal(:), z(:), y(:,:)

    integer(iprec) :: i, nt, nz, nf

    character(len=*), parameter :: pname = "mod_ghf < forward"

    call init_logger

    if(.not.allocated(self%gflux)) then
      call book%record(pname,error,"ERROR - GHF object doesn't include surface flux signal")
      stop
    end if

    nt = size(self%gflux,1)
    nf = size(self%gflux,2)
    nz = size(depth,1)

    tseries = maxval(self%gflux(:,1)) - self%gflux(:,1)
    if(tseries(2).lt.tseries(1)) then
      call book%record(pname,warning,"WARNING - time must go from past to present")
      y = self%gflux(nt:1:-1,:)
      tseries = tseries(nt:1:-1)
    else
      y = self%gflux
    end if


    allocate(self%fflux(nz,nf))

    self%fflux(:,1) = depth
    z = depth

    do i = 2, nf
      signal = y(:,i)
      call fm%new(self%alpha,tseries,signal,z)
      call fm%forward_model()

      self%fflux(:,i) = fm%prof%atemp_original
    end do ! i

    return
  end subroutine forward

  subroutine write_ghf(self,ofile,u2,oo)
    ! Subroutine to save a ghf object into a file
    ! - self :: object to be saved (ghf)
    ! - ofile :: name of the file to save the object (char)
    ! - u2 :: unit to write the file (integer)
    ! - oo :: oo .eq. 0 means do not sort the results
    class(ghf_type), intent(in) :: self
    character(len=*), intent(in) :: ofile
    integer(iprec), optional, intent(in) :: u2
    integer(iprec), optional, intent(in) :: oo

    real(rprec), parameter :: fils = -999.9_rprec
    integer(iprec), parameter :: ifils = -999_iprec

    character(90) :: f ! format for printing
    character(9) :: ntemp_char

    real(rprec), allocatable :: ens(:)

    integer(iprec) :: u, ntemp, nt, nz, o ! Dummy integer
    integer(iprec) :: i, j ! Dummy indices

    call init_logger

    if(present(u2)) then
      u = u2
    else
      u = 90_iprec
    end if

    if(present(oo)) then
      o = oo
    else
      o = 1_iprec
    end if

    open(u,file=trim(ofile),action="write")

    ! First, relevant parameters of the model
    f = "(a,es12.3e3)"
    write(u,trim(f)) '#alpha', self%alpha
    write(u,trim(f)) '#lambda', self%lambda
    write(u,trim(f)) '#ti', self%ti

    write(u,*) ""
    write(u,*) ""

    ! Write flux soluions
    nt = size(self%gflux,1)
    ntemp = size(self%gflux,2)
    write(ntemp_char,'(i0.3)') ntemp

    f = "(a,"//trim(ntemp_char)//"es12.3e3)"
    allocate(ens(ntemp-1))
    do i = 1, nt
      ens(:) = self%gflux(i,2:ntemp)

      if(o.eq.1_iprec) then
        call sort(ntemp-1,ens)
      end if

      write(u,trim(f)) '#gflux', self%gflux(i,1), (ens(j), j=1,ntemp-1)
    end do ! i
    deallocate(ens)


    ! Write flux forward model
    if(allocated(self%fflux)) then
      write(u,*) ""
      write(u,*) ""

      nz = size(self%fflux,1)
      ntemp = size(self%fflux,2)
      write(ntemp_char,'(i0.3)') ntemp

      f = "(a,"//trim(ntemp_char)//"es12.3e3)"
      allocate(ens(ntemp-1))
      do i = 1, nz
        ens(:) = self%fflux(i,2:ntemp)

        if(o.eq.1_iprec) then
          call sort(ntemp-1,ens)
        end if

        write(u,trim(f)) '#fflux', self%fflux(i,1), (ens(j), j=1,ntemp-1)
      end do ! i
    end if

    return
  end subroutine write_ghf

  ! ------------------------
  ! Auxiliar subroutines
  ! ------------------------
  function test_array(x,fils) result(t)
    ! Function to check that 1D arrays does not include strange results.
    ! - x => the array to be checked (real(:))
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

  subroutine inversion2year(logy,x,xtime,y)
    real(rprec), intent(in) :: logy
    real(rprec), intent(in) :: x(:), xtime(:)

    real(rprec), allocatable, intent(out) :: y(:,:)

    integer(iprec) :: n, nt
    integer(iprec) :: i, j

    character(len=*), parameter :: pname="mod_ghf < inversion2year"

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


end module ghf_module
