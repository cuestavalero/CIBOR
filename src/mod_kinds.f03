!
!----------------------------------------------------------------------
! This module provides the kinds of data to be defined in other modules
! and programs:
!   1) SPI: singular precision integer
!   2) DPI: doble precision integer
!   3) SPR: singular precision real
!   4) DPR: doble precision real
!
! Aug-11-2016 (Antigonish)
! Francisco Jose Cuesta Valero
!----------------------------------------------------------------------

module kinds_module
  implicit none
  private

  integer, parameter :: spi = selected_int_kind(6)
  integer, parameter :: dpi = selected_int_kind(15)
  integer, parameter :: spr = selected_real_kind(6)
  integer, parameter :: dpr = selected_real_kind(15)

  integer, parameter, public :: iprec = spi
  integer, parameter, public :: rprec = dpr
 
end module kinds_module


