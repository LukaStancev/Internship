!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for Ganlib support. This module defines the
! interface prototypes of the Ganlib Fortran API, defines TYPE(C_PTR)
! and defines the external functions in the Ganlib API.
!
!Copyright:
! Copyright (C) 2009 Ecole Polytechnique de Montreal
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Lesser General Public
! License as published by the Free Software Foundation; either
! version 2.1 of the License, or (at your option) any later version.
!
!Author(s): A. Hebert
!
!-----------------------------------------------------------------------
!
module GANLIB
   use FILMOD
   use, intrinsic :: iso_c_binding
   interface
      subroutine CUT(name1, name2, ilong)
         use, intrinsic :: iso_c_binding
         character(kind=c_char), dimension(*) :: name1
         character(len=*) :: name2
         integer :: ilong
      end subroutine CUT
   end interface
   interface
      subroutine FIL(name1, name2, ilong)
         use, intrinsic :: iso_c_binding
         character(len=*) :: name1
         character(kind=c_char), dimension(*) :: name2
         integer :: ilong
      end subroutine FIL
   end interface
   interface
      function LCMARA(ilong)
         use, intrinsic :: iso_c_binding
         type(c_ptr) LCMARA
         integer :: ilong
      end function LCMARA
   end interface
   interface
      subroutine LCMOP(iplist, name, imp, medium, impx)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: name
         integer imp, medium, impx
      end subroutine LCMOP
   end interface
   interface
      subroutine LCMPPD(iplist, name, ilong, itype, pt_data)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist, pt_data
         character(len=*) :: name
         integer :: ilong, itype
      end subroutine LCMPPD
   end interface
   interface
      subroutine LCMGPD(iplist, name, pt_data)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist, pt_data
         character(len=*) :: name
      end subroutine LCMGPD
   end interface
   interface
      subroutine LCMLEN(iplist, name, ilong, itylcm)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: name
         integer :: ilong, itylcm
      end subroutine LCMLEN
   end interface
   interface
      subroutine LCMINF(iplist, fnamlcm, fnammy, fempty, ilong, flcml)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: fnamlcm, fnammy
         logical :: fempty, flcml
         integer :: ilong
      end subroutine LCMINF
   end interface
   interface
      subroutine LCMNXT(iplist, name)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: name
      end subroutine LCMNXT
   end interface
   interface
      subroutine LCMVAL(iplist, name)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: name
      end subroutine LCMVAL
   end interface
   interface
      subroutine LCMDEL(iplist, name)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: name
      end subroutine LCMDEL
   end interface
   interface
      function LCMDID(iplist, name)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: LCMDID,iplist
         character(len=*) :: name
      end function LCMDID
   end interface
   interface
      function LCMLID(iplist, name, ilong)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: LCMLID,iplist
         character(len=*) :: name
         integer :: ilong
      end function LCMLID
   end interface
   interface
      function LCMDIL(iplist, ipos)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: LCMDIL,iplist
         integer :: ipos
      end function LCMDIL
   end interface
   interface
      function LCMLIL(iplist, ipos, ilong)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: LCMLIL,iplist
         integer :: ipos, ilong
      end function LCMLIL
   end interface
   interface
      function LCMGID(iplist, name)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: LCMGID,iplist
         character(len=*) :: name
      end function LCMGID
   end interface
   interface
      function LCMGIL(iplist, ipos)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: LCMGIL,iplist
         integer :: ipos
      end function LCMGIL
   end interface
   interface
      subroutine LCMSIX(iplist, name, iact)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(len=*) :: name
         integer :: iact
      end subroutine LCMSIX
   end interface
   interface
      subroutine LCMPPL(iplist, ipos, ilong, itype, pt_data)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist, pt_data
         integer :: ipos, ilong, itype
      end subroutine LCMPPL
   end interface
   interface
      subroutine LCMLEL(iplist, ipos, ilong, itylcm)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer :: ipos, ilong, itylcm
      end subroutine LCMLEL
   end interface
   interface
      subroutine LCMGPL(iplist, ipos, pt_data)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist, pt_data
         integer :: ipos
      end subroutine LCMGPL
   end interface
   interface
      subroutine LCMCL(iplist, iact)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer :: iact
      end subroutine LCMCL
   end interface
   interface
      subroutine LCMEQU(iplis1, iplis2)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplis1, iplis2
      end subroutine LCMEQU
   end interface
   interface
      subroutine LCMLIB(iplist)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
      end subroutine LCMLIB
   end interface
   interface
      subroutine LCMEXP(iplist, impx, nunit, imode, idir)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer, intent(in) :: impx, nunit, imode, idir
      end subroutine LCMEXP
   end interface
   interface
      subroutine LCMEXPV3(iplist, impx, nunit, imode, idir)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer, intent(in) :: impx, nunit, imode, idir
      end subroutine LCMEXPV3
   end interface
   interface
      subroutine XABORT(msg)
         character(len=*) :: msg
      end subroutine XABORT
   end interface
   interface
      subroutine REDGET(ityp, nitma, flott, text, dflot)
         integer :: ityp, nitma
         real :: flott
         character(len=*) :: text
         double precision :: dflot
      end subroutine REDGET
   end interface
   interface
      subroutine REDPUT(ityp, nitma, flott, text, dflot)
         integer :: ityp, nitma
         real :: flott
         character(len=*) :: text
         double precision :: dflot
      end subroutine REDPUT
   end interface
   interface
      subroutine REDOPN(iinp1, iout1, nrec)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iinp1, iout1
         integer :: nrec
      end subroutine REDOPN
   end interface
   interface
      subroutine REDCLS(iinp1, iout1, nrec)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iinp1, iout1
         integer :: nrec
      end subroutine REDCLS
   end interface
   interface
      function KDIOP(name, iactio)
         use, intrinsic :: iso_c_binding
         type(c_ptr) KDIOP
         character(len=*) :: name
         integer :: iactio
      end function KDIOP
   end interface
   interface
      function KDICL(my_file, istatu)
         use, intrinsic :: iso_c_binding
         integer(c_int) KDICL
         type(c_ptr) :: my_file
         integer :: istatu
      end function KDICL
   end interface
   interface
      integer function GANDRV(hmodul, nentry, hentry, ientry, jentry, kentry)
         use, intrinsic :: iso_c_binding
         character(len=*), intent(in) :: hmodul
         integer, intent(in) :: nentry
         character(len=12), dimension(nentry), intent(in) :: hentry
         integer, dimension(nentry), intent(in) :: ientry, jentry
         type(c_ptr), dimension(nentry), intent(in) :: kentry
      end function GANDRV
   end interface
end module GANLIB
