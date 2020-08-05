!
!-----------------------------------------------------------------------
!
!Purpose:
! Fortran-2003 bindings for lcm -- part 2.
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
subroutine LCMPUT(iplist, name, ilong, itype, idata)
   ! store a record in an associative table
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer :: ilong, itype
   integer, target, dimension(*) :: idata
   interface
      subroutine lcmput_c (iplist, namp, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         integer(c_int), value :: ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmput_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmput_c(iplist, name13, ilong, itype, pt_data)
end subroutine LCMPUT
!
subroutine LCMGET(iplist, name, idata)
   ! recover a record from an associative table
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) :: iplist, pt_data
   character(len=*),intent(in) :: name
   character(kind=c_char), dimension(13) :: name13
   integer, target, dimension(*) :: idata
   interface
      subroutine lcmget_c (iplist, namp, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         character(kind=c_char), dimension(*) :: namp
         type(c_ptr), value :: idata
      end subroutine lcmget_c
   end interface
   pt_data=c_loc(idata)
   call STRCUT(name13, name)
   call lcmget_c(iplist, name13, pt_data)
end subroutine LCMGET
!
subroutine LCMPDL(iplist, ipos, ilong, itype, idata)
   ! store a record in an heterogeneous list
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) :: iplist, pt_data
   integer :: ipos, ilong, itype
   integer, target, dimension(*) :: idata
   interface
      subroutine lcmpdl_c (iplist, ipos, ilong, itype, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos, ilong, itype
         type(c_ptr), value :: idata
      end subroutine lcmpdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmpdl_c(iplist, ipos-1, ilong, itype, pt_data)
end subroutine LCMPDL
!
subroutine LCMGDL(iplist, ipos, idata)
   ! recover a record from an heterogeneous list
   use, intrinsic :: iso_c_binding
   use LCMAUX
   type(c_ptr) :: iplist, pt_data
   integer :: ipos
   integer, target, dimension(*) :: idata
   interface
      subroutine lcmgdl_c (iplist, ipos, idata) bind(c)
         use, intrinsic :: iso_c_binding
         type(c_ptr) :: iplist
         integer(c_int), value :: ipos
         type(c_ptr), value :: idata
      end subroutine lcmgdl_c
   end interface
   pt_data=c_loc(idata)
   call lcmgdl_c(iplist, ipos-1, pt_data)
end subroutine LCMGDL
