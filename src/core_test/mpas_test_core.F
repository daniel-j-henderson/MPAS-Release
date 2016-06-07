! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
module test_core

   use mpas_framework
   use mpas_timekeeping

   type (MPAS_Clock_type), pointer :: clock

   contains

   !***********************************************************************
   !
   !  function test_core_init
   !
   !> \brief   MPAS Core Initialization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs the necessary initialization of a core. This can
   !>   involve things like setting up coeffiecients for remapping and advection,
   !>   but more generally can include anything the core needs to initialize.
   !
   !-----------------------------------------------------------------------
   function test_core_init(domain, startTimeStamp) result(iErr)!{{{
   
      use mpas_configure
      use mpas_derived_types
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      character(len=*), intent(out) :: startTimeStamp
   
      type (MPAS_Time_Type) :: startTime
      integer :: iErr

      character (len=StrKIND) :: datestring
      type (MPAS_Time_Type) :: ref_time
      integer :: n
      type (MPAS_timeInterval_type) :: num, den, rem
      integer(kind=8) :: start_time, stop_time, rate

      datestring = "0000-01-01_00:00:00"
      call mpas_set_time(ref_time, dateTimeString=datestring)
      call mpas_set_timeInterval(num, YY=34)
      call mpas_set_timeInterval(den, S=1)


      ! For testing purposes, I've left both the old and new methods, so you can ensure that the 
      ! results are correct and compare execution times. Omit the final argument (or use '1') 
      ! to use the old method, use '2' for the new method.

      call system_clock(start_time, rate)
      call mpas_interval_division(ref_time, num, den, n, rem, 1)
      call system_clock(stop_time)
      call mpas_get_timeInterval(rem, timeString=datestring)
      write(0,*) 'elapsed interval division time: ', real((stop_time - start_time),kind=RKIND) / real(rate,kind=RKIND)
      write(0,*) "n = ", n, " rem = ", datestring

      call system_clock(start_time)
      call mpas_interval_division(ref_time, num, den, n, rem, 2)
      call system_clock(stop_time)
      call mpas_get_timeInterval(rem, timeString=datestring)
      write(0,*) 'New method elapsed interval division time: ', real((stop_time - start_time),kind=RKIND) / real(rate,kind=RKIND)
      write(0,*) "n = ", n, " rem = ", datestring

      iErr = 0

      !
      ! Set "local" clock to point to the clock contained in the domain type
      !
      clock => domain % clock

      !
      ! Set startTimeStamp based on the start time of the simulation clock
      !
      startTime = mpas_get_clock_time(clock, MPAS_START_TIME, iErr)
      call mpas_get_time(startTime, dateTimeString=startTimeStamp) 

   end function test_core_init!}}}


   !***********************************************************************
   !
   !  function test_core_run
   !
   !> \brief   MPAS Core Run Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function performs a run of the MPAS core. This can involve time
   !>  stepping if the core requires time stepping, but more generally includes
   !>  anything a core would typically do after it was initialized.
   !
   !-----------------------------------------------------------------------
   function test_core_run(domain) result(iErr)!{{{
   
      use mpas_derived_types
      use mpas_kind_types
      use mpas_timer
   
      implicit none
   
      type (domain_type), intent(inout) :: domain
      integer :: iErr

      type (mpas_pool_type), pointer :: pool
      type (mpas_pool_iterator_type) :: itr


      iErr = 0

   end function test_core_run!}}}

   !***********************************************************************
   !
   !  function test_core_finalize
   !
   !> \brief   MPAS Core Finalization Function
   !> \author  Doug Jacobsen
   !> \date    04/06/2015
   !> \details 
   !>  This function finalizes the MPAS core. It should at a minimum destroy the
   !>  simulation clock, but can perform any functions necessary to clean up the
   !>  MPAS core.
   !
   !-----------------------------------------------------------------------
   function test_core_finalize(domain) result(iErr)!{{{
   
      use mpas_derived_types
   
      implicit none

      type (domain_type), intent(inout) :: domain 
      integer :: iErr
 

      iErr = 0

      call mpas_destroy_clock(clock, iErr)

   end function test_core_finalize!}}}

end module test_core
