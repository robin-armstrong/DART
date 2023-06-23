! DART software - Copyright UCAR. This open source software is provided
! by UCAR, "as is", without charge, subject to all terms of use at
! http://www.image.ucar.edu/DAReS/DART/DART_download
!
! $Id$

! This file is meant to read a text file containing bottle data from the
! Bermuda Atlantic Time-Series Study (https://bats.bios.asu.edu/).

program text_to_obs

use         types_mod, only : r8, PI, DEG2RAD
use     utilities_mod, only : initialize_utilities, finalize_utilities, &
                              open_file, close_file, &
                              find_namelist_in_file, check_namelist_read, &
                              error_handler, E_ERR, E_MSG, nmlfileunit,   &
                              do_nml_file, do_nml_term
use  time_manager_mod, only : time_type, set_calendar_type, set_date, &
                              operator(>=), increment_time, get_time, &
                              operator(-), NOLEAP, GREGORIAN, operator(+), &
                              print_date
use      location_mod, only : VERTISHEIGHT, VERTISPRESSURE
use  obs_sequence_mod, only : obs_sequence_type, obs_type, read_obs_seq, &
                              static_init_obs_sequence, init_obs, write_obs_seq, & 
                              init_obs_sequence, get_num_obs, & 
                              set_copy_meta_data, set_qc_meta_data

use      obs_kind_mod, only : POLY_ELECTRODE_OXYGEN, TITRATION_ALKALINITY, &
                              CATALYTIC_CARBON, UV_OXY_NITROGEN

implicit none

integer, parameter :: NUM_SCALAR_OBS = 4  ! maximum number of scalar observation variables that will
                                          ! be assimilated at each observation.

! this array defines the order in which observations are read from the file
integer, parameter :: OTYPE_ORDERING(NUM_SCALAR_OBS) &
                      = (/POLY_ELECTRODE_OXYGEN, TITRATION_ALKALINITY, CATALYTIC_CARBON, UV_OXY_NITROGEN/)

! namelist variables, changeable at runtime
character(len=256) :: text_input_file, obs_out_file
integer :: max_lines, read_starting_at_line, date_firstcol, hourminute_firstcol
integer :: lat_cols(2), lon_cols(2), vert_cols(2)
integer :: scalar_obs_cols(2, NUM_SCALAR_OBS)
real(8) :: obs_uncertainties(NUM_SCALAR_OBS)
logical :: debug

namelist /text_to_obs_nml/ text_input_file, max_lines, read_starting_at_line, date_firstcol, &
                           hourminute_firstcol, lat_cols, lon_cols, vert_cols, scalar_obs_cols, &
                           obs_uncertainties, obs_out_file, debug

! local variables
character (len=294) :: input_line

integer :: oday, osec, rcio, iunit, otype, line_number, otype_index
integer :: year, month, day, hour, minute, second, hourminute_raw, date_raw
integer :: num_copies, num_qc, max_obs
           
logical  :: file_exist, first_obs

real(r8) :: temp, terr, qc, wdir, wspeed, werr
real(r8) :: lat, lon, vert, uwnd, uerr, vwnd, verr, ovalue

! the uncertainties corresponding to the observations above

type(obs_sequence_type) :: obs_seq
type(obs_type)          :: obs, prev_obs
type(time_type)         :: ref_day0, time_obs, prev_time


! start of executable code

call initialize_utilities('text_to_obs')

! Read the namelist entries
call find_namelist_in_file("input.nml", "text_to_obs_nml", iunit)
read(iunit, nml = text_to_obs_nml, iostat = rcio)
call check_namelist_read(iunit, rcio, "text_to_obs_nml")

! Record the namelist values used for the run 
if (do_nml_file()) write(nmlfileunit, nml=text_to_obs_nml)
if (do_nml_term()) write(     *     , nml=text_to_obs_nml)

! time setup
call set_calendar_type(GREGORIAN)

! open input text file

iunit = open_file(text_input_file, 'formatted', 'read')
if (debug) print *, 'opened input file ' // trim(text_input_file)

max_obs    = NUM_SCALAR_OBS*max_lines
num_copies = 1
num_qc     = 1

! call the initialization code, and initialize two empty observation types
call static_init_obs_sequence()
call init_obs(obs,      num_copies, num_qc)
call init_obs(prev_obs, num_copies, num_qc)
first_obs = .true.

! create a new, empty obs_seq file.
call init_obs_sequence(obs_seq, num_copies, num_qc, max_obs)

! the first one needs to contain the string 'observation' and the
! second needs the string 'QC'.
call set_copy_meta_data(obs_seq, 1, 'observation')
call set_qc_meta_data(obs_seq, 1, 'Data QC')

! Set the DART data quality control.   0 is good data. 
! increasingly larger QC values are more questionable quality data.
qc = 0.0_r8

line_number = 0 ! counts the number of lines that have been read so far

obsloop: do    ! no end limit - have the loop break when input ends
   ! read in entire text line into a buffer
   read(iunit, "(A)", iostat=rcio) input_line
   line_number = line_number + 1

   if(line_number < read_starting_at_line) then
      cycle obsloop
   end if

   if (rcio /= 0) then 
      if (debug) print *, 'got bad read code from input file at line ',line_number,', rcio = ', rcio
      exit obsloop
   endif

   ! extracting the date when the observation was taken

   read(input_line(date_firstcol:(date_firstcol + 7)), *, iostat=rcio) date_raw
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting date at line ",line_number,", rcio = ",rcio
      exit obsloop

   else if(date_raw == -999) then
      cycle obsloop ! missing date

   end if

   read(input_line(date_firstcol:(date_firstcol + 3)), *, iostat=rcio) year
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting year at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   read(input_line((date_firstcol + 4):(date_firstcol + 5)), *, iostat=rcio) month
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting month at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   read(input_line((date_firstcol + 6):(date_firstcol + 7)), *, iostat=rcio) day
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting day at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   if((month == 02) .and. (day == 29)) then  ! we do not assimilate data taken on leap days
      cycle obsloop
   end if

   ! extracting the time of day when the observation was taken

   read(input_line(hourminute_firstcol:(hourminute_firstcol + 3)), *, iostat=rcio) hourminute_raw
   if(rcio /= 0) then
      if(debug) print *, "got bad read code parsing raw hour-minute at line ",line_number,", rcio = ",rcio
      exit obsloop

   else if(hourminute_raw == -999) then
      cycle obsloop  ! missing timestamp

   else if(hourminute_raw < 60) then
      hour = 0

   else
      read(input_line(hourminute_firstcol:(hourminute_firstcol + 1)), *, iostat=rcio) hour

      if(rcio /= 0) then
         if(debug) print *, "got bad read code getting hour at line ",line_number,", rcio = ",rcio
         exit obsloop

      end if
   end if

   read(input_line((hourminute_firstcol + 2):(hourminute_firstcol + 3)), *, iostat=rcio) minute
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting minute at line ",line_number,", rcio = ",rcio
      exit obsloop
   end if

   ! extracting the observation location

   read(input_line(lat_cols(1):lat_cols(2)), *, iostat=rcio) lat
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting lat at line ",line_number,", rcio = ",rcio
      exit obsloop
      
   else if(abs(lat + 999.0) < 1.0d-8) then
      cycle obsloop ! missing latitude

   end if

   read(input_line(lon_cols(1):lon_cols(2)), *, iostat=rcio) lon
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting lon at line ",line_number,", rcio = ",rcio
      exit obsloop

   else if(abs(lon + 999.0) < 1.0d-8) then
      cycle obsloop ! missing longitude

   end if

   read(input_line(vert_cols(1):vert_cols(2)), *, iostat=rcio) vert
   if(rcio /= 0) then
      if(debug) print *, "got bad read code getting vert at line ",line_number,", rcio = ",rcio
      exit obsloop

   else if(abs(vert + 999.0) < 1.0d-8) then
      cycle obsloop ! missing vertical coordinate

   end if

   ! extracting the observation values

   otype_loop : do otype_index = 1, NUM_SCALAR_OBS
      read(input_line(scalar_obs_cols(1, otype_index):scalar_obs_cols(2, otype_index)), *, iostat=rcio) ovalue
      
      if(rcio /= 0) then
         if(debug) print *, "got bad read code getting observation type ",otype_index," at line ",line_number,", rcio = ",rcio
         exit obsloop

      else if(abs(ovalue + 999.0) < 1.0d-8) then
         cycle otype_loop ! missing observation

      end if

      time_obs = set_date(year, month, day, hours=hour, minutes=minute)
      call get_time(time_obs, osec, days=oday)

      if(debug) then
         call print_date(time_obs, "adding observation taken on")
         print *, " \__ observation type:  ",OTYPE_ORDERING(otype_index)
         print *, "     lat:               ",lat
         print *, "     lon:               ",lon
         print *, "     vert:              ",vert
         print *, "     observation value: ",ovalue
      end if

      call create_3d_obs(lat, lon, vert, VERTISHEIGHT, &
                         ovalue, OTYPE_ORDERING(otype_index), obs_uncertainties(otype_index), &
                         oday, osec, qc, obs)

      call add_obs_to_seq(obs_seq, obs, time_obs, prev_obs, prev_time, first_obs)
   end do otype_loop
end do obsloop

! if we added any obs to the sequence, write it out to a file now.
if ( get_num_obs(obs_seq) > 0 ) then
   if (debug) print *, 'writing obs_seq, obs_count = ', get_num_obs(obs_seq)
   call write_obs_seq(obs_seq, obs_out_file)
endif

! end of main program
call finalize_utilities()

contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   create_3d_obs - subroutine that is used to create an observation
!                   type from observation data.  
!
!       NOTE: assumes the code is using the threed_sphere locations module, 
!             that the observation has a single data value and a single
!             qc value, and that this obs type has no additional required
!             data (e.g. gps and radar obs need additional data per obs)
!
! inputs:
!    lat   - latitude of observation
!    lon   - longitude of observation
!    vval  - vertical coordinate
!    vkind - kind of vertical coordinate (pressure, level, etc)
!    obsv  - observation value
!    otype - observation type
!    oerr  - observation error
!    day   - gregorian day
!    sec   - gregorian second
!    qc    - quality control value
! outputs:
!    obs   - observation type
!
!     created Oct. 2007 Ryan Torn, NCAR/MMM
!     adapted for more generic use 11 Mar 2010, nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine create_3d_obs(lat, lon, vval, vkind, obsv, otype, oerr, day, sec, qc, obs)
use        types_mod, only : r8
use obs_def_mod,      only : obs_def_type, set_obs_def_time, set_obs_def_type_of_obs, &
                             set_obs_def_error_variance, set_obs_def_location
use obs_sequence_mod, only : obs_type, set_obs_values, set_qc, set_obs_def
use time_manager_mod, only : time_type, set_time
use     location_mod, only : set_location

 integer,        intent(in)    :: otype, vkind, day, sec
 real(r8),       intent(in)    :: lat, lon, vval, obsv, oerr, qc
 type(obs_type), intent(inout) :: obs

real(r8)           :: obs_val(1), qc_val(1)
type(obs_def_type) :: obs_def

call set_obs_def_location(obs_def, set_location(lon, lat, vval, vkind))
call set_obs_def_type_of_obs(obs_def, otype)
call set_obs_def_time(obs_def, set_time(sec, day))
call set_obs_def_error_variance(obs_def, oerr * oerr)
call set_obs_def(obs, obs_def)

obs_val(1) = obsv
call set_obs_values(obs, obs_val)
qc_val(1)  = qc
call set_qc(obs, qc_val)

end subroutine create_3d_obs


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!   add_obs_to_seq -- adds an observation to a sequence.  inserts if first
!           obs, inserts with a prev obs to save searching if that's possible.
!
!     seq - observation sequence to add obs to
!     obs - observation, already filled in, ready to add
!     obs_time - time of this observation, in dart time_type format
!     prev_obs - the previous observation that was added to this sequence
!                (will be updated by this routine)
!     prev_time - the time of the previously added observation 
!                (will also be updated by this routine)
!     first_obs - should be initialized to be .true., and then will be
!                updated by this routine to be .false. after the first obs
!                has been added to this sequence.
!
!     created Mar 8, 2010   nancy collins, ncar/image
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine add_obs_to_seq(seq, obs, obs_time, prev_obs, prev_time, first_obs)
 use        types_mod, only : r8
 use obs_sequence_mod, only : obs_sequence_type, obs_type, insert_obs_in_seq
 use time_manager_mod, only : time_type, operator(>=)

  type(obs_sequence_type), intent(inout) :: seq
  type(obs_type),          intent(inout) :: obs, prev_obs
  type(time_type),         intent(in)    :: obs_time
  type(time_type),         intent(inout) :: prev_time
  logical,                 intent(inout) :: first_obs

! insert(seq,obs) always works (i.e. it inserts the obs in
! proper time format) but it can be slow with a long file.
! supplying a previous observation that is older (or the same
! time) as the new one speeds up the searching a lot.

if(first_obs) then    ! for the first observation, no prev_obs
   call insert_obs_in_seq(seq, obs)
   first_obs = .false.
else               
   if(obs_time >= prev_time) then  ! same time or later than previous obs
      call insert_obs_in_seq(seq, obs, prev_obs)
   else                            ! earlier, search from start of seq
      call insert_obs_in_seq(seq, obs)
   endif
endif

! update for next time
prev_obs = obs
prev_time = obs_time

end subroutine add_obs_to_seq

end program text_to_obs

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision$
! $Date$
