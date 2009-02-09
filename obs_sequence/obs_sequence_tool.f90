! Data Assimilation Research Testbed -- DART
! Copyright 2004-2007, Data Assimilation Research Section
! University Corporation for Atmospheric Research
! Licensed under the GPL -- www.gpl.org/licenses/gpl.html

! this latest addition has select by list of obs types.

program obs_sequence_tool

! <next few lines under version control, do not edit>
! $URL$
! $Id$
! $Revision
! $Date

use        types_mod, only : r8, missing_r8
use    utilities_mod, only : timestamp, register_module, initialize_utilities, &
                             find_namelist_in_file, check_namelist_read, &
                             error_handler, E_ERR, E_MSG, nmlfileunit
use     location_mod, only : location_type, get_location, &
                             LocationName !! , vert_is_height   !! set_location2
use      obs_def_mod, only : obs_def_type, get_obs_def_time, get_obs_kind, &
                             get_obs_def_location
use     obs_kind_mod, only : max_obs_kinds, get_obs_kind_name, get_obs_kind_index
use time_manager_mod, only : time_type, operator(>), print_time, set_time, &
                             print_date, set_calendar_type, GREGORIAN
use obs_sequence_mod, only : obs_sequence_type, obs_type, write_obs_seq, &
                             init_obs, assignment(=), get_obs_def, &
                             init_obs_sequence, static_init_obs_sequence, &
                             read_obs_seq_header, read_obs_seq, get_num_obs, &
                             get_first_obs, get_last_obs, get_next_obs, &
                             insert_obs_in_seq, get_num_copies, get_num_qc, &
                             get_copy_meta_data, get_qc_meta_data, &
                             set_copy_meta_data, set_qc_meta_data, &
                             destroy_obs, destroy_obs_sequence, &
                             delete_seq_head, delete_seq_tail, &
                             get_num_key_range, delete_obs_by_typelist, &
                             get_obs_key, &
                             delete_obs_from_seq, get_next_obs_from_key, &
                             delete_obs_by_qc, delete_obs_by_copy

                             !%! select_obs_by_location  ! not in repository yet
implicit none

! <next few lines under version control, do not edit>
character(len=128), parameter :: &
   source   = "$URL$", &
   revision = "$Revision$", &
   revdate  = "$Date$"

type(obs_sequence_type) :: seq_in, seq_out
type(obs_type)          :: obs, prev_obs, next_obs, new_obs
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in, num_copies_in, num_qc_in
integer                 :: size_seq_out, num_copies_out, num_qc_out
integer                 :: num_inserted, iunit, io, i, j, total_num_inserted
integer                 :: max_num_obs, file_id, remaining_obs_count
integer                 :: first_seq
character(len = 129)    :: read_format, meta_data
logical                 :: pre_I_format, all_gone
logical                 :: trim_first, trim_last
character(len = 129)    :: msgstring

! could go into namelist if you wanted more control
integer, parameter      :: print_every = 20000

!----------------------------------------------------------------
! Namelist input with default values

! max_num_input_files : maximum number of input sequence files to be processed
integer, parameter :: max_num_input_files = 50
integer :: num_input_files = 1

! lazy, pick a big number
integer, parameter :: max_obs_input_types = 500
character(len = 32) :: obs_types(max_obs_input_types)
logical :: restrict_by_obs_type
integer :: num_obs_input_types
logical :: restrict_by_location, restrict_by_latlon
logical :: restrict_by_qc, restrict_by_copy, restrict_by_height


character(len = 129) :: filename_seq(max_num_input_files) = 'obs_seq.out'
character(len = 129) :: filename_out  = 'obs_seq.processed'
logical              :: process_file(max_num_input_files)

! Time of first and last observations to be used from obs_sequence
! If negative, these are not used
integer  :: first_obs_days    = -1
integer  :: first_obs_seconds = -1
integer  :: last_obs_days     = -1
integer  :: last_obs_seconds  = -1
! must be location independent...
type(location_type) :: min_loc, max_loc
real(r8) :: min_box(4) = missing_r8, max_box(4) = missing_r8
real(r8) :: min_lat = -90.0_r8
real(r8) :: max_lat =  90.0_r8
real(r8) :: min_lon =   0.0_r8
real(r8) :: max_lon = 360.0_r8
type(time_type) :: first_obs_time, last_obs_time
real(r8) :: min_qc = missing_r8
real(r8) :: max_qc = missing_r8
real(r8) :: min_copy = missing_r8
real(r8) :: max_copy = missing_r8
character(len = 32) :: copy_type = ''
character(len=129) :: qc_metadata = ''
character(len=129) :: copy_metadata = ''
logical  :: keep_types = .true.
logical  :: print_only = .false.
logical  :: gregorian_cal = .true.
real(r8) :: min_gps_height = missing_r8


namelist /obs_sequence_tool_nml/ num_input_files, filename_seq, filename_out, &
         first_obs_days, first_obs_seconds, last_obs_days, last_obs_seconds, &
         obs_types, keep_types, min_box, max_box, print_only, &
         min_lat, max_lat, min_lon, max_lon, min_qc, max_qc, qc_metadata, &
         min_copy, max_copy, copy_metadata, copy_type, gregorian_cal, &
         min_gps_height

!----------------------------------------------------------------
! Start of the routine.
! This routine basically opens the second observation sequence file
! and 'inserts' itself into the first observation sequence file.
! Each observation in the second file is independently inserted 
! or appended into the first file.
! Inserting takes time, appending is much faster when possible.
!----------------------------------------------------------------

call obs_seq_modules_used()

! if you are not using a gregorian cal, set this to false
! if anyone cares, we can add calendar integer to list
if (gregorian_cal) then
   call set_calendar_type(GREGORIAN)
endif

! Initialize input obs_seq filenames and obs_types
do i = 1, max_num_input_files
   filename_seq(i) = ""
enddo
do i = 1, max_obs_input_types
   obs_types(i) = ""
enddo

! Read the namelist entry
call find_namelist_in_file("input.nml", "obs_sequence_tool_nml", iunit)
read(iunit, nml = obs_sequence_tool_nml, iostat = io)
call check_namelist_read(iunit, io, "obs_sequence_tool_nml")

if(num_input_files .gt. max_num_input_files) then
  call error_handler(E_ERR,'obs_sequence_tool', &
     'num_input_files > max_num_input_files. change max_num_input_files in source file', &
     source,revision,revdate)
endif

! Record the namelist values used for the run ...
write(nmlfileunit, nml=obs_sequence_tool_nml)

! See if the user is restricting the obs types to be processed, and set up
! the values if so.
num_obs_input_types = 0
do i = 1, max_obs_input_types
   if ( len(obs_types(i)) .eq. 0 .or. obs_types(i) .eq. "" ) exit
   num_obs_input_types = i
enddo
if (num_obs_input_types == 0) then
   restrict_by_obs_type = .false.
else
   restrict_by_obs_type = .true.
endif

! See if the user is restricting the obs locations to be processed, and set up
! the values if so.
!!if ((minval(min_box).ne.missing_r8) .or. (maxval(min_box).ne.missing_r8) .or. &
!!    (minval(max_box).ne.missing_r8) .or. (maxval(max_box).ne.missing_r8)) then
!!   restrict_by_location = .true.
!!   restrict_by_latlon = .false.
!!   min_loc = set_location2(min_box)
!!   max_loc = set_location2(max_box)
!!else 
if ((min_lat /= -90.0_r8) .or. (max_lat /=  90.0_r8) .or. &
    (min_lon /=   0.0_r8) .or. (max_lon /= 360.0_r8)) then
   ! 3d sphere box check - locations module dependent, but an important one.
   restrict_by_location = .true.
   restrict_by_latlon = .true.
   if (trim(LocationName) /= 'loc3Dsphere') then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'can only use lat/lon box with 3d sphere locations', &
                         source,revision,revdate)
   endif
   ! simple err checks before going on; try to catch radians vs degrees or
   ! just plain junk or typos.
   if (min_lat >= max_lat) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lat must be less than max_lat', &
                         source,revision,revdate)
   endif
   if (min_lat < -90.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lat cannot be less than -90.0 degrees', &
                         source,revision,revdate)
   endif
   if (max_lat >  90.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'max_lat cannot be greater than  90.0 degrees', &
                         source,revision,revdate)
   endif
   ! this is ok - e.g.-180 to 180
   !if (min_lon < 0.0_r8) then
   !   call error_handler(E_ERR,'obs_sequence_tool', &
   !                      'min_lon cannot be less than 0.0 degrees', &
   !                      source,revision,revdate)
   !endif
   if (min_lon > 360.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lon cannot be greater than 360.0 degrees', &
                         source,revision,revdate)
   endif
   if (max_lon > 360.0_r8) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'max_lon cannot be greater than 360.0 degrees', &
                         source,revision,revdate)
   endif

   ! handle wrap in lon; at the end of this block the min is [0,360) and
   ! max is [0, 720) and greater than min.  modulo() must be used here and not
   ! mod() -- the result of modulo() is positive even if the input is negative.
   min_lon = modulo(min_lon,360.0_r8)
   max_lon = modulo(max_lon,360.0_r8)
   if (min_lon > max_lon) max_lon = max_lon + 360.0_r8

   if (min_lon == max_lon) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'min_lon cannot equal max_lon', &
                         source,revision,revdate)
   endif
else
   restrict_by_location = .false.
endif
  
!%! ! SPECIAL: cut off all GPS obs below the given height
!%! if (min_gps_height /= missing_r8) then
!%!    restrict_by_height = .true.
!%! else
   restrict_by_height = .false.
!%! endif


! See if the user is restricting the data values or qc to be processed, 
! and if so, set up
if ((min_qc /= missing_r8) .or. (max_qc /= missing_r8)) then
   if (len(trim(qc_metadata)) == 0) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'must specify the metadata name of a QC field', &
                         source,revision,revdate)
   endif
   restrict_by_qc = .true.
else
   restrict_by_qc = .false.
endif

if ((min_copy /= missing_r8) .or. (max_copy /= missing_r8)) then
   if (len(trim(copy_metadata)) == 0) then
      call error_handler(E_ERR,'obs_sequence_tool', &
                         'must specify the metadata name of a copy field', &
                         source,revision,revdate)
   endif
   restrict_by_copy = .true.
else
   restrict_by_copy = .false.
endif

! Read header information for the sequences to see if we need
! to accomodate additional copies or qc values from subsequent sequences.
! Also, calculate how many observations to be added to the first sequence.
num_copies_out   = 0
num_qc_out       = 0
size_seq_out     = 0

! check to see if we are going to trim the sequence by time
if(first_obs_seconds >= 0 .or. first_obs_days >= 0) then
   first_obs_time = set_time(first_obs_seconds, first_obs_days)
   trim_first = .true.
else
   trim_first = .false.
endif
if(last_obs_seconds >= 0 .or. last_obs_days >= 0) then
   last_obs_time = set_time(last_obs_seconds, last_obs_days)
   trim_last = .true.
else
   trim_last = .false.
endif
if (trim_first .and. trim_last) then
   if (first_obs_time > last_obs_time) then
      call error_handler(E_ERR,'obs_sequence_tool', 'first time cannot be later than last time', &
                         source,revision,revdate)
   endif
endif

! TWO PASS algorithm; open each file, trim it if requested, and count
! the number of actual observations.  then the output file can be
! created with the correct size, and as observations are put into it
! they'll be sorted, and unused obs will be removed.

! pass 1:

first_seq = -1
do i = 1, num_input_files

   if ( len(filename_seq(i)) .eq. 0 .or. filename_seq(i) .eq. "" ) then
      call error_handler(E_ERR,'obs_sequence_tool', &
         'num_input_files and filename_seq mismatch',source,revision,revdate)
   endif

   ! count up the number of observations we are going to eventually have.
   ! if all the observations in a file are not part of the linked list, the
   ! output number of observations might be much smaller than the total size in 
   ! the header.  it is slower, but go ahead and read in the entire sequence
   ! and count up the real number of obs - trim_seq will do the count even if
   ! it is not trimming in time.  this allows us to create an empty obs_seq
   ! output file of exactly the right size.

   call read_obs_seq_header(filename_seq(i), num_copies_in, num_qc_in, &
      size_seq_in, max_num_obs, file_id, read_format, pre_I_format, &
      close_the_file = .true.)
   
   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time, &
                 filename_seq(i), .true., remaining_obs_count)
   call destroy_obs_sequence(seq_in)
   if (remaining_obs_count == 0) then
      process_file(i) = .false.
      cycle
   else
      process_file(i) = .true.
      size_seq_in = remaining_obs_count
   endif

   if ( first_seq < 0 ) then
      first_seq = i
      num_copies_out = num_copies_in
      num_qc_out = num_qc_in
      size_seq_out = size_seq_in
   else
      size_seq_out = size_seq_out + size_seq_in
   endif

enddo

! no valid obs found?  if the index value is still negative, we are
! still waiting to process the first one and never found one.
if (first_seq < 0 .or. size_seq_out == 0) then
   msgstring = 'All input files are empty or all obs excluded by time/type/location'
   call error_handler(E_ERR,'obs_sequence_tool',msgstring,source,revision,revdate)
endif

! pass 2:

! Initialize individual observation variables 
call init_obs(     obs, num_copies_out, num_qc_out)
call init_obs( new_obs, num_copies_out, num_qc_out)
call init_obs(next_obs, num_copies_out, num_qc_out)
call init_obs(prev_obs, num_copies_out, num_qc_out)

total_num_inserted = 0

! Read obs seq to be added, and insert obs from it to the output seq
first_seq = -1
do i = 1, num_input_files

   if (.not. process_file(i)) cycle
 
   write(msgstring,*) 'Starting to process input sequence file ', trim(filename_seq(i))
   call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)

   call read_obs_seq(filename_seq(i), 0, 0, 0, seq_in)

   ! If you get here, there better be observations in this file which
   ! are going to be used (the process_file flag wouldn't be set otherwise.)
   call trim_seq(seq_in, trim_first, first_obs_time, trim_last, last_obs_time, &
                 filename_seq(i), .false., remaining_obs_count)

   ! This would be an error at this point.
   if(remaining_obs_count == 0) then
      call destroy_obs_sequence(seq_in) 
      write(msgstring, *) 'Internal error trying to process file ', trim(filename_seq(i))
      call error_handler(E_ERR,'obs_sequence_tool',msgstring,source,revision,revdate)
   endif

   ! create the output sequence here based on the first input file
   if (first_seq < 0) then
      call init_obs_sequence(seq_out, num_copies_out, num_qc_out, size_seq_out) 
      do j=1, num_copies_out
	 meta_data = get_copy_meta_data(seq_in, j) 
         call set_copy_meta_data(seq_out, j, meta_data)
      enddo 
      do j=1, num_qc_out
	 meta_data = get_qc_meta_data(seq_in, j) 
         call set_qc_meta_data(seq_out, j, meta_data)
      enddo 
      first_seq = i
   else
      ! we have an existing output sequence already.  make sure the next one
      ! is completely compatible.

      ! Compare metadata between the observation sequences.
      ! This routine exits if they do not match.
      call compare_metadata(seq_out, seq_in, filename_seq(first_seq), filename_seq(i))
   endif

   size_seq_out = get_num_key_range(seq_out)   !current size of seq_out
   size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

   if (print_only) then
      call print_obs_seq(seq_in, filename_seq(i))
   endif

   !-------------------------------------------------------------
   ! Start to insert obs from sequence_in into sequence_out
   !
   ! NOTE: insert_obs_in_seq CHANGES the obs passed in.
   !       Must pass a copy of incoming obs to insert_obs_in_seq.
   !--------------------------------------------------------------
   num_inserted = 0
   is_there_one = get_first_obs(seq_in, obs)

   if ( is_there_one )  then

      new_obs      = obs           ! obs records position in seq_out

!#!      call change_variance(new_obs)

      call insert_obs_in_seq(seq_out, new_obs)  ! new_obs linked list info changes

      prev_obs     = new_obs       ! records new position in seq_in
      num_inserted = num_inserted + 1
   
      call get_next_obs(seq_in, obs, next_obs, is_this_last)
      ObsLoop : do while ( .not. is_this_last)

         if (print_every > 0) then
            if (mod(num_inserted,print_every) == 0) then
               print*, 'inserted number ',num_inserted,' of ',size_seq_in
            endif
         endif

         obs     = next_obs   ! essentially records position in seq_out
         new_obs = obs        ! will be modified w/ position in seq_in

!#!         call change_variance(new_obs)

         ! Since the stride through the observation sequence file is always 
         ! guaranteed to be in temporally-ascending order, we can use the
         ! 'previous' observation as the starting point to search for the
         ! correct insertion point.  This speeds up the insert code a lot.

         call insert_obs_in_seq(seq_out, new_obs, prev_obs)

         prev_obs     = new_obs    ! update position in seq_in for next insert
         num_inserted = num_inserted + 1

         call get_next_obs(seq_in, obs, next_obs, is_this_last)

      enddo ObsLoop

      total_num_inserted = total_num_inserted + num_inserted

   else
      write(msgstring,*)'no first observation in ',trim(adjustl(filename_seq(i)))
      call error_handler(E_MSG,'obs_sequence_tool', msgstring, source, revision, revdate)
   endif

   if (.not. print_only) then
      print*, '--------------  Obs seq file # :          ', i
      print*, 'Number of obs in previous seq  :          ', size_seq_out
      print*, 'Number of obs to be  inserted  :          ', size_seq_in
      print*, 'Number of obs really inserted  :          ', num_inserted
      print*, '---------------------------------------------------------'
   endif

   call destroy_obs_sequence(seq_in)

enddo

write(msgstring,*) 'Starting to process output sequence file ', trim(filename_out)
call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)

if (.not. print_only) then
   print*, 'Total number of obs  inserted  :          ', total_num_inserted
   print*, 'Actual number of obs in the new seq file :', get_num_key_range(seq_out)
else
   print*, 'Total number of selected obs in all files :', get_num_key_range(seq_out)
endif

call print_obs_seq(seq_out, filename_out)
if (.not. print_only) then
   call write_obs_seq(seq_out, filename_out)
else
   write(msgstring,*) 'Output sequence file not created; print_only in namelist is .true.'
   print *, trim(msgstring)
   call error_handler(E_MSG,'obs_sequence_tool', msgstring, source, revision, revdate)
endif

! Time to clean up

call destroy_obs_sequence(seq_out)
call destroy_obs(     obs)
call destroy_obs( new_obs)
call destroy_obs(next_obs)
call destroy_obs(prev_obs)

call timestamp(source,revision,revdate,'end')

contains


!---------------------------------------------------------------------
subroutine obs_seq_modules_used()

! Initialize modules used that require it
call initialize_utilities('obs_sequence_tool')
call register_module(source,revision,revdate)
call static_init_obs_sequence()

end subroutine obs_seq_modules_used


!---------------------------------------------------------------------
subroutine compare_metadata(seq1, seq2, fname1, fname2)

!
! This subroutine compares the metadata for two different observation
! sequences and terminates the program if they are not conformable.
! In order to be merged, the two observation sequences must have the same
! number of qc values, the same number of copies ... 
!
! The messages might be a bit misleading 'warning', 'warning', 'error' ...
!

 type(obs_sequence_type), intent(IN) :: seq1, seq2
 character(len=*), optional :: fname1, fname2

integer :: num_copies1, num_qc1
integer :: num_copies2, num_qc2
integer :: num_copies , num_qc, i
character(len=129) :: str1, str2
character(len=255) :: msgstring1, msgstring2

num_qc1     = get_num_qc(    seq1)
num_copies1 = get_num_copies(seq1)

num_qc2     = get_num_qc(    seq2)
num_copies2 = get_num_copies(seq2)

num_copies  = num_copies1
num_qc      = num_qc1

! get this ready in case we have to use it
if (present(fname1) .and. present(fname2)) then
   write(msgstring1,*)'Sequence files ', trim(fname1), ' and ', trim(fname2), &
                      ' are not compatible'
else
  msgstring1 = 'Sequence files cannot be merged because they are not compatible'
endif

if ( num_copies1 /= num_copies2 ) then
   write(msgstring2,*)'Different numbers of data copies found: ', &
                      num_copies1, ' vs ', num_copies2 
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
   num_copies = -1
endif
if ( num_qc1 /= num_qc2 ) then
   write(msgstring2,*)'Different different numbers of QCs found: ', &
                      num_qc1, ' vs ', num_qc2
   call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
   num_qc = -1
endif
if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)
endif

MetaDataLoop : do i=1, num_copies
   str1 = trim(adjustl(get_copy_meta_data(seq1,i)))
   str2 = trim(adjustl(get_copy_meta_data(seq2,i)))

   if( str1 == str2 ) then
      write(msgstring2,*)'metadata ',trim(adjustl(str1)), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
   else
      write(msgstring2,*)'metadata value mismatch. seq1: ', trim(adjustl(str1))
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
      write(msgstring2,*)'metadata value mismatch. seq2: ', trim(adjustl(str2))
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
      call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)
   endif
enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = trim(adjustl(get_qc_meta_data(seq1,i)))
   str2 = trim(adjustl(get_qc_meta_data(seq2,i)))

   if( str1 == str2 ) then
      write(msgstring2,*)'qc metadata ', trim(adjustl(str1)), ' in both.'
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
   else
      write(msgstring2,*)'qc metadata value mismatch. seq1: ', trim(adjustl(str1))
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
      write(msgstring2,*)'qc metadata value mismatch. seq2: ', trim(adjustl(str2))
      call error_handler(E_MSG, 'obs_sequence_tool', msgstring2, source, revision, revdate)
      call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)
   endif
enddo QCMetaData

end subroutine compare_metadata

!---------------------------------------------------------------------
! pass in an already opened sequence and a start/end time.  this routine
! really trims the observations out of the sequence, and returns a count
! of how many remain.
subroutine trim_seq(seq, trim_first, first_time, trim_last, last_time, &
                    seqfilename, print_msg, remaining_obs_count)
 type(obs_sequence_type), intent(inout) :: seq
 logical, intent(in)                    :: trim_first, trim_last
 type(time_type), intent(in)            :: first_time, last_time
 character(len = *), intent(in)         :: seqfilename
 logical, intent(in)                    :: print_msg
 integer, intent(out)                   :: remaining_obs_count

 integer :: i
 logical :: found

   ! Need to find first obs with appropriate time, delete all earlier ones
   if(trim_first) then
      call delete_seq_head(first_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are before first_obs_days:first_obs_seconds'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   
   ! Also get rid of observations past the last_obs_time if requested
   if(trim_last) then
      call delete_seq_tail(last_time, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are after last_obs_days:last_obs_seconds'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   
   ! optionally select only a list of obs types
   if (restrict_by_obs_type) then
      call delete_obs_by_typelist(num_obs_input_types, obs_types, &
                                  keep_types, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            if (keep_types) then
               msgstring = 'Skipping: no obs in ' // trim(seqfilename) // &
                           ' are on the keep obs_types list'
            else
               msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                           ' are on the discard obs_types list'
            endif
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! optionally select obs in a lat/lon box
   ! the more generic is a box in any number of dimensions which matches
   ! the locations module which is determined at link time.  but for now...
   if (restrict_by_location) then
      if (restrict_by_latlon) then
         call select_obs_by_latlon(min_lon, max_lon, min_lat, max_lat, &
                                   seq, all_gone)
      else 
         !%! call select_obs_by_location(min_loc, max_loc, seq, all_gone)
         msgstring = 'Select by general bounding box not implemented yet'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: no obs in ' // trim(seqfilename) // &
                        ' are above the given height'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! optionally select only a range of qc and copies
   if (restrict_by_qc) then
      ! validate the metadata string
      found = .false.
      do i=1, num_qc_in
	 meta_data = get_qc_meta_data(seq_in, i) 
         if (trim(qc_metadata) == trim(meta_data)) then 
            call delete_obs_by_qc(i, min_qc, max_qc, seq, all_gone)
            found = .true.
            exit
         endif
      enddo 
      if (.not. found) then 
         msgstring = 'QC metadata string: ' // trim(qc_metadata) // &
                     ' not found in obs_seq file'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are outside the qc min/max range'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif
   if (restrict_by_copy) then
      ! validate the metadata string
      found = .false.
      do i=1, num_copies_in
	 meta_data = get_copy_meta_data(seq_in, i) 
         if (trim(copy_metadata) == trim(meta_data)) then 
            call delete_obs_by_copy(i, min_copy, max_copy, copy_type, &
                                    seq, all_gone)
            found = .true.
            exit
         endif
      enddo 
      if (.not. found) then 
         msgstring = 'Copy metadata string: ' // trim(copy_metadata) // &
                     ' not found in obs_seq file'
         call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
      endif
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: all obs in ' // trim(seqfilename) // &
                        ' are outside the copy min/max range'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   ! SPECIAL: optionally restrict GPS obs to above a height
   if (restrict_by_height) then
      call select_gps_by_height(min_gps_height, seq, all_gone)
      if(all_gone) then
         if (print_msg) then
            msgstring = 'Skipping: no obs in ' // trim(seqfilename) // &
                        ' are above the GPS height threshold'
            call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
         endif
         remaining_obs_count = 0
         return
      endif
   endif

   remaining_obs_count = get_num_key_range(seq)

end subroutine trim_seq


!---------------------------------------------------------------------
subroutine print_obs_seq(seq_in, filename)

! you can get more info by running the obs_diag program, but this
! prints out a quick table of obs types and counts, overall start and
! stop times, and metadata strings and counts.

type(obs_sequence_type), intent(in) :: seq_in
character(len=*), intent(in)        :: filename

type(obs_type)          :: obs, next_obs
type(obs_def_type)      :: this_obs_def
logical                 :: is_there_one, is_this_last
integer                 :: size_seq_in, num_copies_in, num_qc_in
integer                 :: i, j
integer                 :: this_obs_kind
! max_obs_kinds is a public from obs_kind_mod.f90 and really is
! counting the max number of types, not kinds
integer                 :: type_count(max_obs_kinds), identity_count


! Initialize input obs_types
do i = 1, max_obs_kinds
   type_count(i) = 0
enddo
identity_count = 0

! make sure there are obs left to process before going on.
! num_obs should be ok since we just constructed this seq so it should
! have no unlinked obs.  if it might for some reason, use this instead:
! size_seq_in = get_num_key_range(seq_in)     !current size of seq_in

size_seq_in = get_num_obs(seq_in)
if (size_seq_in == 0) then
   msgstring = 'Obs_seq file is empty.'
   call error_handler(E_MSG,'obs_sequence_tool',msgstring,source,revision,revdate)
   return
endif

! Initialize individual observation variables 
call init_obs(     obs, get_num_copies(seq_in), get_num_qc(seq_in))
call init_obs(next_obs, get_num_copies(seq_in), get_num_qc(seq_in))

! blank line
call error_handler(E_MSG,'',' ')

write(msgstring,*) 'Processing sequence file ', trim(filename)
call error_handler(E_MSG,'',msgstring)

call print_metadata(seq_in, filename)

!-------------------------------------------------------------
! Start to process obs from seq_in
!--------------------------------------------------------------
is_there_one = get_first_obs(seq_in, obs)

if ( .not. is_there_one )  then
   write(msgstring,*)'no first observation in ',trim(adjustl(filename))
   call error_handler(E_MSG,'obs_sequence_tool', msgstring, source, revision, revdate)
endif

! process it here
is_this_last = .false.

call get_obs_def(obs, this_obs_def)
call print_time(get_obs_def_time(this_obs_def), ' First timestamp: ')
if (gregorian_cal) then
   call print_date(get_obs_def_time(this_obs_def), '   Gregorian day: ')
endif

ObsLoop : do while ( .not. is_this_last)

   call get_obs_def(obs, this_obs_def)
   this_obs_kind = get_obs_kind(this_obs_def)
   if (this_obs_kind < 0) then
      identity_count = identity_count + 1
   else
      type_count(this_obs_kind) = type_count(this_obs_kind) + 1
   endif
!   print *, 'obs kind index = ', this_obs_kind
!   if(this_obs_kind > 0)print *, 'obs name = ', get_obs_kind_name(this_obs_kind)

   call get_next_obs(seq_in, obs, next_obs, is_this_last)
   if (.not. is_this_last) then 
      obs = next_obs
   else
      call print_time(get_obs_def_time(this_obs_def), '  Last timestamp: ')
      if (gregorian_cal) then
         call print_date(get_obs_def_time(this_obs_def), '   Gregorian day: ')
      endif
   endif

enddo ObsLoop


write(msgstring, *) 'Number of obs processed  :          ', size_seq_in
call error_handler(E_MSG, '', msgstring)
write(msgstring, *) '---------------------------------------------------------'
call error_handler(E_MSG, '', msgstring)
do i = 1, max_obs_kinds
   if (type_count(i) > 0) then 
      write(msgstring, '(a32,i8,a)') trim(get_obs_kind_name(i)), &
                                     type_count(i), ' obs'
      call error_handler(E_MSG, '', msgstring)
   endif
enddo
if (identity_count > 0) then 
   write(msgstring, '(a32,i8,a)') 'Identity observations', &
                                  identity_count, ' obs'
   call error_handler(E_MSG, '', msgstring)
endif

! another blank line
call error_handler(E_MSG, '', ' ')

! Time to clean up

call destroy_obs(     obs)
call destroy_obs(next_obs)

end subroutine print_obs_seq

!---------------------------------------------------------------------
subroutine print_metadata(seq1, fname1)

!
! print out the metadata strings, trimmed
!

type(obs_sequence_type), intent(in) :: seq1
character(len=*), optional :: fname1

integer :: num_copies , num_qc, i
character(len=129) :: str1
character(len=255) :: msgstring1

num_copies = get_num_copies(seq1)
num_qc     = get_num_qc(    seq1)

if ( num_copies < 0 .or. num_qc < 0 ) then
   call error_handler(E_ERR, 'obs_sequence_tool', msgstring1, source, revision, revdate)
endif

MetaDataLoop : do i=1, num_copies
   str1 = trim(adjustl(get_copy_meta_data(seq1,i)))

   write(msgstring1,*)'Data Metadata: ',trim(adjustl(str1))
   call error_handler(E_MSG, '', msgstring1)

enddo MetaDataLoop

QCMetaData : do i=1, num_qc
   str1 = trim(adjustl(get_qc_meta_data(seq1,i)))

   write(msgstring1,*)'  QC Metadata: ', trim(adjustl(str1))
   call error_handler(E_MSG, '', msgstring1)
   
enddo QCMetaData

end subroutine print_metadata

!---------------------------------------------------------------------
!---------------------------------------------------------------------
! WARNING:
! this routine should be in obs_sequence_mod.f90 but is has one line which
! is dependent on the locations module being the 3d sphere; if it is
! something else, this code will not work.  i have a first cut at a loc-indep
! routine but not all the kinks are worked out yet so put this code here for
! now.  to use merge with any other locations module, comment out
! the calls to get_location() and do not select by bounding box in the nml.

subroutine select_obs_by_latlon(min_lon, max_lon, min_lat, max_lat, &
                                seq, all_gone)

! Delete all observations in the sequence which are outside the bounding box.
! If there are no obs left afterwards return that the sequence is all_gone.

real(r8),                intent(in)    :: min_lon, max_lon, min_lat, max_lat
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

type(obs_def_type)   :: obs_def
type(obs_type)       :: obs, prev_obs
integer              :: i, key
type(location_type)  :: location
logical              :: out_of_range, is_this_last, inside, first_obs
real(r8)             :: ll(3), lat, lon


! Initialize an observation type with appropriate size
call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))

! Iterate entire sequence, deleting obs which are not in the box.
! First, make sure there are obs to delete, and initialize first obs.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
   call destroy_obs(obs)
   call destroy_obs(prev_obs)
   return
endif

first_obs = .true.
prev_obs = obs

! This is going to be O(n), n=num obs in seq
is_this_last = .false.
allobs : do while (.not. is_this_last)

   call get_obs_def(obs, obs_def)
   location = get_obs_def_location(obs_def)

   ! each diff locations mod has a different one of these
   ! FIXME: this next line is what makes this locations dependent.
   ! this also is ignoring the vertical for now.
   ll = get_location(location)
   lon = ll(1)
   lat = ll(2)
   
   ! if wrap in longitude, lon now (0,720] 
   !if (lon < max_lon-360.0_r8) lon = lon + 360.0
   if (max_lon >= 360.0) lon = lon + 360.0

   ! box test.
   if ((lon < min_lon) .or. (lon > max_lon) .or. &
       (lat < min_lat) .or. (lat > max_lat)) then
      inside = .false.
   else
      inside = .true.
   endif
   
   ! same code as delete/keep by obstype; do any code fixes both places
   if (.not. inside) then
      if (first_obs) then
         call delete_obs_from_seq(seq, obs)
         if(.not. get_first_obs(seq, obs)) exit allobs
      else
!print *, 'going to del obs key ', obs%key
!print *, 'prev key is ', prev_obs%key
         call delete_obs_from_seq(seq, obs)
         ! cannot simply use prev_obs; cached copy out of sync with seq one
         key = get_obs_key(prev_obs)
         call get_next_obs_from_key(seq, key, obs, is_this_last)
!print *, 'next obs now is key ', obs%key
      endif
   else
!print *, 'no del, keep this obs key ', obs%key
      first_obs = .false.
     prev_obs = obs
!print *, 'prev obs now is key ', prev_obs%key
!print *, 'obs was key ', obs%key
      call get_next_obs(seq, prev_obs, obs, is_this_last)
!print *, 'obs now is key ', obs%key
   endif
   
end do allobs

! Figure out if there are no more obs left in the sequence.
if(.not. get_first_obs(seq, obs)) then
   all_gone = .true.
else
   all_gone = .false.
endif

! Done.  delete temp storage and return.
call destroy_obs(obs)
call destroy_obs(prev_obs)

end subroutine select_obs_by_latlon


!---------------------------------------------------------------------
subroutine select_gps_by_height(min_height, seq, all_gone)

! CURRENTLY COMMENTED OUT

! Delete all gps observations in the sequence which are below the given ht.
! If there are no obs left afterwards return that the sequence is all_gone.

real(r8),                intent(in)    :: min_height
type(obs_sequence_type), intent(inout) :: seq
logical,                 intent(out)   :: all_gone

all_gone = .false.
return

!%! type(obs_def_type)   :: obs_def
!%! type(obs_type)       :: obs, prev_obs
!%! integer              :: i, key, gps_type_index, this_obs_type
!%! type(location_type)  :: location
!%! logical              :: out_of_range, is_this_last, above, first_obs
!%! real(r8)             :: ll(3), vloc
!%! 
!%! ! figure out what index number is gps
!%! gps_type_index = get_obs_kind_index('GPSRO_REFRACTIVITY')
!%! if (gps_type_index < 0) then
!%!    write(msgstring,*) 'obs_type GPSRO not found'
!%!    call error_handler(E_ERR,'select_gps_by_height', msgstring, &
!%!                       source, revision, revdate)
!%! endif
!%! 
!%! 
!%! ! Initialize an observation type with appropriate size
!%! call init_obs(obs, get_num_copies(seq), get_num_qc(seq))
!%! call init_obs(prev_obs, get_num_copies(seq), get_num_qc(seq))
!%! 
!%! ! Iterate entire sequence, deleting obs which are not in the box.
!%! ! First, make sure there are obs to delete, and initialize first obs.
!%! if(.not. get_first_obs(seq, obs)) then
!%!    all_gone = .true.
!%!    call destroy_obs(obs)
!%!    call destroy_obs(prev_obs)
!%!    return
!%! endif
!%! 
!%! first_obs = .true.
!%! prev_obs = obs
!%! 
!%! ! This is going to be O(n), n=num obs in seq
!%! is_this_last = .false.
!%! allobs : do while (.not. is_this_last)
!%! 
!%!    ! verify GPS obs first, then do height
!%! 
!%!    call get_obs_def(obs, obs_def)
!%!    this_obs_type = get_obs_kind(obs_def)
!%!    if (this_obs_type /= gps_type_index) then
!%!       ! we are going to keep this
!%!       above = .true.   
!%!    else
!%!    
!%!       ! must check height.  at this point, all gps obs are be height only. 
!%!       location = get_obs_def_location(obs_def)
!%!    
!%!       ! this makes the tool sphere_3d dependent.  also assumes height as vert.
!%!       ! should verify. 
!%!       if (.not. vert_is_height(location)) then
!%!          write(msgstring,*) 'obs_type GPSRO vertical location not height'
!%!          call error_handler(E_ERR,'select_gps_by_height', msgstring, &
!%!                             source, revision, revdate)
!%!       endif
!%! 
!%!       ll = get_location(location)
!%!       vloc = ll(3)
!%! 
!%!       if (vloc < min_height) then
!%!          ! delete this one
!%!          above = .false.
!%!       else
!%!          ! will be kept
!%!          above = .true.
!%!       endif
!%!    endif 
!%! 
!%!    ! same code as delete/keep by obstype; do any code fixes both places
!%!    if (.not. above) then
!%!       if (first_obs) then
!%!          call delete_obs_from_seq(seq, obs)
!%!          if(.not. get_first_obs(seq, obs)) exit allobs
!%!       else
!%! !print *, 'going to del obs key ', obs%key
!%! !print *, 'prev key is ', prev_obs%key
!%!          call delete_obs_from_seq(seq, obs)
!%!          ! cannot simply use prev_obs; cached copy out of sync with seq one
!%!          key = get_obs_key(prev_obs)
!%!          call get_next_obs_from_key(seq, key, obs, is_this_last)
!%! !print *, 'next obs now is key ', obs%key
!%!       endif
!%!    else
!%! !print *, 'no del, keep this obs key ', obs%key
!%!       first_obs = .false.
!%!       prev_obs = obs
!%! !print *, 'prev obs now is key ', prev_obs%key
!%! !print *, 'obs was key ', obs%key
!%!       call get_next_obs(seq, prev_obs, obs, is_this_last)
!%! !print *, 'obs now is key ', obs%key
!%!    endif
!%!    
!%! end do allobs
!%! 
!%! ! Figure out if there are no more obs left in the sequence.
!%! if(.not. get_first_obs(seq, obs)) then
!%!    all_gone = .true.
!%! else
!%!    all_gone = .false.
!%! endif
!%! 
!%! ! Done.  delete temp storage and return.
!%! call destroy_obs(obs)
!%! call destroy_obs(prev_obs)

end subroutine select_gps_by_height


!---------------------------------------------------------------------
!#! subroutine change_variance(this_obs)
!#! 
!#! 
!#! use obs_kind_mod
!#! use obs_def_mod
!#! use obs_sequence_mod
!#! 
!#! ! change the variance on specific kinds.
!#! 
!#! type(obs_type), intent(inout) :: this_obs
!#! 
!#! type(obs_def_type) :: this_obs_def
!#! integer            :: this_obs_kind
!#! real(r8)           :: oldvar, newvar
!#! 
!#! 
!#! call get_obs_def(obs, this_obs_def)
!#! this_obs_kind = get_obs_kind(this_obs_def)
!#! 
!#! ! ignore identity obs
!#! if (this_obs_kind < 0) return
!#! 
!#! oldvar = get_obs_def_error_variance(this_obs_def)
!#! 
!#! !print *, 'kind, var = ', this_obs_kind, oldvar
!#! ! SOYOUNG: Change the code here for what you want.
!#! 
!#! ! Set the new variance here based on the type
!#! select case (this_obs_kind)
!#!   case (LAND_SFC_TEMPERATURE)
!#!      newvar = oldvar * 0.5
!#!   case (RADIOSONDE_TEMPERATURE)
!#!      newvar = 1.3
!#!   case (GPSRO_REFRACTIVITY)
!#!      newvar = 1.666
!#!   ! etc
!#!   case default
!#!      newvar = oldvar
!#! end select
!#! 
!#! !print *, 'newvar = ', newvar
!#! if (newvar /= oldvar) then
!#!    call set_obs_def_error_variance(this_obs_def, newvar)
!#!    call set_obs_def(this_obs, this_obs_def)
!#! endif
!#! 
!#! end subroutine change_variance

!---------------------------------------------------------------------
end program obs_sequence_tool