! place problem-specific variables that will go into the probdata namelist
! create a local copy probdata.f90 in your build directorymodule probdata_module

module probdata_module

  implicit none

  private

  ! variables for the probdata namelist
  double precision, save, public :: pres_base = 1.65d6
  double precision, save, public :: dens_base = 1.65d-3
  double precision, save, public :: pert_factor = 8.1d-3
  double precision, save, public :: y_pert_center = 0.7d0
  double precision, save, public :: pert_width = 0.025d0
  logical, save, public :: single = .false.
  logical, save, public :: do_isentropic = .true.

end module probdata_module

subroutine probdata_init(name,namlen)

  use probdata_module

  implicit none

  integer :: namlen
  integer :: name(namlen)

  integer :: un, i, status

  integer, parameter :: maxlen = 256
  character (len=maxlen) :: probin

  namelist /probdata/ pres_base
  namelist /probdata/ dens_base
  namelist /probdata/ pert_factor
  namelist /probdata/ y_pert_center
  namelist /probdata/ pert_width
  namelist /probdata/ single
  namelist /probdata/ do_isentropic

  ! default values
  pres_base = 1.65d6
  dens_base = 1.65d-3
  pert_factor = 0.1d0
  y_pert_center = 0.5d0
  pert_width = 0.1d0
  single = .false.
  do_isentropic = .true.

  ! create the filename
  if (namlen > maxlen) then
     print *, 'probin file name too long'
     stop
  endif

  do i = 1, namlen
     probin(i:i) = char(name(i))
  end do

  ! read in the namelist
  un = 9
  open (unit=un, file=probin(1:namlen), form='formatted', status='old')
  read (unit=un, nml=probdata, iostat=status)

  if (status < 0) then
     ! the namelist does not exist, so we just go with the defaults
     continue

  else if (status > 0) then
     ! some problem in the namelist
     print *, 'ERROR: problem in the probdata namelist'
     stop
  endif

  close (unit=un)

  !$acc update &
  !$acc device(pert_temp_factor, pert_rad_factor, do_small_domain)

end subroutine probdata_init
