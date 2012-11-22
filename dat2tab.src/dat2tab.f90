!!$--------!---------!---------!---------!---------!---------!---------%
!!$
!!$     File Name: DAT2TAB.F90
!!$
!!$     Author:    Carlos Silva Santos, Megajoule Inovação Lda.
!!$
!!$     URL:       www.megajoule.pt
!!$
!!$     Version:   2.0
!!$
!!$     Date:      21.October.2012
!!$
!!$--------!---------!---------!---------!---------!---------!---------%
program dat2tab

  implicit none

  character(len=100) :: filename
  integer :: n,ifreq,nfreq,idir,ndir,ibin,nbin,imeas,nmeas,iseed
  integer :: nmeas_filt,imeas_valid,nmeas_valid,nmeas_filt_vel
  integer :: istat
  real,allocatable,dimension(:) :: fdir,vhm,dirm,vhm_filt,dirm_filt,vhm_valid,dirm_valid
  real,allocatable,dimension(:) :: vh_rubbish,dir_rubbish
  real,allocatable,dimension(:,:) :: fbin_dir,fbin_glo,aux
  real, allocatable, dimension(:) :: lixo
  integer,allocatable,dimension(:) :: ones
  integer(kind=8),allocatable,dimension(:) :: ilixo
  real :: randomdir,randomvel,aperture,flow_dir
  logical, allocatable, dimension(:) :: lvalid

!!$ Get random seed
  integer(4) :: ic4, crate4, cmax4
  call system_clock(count=ic4, count_rate=crate4, count_max=cmax4)

  iseed = ic4

!!$ Test if at least an argument was provided
  n = iargc()  
  if (n.eq.0) then
     write(*,*) '!! ERROR !!'
     write(*,*) 'Usage : dat2tab <file.dat>'
     stop
  end if

  call getarg(1,filename)  
  filename = trim(filename)

!!$ Initialize variables
  ones = 1
  ilixo = 0
  
!!$ Reading DAT file

!!$ read first to find number of measurements for allocation
  open(1,file=filename)
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  nmeas = 1
  do
     read(1,*,iostat=istat)
     if(istat < 0 ) exit 
     nmeas = nmeas + 1
  end do
  close(1)

  nmeas = nmeas - 1

  allocate(dirm(nmeas),vhm(nmeas))
  allocate(dirm_filt(nmeas),vhm_filt(nmeas))
  allocate(dirm_valid(nmeas),vhm_valid(nmeas))
  allocate(ones(nmeas),lixo(nmeas),ilixo(nmeas))
  allocate(lvalid(nmeas))
  allocate(dir_rubbish(nmeas),vh_rubbish(nmeas))

  call read_mjseries(filename, &                      !!$ IN  : input file name
       nmeas, &                                       !!$ IN  : number of records
       dirm(1:nmeas),vhm(1:nmeas), &                  !!$ OUT : dir & vh of measured series
       lixo(1:nmeas),lixo(1:nmeas), &                 !!$ OUT : max & minimum velocities
       lixo(1:nmeas), &                               !!$ OUT : measured series of USTD
       lixo(1:nmeas),lixo(1:nmeas), &                 !!$ OUT : derived quantities u & v
       lixo(1:nmeas), &                               !!$ OUT : derived quantity Iu
       ilixo(1:nmeas))                                !!$ OUT : vector of time stamps

!!$ Filter known direction errors (dir > 360 & vh >= 99.0)
  lvalid = .true.
  where (dirm(1:nmeas).gt.360) lvalid = .false.
  where (vhm(1:nmeas).gt.99.0) lvalid = .false.

!!$ Count valid measurements
  nmeas_valid = count(lvalid)

!!$ Create vectors with valid measurements only
  imeas_valid = 0
  do imeas = 1, nmeas
     if (lvalid(imeas)) then
        imeas_valid = imeas_valid + 1
        vhm_valid(imeas_valid) = vhm(imeas)
        dirm_valid(imeas_valid) = dirm(imeas)
     end if
  end do

!!$ Replace measured data vectors by valid measured data. i.e. invalid measurements are destroyed
  vhm = 0.0
  dirm = 0.0
  nmeas_valid = imeas_valid
  nmeas = nmeas_valid
  vhm(1:nmeas) = vhm_valid(1:nmeas)
  dirm(1:nmeas) = dirm_valid(1:nmeas)

  ndir = 12
  nbin = max(int(maxval(vhm_valid(1:nmeas))),50)

  allocate(fdir(ndir),fbin_dir(ndir,nbin),fbin_glo(ndir,nbin),aux(ndir,nbin))
   
!!$ Filter and calculate frequencies per sector
  do idir = 1, ndir
     dirm_filt = 0.0
     vhm_filt = 0.0
     nmeas_filt = 0
     flow_dir = (idir - 1) * (360 / ndir)
     aperture = (360 / ndir)
     call filter_dir(nmeas, &                       !!$ IN  : number of elements in series
          dirm(1:nmeas), &                          !!$ IN  : dir vector : unfiltered
          vhm(1:nmeas), &                           !!$ IN  : vel vector : unfiltered
          flow_dir,aperture, &                      !!$ IN  : dir filter values
          dirm_filt,vhm_filt, &                     !!$ OUT : dir & vel vectors : FILTERED 
          nmeas_filt)                               !!$ OUT : length of filtered vectors
     
     if (nmeas_filt.ne.0) then
!!$ Calculate sector frequency  
        fdir(idir) = real(nmeas_filt) / real(nmeas_valid)
     else
!!$ Minimum values when no values exist 
        fdir(idir) = 0.0
     end if
     
     do ibin = 1, nbin
        call filter_vel(nmeas_filt, &               !!$ IN  : number of elements in series
             dirm_filt(1:nmeas_filt), &             !!$ IN  : dir vector : filtered by dir
             vhm_filt(1:nmeas_filt), &              !!$ IN  : vel vector : filtered by dir
             real(ibin)-1.0,real(ibin), &           !!$ IN  : min & max velocities
             dir_rubbish,vh_rubbish, &              !!$ OUT : dir & vel vectors : FILTERED 
             nmeas_filt_vel)                        !!$ OUT : length of filtered vectors
!!$ Calculate vbin/sector frequency
        fbin_dir(idir,ibin) = real(nmeas_filt_vel) / real(nmeas_filt)
     end do
  end do

!!$ Write artificial data series in MJ format
  filename = trim(filename)

  call write_tab(len(trim(filename) // '.tab'), &  !!$ IN  : length of filename
       trim(filename) // '.tab', &                 !!$ IN  : input file name
       ndir,nbin, &                                !!$ IN  : number of sector and velocity bins
       fdir(1:ndir), &                             !!$ IN  : sector frequency
       fbin_dir(1:ndir,1:nbin), &                  !!$ IN  : sector freq of velocity bins
       fbin_glo(1:ndir,1:nbin), &                  !!$ IN  : global freq of velocity bins
       len(trim(filename)),trim(filename))         !!$ IN  : label for TAB file

end program dat2tab
!!$ ############################################################################
!!$ ############################################################################
subroutine filter_dir(nmeas,dirms,vhms,dir,aperture,dirms_filt,vhms_filt,nmeas_filt)  

  implicit none

!!$ Passing variables
  integer :: nmeas,nmeas_filt,not_nmeas_filt
  real, dimension(nmeas) :: dirms,vhms,dirms_filt,vhms_filt
  real :: dir,aperture

!!$ Local variables
  integer :: imeas
  real :: dirm,dirp

!!$ Initialize variables
  vhms_filt = 0.0
  dirms_filt = 0.0
  nmeas_filt = 0
  not_nmeas_filt = 0

  dirm = dir - aperture / 2
  dirp = dir + aperture / 2
  if (dirm.lt.0) dirm = dirm + 360.0
  do imeas = 1, nmeas
     if (dirms(imeas).lt.360) then
        if (dirm.gt.dirp) then
           if (dirms(imeas).ge.dirm.or.dirms(imeas).lt.dirp) then
              nmeas_filt = nmeas_filt + 1
              vhms_filt(nmeas_filt) = vhms(imeas)
              dirms_filt(nmeas_filt) = dirms(imeas)     
           else
              not_nmeas_filt = not_nmeas_filt + 1
           end if
        else
           if (dirms(imeas).ge.dirm.and.dirms(imeas).lt.dirp) then
              nmeas_filt = nmeas_filt + 1
              vhms_filt(nmeas_filt) = vhms(imeas)
              dirms_filt(nmeas_filt) = dirms(imeas)     
           else
              not_nmeas_filt = not_nmeas_filt + 1
           end if
        end if
     end if
  end do

end subroutine filter_dir
!!$ ############################################################################
!!$ ############################################################################
subroutine filter_vel(nmeas,dirms,vhms,velm,velp,dirms_filt,vhms_filt,nmeas_filt)  

  implicit none
!!$ Passing variables
  integer :: nmeas,nmeas_filt
  real, dimension(nmeas) :: dirms,vhms,dirms_filt,vhms_filt
  real :: velm,velp

!!$ Local variables
  integer :: imeas
  real :: dirm,dirp

!!$ Initialize variables
  vhms_filt = 0.0
  dirms_filt = 0.0
  nmeas_filt = 0

  do imeas = 1, nmeas
     if (vhms(imeas).ge.velm.and.vhms(imeas).lt.velp) then
        nmeas_filt = nmeas_filt + 1
        vhms_filt(nmeas_filt) = vhms(imeas)
        dirms_filt(nmeas_filt) = dirms(imeas)     
     end if
  end do

end subroutine filter_vel
