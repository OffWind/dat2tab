!########################################################################
!########################################################################
subroutine read_mjseries(filename, &                       !!$ IN  : input file name
             nmeas, &                                      !!$ OUT : number of records
             dir,vh, &                                     !!$ OUT : dir & vh of measured series
             vmax,vmin, &                                  !!$ OUT : max & minimum velocities
             ustd, &                                       !!$ OUT : measured series of USTD
             u,v, &                                        !!$ OUT : derived quantities u & v
             Iu, &                                         !!$ OUT : derived quantity Iu
             time)                                         !!$ OUT : vector of time stamps

  implicit none

!!$ Interface variables
  integer :: nmeas
  character(len=100) :: filename
  character(len=1) :: rubbish
  real,dimension(nmeas) :: dir,vh,vmax,vmin,ustd,u,v,Iu
  integer(kind=8),dimension(nmeas) :: time

!!$ Local variables
  integer :: imeas
  real :: pi

  pi=atan(1.)*4.
  
  open(1,file=trim(filename))
  !!$ skip header lines
  read(1,*)
  read(1,*)
  read(1,*)
  read(1,*)
  do imeas = 1,nmeas
     read(1,*) rubbish,dir(imeas),vh(imeas),vmax(imeas),vmin(imeas), &
          ustd(imeas),time(imeas)
!!$     write(*,*) rubbish,dir(imeas),vh(imeas),vmax(imeas),vmin(imeas), &
!!$          ustd(imeas),time(imeas)
     u(imeas) = -sin(dir(imeas) * pi / 180.0) * vh(imeas)
     v(imeas) = -cos(dir(imeas) * pi / 180.0) * vh(imeas)
     Iu(imeas) = ustd(imeas) / vh(imeas)
  end do
  close(1)
776 format(a5,',',i3.3,',',f5.2,',',f5.2,',',f5.2,',',f5.2,',',i12.12,',')

end subroutine read_mjseries
!########################################################################
!########################################################################
subroutine write_mjseries(filename, &                      !!$ IN  : output file name
             label,labelref, &                             !!$ IN  : point label and reference label
             nmeas, &                                      !!$ IN  : number of records
             dirm_c,vhm, &                                 !!$ IN  : dir & vh of virtual series
             vhmax,vhmin, &                                !!$ IN  : max & minimum velocities (=99.99)
             ustd, &                                       !!$ IN  : virtual series of USTD
             time)                                         !!$ IN  : vector of time stamps

  implicit none

!!$ Interface variables
  integer :: nmeas
  character(len=5) :: label,labelref
  character(len=100) :: filename
  real,dimension(nmeas) :: dirm_c,vhm,vhmax,vhmin,ustd
  integer(kind=8),dimension(nmeas) :: time

!!$ Local variables
  integer :: imeas

  open(1,file=trim(filename))

!!$ Write header
  write(1,656)
  write(1,660) minval(time),maxval(time)
  write(1,666) nmeas
  write(1,670) 
  do imeas = 1,nmeas
     if (nint(dirm_c(imeas)).eq.360) then
        write(1,777) label,0,vhm(imeas), &
             vhmax(imeas),vhmin(imeas),ustd(imeas),time(imeas)
     else
        write(1,777) label,nint(dirm_c(imeas)),vhm(imeas), &
             vhmax(imeas),vhmin(imeas),ustd(imeas),time(imeas)
     end if
  end do

  close(1)

656 format(' :: TAB2DAT Artificial Time Series ::')
660 format('MJ|Inov, Start date ',i12,' :: End date ',i12)
666 format('Registos ',i7)
670 format('Est,dir,Vmed,raj,Vmin,desv.padrao,AAAAMMDDhhmm,')
777 format(a5,',',i3.3,',',f5.2,',',f5.2,',',f5.2,',',f5.2,',',i12.12,',')

end subroutine write_mjseries
!########################################################################
!########################################################################
subroutine read_tab(filename, &              !!$ IN  : input file name
     ndir,nbin, &                            !!$ IN  : number of sector and velocity bins
     fdir, &                                 !!$ OUT : sector frequency
     fbin_dir, &                             !!$ OUT : sector freq of velocity bins
     fbin_glo)                               !!$ OUT : global freq of velocity bins

  implicit none

!!$ Interface variables
  character(len=100) :: filename
  integer :: ndir,nbin
  real,dimension(ndir) :: fdir
  real,dimension(ndir,nbin) :: fbin_dir,fbin_glo

!!$ Local variables
  integer :: ibin,idir
  real :: rubbish

  open(1,file=trim(filename))
  read(1,*)
  read(1,*)
  read(1,*)

  read(1,*) fdir(1:ndir)
  fdir = fdir / 100.0

  do ibin = 1, nbin
     read(1,*) rubbish,fbin_dir(1:ndir,ibin)
  end do
  fbin_dir = fbin_dir / 1000.0

  do idir = 1, ndir
     fbin_glo(idir,1:nbin) = fdir(idir) * fbin_dir(idir,:)
  end do
  close(1)

end subroutine read_tab
!########################################################################
!########################################################################
subroutine write_tab(lenfile, &              !!$ IN  : length of filename 
     filename, &                             !!$ IN  : input file name
     ndir,nbin, &                            !!$ IN  : number of sector and velocity bins
     fdir, &                                 !!$ IN  : sector frequency
     fbin_dir, &                             !!$ IN  : sector freq of velocity bins
     fbin_glo, &                             !!$ IN  : global freq of velocity bins
     lenlabel,label)                         !!$ IN  : file name of DAT file to use as label

  implicit none

!!$ Interface variables
  integer :: ndir,nbin,lenlabel,lenfile
  character(len=lenfile) :: filename
  character(len=lenlabel) :: label
  real,dimension(ndir) :: fdir
  real,dimension(ndir,nbin) :: fbin_dir,fbin_glo

!!$ Local variables
  integer :: ibin,idir
  real :: rubbish

  open(1,file=trim(filename))
  write(1,*) trim(label)
  write(1,10) 0.0, 0.0, nbin * 1.0
  write(1,20) ndir, 1.0, 0.0

  fdir = fdir * 100.0
  write(1,30) fdir(1:ndir)

  do ibin = 1, nbin
     write(1,40) ibin*1.0,fbin_dir(1:ndir,ibin)*1000.0
  end do
  close(1)

10 format(6x,f7.5,5x,f7.5,3x,f5.2)
20 format(5x,i4,3x,f4.2,3x,f4.2)
30 format(9x,12(f7.3,2x))
40 format(1x,f6.2,12(1x,f8.2))
end subroutine write_tab
