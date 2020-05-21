!character(4),parameter::step_string='0587'





subroutine read_data_x( x, n, i, filename) bind(c)
  use iso_c_binding
  
  integer(c_int)::i
  integer(c_int)::n
  real(c_float)::x(3,n)
  integer(c_int)::filename

  !parameters
  real(8),parameter::box=1200. !boxsize
  real(8),parameter::omega0=0.268  
  real(8),parameter::mp=2.77e11*omega0*box**3/3072.**3 !particle mass
  real(8),parameter::pi=4.*atan(1.)

  !work direction
  character(256),parameter::simudir='/data/s1/simu/Jing6620/'
  character(4)::simucode='6620'                                                                     

  !redshift string
  character(4)::step_string="3356"
  character(7),parameter::zstring='0.000'

  !read parameters
  integer(8)::np,ips
  real(4)::ztp,omegat,lambdat,rLbox,xscale,vscale
  real(4)::alpha=1.,dp=0.0288,nstepf=5000.,scalep,scalepf,vfact2
  real(4)::a,H


  integer(8)::pid
  real(4)::t1,t2

  integer(4)::fid,uid
  character(2)::f_string
  character(100)::fname
  real(8)::li

  character(40)::void_char
  integer(8)::loca !the location of file

  character(4)::fn 
  write(fn,'(i4.4)') filename 
  
  fname=trim(simudir)//'pos'//simucode//'.'//fn//'.'//'01'
  write(*,*) 'reading:',trim(fname) 
  open(31,file=fname,status='old',form='unformatted')
  read(31) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
  close(31)
  !write(*,*) 'particle number:',np
  write(*,*) 'ips:',ips
  write(*,*) 'ztp:',ztp
  !write(*,*) 'a:',1./(1.+ztp)
  !write(*,*) 'omegat:',omegat
  !write(*,*) 'lambdat:',lambdat
  !write(*,*) 'boxsize:',rLbox
  !write(*,*) 'xscale:',xscale
  !write(*,*) 'vscale:',vscale
  

  !read position
  write(*,*) "now in fortran"
  call cpu_time(t1)
  fid=i
  write(f_string,'(i2.2)') fid
  fname=trim(simudir)//'pos'//simucode//'.'//fn//'.'//f_string
  write(*,*) 'reading:',trim(fname)
  uid=30+fid
  open(uid,file=fname,status='old',form='unformatted')
  if (fid.eq.1) read(uid) void_char

  
  write(*,*) 'reading particle from',np*(i-1)/4+1,'to',np*i/4
  read(uid) x(1:3,1:np/4)
  close(uid)
  call cpu_time(t2)
  write(*,*) 'time consumed for reading position:',t2-t1,'seconds'

end subroutine read_data_x


subroutine read_data_v( v, n, i, filename) bind(c)
  use iso_c_binding

  integer(c_int)::i
  integer(c_int)::n
  real(c_float)::v(3,n)
  integer(c_int)::filename

  !parameters
  real(8),parameter::box=1200. !boxsize
  real(8),parameter::omega0=0.268
  real(8),parameter::mp=2.77e11*omega0*box**3/3072.**3 !particle mass
  real(8),parameter::pi=4.*atan(1.)

  !work direction
  character(256),parameter::simudir='/data/s1/simu/Jing6620/'
  character(4)::simucode='6620'

  !redshift string
  character(4),parameter::step_string='3356'
  character(7),parameter::zstring='0.000'

  !read parameters
  integer(8)::np,ips
  real(4)::ztp,omegat,lambdat,rLbox,xscale,vscale
  real(4)::alpha=1.,dp=0.0288,nstepf=5000.,scalep,scalepf,vfact2
  real(4)::a,H


  integer(8)::pid
  real(4)::t1,t2

  integer(4)::fid,uid
  character(2)::f_string
  character(100)::fname
  real(8)::li

  character(40)::void_char
  integer(8)::loca !the location of file

  character(4)::fn
  write(fn,'(i4.4)') filename

  fname=trim(simudir)//'pos'//simucode//'.'//fn//'.'//'01'
  write(*,*) 'reading:',trim(fname)
  open(31,file=fname,status='old',form='unformatted')
  read(31) np,ips,ztp,omegat,lambdat,rLbox,xscale,vscale
  close(31)
  !write(*,*) 'particle number:',np
  !write(*,*) 'ips:',ips
  !write(*,*) 'ztp:',ztp
  !write(*,*) 'a:',1./(1.+ztp)
  !write(*,*) 'omegat:',omegat
  !write(*,*) 'lambdat:',lambdat
  !write(*,*) 'boxsize:',rLbox
  !write(*,*) 'xscale:',xscale
  !write(*,*) 'vscale:',vscale

  !read velosity
  call cpu_time(t1)
  fid=i
  write(f_string,'(i2.2)') fid

  fname=trim(simudir)//'vel'//simucode//'.'//fn//'.'//f_string
  write(*,*) 'reading:',trim(fname)
  uid=30+fid
  open(uid,file=fname,status='old',form='unformatted')
  if (fid.eq.1) read(uid) void_char

  read(uid) v(1:3,1:np/4)
  close(uid)
  call cpu_time(t2)
  write(*,*) 'time consumed for reading volesity:',t2-t1,'seconds'

end subroutine read_data_v

                                                                                                                                                                                            96,0-1        Bot

