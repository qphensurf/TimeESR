program background
! Simple tool to remove a linear background from a CW spectrum
  implicit none
  integer, parameter :: q = SELECTED_REAL_KIND(10)
  integer, parameter :: N = 500
  integer :: i
  real (q) :: freq, datum, datini, datfin, dummy 
  real (q) :: ini, fin, slope, absc
  character ( len = 100 ) :: namef

  print *, 'Please, write the initial frequency:'
  read (*,*) ini
  print *, 'Please, write the initial current:'
  read (*,*) datini
  print *, 'Please, write the final frequency:'
  read (*,*) fin
  print *, 'Please, write the final current:'
  read (*,*) datfin
  print *, 'Please, write the name of the file:'
  read (*,*) namef
    open (1, file=namef, status='old')
    open (2, file='corrected.dat')

 slope = (datfin-datini)/(fin - ini)
 absc = datini

  do i = 1, N

     read (1, *) dummy, freq, datum
     write (2, *) freq, datum-(slope*(freq-ini)+absc), dummy

 enddo
 stop
 
end program background
