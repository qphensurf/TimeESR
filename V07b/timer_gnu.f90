module timer
Contains

  subroutine clock(message,i)
  character ( len=* ), intent ( in ) :: message
  real(4), save :: savedTime
  real(4)       :: elapsedTime
  integer       :: i

  if (i == 1) then

  call cpu_time (savedTime)

  endif

  call cpu_time (elapsedTime)

  write (*, *) '====================================================================='
  write (*,'("  ",a,f8.2," s")') "Time: "//message, elapsedTime - savedTime
  write (*, *) '====================================================================='

  savedTime = elapsedTime

  end subroutine

end module timer
