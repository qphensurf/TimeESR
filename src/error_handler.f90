module error_handler

   implicit none

   public :: errore

contains

   subroutine errore( calling_routine, message, ierr )
   ! Similar to Quantum Espresso errore() command, with 
   ! a few options removed. This is a generic error 
   ! handling module.
   ! if ierr <= 0 it does nothing,
   ! if ierr  > 0 it stops.
      implicit none
      
      character(len=*), intent(in) :: calling_routine, message
      integer, intent(in) :: ierr
      character(len=6) :: cerr
      
      if(ierr <= 0 ) return
      
      write(cerr,fmt='(I6)') ierr
      write(*,fmt='(/,1x,78("%"))')
      write(*,fmt='(5X,"Error in routine ",A," (",A,"):")') &
           TRIM(calling_routine), TRIM(ADJUSTL(cerr))
      write(*,fmt='(5x,A)') TRIM(message)
      write(*,fmt='(1x,78("%"),/)')
      write( *, '(5x,"stopping ...")' )
      stop 1
      
      return
      
   end subroutine errore

end module error_handler