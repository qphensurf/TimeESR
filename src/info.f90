MODULE info
   
   USE OpenFiles, ONLY: stdout
   USE art, ONLY: show_art
   
   IMPLICIT NONE
   SAVE
   
   PRIVATE
   PUBLIC :: environment_start, environment_end

CONTAINS

   SUBROUTINE environment_start( code )

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: code
      CHARACTER(LEN=80) :: code_version, version_number

      version_number = '1.0.0'
      code_version = TRIM(code) // " v." // TRIM(version_number)

      CALL opening_message(code_version)

   END SUBROUTINE environment_start


   SUBROUTINE environment_end( )

      IMPLICIT NONE


      CALL closing_message( )
      FLUSH(stdout)

      RETURN

   END SUBROUTINE environment_end


   SUBROUTINE opening_message( code_version )

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN) :: code_version
      CHARACTER(LEN=9)  :: cdate, ctime

      CALL show_art()
      CALL date_and_tim(cdate,ctime)
      WRITE( stdout, '(/5X,A," started on ",A9," at ",A9)' ) &
         TRIM(code_version), cdate, ctime
#ifdef __DEBUG
      WRITE( stdout, '(5X,"(DEBUG MODE)")' )
#endif
      WRITE( stdout, '(5X,"" )' )

   END SUBROUTINE opening_message


   SUBROUTINE closing_message( )
  
      IMPLICIT NONE

      CHARACTER(LEN=9)  :: cdate, ctime
      CHARACTER(LEN=80) :: time_str

      CALL date_and_tim( cdate, ctime )

      time_str = 'Program terminated on:  ' // ctime // ' ' // cdate

      WRITE( stdout,*)
      WRITE( stdout,3334) time_str
      WRITE( stdout,3335)

3334  FORMAT(3X,A60,/)
3335  FORMAT('=',78('-'),'=')

      RETURN
   END SUBROUTINE closing_message


   SUBROUTINE date_and_tim( cdate, ctime )
   
      IMPLICIT NONE
      
      CHARACTER(LEN=9) :: cdate, ctime
      CHARACTER(LEN=3), DIMENSION(12) :: months
      DATA months /'Jan','Feb','Mar','Apr','May','Jun', 'Jul','Aug','Sep','Oct','Nov','Dec'/
      INTEGER date_time(8)
      
      CALL date_and_time(values=date_time)
     
      WRITE (cdate,'(i2,a3,i4)') date_time(3), months(date_time(2)), date_time(1)
      WRITE (ctime,'(i2,":",i2,":",i2)') date_time(5), date_time(6), date_time(7)

   END SUBROUTINE date_and_tim

END MODULE info

