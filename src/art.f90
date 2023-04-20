MODULE art

   USE OpenFiles, ONLY: stdout
   
   IMPLICIT NONE
   SAVE
   
   PRIVATE
   PUBLIC :: show_art
  
CONTAINS
   
   SUBROUTINE show_art()

      IMPLICIT NONE
                      
      WRITE( stdout,'(5X," _____ _                _____ ___________ ",&
                    &/5X,"|_   _(_)              |  ___/  ___| ___ \",&
                    &/5X,"  | |  _ _ __ ___   ___| |__ \ `--.| |_/ /",&
                    &/5X,"  | | | | ._ ` _ \ / _ \  __| `--. \    / ",&
                    &/5X,"  | | | | | | | | |  __/ |___/\__/ / |\ \ ",&
                    &/5X,"  \_/ |_|_| |_| |_|\___\____/\____/\_| \_|")' )
                    
   END SUBROUTINE show_art
   
END MODULE art


