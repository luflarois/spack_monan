program main

   !------------------------------------------------------------------------

   !     -- DESCRIPTION

   !     Master code for SPACK:
   !     Simplified Preprocessor for Atmospheric Chemical Kinetics.

   !     The input files are:
   !     - a file describing the mechanism in a symbolic way,
   !     - a file of chemical species.

   !     The output files are F77 routines describing the chemical
   !     production terms.
   !------------------------------------------------------------------------

   !     -- INPUT VARIABLES

   !     -- INPUT/OUTPUT VARIABLES

   !     -- OUTPUT VARIABLES

   !------------------------------------------------------------------------

   !     -- REMARKS

   !------------------------------------------------------------------------

   !     -- MODIFICATIONS
   !     Rodrigues, L.F., INPE, 2026. luflarois@gmail.com
   !     Converted to Fortran90 and using modules for better code organization 
   !     and maintainability. Useless variables and common blocks have been removed, 
   !     and the code has been  restructured to improve readability and efficiency.

   !------------------------------------------------------------------------

   !     -- AUTHOR(S)
   !     Bruno Sportisse, CEREA, 2003.

   !------------------------------------------------------------------------
   use ModParametre, only: nespmax
   use ModInit, only: lectdata
   use ModMonan, only: create_monan_registry
                
   implicit none
   
   double precision, dimension(nespmax) :: y0
   integer :: neq, indicaq, is_monan


   !     Initialization of data
   !     NEQ has the right dimension.
   call lectdata(y0, neq, indicaq, is_monan)
   if (is_monan == 1) then
      print *, 'MONAN includes will be generated.'
      call create_monan_registry()
   else
      print *, 'MONAN includes will NOT be generated.'
   end if

end program main







