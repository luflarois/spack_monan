module ModLectcinet
    use ModParametre
    implicit none

    contains    

    subroutine lectcinet(ifdth, indicaq, ntuvonline, filename)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Read chemical mechanism.
       !------------------------------------------------------------------------

       !     -- INPUT VARIABLES

       !     IFDTH: input file
       !     FPAR: physical parameters file

       !     -- INPUT/OUTPUT VARIABLES

       !     -- OUTPUT VARIABLES

       !     NESP: number of species.
       !     NR: number of reactions.
       !     XLPHY: physical lumping.
       !     INDPUR: ii=indpur(i,j) true label of j-th species in i-th lumping.
       !     NALG: number of algebraic constraints.

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
       use ModFiccom, only: nequil, nr, nrphot, jer, s
       use ModAuxnom, only: nblanc
       use ModGestion, only: part
       use ModWrite, only: write_header, write_end
       use ModCinet, only: kinreac, kindis, kinhenry, initphase
       use ModGenerator, only: wnonzero
       implicit none

       integer, intent(in out) :: ifdth
       integer, intent(in out) :: indicaq
       integer, intent(in out) :: ntuvonline
       character(len = *) :: filename
       integer, parameter :: nbmot = 100

       !common /nblanc/ nblanc

       character (len = 500) :: chdon
       character (len = 500) :: mot(nbmot)
       integer, dimension(nbmot) :: imot

       integer :: imp, iattente, iarret, nchold, nchplus, nmot

       imp = 6
       iattente = 0
       iarret = 0

       !     Loop for reading the file.
       !     Write the first lines of the output routines
       call write_header(ntuvonline)

       100 read(ifdth, '(a)') chdon
       write(imp, '(a)') chdon

       !     Build the sequence and decompose in words.


       nblanc = nbmot
       call part(chdon, mot, imot, nmot, nblanc)



       !     Case of long sequences (two lines).

       102 continue
       if (mot(nmot) (1:2) == '//') then
          nchold = nmot
          read(ifdth, '(a)') chdon
          write(imp, '(a)') chdon
          nblanc = nbmot - nchold + 1


          call part(chdon, mot(nchold), imot(nchold), nchplus, nblanc)
          nmot = nchold - 1 + nchplus
       end if

       if (mot(nmot) (1:2) == '//') go to 102

       !     Call the different routines according to the keyword.

       if (mot(1) (1:3) == 'KIN') then
          if (iattente == 0) then
             write(*, *) 'ERROR: kinetics before reaction ', nr
             stop
          else if (iattente == 1) then
            call kinreac(mot, imot, nmot, ntuvonline)
          else if (iattente == 2) then
            call kindis(mot, imot, nmot)
          else if (iattente == 3) then
            call kinhenry(mot, imot, nmot)
          end if
          iattente = 0

       else if (mot(1) (1:3) == 'SET') then

          if (mot(2) (1:4) == 'UNIT') then
            call lectunit(mot, imot, nmot)

          else if (mot(2) (1:10) == 'TABULATION') then
            call lect_tabulation(mot, imot, nmot)

          else
             write(*, *) 'ERROR: UNKNOWN SET FUNCTIONS'
             stop
          end if
          !     Symbols for commented lines: %, !,
       else if (mot(1) (1:1) == '%') then
       else if (mot(1) (1:1) == '!') then

          !     END.
       else if (mot(1) (1:3) == 'END') then
          iarret = 1
       else if (iattente /= 0) then
          write(*, *) 'ERROR: I wait for kinetics'
          stop
       else
         call reaction(mot, imot, nmot, iattente)
       end if

       if (iarret == 0) go to 100


       call write_end

       !     Write file non_zero.dat

       call wnonzero(s, nr, jer)

       write(*, *) '########################################'
       write(*, *) '########################################'
       write(6, *) 'Summary for the kinetic scheme'
       write(*, *) 'Total number of reactions =', nr
       write(*, *) 'Number of photolytic reactions =', nrphot
       write(*, *) 'Number of dissociation equilibria =', nequil

       if (indicaq == 0) then
          write(*, *) 'Gas-phase chemistry'
       end if

       if (indicaq == 1) then
          write(*, *) 'Multiphase (gas-phase and aqueous-phase) chemistry'
       end if

       call initphase(ntuvonline, filename)

    end subroutine lectcinet 

   subroutine lectunit(mot, imot, nmot)
      !------------------------------------------------------------------------

      !     -- DESCRIPTION

      !     Read units.
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
      use ModParametre
      use ModFiccom, only: iunitgas, iunitaq
      implicit none

      integer, parameter :: nbmot = 100

      character (len = 500), intent(in out) :: mot(nbmot)
      integer, intent(in out) :: imot(nbmot)
      integer, intent(in out) :: nmot

      iunitgas = 0
      if (mot(2) (1:imot(2)) == 'GAS') then
         if (mot(3) (1:imot(3)) == 'MOLCM3') then
            iunitgas = 0
         else if (mot(3) (1:imot(3)) == 'PPB') then
            iunitgas = 1
         else
            write(*, *) 'ERROR: unknown units for gas-phase kinetics'
            stop
         end if
      else if (mot(2) (1:imot(2)) == 'AQ') then
         if (mot(3) (1:imot(3)) == 'MOLL') then
            iunitaq = 2
         else
            write(*, *) 'ERROR: unknown units for aqueous-phase kinetics'
            stop
         end if
      end if
      return
   end subroutine lectunit

   subroutine lect_tabulation(mot, imot, nmot)
      !------------------------------------------------------------------------

      !     -- DESCRIPTION

      !     Read tabulation for photolysis.
      !     The format is: TABULATION N DEGREES D1 D2 ... DN
      !     N is the number of tabulated angles of values Di (in degrees).
      !     The sequence may be increasing or decreasing.
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
      use ModParametre
      use ModFiccom, only: tabphot, ntabphotmax, ireversetab, ntabphot
      use ModGestion, only: entier, reel
      implicit none

      integer, parameter :: nbmot = 100

      character (len = 500), intent(in out) :: mot(nbmot)
      integer, intent(in out) :: imot(nbmot)
      integer, intent(in out) :: nmot

      double precision, dimension(ntabphotmax) :: t
      integer :: i, j

      ireversetab = 0

      call entier(ntabphot, mot(3), imot(3))
      if (ntabphot == 0) then
         write(*, *) 'ERROR: number of tabulated angles to be checked in'
         write(*, *) 'subroutine lect-tabulation'
         stop
      end if
      if (ntabphot > ntabphotmax) then
         write(*, *) 'ERROR: ntabphot>ntabphotmax'
         stop
      end if

      do i = 1, ntabphot
         call reel(t(i), mot(i + 4), imot(i + 4))
      end do
      !     Check increasing order
      if (t(1) > t(2)) then
         ireversetab = 1
         do j = 1, ntabphot
            tabphot(j) = t(ntabphot + 1 - j)
         end do
      else
         do j = 1, ntabphot
            tabphot(j) = t(j)
         end do
      end if
      do j = 1, ntabphot - 1
         if (tabphot(j) >= tabphot(j + 1)) then
            write(*, *) 'ERROR: the tabulation has to be strictly monotonic'
            stop
         end if
      end do

   end subroutine lect_tabulation

   subroutine reaction(mot, imot, nmot, iattente)

      !     -- DESCRIPTION

      !     Read reactions.
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
      use ModParametre
      use Modcinet, only: creac, cdis, chenry
      implicit none

      integer, parameter :: nbmot = 100

      character (len = 500), intent(in out) :: mot(nbmot)
      integer, intent(in out) :: imot(nbmot)
      integer, intent(in out) :: nmot
      integer, intent(out) :: iattente

      integer :: i

      do i = 1, nmot
         if (mot(i) (1:imot(i)) == '>') then
            if (iattente /= 0) then
               write(*, *) 'ERROR: too many symbols >'
               stop
            else
               iattente = 1
            end if
         else if (mot(i) (1:imot(i)) == '->') then
            if (iattente /= 0) then
               write(*, *) 'ERROR: too many symbols ->'
               stop
            else
               iattente = 1
            end if
         else if (mot(i) (1:imot(i)) == '=') then
            if (iattente /= 0) then
               write(*, *) 'ERROR: too many symbols ='
               stop
            else
               iattente = 2
            end if
         else if (mot(i) (1:imot(i)) == '<H>') then
            if (iattente /= 0) then
               write(*, *) 'ERROR: too many symbols <H>'
               stop
            else
               iattente = 3
            end if
         else if (mot(i) (1:imot(i)) == '=H=') then
            if (iattente /= 0) then
               write(*, *) 'ERROR: too many symbols =H='
               stop
            else
               iattente = 4
            end if
         end if
      end do

      if (iattente == 0) then
         write(*, *) 'ERROR: lack of symbol for reaction'
         stop
      else if (iattente == 1) then
         call creac(mot, imot, nmot)
      else if (iattente == 2) then
         call cdis(mot, imot, nmot)
      else if (iattente == 3) then
         call chenry(mot, imot, nmot, 0)
      else if (iattente == 4) then
         call chenry(mot, imot, nmot, 1)
         iattente = 3
      end if

   end subroutine reaction


end module ModLectcinet