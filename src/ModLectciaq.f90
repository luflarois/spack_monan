module ModLectciaq
    implicit none

contains

    subroutine lectciaq(ifdin)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Read chemical species.
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

       !     Michel Pirre, 28/01/2008  change to create chem1aq_list
       !
       !------------------------------------------------------------------------

       !     CHDON:  string variable corresponding to one line of the file.
       !     MOT: array composed of the words in CHDON.
       !     IMOT: array of sizes for words of MOT.
       !     NMOT: number of words.
       use ModParametre
       use ModFiccom, only: nesp, indaq
       use ModAuxnom, only: nblanc, nom_aux
       use ModGestion, only: part
       implicit none

       integer, intent(in out) :: ifdin
       integer, parameter :: nbmot = 100
       integer :: ie, k, itesaq, imp, nmot

       !common /nblanc/ nblanc

       character (len = 500) :: chdon
       character (len = 500) :: mot(nbmot)
       character(len = 3) :: cie, cieg
       character(len = 2) :: caq
       integer, dimension(nbmot) :: imot
       integer :: ii
       integer, dimension(nespmax) :: igas, inomaq
       character(len = 500) :: nomaq(nespmax)
       character(len = 20) :: cacco(nespmax)
       imp = 6

       !     Loop for reading the input file.

       read(ifdin, *)
       nomaq = ''
       open(unit = 178, file = 'chem1aq_list.f90', status = 'replace')
       write(unit = 178, fmt = '(A)') 'MODULE chem1aq_list'
       write(unit = 178, fmt = '(A)') '  IMPLICIT NONE'
       write(unit = 178, fmt = '(A)') '  '
       write(unit = 178, fmt = '(A)') '  '

       !WRITE(UNIT=178,FMT='(A,A,A,A)') '  CHARACTER(LEN=24),PARAMETER :: chemical_mechanismaq='&
       !          ,char(39),chemical_mechanism(1:len_trim(chemical_mechanism)),char(39)

       write(unit = 178, fmt = '(A)') '  INTEGER,PARAMETER :: maxnspeciesaq= 200'
       write(unit = 178, fmt = '(A,I3.3)') '  INTEGER,PARAMETER :: nspeciesaq=', nesp(3)
       write(unit = 178, fmt = '(A)') '  '
       write(unit = 178, fmt = '(A)') '  '
       write(unit = 178, fmt = '(A)') '  !Name of speciesaq '
       write(unit = 178, fmt = '(A,I3.3,A)') '  CHARACTER(LEN=8),PARAMETER,' // &
       'DIMENSION(nspeciesaq) :: spcaq_name=(/ &'
       caq = 'aq'
       do ie = 1, nesp(3)
          read(ifdin, '(a)') chdon
          print *, ie, chdon
          nblanc = nbmot
          call part(chdon, mot, imot, nmot, nblanc)
          indaq(ie) = 0
          nomaq(ie) (1:imot(1)) = mot(1) (1:imot(1))
          !  nomaq_aux(ie)=nomaq(ie)
          cacco(ie) = mot(2)
          inomaq(ie) = imot(1)
          ii = 4 - len(trim(nomaq(ie)))
          if (ie == 1) then
             write(unit = 178, fmt = '(A,A,A)') '   ' // char(39), trim(nomaq(ie)) // caq // repeat(' ', ii), char(39) // ' & !' 
          else
             write(unit = 178, fmt = '(A,A,A)') '     ,' // char(39), trim(nomaq(ie)) // caq // repeat(' ', ii), char(39) // ' & !' 
          end if
          do k = 1, ie - 1
             if (nomaq(k) (1:inomaq(k)) == nomaq(ie) (1:inomaq(ie))) then
                write(*, *) 'ERROR: aqueous species already initialized: ', nomaq(ie)
                stop
             end if
          end do
          !  END IF
          itesaq = 0
          do k = 1, nesp(2)
             !     IF (nom(k)(1:inom(k)) == nomaq(ie)(1:inomaq(ie))) THEN
             if (nom_aux(k) == nomaq(ie) (1:inomaq(ie))) then
                itesaq = itesaq + 1
                igas(ie) = k
             end if
          end do
          if (itesaq == 0) then
             write(*, *) 'ERROR: aqueous species has no corresponding gas species: ', nomaq(ie)
             stop
          end if
          if (itesaq > 1) then
             write(*, *) 'ERROR: gas species already initialized: ', nomaq(ie)
             stop
          end if
       end do
       write(unit = 178, fmt = '(A)') '   /)'
       write(unit = 178, fmt = '(A)') '  '
       write(unit = 178, fmt = '(A)') '  '

       write(unit = 178, fmt = '(A)') '  !Number of each specie   '
       do ie = 1, nesp(3)
          write(cie, fmt = '(I3.3)') ie
          ii = 4 - len(trim(nomaq(ie)))
          write(unit = 178, fmt = '(A)') '  INTEGER,PARAMETER :: ' // trim(nomaq(ie)) // caq // &
          repeat(' ', ii) // '=' // cie
       end do

       write(unit = 178, fmt = '(A)') '  '
       write(unit = 178, fmt = '(A)') '  '

       write(unit = 178, fmt = '(A)') '!     number of the corresponding gaseous species'
       write(unit = 178, fmt = '(A)') '  INTEGER,PARAMETER,DIMENSION(nspeciesaq) :: ind_gas=(/&'
       do ie = 1, nesp(3)
          write(cie, fmt = '(I3.3)') ie
          write(cieg, fmt = '(I3.3)') igas(ie)
          if (ie < (nesp(3))) then
             write(unit = 178, fmt = '(A)') &
             '    ' // cieg // ' ,   & ! ' // trim(nomaq(ie)) // caq // ' - ' // cie
          else
             write(unit = 178, fmt = '(A)') &
             '    ' // cieg // '     & ! ' // trim(nomaq(ie)) // caq // ' - ' // cie
          end if
       end do
       write(unit = 178, fmt = '(A)') '    /)'
       write(unit = 178, fmt = '(A)') '  '

       write(unit = 178, fmt = '(A)') '! accomodation coefficient'
       write(unit = 178, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspeciesaq) :: acco=(/&'
       do ie = 1, nesp(3)
          write(cie, fmt = '(I3.3)') ie
          if (ie < (nesp(3))) then
             write(unit = 178, fmt = '(A)') &
             '    ' // cacco(ie) // ' ,   & ! ' // trim(nomaq(ie)) // caq // ' - ' // cie
          else
             write(unit = 178, fmt = '(A)') &
             '    ' // cacco(ie) // '     & ! ' // trim(nomaq(ie)) // caq // ' - ' // cie
          end if
       end do
       write(unit = 178, fmt = '(A)') '    /)'

       write(unit = 178, fmt = '(A)') '  '
       write(unit = 178, fmt = '(A)') 'END MODULE chem1aq_list'
       close(unit = 178)

    end subroutine lectciaq

end module ModLectciaq