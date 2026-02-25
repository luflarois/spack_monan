module ModInit
    !Module for initialization of the simulation
    !Author: rodrigues, l.f. 
    !email: luflarois@gmail.com
    !Date: 2026Jan09
    use ModParametre
    implicit none

    character, dimension(nespmax) :: cdry, cwet, coff, cfdd, ctra
    character :: csou(nespmax)
    
    contains


    subroutine lectdata(y0, neq, indicaq, is_monan)
        !------------------------------------------------------------------------

        !     -- DESCRIPTION

        !     Read data.
        !------------------------------------------------------------------------

        !     -- INPUT VARIABLES

        !     -- INPUT/OUTPUT VARIABLES

        !     -- OUTPUT VARIABLES

        !     SLUMP: lumped stoichiometric matrix.
        !     XLPHY: physical lumping.
        !     INDPUR: ii=indpur(i,j) true label of J-th  species in lumping I.
        !     y0: initial conditions.
        !     S: stoichiometrix matrix.
        !     NALG: physical algebraic onstraints.

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
        use ModFiccom, only: nesp, indicaqcom, nalg, nom, ndiff
        use ModNficfort, only: filemeca, filespecies, ipiste
        use ModAuxnom, only: chemical_mechanism
        use ModLectcinet, only: lectcinet
        use ModLectciaq, only: lectciaq
        
        implicit none

        integer, parameter :: nbrem = 0

        double precision, intent(in out) :: y0(nespmax)
        integer, intent(out) :: neq
        integer, intent(out) :: indicaq
        integer, intent(out) :: is_monan
        integer :: ntuvonline
        character (len = 20) :: fd
        character(len = 100) :: filename
        integer :: i, ie, ifdin, ippb, nmaster
        integer :: ifdth

        ipiste = 15
        nmaster = ipiste
        ipiste = ipiste + 1
        !     Initialization of constants and main arrays.
        call initconst
        call initcinet

        nalg = 0

        !     Read master file.
        open(nmaster, file = 'inSPACK', status = 'old')
        write(*, *) '-------------Input files-----------'

        !     Read photolysis rates option
        read(nmaster, *)
        read(nmaster, *) ntuvonline

        read(nmaster, *)
        !iunit = 0
        read(nmaster, *)

        !     Chemical mechanism
        read(nmaster, *) fd
        write(6, *) 'File for chemical mechanism: ', fd
        filemeca = fd
        print *, 'filemeca=', filemeca
        do i = len_trim(filemeca), 1, -1
           if (filemeca(i:i) == '/') exit
        enddo
        chemical_mechanism = trim(filemeca(i + 1:len_trim(filemeca)))
        print *, 'chemical_mechanism=', chemical_mechanism

        call openfic(ipiste, ifdth, fd, 0)
        read(nmaster, *)
        read(nmaster, *)
        read(nmaster, *) fd
        !      read(nmaster,*)
        read(nmaster, *)
        read(nmaster, *)
        read(nmaster, *) filename !To FastJX
        print *, 'FastJx adapter Mechanism: ', filename

      !Incluindo se é ou não para gerar os includes do MONAN
      read(nmaster, *)
      read(nmaster, *)
      read(nmaster, *) is_monan



        close(nmaster)
        !     Chemical species
        write(6, *) 'File for chemical species: ', fd
        filespecies = fd
        call openfic(ipiste, ifdin, fd, 0)

        !     Modif BS multiphase: units molec. cm-3 (ippb=0) ou ppb (1)

        read(ifdin, *)
        read(ifdin, *)
        !     read(Ifdin,*)ippb
        !     read(Ifdin,*)
        !     End modif BS

        ippb = 0
        read(ifdin, *) nesp(2), nesp(3)

        nesp(1) = nesp(2) + nesp(3)
        indicaq = 0
        if (nesp(3) > 0) indicaq = 1
        indicaqcom = indicaq
        write(6, *) 'Number of multiphase species (nesp): ', nesp(1)
        write(6, *) 'Number of gas-phase species: ', nesp(2)
        write(6, *) 'Number of aqueous-phase species: ', nesp(3)
        write(6, *) 'Max number of species (nespmax): ', nespmax
        write(6, *) 'Max number of reactions (nrmax): ', nrmax
        if (nesp(1) > nespmax) then
           write(*, *) 'ERROR: dimension, nesp>nespmax'
           call halte
        end if
        call lectci(ifdin)

        !     ..v.7..x....v....x....v....x....v....x....v....x....v....x....v....x.I
        !     PP 12 02 2002
        !     Check species name
        do ie = 1, nesp(1)
           write(6, 777) ie, nom(ie)
        end do
        777 format(2x, i3, 2x, a10)

        ! change MP 29/01/08 to create chem1_listaq for aqueous species
        if (nesp(3) > 0) then
           call lectciaq(ifdin)
           !DO ie = 1,nesp(1)
           !  WRITE(6,777)ie,nom(ie)
           !END DO
        end if
        !end change MP

        !     ..v.7..x....v....x....v....x....v....x....v....x....v....x....v....x.I
        !     Read chemical mechanism

        write(*, *)
        write(*, *) '-----------Chemical mechanism---------'
        write(*, *)
        call lectcinet(ifdth, indicaq, ntuvonline, filename)
        !     Modif BS for gas-phase only
        !     call convcinet(iunit,indicaq)


        !     Modif BS for gas-phase only
        !     call inition
        !     call initlphy
        !     call initphot
        !     End modif gas-phase only.

        !     Dimension

        neq = ndiff(1) * nbrem

        100 format(a10)

    end subroutine lectdata

    subroutine initconst
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization of constants.
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
       !     Michel Pirre, 28/01/2008 : change to create chem1_lisaq

       !------------------------------------------------------------------------
       use ModFiccom

       av = 6.022d23

    end subroutine initconst


    subroutine initcinet
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization.
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
       use ModFiccom, only: iunitaq, iunitgas, nr, nrphot, nequil, ntabphot, nrp &
                            , nrmol1, nrmol2, nrmol3, naq2, narr1, narr2, narr3 &
                            , narr4, narr5, narr6, narr7, narr8, nthird, rmol &
                            , ihreac, seqion, iphotinv, molec, bp, s, johmoin &
                            , jhplus, iheq, nb, nneq, iemonodi, iedimono, irmonodi &
                            , imolec1, imolec2, imolec3, iarr1, iarr2, iarr3, iarr4 &
                            , iarr5, iarr6, iarr7, iarr8, inaq2, indthird, idifford &
                            , nequil11, nequil12, nequil13, nequil21, nequil22, nequil23 &
                            , nequil3 , nequil41, nequil42, nequil51, nequil52, ispebp &
                            , jer, ittb, iprecalc, jhpoh, jion1, jion2, jaq, ndiff

        implicit none
        integer :: i, ir, j


       iunitaq = -1
       iunitgas = -1
       nr = 0
       nrphot = 0
       nequil = 0
       ntabphot = 0

       do i = 1, nphase
          nrp(i) = 0
          nrmol1(i) = 0
          nrmol2(i) = 0
          nrmol3(i) = 0
          naq2(i) = 0
          narr1(i) = 0
          narr2(i) = 0
          narr3(i) = 0
          narr4(i) = 0
          narr5(i) = 0
          narr6(i) = 0
          narr7(i) = 0
          narr8(i) = 0
          nthird(i) = 0
       end do

       do i = 1, nespmax
          rmol(i) = 0.d0
          ihreac(i) = 0
          do j = 1, nionx
             seqion(i, j) = 0.0d0
          end do
       end do

       do ir = 1, nrmax
          iphotinv(ir) = 0
          molec(ir) = 0
          bp(1, ir) = 0.d0
          bp(2, ir) = 0.d0
          bp(3, ir) = 0.d0
          bp(4, ir) = 0.d0
          bp(5, ir) = 0.d0
          bp(6, ir) = 0.d0
          ispebp(ir) = 0
          jer(1, ir) = 0
          jer(2, ir) = 0
          jer(3, ir) = 0
          ittb(ir) = 0
          iprecalc(ir) = 0

          do i = 1, nphase
             irmonodi(ir, i) = 0
             imolec1(ir, i) = 0
             imolec2(ir, i) = 0
             imolec3(ir, i) = 0
             iarr1(ir, i) = 0
             iarr2(ir, i) = 0
             iarr3(ir, i) = 0
             iarr4(ir, i) = 0
             iarr5(ir, i) = 0
             iarr6(ir, i) = 0
             iarr7(ir, i) = 0
             iarr8(ir, i) = 0
             inaq2(ir, i) = 0
             indthird(ir, i) = 0
          end do

          irmonodi(ir, 1) = ir
          iheq(ir) = 0
          nb(ir) = 0
          nneq(ir) = 0

          do i = 1, nespmax
             s(i, ir) = 0.0d0
          end do

       end do

       do i = 1, nespmax
          iemonodi(i, 1) = i
          iedimono(i, 1) = i
          do j = 2, nphase
             iemonodi(i, j) = 0
             iedimono(i, j) = 0
          end do
       end do

       do i = 1, nphase
          ndiff(i) = 0
       end do

       nequil11 = 0
       nequil12 = 0
       nequil13 = 0
       nequil21 = 0
       nequil22 = 0
       nequil23 = 0
       nequil3 = 0
       nequil41 = 0
       nequil42 = 0
       nequil51 = 0
       nequil52 = 0

       jhplus = 0
       johmoin = 0

       do i = 1, nequilx

          jhpoh(i) = 0
          idifford(i) = 0
          jion1(i) = 0
          jion2(i) = 0
          jaq(i) = 0

       end do

    end subroutine initcinet

    subroutine openfic(ipiste, ifd, fd, inew)
        implicit none

        integer, intent(in out) :: ipiste
        integer, intent(out) :: ifd
        character (len = 20), intent(in) :: fd
        integer, intent(in) :: inew

        ifd = ipiste
        ipiste = ipiste + 1
        if (inew == 0) open(ifd, file = fd, status = 'old')
        if (inew == 1) open(ifd, file = fd, status = 'new')

    end subroutine openfic

    subroutine halte
        stop
    end subroutine halte

    subroutine lectci(ifdin)
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

       !------------------------------------------------------------------------

       !     CHDON:  string variable corresponding to one line of the file.
       !     MOT: array composed of the words in CHDON.
       !     IMOT: array of sizes for words of MOT.
       !     NMOT: number of words.
       use ModFiccom, only: nom, indaq, inom, nesp
       use ModAuxnom, only: nom_aux, chemical_mechanism, nblanc
       use ModGestion, only: part

       integer, intent(in out) :: ifdin
       integer, parameter :: nbmot = 100

       character(len=4) , dimension(nespmax)  :: cwei
       character(len=3) , dimension(nespmax)  :: f0, difrat, mat_eq
       character(len=2) , dimension(nespmax)  :: uv_eq
       character(len=1) , dimension(nespmax)  :: aer_eq
       character(len=20), dimension(nespmax) :: hstar
       character(len=13), dimension(nespmax) :: dvj, dhr
       character(len=12), dimension(nespmax) :: ak0
       character(len=10), dimension(nespmax) :: dak
       character (len = 500) :: chdon
       character (len = 500) :: mot(nbmot)
       character(len = 3) :: cie
       character(len = 3) :: esp
       character(len = 1) :: caer
       integer, dimension(nbmot) :: imot
       integer :: ii, ia, ie, imp, k, nmot
       logical :: donthaveco2 = .true.
 
       imp = 6
       !     Loop for reading the input file.

       read(ifdin, *)
       nom = ''
       open(unit = 78, file = 'chem1_list.f90', status = 'replace')
       write(unit = 78, fmt = '(A)') 'MODULE chem1_list'
       write(unit = 78, fmt = '(A)') '  IMPLICIT NONE'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A,A,A,A)') '  CHARACTER(LEN=24),PARAMETER :: chemical_mechanism=' &
       , char(39), chemical_mechanism(1:len_trim(chemical_mechanism)), char(39)

       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: maxnspecies= 200'
       write(unit = 78, fmt = '(A,I3.3)') '  INTEGER,PARAMETER :: nspecies=', nesp(2)
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  !Name of species '
       write(unit = 78, fmt = '(A,I3.3,A)') '  CHARACTER(LEN=8),PARAMETER,' // &
       'DIMENSION(nspecies) :: spc_name=(/ &'
       do ie = 1, nesp(2)
          read(ifdin, '(a)') chdon
          print *, ie, chdon
          nblanc = nbmot
          call part(chdon, mot, imot, nmot, nblanc)
          indaq(ie) = 0
          nom(ie) (1:imot(1)) = mot(1) (1:imot(1))
          print *, 'LFR: ie,nom(ie)=', ie,nom(ie)
          nom_aux(ie) = nom(ie)
          !LFR
          !LSPR
          cwei(ie) = mot(2)
          csou(ie) = mot(3)
          cdry(ie) = mot(4)
          cwet(ie) = mot(5)
          !cpas(ie)=mot(6)
          cfdd(ie) = mot(6)
          !cfut(ie)=mot(7)
          coff(ie) = mot(7)
          !coff(ie)=mot(8)
          ctra(ie) = mot(8)
          hstar(ie) = mot(9)
          f0(ie) = mot(10)
          difrat(ie) = mot(11)
          uv_eq(ie) = mot(12)
          dvj(ie) = mot(13)
          dhr(ie) = mot(14)
          ak0(ie) = mot(15)
          dak(ie) = mot(16)
          aer_eq(ie) = mot(17) !matrix


          !LFR
          inom(ie) = imot(1)
          ii = 4 - len(trim(nom(ie)))
          if (ie == 1) then
             write(unit = 78, fmt = '(A,A,A)') '   ' // char(39), trim(nom(ie)) // repeat(' ', ii), char(39) // ' & !' 
          else
             write(unit = 78, fmt = '(A,A,A)') '     ,' // char(39), trim(nom(ie)) // repeat(' ', ii), char(39) // ' & !' 
          end if
          do k = 1, ie - 1
             if (nom(k) (1:inom(k)) == nom(ie) (1:inom(ie))) then
                write(*, *) 'ERROR: species already initialized: ', nom(ie)
                stop
             end if
          end do
       end do
       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '  !Number of each specie   '
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          ii = 4 - len(trim(nom(ie)))
          write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: ' // trim(nom(ie)) // &
          repeat(' ', ii) // '=' // cie
          if (trim(nom(ie)) == 'CO2') donthaveco2 = .false.
       end do
       ie = ie + 1
       write(cie, fmt = '(I3.3)') ie
       if (donthaveco2) write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: ' // 'CO2' // &
       repeat(' ', 1) // '=' // cie

       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '  !for memory allocattion: '
       write(unit = 78, fmt = '(A)') '  !This parameters are use for documentation only. '
       write(unit = 78, fmt = '(A)') '  !Use them in a program in substitution of numerical terms.'
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: src = 1 ! source term '
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: ddp = 2 ! dry deposition '
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: wdp = 3 ! wet deposition '
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: fdda = 4! four-dim assimilation '
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: offline = 5! ! off-line emissions: '
       write(unit = 78, fmt = '(A)') '                                  !' // &
       '=1, emission will be read from file'
       write(unit = 78, fmt = '(A)') '                                  ' // &
       '!=0, emission will be calculated during the model simulation (on-line emission)'
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: transport = 6! ! off-line emissions: '
       write(unit = 78, fmt = '(A)') '                                  !' // &
       '=1, yes'
       write(unit = 78, fmt = '(A)') '                                  ' // &
       '!=0, no transport'
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: on = 1'
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER :: off = 0'


       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '  ! spaction(specie,[1=source,2=drydep,3=wetdep,4=fdda,5=offline emission,6=transport]) ]) '
       write(unit = 78, fmt = '(A)') '  INTEGER,PARAMETER,DIMENSION(6,nspecies) :: spc_alloc=' // &
       'RESHAPE((/ &'

       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // csou(ie) // ' , ' // cdry(ie) // ' , ' // cwet(ie) // ' , ' &
             // cfdd(ie) // ' , ' // coff(ie) // ' , ' // ctra(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // csou(ie) // ' , ' // cdry(ie) // ' , ' // cwet(ie) // ' , ' &
             // cfdd(ie) // ' , ' // coff(ie) // ' , ' // ctra(ie) // '    & ! ' // trim(nom(ie)) // ' - ' // cie 
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /),(/6,nspecies/))'

       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A,I3.3,A)') '  INTEGER,PARAMETER,DIMENSION(nspecies) :: spc_uveq=(/ &'
       do ie = 1, nesp(2)
          if (ie /= nesp(2)) then
             write(unit = 78, fmt = '(A)') '   ' // uv_eq(ie) // ', & !' // trim(nom(ie))
          else
             write(unit = 78, fmt = '(A)') '   ' // uv_eq(ie) // '  & !' // trim(nom(ie))
          end if
       end do
       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       !Crossreference between Mechanism and Matrix
       write(unit = 78, fmt = '(A,I3.3,A)') '  INTEGER,PARAMETER,DIMENSION(5) :: spc_matEqv=(/ &'
       mat_eq = '0'

       do ie = 1, nesp(2)
          do ia = 1, 5
             write(caer, fmt = '(I1)') ia
             print *, 'LFR->ia,caer,ie,aer_eqv: ', ia, '|' // caer // '|', ie, '|' // aer_eq(ie) // '|'
             if (trim(aer_eq(ie)) == trim(caer)) then
                write(esp, fmt = '(I3.3)') ie
                print *, 'LFR->ie,esp:', ie, esp
                mat_eq(ia) = esp
                exit
             end if
          end do
       end do

       do ia = 1, 4
          write(unit = 78, fmt = '(A)') '   ' // mat_eq(ia) // ', &'
       end do
       write(unit = 78, fmt = '(A)') '   ' // mat_eq(5) // ' &'
       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '!     HENRYS LAW COEFFICIENTS'
       write(unit = 78, fmt = '(A)') '!     Henrys law coefficient'
       write(unit = 78, fmt = '(A)') '!     [KH298]=mole/(l atm)'
       write(unit = 78, fmt = '(A)') '!     Referencias em R. Sander (1999)'
       write(unit = 78, fmt = '(A)') '!     Compilation of Henry Law Constants '
       write(unit = 78, fmt = '(A)') '!     for Inorganic and Organic Species '
       write(unit = 78, fmt = '(A)') '!     of Potential Importance in '
       write(unit = 78, fmt = '(A)') '!     Environmental Chemistry (Version 3) '
       write(unit = 78, fmt = '(A)') '!     http://www.henrys-law.org '
       write(unit = 78, fmt = '(A)') '!     * indica artigos nao encontrados nesse endereco eletronico'
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: hstar=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // hstar(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // hstar(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /)'

       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '! [1] Noziere B. et al. The uptake of methyl vinyl ketone, methacrolein,'
       write(unit = 78, fmt = '(A)') '! and 2-methyl-3-butene-2-olonto sulfuric acid solutions,Journal of Physical'
       write(unit = 78, fmt = '(A)') '! Chemistry A, Vol.110, No.7, 2387-2395, 2006.'
       write(unit = 78, fmt = '(A)') '! [2] Abraham M. H. et al. Partition of compounds from gas to water and'
       write(unit = 78, fmt = '(A)') '! from gas to physiological saline at 310K: Linear free energy relationships,'
       write(unit = 78, fmt = '(A)') '! elsevier, 2006.'


       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: f0=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // f0(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // f0(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /)'

       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: difrat=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // difrat(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // difrat(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /)'


       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '!     DIFFUSION COEFFICIENTS'
       write(unit = 78, fmt = '(A)') '!     [DV]=cm2/s (assumed: 1/SQRT(molar mass) when not known)'
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: dvj=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // dvj(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // dvj(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /)'



       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '!     -DH/R (for temperature correction)'
       write(unit = 78, fmt = '(A)') '!     [-DH/R]=K'
       write(unit = 78, fmt = '(A)') '!     Referencias em R. Sander (1999)'
       write(unit = 78, fmt = '(A)') '!     Compilation of Henry Law Constants'
       write(unit = 78, fmt = '(A)') '!     for Inorganic and Organic Species '
       write(unit = 78, fmt = '(A)') '!     of Potential Importance in '
       write(unit = 78, fmt = '(A)') '!     Environmental Chemistry (Version 3)'
       write(unit = 78, fmt = '(A)') '!     http://www.henrys-law.org '
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: dhr=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // dhr(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // dhr(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /)'

       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '
       !----------
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: weight=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // cwei(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // cwei(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '
       !----------
       !----------srf
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: init_ajust=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // '1.0' // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // '1.0' // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '
       !----------
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: emiss_ajust=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // '1.0' // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // '1.0' // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '
       !----------



       write(unit = 78, fmt = '(A)') '!    ACID DISSOCIATION CONSTANT AT 298K '
       write(unit = 78, fmt = '(A)') '!     [mole/liter of liquid water]'
       write(unit = 78, fmt = '(A)') '!     Referencias: Barth et al. JGR 112, D13310 2007'
       write(unit = 78, fmt = '(A)') '!     Martell and Smith, 1976, Critical stability'
       write(unit = 78, fmt = '(A)') '!     vol1-4 Plenum Press New York'
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: ak0=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // ak0(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // ak0(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '   /)'
       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '

       write(unit = 78, fmt = '(A)') '!     Temperature correction factor for'
       write(unit = 78, fmt = '(A)') '!     acid dissociation constants'
       write(unit = 78, fmt = '(A)') '!     [K]'
       write(unit = 78, fmt = '(A)') '!     Referencias: Barth et al. JGR 112, D13310 2007'
       write(unit = 78, fmt = '(A)') '  REAL,PARAMETER,DIMENSION(nspecies) :: dak=(/&'
       do ie = 1, nesp(2)
          write(cie, fmt = '(I3.3)') ie
          if (ie < nesp(2)) then
             write(unit = 78, fmt = '(A)') &
             '    ' // dak(ie) // ' ,   & ! ' // trim(nom(ie)) // ' - ' // cie
          else
             write(unit = 78, fmt = '(A)') &
             '    ' // dak(ie) // '     & ! ' // trim(nom(ie)) // ' - ' // cie
          end if
       end do

       write(unit = 78, fmt = '(A)') '    /)'

       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') '  '



       ! change MP 29/01/2008: the following lines can be removed
       !READ(ifdin,*)
       !DO ie=nesp(2)+1,nesp(2)+nesp(3)
       ! READ(ifdin,'(a)')chdon
       !     write(imp,'(a)')chdon
       !  nblanc=nbmot
       !  CALL part(chdon,mot,imot,nmot)
       !  indaq(ie)=1
       !  nom(ie)(1:imot(1))=mot(1)(1:imot(1))
       !  inom(ie)=imot(1)
       !  DO k=1,ie-1
       !    IF (nom(k)(1:inom(k)) == nom(ie)(1:inom(ie))) THEN
       !      WRITE(*,*)'ERROR: species already initialized 1: ',nom(ie)
       !      STOP
       !    END IF
       !  END DO
       !END DO
       ! end change MP
       return
    end subroutine lectci

end module ModInit