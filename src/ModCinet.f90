module Modcinet
    implicit none

contains

    subroutine kinreac(mot, imot, nmot, ntuvonline)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization for kinetics of gas-phase reactions.
       !     The different routines associated to the kinetic laws are
       !     called.
       !------------------------------------------------------------------------

       !     -- INPUT VARIABLES

       !     -- INPUT/OUTPUT VARIABLES

       !     -- OUTPUT VARIABLES

       !     BP(.,NR): coefficient for kinetic rates.

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
       use ModFiccom, only: bp, av, ireversetab, nr, nrphot, ntabphot, tabphot &
                            , nb, ittb, indaqr, molec, iprecalc, ispebp, jer, s
       use ModGestion, only: reel, entier
       use ModGenerator, only: wk1, wk2, wk3, wkc9, wphot, wphot_tuvonline, wtroe &
                              , wtroe7, wtroe10, wcv, wrcfe, wspec, wtb, wfj, wpl
       implicit none

       integer, parameter :: nbmot = 100

       character (len = 500), intent(in out) :: mot(nbmot)
       integer, intent(in out) :: imot(nbmot)
       integer, intent(in out) :: nmot
       integer, intent(in out) :: ntuvonline

       double precision, dimension(ntabphotmax) :: b
       integer :: i, j, jj, it

       !     Arrhenius' law
       !     Gas-phase only
       i = 2

       do 
            if (mot(i) (1:4) == 'ARR1') then
               nb(nr) = 1
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call wk1(nr, bp(1, nr))
            else if (mot(i) (1:4) == 'ARR2') then
               nb(nr) = 2
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call wk2(nr, bp(1, nr), bp(2, nr))
            else if (mot(i) (1:4) == 'ARR3') then
               nb(nr) = 3
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call wk3(nr, bp(1, nr), bp(2, nr), bp(3, nr))
            else if (mot(i) (1:5) == 'ARRC2') then
               nb(nr) = 2
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               iprecalc(nr) = 2
            else if (mot(i) (1:5) == 'ARRC3') then
               nb(nr) = 3
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               iprecalc(nr) = 3

               ! Modification to allow
               ! Arrhenius combinations of the general form
               ! ARR3+ARR2+ARR2+ARR2
            else if (mot(i) (1:5) == 'ARRC9') then
               nb(nr) = 3   !tentativamente, revisar
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call reel(bp(5, nr), mot(i + 5), imot(i + 5))
               call reel(bp(6, nr), mot(i + 6), imot(i + 6))
               call reel(bp(7, nr), mot(i + 7), imot(i + 7))
               call reel(bp(8, nr), mot(i + 8), imot(i + 8))
               call reel(bp(9, nr), mot(i + 9), imot(i + 9))
               call wkc9 (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), bp(5, nr), &
               bp(6, nr), bp(7, nr), bp(8, nr), bp(9, nr))
               !end of modification

               !     PHOT: Photolysis.

            else if (mot(i) (1:4) == 'PHOT') then

               nb(nr) = 10
               nrphot = nrphot + 1
               if (ntabphot == 0) then
                  write(*, *) 'ERROR: tabulated angles not defined for photolysis'
                  write(*, *) 'Needs to be defined by SET TABULATION ...'
                  stop
               end if
               !     Read according increasing or decreasing sequence
               if (ireversetab == 0) then
                  do j = 1, ntabphot
                     call reel(bp(j, nr), mot(i + j), imot(i + j))
                     b(j) = bp(j, nr)
                  end do
               else
                  do j = 1, ntabphot
                     jj = ntabphot + 1 - j
                     call reel(bp(jj, nr), mot(i + j), imot(i + j))
                     b(jj) = bp(jj, nr)
                  end do
               end if


               if (ntuvonline == 0) then
                  it = 0
                  call wphot(nr, ntabphot, b, tabphot, it)
               else
                  call wphot_tuvonline(nr)
               end if
               !     TROE: TROE/Fall-off for the general case.
            else if (mot(i) (1:5) == 'TROE4') then
               nb(nr) = 4
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call wtroe (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), 0.6d0)

            else if (mot(i) (1:5) == 'TROE5') then
               nb(nr) = 4
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call reel(bp(5, nr), mot(i + 5), imot(i + 5))
               call wtroe (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), bp(5, nr))

            else if (mot(i) (1:5) == 'TROE7') then
               nb(nr) = 4
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call reel(bp(5, nr), mot(i + 5), imot(i + 5))
               call reel(bp(6, nr), mot(i + 6), imot(i + 6))
               call reel(bp(7, nr), mot(i + 7), imot(i + 7))
               call wtroe7 (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), bp(5, nr), &
               bp(6, nr), bp(7, nr))

               !     TROE12: TROE/Fall-off for MOCA.
            else if (mot(i) (1:6) == 'TROE10') then
               nb(nr) = 4
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call reel(bp(5, nr), mot(i + 5), imot(i + 5))
               call reel(bp(6, nr), mot(i + 6), imot(i + 6))
               call reel(bp(7, nr), mot(i + 7), imot(i + 7))
               call reel(bp(8, nr), mot(i + 8), imot(i + 8))
               call reel(bp(9, nr), mot(i + 9), imot(i + 9))
               call reel(bp(10, nr), mot(i + 10), imot(i + 10))
               call wtroe10 (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), bp(5, nr), &
               bp(6, nr), bp(7, nr), bp(8, nr), bp(9, nr), bp(10, nr))

               !     CVAR/MOCA: temperature-dependent stoichiometry.
            else if (mot(i) (1:5) == 'CVAR') then
               nb(nr) = 5
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call reel(bp(5, nr), mot(i + 5), imot(i + 5))
               call reel(bp(6, nr), mot(i + 6), imot(i + 6))
               call reel(bp(7, nr), mot(i + 7), imot(i + 7))
               call wcv (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), bp(5, nr), &
               bp(6, nr), bp(7, nr))

               !     RCFE: Reactions Calculated From Equilibria.
            else if (mot(i) (1:4) == 'RCFE') then
               nb(nr) = 8
               call reel(bp(1, nr), mot(i + 1), imot(i + 1))
               call reel(bp(2, nr), mot(i + 2), imot(i + 2))
               call reel(bp(3, nr), mot(i + 3), imot(i + 3))
               call reel(bp(4, nr), mot(i + 4), imot(i + 4))
               call reel(bp(5, nr), mot(i + 5), imot(i + 5))
               call reel(bp(6, nr), mot(i + 6), imot(i + 6))
               call wrcfe (nr, bp(1, nr), bp(2, nr), bp(3, nr), bp(4, nr), bp(5, nr), bp(6, nr))

               !     SPEC: specific reactions.
            else if (mot(i) (1:4) == 'SPEC') then
               nb(nr) = 5
               call entier(ispebp(nr), mot(i + 1), imot(i + 1))
               call wspec (nr, ispebp(nr))

               !     EXTRA: specific reaction with corrected factors
               !     O3 -> 2. OH with corrected photolysis
            else if (mot(i) (1:5) == 'EXTRA') then
               nb(nr) = 10
               if (ntabphot == 0) then
                  write(*, *) 'ERROR: tabulation not given for photolysis'
                  write(*, *) 'Needs to be defined by SET TABULATION ...'
                  stop
               end if
               !     Read according increasing or decreasing sequence
               if (ireversetab == 0) then
                  do j = 1, ntabphot
                     call reel(bp(j, nr), mot(i + j), imot(i + j))
                     b(j) = bp(j, nr)
                  end do
               else
                  do j = 1, ntabphot
                     jj = ntabphot + 1 - j
                     call reel(bp(jj, nr), mot(i + j), imot(i + j))
                     b(jj) = bp(jj, nr)
                  end do
               end if
               it = 1
               call wphot(nr, ntabphot, b, tabphot, 1)

               !     Third body.
            else if (mot(i) (1:2) == 'TB') then
               if (mot(i + 1) (1:1) == 'M') then
                  ittb(nr) = 1
               else if (mot(i + 1) (1:2) == 'O2') then
                  ittb(nr) = 2
               else if (mot(i + 1) (1:2) == 'N2') then
                  ittb(nr) = 3
               else if (mot(i + 1) (1:3) == 'H2O') then
                  ittb(nr) = 4
               else if ((mot(i + 1) (1:2) == 'H2') .and. (imot(i + 1) == 2)) then
                  ittb(nr) = 5
               else
                  write(*, *) 'ERROR: syntax for Third Body'
                  write(*, *) 'M, O2, N2, H20 or H2 expected.'
                  stop
               end if
               i = i + 2
            else
               write(*, *) 'ERROR: unknown syntax for KINETIC definition ', mot(2)
               stop
            end if

            !     Modification BS/KS 21/05/2002
            !     Case of a third body reaction: need for kinetics.

            if ((ittb(nr) /= 0) .and. (nb(nr) == 0)) then
                cycle
            else
                exit
            end if
        end do

       if (ittb(nr) /= 0) call wtb(nr, ittb(nr))

       !     Update the chemical production term and the Jacobian matrix.
       call wfj(s, nr, jer)

       !     Update the production and loss terms (P-Lc formulation)
       call wpl(s, nr, jer)

       !     Conversion mol/l -> molec/cm3.
       if (indaqr(nr) == 1) then
          if (molec(nr) == 2) bp(1, nr) = bp(1, nr) * 1.0d3 / av
          if (molec(nr) == 3) bp(1, nr) = bp(1, nr) * 1.0d6 / av**2
       end if

       return
    end subroutine kinreac

    subroutine kindis(mot, imot, nmot)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization of kinetics for ionic dissociations.

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
       use ModFiccom, only: nequil, xk1, xk2, nom
       use ModGestion, only: reel
       implicit none

       integer, parameter :: nbmot = 100

       character (len = 500), intent(in out) :: mot(nbmot)
       integer, intent(in out) :: imot(nbmot)
       integer, intent(in out) :: nmot

       if (mot(2) (1:5) == 'ARRC2') then
          call reel(xk1(nequil), mot(3), imot(3))
          call reel(xk2(nequil), mot(4), imot(4))
       else
          write(*, *) nom(2) (1:4), ' ', 'ERROR: syntax.'
          stop
       end if

    end subroutine kindis

    subroutine kinhenry(mot, imot, nmot)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization for kinetics of Henry's equilibrium.

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
       use ModFiccom, only: nr, nb, jer, bp, iprecalc, rmol, nom, inom
       use ModGestion, only: reel
       implicit none

       integer, parameter :: nbmot = 100

       character (len = 500), intent(in out) :: mot(nbmot)
       integer, intent(in out) :: imot(nbmot)
       integer, intent(in out) :: nmot

       integer :: i, ilettre, ncomp, nn, ilet

       !     Reaction G ->L
       nr = nr - 1

       nb(nr) = 6
       i = jer(1, nr)

       if (rmol(i) == 0.) then
          ilettre = 1

          do while (ilettre <= inom(i))
             ncomp = 0
             if (nom(i) (ilettre:ilettre) == 'H') then
                ncomp = 1
             else if (nom(i) (ilettre:ilettre) == 'C') then
                ncomp = 12
             else if (nom(i) (ilettre:ilettre) == 'O') then
                ncomp = 16
             else if (nom(i) (ilettre:ilettre) == 'N') then
                ncomp = 14
             else if (nom(i) (ilettre:ilettre) == 'S') then
                ncomp = 32
             else
                write(*, *) 'ERROR: species ', nom(i), ' unknown molar mass.'
                stop
             end if
             nn = 1
             ilet = ilettre + 1
             if (inom(i) >= ilet) then
                if (nom(i) (ilet:ilet) == '2') nn = 2
                if (nom(i) (ilet:ilet) == '3') nn = 3
                if (nom(i) (ilet:ilet) == '4') nn = 4
                if (nom(i) (ilet:ilet) == '5') nn = 5
                if (nom(i) (ilet:ilet) == '6') nn = 6
             end if
             rmol(i) = rmol(i) + ncomp * nn
             if (nn > 1) ilettre = ilettre + 1
             ilettre = ilettre + 1
          end do
          write(*, *) 'Computed molar mass for ', nom(i), '=', rmol(i)
       else
          write(*, *) 'Read molar mass for ', nom(i), '=', rmol(i)
       end if

       if (rmol(i) == 0.) then
          write(*, *) 'ERROR: Species ', nom(i), ' unknown molar mass.'
          stop
       end if

       !     Reaction L->G

       nr = nr + 1
       nb(nr) = 7
       if (mot(2) (1:5) == 'ARRC2') then
          call reel(bp(1, nr), mot(3), imot(3))
          call reel(bp(2, nr), mot(4), imot(4))
          call reel(bp(3, nr), mot(5), imot(5))
          iprecalc(nr) = 2
       else if (mot(2) (1:4) == 'ARR1') then
          call reel(bp(1, nr), mot(3), imot(3))
          bp(2, nr) = 0.d0
       else if (mot(2) (1:4) == 'ARR2') then
          call reel(bp(1, nr), mot(3), imot(3))
          call reel(bp(2, nr), mot(4), imot(4))
       else
          write(*, *) 'ERROR: syntax kinetics Henry ', nr - 1
          stop
       end if

       return
    end subroutine kinhenry

    subroutine creac(mot, imot, nmot)

       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization for gas-phase reactions.

       !------------------------------------------------------------------------

       !     -- INPUT VARIABLES

       !     -- INPUT/OUTPUT VARIABLES

       !     NR: number of reactions.
       !     S: stoichiometric matrix

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
       use ModFiccom, only: nrmax, nom, inom, nesp, indaqr, molec, s, jer, nr &
                            , indaq, indaqr
       use ModGenerator, only: ww, dw
       use ModGestion, only: reel
       implicit none

       integer, parameter :: nbmot = 100

       character (len = 500), intent(inout) :: mot(nbmot)
       integer, intent(in out) :: imot(nbmot)
       integer, intent(in out) :: nmot

       double precision :: stoieir
       character (len = 12) :: nam(5)
       integer :: inam(5)
       integer :: icurseur, je, ie, indic, i

       nr = nr + 1

       if (nr > nrmax) then
          write(*, *) 'ERROR: bad dimension for nr>nrmax.'
          stop
       end if

       !     Check reactants

       molec(nr) = 1
       nam(1) (1:imot(1)) = mot(1) (1:imot(1))
       inam(1) = imot(1)
       icurseur = 1
       do
            icurseur = icurseur + 1
            if (mot(icurseur) (1:imot(icurseur)) == '+') then
               icurseur = icurseur + 1
               molec(nr) = molec(nr) + 1
               nam(molec(nr)) (1:imot(icurseur)) = mot(icurseur) (1:imot(icurseur))
               inam(molec(nr)) = imot(icurseur)
               cycle
            else
               exit
            end if
       end do

       if (molec(nr) > 3) then
          write(*, *) 'ERROR: more than 3 reactants.'
          stop
       end if

       do je = 1, molec(nr)
          indic = 0
          do ie = 1, nesp(1)
             if (inam(je) == inom(ie)) then
                if (nam(je) (1:inam(je)) == nom(ie) (1:inom(ie))) then
                   jer(je, nr) = ie
                   s(ie, nr) = s(ie, nr) -1.d0
                   indic = 1
                end if
             end if
          end do

          if (indic == 0) then
             write(6, *) 'ERROR: the following reactant is unknown ', nam(je) (1:inam(je))
             stop
          end if
       end do

       call ww(nr, molec(nr), jer(1, nr), jer(2, nr), jer(3, nr))
       call dw(nr, molec(nr), jer(1, nr), jer(2, nr), jer(3, nr))

       !     Check products

       if ((mot(icurseur) (1:imot(icurseur)) /= '>') .and. &
       (mot(icurseur) (1:imot(icurseur)) /= '->')) then
          write(*, *) 'ERROR: > or -> expected '
          stop
       end if
       icurseur = icurseur + 1

       200 stoieir = 1.d0
       if (icurseur < nmot) then
          if ((mot(icurseur + 1) (1:imot(icurseur + 1)) /= '+') .and. &
          (mot(icurseur + 1) (1:imot(icurseur + 1)) /= '-')) then
             call reel(stoieir, mot(icurseur), imot(icurseur))
             icurseur = icurseur + 1
          end if
       end if
       nam(1) (1:imot(icurseur)) = mot(icurseur) (1:imot(icurseur))

       indic = 0
       do ie = 1, nesp(1)
          if ((imot(icurseur) == inom(ie)) .and. &
          (nam(1) (1:imot(icurseur)) == nom(ie) (1:inom(ie)))) then
             s(ie, nr) = s(ie, nr) + stoieir
             indic = 1
          end if
       end do
       if (indic == 0) then
          write(6, *) 'WARNING: product unknown ', nam(1) (1:imot(icurseur))
       end if

       icurseur = icurseur + 1
       if (icurseur <= nmot) then
          if ((mot(icurseur) (1:imot(icurseur)) == '+') .or. &
          (mot(icurseur) (1:imot(icurseur)) == '-')) then
             icurseur = icurseur + 1
             go to 200
          else
             write(*, *) 'ERROR: + or - expected.'
             stop
          end if
       end if

       !     Check gas/liquid

       do i = 1, nesp(1)
          if ((s(i, nr) /= 0.) .and. (indaq(i) /= indaqr(nr))) then
             write(*, *) 'ERROR: the phases are not coherent'
             write(*, *) 'reaction ', nr, ' species ', nom(i)
             stop
          end if
       end do

    end subroutine creac

    subroutine cdis(mot, imot, nmot)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization for ionic dissociations.

       !------------------------------------------------------------------------

       !     -- INPUT VARIABLES

       !     -- INPUT/OUTPUT VARIABLES

       !     NR: number of reactions.
       !     S: stoichiometric matrix

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
       use ModFiccom, only: nom, nequil, iesp, nneq, nesp, inom, indaq, seqion
       use ModGestion, only: reel
       implicit none

       integer, parameter :: nbmot = 100

       character (len = 500), intent(in out) :: mot(nbmot)
       integer, intent(in  out) :: imot(nbmot)
       integer, intent(in out) :: nmot

       character (len = 12) :: nam(5)
       integer :: inam(5)
       integer :: icurseur, ie, indic
       double precision :: stoeq

       nequil = nequil + 1
       if (nequil > nionx) then
          write(*, *) 'ERROR: bad dimension for nionx=', nionx
          stop
       end if

       !     Reactant.
       nam(1) (1:imot(1)) = mot(1) (1:imot(1))
       inam(1) = imot(1)

       indic = 0
       do ie = 1, nesp(1)
          if ((inam(1) == inom(ie)) .and. &
          (nam(1) (1:inam(1)) == nom(ie) (1:inom(ie)))) then
             indic = ie
             iesp(nequil, 1) = ie
          end if
       end do
       if (indic == 0) then
          write(6, *) 'ERROR: unknown species ', nam(1)
          stop
       end if
       if (indaq(indic) == 0) then
          write(6, *) 'ERROR: this is a gaseous species ', nam(1)
          stop
       end if

       icurseur = 1
       icurseur = icurseur + 1
       if (mot(icurseur) (1:imot(icurseur)) /= '=') then
          write(*, *) 'ERROR: syntax = expected.'
          stop
       end if

       icurseur = icurseur + 1

       !     Products.
       nneq(nequil) = 2
       inam(nneq(nequil)) = imot(icurseur)
       nam(nneq(nequil)) (1:imot(icurseur)) = mot(icurseur) (1:imot(icurseur))

       200 stoeq = 1.d0

       !     BS 2003: to be checked !
       if (icurseur < nmot) then
          if ((mot(icurseur + 1) (1:imot(icurseur + 1)) /= '+') .and. &
          (mot(icurseur + 1) (1:imot(icurseur + 1)) /= '-')) then
             call reel(stoeq, mot(icurseur), imot(icurseur))
             icurseur = icurseur + 1
          end if
       end if

       nam(nneq(nequil)) (1:imot(icurseur)) = mot(icurseur) (1:imot(icurseur))
       inam(nneq(nequil)) = imot(icurseur)


       write(*, *) 'Ionic species: ', nam(nneq(nequil)) (1:inam(nneq(nequil)))

       indic = 0
       do ie = 1, nesp(1)
          if ((inam(nneq(nequil)) == inom(ie)) .and. &
          (nam(nneq(nequil)) (1:imot(icurseur)) == nom(ie) (1:inom(ie)))) &
          then
             seqion(ie, nequil) = seqion(ie, nequil) + stoeq
             iesp(nequil, nneq(nequil)) = ie
             indic = 1
             if (indaq(ie) == 0) then
                write(*, *) 'ERROR: gaseous species ', nom(ie), ' in dissociation ', nequil
                stop
             end if
          end if
       end do
       if (indic == 0) then
          write(6, *) 'ERROR: unknown species in dissociation ', nam(nneq(nequil))
          stop
       end if

       icurseur = icurseur + 1

       if (icurseur <= nmot) then
          if ((mot(icurseur) (1:imot(icurseur)) == '+') .or. &
          (mot(icurseur) (1:imot(icurseur)) == '-')) then
             icurseur = icurseur + 1
             nneq(nequil) = nneq(nequil) + 1
             go to 200
          else
             write(*, *) 'ERROR: syntax + expected.'
             stop
          end if
       end if

       !     Check lumping for more than 4 ions.
       if (nneq(nequil) > 3) then
          write(*, *) 'ERROR: more than 2 products in dissociation. See routine "lump"'
          stop
       end if

    end subroutine cdis

    subroutine chenry(mot, imot, nmot, ieq)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization for Henry's equilibrium.

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
       use ModFiccom, only: nr, nrmax, nom, molec, inom, nesp, indaq, indaqr, s &
                           , jer, jep, indaq, indaqr, iheq, ihreac
       implicit none

       integer, parameter :: nbmot = 100

       character (len = 500), intent(in) :: mot(nbmot)
       integer, intent(in) :: imot(nbmot)
       integer, intent(in out) :: nmot
       integer, intent(in) :: ieq

       character (len = 12) :: nam(5)

       integer :: inam(5)
       integer :: ie, indic, indp, indr, j1, j2, i

       !     Henry gas --> liq
       !     Creation of Henry liq --> gas.

       nr = nr + 1

       if (nr > nrmax) then
          write(*, *) 'ERROR: Bad dimension for nr>nrmax'
          stop
       end if

       !     Reactant.

       molec(nr) = 1

       nam(1) (1:imot(1)) = mot(1) (1:imot(1))
       inam(1) = imot(1)
       indic = 0
       do ie = 1, nesp(1)
          if (inam(1) == inom(ie)) then
             if (nam(1) (1:inam(1)) == nom(ie) (1:inom(ie))) then
                indr = ie
                jer(1, nr) = ie
                s(ie, nr) = -1.d0
                indic = 1
             end if
          end if
       end do

       if (indic == 0) then
          write(6, *) 'ERROR: unknown species ', nam(1)
          stop
       end if

       !     Product.

       nam(1) (1:imot(3)) = mot(3) (1:imot(3))
       inam(1) = imot(3)
       indic = 0
       do ie = 1, nesp(1)
          if (inam(1) == inom(ie)) then
             if (nam(1) (1:inam(1)) == nom(ie) (1:inom(ie))) then
                indp = ie
                s(ie, nr) = 1.d0
                indic = 1
             end if
          end if
       end do

       if (indic == 0) then
          write(6, *) 'ERROR: unknown species ', nam(1)
          stop
       end if
       indaqr(nr) = indaq(jer(1, nr))

       !     Creation reaction L->G

       nr = nr + 1
       if (nr > nrmax) then
          write(*, *) 'ERROR: Bad dimension for nr>nrmax.'
          stop
       end if

       molec(nr) = 1

       s(indp, nr) = -1.
       s(indr, nr) = 1.d0
       jep(nr) = indr
       jer(1, nr) = indp

       !     Lumping.

       iheq(nr - 1) = ieq
       iheq(nr) = ieq
       j1 = jer(1, nr)
       j2 = jer(1, nr - 1)
       ihreac(j1) = nr
       ihreac(j2) = nr - 1

       if (indaq(jer(1, nr)) /= 1) then
          write(*, *) 'ERROR: phases are not coherent for'
          write(*, *) 'Henry ', nr
          write(*, *) 'Species ', nom(i), ' has to be aqueous.'
          stop
       end if
       if (indaq(jep(nr)) /= 0) then
          write(*, *) 'ERROR phases are not coherent for'
          write(*, *) 'Henry ', nr
          write(*, *) 'Species ', nom(i), ' has to be gaseous.'
          stop
       end if
       indaqr(nr) = 1

       return
    end subroutine chenry

    subroutine initphase(ntuvonline, filename)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Initialization for multiphase description.
       !     IPHASE = 1 : multiphase description.
       !     IPHASE = 2 : gas-phase.
       !     IPHASE = 3 : aqueous-phase.

       !     IRMONODI: index reaction in one phase --> multiphase.

       !     Notice that Henry's reactions are only to be taken into account for
       !     the multiphase case.
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
       use ModFiccom, only: indaqr, nb, nrmol1, nrmol2, nrmol3, narr1, narr2, narr3 &
                          , narr4, narr5, narr6, narr7, narr8, imolec1, imolec2 &
                          , imolec3, iarr1, iarr2, iarr3, iarr4, iarr5, iarr6, iarr7 &
                          , iarr8, molec, nr, nrp, indaq, ittb, irmonodi, naq2, inaq2 &
                          , nthird, indthird
       implicit none

       integer, intent(in) :: ntuvonline
       character(len = *), intent(in) :: filename
       integer :: ir, ip

       do ir = 1, nr

          ip = indaqr(ir) + 2
          nrp(1) = nrp(1) + 1
          if ((nb(ir) /= 6) .and. (nb(ir) /= 7)) then
             nrp(ip) = nrp(ip) + 1
             irmonodi(nrp(ip), ip) = ir
          end if

          if (molec(ir) == 1) then
             nrmol1(1) = nrmol1(1) + 1
             imolec1(nrmol1(1), 1) = nrp(1)
             if ((nb(ir) /= 6) .and. (nb(ir) /= 7)) then
                nrmol1(ip) = nrmol1(ip) + 1
                imolec1(nrmol1(ip), ip) = ir
             end if

          else if (molec(ir) == 2) then
             nrmol2(1) = nrmol2(1) + 1
             imolec2(nrmol2(1), 1) = nrp(1)
             nrmol2(ip) = nrmol2(ip) + 1
             imolec2(nrmol2(ip), ip) = ir

          else if (molec(ir) == 3) then
             nrmol3(1) = nrmol3(1) + 1
             imolec3(nrmol3(1), 1) = nrp(1)
             nrmol3(ip) = nrmol3(ip) + 1
             imolec3(nrmol3(ip), ip) = ir
          end if

          if (nb(ir) == 1) then
             narr1(1) = narr1(1) + 1
             iarr1(narr1(1), 1) = nrp(1)
             narr1(ip) = narr1(ip) + 1
             iarr1(narr1(ip), ip) = ir
          else if (nb(ir) == 2) then
             narr2(1) = narr2(1) + 1
             iarr2(narr2(1), 1) = nrp(1)
             narr2(ip) = narr2(ip) + 1
             iarr2(narr2(ip), ip) = ir
          else if (nb(ir) == 3) then
             narr3(1) = narr3(1) + 1
             iarr3(narr3(1), 1) = nrp(1)
             narr3(ip) = narr3(ip) + 1
             iarr3(narr3(ip), ip) = ir
          else if (nb(ir) == 4) then
             narr4(1) = narr4(1) + 1
             iarr4(narr4(1), 1) = nrp(1)
             narr4(ip) = narr4(ip) + 1
             iarr4(narr4(ip), ip) = ir
          else if (nb(ir) == 5) then
             narr5(1) = narr5(1) + 1
             iarr5(narr5(1), 1) = nrp(1)
             narr5(ip) = narr5(ip) + 1
             iarr5(narr5(ip), ip) = ir
          else if (nb(ir) == 6) then
             narr6(1) = narr6(1) + 1
             iarr6(narr6(1), 1) = nrp(1)
          else if (nb(ir) == 7) then
             narr7(1) = narr7(1) + 1
             iarr7(narr7(1), 1) = nrp(1)
          else if (nb(ir) == 8) then
             narr8(1) = narr8(1) + 1
             iarr8(narr8(1), 1) = nrp(1)
             narr8(ip) = narr8(ip) + 1
             iarr8(narr8(ip), ip) = ir
          end if

          if (indaqr(ir) == 1) then
             if (molec(ir) == 2) then
                naq2(1) = naq2(1) + 1
                inaq2(naq2(1), 1) = ir
                naq2(3) = naq2(3) + 1
                inaq2(naq2(3), 1) = ir
             end if
          end if

          if (ittb(ir) /= 0) then
             nthird(1) = nthird(1) + 1
             indthird(nthird(1), 1) = nrp(1)
             nthird(ip) = nthird(ip) + 1
             indthird(nthird(ip), ip) = ir
          end if

       end do

       !     write(*,*)' Total number of reactions         =',nrp(1)
       write(*, *) ' Number of gas-phase reactions     =', nrp(2)
       write(*, *) ' Number of aqueous-phase reactions =', nrp(3)
       write(*, *) ' Nbr of Henry reversible reactions =', nrp(1) - nrp(2) - nrp(3)
       write(*, *) ' Third body reactions              =', nthird(1)
       write(*, *) '#############################################'

       write(unit = 78, fmt = '(A,I3.3,A)') '  INTEGER,PARAMETER :: nr   =', nrp(2), '!Number of gas-phase reactions'
       write(unit = 78, fmt = '(A,I3.3,A)') '  INTEGER,PARAMETER :: nrt  =', nrp(1), '!Total Number of reactions'
       write(unit = 78, fmt = '(A,I3.3,A)') '  INTEGER,PARAMETER :: nrh2o=', nrp(3), '!Number of aqueous-phase reactions'

       write(unit = 78, fmt = '(A)') '  '

       !IF Fastjx then call the map
       print *, 'Cinet->Ntuvonline=', ntuvonline

       !IF(ntuvonline==1) CALL JxSpack(fileName,nrphot,78)
       call jxspack(filename, ntuvonline, 78)


       write(unit = 78, fmt = '(A)') '  '
       write(unit = 78, fmt = '(A)') 'END MODULE chem1_list'

       close(unit = 78)

    end subroutine initphase

    subroutine jxspack(filename, ntuvonline, noutfile)
       implicit none
       character(len = *), intent(in) :: filename !Input file name
       integer :: nreactions !Total of reactions 
       integer, intent(in) :: noutfile !File to write the output reactions
       integer, intent(in) :: ntuvonline
       integer, parameter :: nfile = 22 !#file to read

       integer, parameter :: maxjcomb = 5 !Maximum of J reactions combinations
       integer, parameter :: sizeofspeciename = 7 !Size of specie name
       integer, parameter :: maxnreactions = 200!Maximum of reactions
       character(len = 80) :: header
       character(len = 20) :: title
       character(len = 3) :: sp, ss
       character(len = 80) :: fmt_txt
       integer :: nofreaction !Number od reaction
       integer, dimension(maxnreactions) :: specprod !Total of species production by reaction
       integer :: i, n

       type tf
          double precision :: factor !Factor of production
          character(len = sizeofspeciename) :: specie !Name of specie    
       end type tf
       type(tf), dimension(maxnreactions, maxjcomb) :: spcdata


       if (ntuvonline > 0) then
          open(unit = nfile, file = filename)
          !Reading the 5 first lines - header
          do i = 1, 5
             read(nfile, fmt = '(A80)') header
          end do
          !Reading the reactions and species' production factors
          spcdata(:, :)%specie = 'NONE'
          do i = 1, maxnreactions
             read(nfile, fmt = *, end = 110) nofreaction, specprod(i), (spcdata(i, n), n = 1, specprod(i))
             !WRITE (*,*) i,nOfReaction,specProd(i),(spcData(i,n),n=1,specProd(i))
             !pause
          end do
       end if
       110 continue

       if (ntuvonline == 0) then
          nreactions = i - 1
          close(unit = 22)
          title = 'FAST-JX'
       elseif(ntuvonline == 0) then
          nreactions = 1
          spcdata(:, :)%specie = 'NONE'
          specprod(:) = 0
          spcdata(:, :)%factor = 0.
          title = 'LUT'
       elseif(ntuvonline == 2) then
          nreactions = 1
          spcdata(:, :)%specie = 'NONE'
          specprod(:) = 0
          spcdata(:, :)%factor = 0.
          title = 'FAST-TUV'
       else
          print *, 'Invalid photojmethod - ntuvonline=', ntuvonline
          stop 100

       endif

       !not oficial bug fix
       select case(ntuvonline)
          case(0)
             title = 'LUT'
          case(1)
             title = 'FAST-JX'
          case(2)
             title = 'FAST-TUV'
       end select

       !Writing the adapt fastJX to Spack routine

       !First part: An array with number of species from Jx which makes a reaction
       write(noutfile, fmt = '(A)') ' !------------------------------------------------------------------------------  '
       write(noutfile, fmt = '(A)') '  '
       write(noutfile, fmt = '(A)') ' ! Photolysis Rate Calculation: method used (LUT=look_up_table, FAST-JX= on-line)'
       write(noutfile, fmt = '(A)') '  '

       write(noutfile, fmt = '(A,A,A,A)') '  CHARACTER(LEN=10),PARAMETER :: PhotojMethod= ', char(39), trim(title), char(39)
       write(noutfile, fmt = '(A,i5,A,i5)') '  INTEGER,PARAMETER :: maxJcomb=', maxjcomb, ', nr_photo=', nreactions
       write(noutfile, fmt = '(A)') '  INTEGER,PARAMETER,DIMENSION(nr_photo) :: nfactors=(/ &'
       do i = 1, nreactions - 1
          write(noutfile, fmt = '(A,I2,A,I3)') '                           ', specprod(i), ', &!', i
       end do
       write(noutfile, fmt = '(A,I2,A,I3)') '                           ', specprod(nreactions), '/)!', i
       write(noutfile, fmt = '(A)') ' '
       write(sp, fmt = '(I2.2)') maxjcomb

       !Second Part: An array with multiplication factor for each species from Jx
       write(noutfile, fmt = '(A,I2.2,A)') '  DOUBLE PRECISION,PARAMETER,DIMENSION(maxJcomb,nr_photo) :: factor=RESHAPE((/ &'
       do i = 1, nreactions - 1
          write(noutfile, fmt = '(A,$)') '                                               '
          do n = 1, maxjcomb
             write(noutfile, fmt = '(D10.3,",",$)') spcdata(i, n)%factor
          end do
          write(noutfile, fmt = '(A)') ' &'
       end do
       write(noutfile, fmt = '(A,$)') '                                               '
       do n = 1, maxjcomb - 1
          write(noutfile, fmt = '(D10.3,",",$)') spcdata(i, n)%factor
       end do
       write(noutfile, fmt = '(D10.3,A)') spcdata(i, maxjcomb)%factor, '/),(/maxJcomb,nr_photo/))'

       write(ss, fmt = '(I2.2)') sizeofspeciename


       write(noutfile, fmt = '(A,I2.2,A)') '  CHARACTER(LEN=' // ss  &
             // '),PARAMETER,DIMENSION(maxJcomb,nr_photo) :: JReactionComp=RESHAPE((/ &'
       do i = 1, nreactions - 1
          write(noutfile, fmt = '(A,$)') '                                               '
          do n = 1, maxjcomb
             write(noutfile, fmt = '(A,",",$)') '"' // spcdata(i, n)%specie // '"'
          end do
          write(noutfile, fmt = '(A)') ' &'
       end do
       write(noutfile, fmt = '(A,$)') '                                               '
       do n = 1, maxjcomb - 1
          write(noutfile, fmt = '(A,",",$)') '"' // spcdata(i, n)%specie // '"'
       end do
       write(noutfile, fmt = '(A,A)') '"' // spcdata(i, maxjcomb)%specie, '"/),(/maxJcomb,nr_photo/))'

       write(noutfile, fmt = '(A)') ''


    end subroutine jxspack


end module ModCinet