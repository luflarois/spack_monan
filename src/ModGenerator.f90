module ModGenerator
    use ModParametre
    implicit none

!------------------------------------------------------------------------
!     Subroutines for automatic generation of fortran files
!     needed for integration of gas-phase chemistry.
!     Authors: Bruno Sportisse and Pierre Plion, CEREA/ENPC
!     Date: February 2003.
!------------------------------------------------------------------------

contains

    subroutine wk1(nr, k)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetic rates case.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick
       implicit none

       integer, intent(in out) :: nr
       double precision, intent(in out) :: k

       write(nfick, 10) nr, k

       10 format(6x, 'rk(ijk,', i3, ') = ', d23.16)

    end subroutine wk1

    subroutine wk2(nr, k, tact)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetic rates case ARR2:
       !     k(T) = b1 * exp(-b2/T)

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick, petit
       implicit none

       integer, intent(in out) :: nr
       double precision, intent(out) :: k
       double precision, intent(in out) :: tact


       !Modification, added (ABS(k) < petit)
       if ((abs(tact) < petit) .or. (abs(k) < petit)) then
          !end of modification
          call wk1(nr, k)
       else
          k = log(k)
          write(nfick, 10) nr, k, tact
       end if

       10 format(6x, 'rk(ijk,', i3, ') =  DEXP(', d23.16, 3x, '&', / 8x, &
       ' - ( ', d23.16, ' )/temp(ijk))')

    end subroutine wk2
    !------------------------------------------------------------------------

    subroutine wk3(nr, k, expt, tact)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetic rates case.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick, petit
       implicit none

       integer, intent(in out) :: nr
       double precision, intent(out) :: k
       double precision, intent(in out) :: expt
       double precision, intent(in out) :: tact

       if (abs(expt) < petit) then
          call wk2(nr, k, tact)
       else
          k = log(k)
          if (tact > petit) then
             write(nfick, 10) nr, k, expt, tact
          else if (tact < -petit) then
             write(nfick, 11) nr, k, expt, (-tact)
          else
             write(nfick, 12) nr, k, expt
          end if
       end if


       10 format(6x, 'rk(ijk,', i3, ') =  DEXP(', d23.16, 3x, '&',  /, &
       8x, '+ (', d23.16, ' * LOG(temp(ijk)))', 3x, '&',  /, 8x, '- ', d23.16, '/temp(ijk))')


       11 format(6x, 'rk(ijk,', i3, ') =  DEXP(', d23.16, 3x, '&',  /, &
       8x, '+ (', d23.16, ' * LOG(temp(ijk)))', 3x, '&',  /, 8x, '+ ', d23.16, '/temp(ijk))')
       12 format(6x, 'rk(ijk,', i3, ') =  DEXP(', d23.16, 3x, '&',  /, &
       8x, '+ (', d23.16, ' * LOG(temp(ijk))) ) ')

    end subroutine wk3

    subroutine wkc9(nr, k1, expt, tact1, k2, tact2, k3, tact3, k4, tact4)
       !------------------------------------------------------------------------
       !     -- DESCRIPTION
       ! Arrhenius combinations
       ! of the general form ARR3+ARR2+ARR2+ARR2
       ! 9 parameters, 3 for the ARR3 and 2 for each following ARR2 (3*2=6):
       ! k1*(T^expt)*exp(-tact1/T)+k2*exp(-tact2/T)+k3*exp(-tact3/T)+k4*exp(-tact4/T)
       ! 6 combinations implemented (formats are more general):
       !-----------------------------
       ! arr2+0+0+arr1 -> k1,0,tact1,0,0,0,0,k4,0
       ! arr2+arr2+0+0 -> k1,0,tact1,k2,tact2,0,0,0,0
       ! arr2+arr2+0+arr1 -> k1,0,tact1,k2,tact2,0,0,k4,0
       ! arr2+arr2+arr2+0 -> k1,0,tact1,k2,tact2,k3,tact3,0,0
       ! arr2+arr2+arr2+arr2 -> k1,0,tact1,k2,tact2,k3,tact3,k4,tact4
       ! arr3+arr2+arr2+arr1 -> k1,expt,tact1,k2,tact2,k3,tact3,k4,tact4
       ! ARR3 = ARR2 with T^0; ARR2 = ARR1 with e^(-0/T)
       ! olhio! aindã sem chequeo de erros!!!!
       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Madeleine Sánchez, 2009.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick, petit
       implicit none


       integer, intent(in out) :: nr
       double precision, intent(in out) :: k1
       double precision, intent(in out) :: k2
       double precision, intent(in out) :: k3
       double precision, intent(in out) :: k4
       double precision, intent(in out) :: expt
       double precision, intent(in out) :: tact1
       double precision, intent(in out) :: tact2
       double precision, intent(in out) :: tact3
       double precision, intent(in out) :: tact4
       if (abs(expt) < petit) then !arr2 + ..
          if ((abs(tact2) < petit) .and. (abs(k2) < petit)) then !arr2 +..+
             if ((abs(tact3) < petit) .and. (abs(k3) < petit)) then !arr2 +..+..+
                if ((abs(tact4) < petit) .and. (abs(k4) > petit)) then !arr2+arr1 *
                   !        WRITE(nfick,*) nr, k1, tact1, k4
                   write(nfick, 10) nr, k1, tact1, k4
                end if
             end if
          else !arr2+arr2 +..
             if ((abs(tact3) < petit) .and. (abs(k3) < petit)) then !arr2+arr2+..+..
                if ((abs(tact4) < petit) .and. (abs(k4) < petit)) then !arr2+arr2  *
                   write(nfick, 11) nr, k1, tact1, k2, tact2
                else if ((abs(tact4) < petit) .and. (abs(k4) > petit)) then !arr2+arr2+arr1 *
                   write(nfick, 12) nr, k1, tact1, k2, tact2, k4
                end if
             else !arr2+arr2+arr2+..
                if ((abs(tact4) < petit) .and. (abs(k4) < petit)) then !arr2+arr2+arr2 *
                   write(nfick, 13) nr, k1, tact1, k2, tact2, k3, tact3
                else if ((abs(tact4) > petit)) then !arr2+arr2+arr2+arr2 *
                   write(nfick, 14) nr, k1, tact1, k2, tact2, k3, tact3, k4, tact4
                end if
             end if
          end if
       else !arr3+arr2+arr2+arr1 *
          if (tact1 > petit) then
             write(nfick, 15) nr, k1, expt, tact1, k2, tact2, k3, tact3, k4
          else if (tact1 < -petit) then
             write(nfick, 16) nr, k1, expt, (-tact1), k2, tact2, k3, tact3, k4
          else
             write(nfick, 17) nr, k1, expt, k2, tact2, k3, tact3, k4
          end if
       end if

       ! arr2+arr1
       10 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, '* DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16)

       ! arr2+arr2
       11 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk))')

       ! arr2+arr2+arr1
       12 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16)

       ! arr2+arr2+arr2
       13 format(6x, 'rk(ijk,', i3, ') =  ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk))')

       ! arr2+arr2+arr2+arr2
       14 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk))')

       ! arr3+arr2+arr2+arr1
       ! tact1 > 0
       15 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(', d23.16, ' * LOG(temp(ijk))', 3x, '&',  /, &
       8x, '- ', d23.16, '/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16)
       ! k1*(T^expt)*exp(-tact1/T)+k2*exp(-tact2/T)+k3*exp(-tact3/T)+k4*exp(-tact4/T)
       ! tact1 < 0
       16 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(', d23.16, ' *   LOG(temp(ijk))', 3x, '&',  /, &
       8x, '+ ', d23.16, '/temp(ijk)) + &',  /, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16)
       ! tact1 = 0
       17 format(6x, 'rk(ijk,', i3, ') = ', d23.16, 3x, '&',  /, &
       8x, ' * DEXP(', d23.16, ' * LOG(temp(ijk))) + &',  /, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16, 3x, '&',  /, &
       8x, ' * DEXP(-( ', d23.16, ' )/temp(ijk)) + & ', / 8x, &
       8x, d23.16)

    end subroutine wkc9

    subroutine wphot(nr, ntabphot, b, tabphot, it)

       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetics for tabulated photolysis at angles TABPHOT(I)
       !     for 1<=I<=NTABPHOT.
       !     The corresponding values for photolysis rates are BP(I,NR)

       !     at 90: no photolysis.
       !     at 0 and 90: first derivative = 0.

       !     Interpolation with standard cubic spline (second derivative =0
       !     at boundaries).

       !     Correction factor according to IT (IT =1 : O3 >2 OH).

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Bruno Sportisse, 2003.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick
       implicit none

       integer, intent(in out) :: nr
       integer, intent(in) :: ntabphot
       double precision, intent(in) :: b(ntabphotmax)
       double precision, intent(in) :: tabphot(ntabphotmax)
       integer, intent(in) :: it
       integer :: i, j
       double precision :: a(ntabphotmax), c(4, nintphotmax)

       !      print*,'nr= ',nr
       !      print*,'ntabphot= ',ntabphot
       !      print*,'b= ',b
       !      print*,'tabphot= ',tabphot
       !      print*,'it= ',it

       do i = 1, ntabphot
          a(i) = tabphot(i)
       end do

       call spl3(ntabphot, a, b, c)

       do i = 1, ntabphot - 1
          write(6, 100) b(i), (c(j, i), j = 1, 4), b(i + 1)
       end do

       write(nfick, 12)

       do i = 1, ntabphot - 1
          write(nfick, 11) a(i), a(i + 1)
          write(nfick, 14) nr, c(4, i)
          write(nfick, 15) nr, c(3, i), a(i), nr
          write(nfick, 15) nr, c(2, i), a(i), nr
          write(nfick, 15) nr, c(1, i), a(i), nr


       end do
       write(nfick, 16) nr, b(ntabphot)
       write(nfick, 13)
       write(nfick, 10) nr, nr


       if (it == 1) then
          write(nfick, 103) nr, nr
       end if

       !srf
       !10   FORMAT(6X,'if(att.lt.0.99999) rk(ijk,',i3,') = rk(ijk,',i3,') * att',/'!')
       10 format(6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * att(ijk)', /'!')
       11 format(6x, 'elseif(azi(ijk).ge.', d9.2, ' .and. azi(ijk).lt.', d9.2, ') then')
       14 format(7x, 'rk(ijk,', i3, ')=', d23.16)
       15 format(7x, 'rk(ijk,', i3, ')=', d23.16, '+(azi(ijk)-', d9.2, ') * rk(ijk,', i3, ')')
       12 format('!', / 6x, 'if(azi(ijk).lt.0.D0)then',  /, 7x, 'STOP')
       16 format(6x, 'elseif(azi(ijk).ge.90.D0)then',  /, 7x, 'rk(ijk,', i3, ')=', d23.16)
       13 format(6x, 'endif')
       100 format(6(2x, d23.16))

       103 format(6x, 'VO2  = 3.2D-11 * exp(70.D0 /temp(ijk)) * (0.21D0*SumM(ijk))',  / &
       6x, 'VN2  = 1.8D-11 * exp(110.D0/temp(ijk)) * (0.79D0*SumM(ijk))',  / &
       6x, 'VH2O = 2.2D-10 * YlH2O(ijk)' / &
       6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * VH2O / (VH2O + VN2 + VO2) ')

    end subroutine wphot

    subroutine spl3(n1, a, b, c)
       !------------------------------------------------------------------------
       !     Determination des coefficients des arcs de cubiques permettant
       !     d'interpoler b(a) en respectant la condition de derivee premiere
       !     nulle au bord

       !     n1 nombre de points, n0 nombre de segments

       !     Pour a(i) < aa < a(i+1)
       !     On ecrira
       !     bb = c(1,i) + c(2,i)*aa + c(3,i)*aa**2 + c(4,i)*aa**3
       !     ou aa = angle -a(i)

       !     c(1,i) = b(i)

       !     bp = c(2,i) + 2*c(3,i)*aa + 3*c(4,i)*aa**2
       !     bs = 2*c(3,i) + 6*c(4,i)*aa

       !     l'annulation de la derivee premiere au premier point fournit
       !     c(2,1) = 0

       !     1<i<n0
       !     la continuite de la fonction s'ecrit
       !     bb(i,a(i+1) = bb(i+1,a(i+1)) = b(i+1)
       !     c(1,i)+c(2,i)*d(i)+c(3,i)*d(i)**2+c(4,i)*d(i)**3 = b(i+1)
       !     ou d(i) = a(i+1)-a(i)

       !     la continuite des derivees premieres s'ecrit

       !     1<i<(n0-1)
       !     bp(i,a(i+1)) = bp(i+1,a(i+1))
       !     c(2,i) +2*c(3,i)*d(i)+3*c(4,i)*d(i)**2 = c2(i+1)

       !     la continuite des derivees secondes s'ecrit

       !     1<i<(n0-1)
       !     bs(i,a(i+1)) = bs(i+1,a(i+1))
       !     c(3,i)+3.*c(4,i)*d(i) = c(3,i+1)

       !     L'annulation de la derivee premiere au dernier point s'ecrit

       !     bp(n0,a(n1)) = 0
       !     c(2,n0)+2.*c(3,n0)*d(n0)+3.*c(4,n0)*d(n0)**2 = 0

       !     i/ La continuite de la fonction sert de relation de recurrence
       !     pour calculer c(4,i) en fonction de c(3,i)
       !     ii/ La continuite de la derivee premiere sert de relation de recurrence
       !     pour calculer c(2,i+1) en fonction de c(2,i),c(3,i),c(4,i)
       !     iii/ La continuite de la derivee seconde sert de relation de recurrence
       !     pour calculer c(3,i+1) en fonction de c(3,i),c(4,i)

       !     On utilise donc comme inconnue auxilliaire c(3,1) qui sera determine
       !     par la condition sur la derivee premiere au dernier point
       !------------------------------------------------------------------------
       implicit none

       integer, intent(in) :: n1
       double precision, intent(in) :: a(ntabphotmax)
       double precision, intent(in) :: b(ntabphotmax)
       double precision, intent(out) :: c(4, nintphotmax)
       integer :: i, n0

       double precision :: d(ntabphotmax)
       double precision :: c20(nintphotmax), c21(nintphotmax)
       double precision :: c30(nintphotmax)
       double precision :: c31(nintphotmax), c40(nintphotmax)
       double precision :: c41(nintphotmax), cc31

       n0 = n1 - 1
       do i = 1, n0
          c(1, i) = b(i)
          d(i) = a(i + 1) - a(i)
       end do
       c20(1) = 0.d0
       c21(1) = 0.d0
       c30(1) = 0.d0
       c31(1) = 1.d0
       c40(1) = (b(2) - c(1, 1) - c20(1) * d(1) - c30(1) * d(1)**2) / d(1)**3
       c41(1) = (-c21(1) * d(1) - c31(1) * d(1)**2) / d(1)**3
       do i = 2, n0
          !     c(2,i) +2*c(3,i)*d(i)+3*c(4,i)*d(i)**2 = c2(i+1)
          c20(i) = c20(i - 1) + 2.d0 * c30(i - 1) * d(i - 1) + 3.d0 * c40(i - 1) * d(i - 1)**2
          c21(i) = c21(i - 1) + 2.d0 * c31(i - 1) * d(i - 1) + 3.d0 * c41(i - 1) * d(i - 1)**2
          !     c(3,i)+3.*c(4,i)*d(i) = c(3,i+1)
          c30(i) = c30(i - 1) +3.d0 * c40(i - 1) * d(i - 1)
          c31(i) = c31(i - 1) +3.d0 * c41(i - 1) * d(i - 1)
          !     c(1,i)+c(2,i)*d(i)+c(3,i)*d(i)**2+c(4,i)*d(i)**3 = b(i+1)
          !     c4(i) = (b(i+1)-c(1,i)-c(2,i)*d(i)-c(3,i)*d(i)**2)/d(i)**3
          c40(i) = (b(i + 1) - c(1, i) - c20(i) * d(i) - c30(i) * d(i)**2) / d(i)**3
          c41(i) = (-c21(i) * d(i) - c31(i) * d(i)**2) / d(i)**3
       end do
       !     c(2,n0)+2.*c(3,n0)*d(n0)+3.*c(4,n0)*d(n0)**2 = 0
       !     (c20(n0)+c21(n0)*cc31) + 2.*(c30(n0)+c31(n0)*cc31)....
       cc31 = c20(n0) + 2.d0 * c30(n0) * d(n0) + 3.d0 * c40(n0) * d(n0)**2
       cc31 = -cc31 / (c21(n0) + 2.d0 * c31(n0) * d(n0) + 3.d0 * c41(n0) * d(n0)**2)

       do i = 1, n0
          c(2, i) = c20(i) + cc31 * c21(i)
          c(3, i) = c30(i) + cc31 * c31(i)
          c(4, i) = c40(i) + cc31 * c41(i)
       end do

    end subroutine spl3

    subroutine wphot_tuvonline(nr)
       use ModNficfort, only: nfick
       implicit none

       integer, intent(in) :: nr

       write(nfick, 10) nr, nr

       10 format(7x, 'rk(ijk,', i3, ') = Jphoto(ijk,', i3, ')')

    end subroutine wphot_tuvonline

    subroutine wtroe(nr, b1, b2, b3, b4, b5)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetic rates case for Troe reactions:

       !     k(T)=k0*M/(1+r)  * b5**(1/(1+[log10 r]**2))
       !     with r= k0*M/kinf, k0=b1*(T/300)**(-b2), kinf=b3*(T/300)**(-b4)
       !     b5=0.6 if not specified (case of KINETIC TROE4 ...).
       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion and Bruno Sportisse, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick
       implicit none

       integer, intent(in) :: nr
       double precision, intent(in) :: b1
       double precision, intent(in) :: b2
       double precision, intent(in) :: b3
       double precision, intent(in) :: b4
       double precision, intent(in) :: b5

       write(nfick, 10) b1, b2
       write(nfick, 11) b3, b4
       write(nfick, 12) nr, b5

       !     Modification BS/KS 06/2002
       !     Replace LOG10 by LOG and Effko*Rapk by Effko/Rapk

       10 format(6x, 'Effko(ijk) = ', d23.16, '* (temp(ijk) / 3.d2)', 3x, '&',  / &
       8x, '           **(- (', d23.16, '))')
       11 format(6x, 'Rapk(ijk) = ', d23.16, '* (temp(ijk) / 3.d2)', 3x, '&',  / &
       8x, '            **(- (', d23.16, '))')
       12 format(6x, 'rk(ijk,', i3, ') = (Effko(ijk) * SumM(ijk) / ', 3x, '&',  / &
       8x, '            ( 1.0d0 + Effko(ijk) * SumM(ijk) / Rapk(ijk))) *', 3x, '&',  / &
       8x, '            ', d10.4, '** (1.0d0 / (1.0d0 + ', 3x, '&',  / &
       8x, '             (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))')


    end subroutine wtroe

    subroutine wtroe7(nr, b1, b2, b3, b4, b5, b6, b7)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetics for reactions computed from equilibria.
       !     k(T)=k0*M/(1+r)  * b7**(1/(1+[log10 r]**2))
       !     with r= k0*M/kinf
       !     k0=b1*exp(-b2/T)*(T/300)**(-b3)
       !     kinf=b4*exp(-b5/T)*(T/300)**(-b6)

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion and Bruno Sportisse 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick
       implicit none

       integer, intent(in) :: nr
       double precision, intent(in) :: b1
       double precision, intent(in) :: b2
       double precision, intent(in) :: b3
       double precision, intent(in) :: b4
       double precision, intent(in) :: b5
       double precision, intent(in) :: b6
       double precision, intent(in) :: b7

       write(nfick, 10) b1, b3, b2
       write(nfick, 11) b4, b6, b5
       write(nfick, 12) nr, b7

       10 format(6x, 'Effko(ijk) = ', d23.16, '* (temp(ijk) / 3.d2)', 3x, '&',  / &
       8x, '           **(- (', d23.16, '))*', 3x, '&',  / &
       8x, '           dexp(-', d23.16, '/temp(ijk))')
       11 format(6x, 'Rapk(ijk) = ', d23.16, '* (temp(ijk) / 3.d2)', 3x, '&',  / &
       8x, '            **(- (', d23.16, '))*', 3x, '&',  / &
       8x, '           dexp(-', d23.16, '/temp(ijk))')
       12 format(6x, 'rk(ijk,', i3, ') = (Effko(ijk)*SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) / ', 3x, '&',  / &
       8x, '         Rapk(ijk))) * ', d10.4, '** (1.0d0 / (1.0d0 + ', 3x, '&',  / &
       8x, '          (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))')
       !     13    format(6x,'rk(',i3,') = facteur * (',D23.16,') * '/
       !     &       5x,'&      dexp((',D23.16,')/temp)')

    end subroutine wtroe7

    subroutine wtroe10(nr, b1, b2, b3, b4, b5, b6, b7, b8, b9, b10)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetic rates case for Troe reactions/MOCA.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion and Bruno Sportisse, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick, petit
       implicit none

       integer, intent(in) :: nr
       double precision, intent(in) :: b1
       double precision, intent(in) :: b2
       double precision, intent(in) :: b3
       double precision, intent(in) :: b4
       double precision, intent(in) :: b5
       double precision, intent(in) :: b6
       double precision, intent(in) :: b7
       double precision, intent(in) :: b8
       double precision, intent(in) :: b9
       double precision, intent(in) :: b10



       !     ! ! ! MOCA computes QLN = log(k)
       !     Modification BS/KS: 22/03/2002
       !     Error for b5,b6 -> b6,b5

       if (b6 > petit) then
          write(nfick, 10) b4, b6, b5
       else if (b6 < -petit) then
          write(nfick, 20) b4, (-b6), (b5)
       else
          write(nfick, 30) b4, b5
       end if

       !     End modification.

       !     Modification BS/KS: 20/03/2002
       !     Error b1 --> log(b1)

       if (b3 > petit) then
          write(nfick, 11) log(b1), b2, b3
       else if (b3 < -petit) then
          write(nfick, 21) log(b1), b2, (-b3)
       else
          write(nfick, 31) log(b1), b2
       end if

       !     End modification.
       write(nfick, 12)
       if (abs(b8) < petit) then
          write(nfick, 14) nr, b7
       else
          write(nfick, 15) (1.d0 - b7), b8, b7, b9, b10
          if (abs(b10) < petit) write(nfick, 16)
          write(nfick, 17) nr
       end if

       10 format(6x, 'Effko(ijk) = ', d23.16, '* SumM(ijk) * ', 3x, '&',  / &
       8x, ' DEXP(-(', d23.16, ')/temp(ijk) + ', 3x, '&',  / &
       8x, '(', d23.16, ') * LOG(temp(ijk)/3.d2))')

       20 format(6x, 'Effko(ijk) = ', d23.16, '*DEXP(', d23.16, '/temp(ijk) + ', 3x, '&',  / &
       8x, '(', d23.16, ') * LOG(temp(ijk)/3.d2)) * SumM(ijk)')

       30 format(6x, 'Effko(ijk) = ', d23.16, '*DEXP(',  / &
       5x, 's   ', d23.16, ' * LOG(temp(ijk)/3.d2)) * SumM(ijk)')

       11 format(6x, 'Rapk(ijk) = Effko(ijk) / DEXP(', d23.16, ' + ', 3x, '&',  / &
       8x, '     (', d23.16, ')* LOG(temp(ijk)/3.d2) - ', 3x, '&',  / &
       8x, '     (', d23.16, ')/temp(ijk))')

       21 format(6x, 'Rapk(ijk) = Effko(ijk) / DEXP(', d23.16, ' + ', 3x, '&',  / &
       8x, '     (', d23.16, ')* LOG(temp(ijk)/3.d2) + ', 3x, '&',  / &
       8x, '     (', d23.16, ')/temp(ijk))')

       31 format(6x, 'Rapk(ijk) = Effko(ijk) / DEXP(', d23.16, ' + ', 3x, '&',  /, &
       8x, '     (', d23.16, ')* LOG(temp(ijk)/3.d2) ) ')

       12 format(6x, 'facteur(ijk) = 1.d0/(1.d0+LOG10(Rapk(ijk))**2)')
       !     12    format(6x,'facteur = 0.d0',/
       !     &   6x,'if(Rapk.gt.0.d0) facteur = 1.d0/(1.d0+LOG10(Rapk)**2)')
       14 format(6x, 'rk(ijk,', i3, ') = DEXP( LOG(Effko(ijk)/(1.d0+Rapk(ijk))) + ', 3x, '&',  / &
       8x, '             LOG(', d23.16, ') * facteur(ijk))')

       15 format(6x, 'Fcent(ijk) = ', d23.16, ' * DEXP(-temp(ijk)/', 3x, '&',  / &
       8x, '(', d23.16, '))',  / 8x, '      + (', d23.16, ') * DEXP(-temp(ijk)/', 3x, '&', / &
       8x, '(', d23.16, '))', 3x, '&',  / 8x, '      +  DEXP (-(', d23.16, ')/temp(ijk))')
       16 format(6x, 'Fcent(ijk) = Fcent(ijk) -1.d0')
       17 format(6x, 'rk(ijk,', i3, ') = DEXP ( LOG(Effko(ijk)/(1.d0+Rapk(ijk))) ', 3x, '&',  / &
       8x, '           + facteur(ijk)*LOG(Fcent(ijk)))')

    end subroutine wtroe10

    subroutine wcv(nr, b1, b2, b3, b4, b5, b6, b7)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetics for variable stoichiometry for MOCA.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick, petit
       implicit none

       integer, intent(inout) :: nr
       double precision, intent(inout) :: b1
       double precision, intent(inout) :: b2
       double precision, intent(inout) :: b3
       double precision, intent(in) :: b4
       double precision, intent(in) :: b5
       double precision, intent(in) :: b6
       double precision, intent(in) :: b7

       double precision :: d54, d65, d76

       if (abs(b2) < petit) then
          call wk2(nr, b1, b3)
       else
          call wk3(nr, b1, b2, b3)
       end if
       !     write(nficK,10)nr,b4,b5,b6,b7
       d54 = (b5 - b4) *.05d0
       d65 = (b6 - b5) *.05d0
       d76 = (b7 - b6) *.05d0
       write(nfick, 12)
       if (abs(b4) > petit) then
          write(nfick, 13) nr, nr, b4
       else
          write(nfick, 18) nr
       end if
       write(nfick, 14) 260. , 280.
       if (abs(b4) > petit .or. abs(d54) > petit) then
          write(nfick, 15) nr, nr, b4, 260. , d54
       else
          write(nfick, 18) nr
       end if
       write(nfick, 14) 280. , 300.
       if (abs(b5) > petit .or. abs(d65) > petit) then
          write(nfick, 15) nr, nr, b5, 280. , d65
       else
          write(nfick, 18) nr
       end if
       write(nfick, 14) 300. , 320.
       if (abs(b6) > petit .or. abs(d76) > petit) then
          write(nfick, 15) nr, nr, b6, 300. , d76
       else
          write(nfick, 18) nr
       end if
       write(nfick, 16)
       if (abs(b7) > petit) then
          write(nfick, 13) nr, nr, b7
       else
          write(nfick, 18) nr
       end if
       write(nfick, 17)

       10 format('!', 6x, i3, 4(2x, d23.16))
       12 format(6x, 'if (temp(ijk).le.260.d0) then ')
       13 format(7x, ' rk(ijk,', i3, ') = rk(ijk,', i3, ') *  ', d23.16)
       14 format(6x, 'elseif(temp(ijk).gt.', d10.3, '.and.temp(ijk).le.', d10.3, ') then')
       15 format(7x, ' rk(ijk,', i3, ') = rk(ijk,', i3, ') * (', &
       d23.16, 3x, '&', / 8x, ' + (temp(ijk)-', d10.3, ') * (', d23.16, '))')
       16 format(6x, 'else')
       17 format(6x, 'endif', /'!')
       18 format(7x, 'rk(ijk,', i3, ') = 0.d0')

    end subroutine wcv

    subroutine wrcfe(nr, b1, b2, b3, b4, b5, b6)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetics for reactions computed from equilibria.
       !     k(T)=k0*M/(1+r)  * 0.6**(1/(1+[log10 r]**2)) * b5* exp(-b6/T)
       !     with r= k0*M/kinf, k0=b1*(T/300)**(-b2), kinf=b3*(T/300)**(-b4)

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion and Bruno Sportisse, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick, petit
       implicit none

       integer, intent(in out) :: nr
       double precision, intent(in) :: b1
       double precision, intent(in) :: b2
       double precision, intent(in) :: b3
       double precision, intent(in) :: b4
       double precision, intent(in out) :: b5
       double precision, intent(in out) :: b6

       call wk2(nr, b5, b6)

       write(nfick, 10) b1, b2
       write(nfick, 11) b3, b4
       write(nfick, 12)
       write(nfick, 13) nr, nr

       10 format(6x, 'Effko(ijk) = ', d23.16, '* (temp(ijk) / 3.d2)', 3x, '&',  / &
       8x, '           **(- (', d23.16, '))')
       11 format(6x, 'Rapk(ijk) = ', d23.16, '* (temp(ijk) / 3.d2)', 3x, '&',  / &
       8x, '            **(- (', d23.16, '))')
       12 format(6x, 'facteur(ijk) = (Effko(ijk) * SumM(ijk) / ( 1.0d0 + Effko(ijk) * SumM(ijk) / ', &
       3x, '&',  / 8x, '         Rapk(ijk))) * 0.6d0 ** (1.0d0 / (1.0d0 + ', 3x, '&',  / &
       8x, '          (LOG10(Effko(ijk) * SumM(ijk) / Rapk(ijk)))**2))')
       13 format(6x, 'rk(ijk,', i3, ') = facteur(ijk) * rk(ijk,', i3, ')')
       !     13    format(6x,'rk(',i3,') = facteur * (',D23.16,') * '/
       !     &       5x,'&      dexp((',D23.16,')/temp)')


    end subroutine wrcfe

    subroutine wspec (nr, ispebp)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write kinetic rates case for specific reactions.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion and Bruno Sportisse, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick
       implicit none


       integer, intent(in out) :: nr
       integer, intent(in out) :: ispebp

       ! RACM
       if (ispebp == -1) then
          write(nfick, 11) nr
       else if (ispebp == -2) then
          write(nfick, 12) nr
       else if (ispebp == -3) then
          write(nfick, 13) nr
       else if (ispebp == -4) then
          write(nfick, 14) nr
       else if (ispebp == -5) then
          write(nfick, 15) nr
       else if (ispebp == -6) then
          write(nfick, 16) nr
       else if (ispebp == -7) then
          write(nfick, 17) nr
          !CB2002
       else if (ispebp == -8) then
          write(nfick, 18) nr
       else if (ispebp == -9) then
          write(nfick, 19) nr
       else if (ispebp == -10) then
          write(nfick, 110) nr
       else if (ispebp == -11) then
          write(nfick, 111) nr
       else if (ispebp == -12) then
          write(nfick, 112) nr
       else
          write(*, *) 'ERROR: unknown specific reaction ', ispebp
          stop
       end if

       11 format(6x, 'rk(ijk,', i3, ') = SumM(ijk) * 6.0d-34 * (temp(ijk)/3.d2) ** (-2.3d0)')
       !     BS 05/02/2003 values given by RACM
       !     12   format(6x,'rk(',i3,') = 2.3d-13 * dexp(600.0d0 / temp)',/
       !     &     5x,'&            + 1.73d-33* SumM * dexp(1000.0d0 / temp)')
       !     13   format(6x,'rk(',i3,') = 3.22d-34 * dexp(2800.0d0 / temp) + ',/
       !     &     5x,'&            2.38d-54 * SumM * dexp(3200.0d0 / temp)')
       !     MODIF BS 06/06/2003 on the basis of CMAQ
       12 format(6x, 'rk(ijk,', i3, ') = 2.2d-13 * dexp(620.0d0 / temp(ijk))', 3x, '&',  / &
       8x, '            + 1.9d-33* SumM(ijk) * dexp(980.0d0 / temp(ijk))')
       13 format(6x, 'rk(ijk,', i3, ') = 3.08d-34 * dexp(2820.0d0 / temp(ijk)) + ', 3x, '&',  / &
       8x, '            2.66d-54 * SumM(ijk) * dexp(3180.0d0 / temp(ijk))')
       !     END MODIF
       14 format(6x, 'Effko(ijk) = 7.2d-15 * dexp(785.0d0 / temp(ijk))',  / &
       6x, 'Rapk(ijk) = 4.1d-16 * dexp(1440.0d0 / temp(ijk))',  / &
       6x, 'facteur(ijk) =1.9d-33 * dexp(725.0d0 / temp(ijk)) * SumM(ijk)',  / &
       6x, 'rk(ijk,', i3, ') = Effko(ijk) + facteur(ijk)/(1.0d0 + facteur(ijk) / Rapk(ijk))')
       15 format(6x, 'rk(ijk,', i3, ') = 1.5d-13 * (1.0d0 + 2.439d-20 * SumM(ijk))')
       16 format(6x, 'Rapk(ijk) = 3.4d-30 * (300./temp(ijk))**(3.2D0)*SumM(ijk)',  / &
       6x, 'Effko(ijk) = Rapk(ijk)/(4.77D-11*(300.D0/temp(ijk))**1.4D0)',  / &
       6x, 'rk(ijk,', i3, ')=(Rapk(ijk)/(1.D0+Effko(ijk)))*0.3D0**', 3x, '&',  / &
       8x, '(1.0d0/(1.0d0+ ((LOG10(Effko(ijk))-0.12D0)/1.2D0)**2))')
       17 format(6x, 'rk(ijk,', i3, ') = 2.0d-39 * YlH2O(ijk) * YlH2O(ijk)')
       !
       ! CB2002
       !     O + O2 -> O3 in NASA 2006
       18 format(6x, 'rk(ijk,', i3, ') = SumM(ijk) * 6.0d-34 * (temp(ijk)/3.d2) ** (-2.4d0)')

       !     N2O5 + H2O + H2O -> 2 HNO3 + H2O in CB2002
       19 format(6x, 'rk(ijk,', i3, ') = 1.80d-39 * YlH2O(ijk) * YlH2O(ijk)')
       !     HO2 + HO2 + H2O -> H2O2 + O2 + H2O in NASA 2006
       110 format(6x, 'rk(ijk,', i3, ') = 1.7d-33 * SumM(ijk) * ', 3x, '&',  / &
       8x, '            dexp(1000.0d0 / temp(ijk)) * (1 + 1.4D-21 * ', 3x, '&',  / &
       8x, '            YlH2O(ijk) * dexp(2200.0d0 / temp(ijk)))')
       !srf- 2 parenteses missing
       !    8X,'            YlH2O(ijk) * dexp(2200.0d0 / temp(ijk)')

       !     OH + HNO3 -> NO3 + H2O in NASA 2006
       111 format(6x, 'Effko(ijk) = 2.4d-14 * dexp(460.0d0 / temp(ijk))',  / &
       6x, 'Rapk(ijk) = 2.7d-17 * dexp(2199.0d0 / temp(ijk))',  / &
       6x, 'facteur(ijk) = 6.5d-34 * dexp(1335.0d0 / temp(ijk)) * SumM(ijk)',  / &
       6x, 'rk(ijk,', i3, ') = Effko(ijk) + facteur(ijk)/(1.0d0 + facteur(ijk) / Rapk(ijk))')
       !     CO + OH -> H + CO2 then H + O2 -> HO2  in NASA 2006
       112 format(6x, 'Effko(ijk) = 1.5d-13  * (temp(ijk)/300.D0)**(0.6D0)',  / &

       6x, 'Rapk(ijk) = 2.1d+09 * (temp(ijk)/300.D0)**(6.1D0)',  / &
       !srf- paremteses missing
       !    6X,'Rapk(ijk) = 2.1d+09 * (temp(ijk)/300.)**6.1)',/  &

       6x, 'rk(ijk,', i3, ')=(Effko(ijk)/(1.D0+Effko(ijk)/(Rapk(ijk)/SumM(ijk)))) * ', 3x, '&',  / &
       8x, '0.6D0**(1.0d0+(DLOG10(Effko(ijk)/(Rapk(ijk)/SumM(ijk)))**2.0D0)**(-1.0D0))')

    end subroutine wspec

    subroutine wtb(nr, ittb)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write third body kinetics.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion and Bruno Sportisse, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nfick
       implicit none


       integer, intent(in out) :: nr
       integer, intent(in out) :: ittb


       if (ittb == 1) then
          write(nfick, 11) nr, nr
       else if (ittb == 2) then
          write(nfick, 12) nr, nr
       else if (ittb == 3) then
          write(nfick, 13) nr, nr
       else if (ittb == 4) then
          write(nfick, 14) nr, nr
       else if (ittb == 5) then
          write(nfick, 15) nr, nr
       end if

       11 format(6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * SumM(ijk)')
       !     Seinfeld pp 22: N2 0.78; O2 0.21
       12 format(6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * SumM(ijk) * 0.2d0')
       13 format(6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * SumM(ijk) * 0.8d0')
       14 format(6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * YlH2O(ijk)')
       15 format(6x, 'rk(ijk,', i3, ') = rk(ijk,', i3, ') * SumM(ijk) * 5.8d-7')

    end subroutine wtb

    subroutine ww(nr, ne, ie1, ie2, ie3)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write reaction rates.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nficw
       implicit none


       integer, intent(in out) :: nr
       integer, intent(in out) :: ne
       integer, intent(in out) :: ie1
       integer, intent(in out) :: ie2
       integer, intent(in out) :: ie3


       if (ne == 1) then
          write(nficw, 10) nr, nr, ie1
          !     write(nficFF,100)nr,ie1
       else if (ne == 2) then
          write(nficw, 11) nr, nr, ie1, ie2
          !     write(nficFF,110)nr,ie1,ie2
       else if (ne == 3) then
          write(nficw, 12) nr, nr, ie1, ie2, ie3
          !     write(nficFF,120)nr,ie1,ie2,ie3
       end if

       10 format(6x, 'w(ijk,', i3, ') =  rk(ijk,', i3, ') * Y(ijk,', i3, ')')
       11 format(6x, 'w(ijk,', i3, ') =  rk(ijk,', i3, ') * Y(ijk,', i3, ')', ' * Y(ijk,', i3, ')')
       12 format(6x, 'w(ijk,', i3, ') =  rk(ijk,', i3, ') * Y(ijk,', i3, ')', ' * Y(ijk,', i3, ')', &
       ' * Y(ijk,', i3, ')')

    end subroutine ww
    !------------------------------------------------------------------------

    subroutine dw(nr, ne, ie1, ie2, ie3)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Write derivative of reaction rates.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------
       use ModNficfort, only: nficdw
       implicit none


       integer, intent(in out) :: nr
       integer, intent(in out) :: ne
       integer, intent(in out) :: ie1
       integer, intent(in out) :: ie2
       integer, intent(in out) :: ie3


       if (ne == 1) then
          write(nficdw, 10) nr, ie1, nr
          !     write(nficJJ,100)nr
       else if (ne == 2) then
          write(nficdw, 11) nr, ie1, nr, ie2
          !     write(nficJJ,110)nr,ie2
          write(nficdw, 11) nr, ie2, nr, ie1
          !     write(nficJJ,111)nr,ie1
       else if (ne == 3) then
          write(nficdw, 12) nr, ie1, nr, ie2, ie3
          !     write(nficJJ,120)nr,ie2,ie3
          write(nficdw, 12) nr, ie2, nr, ie1, ie3
          !     write(nficJJ,121)nr,ie1,ie3
          write(nficdw, 12) nr, ie3, nr, ie1, ie2
          !     write(nficJJ,122)nr,ie1,ie2
       end if

       10 format(6x, 'dw(ijk,', i3, ',', i3, ') =  rk(ijk,', i3, ')')
       11 format(6x, 'dw(ijk,', i3, ',', i3, ') =  rk(ijk,', i3, ') * Y(ijk,', i3, ')')
       12 format(6x, 'dw(ijk,', i3, ',', i3, ') =  rk(ijk,', i3, ') * Y(ijk,', i3, ')', ' * Y(ijk,', i3, ')')

    end subroutine dw

    subroutine wnonzero(s, nrtot, jer)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Create file for nonzero entries of Jacobian matrix:
       !     non_zero.dat

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Bruno Sportisse, 2003.

       !------------------------------------------------------------------------
       use ModNficfort, only: nficnz, petit
       use ModParametre
       implicit none

       double precision, intent(in out) :: s(nespmax, nrmax)
       integer, intent(in) :: nrtot
       integer, intent(in) :: jer(3, nrmax)
       double precision :: jnz(nespmax, nespmax)
       integer :: nr
       integer :: ie, je, nztot

       do ie = 1, nespmax
          do je = 1, nespmax
             jnz(ie, je) = 0
          end do
       end do

       do nr = 1, nrtot
          do ie = 1, nespmax
             if (abs(s(ie, nr)) > petit) then
                !     Products
                if (s(ie, nr) > 0.d0) then
                   do je = 1, 3
                      if (jer(je, nr) > 0) then
                         jnz(ie, jer(je, nr)) = 1
                      end if
                   end do
                else
                   !     Reactants
                   do je = 1, 3
                      if (jer(je, nr) > 0) then
                         jnz(ie, jer(je, nr)) = 1
                      end if
                   end do
                end if
             end if
          end do
       end do

       nztot = 0
       do ie = 1, nespmax
          do je = 1, nespmax
             if (jnz(ie, je) == 1) nztot = nztot + 1
          end do
       end do

       write(nficnz, *) nztot
       do ie = 1, nespmax
          do je = 1, nespmax
             if (jnz(ie, je) == 1) write(nficnz, *) ie, je
          end do
       end do

    end subroutine wnonzero

    subroutine wfj(s, nr, jer)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Update chemical production term and Jacobian Matrix.

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Pierre Plion, 2002.

       !------------------------------------------------------------------------

       use ModNficfort, only: nficf, nficj, petit
       use ModParametre
       implicit none

       double precision, intent(in out) :: s(nespmax, nrmax)
       integer, intent(in out) :: nr
       integer, intent(in out) :: jer(3, nrmax)


       integer :: ie, je

       do ie = 1, nespmax
          if (abs(s(ie, nr)) > petit) then
             !     Products
             if (s(ie, nr) > 0.d0) then
                !     Products with  stoichiometry = 1
                if (abs(s(ie, nr) - 1.d0) < petit) then
                   write(nficf, 30) ie, ie, nr
                   do je = 1, 3
                      if (jer(je, nr) > 0) then
                         write(nficj, 31) ie, jer(je, nr), ie, jer(je, nr), nr, jer(je, nr)
                      end if
                   end do
                else
                   write(nficf, 10) ie, ie, s(ie, nr), nr
                   do je = 1, 3
                      if (jer(je, nr) > 0) write(nficj, 11) ie, jer(je, nr), ie, jer(je, nr), &
                      s(ie, nr), nr, jer(je, nr)
                   end do
                end if
             else
                !     Reactants
                if (abs(s(ie, nr) + 1.d0) < petit) then
                   !     Reactants with stochiometry = 1
                   write(nficf, 40) ie, ie, nr
                   do je = 1, 3
                      if (jer(je, nr) > 0) write(nficj, 41) &
                      ie, jer(je, nr), ie, jer(je, nr), nr, jer(je, nr)
                   end do
                else
                   write(nficf, 20) ie, ie, (-s(ie, nr)), nr
                   do je = 1, 3
                      if (jer(je, nr) > 0) write(nficj, 21) ie, jer(je, nr), ie, jer(je, nr), &
                      (-s(ie, nr)), nr, jer(je, nr)
                   end do
                end if
             end if
          end if
       end do

       10 format(6x, 'chem(ijk,', i3, ') = chem(ijk,', i3, ') + ', d23.16, ' * w(ijk,', i3, ')')
       11 format(6x, 'JacC(ijk,', i3, ',', i3, ') = JacC(ijk,', i3, ',', i3, ')+', &
       d23.16, '*dw(ijk,', i3, ',', i3, ')')

       20 format(6x, 'chem(ijk,', i3, ') = chem(ijk,', i3, ') - ', d23.16, ' * w(ijk,', i3, ')')
       21 format(6x, 'JacC(ijk,', i3, ',', i3, ') = JacC(ijk,', i3, ',', i3, ')-', &
       d23.16, '*dw(ijk,', i3, ',', i3, ')')

       30 format(6x, 'chem(ijk,', i3, ') = chem(ijk,', i3, ') + w(ijk,', i3, ')')
       40 format(6x, 'chem(ijk,', i3, ') = chem(ijk,', i3, ') - w(ijk,', i3, ')')
       31 format(6x, 'JacC(ijk,', i3, ',', i3, ') = JacC(ijk,', i3, ',', i3, &
       ') + dw(ijk,', i3, ',', i3, ')')
       41 format(6x, 'JacC(ijk,', i3, ',', i3, ') = JacC(ijk,', i3, ',', i3, &
       ') - dw(ijk,', i3, ',', i3, ')')

       !~ 10   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') + ',d23.16,' * w(ijk,',i3,')')
       !~ 11   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,')+',  &
       !~ d23.16,'*dw(ijk,',i3,',',i3,')')

       !~ 20   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') - ',d23.16,' * w(ijk,',i3,')')
       !~ 21   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,')-',  &
       !~ d23.16,'*dw(ijk,',i3,',',i3,')')

       !~ 30   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') + w(ijk,',i3,')')
       !~ 40   FORMAT(6X,'chem(ijk,',A,') = chem(ijk,',A,') - w(ijk,',i3,')')
       !~ 31   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,  &
       !~ ') + dw(ijk,',i3,',',i3,')')
       !~ 41   FORMAT(6X,'JacC(',A,',',i3,') = JacC(',A,',',i3,  &
       !~ ') - dw(ijk,',i3,',',i3,')')


    end subroutine wfj

    subroutine wpl(s, nr, jer)
       !------------------------------------------------------------------------

       !     -- DESCRIPTION

       !     Update production and loss terms (P-Lc formulation).

       !------------------------------------------------------------------------

       !     -- AUTHOR(S)

       !     Bruno Sportisse, 2003.

       !------------------------------------------------------------------------
       use ModNficfort, only: nficprod, nficloss, petit
       use ModParametre
       use ModAuxnom, only: nom_aux
       implicit none

       double precision, intent(in out) :: s(nespmax, nrmax)
       integer, intent(in out) :: nr
       integer, intent(in out) :: jer(3, nrmax)

       integer :: ie

       do ie = 1, nespmax
          if (abs(s(ie, nr)) > petit) then
             !     Products
             if (s(ie, nr) > 0.d0) then
                !     Products with  stoichiometry = 1
                if (abs(s(ie, nr) - 1.d0) < petit) then
                   write(nficprod, 30) ie, ie, nr
                   !        WRITE(nficprod,30)nom_aux(ie),nom_aux(ie),nr
                else
                   write(nficprod, 10) ie, ie, s(ie, nr), nr
                   !        WRITE(nficprod,10)nom_aux(ie),nom_aux(ie),s(ie,nr),nr
                end if
             else
                !     Reactants
                if (abs(s(ie, nr) + 1.d0) < petit) then
                   !     Reactants with stochiometry = 1
                   write(nficloss, 40) ie, ie, nr, ie
                   !        WRITE(nficloss,40)nom_aux(ie),nom_aux(ie),nr,nom_aux(ie)
                else
                   write(nficloss, 20) ie, ie, (-s(ie, nr)), nr, ie
                   !        WRITE(nficloss,20)nom_aux(ie),nom_aux(ie),(-s(ie,nr)),nr,nom_aux(ie)
                end if
                !     Correction if bimolecular reactant: none
                !     Because d(c*c)/dc=c and not 2*c as computed !
                !     nm=0
                !     do je=1,3
                !     if (jer(je,nr).eq.ie) nm=nm+1
                !     enddo
                !     if (nm.gt.1) write(nficloss,50)ie,ie,nm*1.
             end if
          end if
       end do

       10 format(6x, 'prod(ijk,', i3, ') = prod(ijk,', i3, ') + ', d23.16, ' * w(ijk,', i3, ')')
       30 format(6x, 'prod(ijk,', i3, ') = prod(ijk,', i3, ') + w(ijk,', i3, ')')

       20 format(6x, 'loss(ijk,', i3, ') = loss(ijk,', i3, ') + ', d23.16, ' * dw(ijk,', i3, ',', &
       i3, ')')
       40 format(6x, 'loss(ijk,', i3, ') = loss(ijk,', i3, ') + dw(ijk,', i3, ',', i3, ')')

       50 format(6x, 'loss(ijk,', i3, ') = loss(ijk,', i3, ')/', d23.16)

       !~ 10   FORMAT(6X,'prod(ijk,',A,') = prod(ijk,',A,') + ',d23.16,' * w(ijk,',i3,')')
       !~ 30   FORMAT(6X,'prod(ijk,',A,') = prod(ijk,',A,') + w(ijk,',i3,')')

       !~ 20   FORMAT(6X,'loss(ijk,',A,') = loss(ijk,',A,') + ',d23.16,' * dw(ijk,',i3,',',  &
       !~ A,')')
       !~ 40   FORMAT(6X,'loss(ijk,',A,') = loss(ijk,',A,') + dw(ijk,',i3,',',A,')')

       !~ 50   FORMAT(6X,'loss(ijk,',i3,') = loss(ijk,',i3,')/',d23.16)

    end subroutine wpl

end module ModGenerator