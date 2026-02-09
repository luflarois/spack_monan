module ModGestion
    implicit none

contains

    subroutine part(chdon, mot, imot, nmot, nblanc)
       implicit none
       integer, parameter :: nbmot = 100

       character (len = 500), intent(in) :: chdon
       character (len = 500), intent(out) :: mot(nbmot)
       integer, intent(out) :: imot(nbmot)
       integer, intent(out) :: nmot
       integer, intent(in) :: nblanc
       !     Decomposition of a sequence (CHDON) in NMOT words (MOT)

       integer :: i0, iblanc, ilettre, nl

       !common /nblanc/ nblanc

       do i0 = 1, nblanc
          mot(i0) = '                                                    '
       end do
       nmot = 0
       iblanc = 1
       ilettre = 0
       do nl = 1, 500
          if (chdon(nl:nl) /= ' ') then
             if (iblanc == 1) then
                iblanc = 0
                nmot = nmot + 1
                ilettre = 1
                mot(nmot) (ilettre:ilettre) = chdon(nl:nl)
             else
                ilettre = ilettre + 1
                mot(nmot) (ilettre:ilettre) = chdon(nl:nl)
             end if
          else
             if (iblanc == 0) imot(nmot) = ilettre
             iblanc = 1
          end if
       end do

    end subroutine part


   subroutine reel(r, ch, ich)
      implicit none

      double precision, intent(in out) :: r
      character (len = 500), intent(in out) :: ch
      integer, intent(in out) :: ich

      !     Read a real in a string
      character (len = 2) :: chtemp
      character (len = 20) :: chtemp2
      write(chtemp, '(i2)') ich

      chtemp2 = '(e' // chtemp // '.0)'
      read(ch, chtemp2) r

   end subroutine reel

   subroutine entier(i, ch, ich)
      implicit none

      integer, intent(in out) :: i
      character (len = 500), intent(in out) :: ch
      integer, intent(in out) :: ich

      !     Read an integer in a string


      character (len = 1) :: chtemp
      character (len = 20) :: chtemp2
      write(chtemp, '(i1)') ich

      chtemp2 = '(i' // chtemp // ')'
      read(ch, chtemp2) i

   end subroutine entier

end module ModGestion