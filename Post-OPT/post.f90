program post
   !Post-processor for spack
   !Author: Luiz Flavio Rodrigues
   !        23/10/2007
   !Read some files from Spack and then serializes the arrays
   !
   character(len = *), parameter :: file1 = "chem_spack_fexchem.f90"
   character(len = *), parameter :: file1_s = "chem_spack_fexchem_s.f90"
   character(len = *), parameter :: tipo1 = "chem"
   character(len = *), parameter :: file2 = "chem_spack_fexloss.f90"
   character(len = *), parameter :: file2_s = "chem_spack_fexloss_s.f90"
   character(len = *), parameter :: tipo2 = "loss"
   character(len = *), parameter :: file3 = "chem_spack_fexprod.f90"
   character(len = *), parameter :: file3_s = "chem_spack_fexprod_s.f90"
   character(len = *), parameter :: tipo3 = "prod"
   character(len = *), parameter :: file4 = "chem_spack_jacdchemdc.f90"
   character(len = *), parameter :: file4_s = "chem_spack_jacdchemdc_s.f90"

   call post_jacc(file4, file4_s)
   call post_convert(file1, file1_s, tipo1)
   call post_convert(file2, file2_s, tipo2)
   call post_convert(file3, file3_s, tipo3)

end program post

subroutine post_jacc(file1, file1_s)
   use chem1_list, only: nspecies, nr
   implicit none
   character(len = *), intent(in) :: file1
   character(len = *), intent(in) :: file1_s

   character(len = 45) :: aux
   character(len = 150) :: linha
   character(len = 150), dimension(10000) :: buffer
   logical, dimension(10000) :: item_add
   character(len = 17) :: item
   character(len = 15) :: lixo1
   character(len = 18) :: lixo2
   character(len = 3) :: num
   integer, dimension(nspecies) :: fexprod_count
   integer :: istat, ngas, cl, buf_count, i, j
   logical :: laco = .false.

   fexprod_count = 0
   cl = 0

   open(unit = 22, file = file1)
   open(unit = 55, file = file1_s)

   do

      if (.not. laco) then
         !PRINT *,'Reading jacc.....',cl
         cl = cl + 1
         read(22, fmt = '(A150)', iostat = istat) linha
         if (istat /= 0) exit
         if (linha(1:17) == '        JacC=0.d0') cycle
         write (55, fmt = '(A150)') linha
      end if
      if (linha(7:26) == 'DO ijk=ijkbeg,ijkend' .and. .not. laco) then
         print *, 'Laco encontrado Jacc.....'
         laco = .true.
         buf_count = 0
         cycle
      end if

      if (laco) then
         !PRINT *,'Reading Jacc.....',cl
         cl = cl + 1
         buf_count = buf_count + 1
         read(22, fmt = '(A150)') buffer(buf_count)
         if (buffer(buf_count) (8:13) == 'END DO' .and. laco) then
            item_add = .false.
            !=====================================
            do i = 1, buf_count - 1
               if (item_add(i)) cycle !Já computado
               item = buffer(i) (7:23)
               item_add(i) = .true.
               write (55, fmt = '(A,$)') buffer(i) (1:26) // trim(buffer(i) (44:))
               do j = i + 1, buf_count - 1
                  if (item == buffer(j) (7:23) .and. .not. item_add(j)) then
                     item_add(j) = .true.
                     write (55, fmt = '(A)') ' &'
                     write (55, fmt = '(A, $)') repeat(' ', 26) // trim(buffer(j) (44:))
                  end if
               end do
               write(55, fmt = '(A)') ''
            end do
            !======================================

            write (55, fmt = '(A150)') buffer(buf_count)
            laco = .false.
            cycle
         end if
      end if

   end do

   close(unit = 22)
   close(unit = 55)

end subroutine post_jacc

subroutine post_convert(file1, file1_s, tipo)
   use chem1_list, only: nspecies, nr
   implicit none
   character(len = *), intent(in) :: file1
   character(len = *), intent(in) :: file1_s
   character(len = *), intent(in) :: tipo
   character(len = 45), dimension(nspecies, nr) :: fexprod
   character(len = 45) :: aux
   character(len = 80) :: linha
   character(len = 15) :: lixo1
   character(len = 18) :: lixo2
   character(len = 3) :: num
   integer, dimension(nspecies) :: fexprod_count
   integer :: istat, ngas, cl
   logical :: laco

   fexprod_count = 0
   cl = 0
   fexprod(:, :) = 'empty'
   laco = .false.


   open(unit = 22, file = file1)
   open(unit = 55, file = file1_s)
   print *, 'Arquivo:' // trim(file1)
   do

      if (.not. laco) then
         print *, 'Reading (Fora do laco).....', cl
         cl = cl + 1
         read(22, fmt = '(A80)', iostat = istat) linha
         if (istat /= 0) exit
         if (linha(1:17) == '       ' // tipo // '=0.d0') cycle
         write (55, fmt = '(A80)') linha
      end if

      if (linha(1:29) == '      DO ijk=ijkbegin,ijkend' .and. .not. laco) then
         print *, 'Laco encontrado.....'
         laco = .true.
         cycle
      end if

      if (linha(1:12) == '      END DO' .and. laco) then
         print *, 'Encerrando compactacao....'
         call fecha_laco(nspecies, nr, fexprod, fexprod_count)
         write (55, fmt = '(A80)') linha
         laco = .false.
         cycle
      end if

      if (laco) then
         print *, 'Reading.....', cl
         cl = cl + 1
         read(22, fmt = '(A15,I3,A18,A45)') lixo1, ngas, lixo2, aux
         if (lixo1(1:12) == '      END DO') then
            print *, 'Encerrando compactacao....'
            call fecha_laco(nspecies, nr, fexprod, fexprod_count, tipo)
            write (55, fmt = '(A15)') lixo1
            laco = .false.
            cycle
         end if
         fexprod_count(ngas) = fexprod_count(ngas) + 1
         fexprod(ngas, fexprod_count(ngas)) = aux
      end if

   end do

   close(unit = 22)
   close(unit = 55)

end subroutine post_convert

subroutine fecha_laco(nspecies, nr, fexprod, fexprod_count, tipo)
   implicit none
   integer, intent(in) :: nspecies, nr
   character(len = *), intent(in) :: tipo
   character(len = 45), dimension(nspecies, nr), intent(in) :: fexprod
   integer, dimension(nspecies), intent(in) :: fexprod_count
   integer :: i, j


   do i = 1, nspecies
      if (trim(fexprod(i, 1)) == 'empty') cycle
      write(55, fmt = '(A15,I3,A4,A,$)') '      ' // tipo // '(ijk,', i, ') = ', trim(fexprod(i, 1))
      !print*,' xx',i,trim(fexprod(i,1))
      if (fexprod_count(i) > 1) then
         write(55, fmt = '(A)') ' &'
      else
         write(55, fmt = '(A)') ''
         cycle
      end if
      do j = 2, fexprod_count(i) - 1
         write(55, fmt = '(A)') '          ' // trim(fexprod(i, j)) // ' &'
      end do
      write(55, fmt = '(A)') '          ' // trim(fexprod(i, fexprod_count(i)))
   end do

end subroutine fecha_laco



