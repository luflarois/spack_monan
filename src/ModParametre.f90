module ModParametre
  ! Module defining global parameters for the simulation
  ! Author: rodrigues, l.f. 
  ! email: luflarois@gmail.com
  ! Date: 2026Jan05
  implicit none

    integer, parameter :: nrmax = 350
    integer, parameter :: nphotmax = 2
    integer, parameter :: ntabphotmax = 12
    integer, parameter :: nintphotmax = ntabphotmax - 1
    integer, parameter :: nbpmax = 12 + ntabphotmax
    integer, parameter :: nespmax = 135
    integer, parameter :: nphase = 3
    integer, parameter :: nequilx = 25
    integer, parameter :: nionx = 25
    integer, parameter :: nangl = 20
    integer, parameter :: nlo = 3 * nangl
    integer, parameter :: ntab1 = 11
    integer, parameter :: neqmax = nionx + nequilx
    ! Parâmetros dimensionais
    integer, parameter :: nequilx_plus_nionx = nequilx + nionx

end module ModParametre