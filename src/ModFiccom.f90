module ModFiccom
  use ModParametre
  implicit none
  
  ! Variáveis do common /cinetique1/
  double precision, dimension(nbpmax, nrmax) :: bp
  double precision, dimension(5) :: debug
  double precision, dimension(nespmax, nrmax) :: s, slump, seff
  integer :: nr
  integer, dimension(nphase) :: nrp
  integer, dimension(nrmax) :: nb, molec, ispebp, indaqr, ittb, iprecalc, jep
  integer, dimension(3, nrmax) :: jer
  integer, dimension(nrmax, nphase) :: imolec1, imolec2, imolec3, iarr1, iarr2, iarr3, iarr4, iarr5, iarr6, iarr7, iarr8
  integer, dimension(nphase) :: nrmol1, nrmol2, nrmol3, narr1, narr2, narr3, narr4, narr5, narr6, narr7, narr8, naq2
  integer, dimension(nrmax, nphase) :: inaq2, indthird
  integer, dimension(nphase) :: nthird
  integer :: iunitaq, iunitgas
  integer, dimension(nespmax, nphase) :: iemonodi, iedimono
  integer, dimension(nrmax, nphase) :: irmonodi
  integer :: iphasecom
  
  ! Variáveis do common /photolyse/
  double precision, dimension(nphotmax, nangl) :: xtab, ytab
  double precision, dimension(nphotmax, nlo) :: cpolg
  double precision, dimension(ntabphotmax) :: tabphot
  integer, dimension(nrmax) :: iphot, iphotinv
  integer, dimension(nphotmax) :: ntab
  integer :: ntabphot, ireversetab, nrphot
  
  ! Variáveis do common /especes/
  character(len=12), dimension(nespmax) :: nom
  integer, dimension(nespmax) :: inom
  
  ! Variáveis do common /dimens/
  integer, dimension(nphase) :: nesp, ndiff
  integer :: nalg
  
  ! Variáveis do common /reducphys/
  double precision, dimension(nespmax, 5) :: xlphy
  integer, dimension(nespmax, 5) :: indpur
  integer, dimension(nespmax) :: idlump, indlump, idifford
  
  ! Variáveis do common /indicateurs/
  integer :: indicaqcom, ireductcom
  
  ! Variáveis do common /indic1/
  integer :: ixl
  
  ! Variáveis do common /aqueous/
  double precision, dimension(nespmax) :: alpha, dg, rmol
  double precision :: theta, diam, xl0, xlmax, t0, t1cloud
  integer :: indmod, iswich
  integer, dimension(nespmax) :: indaq, ihreac
  integer, dimension(nrmax) :: iheq, ieqhion
  integer :: jhplus, johmoin
  
  ! Variáveis do common /comph/
  double precision :: ph
  integer :: iph, nion, nitph
  integer, dimension(nespmax) :: iion, ival
  
  ! Variáveis do common /const/
  double precision :: av
  
  ! Variáveis do common /equil/
  double precision, dimension(nionx) :: xk1, xk2
  double precision, dimension(nespmax, nionx) :: seqion
  integer, dimension(nequilx, 10) :: iesp
  integer, dimension(nrmax) :: nneq
  integer, dimension(nequilx) :: jion1, jion2, nr1, nr2
  integer, dimension(nequilx_plus_nionx) :: jaq, jgaz
  integer :: nequil, nequil11, nequil12, nequil13, nequil21, nequil22, nequil23, nequil3
  integer :: nequil41, nequil42, nequil51, nequil52
  integer, dimension(neqmax) :: iequil11, iequil12, iequil13, iequil21, iequil22, iequil23, iequil3
  integer, dimension(neqmax, 2) :: iequil41, iequil42
  integer, dimension(neqmax, 5) :: iequil51, iequil52
  integer, dimension(nionx) :: jhpoh
  

contains
  
end module ModFiccom