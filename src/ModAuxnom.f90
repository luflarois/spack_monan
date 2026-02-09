module ModAuxnom
    use ModParametre, only: nespmax
    implicit none

    character(len=12),dimension(nespmax)  :: nom_aux
    character (len=240) :: chemical_mechanism
    integer :: nblanc

    !common /nomesp/ nom_aux, chemical_mechanism
    !common /nblanc/ nblanc

end module ModAuxnom