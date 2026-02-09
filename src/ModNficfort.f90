module ModNficfort
    double precision, parameter :: petit = 1.d-30
    
    ! Labels for files kinetic.f, fexchem.f and jacdchemdc.f
    ! Label for non_zero.dat
    ! Labels for fexprod.f, fexloss.f,rates.f,dratedc.f
    
    integer :: nfick, nficf, nficj, nficnz, nficloss, nficprod
    integer :: nficw, nficdw
    !common /labfiles/ nfick, nficf, nficj, nficnz, nficloss, nficprod, nficw, nficdw
    
    !   Files for species and mechanism
    character (len = 20) :: filemeca, filespecies
    !common /namefiles/ filemeca, filespecies
    
    ! Current free label
    integer :: ipiste
    !common /piste/ ipiste
End module ModNficfort

