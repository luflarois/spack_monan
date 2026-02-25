module ModMonan
    use ModFiccom, only: nom, nesp
    use ModInit, only: csou, cdry, cwet, coff, cfdd, ctra
    implicit none
    private
    public :: create_monan_registry

contains

    subroutine create_monan_registry()
        implicit none
        integer :: ns
    
            !------------------------------------------------------------------------        

        ! This subroutine would contain the logic to create the registry of MONAN includes.
        ! The actual implementation would depend on the specific requirements and structure of the MONAN includes.
        print *, 'Creating MONAN registry...'
        open (unit=10, file='chem_registry.xml', status='replace')

        write(10,*) "#ifdef DO_CHEMISTRY"
        write(10,*)     '        <var_array name="chem" type="real" dimensions="nVertLevels nCells Time">'
        ! Laço nas variáveis químicas para escrever as linhas correspondentes no arquivo XML. 
        do ns = 1, nesp(2)
            print *, "LFR-MONAN: ns, nom(ns)=", ns, nom(ns)
            write(10,*) '                <var_name="'//trim(nom(ns))//'" array_group="chem_con" units="kg m^{-3}"'
            write(10,*) '                         description="species concentration for specie '//trim(nom(ns))//'" />'           
        end do
        write(10,*) " "
        write(10,*)      '        </var_array>'
        write(10,*) "#endif"

        close(10)
        print *, 'MONAN registry created successfully.'

    end subroutine create_monan_registry

end module ModMonan