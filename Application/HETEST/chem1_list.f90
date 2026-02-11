MODULE chem1_list
  IMPLICIT NONE
  
  
  CHARACTER(LEN=24),PARAMETER :: chemical_mechanism='HETEST'
  INTEGER,PARAMETER :: maxnspecies= 200
  INTEGER,PARAMETER :: nspecies=031
  
  
  !Name of species 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspecies) :: spc_name=(/ &
   'O3  ' & !
     ,'H2  ' & !
     ,'H2O2' & !
     ,'NO  ' & !
     ,'NO2 ' & !
     ,'NO3 ' & !
     ,'N2O5' & !
     ,'HONO' & !
     ,'HNO3' & !
     ,'HNO4' & !
     ,'CO  ' & !
     ,'O3P ' & !
     ,'O1D ' & !
     ,'HO  ' & !
     ,'HO2 ' & !
     ,'CH4 ' & !
     ,'HCHO' & !
     ,'OP1 ' & !
     ,'MO2 ' & !
     ,'CL  ' & !
     ,'HCL ' & !
     ,'CLO ' & !
     ,'CLN3' & !
     ,'HOCL' & !
     ,'CLO2' & !
     ,'CL2 ' & !
     ,'BRO ' & !
     ,'BR  ' & !
     ,'HBR ' & !
     ,'BRN3' & !
     ,'HOBR' & !
   /)
  
  
  !Number of each specie   
  INTEGER,PARAMETER :: O3  =001
  INTEGER,PARAMETER :: H2  =002
  INTEGER,PARAMETER :: H2O2=003
  INTEGER,PARAMETER :: NO  =004
  INTEGER,PARAMETER :: NO2 =005
  INTEGER,PARAMETER :: NO3 =006
  INTEGER,PARAMETER :: N2O5=007
  INTEGER,PARAMETER :: HONO=008
  INTEGER,PARAMETER :: HNO3=009
  INTEGER,PARAMETER :: HNO4=010
  INTEGER,PARAMETER :: CO  =011
  INTEGER,PARAMETER :: O3P =012
  INTEGER,PARAMETER :: O1D =013
  INTEGER,PARAMETER :: HO  =014
  INTEGER,PARAMETER :: HO2 =015
  INTEGER,PARAMETER :: CH4 =016
  INTEGER,PARAMETER :: HCHO=017
  INTEGER,PARAMETER :: OP1 =018
  INTEGER,PARAMETER :: MO2 =019
  INTEGER,PARAMETER :: CL  =020
  INTEGER,PARAMETER :: HCL =021
  INTEGER,PARAMETER :: CLO =022
  INTEGER,PARAMETER :: CLN3=023
  INTEGER,PARAMETER :: HOCL=024
  INTEGER,PARAMETER :: CLO2=025
  INTEGER,PARAMETER :: CL2 =026
  INTEGER,PARAMETER :: BRO =027
  INTEGER,PARAMETER :: BR  =028
  INTEGER,PARAMETER :: HBR =029
  INTEGER,PARAMETER :: BRN3=030
  INTEGER,PARAMETER :: HOBR=031
  
  
  !for memory allocattion: 
  !This parameters are use for documentation only. 
  !Use them in a program in substitution of numerical terms.
  INTEGER,PARAMETER :: src = 1 ! source term 
  INTEGER,PARAMETER :: ddp = 2 ! dry deposition 
  INTEGER,PARAMETER :: wdp = 3 ! wet deposition 
  INTEGER,PARAMETER :: fdda = 4! four-dim assimilation 
  INTEGER,PARAMETER :: offline = 5! ! off-line emissions: 
                                  !=1, emission will be read from file
                                  !=0, emission will be calculated during the model simulation (on-line emission)
  INTEGER,PARAMETER :: transport = 6! ! off-line emissions: 
                                  !=1, yes
                                  !=0, no transport
  INTEGER,PARAMETER :: on = 1
  INTEGER,PARAMETER :: off = 0
  
  
  ! spaction(specie,[1=source,2=drydep,3=wetdep,4=fdda,5=offline emission,6=transport]) ]) 
  INTEGER,PARAMETER,DIMENSION(6,nspecies) :: spc_alloc=RESHAPE((/ &
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! O3 - 001
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! H2 - 002
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! H2O2 - 003
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! NO - 004
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! NO2 - 005
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! NO3 - 006
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! N2O5 - 007
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! HONO - 008
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! HNO3 - 009
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! HNO4 - 010
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! CO - 011
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! O3P - 012
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! O1D - 013
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! HO - 014
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! HO2 - 015
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! CH4 - 016
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! HCHO - 017
    0 , 1 , 1 , 0 , 0 , 1 ,   & ! OP1 - 018
    0 , 1 , 0 , 0 , 0 , 1 ,   & ! MO2 - 019
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! CL - 020
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! HCL - 021
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! CLO - 022
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! CLN3 - 023
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! HOCL - 024
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! CLO2 - 025
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! CL2 - 026
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! BRO - 027
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! BR - 028
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! HBR - 029
    0 , 0 , 0 , 0 , 0 , 1 ,   & ! BRN3 - 030
    0 , 0 , 0 , 0 , 0 , 1    & ! HOBR - 031
    /),(/6,nspecies/))
  
  
  INTEGER,PARAMETER,DIMENSION(nspecies) :: spc_uveq=(/ &
   02, & !O3
   00, & !H2
   07, & !H2O2
   04, & !NO
   09, & !NO2
   10, & !NO3
   11, & !N2O5
   00, & !HONO
   13, & !HNO3
   14, & !HNO4
   00, & !CO
   00, & !O3P
   00, & !O1D
   00, & !HO
   00, & !HO2
   00, & !CH4
   00, & !HCHO
   00, & !OP1
   00, & !MO2
   00, & !CL
   00, & !HCL
   00, & !CLO
   00, & !CLN3
   00, & !HOCL
   00, & !CLO2
   00, & !CL2
   00, & !BRO
   00, & !BR
   00, & !HBR
   00, & !BRN3
   00  & !HOBR
   /)
  
  
!     HENRYS LAW COEFFICIENTS
!     Henrys law coefficient
!     [KH298]=mole/(l atm)
!     Referencias em http://www.mpch-mainz.mpg.de/~sander/res/henry.html
!     * indica artigos nao encontrados nesse endereço eletronico
  REAL,PARAMETER,DIMENSION(nspecies) :: hstar=(/&
    1.10E-2              ,   & ! O3 - 001
    7.80E-4              ,   & ! H2 - 002
    8.30E+4              ,   & ! H2O2 - 003
    1.90E-3              ,   & ! NO - 004
    1.20E-2              ,   & ! NO2 - 005
    1.8E+00              ,   & ! NO3 - 006
    1.0E+10              ,   & ! N2O5 - 007
    5.00E+1              ,   & ! HONO - 008
    2.10E+5              ,   & ! HNO3 - 009
    1.20E+4              ,   & ! HNO4 - 010
    9.90E-4              ,   & ! CO - 011
    0.00E+0              ,   & ! O3P - 012
    0.00E+0              ,   & ! O1D - 013
    3.00E+1              ,   & ! HO - 014
    5.70E+3              ,   & ! HO2 - 015
    1.40E-3              ,   & ! CH4 - 016
    3.20E+3              ,   & ! HCHO - 017
    3.10E+2              ,   & ! OP1 - 018
    2.00E+3              ,   & ! MO2 - 019
    0.0E+00              ,   & ! CL - 020
    2.10E+4              ,   & ! HCL - 021
    0.0E+00              ,   & ! CLO - 022
    0.0E+00              ,   & ! CLN3 - 023
    1.20E+4              ,   & ! HOCL - 024
    0.0E+00              ,   & ! CLO2 - 025
    0.0E+00              ,   & ! CL2 - 026
    0.0E+00              ,   & ! BRO - 027
    0.0E+00              ,   & ! BR - 028
    2.10E+4              ,   & ! HBR - 029
    0.0E+00              ,   & ! BRN3 - 030
    1.20E+4                  & ! HOBR - 031
    /)
  
! [1] Noziere B. et al. The uptake of methyl vinyl ketone, methacrolein,
! and 2-methyl-3-butene-2-olonto sulfuric acid solutions,Journal of Physical
! Chemistry A, Vol.110, No.7, 2387-2395, 2006.
! [2] Abraham M. H. et al. Partition of compounds from gas to water and
! from gas to physiological saline at 310K: Linear free energy relationships,
! elsevier, 2006.
  REAL,PARAMETER,DIMENSION(nspecies) :: f0=(/&
    1.0 ,   & ! O3 - 001
    0.0 ,   & ! H2 - 002
    1.0 ,   & ! H2O2 - 003
    0.0 ,   & ! NO - 004
    0.1 ,   & ! NO2 - 005
    0.0 ,   & ! NO3 - 006
    0.0 ,   & ! N2O5 - 007
    0.0 ,   & ! HONO - 008
    0.0 ,   & ! HNO3 - 009
    0.0 ,   & ! HNO4 - 010
    0.0 ,   & ! CO - 011
    0.0 ,   & ! O3P - 012
    0.0 ,   & ! O1D - 013
    0.0 ,   & ! HO - 014
    0.0 ,   & ! HO2 - 015
    0.0 ,   & ! CH4 - 016
    0.0 ,   & ! HCHO - 017
    0.1 ,   & ! OP1 - 018
    0.0 ,   & ! MO2 - 019
    0.0 ,   & ! CL - 020
    0.0 ,   & ! HCL - 021
    0.0 ,   & ! CLO - 022
    0.0 ,   & ! CLN3 - 023
    0.0 ,   & ! HOCL - 024
    0.0 ,   & ! CLO2 - 025
    0.0 ,   & ! CL2 - 026
    0.0 ,   & ! BRO - 027
    0.0 ,   & ! BR - 028
    0.0 ,   & ! HBR - 029
    0.0 ,   & ! BRN3 - 030
    0.0     & ! HOBR - 031
    /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: difrat=(/&
    1.6 ,   & ! O3 - 001
    0.0 ,   & ! H2 - 002
    1.4 ,   & ! H2O2 - 003
    1.3 ,   & ! NO - 004
    1.6 ,   & ! NO2 - 005
    0.0 ,   & ! NO3 - 006
    0.0 ,   & ! N2O5 - 007
    0.0 ,   & ! HONO - 008
    1.9 ,   & ! HNO3 - 009
    0.0 ,   & ! HNO4 - 010
    0.0 ,   & ! CO - 011
    0.0 ,   & ! O3P - 012
    0.0 ,   & ! O1D - 013
    0.0 ,   & ! HO - 014
    0.0 ,   & ! HO2 - 015
    0.0 ,   & ! CH4 - 016
    1.3 ,   & ! HCHO - 017
    1.6 ,   & ! OP1 - 018
    0.0 ,   & ! MO2 - 019
    0.0 ,   & ! CL - 020
    0.0 ,   & ! HCL - 021
    0.0 ,   & ! CLO - 022
    0.0 ,   & ! CLN3 - 023
    0.0 ,   & ! HOCL - 024
    0.0 ,   & ! CLO2 - 025
    0.0 ,   & ! CL2 - 026
    0.0 ,   & ! BRO - 027
    0.0 ,   & ! BR - 028
    0.0 ,   & ! HBR - 029
    0.0 ,   & ! BRN3 - 030
    0.0     & ! HOBR - 031
    /)
  
  
!     DIFFUSION COEFFICIENTS
!     [DV]=cm2/s (assumed: 1/SQRT(molar mass) when not known)
  REAL,PARAMETER,DIMENSION(nspecies) :: dvj=(/&
    0.1750000	    ,   & ! O3 - 001
    0.7070000     ,   & ! H2 - 002
    0.1710000	    ,   & ! H2O2 - 003
    0.1830000	    ,   & ! NO - 004
    0.1470000	    ,   & ! NO2 - 005
    0.1270000	    ,   & ! NO3 - 006
    0.1100000	    ,   & ! N2O5 - 007
    0.1530000	    ,   & ! HONO - 008
    0.1260000	    ,   & ! HNO3 - 009
    0.1130000	    ,   & ! HNO4 - 010
    0.1890000	    ,   & ! CO - 011
    0.2500000	    ,   & ! O3P - 012
    0.2500000	    ,   & ! O1D - 013
    0.2430000	    ,   & ! HO - 014
    0.1740000	    ,   & ! HO2 - 015
    0.2500000	    ,   & ! CH4 - 016
    0.1830000	    ,   & ! HCHO - 017
    0.1440000	    ,   & ! OP1 - 018
    0.1458650	    ,   & ! MO2 - 019
    0.1680000     ,   & ! CL - 020
    0.1660000     ,   & ! HCL - 021
    0.1390000     ,   & ! CLO - 022
    0.1010000     ,   & ! CLN3 - 023
    0.1380000     ,   & ! HOCL - 024
    0.1040000     ,   & ! CLO2 - 025
    0.1190000     ,   & ! CL2 - 026
    0.1020000     ,   & ! BRO - 027
    0.1120000     ,   & ! BR - 028
    0.1110000     ,   & ! HBR - 029
    0.0839000     ,   & ! BRN3 - 030
    0.1020000         & ! HOBR - 031
    /)
  
  
!     -DH/R (for temperature correction)
!     [-DH/R]=K
  REAL,PARAMETER,DIMENSION(nspecies) :: dhr=(/&
    2400.         ,   & ! O3 - 001
    500.          ,   & ! H2 - 002
    7400.         ,   & ! H2O2 - 003
    1400.         ,   & ! NO - 004
    2500.         ,   & ! NO2 - 005
    2000.         ,   & ! NO3 - 006
    3400.         ,   & ! N2O5 - 007
    4900.         ,   & ! HONO - 008
    8700.         ,   & ! HNO3 - 009
    6900.         ,   & ! HNO4 - 010
    1300.         ,   & ! CO - 011
    0.            ,   & ! O3P - 012
    0.            ,   & ! O1D - 013
    4500.         ,   & ! HO - 014
    5900.         ,   & ! HO2 - 015
    1600.         ,   & ! CH4 - 016
    6800.         ,   & ! HCHO - 017
    5200.         ,   & ! OP1 - 018
    6600.         ,   & ! MO2 - 019
    0.            ,   & ! CL - 020
    8700.         ,   & ! HCL - 021
    0.0           ,   & ! CLO - 022
    0.0           ,   & ! CLN3 - 023
    6900.         ,   & ! HOCL - 024
    0.0           ,   & ! CLO2 - 025
    0.0           ,   & ! CL2 - 026
    0.0           ,   & ! BRO - 027
    0.0           ,   & ! BR - 028
    8700.         ,   & ! HBR - 029
    0.0           ,   & ! BRN3 - 030
    6900.             & ! HOBR - 031
    /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: weight=(/&
    48.  ,   & ! O3 - 001
    2.   ,   & ! H2 - 002
    34.  ,   & ! H2O2 - 003
    30.  ,   & ! NO - 004
    46.  ,   & ! NO2 - 005
    62.  ,   & ! NO3 - 006
    108. ,   & ! N2O5 - 007
    47.  ,   & ! HONO - 008
    63.  ,   & ! HNO3 - 009
    79.  ,   & ! HNO4 - 010
    28.  ,   & ! CO - 011
    16.  ,   & ! O3P - 012
    16.  ,   & ! O1D - 013
    17.  ,   & ! HO - 014
    33.  ,   & ! HO2 - 015
    16.  ,   & ! CH4 - 016
    30.  ,   & ! HCHO - 017
    48.  ,   & ! OP1 - 018
    47.  ,   & ! MO2 - 019
    35.5 ,   & ! CL - 020
    36.5 ,   & ! HCL - 021
    51.5 ,   & ! CLO - 022
    97.5 ,   & ! CLN3 - 023
    52.5 ,   & ! HOCL - 024
    92.  ,   & ! CLO2 - 025
    71.  ,   & ! CL2 - 026
    96.  ,   & ! BRO - 027
    80.  ,   & ! BR - 028
    81.  ,   & ! HBR - 029
    142. ,   & ! BRN3 - 030
    97.      & ! HOBR - 031
   /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: init_ajust=(/&
    1.0 ,   & ! O3 - 001
    1.0 ,   & ! H2 - 002
    1.0 ,   & ! H2O2 - 003
    1.0 ,   & ! NO - 004
    1.0 ,   & ! NO2 - 005
    1.0 ,   & ! NO3 - 006
    1.0 ,   & ! N2O5 - 007
    1.0 ,   & ! HONO - 008
    1.0 ,   & ! HNO3 - 009
    1.0 ,   & ! HNO4 - 010
    1.0 ,   & ! CO - 011
    1.0 ,   & ! O3P - 012
    1.0 ,   & ! O1D - 013
    1.0 ,   & ! HO - 014
    1.0 ,   & ! HO2 - 015
    1.0 ,   & ! CH4 - 016
    1.0 ,   & ! HCHO - 017
    1.0 ,   & ! OP1 - 018
    1.0 ,   & ! MO2 - 019
    1.0 ,   & ! CL - 020
    1.0 ,   & ! HCL - 021
    1.0 ,   & ! CLO - 022
    1.0 ,   & ! CLN3 - 023
    1.0 ,   & ! HOCL - 024
    1.0 ,   & ! CLO2 - 025
    1.0 ,   & ! CL2 - 026
    1.0 ,   & ! BRO - 027
    1.0 ,   & ! BR - 028
    1.0 ,   & ! HBR - 029
    1.0 ,   & ! BRN3 - 030
    1.0     & ! HOBR - 031
   /)
  
  
  REAL,PARAMETER,DIMENSION(nspecies) :: emiss_ajust=(/&
    1.0 ,   & ! O3 - 001
    1.0 ,   & ! H2 - 002
    1.0 ,   & ! H2O2 - 003
    1.0 ,   & ! NO - 004
    1.0 ,   & ! NO2 - 005
    1.0 ,   & ! NO3 - 006
    1.0 ,   & ! N2O5 - 007
    1.0 ,   & ! HONO - 008
    1.0 ,   & ! HNO3 - 009
    1.0 ,   & ! HNO4 - 010
    1.0 ,   & ! CO - 011
    1.0 ,   & ! O3P - 012
    1.0 ,   & ! O1D - 013
    1.0 ,   & ! HO - 014
    1.0 ,   & ! HO2 - 015
    1.0 ,   & ! CH4 - 016
    1.0 ,   & ! HCHO - 017
    1.0 ,   & ! OP1 - 018
    1.0 ,   & ! MO2 - 019
    1.0 ,   & ! CL - 020
    1.0 ,   & ! HCL - 021
    1.0 ,   & ! CLO - 022
    1.0 ,   & ! CLN3 - 023
    1.0 ,   & ! HOCL - 024
    1.0 ,   & ! CLO2 - 025
    1.0 ,   & ! CL2 - 026
    1.0 ,   & ! BRO - 027
    1.0 ,   & ! BR - 028
    1.0 ,   & ! HBR - 029
    1.0 ,   & ! BRN3 - 030
    1.0     & ! HOBR - 031
   /)
  
  
!    ACID DISSOCIATION CONSTANT AT 298K 
!     [mole/liter of liquid water]
!     Referencias: Barth et al. JGR 112, D13310 2007
!     Martell and Smith, 1976, Critical stability
!     vol1-4 Plenum Press New York
  REAL,PARAMETER,DIMENSION(nspecies) :: ak0=(/&
    0.00E+00	    ,   & ! O3 - 001
    0.00E+00     ,   & ! H2 - 002
    2.20E-12     ,   & ! H2O2 - 003
    0.00E+00	    ,   & ! NO - 004
    0.00E+00	    ,   & ! NO2 - 005
    0.00E+00	    ,   & ! NO3 - 006
    0.00E+00	    ,   & ! N2O5 - 007
    7.10E-04	    ,   & ! HONO - 008
    1.54E+01	    ,   & ! HNO3 - 009
    0.00E+00	    ,   & ! HNO4 - 010
    0.00E+00	    ,   & ! CO - 011
    0.00E+00	    ,   & ! O3P - 012
    0.00E+00	    ,   & ! O1D - 013
    0.00E+00	    ,   & ! HO - 014
    3.50E-05	    ,   & ! HO2 - 015
    0.00E+00	    ,   & ! CH4 - 016
    0.00E+00	    ,   & ! HCHO - 017
    0.00E+00	    ,   & ! OP1 - 018
    0.00E+00	    ,   & ! MO2 - 019
    0.00E+00     ,   & ! CL - 020
    1.54E+01     ,   & ! HCL - 021
    0.00E+00     ,   & ! CLO - 022
    0.00E+00     ,   & ! CLN3 - 023
    0.00E+00     ,   & ! HOCL - 024
    0.00E+00     ,   & ! CLO2 - 025
    0.00E+00     ,   & ! CL2 - 026
    0.00E+00     ,   & ! BRO - 027
    0.00E+00     ,   & ! BR - 028
    1.54E+01     ,   & ! HBR - 029
    0.00E+00     ,   & ! BRN3 - 030
    0.00E+00         & ! HOBR - 031
   /)
  
  
!     Temperature correction factor for
!     acid dissociation constants
!     [K]
!     Referencias: Barth et al. JGR 112, D13310 2007
  REAL,PARAMETER,DIMENSION(nspecies) :: dak=(/&
    0.         ,   & ! O3 - 001
    0.         ,   & ! H2 - 002
    -3700.     ,   & ! H2O2 - 003
    0.         ,   & ! NO - 004
    0.         ,   & ! NO2 - 005
    0.         ,   & ! NO3 - 006
    0.         ,   & ! N2O5 - 007
    0.         ,   & ! HONO - 008
    0.         ,   & ! HNO3 - 009
    0.         ,   & ! HNO4 - 010
    0.         ,   & ! CO - 011
    0.         ,   & ! O3P - 012
    0.         ,   & ! O1D - 013
    0.         ,   & ! HO - 014
    0.         ,   & ! HO2 - 015
    0.         ,   & ! CH4 - 016
    0.         ,   & ! HCHO - 017
    0.         ,   & ! OP1 - 018
    0.         ,   & ! MO2 - 019
    0.         ,   & ! CL - 020
    0.         ,   & ! HCL - 021
    0.         ,   & ! CLO - 022
    0.         ,   & ! CLN3 - 023
    0.         ,   & ! HOCL - 024
    0.         ,   & ! CLO2 - 025
    0.         ,   & ! CL2 - 026
    0.         ,   & ! BRO - 027
    0.         ,   & ! BR - 028
    0.         ,   & ! HBR - 029
    0.         ,   & ! BRN3 - 030
    0.             & ! HOBR - 031
    /)
  
  
  INTEGER,PARAMETER :: nr   =096!Number of gas-phase reactions
  INTEGER,PARAMETER :: nrt  =096!Total Number of reactions
  INTEGER,PARAMETER :: nrh2o=000!Number of aqueous-phase reactions
  
 !------------------------------------------------------------------------------  
  
 ! Photolysis Rate Calculation: method used (LUT=look_up_table, FAST-JX= on-line)
  
  CHARACTER(LEN=10),PARAMETER :: PhotojMethod= 'FAST-TUV'
  INTEGER,PARAMETER :: maxJcomb=    5, nr_photo=   19
  INTEGER,PARAMETER,DIMENSION(nr_photo) :: nfactors=(/ &
                            1, &!  1
                            1, &!  2
                            1, &!  3
                            1, &!  4
                            1, &!  5
                            1, &!  6
                            1, &!  7
                            1, &!  8
                            1, &!  9
                            1, &! 10
                            1, &! 11
                            1, &! 12
                            1, &! 13
                            1, &! 14
                            1, &! 15
                            1, &! 16
                            1, &! 17
                            1, &! 18
                            1/)! 19
 
  DOUBLE PRECISION,PARAMETER,DIMENSION(maxJcomb,nr_photo) :: factor=RESHAPE((/ &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.900D+00, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00, &
                                                0.100D+01, 0.000D+00, 0.000D+00, 0.000D+00, 0.000D+00/),(/maxJcomb,nr_photo/))
  CHARACTER(LEN=07 ),PARAMETER,DIMENSION(maxJcomb,nr_photo) :: JReactionComp=RESHAPE((/ &
                                               "NO2    ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "O3(1D) ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "O3     ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "HNO2   ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "HNO3   ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "HNO4   ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "NO3    ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "NO3    ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "H2O2   ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "H2COb  ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "H2COa  ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "CH3OOH ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "ClNO3a ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "ClNO3b ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "Cl2O2  ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "Cl2    ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "BrO    ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "BrNO3  ","NONE   ","NONE   ","NONE   ","NONE   ", &
                                               "HOBr   ","NONE   ","NONE   ","NONE   ","NONE   "/),(/maxJcomb,nr_photo/))

  
END MODULE chem1_list
