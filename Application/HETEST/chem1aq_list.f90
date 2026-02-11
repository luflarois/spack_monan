MODULE chem1aq_list
  IMPLICIT NONE
  
  
  INTEGER,PARAMETER :: maxnspeciesaq= 200
  INTEGER,PARAMETER :: nspeciesaq=014
  
  
  !Name of speciesaq 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspeciesaq) :: spcaq_name=(/ &
   'O3aq  ' & !
     ,'H2O2aq' & !
     ,'HONOaq' & !
     ,'HNO3aq' & !
     ,'HNO4aq' & !
     ,'HOaq  ' & !
     ,'HO2aq ' & !
     ,'HCHOaq' & !
     ,'OP1aq ' & !
     ,'MO2aq ' & !
     ,'HCLaq ' & !
     ,'HOCLaq' & !
     ,'HBRaq ' & !
     ,'HOBRaq' & !
   /)
  
  
  !Number of each specie   
  INTEGER,PARAMETER :: O3aq  =001
  INTEGER,PARAMETER :: H2O2aq=002
  INTEGER,PARAMETER :: HONOaq=003
  INTEGER,PARAMETER :: HNO3aq=004
  INTEGER,PARAMETER :: HNO4aq=005
  INTEGER,PARAMETER :: HOaq  =006
  INTEGER,PARAMETER :: HO2aq =007
  INTEGER,PARAMETER :: HCHOaq=008
  INTEGER,PARAMETER :: OP1aq =009
  INTEGER,PARAMETER :: MO2aq =010
  INTEGER,PARAMETER :: HCLaq =011
  INTEGER,PARAMETER :: HOCLaq=012
  INTEGER,PARAMETER :: HBRaq =013
  INTEGER,PARAMETER :: HOBRaq=014
  
  
!     number of the corresponding gaseous species
  INTEGER,PARAMETER,DIMENSION(nspeciesaq) :: ind_gas=(/&
    001 ,   & ! O3aq - 001
    003 ,   & ! H2O2aq - 002
    008 ,   & ! HONOaq - 003
    009 ,   & ! HNO3aq - 004
    010 ,   & ! HNO4aq - 005
    014 ,   & ! HOaq - 006
    015 ,   & ! HO2aq - 007
    017 ,   & ! HCHOaq - 008
    018 ,   & ! OP1aq - 009
    019 ,   & ! MO2aq - 010
    021 ,   & ! HCLaq - 011
    024 ,   & ! HOCLaq - 012
    029 ,   & ! HBRaq - 013
    031     & ! HOBRaq - 014
    /)
  
! accomodation coefficient
  REAL,PARAMETER,DIMENSION(nspeciesaq) :: acco=(/&
    5.3E-4               ,   & ! O3aq - 001
    2.0E-1               ,   & ! H2O2aq - 002
    5.0E-2               ,   & ! HONOaq - 003
    2.0E-1               ,   & ! HNO3aq - 004
    2.0E-1               ,   & ! HNO4aq - 005
    5.0E-2               ,   & ! HOaq - 006
    2.0E-1               ,   & ! HO2aq - 007
    2.0E-1               ,   & ! HCHOaq - 008
    5.0E-2               ,   & ! OP1aq - 009
    2.0E-1               ,   & ! MO2aq - 010
    2.0E-1               ,   & ! HCLaq - 011
    2.0E-1               ,   & ! HOCLaq - 012
    2.0E-1               ,   & ! HBRaq - 013
    2.0E-1                   & ! HOBRaq - 014
    /)
  
END MODULE chem1aq_list
