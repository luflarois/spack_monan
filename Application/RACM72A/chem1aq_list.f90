MODULE chem1aq_list
  IMPLICIT NONE
  
  
  INTEGER,PARAMETER :: maxnspeciesaq= 200
  INTEGER,PARAMETER :: nspeciesaq=020
  
  
  !Name of speciesaq 
  CHARACTER(LEN=8),PARAMETER,DIMENSION(nspeciesaq) :: spcaq_name=(/ &
   'O3aq  ' & !
     ,'H2O2aq' & !
     ,'HONOaq' & !
     ,'HNO3aq' & !
     ,'HNO4aq' & !
     ,'SO2aq ' & !
     ,'SULFaq' & !
     ,'HOaq  ' & !
     ,'HO2aq ' & !
     ,'CSLaq ' & !
     ,'HCHOaq' & !
     ,'GLYaq ' & !
     ,'MGLYaq' & !
     ,'DCBaq ' & !
     ,'OP1aq ' & !
     ,'OP2aq ' & !
     ,'PAAaq ' & !
     ,'ORA1aq' & !
     ,'ORA2aq' & !
     ,'MO2aq ' & !
   /)
  
  
  !Number of each specie   
  INTEGER,PARAMETER :: O3aq  =001
  INTEGER,PARAMETER :: H2O2aq=002
  INTEGER,PARAMETER :: HONOaq=003
  INTEGER,PARAMETER :: HNO3aq=004
  INTEGER,PARAMETER :: HNO4aq=005
  INTEGER,PARAMETER :: SO2aq =006
  INTEGER,PARAMETER :: SULFaq=007
  INTEGER,PARAMETER :: HOaq  =008
  INTEGER,PARAMETER :: HO2aq =009
  INTEGER,PARAMETER :: CSLaq =010
  INTEGER,PARAMETER :: HCHOaq=011
  INTEGER,PARAMETER :: GLYaq =012
  INTEGER,PARAMETER :: MGLYaq=013
  INTEGER,PARAMETER :: DCBaq =014
  INTEGER,PARAMETER :: OP1aq =015
  INTEGER,PARAMETER :: OP2aq =016
  INTEGER,PARAMETER :: PAAaq =017
  INTEGER,PARAMETER :: ORA1aq=018
  INTEGER,PARAMETER :: ORA2aq=019
  INTEGER,PARAMETER :: MO2aq =020
  
  
!     number of the corresponding gaseous species
  INTEGER,PARAMETER,DIMENSION(nspeciesaq) :: ind_gas=(/&
    001 ,   & ! O3aq - 001
    002 ,   & ! H2O2aq - 002
    007 ,   & ! HONOaq - 003
    008 ,   & ! HNO3aq - 004
    009 ,   & ! HNO4aq - 005
    010 ,   & ! SO2aq - 006
    011 ,   & ! SULFaq - 007
    015 ,   & ! HOaq - 008
    016 ,   & ! HO2aq - 009
    031 ,   & ! CSLaq - 010
    032 ,   & ! HCHOaq - 011
    035 ,   & ! GLYaq - 012
    036 ,   & ! MGLYaq - 013
    037 ,   & ! DCBaq - 014
    044 ,   & ! OP1aq - 015
    045 ,   & ! OP2aq - 016
    046 ,   & ! PAAaq - 017
    047 ,   & ! ORA1aq - 018
    048 ,   & ! ORA2aq - 019
    049     & ! MO2aq - 020
    /)
  
! accomodation coefficient
  REAL,PARAMETER,DIMENSION(nspeciesaq) :: acco=(/&
    5.3E-4               ,   & ! O3aq - 001
    2.0E-1               ,   & ! H2O2aq - 002
    5.0E-2               ,   & ! HONOaq - 003
    2.0E-1               ,   & ! HNO3aq - 004
    2.0E-1               ,   & ! HNO4aq - 005
    5.0E-2               ,   & ! SO2aq - 006
    2.0E-1               ,   & ! SULFaq - 007
    5.0E-2               ,   & ! HOaq - 008
    2.0E-1               ,   & ! HO2aq - 009
    1.0E-1               ,   & ! CSLaq - 010
    2.0E-1               ,   & ! HCHOaq - 011
    2.0E-1               ,   & ! GLYaq - 012
    2.0E-1               ,   & ! MGLYaq - 013
    2.0E-1               ,   & ! DCBaq - 014
    5.0E-2               ,   & ! OP1aq - 015
    2.0E-1               ,   & ! OP2aq - 016
    1.0E-1               ,   & ! PAAaq - 017
    2.0E-1               ,   & ! ORA1aq - 018
    2.0E-2               ,   & ! ORA2aq - 019
    2.0E-1                   & ! MO2aq - 020
    /)
  
END MODULE chem1aq_list
