PROGRAM post
 !Post-processor for spack
  !Author: Luiz Flavio Rodrigues
  !        23/10/2007
  !Read some files from Spack and then serializes the arrays
  !
  CHARACTER(LEN=*),PARAMETER :: file1="chem_spack_fexchem.f90"
  CHARACTER(LEN=*),PARAMETER :: file1_s="chem_spack_fexchem_s.f90"
  CHARACTER(LEN=*),PARAMETER :: tipo1="chem"
  CHARACTER(LEN=*),PARAMETER :: file2="chem_spack_fexloss.f90"
  CHARACTER(LEN=*),PARAMETER :: file2_s="chem_spack_fexloss_s.f90"
  CHARACTER(LEN=*),PARAMETER :: tipo2="loss"
  CHARACTER(LEN=*),PARAMETER :: file3="chem_spack_fexprod.f90"
  CHARACTER(LEN=*),PARAMETER :: file3_s="chem_spack_fexprod_s.f90"
  CHARACTER(LEN=*),PARAMETER :: tipo3="prod"
  CHARACTER(LEN=*),PARAMETER :: file4="chem_spack_jacdchemdc.f90"
  CHARACTER(LEN=*),PARAMETER :: file4_s="chem_spack_jacdchemdc_s.f90"

  CALL post_jacc(file4,file4_s) 
  CALL post_convert(file1,file1_s,tipo1)
  CALL post_convert(file2,file2_s,tipo2)
  CALL post_convert(file3,file3_s,tipo3)

END PROGRAM Post
 
SUBROUTINE post_jacc(file1,file1_s)
  USE chem1_list, ONLY: nspecies,nr
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: file1
  CHARACTER(LEN=*),INTENT(IN) :: file1_s
 
  CHARACTER(LEN=45) :: aux
  CHARACTER(LEN=150) :: linha
  CHARACTER(LEN=150),DIMENSION(10000) :: buffer
  LOGICAL,DIMENSION(10000) :: item_add
  CHARACTER(LEN=17) :: item
  CHARACTER(LEN=15) :: lixo1
  CHARACTER(LEN=18) :: lixo2
  CHARACTER(LEN=3) :: num
  INTEGER,DIMENSION(nspecies) :: fexprod_count
  INTEGER :: istat,ngas,cl,buf_count,i,j
  LOGICAL :: laco=.false.

  fexprod_count=0
  cl=0

  OPEN(unit=22,file=file1)
  OPEN(unit=55,file=file1_s)
  
  DO

    IF(.not. laco) THEN
      !PRINT *,'Reading jacc.....',cl
      cl=cl+1
      READ(22,FMT='(A150)',IOSTAT=istat) linha
      IF(istat/=0) EXIT
      IF(linha(1:17)=='        JacC=0.d0') CYCLE
      WRITE (55,FMT='(A150)') linha
    END IF
    IF(linha(7:26)=='DO ijk=ijkbeg,ijkend' .and. .not. laco) THEN
       print *,'Laco encontrado Jacc.....'
       laco=.true.
       buf_count=0
       CYCLE
    END IF
    
    IF(laco) THEN
      !PRINT *,'Reading Jacc.....',cl
      cl=cl+1
      buf_count=buf_count+1
      READ(22,FMT='(A150)') buffer(buf_count)
      IF(buffer(buf_count)(8:13)=='END DO' .and. laco) THEN
	item_add=.false.	
        !=====================================
        DO i=1,buf_count-1
          IF(item_add(i)) CYCLE !Já computado
          item=buffer(i)(7:23)
          item_add(i)=.true.
          WRITE (55,FMT='(A,$)') buffer(i)(1:26)//trim(buffer(i)(44:))      
          DO j=i+1,buf_count-1
            IF(item==buffer(j)(7:23) .and. .not. item_add(j)) THEN
              item_add(j)=.true.
              WRITE (55,FMT='(A)') ' &'
              WRITE (55,FMT='(A, $)') repeat(' ',26)//trim(buffer(j)(44:)) 
            END IF
          END DO
          WRITE(55,FMT='(A)') ''
        END DO     
        !======================================       
              
        WRITE (55,FMT='(A150)') buffer(buf_count)
        laco=.false.
        CYCLE
      END IF
    END IF
    
  END DO

  CLOSE(unit=22)
  CLOSE(unit=55)

END SUBROUTINE post_jacc 

SUBROUTINE post_convert(file1,file1_s,tipo)
  USE chem1_list, ONLY: nspecies,nr
  IMPLICIT NONE
  CHARACTER(LEN=*),INTENT(IN) :: file1
  CHARACTER(LEN=*),INTENT(IN) :: file1_s
  CHARACTER(LEN=*),INTENT(IN) :: tipo
  CHARACTER(LEN=45),DIMENSION(nspecies,nr) :: fexprod
  CHARACTER(LEN=45) :: aux
  CHARACTER(LEN=80) :: linha
  CHARACTER(LEN=15) :: lixo1
  CHARACTER(LEN=18) :: lixo2
  CHARACTER(LEN=3) :: num
  INTEGER,DIMENSION(nspecies) :: fexprod_count
  INTEGER :: istat,ngas,cl
  LOGICAL :: laco

  fexprod_count=0
  cl=0
  fexprod(:,:)='empty'
  laco=.false.
  

  OPEN(unit=22,file=file1)
  OPEN(unit=55,file=file1_s)
  PRINT *,'Arquivo:'//trim(file1)
  DO

    IF(.not. laco) THEN
      PRINT *,'Reading (Fora do laco).....',cl
      cl=cl+1
      READ(22,FMT='(A80)',IOSTAT=istat) linha
      IF(istat/=0) EXIT
      IF(linha(1:17)=='       '//tipo//'=0.d0') CYCLE
      WRITE (55,FMT='(A80)') linha
    END IF

    IF(linha(1:29)=='      DO ijk=ijkbegin,ijkend' .and. .not. laco) THEN
       print *,'Laco encontrado.....'
       laco=.true.
       CYCLE
    END IF
    
    IF(linha(1:12)=='      END DO' .and. laco) THEN
      print *,'Encerrando compactacao....'
      CALL fecha_laco(nspecies,nr,fexprod,fexprod_count)
      WRITE (55,FMT='(A80)') linha
      laco=.false.
      CYCLE
    END IF

    IF(laco) THEN
      PRINT *,'Reading.....',cl
      cl=cl+1
      READ(22,FMT='(A15,I3,A18,A45)') lixo1,ngas,lixo2,aux
      IF(lixo1(1:12)=='      END DO') THEN
        print *,'Encerrando compactacao....'
        CALL fecha_laco(nspecies,nr,fexprod,fexprod_count,tipo)
        WRITE (55,FMT='(A15)') lixo1
        laco=.false.
        CYCLE
      END IF        
      fexprod_count(ngas)=fexprod_count(ngas)+1
      fexprod(ngas,fexprod_count(ngas))=aux
    END IF           
    
  END DO

  CLOSE(unit=22)
  CLOSE(unit=55)

END SUBROUTINE post_convert

SUBROUTINE fecha_laco(nspecies,nr,fexprod,fexprod_count,tipo)
  IMPLICIT NONE
  INTEGER,INTENT(IN) :: nspecies,nr
  CHARACTER(LEN=*),INTENT(IN) :: tipo
  CHARACTER(LEN=45),DIMENSION(nspecies,nr),INTENT(IN) :: fexprod
  INTEGER,DIMENSION(nspecies),INTENT(IN) :: fexprod_count
  INTEGER :: i,j


  DO i=1,nspecies
     if(trim(fexprod(i,1)) == 'empty') cycle
     WRITE(55,FMT='(A15,I3,A4,A,$)') '      '//tipo//'(ijk,',i,') = ',trim(fexprod(i,1))
     !print*,' xx',i,trim(fexprod(i,1))
     IF(fexprod_count(i)>1) THEN 
         WRITE(55,FMT='(A)') ' &'
     ELSE
       WRITE(55,FMT='(A)') ''
       CYCLE
     END IF
     DO j=2,fexprod_count(i)-1
       WRITE(55,FMT='(A)') '          '//trim(fexprod(i,j))//' &'
     END DO
     WRITE(55,FMT='(A)') '          '//trim(fexprod(i,fexprod_count(i)))
  END DO

END SUBROUTINE fecha_laco

        

