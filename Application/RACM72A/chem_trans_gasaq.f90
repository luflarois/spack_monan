MODULE mod_chem_trans_gasaq

  IMPLICIT NONE

  PRIVATE

  PUBLIC :: trans_gaq ! Subroutine

CONTAINS


  !======================================================================
  SUBROUTINE trans_gaq(m1,m2,m3,ia,iz,ja,jz)

    !====================================================================
    !
    !  author M. Pirre 18/02/2008 
    !
    !==================================================================
    !
    USE mem_chem1   , ONLY:chem1_g
    USE mem_chem1aq , ONLY:chem1aq_g
    USE mem_grid    , ONLY:ngrids,ngrid
    USE mem_basic   , ONLY:basic_g
    USE mem_micro   , ONLY:micro_g
    USE chem1aq_list, ONLY:nspeciesaq,ind_gas
    IMPLICIT NONE
    INTEGER ngr,iaq
    INTEGER i,j,k,ijk
    INTEGER, INTENT(IN) ::m1,m2,m3,ia,iz,ja,jz
    REAL, dimension(m1,m2,m3)::temp,volc,volr
    !
    ngr = ngrid

    call trans_para(m1,m2,m3,ia,iz,ja,jz  &
         ! parameters to compute temperature and pressure
         ,basic_g(ngr)%theta     (1,1,1)     &
         ,basic_g(ngr)%pi0       (1,1,1)     &
         ,basic_g(ngr)%pp        (1,1,1)     &
         ,micro_g(ngr)%rrp       (1,1,1)     &! rain mixing ratio (kg(water)/kg (air))          
         ,micro_g(ngr)%rcp       (1,1,1)     &! cloud mixing ratio (kg(water)/kg (air))
         ,temp,volc,volr)                     ! output: temperature, cloud water and rain mixing ratio

    !        
    DO iaq = 1,nspeciesaq
      call trans_tend(m1,m2,m3,ia,iz,ja,jz,iaq   &
            ,temp,volc,volr                      &
            ,chem1_g(ind_gas(iaq),ngr)%sc_p      &! gas species mixing ratio (kg(species)/kg(air)  
            ,chem1aq_g(iaq,ngr)%sc_pc            &! mixing ratio of species in cloud and rain drops
            ,chem1aq_g(iaq,ngr)%sc_pr            &! (kg(species)/kg(air)
            ,chem1_g(ind_gas(iaq),ngr)%sc_t      &! tendency of gas species (kg(species)/kg(air).s^-1)    
            ,chem1aq_g(iaq,ngr)%sc_tc            &! tendencies of mixing ratio of species in cloud and rain drop
            ,chem1aq_g(iaq,ngr)%sc_tr            &! (kg(species)/kg(air).s^-1)
            )
    ENDDO

  end subroutine trans_gaq

  !======================================================================
  SUBROUTINE trans_tend(m1,m2,m3,ia,iz,ja,jz,iaq,       &
       temp,volc,volr,                            &
       sc_p, sc_pc,sc_pr, sc_t, sc_tc, sc_tr)

    !====================================================================
    ! author: M. Pirre  18/02/2008
    !
    ! Computation of the tendencies of the mixing ratio of the chemical species
    ! in the gas phase: sc_t (kg(species).kg(air)^-1.s^-1), 
    ! in the liquid cloud phase: sc_tc (kg(species).kg(air)^-1.s^-1)
    ! and in the liquid rain phase: sc_tr (kg(species).kg(air)^1.s^-1)
    ! due to mass tranfer between the gas phase and the liquid phase

    ! reference: Barth et al., JGR, vol 106, D12, p 12,381-12,400, 2001
    !
    !=====================================================================
    !
    USE chem1_list  , ONLY:nspecies,weight,hstar,dvj,dhr,ak0,dak
    USE chem1aq_list, ONLY:spcaq_name,ind_gas,acco
    USE mem_grid    , ONLY:grid_g,ngrid,time,dtlt,itime1,zt
    USE rconstants  , ONLY:rgas

    IMPLICIT NONE
    INTEGER m1,m2,m3,ia,iz,ja,jz
    REAL pi,weightair,patm,rmol,ratm,rcloud,rrain,temp0i,vit,tauic  &
         ,tauir,tcorr,corrhc,corrhr,taudc,taudr,tend_c,tend_g ,tend_r
    REAL hplusc,hplusr,henry_fc,henry_fr,henry_t,ktc,ktr,ctot,cgash  &
         ,cliqh,acliqh,acgash,accloudh,acrainh,voltot
    REAL, dimension(m1,m2,m3)  :: temp, volc, volr   &
         , sc_p, sc_pc, sc_pr, sc_t, sc_tc, sc_tr
    INTEGER i,j,k,iaq
    INTEGER itest 
    ! define some constants
    parameter (pi  = 3.14159265)         !pi number
    parameter (weightair = 28.97E-3)     !air molar weight in kg
    parameter (patm = 1.01325E+5)        !1 atmosphere in Pa
    parameter (rmol = rgas*weightair)    !universal gas constant: 8.314
    parameter (ratm = rmol*1000./patm)   !universal gas constant in dm3.atm.K-1.mole-1: 0.082057
    parameter (rcloud = 1.5e-3)          !mean cloud radius in cm
    parameter (rrain  = 1.e-1)           !mean raindrop radius in cm
    parameter (temp0i=1./298.15)         !inverse of standard temperature in K-1


    ! gas diffusion timescale toward cloud (taudc) and rain (taudr) drop  
    taudc=(rcloud*rcloud)/ (3.*dvj(ind_gas(iaq)))
    taudr=(rrain*rrain)  / (3.*dvj(ind_gas(iaq)))
    !      
    hplusc=1.175E-4 !pH=3.93 RAMS chimie 
    hplusr= 7.95E-5 !pH=4.10 RAMS chimie
    itest=0
    DO j = ja , jz
       DO i = ia , iz
          DO k =  2,m1-1
             if(volc(k,i,j) <= 0. .and. volr(k,i,j) <= 0.) cycle 
             ! timescale for transfer across the interface of the liquid cloud drop (tauic) 
             ! and the liquid rain drop (tauir)
             vit=100.*sqrt(8000.*rmol*temp(k,i,j)/(pi*weight(ind_gas(iaq))))
             tauic=4.*rcloud/(3.*vit*acco(iaq))
             tauir=4.*rrain /(3.*vit*acco(iaq))
             !              
             ! rate constant (cm^3(air).cm^-3(water).s^-1 for transfer through the gas phase
             ! and across the interface of cloud drop (ktc) and rain drop (ktr)              
             ktc  =1./(taudc+tauic)
             ktr  =1./(taudr+tauir)
             !              
             ! effective henry constant for cloud and rain * rmol* temp: henry_fc , henry_fr
             ! in cm^3(air).cm^-3(water)             
             tcorr=1./temp(k,i,j)-temp0i
             corrhc=1.+ak0(ind_gas(iaq))*exp(dak(ind_gas(iaq))*tcorr)/hplusc
             henry_fc=hstar(ind_gas(iaq))*exp(dhr(ind_gas(iaq))*tcorr)*corrhc&
                  *ratm*temp(k,i,j)
             corrhr=1.+ak0(ind_gas(iaq))*exp(dak(ind_gas(iaq))*tcorr)/hplusr
             henry_fr=hstar(ind_gas(iaq))*exp(dhr(ind_gas(iaq))*tcorr)*corrhr&
                  *ratm*temp(k,i,j)
             !                      
             !              
             ! tendencies for species in gas, cloud and rain drops             
             tend_c=ktc*volc(k,i,j)*sc_p(k,i,j)-ktc*sc_pc(k,i,j)/henry_fc
             tend_r=ktr*volr(k,i,j)*sc_p(k,i,j)-ktr*sc_pr(k,i,j)/henry_fr
             tend_g=-tend_c-tend_r
             !              
             ! mixing ratio of gas species at time t + dtlt - mixing ratio at time t
             ! assuming Henry equilibrium              
             ctot=sc_p(k,i,j)+sc_pc(k,i,j)+sc_pr(k,i,j)
             henry_t=henry_fc*volc(k,i,j)+henry_fr*volr(k,i,j)
             cgash=ctot/(1+henry_t)
             acgash=cgash-sc_p(k,i,j)
             cliqh=henry_t*ctot/(1+henry_t)
             acliqh=cliqh-sc_pc(k,i,j)-sc_pr(k,i,j)
             if(abs(tend_g)*dtlt .le. abs(acliqh))then
                ! cumulative tendencies assuming the mass tranfer equations                      
                sc_t(k,i,j)=sc_t(k,i,j)+tend_g
                sc_tc(k,i,j)=sc_tc(k,i,j)+tend_c
                sc_tr(k,i,j)=sc_tr(k,i,j)+tend_r
             !   if((k.eq.14).and.(i.eq.13).and.(j.eq.26)) print*,'dum tend',iaq,k,i,j,sc_t(k,i,j),sc_tc(k,i,j),sc_tr(k,i,j) 
             !   if(itest.lt.10)then
             !      if(k.eq.20) print*,'lin tend',iaq,k,i,j,spcaq_name(iaq),sc_t(k,i,j),sc_tc(k,i,j),sc_tr(k,i,j)
             !      if(k.eq.20) print*,'lin rap',iaq,k,i,j,spcaq_name(iaq),sc_p(k,i,j),sc_pc(k,i,j),sc_pr(k,i,j)
             !      itest=itest+1
             !   endif
             else
                ! cumulative tendencies assuming Henry equilibrium                      
                sc_t(k,i,j)=sc_t(k,i,j)+acgash/dtlt
                voltot=volc(k,i,j)+volr(k,i,j)
                accloudh=(cliqh*volc(k,i,j)/voltot)-sc_pc(k,i,j)
                acrainh=(cliqh*volr(k,i,j)/voltot)-sc_pr(k,i,j)
                sc_tc(k,i,j)=sc_tc(k,i,j)+accloudh/dtlt
                sc_tr(k,i,j)=sc_tr(k,i,j)+acrainh/dtlt
               ! if((k.eq.14).and.(i.eq.13).and.(j.eq.26)) print*,'dum henry',iaq,k,i,j,sc_t(k,i,j),sc_tc(k,i,j),sc_tr(k,i,j) 
               ! if(itest.lt.10)then
               !    if(k.eq.20)print*,'henry tend',iaq,k,i,j,spcaq_name(iaq),sc_t(k,i,j),sc_tc(k,i,j),sc_tr(k,i,j)
               !    if(k.eq.20)print*,'henry rap',iaq,k,i,j,spcaq_name(iaq),sc_p(k,i,j),sc_pc(k,i,j),sc_pr(k,i,j)
               !    itest=itest+1
               ! endif
             end if
          ENDDO
       ENDDO
    ENDDO

  end subroutine trans_tend
  !    
  !===================================================================
  SUBROUTINE trans_para(m1,m2,m3,ia,iz,ja,jz &
       ,theta,pi0,pp,rrp,rcp     &
       ,temp,volc,volr)

    !====================================================================
    !
    USE rconstants , ONLY:rgas,cp,cpor,p00

    IMPLICIT NONE
    INTEGER m1,m2,m3,ia,iz,ja,jz
    REAL pred, press
    REAL, dimension(m1,m2,m3) :: theta, pi0, pp, rcp, rrp, temp, volc, volr
    INTEGER i,j,k
    INTEGER itest
    !
    !itest=0
    !
    DO j = ja , jz
       DO i = ia , iz
          DO k = 1 , m1
             ! temperature (K) : temp
             pred=(pp(k,i,j)+pi0(k,i,j))/cp
             temp(k,i,j)=theta(k,i,j)*pred
             !
             ! pressure (Pa) : press
             press=p00*pred**cpor
             !
! volume of cloud and rain water in cm^3(water)/cm^-3(air)
           
          volc(k,i,j)=0.
          if(rcp(k,i,j)>1.e-9)then
            volc(k,i,j)=1.e-3*rcp(k,i,j)*press/(rgas*temp(k,i,j))
          endif
          volr(k,i,j)=0.
          if(rrp(k,i,j)>1.e-9)then
            volr(k,i,j)=1.e-3*rrp(k,i,j)*press/(rgas*temp(k,i,j))
          endif
             !
             !             
             !         if(volc(k,i,j).gt.0..and.volr(k,i,j).gt.0..and.itest.lt.5) then
             !             print*,'k,temp,press,rcp,volc,rrp,volr ',k,temp(k,i,j),press,rcp(k,i,j),volc(k,i,j) &
             !                    , rrp(k,i,j),volr(k,i,j)
             !             itest=itest+1
             !         endif
             !
          ENDDO
       ENDDO
    ENDDO
    !    
  end subroutine trans_para
  !===============================================================


END MODULE mod_chem_trans_gasaq
