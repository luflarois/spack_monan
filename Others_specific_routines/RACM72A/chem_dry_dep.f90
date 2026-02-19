   
MODULE module_dry_dep
!
      use chem1_list, only : nspecies_chem =>nspecies
      use aer1_list, only:  nspecies_aer=> nspecies,nmodes
!      
      integer, parameter :: dep_seasons = 5 ,   nlu = 25,                   &
        nseason = 1, nseasons = 2
!
! following currently hardwired to leaf-3 classes
      integer, parameter :: isice=2,iswater=1
!
      integer :: ixxxlu(nlu)
      real ::  kpart(nlu),&
               rac(nlu,dep_seasons), rclo(nlu,dep_seasons), rcls(nlu,dep_seasons), &
              rgso(nlu,dep_seasons), rgss(nlu,dep_seasons),                        &
                ri(nlu,dep_seasons),  rlu(nlu,dep_seasons)

      real, dimension(nspecies_chem) :: dratio

!- aerosol section
      integer, parameter :: ntotal=nspecies_aer*nmodes
      integer :: nspecies
      integer :: naer_a,naer_z
      integer, dimension(ntotal)   :: ind_aer
      integer, dimension(ntotal)   :: ind_mode
 
      logical           :: aer_alloc = .false.

      type sedim_type
      
            real,pointer,dimension  (:,:,:)  :: v_sed_part   
            real,pointer,dimension  (:,:,:)  :: r_lsl_part   
            real,pointer,dimension  (:,:  )  :: v_dep_part   
       
      end type
      
      type(sedim_type), ALLOCATABLE,DIMENSION(:,:)   :: dd_sedim
      
      integer :: NAER_TRANSPORTED


   CONTAINS
!========================================================================
subroutine drydep_driver(m1,m2,m3,ia,iz,ja,jz)

 
    USE mem_grid         ,  ONLY :  npatch,naddsc
    USE mem_grid         ,  ONLY :  MAXGRDS,nzpmax,dzt,zt,ngrids  &
                                   ,ngrid,nzg,nzs,jdim,dtlt,grid_g
    USE mem_basic        ,  ONLY :  basic_g
    USE mem_turb         ,  ONLY :  turb_g
    USE mem_micro        ,  ONLY :  micro_g
    USE mem_grid         ,  ONLY :  grid_g
    USE mem_scratch      ,  ONLY :  scratch
    USE mem_leaf         ,  ONLY :  leaf_g 
    USE extras           ,  ONLY :  extra2d 
    USE mem_radiate      ,  ONLY :  radiate_g 
    use node_mod         ,  ONLY :  mynum
    use mem_chem1, only: CHEMISTRY
    use chem1_list, nspecies_chem=> nspecies
    use aer1_list, nspecies_aer=> nspecies
    use mem_aer1, only: AEROSOL

  IMPLICIT NONE
  INTEGER,INTENT(IN) :: ia,iz,ja,jz,m1,m2,m3

!- no tracers, aerosols or chemistry, go back
if(CHEMISTRY < 0) return

      
    call dry_dep(ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,jdim,mynum,dtlt &
	    	,leaf_g(ngrid)%r_aer            &
            	,basic_g(ngrid)%theta	        &
            	,basic_g(ngrid)%rv	        &
            	,basic_g(ngrid)%pp	        &
            	,basic_g(ngrid)%dn0	        &
            	,basic_g(ngrid)%pi0	        &
            	,basic_g(ngrid)%up	        &
            	,basic_g(ngrid)%vp	        &
            	,turb_g(ngrid)%tkep	        &
            	,turb_g(ngrid)%sflux_t          &
	        ,turb_g(ngrid)%sflux_r          &
	        ,turb_g(ngrid)%sflux_u          &
	        ,turb_g(ngrid)%sflux_v          &
            	,radiate_g(ngrid)%rshort        &
		,micro_g(ngrid)%rcp	        &
            	,grid_g(ngrid)%rtgt	        &
            	,scratch%vt2da  	        &
            	,scratch%vt2db  	        &
            	,scratch%vt2dc  	        &
            	,scratch%vt2dd  	        &
            	,scratch%vt2de  	        &
            	,scratch%vt2df  	        &
		,scratch%vt3da                  &
		,scratch%vt3db                  &
                ,scratch%vt3dc                  &
		,scratch%vt3dd                  &
		,scratch%vt3de                  &
	    	,leaf_g(ngrid)%ustar	        &
	    	,leaf_g(ngrid)%tstar	        &
	    	,leaf_g(ngrid)%patch_area       &
	    	,leaf_g(ngrid)%leaf_class       &
	    	,leaf_g(ngrid)%patch_rough      &
	    	,MAXGRDS,ngrids		        &
	    	,dzt,zt,nzpmax)
	    
return
end subroutine drydep_driver
!========================================================================

subroutine dry_dep(ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,jdim,mynum,dt &
                  ,r_aer                                     &
		  ,theta,rv,pp,dn0,pi0,up,vp                 &
		  ,tke,sflux_t,sflux_r,sflux_u,sflux_v       &  
		  ,rshort,rcp,rtgt                           &
		  ,temps,prss,dens,vels,rvs,Zi               &
                  ,vt3da,vt3db,vt3dc,temp3d,air_dens3d       &
		  ,ustar,tstar,patch_area,veg,Z0m	     &
	          ,maxgrds,ngrids,dzt,zt,nzpmax  )

    USE rconstants         ,  ONLY :  cpi,cpor,p00,g,vonk
    USE mem_scalar         ,  ONLY :  scalar_g
    use mem_cuparm         ,  ONLY :  nnqparm,cuparm_g
    USE mem_micro          ,  ONLY :  micro_g
    use micphys            ,  only :  level
    use chem1_list, only:   nspecies_chem  =>nspecies    &
                          ,spc_alloc_chem  =>spc_alloc &
                          , spc_name_chem  =>spc_name    &
                          ,ddp,transport,off 
    use mem_chem1, only: chem1_g,CHEMISTRY
    
    use aer1_list, only: nspecies_aer   =>nspecies &
                        ,spc_alloc_aer  =>spc_alloc 
    use mem_aer1 , only: aer1_g ,AEROSOL
!


    IMPLICIT NONE
    INTEGER, INTENT(IN)  ::  ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,jdim,mynum
    INTEGER, INTENT(IN)  ::  maxgrds,nzpmax,ngrids
    REAL, INTENT (IN)    ::  dt
    REAL, PARAMETER      ::  ubmin = 0.25
   
    
    REAL, INTENT (IN)    ::  r_aer(m2,m3,npatch)
    REAL, INTENT (IN)    :: theta(m1,m2,m3)
    REAL, INTENT (IN)    ::    rv(m1,m2,m3)
    REAL, INTENT (IN)    ::    pp(m1,m2,m3)
    REAL, INTENT (IN)    ::   dn0(m1,m2,m3)
    REAL, INTENT (IN)    ::   pi0(m1,m2,m3)
    REAL, INTENT (IN)    ::    up(m1,m2,m3)
    REAL, INTENT (IN)    ::    vp(m1,m2,m3)
    REAL, INTENT (IN)    ::   tke(m1,m2,m3)
    REAL, INTENT (IN)    ::   rcp(m1,m2,m3)
    REAL, INTENT (INOUT) ::   vt3da(m1,m2,m3)
    REAL, INTENT (INOUT) ::   vt3db(m1,m2,m3)
    REAL, INTENT (INOUT) ::   vt3dc(m1,m2,m3)
    REAL, INTENT (INOUT) ::   temp3d(m1,m2,m3)
    REAL, INTENT (INOUT) ::   air_dens3d(m1,m2,m3)

    REAL, INTENT (IN)    ::   sflux_t(m2,m3)
    REAL, INTENT (IN)    ::   sflux_r(m2,m3)
    REAL, INTENT (IN)    ::   sflux_u(m2,m3)
    REAL, INTENT (IN)    ::   sflux_v(m2,m3)
    REAL, INTENT (IN)    ::      rtgt(m2,m3)
    REAL, INTENT (INOUT) ::     temps(m2,m3)   
    REAL, INTENT (INOUT) ::      prss(m2,m3)    
    REAL, INTENT (INOUT) ::      dens(m2,m3)    
    REAL, INTENT (INOUT) ::      vels(m2,m3)    
    REAL, INTENT (INOUT) ::       rvs(m2,m3)     
    REAL, INTENT (IN)    ::        Zi(m2,m3)
    REAL, INTENT (IN)    ::    rshort(m2,m3)
    
    REAL, INTENT (IN)    ::       ustar(m2,m3,npatch)
    REAL, INTENT (IN)    ::       tstar(m2,m3,npatch)
    REAL, INTENT (IN)    ::  patch_area(m2,m3,npatch)
    REAL, INTENT (IN)    ::         veg(m2,m3,npatch)
    REAL, INTENT (IN)    ::         Z0m(m2,m3,npatch)
    
    REAL, INTENT (IN)    :: dzt(nzpmax)
    REAL, INTENT (IN)    ::  zt(nzpmax)

    integer :: i,j,ipatch,idry_part,ispc,k
    real ups,vps,pis,sbf

    REAL   ::  v_dep_gas(nspecies_chem,m2,m3)
    real   :: check_rain(m2,m3              )
    REAL   ::	    rmol(m2,m3              )!1./Monin-Obukhob
    REAL   ::    rhchem(m2,m3)     
    real,external :: rs


    v_dep_gas = 0.0   
    check_rain = 0.
    rmol = 0.
!-aux variables

do j = ja,jz
 do i = ia,iz
    rvs  (i,j) = rv(2,i,j)
    pis        = ( pp(1,i,j) + pp(2,i,j) + pi0(1,i,j) + pi0(2,i,j) ) * .5 * cpi  
    prss (i,j) = pis ** cpor * p00                                               
    dens (i,j) = ( dn0(1,i,j) + dn0(2,i,j) ) * .5
    temps(i,j) = theta(2,i,j) * pis        ! temps=theta*Exner/CP
    rhchem(i,j)=100.*min(1.,max(0.05,rvs(i,j)/rs(prss(i,j),temps(i,j))))
    ups        = (up(2,i,j) + up(2,i-1,j   )) * .5
    vps        = (vp(2,i,j) + vp(2,i,j-jdim)) * .5
    vels (i,j) = max(ubmin,sqrt(ups** 2 + vps** 2))
    !- compute surface buoyancy flux [sbf] and 1 / Monin Obukhov height [1/zl]	    
    sbf = (g/theta(2,i,j) * (1. + .61 * rv(2,i,j))) &
        * ( sflux_t(i,j) * (1. + .61 * rv(2,i,j)) &
	+   sflux_r(i,j) * .61 * theta(2,i,j)   )
    !-  reciprocal of the Monin-Obukhov length (1/m)
    rmol (i,j) = - vonk*sbf / max(0.1,sqrt(sqrt(sflux_u(i,j)**2 + sflux_v(i,j)**2)))
    if(abs(rmol (i,j)) < 1.e-7) rmol (i,j) = 0. 
!    print*,'rmol=',rmol(i,j)
 enddo
enddo
if (nnqparm(ngrid) > 0) then
   check_rain(ia:iz,ja:jz)  = cuparm_g(ngrid)%conprr(ia:iz,ja:jz)    
endif
if (level >= 3) then
   check_rain(ia:iz,ja:jz)  = check_rain(ia:iz,ja:jz)  +  micro_g(ngrid)%pcpg(ia:iz,ja:jz)
endif

call define_PBL_height(m1,m2,m3,npatch,ia,iz,ja,jz,zt     &
		      ,tke,sflux_t,rcp,rtgt &
		      ,zi)


!---  dry deposition for gases
call dry_dep_gases(m1,m2,m3,nspecies_chem,npatch,ia,iz,ja,jz  &
                      ,v_dep_gas,r_aer  		           &
		      ,prss,temps,dens,vels,rvs,rcp,Zi		   &
                      ,ustar,tstar,patch_area,veg,Z0m		   &
		      ,rshort,rtgt,dzt,check_rain,rmol,rhchem)
		      
!-apply dry deposition on the tracers concentration and get the deposited mass on surface
do ispc=1,nspecies_chem
    if(spc_alloc_chem(ddp,ispc) == off .or. spc_alloc_chem(transport,ispc) == off) cycle
   
       call apply_drydep(m1,m2,m3,ia,iz,ja,jz		&
        		,v_dep_gas(ispc,1:m2,1:m3)	& ! don't change this
        	        ,chem1_g  (ispc,ngrid)%sc_t     &! tendency array
        		,chem1_g  (ispc,ngrid)%sc_dd    &! deposited mass
        		,chem1_g  (ispc,ngrid)%sc_p     &! mixing ratio
        		,dens,rtgt,dzt,dt)       
enddo

!--- aerosol dry deposition and sedimentation 

If(AEROSOL > 0) then
       
 do j = ja,jz
  do i = ia,iz
   do k = 1, m1-1
 	 temp3d(k,i,j) = 0.5*( theta(k+1,i,j) * ( pp(k+1,i,j)+pi0(k+1,i,j) ) * cpi + &
   			       theta(k  ,i,j) * ( pp(k  ,i,j)+pi0(k  ,i,j) ) * cpi   )
     air_dens3d(k,i,j) = 0.5 *(  dn0(k+1,i,j) +  dn0(k  ,i,j) )
    enddo    
 	 temp3d(m1,i,j) =	 temp3d(m1-1,i,j) 
   					      
     air_dens3d(m1,i,j) =    air_dens3d(m1-1,i,j)
   enddo
 enddo   
    !print*,'temp=',temp3d(1:2,int(iz/2),int(jz/2))
    
 !--  dry deposition/sedimentation for particles
  call dry_dep_sedim_particles(ngrids,ngrid,m1,m2,m3,nspecies_aer,npatch,ia,iz,ja,jz  &
  			    ,r_aer,temp3d,air_dens3d,temps,dens,vels,rvs,Zi	      &
  			    ,ustar,tstar,patch_area,veg,Z0m		 )
    
  !- loop over the transported aerosols particles
  do ispc=naer_a,naer_z
  
 
   if(spc_alloc_aer(ddp      ,ind_mode(ispc),ind_aer(ispc)) == off .or. &
      spc_alloc_aer(transport,ind_mode(ispc),ind_aer(ispc)) == off) cycle

 
      !print*,'dd=',dd_sedim(ispc,ngrid)%v_dep_part(int(iz/2),int(jz/2))*100.
       
      call apply_drydep(m1,m2,m3,ia,iz,ja,jz			      &
     		,dd_sedim(ispc,ngrid)%v_dep_part           	      &
                ,aer1_g  (ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t    &! tendency array
      	        ,aer1_g  (ind_mode(ispc),ind_aer(ispc),ngrid)%sc_dd   &! deposited mass
                ,aer1_g  (ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p    &! mixing ratio
                ,dens,rtgt,dzt,dt)
       !call apply_drydep(m1,m2,m3,ia,iz,ja,jz  			  &
       !	   ,dd_sedim(ispc,ngrid)%v_dep_part(1:,1:)		 &
       ! 	   ,aer1_g  (ind_mode(ispc),ind_aer(ispc),ngrid)%sc_t (1:    )&! tendency array
       !	   ,aer1_g  (ind_mode(ispc),ind_aer(ispc),ngrid)%sc_dd(1:,1:  )&! deposited mass
       !	   ,aer1_g  (ind_mode(ispc),ind_aer(ispc),ngrid)%sc_p (1:,1:,1:)&! mixing ratio
       !	   ,dens,rtgt,dzt,dt)
       
  enddo    
endif  ! endif of aerosol

end subroutine dry_dep
!========================================================================

subroutine dry_dep_gases(m1,m2,m3,nspecies_chem,npatch,ia,iz,ja,jz &
                       ,V_dep,r_aer,prss,temps,dens,vels,rvs,rcp,Zi&
		       ,ustar,tstar,patch_area,veg,Z0m,rshort,rtgt,dzt&
		       ,check_rain,rmol,rhchem)
use chem1_list, only : O3, SULF                 &
                , spc_alloc_chem  =>spc_alloc   &
                , spc_name_chem   =>spc_name    &
                , transport,ON 
		
use mem_grid,   only: time,iyear1,imonth1,idate1,itime1
implicit none
integer :: m1,m2,m3,nspecies_chem,npatch,ia,iz,ja,jz,i,j,ipatch
!
real, dimension(m1,m2,m3)     :: rcp
real, dimension(m2,m3)        :: prss,temps,dens,vels,rvs,Zi,rshort,rtgt,check_rain,rmol,rhchem
real, dimension(m2,m3,npatch) :: ustar,tstar,patch_area,veg,Z0m,r_aer
real, dimension(m1) :: dzt  

real, dimension(nspecies_chem,m2,m3)     :: V_dep
integer,external :: julday


!- local var
real,  dimension(nspecies_chem) :: vgs0d  !,srfres
integer iveg,iseason,jday
logical :: highnh3, rainflag, vegflag, wetflag
real dvpart,dvfog,clwchem,z1,ustarw,vgpart
!real , dimension(m2,m3) :: aer_res
!real   ::  r_lsl(m2,m3,npatch,nspecies_chem),rcx(npatch,nspecies_chem)
real   ::  r_lsl(nspecies_chem),rcx(nspecies_chem)

!   Set the reference height (10.0 m)
real, parameter ::  zr = 10.0


  jday=julday(imonth1,idate1,iyear1)
  iseason = 1
  if( jday.lt.90 .or. jday.gt.270)  iseason=2

do j = ja,jz
  do i = ia,iz
      ustarw = 0.
!
!     Set logical default values
      rainflag = .FALSE.
      if (check_rain(i,j) > 0.  ) rainflag = .true.
      
      wetflag  = .FALSE.
      if (rhchem(i,j) >= 95.) wetflag  = .true.

      clwchem = rcp(2,i,j)
      
!-not using
      highnh3 = .FALSE.
!      if(chem(i,kts,j,p_nh3).gt.2.*chem(i,kts,j,p_so2))highnh3 = .true.
!----


      !zntt = znt(i,j)  ! Surface roughness height (m)
      
      z1 =   rtgt(i,j)/dzt(2)! height of the fisrt cell

      do ipatch = 1,npatch
         if (patch_area(i,j,ipatch) .ge. .009) then 
            		   
	    if(ipatch == 1) then
	       iveg = 1           ! rc routine needs water as 1
	    else
	       iveg = nint(veg(i,j,ipatch))     
	    endif
	     
	    ustarw=ustarw+ustar(i,j,ipatch)*patch_area(i,j,ipatch)	
	    !
            !- calculates surface resistance (rcx) 
            !
	    call rc(rcx,temps(i,j),rshort(i,j),rhchem(i,j),iveg,iseason      &
                   ,nspecies_chem,wetflag,rainflag,highnh3)
	    !- get R_b 
            call get_rb(nspecies_chem,r_lsl,temps(i,j) &
	                                   ,prss(i,j),ustar(i,j,ipatch))
            !
            !- deposition velocity for each patch 
            vgs0d(1:nspecies_chem) = 1./( r_aer(i,j,ipatch     )  + &
	                                  r_lsl(1:nspecies_chem)  + &
                                          rcx  (1:nspecies_chem)    )
	    !			  
            !- special treatment for SULF
	    if( spc_alloc_chem(transport,SULF) == ON) then
             CALL deppart(rmol(i,j),ustar(i,j,ipatch),rhchem(i,j),clwchem,iveg,dvpart,dvfog)
             CALL depvel(nspecies_chem,rmol(i,j),zr,Z0m(i,j,ipatch)   &
	               ,ustar(i,j,ipatch),vgpart)
             vgs0d(sulf)=1.0/((1.0/vgpart)+(1.0/dvpart))
	    endif
	    !- end SULF treatment
	    !
	    !
            !- effective deposition velocity for all patches
            V_dep(1:nspecies_chem,i,j) = V_dep(1:nspecies_chem,i,j) + &
	                                 patch_area(i,j,ipatch)*vgs0d(1:nspecies_chem)

         endif
     enddo                     
!	print*,'O3-1 cm/s=',V_dep(O3,i,j)*100; call flush(6)
     !- get the effective deposition velocity
     call cellvg(V_dep(1:nspecies_chem,i,j),ustarw,z1,zr,rmol(i,j),nspecies_chem)
!	print*,'O3-2 cm/s=',V_dep(O3,i,j)*100; call flush(6)

  enddo                     
enddo                     
return
end subroutine dry_dep_gases
!========================================================================

subroutine get_rb(nspecies_chem,r_lsl,temps,prss,ustar)
use chem1_list, only : dvj

implicit none
integer nspecies_chem
real temps,ustar,prss
real r_lsl(nspecies_chem)
!local var
real, parameter :: i23=2./3., STDTEMP = 273.15, STDATMPA = 101325.0 
integer ispc
real dvj_tmp,sc,scpr23,eta


!-  kinematic viscosity of air at surface (cm^2/s)
eta = -1.1555E-10*temps**3 + 9.5728E-07*temps**2 + 3.7604E-04*temps - 3.4484E-02  

do ispc=1,nspecies_chem

!         dvj(ispc) = dvj(ispc)*(293.15/298.15)**1.75
!         dvj_tmp = dvj(ispc)*(temps/(temps+5.))**1.75
          dvj_tmp = dvj(ispc)* (( temps / STDTEMP ) ** 1.75)&
                   * (STDATMPA / prss ) ! cm^2/s
!          dratio(ispc) = 0.242/dvj(ispc) 

!         sc = 0.15/dvj_tmp ! Schmidt Number at 20°C   
          sc =  eta/dvj_tmp ! Schmidt Number at any temps 

	  scpr23 = (sc/0.72)**(i23) ! (Schmidt # / Prandtl #)**
          r_lsl(ispc) = 5.*scpr23/ustar
enddo
end subroutine get_rb

! **********************************************************************
! **************************  SUBROUTINE RC  ***************************
! **********************************************************************
      SUBROUTINE rc(rcx,t,rad,rh,iland,iseason,nspecies_chem,             &
                    wetflag,rainflag,highnh3)
!     THIS SUBROUTINE CALCULATES SURFACE RESISTENCES ACCORDING
!     TO THE MODEL OF
!     M. L. WESELY,
!     ATMOSPHERIC ENVIRONMENT 23 (1989), 1293-1304
!     WITH SOME ADDITIONS ACCORDING TO
!     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
!     WRITTEN BY  WINFRIED SEIDL, APRIL 1997
!     MODYFIED BY WINFRIED SEIDL, MARCH 2000
!                    FOR MM5 VERSION 3
!----------------------------------------------------------------------
      use chem1_list, spc_alloc_chem  =>spc_alloc   &
                     , spc_name_chem   =>spc_name  
     implicit none

! .. Scalar Arguments ..
        REAL :: rad, rh, t, temp_i
        INTEGER :: iland, iseason, nspecies_chem
        LOGICAL :: highnh3, rainflag, wetflag
        real :: rcx(nspecies_chem)
! ..
! .. Array Arguments ..
        INTEGER :: iprt
! ..
! .. Local Scalars ..
        REAL :: rclx, rdc, resice, rgsx, rluo1, rluo2, rlux, rmx, rs, rsmx, &
          tc, rdtheta, z, hplus, corrh
        INTEGER :: n
! ..
! .. Local Arrays ..
        REAL :: hstary(nspecies_chem)
! ..
! .. Intrinsic Functions ..
        INTRINSIC exp
!srf- declaring NH3, because it is not in the chem1-list for RACM
        integer :: NH3 ! see return comand before the calculation of RC for NH3
!srf---
! ..
        DO n = 1, nspecies_chem
          rcx(n) = 1.
        END DO

        tc = t - 273.15
        rdtheta = 0.
	temp_i = 1./t-1./298.

        z = 200./(rad+0.1)

!!!  HARDWIRE VALUES FOR TESTING
!       z=0.4727409
!       tc=22.76083
!       t=tc+273.15
!       rad = 412.8426
!       rainflag=.false.
!       wetflag=.false.

        IF ((tc<=0.) .OR. (tc>=40.)) THEN
          rs = 9999.
        ELSE
          rs = ri(iland,iseason)*(1+z*z)*(400./(tc*(40.-tc)))
        END IF
        rdc = 100*(1.+1000./(rad+10))/(1+1000.*rdtheta)
        rluo1 = 1./(1./3000.+1./3./rlu(iland,iseason))
        rluo2 = 1./(1./1000.+1./3./rlu(iland,iseason))
        resice = 1000.*exp(-tc-4.)

! change MP 11/12/07 taking into account the acid dissociation constant
! and assuming pH=7 in the vegetation
!srf-opt  hstary(n) = hstar(n)*exp(dhr(n)*(1./t-1./298.))
!                    *(1+ak0(n)*exp(dak(n)*(1./t-1/298.))/hplus)
!          pH=7  hplus=10**(-pH) 
          hplus=1.E-7
        DO n = 1, nspecies_chem
!srf      IF (hstar(n)==0.) GO TO 10
          IF (hstar(n)==0.) cycle
!srf
          corrh=1+ak0(n)*exp(dak(n)*(temp_i))/hplus
          hstary(n) = hstar(n)*exp(dhr(n)*(temp_i))*corrh
!
!         if ((n == 8).or.(n == 52))then
!           print*, n, hstary(n),t,temp_i,ak0(n), dak(n),hplus,corrh &
!                , hstar(n), dhr(n)
!         endif
! end change MP
          rmx = 1./(hstary(n)/3000.+100.*f0(n))
          rsmx = rs*dratio(n) + rmx
          rclx = 1./(hstary(n)/1.E+5/rcls(iland,iseason)+f0(n)/rclo(iland, &
            iseason)) + resice
          rgsx = 1./(hstary(n)/1.E+5/rgss(iland,iseason)+f0(n)/rgso(iland, &
            iseason)) + resice
          rlux = rlu(iland,iseason)/(1.E-5*hstary(n)+f0(n)) + resice
          IF (wetflag) THEN
            rlux = 1./(1./3./rlu(iland,iseason)+1.E-7*hstary(n)+f0(n)/rluo1)
          END IF
        IF (rainflag) THEN
            rlux = 1./(1./3./rlu(iland,iseason)+1.E-7*hstary(n)+f0(n)/rluo2)
          END IF
          rcx(n) = 1./(1./rsmx+1./rlux+1./(rdc+rclx)+1./(rac(iland, &
            iseason)+rgsx))
          IF (rcx(n)<1.) rcx(n) = 1.
10      END DO
!

       if( spc_alloc_chem(transport,O3) == ON) then

!--------------------------------- O3 
!        SPECIAL TREATMENT FOR OZONE
!srf       hstary(O3) = hstar(O3)*exp(dhr(O3)*(1./t-1./298.))
           hstary(O3) = hstar(O3)*exp(dhr(O3)*(temp_i))
           rmx = 1./(hstary(O3)/3000.+100.*f0(O3))
           rsmx = rs*dratio(O3) + rmx
           rlux = rlu(iland,iseason)/(1.E-5*hstary(O3)+f0(O3)) + resice
           rclx = rclo(iland,iseason) + resice
           rgsx = rgso(iland,iseason) + resice
           IF (wetflag) rlux = rluo1
           IF (rainflag) rlux = rluo2
           rcx(O3) = 1./(1./rsmx+1./rlux+1./(rdc+rclx)+1./(rac(iland, &
             iseason)+rgsx))
           rcx(O3)=max(rcx(O3),1.) !- old: IF (rcx(O3)<1.) rcx(O3) = 1.
      endif


      if( spc_alloc_chem(transport,SO2) == ON) then

!    	      SPECIAL TREATMENT FOR SO2 (Wesely)
!    		HSTARY(SO2)=HSTAR(SO2)*EXP(DHR(SO2)*(1./T-1./298.))
!    		RMX=1./(HSTARY(SO2)/3000.+100.*F0(SO2))
!    		RSMX=RS*DRATIO(SO2)+RMX
!    		RLUX=RLU(ILAND,ISEASON)/(1.E-5*HSTARY(SO2)+F0(SO2))
!    	     &       +RESICE
!    		RCLX=RCLS(ILAND,ISEASON)+RESICE
!    		RGSX=RGSS(ILAND,ISEASON)+RESICE
!    		IF ((wetflag).OR.(RAINFLAG)) THEN
!    		  IF (ILAND.EQ.1) THEN
!    		    RLUX=50.
!    		  ELSE
!    		    RLUX=100.
!    		  END IF
!    		END IF
!    		RCX(SO2)=1./(1./RSMX+1./RLUX+1./(RDC+RCLX)
!    	     &  	      +1./(RAC(ILAND,ISEASON)+RGSX))
!    		IF (RCX(SO2).LT.1.) RCX(SO2)=1.
!
!
!----------------------------------------- SO2 
!    	      SO2 according to Erisman et al. 1994
!    		R_STOM
     		rsmx = rs*dratio(SO2)
!    		R_EXT
     		IF (tc>(-1.)) THEN
     		  IF (rh<81.3) THEN
     		    rlux = 25000.*exp(-0.0693*rh)
     		  ELSE
     		    rlux = 0.58E12*exp(-0.278*rh)
     		  END IF
     		END IF
     		IF (((wetflag) .OR. (rainflag)) .AND. (tc>(-1.))) THEN
     		  rlux = 1.
     		END IF
     		IF ((tc>=(-5.)) .AND. (tc<=(-1.))) THEN
     		  rlux = 200.
     		END IF
     		IF (tc<(-5.)) THEN
     		  rlux = 500.
     		END IF
!    		INSTEAD OF R_INC R_CL and R_DC of Wesely are used
     		rclx = rcls(iland,iseason)
!    		DRY SURFACE
     		rgsx = 1000.
!    		WET SURFACE
     		IF ((wetflag) .OR. (rainflag)) THEN
     		  IF (highnh3) THEN
     		    rgsx = 0.
     		  ELSE
     		    rgsx = 500.
     		  END IF
     		END IF
!    		WATER
     		IF (iland==iswater) THEN
     		  rgsx = 0.
     		END IF
!    		SNOW
     		IF ((iseason==4) .OR. (iland==isice)) THEN
     		  IF (tc>2.) THEN
     		    rgsx = 0.
     		  END IF
     		  IF ((tc>=(-1.)) .AND. (tc<=2.)) THEN
     		    rgsx = 70.*(2.-tc)
     		  END IF
     		  IF (tc<(-1.)) THEN
     		    rgsx = 500.
     		  END IF
     		END IF
!    		TOTAL SURFACE RESISTENCE
     		IF ((iseason/=4) .AND. (ixxxlu(iland)/=1) .AND. (iland/=iswater) .AND. &
     		    (iland/=isice)) THEN
     		  rcx(SO2) = 1./(1./rsmx+1./rlux+1./(rclx+rdc+rgsx))
     		ELSE
     		  rcx(SO2) = rgsx
     		END IF
     		IF (rcx(SO2)<1.) rcx(SO2) = 1.
!
      endif


!------------ NOT IN USE BELOW -----------------------------------------------
!srf- dry dep for NH3 NOT in use :
RETURN 
!srf ---

!     NH3 according to Erisman et al. 1994
!       R_STOM
        rsmx = rs*dratio(NH3)
!       GRASSLAND (PASTURE DURING GRAZING)
        IF (ixxxlu(iland)==3) THEN
          IF (iseason==1) THEN
!           SUMMER
            rcx(NH3) = 1000.
          END IF
          IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
!           WINTER, NO SNOW
            IF (tc>-1.) THEN
              IF (rad/=0.) THEN
                rcx(NH3) = 50.
              ELSE
                rcx(NH3) = 100.
              END IF
              IF ((wetflag) .OR. (rainflag)) THEN
                rcx(NH3) = 20.
              END IF
            END IF
            IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
              rcx(NH3) = 200.
            END IF
            IF (tc<(-5.)) THEN
              rcx(NH3) = 500.
            END IF
          END IF
        END IF
!       AGRICULTURAL LAND (CROPS AND UNGRAZED PASTURE)
        IF (ixxxlu(iland)==2) THEN
          IF (iseason==1) THEN
!           SUMMER
            IF (rad/=0.) THEN
              rcx(NH3) = rsmx
            ELSE
              rcx(NH3) = 200.
            END IF
            IF ((wetflag) .OR. (rainflag)) THEN
              rcx(NH3) = 50.
            END IF
          END IF
          IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
!           WINTER, NO SNOW
            IF (tc>-1.) THEN
              IF (rad/=0.) THEN
                rcx(NH3) = rsmx
              ELSE
                rcx(NH3) = 300.
              END IF
              IF ((wetflag) .OR. (rainflag)) THEN
                rcx(NH3) = 100.
              END IF
            END IF
            IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
              rcx(NH3) = 200.
            END IF
            IF (tc<(-5.)) THEN
              rcx(NH3) = 500.
            END IF
          END IF
        END IF
!       SEMI-NATURAL ECOSYSTEMS AND FORESTS
        IF ((ixxxlu(iland)==4) .OR. (ixxxlu(iland)==5) .OR. (ixxxlu( &
            iland)==6)) THEN
          IF (rad/=0.) THEN
            rcx(NH3) = 500.
          ELSE
            rcx(NH3) = 1000.
          END IF
          IF ((wetflag) .OR. (rainflag)) THEN
            IF (highnh3) THEN
              rcx(NH3) = 100.
            ELSE
              rcx(NH3) = 0.
            END IF
          END IF
          IF ((iseason==2) .OR. (iseason==3) .OR. (iseason==5)) THEN
!           WINTER, NO SNOW
            IF ((tc>=(-5.)) .AND. (tc<=-1.)) THEN
              rcx(NH3) = 200.
            END IF
            IF (tc<(-5.)) THEN
              rcx(NH3) = 500.
            END IF
          END IF
        END IF
!       WATER
        IF (iland==iswater) THEN
          rcx(NH3) = 0.
        END IF
!       URBAN AND DESERT (SOIL SURFACES)
        IF (ixxxlu(iland)==1) THEN
          IF ( .NOT. wetflag) THEN
            rcx(NH3) = 50.
          ELSE
            rcx(NH3) = 0.
          END IF
        END IF
!       SNOW COVERED SURFACES OR PERMANENT ICE
        IF ((iseason==4) .OR. (iland==isice)) THEN
          IF (tc>2.) THEN
            rcx(NH3) = 0.
          END IF
          IF ((tc>=(-1.)) .AND. (tc<=2.)) THEN
            rcx(NH3) = 70.*(2.-tc)
          END IF
          IF (tc<(-1.)) THEN
            rcx(NH3) = 500.
          END IF
        END IF
        IF (rcx(NH3)<1.) rcx(NH3) = 1.

      END SUBROUTINE rc
! **********************************************************************
      SUBROUTINE cellvg(vgtemp,ustar,dz,zr,rmol,nspec)
!     THIS PROGRAM HAS BEEN DESIGNED TO CALCULATE THE CELL AVERAGE
!     DEPOSITION VELOCITY GIVEN THE VALUE OF VG AT SOME REFERENCE
!     HEIGHT ZR WHICH IS MUCH SMALLER THAN THE CELL HEIGHT DZ.
!       PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner    (February 1991)
!.....PROGRAM VARIABLES...
!     VgTemp   - DEPOSITION VELOCITY AT THE REFERENCE HEIGHT
!     USTAR    - FRICTION VELOCITY
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     DZ       - CELL HEIGHT
!     CELLVG   - CELL AVERAGE DEPOSITION VELOCITY
!     VK       - VON KARMAN CONSTANT
!     Local Variables
! .. Scalar Arguments ..
        REAL :: dz, rmol, ustar, zr
        INTEGER :: nspec
! ..
! .. Array Arguments ..
        REAL :: vgtemp(nspec)
! ..
! .. Local Scalars ..
        REAL :: a, fac, pdz, pzr,fx
        INTEGER :: nss
! ..
! .. Intrinsic Functions ..
        INTRINSIC alog, sqrt
! ..
!     Set the von Karman constant
       REAL, parameter :: vk = 0.4

!     DETERMINE THE STABILITY BASED ON THE CONDITIONS
!             1/L < 0 UNSTABLE
!             1/L = 0 NEUTRAL
!             1/L > 0 STABLE


       IF (rmol<0) THEN
         pdz = sqrt(1.0-9.0*dz*rmol)
         pzr = sqrt(1.0-9.0*zr*rmol)
         fac = ((pdz-1.0)/(pzr-1.0))*((pzr+1.0)/(pdz+1.0))
         a = 0.74*dz*alog(fac) + (0.164/rmol)*(pdz-pzr)
       ELSE IF (rmol==0) THEN
         a = 0.74*(dz*alog(dz/zr)-dz+zr)
       ELSE
         a = 0.74*(dz*alog(dz/zr)-dz+zr) + (2.35*rmol)*(dz-zr)**2
       END IF
      
       fx = a/(vk*ustar*(dz-zr))
      
!      CALCULATE THE DEPOSITION VELOCITIY
       DO nss = 1, nspec
           vgtemp(nss) = vgtemp(nss)/(1.0+vgtemp(nss)*fx)
       END DO
      
        RETURN
      END SUBROUTINE cellvg


!========================================================================

subroutine apply_drydep(m1,m2,m3,ia,iz,ja,jz,V_dep,sclt,M_dep,sclp &
                       ,dens,rtgt,dzt,dt)
implicit none
integer :: m1,m2,m3,ia,iz,ja,jz,i,j

real, dimension(m2,m3)    :: V_dep &! dry deposition velocity (m/s)
                            ,M_dep  ! accumulated mass on surface due dry dep 
				    ! process (kg m^-2)
real, dimension(m1,m2,m3) :: sclp & ! tracer mixing ratio (kg/kg)
                            ,sclt   ! tendency (kg[tracer] kg[air]^-1 s^-1)
real, dimension(m2,m3) :: rtgt,dens  
real, dimension(m1) :: dzt  
real dt,dz,tend_drydep			    

do j=ja,jz
   do i=ia,iz
!-1st vertical layer thickness
        dz = rtgt(i,j)/dzt(2) ! dzt=1/(z(k)-z(k-1))
!- tendency to dry deposition process	
	tend_drydep  = - V_dep(i,j)*sclp(2,i,j)/(dz)
!- update total tendency kg[tracer]/kg[air]/s
        sclt(2,i,j) = sclt(2,i,j) + tend_drydep 
!- accumulate the surface deposited mass of the tracer by this process
!	M_dep(i,j) = M_dep(i,j) + dt*tend_drydep*dens(i,j) ! kg[tracer] m^-2
	M_dep(i,j) = M_dep(i,j) - dt*tend_drydep*dens(i,j)*dz ! kg[tracer] m^-2
   enddo
enddo
return
end subroutine apply_drydep
!========================================================================

subroutine define_PBL_height(m1,m2,m3,npatch,ia,iz,ja,jz,zt,tke,shf,rcp,rtgt,Zi)


implicit none
integer :: m1,m2,m3,npatch,ia,iz,ja,jz,i,j,k
real, dimension(m1)    ::    zt
real, dimension(m2,m3) ::    shf,rtgt,Zi
real, dimension(m1,m2,m3) :: tke,rcp
real pblht
REAL,PARAMETER :: tkethrsh=0.001       !   tke threshold for PBL height in m2/s2
                                         !   tkmin    = 5.e-4   minimum TKE in RAMS 
REAL,PARAMETER :: rcmin=1.e-4          !   min liq water = 0.1 g/kg

   do j=ja,jz
      do i=ia,iz

         Zi(i,j) = 0.
!- convective layer
	 if(shf(i,j) >= 1.e-8) then 
      
          pblht=0.
          do k=2,m1-1
            pblht=zt(k)*rtgt(i,j) 
!                  if(i.ge.10.and.i.le.25.and.j.ge.13.and.j.le.25)
!     &               print*,'i,j,k,z,pbl=',i,j,k,ztn(k,ngrd),pblht
	    if( rcp(k,i,j) .gt. rcmin     )goto 10 ! dry convective layer
            if( tke(k,i,j) .le. tkethrsh  )goto 10 
          enddo
 10       continue

          Zi(i,j)=pblht
	 endif
	 !print*,'PBLh',Zi(i,j)

      enddo
   enddo
return
end subroutine define_PBL_height
!========================================================================
!========================================================================
!========================================================================
!========================================================================
      SUBROUTINE dep_init
      use chem1_list, nspecies_chem =>nspecies
      implicit none
! ..
! ..
! .. Local Scalars ..
        REAL :: sc
        INTEGER :: iland, iseason, l
        integer :: iprt
! ..
! .. Local Arrays ..
        REAL :: dat1(nlu,dep_seasons), dat2(nlu,dep_seasons),         &
                dat3(nlu,dep_seasons), dat4(nlu,dep_seasons),         &
                dat5(nlu,dep_seasons), dat6(nlu,dep_seasons),         &
                dat7(nlu,dep_seasons)
! ..
! .. Data Statements ..
!     RI for stomatal resistance
!      data ((ri(ILAND,ISEASON),ILAND=1,nlu),ISEASON=1,dep_seasons)/0.10E+11, &
        DATA ((dat1(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.60E+02, 0.60E+02, 0.60E+02, 0.60E+02, 0.70E+02, 0.12E+03, &
          0.12E+03, 0.12E+03, 0.12E+03, 0.70E+02, 0.13E+03, 0.70E+02, &
          0.13E+03, 0.10E+03, 0.10E+11, 0.80E+02, 0.10E+03, 0.10E+11, &
          0.80E+02, 0.10E+03, 0.10E+03, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.10E+11, 0.10E+11, &
          0.70E+02, 0.25E+03, 0.50E+03, 0.10E+11, 0.10E+11, 0.50E+03, &
          0.10E+11, 0.10E+11, 0.50E+03, 0.50E+03, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.10E+11, &
          0.10E+11, 0.70E+02, 0.25E+03, 0.50E+03, 0.10E+11, 0.10E+11, &
          0.50E+03, 0.10E+11, 0.10E+11, 0.50E+03, 0.50E+03, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.70E+02, 0.40E+03, 0.80E+03, 0.10E+11, &
          0.10E+11, 0.80E+03, 0.10E+11, 0.10E+11, 0.80E+03, 0.80E+03, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.12E+03, 0.12E+03, &
          0.12E+03, 0.12E+03, 0.14E+03, 0.24E+03, 0.24E+03, 0.24E+03, &
          0.12E+03, 0.14E+03, 0.25E+03, 0.70E+02, 0.25E+03, 0.19E+03, &
          0.10E+11, 0.16E+03, 0.19E+03, 0.10E+11, 0.16E+03, 0.19E+03, &
          0.19E+03, 0.10E+11, 0.10E+11, 0.10E+11/
! ..
        IF (nlu/=25) THEN
          PRINT *, 'number of land use classifications not correct '
          STOP
        END IF
        IF (dep_seasons/=5) THEN
          PRINT *, 'number of dep_seasons not correct '
          STOP
        END IF

!     SURFACE RESISTANCE DATA FOR DEPOSITION MODEL OF
!     M. L. WESELY, ATMOSPHERIC ENVIRONMENT 23 (1989) 1293-1304

!     Seasonal categories:
!     1: midsummer with lush vegetation
!     2: autumn with unharvested cropland
!     3: late autumn with frost, no snow
!     4: winter, snow on ground and subfreezing
!     5: transitional spring with partially green short annuals

!     Land use types:
!     USGS type                                Wesely type
!      1: Urban and built-up land              1
!      2: Dryland cropland and pasture         2
!      3: Irrigated cropland and pasture       2
!      4: Mix. dry/irrg. cropland and pasture  2
!      5: Cropland/grassland mosaic            2
!      6: Cropland/woodland mosaic             4
!      7: Grassland                            3
!      8: Shrubland                            3
!      9: Mixed shrubland/grassland            3
!     10: Savanna                              3, always summer
!     11: Deciduous broadleaf forest           4
!     12: Deciduous needleleaf forest          5, autumn and winter modi
!     13: Evergreen broadleaf forest           4, always summer
!     14: Evergreen needleleaf forest          5
!     15: Mixed Forest                         6
!     16: Water Bodies                         7
!     17: Herbaceous wetland                   9
!     18: Wooded wetland                       6
!     19: Barren or sparsely vegetated         8
!     20: Herbaceous Tundra                    9
!     21: Wooded Tundra                        6
!     22: Mixed Tundra                         6
!     23: Bare Ground Tundra                   8
!     24: Snow or Ice                          -, always winter
!     25: No data                              8


!     Order of data:
!      |
!      |   seasonal category
!     \|/
!     ---> landuse type
!     1       2       3       4       5       6       7       8       9
!     RLU for outer surfaces in the upper canopy
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            ri(iland,iseason) = dat1(iland,iseason)
          END DO
        END DO
!      data ((rlu(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
        DATA ((dat2(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.10E+11, 0.25E+04, 0.20E+04, 0.10E+11, &
          0.25E+04, 0.20E+04, 0.20E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, 0.90E+04, &
          0.20E+04, 0.40E+04, 0.80E+04, 0.10E+11, 0.90E+04, 0.80E+04, &
          0.10E+11, 0.90E+04, 0.80E+04, 0.80E+04, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, &
          0.90E+04, 0.20E+04, 0.40E+04, 0.80E+04, 0.10E+11, 0.90E+04, &
          0.80E+04, 0.10E+11, 0.90E+04, 0.80E+04, 0.80E+04, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.20E+04, 0.60E+04, 0.90E+04, 0.10E+11, &
          0.90E+04, 0.90E+04, 0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.40E+04, 0.40E+04, &
          0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, &
          0.20E+04, 0.40E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.30E+04, &
          0.10E+11, 0.40E+04, 0.30E+04, 0.10E+11, 0.40E+04, 0.30E+04, &
          0.30E+04, 0.10E+11, 0.10E+11, 0.10E+11/
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            rlu(iland,iseason) = dat2(iland,iseason)
          END DO
        END DO
!     RAC for transfer that depends on canopy height and density
!      data ((rac(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+03, &
        DATA ((dat3(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+04, 0.10E+03, &
          0.10E+03, 0.10E+03, 0.10E+03, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.00E+00, 0.30E+03, 0.20E+04, 0.00E+00, &
          0.30E+03, 0.20E+04, 0.20E+04, 0.00E+00, 0.00E+00, 0.00E+00, &
          0.10E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+04, &
          0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.15E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.17E+04, 0.00E+00, 0.20E+03, 0.17E+04, &
          0.00E+00, 0.20E+03, 0.17E+04, 0.17E+04, 0.00E+00, 0.00E+00, &
          0.00E+00, 0.10E+03, 0.10E+02, 0.10E+02, 0.10E+02, 0.10E+02, &
          0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+04, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, 0.00E+00, 0.10E+03, &
          0.15E+04, 0.00E+00, 0.10E+03, 0.15E+04, 0.15E+04, 0.00E+00, &
          0.00E+00, 0.00E+00, 0.10E+03, 0.10E+02, 0.10E+02, 0.10E+02, &
          0.10E+02, 0.10E+04, 0.10E+02, 0.10E+02, 0.10E+02, 0.10E+02, &
          0.10E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, 0.00E+00, &
          0.50E+02, 0.15E+04, 0.00E+00, 0.50E+02, 0.15E+04, 0.15E+04, &
          0.00E+00, 0.00E+00, 0.00E+00, 0.10E+03, 0.50E+02, 0.50E+02, &
          0.50E+02, 0.50E+02, 0.12E+04, 0.80E+02, 0.80E+02, 0.80E+02, &
          0.10E+03, 0.12E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.15E+04, &
          0.00E+00, 0.20E+03, 0.15E+04, 0.00E+00, 0.20E+03, 0.15E+04, &
          0.15E+04, 0.00E+00, 0.00E+00, 0.00E+00/
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            rac(iland,iseason) = dat3(iland,iseason)
          END DO
        END DO
!     RGSS for ground surface  SO2
!      data ((rgss(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.40E+03, &
        DATA ((dat4(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.40E+03, &
          0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.50E+03, 0.35E+03, &
          0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, 0.50E+03, 0.50E+03, &
          0.50E+03, 0.10E+03, 0.10E+01, 0.10E+01, 0.10E+03, 0.10E+04, &
          0.10E+01, 0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, 0.10E+04, &
          0.40E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.50E+03, &
          0.35E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, 0.50E+03, &
          0.50E+03, 0.50E+03, 0.10E+03, 0.10E+01, 0.10E+01, 0.10E+03, &
          0.10E+04, 0.10E+01, 0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, &
          0.10E+04, 0.40E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, &
          0.50E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.35E+03, 0.50E+03, &
          0.50E+03, 0.50E+03, 0.50E+03, 0.20E+03, 0.10E+01, 0.10E+01, &
          0.20E+03, 0.10E+04, 0.10E+01, 0.20E+03, 0.20E+03, 0.10E+04, &
          0.10E+03, 0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, &
          0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, 0.10E+03, &
          0.10E+03, 0.10E+03, 0.50E+03, 0.10E+03, 0.10E+03, 0.10E+01, &
          0.10E+03, 0.10E+03, 0.10E+04, 0.10E+03, 0.10E+03, 0.10E+03, &
          0.10E+04, 0.10E+03, 0.10E+04, 0.50E+03, 0.15E+03, 0.15E+03, &
          0.15E+03, 0.15E+03, 0.50E+03, 0.35E+03, 0.35E+03, 0.35E+03, &
          0.35E+03, 0.50E+03, 0.50E+03, 0.50E+03, 0.50E+03, 0.20E+03, &
          0.10E+01, 0.10E+01, 0.20E+03, 0.10E+04, 0.10E+01, 0.20E+03, &
          0.20E+03, 0.10E+04, 0.10E+03, 0.10E+04/
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            rgss(iland,iseason) = dat4(iland,iseason)
          END DO
        END DO
!     RGSO for ground surface  O3
!      data ((rgso(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.30E+03, &
        DATA ((dat5(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.30E+03, &
          0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.30E+03, 0.20E+04, 0.10E+04, 0.30E+03, 0.40E+03, &
          0.10E+04, 0.30E+03, 0.30E+03, 0.40E+03, 0.35E+04, 0.40E+03, &
          0.30E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.30E+03, 0.20E+04, 0.80E+03, 0.30E+03, &
          0.40E+03, 0.80E+03, 0.30E+03, 0.30E+03, 0.40E+03, 0.35E+04, &
          0.40E+03, 0.30E+03, 0.15E+03, 0.15E+03, 0.15E+03, 0.15E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.30E+03, 0.20E+04, 0.10E+04, &
          0.30E+03, 0.40E+03, 0.10E+04, 0.30E+03, 0.30E+03, 0.40E+03, &
          0.35E+04, 0.40E+03, 0.60E+03, 0.35E+04, 0.35E+04, 0.35E+04, &
          0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, 0.35E+04, &
          0.35E+04, 0.35E+04, 0.20E+03, 0.35E+04, 0.35E+04, 0.20E+04, &
          0.35E+04, 0.35E+04, 0.40E+03, 0.35E+04, 0.35E+04, 0.35E+04, &
          0.40E+03, 0.35E+04, 0.40E+03, 0.30E+03, 0.15E+03, 0.15E+03, &
          0.15E+03, 0.15E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, &
          0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.20E+03, 0.30E+03, &
          0.20E+04, 0.10E+04, 0.30E+03, 0.40E+03, 0.10E+04, 0.30E+03, &
          0.30E+03, 0.40E+03, 0.35E+04, 0.40E+03/
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            rgso(iland,iseason) = dat5(iland,iseason)
          END DO
        END DO
!     RCLS for exposed surfaces in the lower canopy  SO2
!      data ((rcls(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
        DATA ((dat6(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.20E+04, &
          0.20E+04, 0.20E+04, 0.10E+11, 0.25E+04, 0.20E+04, 0.10E+11, &
          0.25E+04, 0.20E+04, 0.20E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, 0.90E+04, &
          0.20E+04, 0.20E+04, 0.40E+04, 0.10E+11, 0.90E+04, 0.40E+04, &
          0.10E+11, 0.90E+04, 0.40E+04, 0.40E+04, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.90E+04, 0.90E+04, 0.90E+04, 0.90E+04, 0.20E+04, 0.90E+04, &
          0.90E+04, 0.20E+04, 0.30E+04, 0.60E+04, 0.10E+11, 0.90E+04, &
          0.60E+04, 0.10E+11, 0.90E+04, 0.60E+04, 0.60E+04, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.90E+04, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.90E+04, 0.90E+04, 0.20E+04, 0.20E+03, 0.40E+03, 0.10E+11, &
          0.90E+04, 0.40E+03, 0.10E+11, 0.90E+04, 0.40E+03, 0.40E+03, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.40E+04, 0.40E+04, &
          0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, 0.40E+04, &
          0.20E+04, 0.40E+04, 0.20E+04, 0.20E+04, 0.20E+04, 0.30E+04, &
          0.10E+11, 0.40E+04, 0.30E+04, 0.10E+11, 0.40E+04, 0.30E+04, &
          0.30E+04, 0.10E+11, 0.10E+11, 0.10E+11/
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            rcls(iland,iseason) = dat6(iland,iseason)
          END DO
        END DO
!     RCLO for exposed surfaces in the lower canopy  O3
!      data ((rclo(ILAND,ISEASON),ILAND=1,25),ISEASON=1,5)/0.10E+11, &
        DATA ((dat7(iland,iseason),iland=1,nlu),iseason=1,dep_seasons)/0.10E+11, &
          0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.10E+04, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+11, &
          0.10E+04, 0.10E+04, 0.10E+04, 0.10E+11, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, &
          0.40E+03, 0.40E+03, 0.40E+03, 0.10E+04, 0.40E+03, 0.40E+03, &
          0.10E+04, 0.10E+04, 0.60E+03, 0.10E+11, 0.40E+03, 0.60E+03, &
          0.10E+11, 0.40E+03, 0.60E+03, 0.60E+03, 0.10E+11, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.40E+03, 0.40E+03, 0.40E+03, 0.40E+03, 0.10E+04, 0.40E+03, &
          0.40E+03, 0.10E+04, 0.10E+04, 0.60E+03, 0.10E+11, 0.80E+03, &
          0.60E+03, 0.10E+11, 0.80E+03, 0.60E+03, 0.60E+03, 0.10E+11, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.40E+03, 0.10E+04, 0.10E+04, 0.10E+04, 0.10E+04, &
          0.40E+03, 0.40E+03, 0.10E+04, 0.15E+04, 0.60E+03, 0.10E+11, &
          0.80E+03, 0.60E+03, 0.10E+11, 0.80E+03, 0.60E+03, 0.60E+03, &
          0.10E+11, 0.10E+11, 0.10E+11, 0.10E+11, 0.10E+04, 0.10E+04, &
          0.10E+04, 0.10E+04, 0.50E+03, 0.50E+03, 0.50E+03, 0.50E+03, &
          0.10E+04, 0.50E+03, 0.15E+04, 0.10E+04, 0.15E+04, 0.70E+03, &
          0.10E+11, 0.60E+03, 0.70E+03, 0.10E+11, 0.60E+03, 0.70E+03, &
          0.70E+03, 0.10E+11, 0.10E+11, 0.10E+11/
        DO iseason = 1, dep_seasons
          DO iland = 1, nlu
            rclo(iland,iseason) = dat7(iland,iseason)
          END DO
        END DO
        
	DO l = 1, nspecies_chem
	  ! if(dvj(l) < 0.) dvj(l) = 1./(sqrt(weight(l)))
!          hstar4(l) = hstar(l) ! preliminary              
! Correction of diff. coeff
!          dvj(l) = dvj(l)*(293.15/298.15)**1.75
!          sc = 0.15/dvj(l) ! Schmidt Number at 20°C   
          dratio(l) = 0.242/dvj(l) !    ! of water vapor and gas at
! Ratio of diffusion coeffi
!          scpr23(l) = (sc/0.72)**(2./3.) ! (Schmidt # / Prandtl #)**
        END DO


!     DATA FOR AEROSOL PARTICLE DEPOSITION FOR THE MODEL OF
!     J. W. ERISMAN, A. VAN PUL AND P. WYERS
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607

!     vd = (u* / k) * CORRECTION FACTORS

!     CONSTANT K FOR LANDUSE TYPES:
! urban and built-up land                  
        kpart(19) = 500.
! dryland cropland and pasture             
        kpart(2) = 500.
! irrigated cropland and pasture           
        kpart(16) = 500.
! mixed dryland/irrigated cropland and past
        kpart(16) = 500.
! cropland/grassland mosaic                
        kpart(15) = 500.
! cropland/woodland mosaic                 
        kpart(14) = 100.
! grassland                                
        kpart(8) = 500.
! shrubland                                
        kpart(13) = 500.
! mixed shrubland/grassland                
        kpart(9) = 500.
! savanna                                  
        kpart(9) = 500.
! deciduous broadleaf forest               
        kpart(5) = 100.
! deciduous needleleaf forest              
        kpart(6) = 100.
! evergreen broadleaf forest               
        kpart(7) = 100.
! evergreen needleleaf forest              
        kpart(4) = 100.
! mixed forest                             
        kpart(14) = 100.
! water bodies                             
        kpart(1) = 500.
! herbaceous wetland                       
        kpart(17) = 500.
! wooded wetland                           
        kpart(18) = 500.
! barren or sparsely vegetated             
        kpart(10) = 500.
! herbaceous tundra                        
        kpart(11) = 500.
! wooded tundra                            
        kpart(11) = 100.
! mixed tundra                             
        kpart(11) = 500.
! bare ground tundra                       
        kpart(3) = 500.
! snow or ice                              
        kpart(2) = 500.
!     Comments:
        kpart(25) = 500.
!     Erisman et al. (1994) give
!     k = 500 for low vegetation and k = 100 for forests.

!     For desert k = 500 is taken according to measurements
!     on bare soil by
!     J. Fontan, A. Lopez, E. Lamaud and A. Druilhet (1997)
!     "Vertical Flux Measurements of the Submicronic Aerosol Particles
!     and Parametrisation of the Dry Deposition Velocity"
!     in: "Biosphere-Atmosphere Exchange of Pollutants
!     and Trace Substances"
!     Editor: S. Slanina. Springer-Verlag Berlin, Heidelberg, 1997
!     pp. 381-390

!     For coniferous forest the Erisman value of  k = 100 is taken.
!     Measurements of Erisman et al. (1997) in a coniferous forest
!     in the Netherlands, lead to values of k between 20 and 38
!     (Atmospheric Environment 31 (1997), 321-332).
!     However, these high values of vd may be reached during
!     instable cases. The eddy correlation measurements
!     of Gallagher et al. (1997) made during the same experiment
!     show for stable cases (L>0) values of k between 200 and 250
!     at minimum (Atmospheric Environment 31 (1997), 359-373).
!     Fontan et al. (1997) found k = 250 in a forest
!     of maritime pine in southwestern France.

!     For gras, model calculations of Davidson et al. support
!     the value of 500.
!     C. I. Davidson, J. M. Miller and M. A. Pleskov
!     "The Influence of Surface Structure on Predicted Particles
!     Dry Deposition to Natural Gras Canopies"
!     Water, Air, and Soil Pollution 18 (1982) 25-43

!     Snow covered surface: The experiment of Ibrahim et al. (1983)
!     gives k = 436 for 0.7 um diameter particles.
!     The deposition velocity of Milford and Davidson (1987)
!     gives k = 154 for continental sulfate aerosol.
!     M. Ibrahim, L. A. Barrie and F. Fanaki
!     Atmospheric Environment 17 (1983), 781-788

!     J. B. Milford and C. I. Davidson
!     "The Sizes of Particulate Sulfate and Nitrate in the Atmosphere
!     - A Review"
!     JAPCA 37 (1987), 125-134
!
!------------!-------------!----------------------------------------------------------------
!  LEAF-3 CLASS (20)                      !  Wesely type !   Land use types:
!---AND DESCRIPTION 			  !		 !   USGS type  			 
!  0  Ocean                               !   7 	 !    1: Urban and built-up land	      
!  1  Lakes, rivers, streams              !   7 	 !    2: Dryland cropland and pasture	     
!  2  Ice cap/glacier                     !   - 	 !    3: Irrigated cropland and pasture      
!  3  Desert, bare soil                   !   8 	 !    4: Mix. dry/irrg. cropland and pasture 
!  4  Evergreen needleleaf tree           !   5 	 !    5: Cropland/grassland mosaic	     
!  5  Deciduous needleleaf tree           !   5 	 !    6: Cropland/woodland mosaic	     
!  6  Deciduous broadleaf tree            !   4 	 !    7: Grassland			     
!  7  Evergreen broadleaf tree            !   4 	 !    8: Shrubland			     
!  8  Short grass                         !   2 	 !    9: Mixed shrubland/grassland				
!  9  Tall grass                          !   3 	 !   10: Savanna			     
! 10  Semi-desert                         !   8 	 !   11: Deciduous broadleaf forest	     
! 11  Tundra                              !   8 	 !   12: Deciduous needleleaf forest	     
! 12  Evergreen shrub                     !   3 	 !   13: Evergreen broadleaf forest	     
! 13  Deciduous shrub                     !   3 	 !   14: Evergreen needleleaf forest	     
! 14  Mixed woodland                      !   6 	 !   15: Mixed Forest			     
! 15  Crop/mixed farming, C3 grassland    !   2 	 !   16: Water Bodies			     
! 16  Irrigated crop                      !   2 	 !   17: Herbaceous wetland		     
! 17  Bog or marsh                        !   9 	 !   18: Wooded wetland 		     
! 18  Wooded grassland                    !   6 	 !   19: Barren or sparsely vegetated	     
! 19  Urban and built up                  !   1 	 !   20: Herbaceous Tundra		     
! 20  Wetland evergreen broadleaf tree    !   6 	 !   21: Wooded Tundra  		     
                                          !		 !   22: Mixed Tundra			     
					  !		 !   23: Bare Ground Tundra		     
					  !	 	 !   24: Snow or Ice			       
!--------------------------------------------------------------------------------------------
          ixxxlu(1)  =7 
          ixxxlu(2)  =0  
          ixxxlu(3)  =8 
          ixxxlu(4)  =5 
          ixxxlu(5)  =5 
          ixxxlu(6)  =4 
          ixxxlu(7)  =4 
          ixxxlu(8)  =2 
          ixxxlu(9)  =3 
          ixxxlu(10) =8 
          ixxxlu(11) =8 
          ixxxlu(12) =3 
          ixxxlu(13) =3 
          ixxxlu(14) =6 
          ixxxlu(15) =2 
          ixxxlu(16) =2 
          ixxxlu(17) =9 
          ixxxlu(18) =6 
          ixxxlu(19) =1 
          ixxxlu(20) =6 


   end      SUBROUTINE dep_init

! **********************************************************************
      SUBROUTINE deppart(rmol,ustar,rh,clw,iland,dvpart,dvfog)
!     THIS SUBROUTINE CALCULATES SURFACE DEPOSITION VELOCITIES
!     FOR FINE AEROSOL PARTICLES ACCORDING TO THE MODEL OF
!     J. W. ERISMAN, A. VAN PUL, AND P. WYERS,
!     ATMOSPHERIC ENVIRONMENT 28 (1994), 2595-2607
!     WRITTEN BY WINFRIED SEIDL, APRIL 1997
!     MODIFIED BY WINFRIED SEIDL, MARCH 2000
!            FOR MM5 VERSION 3
! ----------------------------------------------------------------------
        implicit none
! .. Scalar Arguments ..
        REAL :: clw, dvfog, dvpart, rh, rmol, ustar
        INTEGER :: iland
! ..
! .. Intrinsic Functions ..
        INTRINSIC exp
! ..
        dvpart = ustar/kpart(iland)
	!print*,ustar,iland,kpart(iland)
	
        IF (rmol<0.) THEN
!         INSTABLE LAYERING CORRECTION
          dvpart = dvpart*(1.+(-300.*rmol)**0.66667)
        END IF
        IF (rh>80.) THEN
!         HIGH RELATIVE HUMIDITY CORRECTION
!         ACCORDING TO J. W. ERISMAN ET AL.
!         ATMOSPHERIC ENVIRONMENT 31 (1997), 321-332
          dvpart = dvpart*(1.+0.37*exp((rh-80.)/20.))
        END IF
        !print*,'dvpart=', dvpart
!srf - not using fog deposition (in case on, set clw,iland=iveg)
        return

!       SEDIMENTATION VELOCITY OF FOG WATER ACCORDING TO
!       R. FORKEL, W. SEIDL, R. DLUGI AND E. DEIGELE
!       J. GEOPHYS. RES. 95D (1990), 18501-18515
        dvfog = 0.06*clw
        IF (ixxxlu(iland)==5) THEN
!         TURBULENT DEPOSITION OF FOG WATER IN CONIFEROUS FOREST ACCORDI
!         A. T. VERMEULEN ET AL.
!         ATMOSPHERIC ENVIRONMENT 31 (1997), 375-386
          dvfog = dvfog + 0.195*ustar*ustar
        END IF

      END SUBROUTINE deppart
! **********************************************************************
      SUBROUTINE depvel(numchem,rmol,zr,z0,ustar,vgpart)!,polint)
!     THIS FUNCTION HAS BEEN DESIGNED TO EVALUATE AN UPPER LIMIT
!     FOR THE POLLUTANT DEPOSITION VELOCITY AS A FUNCTION OF THE
!     SURFACE ROUGHNESS AND METEOROLOGICAL CONDITIONS.
!     PROGRAM WRITTEN BY GREGORY J.MCRAE (NOVEMBER 1977)
!         Modified by Darrell A. Winner  (Feb. 1991)
!                  by Winfried Seidl     (Aug. 1997)
!.....PROGRAM VARIABLES...
!     RMOL     - RECIPROCAL OF THE MONIN-OBUKHOV LENGTH
!     ZR       - REFERENCE HEIGHT
!     Z0       - SURFACE ROUGHNESS HEIGHT
!     SCPR23   - (Schmidt #/Prandtl #)**(2/3) Diffusion correction fact
!     UBAR     - ABSOLUTE VALUE OF SURFACE WIND SPEED
!     DEPVEL   - POLLUTANT DEPOSITION VELOCITY
!     Vk       - VON KARMAN CONSTANT
!     USTAR    - FRICTION VELOCITY U*
!     POLINT   - POLLUTANT INTEGRAL
!.....REFERENCES...
!     MCRAE, G.J. ET AL. (1983) MATHEMATICAL MODELING OF PHOTOCHEMICAL
!       AIR POLLUTION, ENVIRONMENTAL QUALITY LABORATORY REPORT 18,
!       CALIFORNIA INSTITUTE OF TECHNOLOGY, PASADENA, CALIFORNIA.
!.....RESTRICTIONS...
!     1. THE MODEL EDDY DIFFUSIVITIES ARE BASED ON MONIN-OBUKHOV
!        SIMILARITY THEORY AND SO ARE ONLY APPLICABLE IN THE
!        SURFACE LAYER, A HEIGHT OF O(30M).
!     2. ALL INPUT UNITS MUST BE CONSISTENT
!     3. THE PHI FUNCTIONS USED TO CALCULATE THE FRICTION
!        VELOCITY U* AND THE POLLUTANT INTEGRALS ARE BASED
!        ON THE WORK OF BUSINGER ET AL.(1971).
!     4. THE MOMENTUM AND POLLUTANT DIFFUSIVITIES ARE NOT
!        THE SAME FOR THE CASES L<0 AND L>0.
       implicit none
!     Local Variables
! .. Scalar Arguments ..
        REAL :: rmol, ustar, vgpart, z0, zr
        INTEGER :: numchem
! ..
! .. Array Arguments ..
!        REAL :: depv(numchem)
! ..
! .. Local Scalars ..
        REAL :: ao, ar, polint
        INTEGER :: l
! ..
! .. Intrinsic Functions ..
        INTRINSIC alog
! ..
!     Set the von Karman constant
        real, parameter :: vk = 0.4

!     DETERMINE THE STABILITY BASED ON THE CONDITIONS
!             1/L < 0 UNSTABLE
!             1/L = 0 NEUTRAL
!             1/L > 0 STABLE

        IF (rmol<0) THEN
          ar = ((1.0-9.0*zr*rmol)**(0.25)+0.001)**2
          ao = ((1.0-9.0*z0*rmol)**(0.25)+0.001)**2
          polint = 0.74*(alog((ar-1.0)/(ar+1.0))-alog((ao-1.0)/(ao+1.0)))
        ELSE IF (rmol==0) THEN
          polint = 0.74*alog(zr/z0)
        ELSE
          polint = 0.74*alog(zr/z0) + 4.7*rmol*(zr-z0)
        END IF

!     CALCULATE THE Maximum DEPOSITION VELOCITY

        !DO l = 1, numchem
        !  depv(l) = ustar*vk/(2.0*scpr23(l)+polint)
        !END DO
        vgpart = ustar*vk/polint
        !print*,polint,ustar,vgpart,z0,ustar*vk/polint
        RETURN
      END SUBROUTINE depvel

!========================================================================

! -----------------   AEROSOL DRY DEP AND SEDIMENTATION -----------------

!========================================================================

subroutine dry_dep_sedim_particles(ngrids,ngrid,m1,m2,m3,naddsc,npatch,ia,iz,ja,jz &
                            ,r_aer               &
			    ,temp3d,air_dens3d,temps,dens,vels,rvs,Zi                  &
			    ,ustar,tstar                             &
			    ,patch_area,veg,Z0m)

implicit none
integer :: m1,m2,m3,naddsc,npatch,ia,iz,ja,jz,i,j,ipatch,ngrids,ngrid,ispc
real vdtmp
real, dimension(m2,m3)        :: temps,dens,vels,rvs,Zi
real, dimension(m2,m3,npatch) :: ustar,tstar,patch_area,veg,Z0m

real, dimension(m2,m3,npatch) :: r_aer
real, dimension(m1,m2,m3)  :: temp3d,air_dens3d

 
 
if(.not. aer_alloc) then
   call alloc_aer_sedim(m1,m2,m3,npatch,ngrids) 
end if
 
if( NAER_TRANSPORTED == 0 ) return

!- sedimentation  parameterization 
 call sedim_particles_3d(ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,temp3d,air_dens3d)


!- dry deposition parameterization

!- laminar sub-layer resistance
 call lsl_particles(ngrid,m2,m3,npatch,ia,iz,ja,jz &
                  ,temps,dens,vels,rvs,Zi,ustar,tstar,patch_area,veg,Z0m)

!- particles deposition velocity (m/s)

do j = ja,jz
  do i = ia,iz
    do ispc=naer_a,naer_z

      dd_sedim(ispc,ngrid)%v_dep_part(i,j) = 0.
      
      do ipatch = 1,npatch

        if (patch_area(i,j,ipatch) .ge. .009) then


   	 vdtmp = dd_sedim(ispc,ngrid)%v_sed_part(1,i,j)				      + &
   	 	 1./( r_aer(i,j,ipatch) + dd_sedim(ispc,ngrid)%r_lsl_part(i,j,ipatch) + &
   		      r_aer(i,j,ipatch) * dd_sedim(ispc,ngrid)%r_lsl_part(i,j,ipatch) * &
   		      dd_sedim(ispc,ngrid)%v_sed_part(1,i,j)				  ) 

   	 dd_sedim(ispc,ngrid)%v_dep_part(i,j) = dd_sedim(ispc,ngrid)%v_dep_part(i,j)  + &
   					        patch_area(i,j,ipatch)*vdtmp

        endif
      enddo		    

     !print*,'V_dep (cm/s)',i,j,dd_sedim(ispc,ngrid)%v_dep_part(i,j)*100.
     !print*,'r_lsl-v_sed =', dd_sedim(ispc,ngrid)%r_lsl_part(i,j,1:npatch), &
     !	dd_sedim(ispc,ngrid)%v_sed_part(1,i,j) !,r_aer(i,j,ipatch)

    enddo		   
 enddo
enddo

return

!- replace the sedimentation velocity at surface with the deposition velocity
!- (which includes v_sed_part plus the turbulent flux in laminar sub-layer)

do j = ja,jz
  do i = ia,iz
    do ispc=naer_a,naer_z
        dd_sedim(ispc,ngrid)%v_sed_part(1,i,j) = dd_sedim(ispc,ngrid)%v_dep_part(i,j)
    enddo
  enddo
enddo  


end subroutine dry_dep_sedim_particles
!========================================================================

subroutine lsl_particles(ng,m2,m3,npatch,ia,iz,ja,jz&
                        ,temps,dens,vels,rvs,Zi,ustar,tstar,patch_area,veg,Z0m)

  use aer1_list, only : part_radius,part_dens
  implicit none
  REAL,PARAMETER :: kB = 1.3807e-23      ! const Boltzmann - kg m^2 s^-2 K^-1 molecule^-1
  REAL,PARAMETER :: ASP = 1.257          ! 1.249
  REAL,PARAMETER :: BSP = 0.4            ! 0.42
  REAL,PARAMETER :: CSP = 1.1            ! 0.87 

  REAL,PARAMETER :: M_AVEG = 4.8096e-26  ! average mass of one molecure - kg  molecule^-1
  REAL,PARAMETER :: vonK = 0.40          ! von Karman constant
  REAL,PARAMETER :: Cpd = 1004.67        ! specific heat of dry air [J/kg/K]
  REAL,PARAMETER :: em23 = -2./3., ep23 = +2./3., ep13=1./3.    ! exponents 2./3. 1/3
  REAL,PARAMETER :: pi = 3.1415927, g = 9.80


integer :: ng,m2,m3,npatch,ia,iz,ja,jz,i,j,ipatch,ispc
real, dimension(m2,m3) :: temps,dens,vels,rvs,Zi
real, dimension(m2,m3,npatch) :: ustar,tstar,patch_area,veg,Z0m
real wptp,wstar
!
! dd_sedim(ispc,ng)%v_sed      == particle sedimentation velocity   (m/s)
! dd_sedim(ispc,ng)%r_lsl_part == resistance to molecular diffusion (s/m)
!
real Kn,n_air,Gi,v_air  ,mfp,D,nu,Sc,St,Kd,Dh,Z0h,Pr    

real,parameter :: limite = -30. 

do j = ja,jz
  do i = ia,iz

   !- several particle/environment properties

   !- mean speed of air molecules (m/s)
   !  v_air = sqrt( 8. * kB   * temps(i,j) / (pi * M_AVEG) )
      v_air = sqrt( 7.3102e+2 * temps(i,j)			  )

   !-dynamic viscosity of air (kg m^-1 s^-1)
   !  n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16)**1.5
   !optimized version
      n_air = 1.8325e-5*(416.16/(temps(i,j)+120.))*(temps(i,j)/296.16) * & 
 	      sqrt(temps(i,j)/296.16)
      
   !- mean free path of an air molecule (m)
      mfp = 2.* n_air /(dens(i,j)*v_air)

   !-  kinematic viscosity of air ()
      nu = n_air/dens(i,j)
      
      do ispc=naer_a,naer_z
  
       !- Knudsen number
       Kn = mfp/part_radius(ind_mode(ispc),ind_aer(ispc))


       !- Slip correction factor (Gi)
       Gi = 1. + Kn*( ASP + BSP*exp(-CSP/Kn) )

       !- Schmidt number determination (Jacobson)
       !-- molecular diffusivity (Brownian diffusivity coeficient)   
       !    D =  (kB/M_AVEG) * temps(i,j) * Gi / (6.*pi*part_radius*n_air)

       D =  kB * temps(i,j) * Gi / (6.*pi*part_radius(ind_mode(ispc),ind_aer(ispc))*n_air)


       !-  Schmidt number
       Sc = nu/D

         !- laminar sub-layer resistance for particles (s m^-1)
         do ipatch=1,npatch

          !-  Stokes number determination (Slinn, 1980)
          !-  St =     ustar^2 * V_sedim /( g * kinematic viscosity of air)    
          St = max(ustar(i,j,ipatch)**2. * dd_sedim(ispc,ng)%v_sed_part(1,i,j) &
	       / (g * nu), 0.01) 
    
          !-  laminar sub-layer resistance for particles (s m^-1) 

          if(ipatch == 1) then !- water

               !-from Slinn 1980 (Atmos. Env.)
               !dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = ( vonK * vels(i,j)/ustar(i,j,ipatch)**2. ) / &
               !	                  ( 1./sqrt(Sc) + 10.**(-3./St) )
                dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = ( vonK * vels(i,j)/ustar(i,j,ipatch)**2. ) / &
	                            ( 1./sqrt(Sc) + 10.**max((-3./St),limite) )
			  
			  
          else
               if(patch_area(i,j,ipatch) > 0.009) then
          
               !-- for smooth land surfaces !- bare ground (3), desert(3) or ice (2)
               !-- Seinfeld & Pandis (1998)
	         if(nint(veg(i,j,ipatch)) == 3 ) then 
            
                    ! dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = 1./(ustar(i,j,ipatch) * (Sc**em23 + 10.**(-3./St)))
                      dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = 1./(ustar(i,j,ipatch) * &
                                     (Sc**em23 + 10.**max((-3./St),limite)))
	                      
                                       
                 else
                    !- lsl resistance according  Jacobson(1999)
                    !
                    !- thermal conductivity of dry air (Kd)
	            !Kd = 0.023807 + 7.1128e-5*(temps(i,j) - 273.15) !- Eq.(2.3) 
                    !- Prandt number
                    !Pr =  n_air*Cpd*(1.+0.859*rvs(i,j))/Kd           !- Eq.(17.32)  
	   
                    !- energy moisture roughness lengths (Z0h)                 !- Eq.(8.10)
                    !-- molecular thermal diffusion coeff. (m^2 s^-1)
                    !Dh=Kd/(dens(i,j)*Cpd)
                    !- Z0h	   
                    !Z0h=Dh/(vonK*ustar(i,j,ipatch))
                    !- lsl resistance according Jacobson Eq. (20.14)
                    ! dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) =  log( Z0m(i,j,ipatch)/Z0h )* ( (Sc/Pr)**ep23 ) &
                    !	                     /  ( vonK*ustar(i,j,ipatch) )

                    !---------  		 
		    !- lsl resistance according :
                    !- from Wesely et al. (1985) following Slinn (1982)
                    !- also Binkowski & Shankar, JGR, 100,D12,26191-26209, 1995
                    !- Rb= (u* (1+ (w*/u*)^2)) (Sc^2./3 + 10^(-3/St))	 

                    wptp  = - ustar(i,j,ipatch) * tstar(i,j,ipatch) ! sensible heat flux
	            wstar = ( max (0., g* Zi(i,j)*  wptp/temps(i,j) ) )**ep13  
                    !	     
		    !dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = 1./(  					  &
                    !				     ustar(i,j,ipatch)* 			 &
                    !				     (1. + 0.24*(wstar/ustar(i,j,ipatch))**2.)*  &
                    !				     ( (Sc**em23 + 10.**(-3./St)) )		 &
                    !				    )
		     dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch) = 1./(						&
                        		    ustar(i,j,ipatch)*  			&
                        		    (1. + 0.24*(wstar/ustar(i,j,ipatch))**2.)*  &
                        		    ( (Sc**em23 + 10.**max((-3./St),limite) ) ) )

!       print*,'================'
!       print*,'LSL_WE',i,j,ipatch,nint(veg(i,j,ipatch)),patch_area(i,j,ipatch)
!	PRINT*,wstar,ustar(i,j,ipatch),r_lsl(i,j,ipatch)

                 endif
               endif
        endif


        !print*,'laminar sub-layer resistance'
        !print*,'LSL_J',i,j,ipatch,dd_sedim(ispc,ng)%r_lsl_part(i,j,ipatch)
        !print*,'ZOM',Z0m(i,j,ipatch)
        !print*,'Z0H',Z0h
        !print*,'SC PR',Sc,Pr
        !print*,'u*',ustar(i,j,ipatch)
	
      enddo
    enddo

  enddo
enddo

end subroutine lsl_particles

!=========================================================================
  subroutine alloc_aer_sedim(m1,m2,m3,npatch,ngrids)
    use aer1_list, only:  nspecies_aer=> nspecies &
                         ,spc_alloc_aer=> spc_alloc &
			 ,mode_alloc, nucle, accum, coarse &
			 ,transport,on,off
    use mem_aer1, only : AEROSOL
    use node_mod, only : mmxp,mmyp,mmzp
    
    implicit none
    integer, intent (IN):: m1,m2,m3,ngrids,npatch
    integer :: ispc,imode,ng


    if(aer_alloc) then
       print *,'ERROR: aer_alloc already allocated'
       print *,'Routine: aer_alloc File: chem_dry_dep.f90'
       stop
    end if

    NAER_TRANSPORTED = 0
    
    ! aerosol section : mapping  
    naer_a = 1
    naer_z = 0
      
      do ispc = 1,nspecies_aer
         do imode=1,nmodes
	 
           if(mode_alloc   (imode,ispc) == on ) then 
	      
	      naer_z = naer_z + 1
	      ind_aer (naer_z) = ispc
	      ind_mode(naer_z) = imode
	       
            endif
         enddo
      enddo 
     
    ! total number of species (aer) to be transported 
    NAER_TRANSPORTED  =  naer_z
    if(NAER_TRANSPORTED == 0) return

    allocate (dd_sedim(NAER_TRANSPORTED,ngrids))
    do ng=1,ngrids
     do ispc=1,NAER_TRANSPORTED
         allocate(dd_sedim(ispc,ng)%v_sed_part (mmzp(ng),mmxp(ng),mmyp(ng))     ); dd_sedim(ispc,ng)%v_sed_part = 0.
         allocate(dd_sedim(ispc,ng)%r_lsl_part (mmxp(ng),mmyp(ng),npatch) ); dd_sedim(ispc,ng)%r_lsl_part = 0.
         allocate(dd_sedim(ispc,ng)%v_dep_part (mmxp(ng),mmyp(ng)       ) ); dd_sedim(ispc,ng)%v_dep_part = 0.
     enddo
    enddo	    
	     

    aer_alloc=.true.

  end subroutine alloc_aer_sedim
!------------------------------------------------------------------------

subroutine sedim_particles_3d(ngrid,m1,m2,m3,npatch,ia,iz,ja,jz,temp3d,air_dens3d)
  
  use aer1_list, only : part_radius,part_dens
  use rconstants,  only :  g
  implicit none

  integer :: ngrid,m1,m2,m3,naddsc,npatch,ia,iz,ja,jz,i,j,ipatch,k,ispc
  real, dimension(m1,m2,m3) :: temp3d,air_dens3d
  REAL,PARAMETER :: ASP = 1.257          ! 1.249
  REAL,PARAMETER :: BSP = 0.4            ! 0.42
  REAL,PARAMETER :: CSP = 1.1            ! 0.87 
  REAL,PARAMETER :: onesix = 1./6.         
  
!
!For small particles diameter (<20 micrometer) (low Reynolds number)
!the fall velocity is given by (Seinfeld & Pandis, Jacobson)
! V_s = 2 r**2 (rho_p - rho_air) * g * Gi / 9 n_air
! where:
! r == radius of particle (m)
! rho_p,a = density of particle, air (kg/m^3)
! g = gravity accel (9.8 m/s^2)
! Gi = Kn(A' + B' + C' exp(-C'/Kn) , Kn = Knudsen number
! n_air = dynamic viscosity of air
!- constantes
!parameter (ASP=1.257,BSP=0.4,CSP=1.1) ! Seinfeld & Pandis 
!parameter (kB = 1.3807e-23)     ! const Boltzmann - kg m^2 s^-2 K^-1 molecule^-1
!parameter (M_AVEG = 4.8096e-26) ! average mass of one molecure - kg  molecule^-1
!
!real, dimension(m1,m2,m3,nmodes,aer_nspecies)    :: v_sed

real A,B,C,Kn,n_air,Gi,v_air,mfp,nu,re,z,b0,x,y


do j = ja,jz
  do i = ia,iz
    do k=1,m1  

     !- several particle/environment properties
     
     !- mean speed of air molecules (m/s)
     !  v_air = sqrt( 8. * kB	* temp3d(k,i,j) / (pi * M_AVEG) )
 	v_air = sqrt( 7.3102e+2 * temp3d(k,i,j)  		      )

     !-dynamic viscosity of air (kg m^-1 s^-1)
     !  n_air = 1.8325e-5*(416.16/(temp3d(k,i,j)+120.))*(temp3d(k,i,j)/296.16)**1.5
     !optimized version
 	n_air = 1.8325e-5*(416.16/(temp3d(k,i,j)+120.))*(temp3d(k,i,j)/296.16) * & 
 		sqrt(temp3d(k,i,j)/296.16)   
     !- kinematic viscosity of air ()
 	nu = n_air/air_dens3d(k,i,j)
     !- mean free path of an air molecule (m)
 	mfp = 2.* n_air /(air_dens3d(k,i,j)*v_air)
      
     do ispc=naer_a,naer_z
 
       !print*,'particle=',ispc
    
 
       !- Knudsen number
       Kn = mfp/part_radius(ind_mode(ispc),ind_aer(ispc))

       !- Slip correction factor (Gi)
       Gi = 1. + Kn*( ASP + BSP*exp(-CSP/Kn) )


       !- This is regime 1 (Pruppacher and Klett (chap. 10)) , first guess for terminal velocity
       !- 0.5 < part_radius < 10 micrometers  / 1.e-6 < Re < 1e-2
       !  
       !- particle sedimentation velocity (m/s)
       !  v_sed(k,i,j) = (2./9.)*g*part_radius**2  *(part_dens-air_dens3d(k,i,j))*Gi/n_air
       !- opt
          dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) = 2.18*part_radius(ind_mode(ispc),ind_aer(ispc))**2.  &
	                      *(part_dens(ind_mode(ispc),ind_aer(ispc)))*Gi/n_air !part_dens >>dens_air
       
       
       ! Reynolds number
       !   re = 2*part_radius*v_sed_part(i,j)/nu
           re = 2. * air_dens3d(k,i,j) * part_radius(ind_mode(ispc),ind_aer(ispc)) * &
	        dd_sedim(ispc,ngrid)%v_sed_part(k,i,j)/n_air
       

       if( re .ge. 0.01 .and. re .le. 300.) then
       !  This is "regime 2" in Pruppacher and Klett (chap. 10).

            x = log(24.*re/Gi)
            y = -3.18657 + x*(0.992696    - x*(.00153193 &
             		 + x*(0.000987059 + x*(.000578878&
             		 - x*(8.55176E-05 - x* 3.27815E-06 )))))
            !if( y .lt. -675. ) y = -675.
            !if( y .ge.  741. ) y =  741.
            y=min(max(y,-675.),741.)
	    
            re = exp(y)*Gi

            dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) = re * n_air / &
	                       (2.*part_radius(ind_mode(ispc),ind_aer(ispc))* air_dens3d(k,i,j))
       endif
       
       if( re > 300. )then
        !  This is "regime 3" in Pruppacher and Klett (chap. 10).
            z  = ((1.e6*air_dens3d(k,i,j)**2) &
                /  (g * part_dens(ind_mode(ispc),ind_aer(ispc)) * n_air**4) )**(onesix)
              b0 = (24.*dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) *n_air)/100.
              x  = log(z*b0)
              y  = -5.00015 + x*(5.23778   - x*(2.04914 - x*(0.475294 &
                           - x*(0.0542819 - x* 0.00238449 ))))
              !if( y .lt. -675. )  y = -675.0
              !if( y .ge.  741. )  y =  741.0
              y=min(max(y,-675.),741.)
              re = z*exp(y)*Gi
 
              dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) = re * n_air / &
	                         (2.*part_radius(ind_mode(ispc),ind_aer(ispc)) * air_dens3d(k,i,j))

        endif
	!tmp
	!for testing 
	!dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) = float(ispc)/500.
	!if(k==1 .or. K==2)print*,'v-sed=',k,dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) 
	!tmp	
	
       
      enddo ! ispc loop
     enddo  !   k loop    
  enddo
enddo
return
end subroutine sedim_particles_3d
!------------------------------------------------------------------------
subroutine fa_preptc_with_sedim(m1,m2,m3,vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
                                 ,vt3dh,vt3di,vt3dj,vt3dk,dn0,dn0u,dn0v  &
                                 ,rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t  &
                                 ,dxu,dyv,dxt,dyt,mynum &
				 ,dtlt,N_current_scalar,num_scalar_aer_1st, wp,wc,vt3dp)

use mem_grid, only : ngrid,hw4,dzm,dzt

implicit none

integer :: m1,m2,m3,j,i,k,im,ip,jm,jp,mynum

real :: c1,c2,c3,c4,rtgti

real, dimension(m1,m2,m3) :: vt3da,vt3db,vt3dc,vt3dd,vt3de,vt3df  &
   ,vt3dh,vt3di,vt3dj,vt3dk,dn0,dn0u,dn0v

real, dimension(m2,m3) :: rtgt,rtgu,rtgv,fmapt,fmapui,fmapvi,f13t,f23t  &
   ,dxu,dyv,dxt,dyt

! VT3DA, VT3DB, and VT3DC are input as the velocity components (averaged
! between past and current time levels) times dtlt.

! Add contribution to VT3DC from horiz winds crossing sloping sigma surfaces,
!    and include 1/rtgt factor in VT3DC
! Compute half Courant numbers: VT3DD, VT3DE, and VT3DF
! Compute weight at scalar point: VT3DH
! Compute advective weights for the linear term: VT3DI, VCTR1, and VCTR2


!srf- aerosol section 
integer, intent(in) :: num_scalar_aer_1st,N_current_scalar
real, intent(IN) :: dtlt
real, dimension(m1,m2,m3) :: wp,wc,vt3dp
integer ISPC

!------------------ aerosol mapping
! num_scalar_aer_1st                  corresponds do aerosol  1 (=naer_a)
! num_scalar_aer_1st +1               corresponds do aerosol  2
! ...
! num_scalar_aer_1st +NAER_TRANSPORTED corresponds do aerosol naer_z
! in this case, to access the correct V_SED use "ISPC"
ISPC = N_current_scalar - num_scalar_aer_1st + 1

! sedim is included at vertical direction (vt3dc)
do j = 1,m3
   do i = 1,m2
      do k = 1,m1
         !if(j==20 .and. i==20) print*,'particle',ISPC,sedim(ispc,ngrid)%v_sed_part(k,i,j),wp(k,i,j)
	 vt3dc(k,i,j) = ( ( wp(k,i,j) + wc(k,i,j) )* 0.5 - dd_sedim(ispc,ngrid)%v_sed_part(k,i,j) )* dtlt
!-test	 vt3dc(k,i,j) = ( - sedim(ispc,ngrid)%v_sed_part(k,i,j) )* dtlt
      enddo
   enddo
enddo
!print*,'W=',wp(1:2,int(m2/2),int(m3/2)),wc(1:2,int(m2/2),int(m3/2))
!!!! (be sure the lines below are executed each timestep)!!!!
!- only necessary one time 
IF(ISPC==1) then
      do j = 1,m3
 	 jm = max(1,j-1)
 	 !jp = min(m3,j+1)
 	 do i = 1,m2
 	    im = max(1,i-1)
 	    !ip = min(m2,i+1)
 	    !rtgti = 1. / rtgt(i,j)
 	    do k = 1,m1-1
 	       vt3dp(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
 		  + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
 		  + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
 		  + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)
	       vt3dk(k,i,j) = dzt(k) / dn0(k,i,j)
 	    enddo

 	 enddo
       enddo
endif


!-  then the fluxes are re-evaluated
do j = 1,m3
   jm = max(1,j-1)
   !jp = min(m3,j+1)
   do i = 1,m2
      im = max(1,i-1)
      !ip = min(m2,i+1)
      rtgti = 1. / rtgt(i,j)
      do k = 1,m1-1
!-orig way
!        vt3dc(k,i,j) = ((vt3da(k,i,j) + vt3da(k+1,i,j)  &
!           + vt3da(k,im,j) + vt3da(k+1,im,j)) * f13t(i,j)  &
!           + (vt3db(k,i,j) + vt3db(k+1,i,j) + vt3db(k,i,jm)  &
!           + vt3db(k+1,i,jm)) * f23t(i,j)) * hw4(k)  &
!           + vt3dc(k,i,j) * rtgti
!opt way
         vt3dc(k,i,j) = vt3dp(k,i,j) &
            + vt3dc(k,i,j) * rtgti


	 vt3df(k,i,j) = .5 * vt3dc(k,i,j) * dzm(k)
	 !vt3dk(k,i,j) = dzt(k) / dn0(k,i,j)

      enddo

   enddo
enddo
!do k = 1,m1-1
!   vctr1(k) = (zt(k+1) - zm(k)) * dzm(k)
!   vctr2(k) =  (zm(k) - zt(k)) * dzm(k)
!enddo
!print*,'W=',wp(1:2,int(m2/2),int(m3/2)),wc(1:2,int(m2/2),int(m3/2))

!convert velocity components * dtlt (VT3DA, VT3DB, VT3DC)
! into mass fluxes times dtlt.

do j = 1,m3
   do i = 1,m2
      !c1 = fmapui(i,j) * rtgu(i,j)
      !c2 = fmapvi(i,j) * rtgv(i,j)
      do k = 1,m1-1
         vt3dc(k,i,j) = vt3dc(k,i,j) * .5  &
            * (dn0(k,i,j) + dn0(k+1,i,j))      !air_dens3d ????
      enddo
   enddo
enddo
end subroutine fa_preptc_with_sedim
!--------------------------------------------------------------------------

END MODULE module_dry_dep
