      subroutine vumat(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
      character*80 cmname
C      
      integer readflag
      real*4  edgedim(350000,3)
      common /integerbuf/readflag 
      common / realbuf/edgedim

        
        IF ( cmname(1:8).eq.'MATFIBRE' ) THEN
C		
            CALL  vumatfibre( 
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
	 
		ELSEIF ( cmname(1:6).eq.'MATCOH' ) THEN
 
            CALL vumatcoh(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
     
		ENDIF
C         
      RETURN
      END       

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                    C
C  subroutine: vumatfibre                                            C
C  function: Fibre failure activation using maximum stress criteria  C
C            and evolution                                           C 
C                                                                    C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC 
C SDV1  - Strain 11
C SDV2  - Strain 22
C SDV3  - Strain 33
C SDV4  - Strain 12, engineering shear strain gamma
C SDV5  - Original strain of the unloading part 
C SDV6  - The sign of the increment of strain gamma_12 
C SDV7  - The flow of loading; Loading: flow=1; Unloading: flow=-1
C SDV8  - Origin strain of the loading stage (origin of the nonlinear shear)
C SDV9  - Maximum strain experienced during mono loading stage (the sign of the stress does not change)
C SDV10 - Maximum stress ever experienced during one cycle
C SDV11 - Sign change of shear stress: sign_change
C SDV12 - Shear failure index under tension SHR_FIt
C SDV13 - Shear failure index under compression SHR_FIc
C SDV14 - Shear damage variable D_SHR
C SDV15 - Element deletion flag
C SDV16 - Stress 12 at shear damage initiation under tension: sig12_f0_T
C SDV17 - Strain 12 at shear damage initiation under tension: eps12_f0_T
C SDV18 - Strain 12 at shear damage completion under tension: eps12_ff_T
C SDV19 - Charcteristic length L
C SDV20 - Shear failure flag
C SDV21 - Shear damage variable under tension, SHR_DPlus
C SDV22 - Shear damage variable under compression, SHR_DMinus
C SDV23 - Stress 12 at shear damage initiation under compression: sig12_f0_C
C SDV24 - Strain 12 at shear damage initiation under compression: eps12_f0_C
C SDV25 - Strain 12 at shear damage completion under compression: eps12_ff_C
C SDV26 - Shear failure flag under tension
C SDV27 - Shear failure flag under compression 

		 subroutine vumatfibre(
C Read only (unmodifiable)variables -
     1  nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     2  stepTime, totalTime, dt, cmname, coordMp, charLength,
     3  props, density, strainInc, relSpinInc,
     4  tempOld, stretchOld, defgradOld, fieldOld,
     5  stressOld, stateOld, enerInternOld, enerInelasOld,
     6  tempNew, stretchNew, defgradNew, fieldNew,
C Write only (modifiable) variables -
     7  stressNew, stateNew, enerInternNew, enerInelasNew )
C
      include 'vaba_param.inc'
C
      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     8  defgradNew(nblock,ndir+nshr+nshr),
     9  fieldNew(nblock,nfieldv),
     1  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     2  enerInternNew(nblock), enerInelasNew(nblock)
C
       parameter(zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, 
     * third = one/three, half = 0.5d0, twothird = two/three, tol=0.1d0,
     * threehalf = 1.5d0)              
C     
       character*80 cmname,cpname
       character*256 outdir,fullpath
       integer intnum,locnum,jrcd,lenoutdir
       double precision E11,E22,E33,v12,v21,v13,v31,v23,v32,G12,G13,G23
       double precision C11 ,C12 ,C13 ,C14 ,C21 ,C22,
     *   C23 , C24,C31 ,C32 ,C33 ,C34,
     *   C41 , C42 ,C43 ,C44
	   double precision n,m,delta
	   double precision Xt,Xc,Yc,Yt
	   double precision Gft,Gfc
	   
	!  Nonlinear shear variables 
       double precision NLS, A, B, G_ply, GSHRt, GSHRc
	   double precision sgn_d_strain, flow, s12_max, sign_change
	   double precision SHR_FIt, SHR_FIc, SHR_FIt_old, SHR_FIc_old
	   double precision SHR_FIt_max, SHR_FIc_max
	   double precision D_SHR, SHR_DPlus, SHR_DMinus
	   double precision sig12_f0_T, eps12_f0_T, eps12_ff_T
	   double precision sig12_f0_C, eps12_f0_C, eps12_ff_C
	   
       double precision FIt, FIc, FIt_old, FIc_old, FIt_max, FIc_max
       double precision sig_f0, eps_f0, eps_f
	   double precision sig_f0_T, eps_f0_T, eps_ff_T, eps_f_T
	   double precision sig_f0_C, eps_f0_C, eps_ff_C, eps_f_C
	   double precision D, d1Plus, d1Minus
       real*8 Stress(ndir+nshr), Strain(ndir+nshr)
	   real*8 NLS
C	   
	   integer i,j
C
C INITIALISATION: MATERIAL CARD PARAMETERS
C =========================================
C
      E11 = props(1) 
      E22 = props(2) 
      E33 = props(3) 
      v12 = props(4) 
      v13 = props(5) 
      v23 = props(6) 
      G12 = props(7) 
      G13 = props(8) 
      G23 = props(9) 
	  NLS = props(10) ! Nonlinear shear flag
	  A   = props(11)
	  B   = props(12)
	  G_ply = props(13) ! Shear modulus of single ply (important for nonlinear shear model!)
	  GSHRt = props(14)
	  GSHRc = props(15) 

C
C AXISYMMETRIC MODEL STIFFNESS MATRIX 
C ====================================
C
      v21 = v12*E22/E11
      v31 = v13*E33/E11
      v32 = v12
	  
      n = E11/E22
	  m = G12/E22
C
	  delta = (one + v13)*(one - v13 - two*n*v21**2)
C
	  C11 = (E22/delta) * n * (one - n*v21**2)
	  C12 = (E22/delta) * n * v21 * (one + v13)
	  C13 = (E22/delta) * n * (v13 + n* v21**2)
	  C14 = zero
	  C21 = C12
	  C22 = (E22/delta) * (one - v13**2)
	  C23 = (E22/delta) * n * v21 * (one + v13)
	  C24 = zero
	  C31 = C13
	  C32 = C23
	  C33 = (E22/delta) * n * (one - n*v21**2)
	  C34 = zero
	  C41 = zero
	  C42 = zero
	  C43 = zero
	  C44 = G12 
	  
C	  
C START CONTINNUM DAMAGE MODEL
C ============================
C
C INITIAL ELASTIC STEP FOR ABAQUS TESTS
C -------------------------------------	  
C        
      IF (totalTime.le.dt) THEN
      
          DO i = 1, nblock 
		  
		      stressNew(i,1)= stressOld(i,1)+ C11*straininc(i,1)
     *           + C12*straininc(i,2) + C13*straininc(i,3)
        
		      stressNew(i,2)= stressOld(i,2)+ C21*straininc(i,1)
     *           + C22*straininc(i,2) + C23*straininc(i,3) 

              stressNew(i,3)=stressOld(i,3) + C31*straininc(i,1)
     *           + C32*straininc(i,2) + C33*straininc(i,3)
	 
              stressNew(i,4) = stressOld(i,4) + two*C44*straininc(i,4)  
			  
			  sign_change = one
			  stateNew(i,11) = sign_change
			  
			  D_SHR = stateNew(i,14)
			  stateNew(i,20) = zero    ! SDV20 - Shear failure flag: 0-elastic, 1-softening
			  stateNew(i,26) = zero    ! SDV26 - Shear failure flag under tension
			  stateNew(i,27) = zero    ! SDV27 - Shear failure flag under compression
			  stateNew(i,15) = one     ! SDV15 - Element delete flag; 0-delete
			  stateNew(i,19） = charLength(i)
			  
          ENDDO
C
C IF NOT INITIAL STEP, CAL. STRESSES ACCRODING TO CDM MODEL
C =========================================================
C      
      ELSE   
C         
          DO i = 1, nblock
C			  
		      Do j = 1, nstatev
			      stateNew(i,j)=stateOld(i,j)
		      ENDDO
			  
			  D_SHR = stateNew(i,14)
			  
			  IF ( stateNew(i,15) .eq. one) THEN  ! If the element is not deleted
			  
				  stateNew(i,1)= stateOld(i,1)+ straininc(i,1)
				  stateNew(i,2)= stateOld(i,2)+ straininc(i,2)
				  stateNew(i,3)= stateOld(i,3)+ straininc(i,3)
				  stateNew(i,4)= stateOld(i,4)+ two*straininc(i,4)  
C       
				  stressNew(i,1)= C11*stateNew(i,1)+C12*stateNew(i,2)
     *                  +C13*stateNew(i,3)
				  stressNew(i,2)= C21*stateNew(i,1)+C22*stateNew(i,2)
     *                 + C23*stateNew(i,3)
				  stressNew(i,3)= C31*stateNew(i,1)+C32*stateNew(i,2)
     *                 + C33*stateNew(i,3)  

C
C CALCULATE ELASTIC SHEAR STRESS USING NONLINEAR SHEAR MODEL
C =========================================================
C  
				  IF (NLS .EQ. ZERO) THEN
				  
					  stressNew(i,4) = C44*stateNew(i,4)
				  
				  ELSE IF (NLS. EQ. ONE) THEN 
C		  
					  sgn_d_strain = Sign(one, strainInc(i,4))
					  stateNew(i,6) = sgn_d_strain
				  
					  flow = Sign(one, stateNew(i,6)*stressOld(i,4))
					  stateNew(i,7) = flow
				  
					  IF (stateNew(i,7).EQ.one) THEN                              ! flow=1: loading; flow=-1: unloading					  
					  
						  IF ((stateNew(i,7) .NE. stateOld(i,7)) .AND.
     *					     ((stateNew(i,11) .EQ. one))) THEN         ! Record the shear strain at the point that flow changes
						  
							  stateNew(i,8) = stateOld(i,4)              ! SDV8: Origin of the nonlinear part/loading part					       						  
						  
						  ENDIF					  	
					   
						  IF ( stateNew(i,6)*stateNew(i,4) .GE. 
     *					       stateNew(i,6)*stateNew(i,9) ) THEN ! SDV9: max strain experienced
						  
							  stateNew(i,9) = stateNew(i,4)
							  stressNew(i,4) = stateNew(i,6)*abs(A*(one-
     *					     	exp(-B*abs(stateNew(i,4)-stateNew(i,8)))))
	 
							  s12_max = stressNew(i,4) 
							  stateNew(i,10) = s12_max
						  
						  ELSE
					  
							  stateNew(i,9) = stateOld(i,9)
							  stressNew(i,4) = stateNew(i,6)*abs(
     *					      abs(G_ply*(abs(stateNew(i,9)- stateNew(i,4))))
     *	                      -abs(stateNew(i,10)))

						  ENDIF
					  		  
					  ELSE IF (stateNew(i,7).NE.one) THEN
					  
						  IF (stateNew(i,7) .NE. stateOld(i,7)) THEN
							  stateNew(i,5) = stateNew(i,9)              ! SDV5: Origin of the unloading part
						  ELSE
							  stateNew(i,5) = stateOld(i,5)
						  ENDIF                                          !      /Maximum strain in the loading part
					  
						  stressNew(i,4) = stateNew(i,6)*(
     *					      G_ply*(abs(stateNew(i,4)- stateNew(i,5)))
     *	                     -abs(stateNew(i,10)))  

					  ENDIF
				  
					  IF (Sign(one,stressNew(i,4)).EQ.Sign(one,stressOld(i,4))) THEN
						  sign_change = zero 
					  ELSE 
						  sign_change = one
					  ENDIF
					  stateNew(i,11) = sign_change
				  
				  ENDIF ! for nls flag
			  
				  Stress(1:4) = stressNew(i,1:4) ! Effective stresses in the undamaged material
			  
C
C CALCULATE DAMAGED SHEAR STRESS
C =========================================================
C  			  
			      IF (stateNew(i,20) .eq. zero ) THEN  ! In the elastic region, calculate failure indices
					  SHR_FIt = zero
					  SHR_FIc = zero

C
C CALCULATE FAILURE INDICES
C -------------------------
					  
					  IF (Stress(4) >= zero ) THEN
						  SHR_FIt = Stress(4)/A
					  ELSEIF (Stress(4) < zero ) THEN
						  SHR_FIc = -Stress(4)/A
					  ENDIF
					  
					  SHR_FIt_old = stateOld(i,12)
					  SHR_FIc_old = stateOld(i,13)
					  
					  SHR_FIt_max = max(SHR_FIt, SHR_FIt_old)
					  SHR_FIc_max = max(SHR_FIc, SHR_FIc_old)
					  
					  stateNew(i,12) = SHR_FIt_max
					  stateNew(i,13) = SHR_FIc_max

C
C CALCULATE DAMAGE VARIABLE 
C -------------------------	
C 
					  SHR_DPlus = zero
					  SHR_DMinus = zero
					  
					  IF ( stateNew(i,12) >= one .or. stateNew(i,13) >= one ) THEN
C					      
						  stateNew(i,20) = one ! Shear soften				
C						  
					      IF ( stateNew(i,4) > zero ) THEN
						  
							  stateNew(i,26) = one    ! Shear failure flag under tension
							  
							  sig12_f0_T = Stress(4)
							  stateNew(i,16) = sig12_f0_T
							  
							  eps12_f0_T = stateNew(i,4)
							  stateNew(i,17) = eps12_f0_T
							  
						      eps12_ff_T = 2.d0*GSHRt/(A*charLength(i)) ! Strain11 at final damage
						      stateNew(i,18) = eps12_ff_T			
							  
						  ELSEIF ( stateNew(i,4) < zero ) THEN
						  
							  stateNew(i,27) = one    ! Shear failure flag under compression
							  
							  sig12_f0_C = Stress(4)
							  stateNew(i,23) = sig12_f0_C
							  
							  eps12_f0_C = stateNew(i,4)
							  stateNew(i,24) = eps12_f0_C
							  
							  eps12_ff_C = -2.d0*GSHRc/(A*charLength(i)) ! eps12_ff_C is minus
							  stateNew(i,25) = eps12_ff_C			
							  
						  ENDIF
C						  
					  ENDIF
C					  
				  ELSEIF ( stateNew(i,20) .eq. one ) THEN ! In the soften region, calculate damage variables
C
					  eps12_f = stateNew(i,4)
C
					  IF ( stateNew(i,4) .GT. zero ) THEN
C					  
						  IF ( stateNew(i,26) .EQ. zero ) THEN

						      SHR_FIt = Stress(4)/A
							  SHR_FIt_old = stateOld(i,12)
							  SHR_FIt_max = max(SHR_FIt, SHR_FIt_old)
							  stateNew(i,12) = SHR_FIt_max
C
							  IF ( stateNew(i,12) .GT. one ) THEN  !If tension damage initiates
							  
							      stateNew(i,26) = one				
								  
								  sig12_f0_T = Stress(4)
								  stateNew(i,16) = sig12_f0_T
								  
								  eps12_f0_T = stateNew(i,4)
								  stateNew(i,17) = eps12_f0_T
								  
								  eps12_ff_T = 2.d0*GSHRt/(A*charLength(i)) ! Strain11 at final damage
								  stateNew(i,18) = eps12_ff_T
								  
							  ENDIF
							  
						  ELSEIF ( stateNew(i,26) .EQ. one ) THEN
						  
							  eps12_f_T  = stateNew(i,4)
							  eps12_f0_T = stateNew(i,17)
							  eps12_ff_T = stateNew(i,18)
							  
							  SHR_DPlus = (eps12_ff_T * (eps12_f_T - eps12_f0_T))/
     *						           (eps12_f_T * (eps12_ff_T - eps12_f0_T))
							  SHR_DPlus = min(one, max(zero,SHR_DPlus))
							  SHR_DPlus = max(stateNew(i,21), SHR_DPlus)
							  stateNew(i,21) = SHR_DPlus
						  
						  ENDIF
C						  
					  ELSEIF ( stateNew(i,4) .LT. zero ) THEN
					  
						  IF ( stateNew(i,27) .EQ. zero ) THEN
						  
						      SHR_FIc = abs(Stress(4)/A)
							  SHR_FIc_old = stateOld(i,13)
							  SHR_FIc_max = max(SHR_FIc, SHR_FIc_old)
							  stateNew(i,13) = SHR_FIc_max

							  IF ( SHR_FIc_max .GT. one ) THEN
								  stateNew(i,27) = one
								  
								  sig12_f0_C = Stress(4)
								  stateNew(i,23) = sig12_f0_C
								  
								  eps12_f0_C = stateNew(i,4)
								  stateNew(i,24) = eps12_f0_C
								  
								  eps12_ff_C = -2.d0*Gfc/(Xc*charLength(i)) ! eps12_ff_C is minus
								  stateNew(i,25) = eps12_ff_C
								  
							  ENDIF
C							  
						  ELSEIF ( stateNew(i,27) .EQ. one ) THEN
							  eps12_f_C = stateNew(i,4)
							  eps12_f0_C = stateNew(i,24)
							  eps12_ff_C = stateNew(i,25)
							  SHR_DMinus = (eps12_ff_C * (eps12_f_C - eps12_f0_C))/
     *			        			    (eps12_f_C * (eps12_ff_C - eps12_f0_C))
							  SHR_DMinus = min(one, max(zero,SHR_DMinus))
							  SHR_DMinus = max(stateNew(i,22), SHR_DMinus)
							  stateNew(i,22) = SHR_DMinus
						  ENDIF
C
					  ENDIF
C
C CRACK CLOSURE UNDER LOAD REVERSAL
C ---------------------------------					  
C				  
					  D_SHR = max(SHR_DPlus, SHR_DMinus)
					  D_SHR = min(one, max(zero, D_SHR))
					  D_SHR = max(stateNew(i,14), D_SHR)
					  stateNew(i,14) = D_SHR
					  
C 
C STRESS DEGRADATION
C ------------------
					  Stress(4) = (1-D_SHR)*Stress(4)
				  
				  ENDIF
C
C TRANSFER STRESS VECTOR BACK TO ABAQUS
C -------------------------------------
				  stressNew(i,1:4) = Stress(1:4)
				  
			  ENDIF ! FOR IF OF ELEMENT DELETION FLAG (SDV15)
	      ENDDO
C	  
	  ENDIF
C	
	  RETURN 	  
	  END

!-------------------------------------------------------------------!
!      Subroutine vumatcoh:                                         !
!-------------------------------------------------------------------!       

#include 'Abaqus_Definitions.f'
#include 'Fatigue_Globals.f'
#include 'Mesh_Utils.f'
#include 'Fatigue_ANN.f'
#include 'Fatigue.f'
#include 'CZM.f'

      !> VUMAT Main entry point from Abaqus
      subroutine vumatcoh(
     *     jblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew )

        use CZM, only: vumat_cohesive_fatigue
        include 'vaba_param.inc'

        dimension jblock(*), props(nprops),density(*), coordMp(*),
     *     charLength(*), strainInc(*), relSpinInc(*), tempOld(*),
     *     stretchOld(*), defgradOld(*), fieldOld(*), stressOld(*),
     *     stateOld(*), enerInternOld(*),enerInelasOld(*),
     *     tempNew(*), stretchNew(*), defgradNew(*), fieldNew(*),
     *     stressNew(*), stateNew(*), enerInternNew(*), enerInelasNew(*)

        character*80 cmname

        call vumat_cohesive_fatigue ( jblock(1),
     *     ndir, nshr, nstatev, nfieldv, nprops, lanneal,
     *     stepTime, totalTime, dt, cmname, coordMp, charLength,
     *     props, density, strainInc, relSpinInc,
     *     tempOld, stretchOld, defgradOld, fieldOld,
     *     stressOld, stateOld, enerInternOld, enerInelasOld,
     *     tempNew, stretchNew, defgradNew, fieldNew,
     *     stressNew, stateNew, enerInternNew, enerInelasNew,
     *     jblock(2), jblock(3), jblock(4), jblock(5))

      end subroutine vumatcoh

      
      !> Vexternaldb entry from Abaqus
      subroutine vexternaldb(lOp, i_Array, niArray, r_Array, nrArray)

        use Abaqus_Definitions, only: j_int_StartAnalysis
        use CZM, only: CZM_initialisation
        include 'vaba_param.inc'

        dimension i_Array(niArray), r_Array(nrArray)

        ! Initialisation at start of the analysis
        if (lOp .eq. j_int_StartAnalysis) then
          call CZM_initialisation()
        endif

      end subroutine vexternaldb			
