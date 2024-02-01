!DIR$ FREEFORM

!> Module for cohesive zone modelling
module CZM
  use Abaqus_Definitions, only: wp=>abaqus_real_kind
  implicit none

contains

!> Static/fatigue cohesive law vumat
  subroutine vumat_cohesive_fatigue ( &
         nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
         stepTime, totalTime, dt, cmname, coordMp, charLength, &
         props, density, strainInc, relSpinInc, &
         tempOld, stretchOld, defgradOld, fieldOld, &
         stressOld, hsvOld, enerInternOld, enerInelasOld, &
         tempNew, stretchNew, defgradNew, fieldNew, &
         stressNew, hsvNew, enerInternNew, enerInelasNew, &
         nMatPoint, nLayer, nSecPoint, nElement)

    use Fatigue_Globals
    use Fatigue, only: calculateParisLawParams, zhangsMethod_DamageInitiated, &
                       mayTaoMethod_DamageInitiated, checkFatigueSoftening, &
                       fatigueSoftening, storeNeighboursAndDimensions
    
    integer, intent(in) :: nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal
    real(wp), intent(in) :: stepTime, totalTime, dt
    character(*), intent(in) :: cmname
    real(wp), intent(in) :: coordMp(nblock,*), charLength(nblock), props(nprops), &
                            density(nblock), strainInc(nblock,ndir+nshr), &
                            relSpinInc(nblock,nshr), tempOld(nblock), &
                            stretchOld(nblock,ndir+nshr), &
                            defgradOld(nblock,ndir+nshr+nshr), &
                            fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr), &
                            hsvOld(nblock,nstatev), enerInternOld(nblock), &
                            enerInelasOld(nblock), tempNew(nblock), &
                            stretchNew(nblock,ndir+nshr), &
                            defgradNew(nblock,ndir+nshr+nshr), &
                            fieldNew(nblock,nfieldv)
    real(wp), intent(out) :: stressNew(nblock,ndir+nshr), hsvNew(nblock,nstatev)
    real(wp), intent(inout) :: enerInternNew(nblock), enerInelasNew(nblock)
                 
    integer, intent(in) :: nElement(nblock),nMatPoint(nblock), &
                           nLayer(nblock),nSecPoint(nblock)


    character*80 cpname
    real(wp), parameter :: zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, &
              third = one/three, four = 4.d0,  half = 0.5d0, sf=1.e-20

    real*8 GIc,GIIc,GIIc_n,alpha,FIlimit,FIIlimit,FIIlimit_n,eta_F, &
       eta_G,EI,EII,fatinit,frq,CI,CII,expI,expII,a_SNI,a_SNII,R,FCL, &
       FI,FII,Fmixed,uI,uII,u,cosI,cosII,u0,FIyield,FIIyield,umax, &
       Fyield,Gc,Gc_IP,Ds_max,uratio,fatcycle,D_rate,D_inc, &
       cMM,expMM,dadN,failure_time,Gmax,area, &
       L,W,T,ulocal(3),score(4),l1x,l1y,l1z,l2x,l2y,l2z,l1,l2,e1,e2, &
       em,temp,tempI,tempII,sig_SNI, &
       sig_SNII,FIType,FLawType,FCLType, &
       MMPPType,Gratio,FTType,FDRType,IfNoTrack,CLType,BK,TTCELaw,cm, &
       expm,Di_I,Di_II,C,b0I,h,Alpha_0,IfFatigue,FIlimit_0, &
       FIIlimit_0,ElemType,IPnum4Tip, Ds, Dtotal, &
       IfElemTip,E11,E33,G13,v13,u_sum,Phi(4),Z !ANN parameters
    integer currentelem(9),locnum,threadID,tempcount,locIPnum, &
            ElemNum1CubR,i_33,i_23,i_13,i_12,i_22,NumIPs, i, j, jrcd
    character*256 outdir,filepath,jobname
    integer lenoutdir,lenjobname,keyfile


    !**** Start of the main programme from here

    ! Read material properties related to static law
    GIc=props(1)
    GIIc=props(2)
    FIlimit_0=props(3)
    FIIlimit_0=props(4)
    EI=props(5)
    EII=props(6)
    alpha=props(7)
    TTCELaw=props(8)
    eta_F=props(9)
    eta_G=props(10)
    
    if (FIIlimit_0.lt.sf) then
      FIIlimit_0=FIlimit_0
    end if
    
    if (EII.lt.sf) then
      EII=EI
    end if
  
    ! Read material properties related to fatigue law
    IfFatigue=props(11)
    fatinit=props(12)
    frq=props(13)
    R=props(14)
    FIType=props(15)
    a_SNI=props(16)
    a_SNII=props(17)
    FLawType=props(18)
    if ((FLawType.ne.one).and.(FLawType.ne.two).and. &
       (FLawType.ne.three).and.(FLawType.ne.four)) then
    FLawType=zero
    end if
    CI=props(19)
    expI=props(20)
    CII=props(21)
    expII=props(22)
    ! when Allegris's law is used for fatigue crack growth rate
    if (FLawType.eq.zero) then
      C=CI; b0I=expI; h=CII
    elseif (FLawType.eq.four) then ! when ANN used
      E11=CI; E33=expI; G13=CII; v13=expII
    end if
    MMPPType=props(23)
    cm=props(24)
    expm=props(25)
    FDRType=props(26)
    IfNoTrack=props(27)
    FCLType=props(28)
    FTType=props(29)
    ElemType=props(30)
    if (ElemType.eq.one) then 
      NumIPs=1 ! C3D8R
      i_33=3
      i_23=5
      i_13=6
    else
      NumIPs=4 ! COH3D8
      i_22=1 
	  i_12=2 
!      i_23=3 
	
      
    end if
    IfElemTip=props(31)
    IPnum4Tip=props(32)
    if ((IPnum4Tip.ne.one).and.(IPnum4Tip.ne.two).and. &
      (IPnum4Tip.ne.three)) then
    IPnum4Tip=four
  end if

  do i = 1, nblock

    ! Save last dt's history variable into current dt's
    do j=1,nstatev
      hsvNew(i,j)=hsvOld(i,j)
    enddo

    T = charLength(i)

    if (totalTime.le.dt) then

      ! Save constitutive thickness by assuming it is COH3D8 element
      hsvNew(i,44) = T

    end if

    if (IfFatigue.ne.zero) then

      if (totalTime.le.dt) then

        ! Retrieve original element ID
        call vgetpartinfo(nElement(i),1,cpname,locnum,jrcd)
        hsvNew(i,42)=locnum

        ! Save local IP number by its call sequence by the VUMAT
        GV4Elem(locnum,7)=GV4Elem(locnum,7)+1
        hsvNew(i,41)=GV4Elem(locnum,7)
        if (GV4Elem(locnum,7).ge.NumIPs) then
          GV4Elem(locnum,7)=zero
        end if

        ! Store neighbours and geometric dimensions into history variables
        call storeNeighboursAndDimensions(hsvNew, nblock, nstatev, &
                      i, locnum, NumIPs, T, ElemDimens, NodeCoords, &
                      IPEdgNbrElem, CohElemNum, CohElem, &
                      IP1NbrElem, IP2NbrElem, IP3NbrElem, IP4NbrElem)


        ! if C3D8R elem used, save geometric thickness as constitutive thickness
        if (ElemType.eq.one)then
          T = max(T,sf)
          hsvNew(i,44) = T
        end if

      else ! if not first increment

        T = hsvNew(i,44)  ! Check out element constitutive thickness
        locnum = hsvNew(i,42)  ! Check out original element ID

      end if ! (totalTime.le.dt)

    end if ! (IfFatigue.ne.zero)

	i_22 = 1
	i_12 = 2
    ! mode I and II strains and displacements
    hsvNew(i,16) = hsvOld(i,16) + straininc(i,i_22)  ! mode I strain
    hsvNew(i,17) = hsvOld(i,17) + straininc(i,i_12)  ! mode II 13 shear strain

    ulocal(1) = T*hsvNew(i,16)       ! mode I displacement
    ulocal(2) = T*two*hsvNew(i,17)   ! mode II 13 shear displacement

    hsvNew(i,1)=ulocal(1)
    hsvNew(i,2)=ulocal(2)

    if (totalTime.le.dt) then

      ! Purely elastic definition at the beginning to avoid unstable solution
      stressNew(i,i_22)=  EI*ulocal(1)
      stressNew(i,i_12)=  EII*ulocal(2)

      ! Store original mode I and II strengths
      hsvNew(i,20) = FIlimit_0
      hsvNew(i,10) = FIIlimit_0

    else ! if not first increment

      ! Static and fatigue cohesive law starts from the second dt

      ! Displacments
      uI=max(zero,ulocal(1))
      hsvNew(i,4)=uI
      uII=dsqrt(ulocal(2)**2) ! resultant shear disp. remove ulocal(3)**2!
      hsvNew(i,5)=uII
      u=dsqrt(uI**2+uII**2) !  mixed-mode disp.

      ! Direction cosines
      if(u.eq.zero)then
        cosI=zero
      else
        cosI=uI/u
      end if
      cosII=dsqrt(dabs(one-cosI**2))

      ! Check out mode I strength, possibly degraded due to fatigue initiation
      FIlimit=hsvNew(i,20)
      call calcMIIStrengthAndFractureToughness(hsvNew, nblock, &
                        nstatev, &
                        FIIlimit_n, FIIlimit_0, FIType, totalTime, fatinit, &
                        FIILimit, TTCELaw, GIIc_n, GIIc, eta_G, eta_F, EI, &
                        ulocal, i, IfFatigue)

      ! Calculating damage initiation displacemnt
      u0= dsqrt(one/((EI*cosI/FIlimit)**2+(EII*cosII/FIIlimit)**2))
	  !WRITE(*,*) 'u0=', u0

      ! Mode I, II in-situ damage onset strengths
      FIyield=EI*u0*cosI
      FIIyield=EII*u0*cosII

      ! Mode-mixity ratio by GII/(GI+GII) at full failure
      Gratio=EII*CosII*COSII/(EI*cosI*cosI+EII*cosII*cosII)

      ! Mixed-mode failure displacement by Power law or B-K law
      if (alpha.gt.zero) then  ! power law
		! umax: u^f displacement at full failure
        umax=one/((half*FIyield*cosI/GIc)**alpha + &
                      (half*FIIyield*cosII/GIIc_n)**alpha)
        umax=umax**(1/alpha)
      else   ! B-K law
        Gc_IP=GIc+(GIIc_n-GIc)*(Gratio**(-alpha))
        umax=two*Gc_IP/(FIyield*cosI+FIIyield*cosII)
      end if	
		
      !  Mixed mode yield stress and ERR for complete decohesion
      Fyield= dsqrt((FIyield**2)+(FIIyield**2))
      Gc_IP=half*(FIyield*cosI+FIIyield*cosII)*umax    ! compute Gc=GI+GII
      hsvNew(i,25)=Gc_IP

      ! Static damage component Ds
      Ds = (u-u0)/(umax-u0)
      Ds_max=max(hsvNew(i,15),Ds)  ! unloading not affecting Ds
      hsvNew(i,15)=Ds_max    ! hsv 15 Ds_max

      ! If fatigue, compute Paris parameters, average G,Gc,c,m and Ds if needed
      if (IfFatigue.ne.zero) then ! Only execute for fatigue
        ! Paris law parameters for the current IP If ANN not used
        call calculateParisLawParams(FLawType, MMPPType, cMM, cI, &
                      cm, Gratio, cII, expMM, expI, expII, expm, C, &
                      Gc_IP, GIc, GIIc_n, h, b0I, R, &
                      hsvNew, nblock, nstatev, IfElemTip, &
                      locnum, NumIPs, i)

      end if

      ! Global elapsed number of cycles
      if ((IfFatigue.ne.zero).and.(totalTime.ge.fatinit))then
        fatcycle=hsvNew(i,33) + frq*dt  ! number of elapsed cycles
        hsvNew(i,33)=fatcycle
      end if

      !****Traction update depending on which stage the elemnent/IP is in***
      !*Distinguish between elastic loading, already full failure, full failure
      !*in the current dt, static sofening or fatigue softening

      if ((Ds_max.le.zero).and.(hsvNew(i,14).le.one)) then

        !>>Elastic loading

        uratio=u/(u0+sf)         !Jiang-EQ(16) with D=0
        FI=FIyield*uratio        !Jiang-EQ(14)
        FII=FIIyield*uratio      !Jiang-EQ(15)
        Fmixed=Fyield*uratio     
        hsvNew(i,14)=1 ! hsv 14 <= 1: elastic loading

        !Check if damage is initiated due to fatigue at the elastic stage
        if ((IfFatigue.ne.zero).and.(totalTime.gt.fatinit)) then

          if (FIType.eq.one) then

            call zhangsMethod_DamageInitiated(hsvNew, nblock, nstatev, &
                              EI, u, cosI, FILimit, a_SNI, a_SNII, frq, dt, EII, cosII, &
                              FIILimit, uratio, i, locnum)

          elseif (FIType.eq.two) then ! Mukhopadhyay's method

            hsvNew(i,20)=FIlimit_0*(1.d0-a_SNI*log10(fatcycle))
            hsvNew(i,10)=FIIlimit_n*(1.d0-a_SNII*log10(fatcycle))

          else
            call mayTaoMethod_DamageInitiated(hsvNew, nblock, nstatev, &
                              EI, u, cosI, FILimit, a_SNI, a_SNII, frq, dt, EII, cosII, &
                              FIILimit, uratio, i, locnum)
          end if ! (FIType.eq.one)

        end if ! end of 'if ((IfFatigue.ne.zero).and.(totalTime.gt.fatinit))'

      elseif (hsvNew(i,14).eq.four) then

        !>>Already fully failed

        FI=zero
        FII=zero
        Fmixed=zero
        hsvNew(i,14)=4
        hsvNew(i,28)=zero

      elseif ((Ds_max.ge.one).and.(hsvNew(i,14).ne.four)) then

        !>>Fully failed in the current step

        FI=zero
        FII=zero
        Fmixed=zero
        hsvNew(i,14)=4
        hsvNew(i,28)=zero
        hsvNew(i,32)=totalTime   ! recording failure time
        if (IfFatigue.ne.zero) then ! Only execute for fatigue
          GV4Elem(locnum,1)=GV4Elem(locnum,1)+one ! update the no. of failed IPs
        end if


      else !>>Static or Fatigue softening

        ! Check if fatigue softening
        if ((IfFatigue.ne.zero).and.(totalTime.ge.fatinit)) then

          ! If not previously checked out
          if (hsvNew(i,14).eq.2.0d0) then
            call checkFatigueSoftening(hsvNew, nblock, nstatev, locnum, i, &
                                                IfNoTrack, IfElemTip, IPnum4Tip)
          end if

        end if ! end of 'if ((IfFatigue.ne.zero).and.(totalTime.ge.fatinit))'

        if (hsvNew(i,14).eq.three) then

          !>Fatigue softening stage

          call fatigueSoftening(hsvNew, nblock, nstatev, i, FDRType, &
                        IfElemTip, locnum, Gmax, Gc, R, FLawType, &
                        cMM, expMM, dadN, E11, E33, G13, v13, FIlimit, &
                        FIIlimit, GIC, GIIC_n, NumIPs, Z, &
                        FCLType, ElemType, FCL, frq, dt, Ds_max, &
                        D_inc, Dtotal, FI, FII, Fmixed, FTType, FIyield, &
                        FIIyield, Fyield, totalTime, u0, umax, u)

        else

          !>Static softening stage

          call staticSoftening(hsvNew, nblock, nstatev, FI, FII, &
                            FIyield, FIIyield, Fmixed, Fyield, Ds, Ds_max, &
                            u, u0, umax, IfFatigue, totalTime, fatinit, &
                            FILimit, FIILimit, FIlimit_0, FIIlimit_0, &
                            a_SNI, a_SNII, cosI, cosII, frq, dt, i, FIType)

        end if   ! end of if(hsvNew(i,14).eq.three

      end if ! end of if((Ds_max.le.zero).and.(hsvNew(i,14).eq.zero))then

      ! Assign tractions to stressNew components
      if(hsvNew(i,1).gt.zero) then
        stressNew(i,i_22)=FI
      else
        stressNew(i,i_22)=EI*hsvNew(i,1)
      end if
      hsvNew(i,13)=stressNew(i,i_22)

      if(uII.eq.zero)then
        stressNew(i,i_12)=zero
      else
        stressNew(i,i_12)=FII*hsvNew(i,2)/(uII+sf)
      end if
      hsvNew(i,11)=stressNew(i,i_12)


      ! Assign zero to other two stress components (S11,S33)
      do j=1,ndir+nshr
        if ((j.ne.i_22).and.(j.ne.i_12)) then
          stressNew(i,j)=zero
        end if
      enddo

      ! Save energies and mode mixity
      if(hsvNew(i,1).gt.0.0)then
        if(hsvNew(i,1).gt.hsvNew(i,19))then
          hsvNew(i,6)=hsvNew(i,6)+ &
                            half*(hsvNew(i,22)+FI)*(hsvNew(i,1)-hsvNew(i,19))  ! GI
          hsvNew(i,22)=FI
        end if
      end if
      hsvNew(i,19)=max(hsvNew(i,19),hsvNew(i,1)) ! max of mode I disp over time

      if(hsvNew(i,5).gt.hsvNew(i,21))then
        hsvNew(i,7)=hsvNew(i,7)+ &
                        half*(hsvNew(i,23)+FII)*(hsvNew(i,5)-hsvNew(i,21))  ! GII
        hsvNew(i,23)=FII
      end if
      hsvNew(i,21)=max(hsvNew(i,21),hsvNew(i,5)) ! max of mode II disp over time

      hsvNew(i,8)=hsvNew(i,6)+hsvNew(i,7)        ! GI+GII

      if(hsvNew(i,14).gt.one)then
        hsvNew(i,9)=hsvNew(i,7)/(hsvNew(i,8)+sf) ! mode mixity
      end if

      !  Gratio calculation (G/Gc)
      if(hsvNew(i,25).eq.0.d0)then
        hsvNew(i,26)=0.d0
      else
        hsvNew(i,26)=hsvNew(i,8)/hsvNew(i,25)
      end if


      if (IfFatigue.ne.zero) then

        hsvNew(i,43)=GV4Elem(locnum,1)   ! Save the number of failed IPs

        ! Save ERR,Gc,c,m,Ds,failure flag, GIIC_n, phi to common blocks if fatigue
        GV4IP(locnum,5*hsvNew(i,41)-4)=hsvNew(i,8)   ! ERR
        GV4IP(locnum,5*hsvNew(i,41)-3)=hsvNew(i,25)  ! Gc
        GV4IP(locnum,5*hsvNew(i,41)-2)=hsvNew(i,34)  ! c
        GV4IP(locnum,5*hsvNew(i,41)-1)=hsvNew(i,35)  ! m
        GV4IP(locnum,5*hsvNew(i,41)-0)=hsvNew(i,15)  ! Ds
        GV4IP(locnum,hsvNew(i,41)+20)=hsvNew(i,14)   ! failure flag
        GV4IP(locnum,hsvNew(i,41)+24)=hsvNew(i,24)   ! GII_n
        GV4IP(locnum,hsvNew(i,41)+28)=hsvNew(i,9)    ! Mode mixity
        !GV4Elem(locnum,6)=max(GV4Elem(locnum,6),hsvNew(i,15)) ! max Ds

      end if

    end if ! (totalTime.le.dt)

  end do ! nblock

  end subroutine vumat_cohesive_fatigue


  !> 
  subroutine calcMIIStrengthAndFractureToughness(hsvNew, nblock, nstatev, &
                    FIIlimit_n, FIIlimit_0, FIType, totalTime, fatinit, &
                    FIILimit, TTCELaw, GIIc_n, GIIc, eta_G, eta_F, EI, &
                    ulocal, i, IfFatigue)
    real(wp) ::  hsvNew(nblock, nstatev), totalTime
    real*8 FIILimit, FIIlimit_0, FIIlimit_n, FIType, fatinit, &
          TTCELaw, GIIc, GIIc_n, eta_G, eta_F, EI, ulocal(3), IfFatigue
    integer nblock, nstatev, i

    if (hsvNew(i,14).le.1.0d0) then  ! elastic stage
      FIIlimit_n=FIIlimit_0-eta_F*EI*min(0.d0,hsvNew(i,1)) ! TTCE
      if ((IfFatigue.ne.0.0d0).and.(FIType.eq.2.0d0).and. &
                            (totalTime.gt.fatinit)) then
        FIIlimit=hsvNew(i,10) !read reduced strength if Mukhopadhyay's used
      else	! only affected by TTCE at the elastic stage
        FIIlimit=FIIlimit_n
        hsvNew(i,10)=FIIlimit
      end if
    else ! softening stage
      FIIlimit=hsvNew(i,10) !possibly reduced due to fatigue initiation
    end if

    ! Mode II fracture toughness, possibly enhanced due to TTCE
    if (hsvNew(i,14).le.1.0d0) then ! only consider TTCE at elastic stage
      if (TTCELaw.eq.2.0d0) then  ! law C
        GIIc_n=GIIc*(FIIlimit/FIIlimit_0)**2
      elseif (TTCELaw.eq.1.0d0) then  ! law B
        GIIc_n=GIIc*(FIIlimit/FIIlimit_0)
      else  ! law A
        GIIc_n=GIIC*(1-eta_G*EI*min(0.0d0,ulocal(1)))
      end if
      GIIc_n=max(GIIc_n,GIIc) ! NO effect of fatigue initiation on GIIC
    else
      GIIc_n=hsvNew(i,24)
    end if
    hsvNew(i,24)=GIIc_n

  end subroutine calcMIIStrengthAndFractureToughness
  
  
  !> 
  subroutine staticSoftening(hsvNew, nblock, nstatev, FI, FII, &
                      FIyield, FIIyield, Fmixed, Fyield, Ds, Ds_max, &
                      u, u0, umax, IfFatigue, totalTime, fatinit, &
                      FILimit, FIILimit, FIlimit_0, FIIlimit_0, &
                      a_SNI, a_SNII, cosI, cosII, frq, dt, i, FIType)

    real(wp) :: hsvNew(nblock, nstatev), totalTime, dt

    real*8 FI, FII, FIYield, FIIyield, Fmixed, Fyield, Ds, Ds_max, &
          u, u0, umax, IfFatigue, fatinit, FILimit, &
          FIILimit, FIlimit_0, FIIlimit_0, a_SNI, a_SNII, &
          cosI, cosII, frq, FIType

    real*8 temp, tempI, tempII, uratio, up

    integer nblock, nstatev, i

    if (Ds.eq.Ds_max) then    ! loading

        FI=FIyield*(1.0d0-Ds)
        FII=FIIyield*(1.0d0-Ds)
        Fmixed=Fyield*(1.0d0-Ds)

    elseif (Ds.lt.Ds_max) then ! elastic unloading/reloading

        up=u0+Ds_max*(umax-u0)
        uratio=u/(up+1.e-20)
        FI=FIyield*uratio*(1.0d0-Ds_max)
        FII=FIIyield*uratio*(1.0d0-Ds_max)
        Fmixed=Fyield*uratio*(1.0d0-Ds_max)

    endif

    hsvNew(i,14)=2

    ! If Tao's initiaiton method is used, it futher reduces strength in a
    ! fatigue damage initiating element if it not a fatigue crack tip
    if ((IfFatigue.ne.0.0d0).and.(totalTime.gt.fatinit)) then
        if (FIType.eq.3.0d0) then ! Tao's method
            if (hsvNew(i,36).ge.1.0d0) then
                temp=1.d0-dsqrt((FI/FIlimit_0)**2+(FII/FIIlimit_0)**2)
                temp=10.d0**(temp/(a_SNI*cosI**2.0d0+a_SNII*cosII**2.0d0+1.e-20))
                temp=frq*dt/temp
                tempI=FIlimit*(1.0d0-temp)
                tempII=FIIlimit*(1.0d0-temp)
                hsvNew(i,20)=max(0.0d0,tempI)
                hsvNew(i,10)=max(0.0d0,tempII)
            endif
        endif
    endif

  end subroutine staticSoftening
     

  !> CZM initialisation procedure
  !> Run on startup from vexternaldb (lOp .eq. j_int_StartAnalysis)
  subroutine CZM_initialisation()
    use Fatigue_Globals
    use Mesh_Utils
    use Fatigue, only: JudgeIfFatigue
    use Fatigue_ANN, only: ReadNeuralNetworksParam

    integer :: kProcessNum, threadID
    integer :: lenoutdir, lenjobname
    integer :: get_thread_id
    character*256 outdir,filepath,jobname
    double precision :: Alpha_0, isFatigue

    call vgetrank(kprocessnum)
    threadID=get_thread_id()

    if (kprocessnum==0 .and. threadID==0) then !the first processor/thread

      call vgetoutdir(outdir,lenoutdir)
      call vgetjobname(jobname,lenjobname)
      filepath=outdir(1:lenoutdir)//'/'//jobname(1:lenjobname)//'.inp'

      ! Read fatigue/static switch of the material card in the model .inp file
      call JudgeIfFatigue(filepath,outdir,lenoutdir,isFatigue)

      if (isFatigue.ne.0.d0) then
        write(*,*) 'This model is judged to be a fatigue analysis'
      else
        write(*,*) 'This model is judged to be a static analysis'
      endif

      if (isFatigue.ne.0.d0) then

        ! searching neighbours and computing element dimensions
        Alpha_0=400.d0 !Shreshold angle used to distinguish between inter
                  !and intralaminate elements; Assign a larger than 180
                        !value to it so distinguishing is disabled.

        call SNBCD(filepath,outdir,lenoutdir,maxCohElemNum, &
                    maxNodeNum,Alpha_0,CohElem,IPEdgNbrElem, &
                    IP1NbrElem,IP2NbrElem,IP3NbrElem,IP4NbrElem, &
                    CohElemNum,ElemDimens,NodeCoords)
      !Read Neural networks parameters
        call ReadNeuralNetworksParam(outdir,lenoutdir,NNW,NNB,NNC, &
                                    NNMinMax,NNdSRGthP,NNN)
      end if ! is fatigue

    end if ! first processor/thread

  end subroutine CZM_initialisation

end module CZM
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132
