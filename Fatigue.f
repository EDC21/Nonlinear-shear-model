#include 'Fatigue_CharacLength.f'

!DIR$ FREEFORM

!> Module for Fatigue subroutines
module Fatigue
  use Abaqus_Definitions, only: wp=>abaqus_real_kind
  use Mesh_Utils, only: to_upper
  implicit none

contains

  !> Compute Paris parameters for current IP, average G,Gc,c,m and Ds if needed
  subroutine calculateParisLawParams(FLawType, MMPPType, cMM, cI, &
                   cm, Gratio, cII, expMM, expI, expII, expm, C, &
                   Gc_IP, GIc, GIIc_n, h, b0I, R, &
                   hsvNew, nblock, nstatev, IfElemTip, &
                   locnum, NumIPs, i)
    use Fatigue_Globals
    implicit none
    real(wp) hsvNew(nblock, nstatev)
    real*8 FLawType, MMPPType, cMM, cI, cII, cm, expMM, expI, expm, &
              expII, Gratio, C, Gc_IP, GIc, GIIc_n, h, b0I, R, &
             IfElemTip
    integer nblock, nstatev, i, m, locnum, NumIPs
    real*8 G_Sum, Gc_Sum, C_Sum, exp_Sum, Ds_Sum


    if (FLawType.ne.4.d0) then ! If ANN not used
      if (FLawType.ne.0.d0) then ! If Allegri's Paris law not used
        ! Calculate cMM and expMM
        if (MMPPType.eq.1.d0) then
          call BlancoMMPP(Gratio, cI, cII, expI, expII, cm, expm, &
                                        cMM, expMM)
        elseif (MMPPType.eq.2.d0) then
          call KeaneBenzeggaghMMPP(Gratio, cI, cII, expI, expII, &
                                                 cm, expm, cMM, expMM)
        else
          call RussellStreetMMPP(Gratio, cI, cII, expI, expII, &
                                               cMM, expMM)
        end if
      else
        call AllegriParisLaw(Gratio, Gc_IP, GIc, GIIC_n, b0I, h, &
                                       C, R, cMM, expMM)
      end if

      hsvNew(i,34) = cMM         !C
      hsvNew(i,35) = expMM       !m
    end if

    ! If element level tip tracking used, average ERR, Gc, Paris parameters
    ! and Ds when IP 1 is called
    if (IfElemTip.ne.0.d0) then
      if (hsvNew(i,41).eq.1.d0) then
        G_Sum = 0;  Gc_Sum = 0;  C_Sum = 0;  exp_Sum = 0;  Ds_Sum = 0
        do m=1,NumIPs
          G_Sum = G_Sum+GV4IP(locnum,(m-1)*5+1)
          Gc_Sum = Gc_Sum+GV4IP(locnum,(m-1)*5+2)
          C_Sum = C_Sum+GV4IP(locnum,(m-1)*5+3)
          exp_Sum = exp_Sum+GV4IP(locnum,(m-1)*5+4)
          Ds_Sum = Ds_Sum+GV4IP(locnum,m*5)
        end do
        GV4Elem(locnum,2) = G_Sum/NumIPs   ! ERR
        GV4Elem(locnum,3) = Gc_Sum/NumIPs  ! Gc
        GV4Elem(locnum,4) = C_Sum/NumIPs   ! Paris constant
        GV4Elem(locnum,5) = exp_Sum/NumIPs ! Paris power
        GV4Elem(locnum,6) = Ds_Sum/NumIPs  ! Static damage variable
      end if
      ! Save the averaged G value obtained in the last time incremet
      hsvNew(i,40) = GV4Elem(locnum,2)
    end if

  end subroutine calculateParisLawParams


  subroutine RussellStreetMMPP(Gratio, cI, cII, expI, expII, cMM, expMM)
    double precision, intent (in) :: Gratio, cI, cII, expI, expII
    double precision, intent (out) :: cMM, expMM

    cMM   = (1-Gratio)*cI   + Gratio*cII
    expMM = (1-Gratio)*expI + Gratio*expII

  end subroutine RussellStreetMMPP


  subroutine KeaneBenzeggaghMMPP(Gratio, cI, cII, expI, expII, cm, expm, &
                                   cMM, expMM)
    double precision, intent (in) :: Gratio, cI, cII, expI, expII, cm, expm
    double precision, intent (out) :: cMM, expMM

    cMM = exp( log(cII) + (log(cI) - log(cII)) * ((1-Gratio)**cm) )
    expMM = expI + (expII-expI)*(Gratio**expm)

  end subroutine KeaneBenzeggaghMMPP


  subroutine BlancoMMPP(Gratio, cI, cII, expI, expII, cm, expm, &
                          cMM, expMM)
    double precision, intent (in) :: Gratio, cI, cII, expI, expII, cm, expm
    double precision, intent (out) :: cMM, expMM

    cMM = log10(cI) + log10(cm)*Gratio + log10(cII/cI/cm)*(Gratio**2)
    cMM = 10**cMM
    expMM = expI + expm*Gratio + (expII-expI-expm)*(Gratio**2)

  end subroutine BlancoMMPP


  subroutine AllegriParisLaw(Gratio, Gc_IP, GIc, GIIC_n, b0I, h, C, R, &
                               cMM, expMM)
    double precision, intent (in) :: Gratio, Gc_IP, GIc, GIIC_n, &
                                        b0I, h, C, R
    double precision, intent (out) :: cMM, expMM
    double precision :: alpha_phi

    cMM = C
    alpha_phi = (Gc_IP-GIc)/(GIIc_n-GIc)
    expMM = b0I*exp(-h*Gratio)/((1-R)**(1+alpha_phi))

  end subroutine AllegriParisLaw


  !> Check if damage is initiated due to fatigue using Zhang's Method
  subroutine zhangsMethod_DamageInitiated(hsvNew, nblock, nstatev, &
                      EI, u, cosI, FILimit, a_SNI, a_SNII, frq, dt, EII, cosII, &
                      FIILimit, uratio, i, locnum)
    use Fatigue_Globals
    implicit none
    real(wp) hsvNew(nblock, nstatev), dt
    real*8 EI, u, cosI, FILimit, a_SNI, a_SNII, frq, EII, cosII, &
          FIILimit, uratio, Di, temp
    integer i, locnum, nblock, nstatev

    temp = 10.d0**((1-EI*u*cosI/FIlimit)/a_SNI)
    hsvNew(i,36) = hsvNew(i,36)+frq*dt/temp
    temp = 10.d0**((1-EII*u*cosII/FIIlimit)/a_SNII)
    hsvNew(i,37) = hsvNew(i,37)+frq*dt/temp
    Di = hsvNew(i,36)+hsvNew(i,37)
    if (Di.ge.1.d0) then
      if (GV4Elem(locnum,8).eq.0.d0) then
        GV4Elem(locnum,8) = hsvNew(i,33)
      end if
      hsvNew(i,20) = FIlimit*uratio
      hsvNew(i,10) = FIIlimit*uratio
      hsvNew(i,14) = 2
    end if
  
  end subroutine zhangsMethod_DamageInitiated


  !> Check if damage is initiated due to fatigue using May Tao Method
  subroutine mayTaoMethod_DamageInitiated(hsvNew, nblock, nstatev, &
                  EI, u, cosI, FILimit, a_SNI, a_SNII, frq, dt, EII, cosII, &
                  FIILimit, uratio, i, locnum)
    use Fatigue_Globals
    implicit none

    real(wp) hsvNew(nblock, nstatev), dt
    real*8 EI, u, cosI, FILimit, a_SNI, a_SNII, frq, EII, cosII, &
          FIILimit, uratio, temp
    integer i, locnum, nblock, nstatev

    temp = 1.0d0-dsqrt((EI*u*cosI/FIlimit)**2 + &
                    (EII*u*cosII/FIIlimit)**2)
    temp = 10.d0**(temp/(a_SNI*cosI**2.0d0+a_SNII*cosII**2.0d0+1.e-20))
    hsvNew(i,36) = hsvNew(i,36)+frq*dt/temp
    if (hsvNew(i,36).ge.1.0d0) then
        hsvNew(i,36) = 1.0d0
        if (GV4Elem(locnum,8).eq.0.0d0) then
            GV4Elem(locnum,8) = hsvNew(i,33)
        end if
        hsvNew(i,20) = FIlimit*uratio
        hsvNew(i,10) = FIIlimit*uratio
        hsvNew(i,14) = 2
    end if

  end subroutine mayTaoMethod_DamageInitiated


  !> Check for fatigue softening 
  subroutine checkFatigueSoftening(hsvNew, nblock, nstatev, locnum, i, &
                  IfNoTrack, IfElemTip, IPnum4Tip)

    use Fatigue_Globals
    implicit none

    real(wp) hsvNew(nblock, nstatev)
    real*8 IfNoTrack, IfElemTip, IPnum4Tip, G_e, &
             edge_1, edge_2, edge_3, edge_4, edge_5, edge_6, edge_7, edge_8
    integer nblock, nstatev, locnum, i, j, k

    edge_1 = hsvNew(i, 50); edge_2 = hsvNew(i, 51); edge_3 = hsvNew(i, 52)
    edge_4 = hsvNew(i, 53); edge_5 = hsvNew(i, 54); edge_6 = hsvNew(i, 55)
    edge_7 = hsvNew(i, 56); edge_8 = hsvNew(i, 57)
    G_e = GV4Elem(locnum, 2)


    if (IfNoTrack.eq.1) then ! Tip tracking not used
      hsvNew(i,14) = 3
      return
    else
      if (IfElemTip.ne.0.0d0) then ! Element level tracking
        ! First, check for failed neighbours
        do j=1,4
          ! bottom facet connected neighbours
          if (hsvNew(i,j+49).gt.0) then
            if (GV4Elem(hsvNew(i,j+49),1).eq.IPnum4Tip) then
              hsvNew(i,14) = 3
              return
            end if
          end if
          ! top facet connected neighbours
          if (hsvNew(i,j+53).gt.0) then
            if (GV4Elem(hsvNew(i,j+53),1).eq.IPnum4Tip) then
              hsvNew(i,14) = 3
              return
            end if
          end if
        end do

        ! Second, check if the element at edge
        if ((edge_1.eq.0).and.(edge_5.eq.0)) then !edge 1-2 and 5-6
          if (edge_2.gt.0) then
            if (G_e.gt.GV4Elem(edge_2,2)) then ! G_e > G_n2
              hsvNew(i,14) = 3
              return
            end if
          end if
          if (edge_6.gt.0) then
            if (G_e.gt.GV4Elem(edge_6,2)) then ! G_e > G_n6
              hsvNew(i,14) = 3
              return
            end if
          end if
        end if

        if ((edge_2.eq.0).and.(edge_6.eq.0)) then !edge 3-4 and 7-8
          if (edge_1.gt.0) then
            if (G_e.gt.GV4Elem(edge_1,2)) then ! G_e > G_n1
              hsvNew(i,14) = 3
              return
            end if
          end if
          if (edge_5.gt.0) then
            if (G_e.gt.GV4Elem(edge_5,2)) then ! G_e > G_n5
              hsvNew(i,14) = 3
              return
            end if
          end if
        end if

        if ((edge_3.eq.0).and.(edge_7.eq.0)) then !edge 2-3 and 6-7
          if (edge_4.gt.0) then
            if (G_e.gt.GV4Elem(edge_4,2)) then ! G_e > G_n4
              hsvNew(i,14) = 3
              return
            end if
          end if
          if (edge_8.gt.0) then
            if (G_e.gt.GV4Elem(edge_8,2)) then ! G_e > G_n8
              hsvNew(i,14) = 3
              return
            end if
          end if
        end if

        if ((edge_4.eq.0).and.(edge_8.eq.0)) then !edge 1-4 and 5-8
          if (hsvNew(i,52).gt.0) then
            if (G_e.gt.GV4Elem(edge_3,2)) then ! G_e > G_n3
              hsvNew(i,14) = 3
              return
            end if
          end if
          if (edge_7.gt.0) then
            if (G_e.gt.GV4Elem(edge_7,2)) then ! G_e > G_n7
              hsvNew(i,14) = 3
              return
            end if
          end if
        end if

      else  ! IP-based crack tip judgement

        ! Check sister IPs' status
        do j=1,4
          if (j.ne.hsvNew(i,41)) then ! not comparing with itself
            if (GV4IP(locnum,j+20).eq.4.0d0) then
              hsvNew(i,14) = 3 ! sister IP faled
              return
            end if
          end if
        end do
        ! Check IPs of neighbour elements
        do j=1,11
          if (hsvNew(i,j+69).gt.0.0d0) then ! neighbour element exists
            do k=1,4
              if (GV4IP(hsvNew(i,j+69),k+20).eq.4.0d0) then
                hsvNew(i,14) = 3 ! neighbour IP faled
                return
              end if
            end do
          else
            exit
          end if
        end do

        ! check if the IP at edge or discontiuous location
        if (hsvNew(i,81).eq.1.0d0) then
          hsvNew(i,14) = 3
          return
        end if

      end if	! if (IfElemTip.ne.zero) then
    end if ! end of 'if (IfNoTrack.eq.1)'

  end subroutine checkFatigueSoftening


  !> Fatigue softening stage
  subroutine fatigueSoftening(hsvNew, nblock, nstatev, i, FDRType, &
                  IfElemTip, locnum, Gmax, Gc, R, FLawType, &
                  cMM, expMM, dadN, E11, E33, G13, v13, FIlimit, &
                  FIIlimit, GIC, GIIC_n, NumIPs, Z, &
                  FCLType, ElemType, FCL, frq, dt, Ds_max, &
                  D_inc, Dtotal, FI, FII, Fmixed, FTType, FIyield, &
                  FIIyield, Fyield, totalTime, u0, umax, u)
    use Fatigue_Globals
    use Fatigue_CharacLength, only: calcCharacteristicLength
    implicit none

    real(wp) hsvNew(nblock, nstatev), dt, totalTime

    integer nblock, nstatev, i, locnum, NumIPs

    real*8 FDRType, IfElemTip, Gmax, Gc, R, FLawType, cMM, &
             expMM, dadN, E11, E33, G13, v13, FILimit, FIILimit, GIC, &
             GIIC_n, Z, FCLType, ElemType, FCL, frq, Ds_max, &
             D_inc, FI, FII, Fmixed, FTType, FIyield, &
             FIIyield, Fyield, u0, umax, u, Dtotal
    real*8 up, uratio, deltaG, Df

    ! Set fatigue softening flag and start cyclc number
    if(hsvNew(i,27).eq.0.d0)then
      hsvNew(i,27) = 1.d0
      hsvNew(i,38) = hsvNew(i,33) ! save fatigue softening start cycle number
    end if

    ! Achieve Gmax and Gc
    call calculateGMaxAndGC(hsvNew, nblock, nstatev, FDRType, IfElemTip, &
                              i, locnum, Gmax, Gc)

    ! Compute dG
    deltaG = Gmax*(1.0d0-R*abs(R))

    ! Retrieve Pairs parameters if ANN not used
    if (FLawType.ne.4.0d0) then
      if (IfElemTip.ne.0.0d0) then
        cMM = GV4Elem(locnum,4)
        expMM = GV4Elem(locnum,5)
      end if
    end if

    ! Compute da/dN
    call calculate_dadN(FLawType, cMM, deltaG, Gc, expMM, locnum, &
                  hsvNew, nblock, nstatev, i, dadN, &
                  E11, E33, G13, v13, FILimit, FIILimit, GIC, GIIC_n, &
                  R, Gmax, NumIPs, IfElemTip, Z)

    ! Characteristic length, value returned in FCL
    call calcCharacteristicLength(hsvNew, FCLType, ElemType, FCL, &
                                    nblock, nstatev, locnum, i)


    ! Half the length if IP based modelling
    if (IfElemTip.ne.0.0d0) then
      hsvNew(i,49) = FCL
    else
      hsvNew(i,49) = FCL/2.0d0
    end if

    ! Compute fatigue damage increment
    call calculateFatigueDamageInc(FCL, dadN, FDRType, &
                  frq, dt, Ds_max, D_inc, hsvNew, nblock, nstatev, i)

    ! Compute fatigue damage variable and total damage variable
    Df = hsvNew(i,30)+ D_inc
    hsvNew(i,30) = Df
    Dtotal = Ds_max+Df
    hsvNew(i,31) = Dtotal

    ! Compute mode I, II tractions
    if(Dtotal.ge.1.0d0)then !Fail in this dt
      FI = 0.0d0
      FII = 0.0d0
      Fmixed = 0.0d0
      hsvNew(i,14) = 4
      hsvNew(i,28) = 0.0d0
      hsvNew(i,32) = totalTime  ! recording failure time
      GV4Elem(locnum,1) = GV4Elem(locnum,1)+1.0d0 !update the no. of failed IPs

    else ! Not failed
      if (FTType.eq.1.0d0) then  ! by the yield strengths
        FI = FIyield*(1.0d0-Dtotal)
        FII = FIIyield*(1.0d0-Dtotal)
        Fmixed = Fyield*(1.0d0-Dtotal)
      else   ! by the unloading stiffness
        up = u0+Dtotal*(umax-u0)
        uratio = u/(up+1.0e-20)
        FI = FIyield*uratio*(1.0d0-Dtotal)
        FII = FIIyield*uratio*(1.0d0-Dtotal)
        Fmixed = Fyield*uratio*(1.0d0-Dtotal)
      end if
      hsvNew(i,14) = 3
    end if

  end subroutine fatigueSoftening

  
  !> Compute fatigue damage increment
  subroutine calculateFatigueDamageInc(FCL, dadN, FDRType, &
    frq, dt, Ds_max, D_inc, hsvNew, nblock, nstatev, i)
    implicit none
    real(wp), intent(inout) :: hsvNew(nblock, nstatev)
    real(wp), intent(in) :: dt
    double precision, intent (in) :: FCL, dadN, FDRType, frq, Ds_max
    double precision, intent (inout) :: D_inc
    integer, intent (in) :: nblock, nstatev, i
    double precision :: failure_cycles, remaining_cycles, failure_dtNum, rd

    failure_cycles=FCL/dadN
    if (FDRType.eq.2.0d0) then ! Tao's method
      remaining_cycles=failure_cycles-(hsvNew(i,33)-hsvNew(i,38))
      if (remaining_cycles.gt.0.0d0) then
        failure_dtNum=remaining_cycles/frq/dt
        rd=1.0d0-(0.01/(1-Ds_max-hsvNew(i,30)))**(1/failure_dtNum)
        rd=min(max(rd, 0.0d0), 1.0d0)
        D_inc=1-Ds_max-hsvNew(i,30)
        if (D_inc.gt.0.0d0) then
          D_inc=D_inc*rd
        else
          D_inc=0.5d0*(D_inc+hsvNew(i,31)-hsvNew(i,30)-Ds_max) ! why?
        end if
      else
        rd=0
        D_inc=1.0d0-Ds_max-hsvNew(i,30)
      end if
      hsvNew(i,39)=rd
    else  ! Kawashita's method
      D_inc=((1.0d0-Ds_max)/failure_cycles)*frq*dt
    end if

    return
  end subroutine calculateFatigueDamageInc


  !> Calculate fatigue crack propagation rate
  subroutine calculate_dadN(FLawType, cMM, deltaG, Gc, expMM, locnum, &
    hsvNew, nblock, nstatev, i, dadN, &
    E11, E33, G13, v13, FILimit, FIILimit, GIC, GIIC_n, &
    R, Gmax, NumIPs, IfElemTip, Z)
    use Fatigue_Globals
    use Fatigue_ANN, only: NeuralNetworksDaDN
    implicit none

    real(wp) hsvNew(nblock, nstatev)

    real*8 FLawType, cMM, deltaG, expMM, dadN, E11, E33, G13, v13, &
      FILimit, FIILimit, GIC, GIIC_n, R, Gc, Gmax, IfElemTip, Z
    integer nblock, nstatev, i, locnum, NumIPs
    real*8 Phi(4), chif, chis, Gmin

    if(FLawType.eq.2.0d0) then
      dadN = cMM*(deltaG/Gc)**expMM
    elseif(FLawType.eq.3.0d0) then
      dadN = cMM*(deltaG)**expMM
    elseif(FLawType.eq.4.0d0) then
      Phi(1)=GV4IP(locnum,29)
      if (NumIPs.gt.1.0d0) then
        Phi(2)=GV4IP(locnum,30);Phi(3)=GV4IP(locnum,31)
        Phi(4)=GV4IP(locnum,32)
      end if

      call NeuralNetworksDaDN(dadN,E11,E33,G13,v13, &
        FILimit,FIILimit,GIC,GIIC_n,R,Gc,Gmax,Gmin,NumIPs,Phi, &
        IfElemTip,chif,chis,Z,NNW,NNB,NNC,NNMinMax,NNdSRGthP,NNN)
      hsvNew(i,39)=chif
      hsvNew(i,34)=chis
      hsvNew(i,35)=Z
    else
      dadN = cMM*(Gmax/Gc)**expMM
    end if
    hsvNew(i,29)=dadN

  end subroutine calculate_dadN


  subroutine calculateGMaxAndGC(hsvNew, nblock, nstatev, FDRType, IfElemTip, &
                                      i, locnum, Gmax, Gc)

    use Fatigue_Globals
    implicit none
    real(wp) :: hsvNew(nblock, nstatev)
    integer nblock, nstatev, m, n, j, i, k, locnum
    real*8 FDRType, IfElemTip, Gmax, Gc

    if (FDRType.eq.1) then !Non-local, maximum Gmax of itself and neighbours
      if (IfElemTip.ne.0.d0) then !Element level tracking
        Gmax=GV4Elem(locnum,2)
        Gc=GV4Elem(locnum,3)
        do j=1,4
          m=int(hsvNew(i,j+49))
          n=int(hsvNew(i,j+53))
          if(m.gt.0)then
            if (Gmax.lt.GV4Elem(m,2)) then
              Gmax=GV4Elem(m,2)
              Gc=GV4Elem(m,3)
            end if
          end if
          if(n.gt.0) then
            if (Gmax.lt.GV4Elem(n,2)) then
              Gmax=GV4Elem(n,2)
              Gc=GV4Elem(n,3)
            end if
          end if
        end do
      else !IP level tracking
        Gmax=hsvNew(i,8)
        Gc=hsvNew(i,25)
        do j=1,4 ! Compare with sister IPs
          if (j.ne.hsvNew(i,41)) then ! not comparing with itself
            if (Gmax.lt.GV4IP(locnum,j*5-4)) then
              Gmax=GV4IP(locnum,j*5-4)
              Gc=GV4IP(locnum,j*5-3)
            end if
          end if
        end do
        do j=1,11 ! Compare with IPs of neighbour elements
          if(hsvNew(i,j+69).gt.0.d0)then ! neighbour element exists
            do k=1,4
              if (Gmax.lt.GV4IP(hsvNew(i,j+69),k*5-4)) then
                Gmax=GV4IP(hsvNew(i,j+69),k*5-4)
                Gc=GV4IP(hsvNew(i,j+69),k*5-3)
              end if
            end do
          end if
        end do
      end if
    else ! local, its own values
      if (IfElemTip.ne.0.d0) then !Element level tracking
        Gmax=GV4Elem(locnum,2)
        Gc=GV4Elem(locnum,3)
      else !IP level tracking
        Gmax=hsvNew(i,8)
        Gc=hsvNew(i,25)
      end if
    end if

  end subroutine calculateGMaxAndGC


  subroutine storeNeighboursAndDimensions(hsvNew, nblock, nstatev, &
                  i, locnum, NumIPs, T, ElemDimens, NodeCoords, &
                  IPEdgNbrElem, CohElemNum, CohElem, &
                  IP1NbrElem, IP2NbrElem, IP3NbrElem, IP4NbrElem)

    use Fatigue_Globals, only: maxCohElemNum, maxNodeNum
    implicit none
    real(wp) hsvNew(nblock, nstatev)
    integer i, j, k, locnum, nblock, nstatev, NumIPs
    real*8 T

    real*8 ElemDimens(maxCohElemNum, 4), NodeCoords(maxNodeNum, 3)
    integer CohElem(maxCohElemNum,9), IPEdgNbrElem(maxCohElemNum,9), &
              IP1NbrElem(maxCohElemNum,12), IP2NbrElem(maxCohElemNum,12), &
              IP3NbrElem (maxCohElemNum,12), IP4NbrElem(maxCohElemNum,12), &
              CohElemNum


    do j = 1, CohElemNum
      if (IPEdgNbrElem(j,1).eq.locnum) then !Local ID in the SNBCD com blocks
        ! All in-plane edge neighbours of the IP host element
        do k=1,4
          hsvNew(i,49+k)=IPEdgNbrElem(j,k+1)
          hsvNew(i,53+k)=IPEdgNbrElem(j,k+5)
        end do

        !Element geometric thickness, length, width and mid-plane area
        T=ElemDimens(j,1) !thickness
        hsvNew(i,45)=ElemDimens(j,2) !length
        hsvNew(i,46)=ElemDimens(j,3) !width
        hsvNew(i,47)=ElemDimens(j,4) !mid-plane area

        ! Mid-plane corner coodinates, for dynamic char. length if needed
        do k=1,3
          hsvNew(i,57+k)=0.5d0*(NodeCoords(CohElem(j,2),k)+ &
                                             NodeCoords(CohElem(j,6),k))
          hsvNew(i,60+k)=0.5d0*(NodeCoords(CohElem(j,3),k)+ &
                                             NodeCoords(CohElem(j,7),k))
          hsvNew(i,63+k)=0.5d0*(NodeCoords(CohElem(j,4),k)+ &
                                             NodeCoords(CohElem(j,8),k))
          hsvNew(i,66+k)=0.5d0*(NodeCoords(CohElem(j,5),k)+ &
                                             NodeCoords(CohElem(j,9),k))
        end do

        ! IP host vertical edge neighbour elements if COH3D8 element used
        if (NumIPs.eq.4.0d0) then
          if (hsvNew(i,41).eq.1.0d0) then
            hsvNew(i,70:81)=IP1NbrElem(j,1:12)
          elseif (hsvNew(i,41).eq.2.0d0) then
            hsvNew(i,70:81)=IP2NbrElem(j,1:12)
          elseif (hsvNew(i,41).eq.3.0d0) then
            hsvNew(i,70:81)=IP3NbrElem(j,1:12)
          elseif (hsvNew(i,41).eq.4.0d0) then
            hsvNew(i,70:81)=IP4NbrElem(j,1:12)
          end if
        end if

        exit
      end if
    end do

  end subroutine storeNeighboursAndDimensions


  !> Checking if model is fatigue or static by reading .inp file
  subroutine JudgeIfFatigue(filepath,outdir,lenoutdir,Fatigue)

    real*8, parameter :: zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0, &
              three = 3.d0, four = 4.d0, seven = 7.d0
    character*256 outdir,filepath,line,capitalline
    integer lenoutdir,keyfile,a(2),i,j,k,numdata
    real*8 Fatigue,b,c
    character(len=8) :: s1(2)=(/"*MATERIA", "COHESIVE"/)

    Fatigue=0; i=0; j=0
    open(newunit=keyfile,file=filepath,status='unknown')
1500 read(keyfile,1550, end=1800) line
1550 format(a)
    if (line(1:1).eq.'*') then
      capitalline=to_upper(line)
      a(1:2)=index(capitalline,s1)
      if ((a(1).ne.zero).and.(a(2).ne.zero)) then !'*Material, name=COHESIVE'
1600    read(keyfile,*,end=1800) line
        i=i+1
        capitalline=to_upper(line)
        if ((capitalline(1:5).eq.'*USER').and.(i.le.seven)) then !'User'
1700      read(keyfile,1550,end=1800) line
          if (line(1:1).ne.'*') then ! material card line
            j=j+1
            if (j.eq.two)then
              numdata = count([(line(k:k),k=1,len(line))].eq.',')
              If (numdata.le.3) then !2 variables with 2 or 3 commas, then static
                go to 1800
              else  ! more than 2 variables, then read the value of 'IfFatigue'
                backspace(keyfile)
                read(keyfile,*,end=1800) b,c,Fatigue
                go to 1800
              end if
            end if
          end if !End of 'material card line'
          go to 1700
        end if !End of the line with 'User'
        go to 1600
      end if !End of the line of '*Material, name=COHESIVE'
      go to 1500
    end if
    goto 1500  ! End of Keep reading lines
1800 continue
    close(keyfile)

  end subroutine JudgeIfFatigue

end module Fatigue
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132
