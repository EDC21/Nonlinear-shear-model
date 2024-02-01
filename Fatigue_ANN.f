!DIR$ FREEFORM

!>  Subroutines to get fatigue crack propagation rate da/dN with Artificial Neural Network
module Fatigue_ANN
  use Abaqus_Definitions, only: wp=>abaqus_real_kind
  use Mesh_Utils, only: to_upper
  implicit none

contains

  !> Read neural network weights from file
  subroutine ReadNeuralNetworksParam(outdir,lenoutdir,NNW,NNB,NNC, &
                                       NNMinMax,NNdSRGthP,NNN)
    integer N,i,j
    character*256 outdir,filepath,line,capitalline
    integer lenoutdir,keyfile
    real*8 NNW(500,3),NNB(500),NNC(500),NNMinMax(4,2),NNdSRGthP(5,5)
    integer NNN

    ! Initialise N for if NN is not used
    N = 0

    ! Read Nerual networks relevant parameters
    filepath=outdir(1:lenoutdir)//'/'//'NNP.inp'
    open(newunit=keyfile,file=filepath,status='unknown')
9000 read(keyfile,9050, end=9600) line
9050 format(a)
    if (line(1:1).EQ.'*') then
      capitalline=to_upper(line)
      if(capitalline(1:4).EQ.'*INN') then ! Inner weights
        N=0 ! number of neurons
9100    read(keyfile,*) line
        if (line(1:1).NE.'*') then
          N=N+1
          backspace(keyfile)
          read(keyfile,*)  NNW(N,1),NNW(N,2),NNW(N,3)
          write(*,*) NNW(N,1),NNW(N,2),NNW(N,3)
          go to 9100  ! To read out all dats under the current keyword
        else
          backspace(keyfile)
          go to 9000 !  Go to read the next keyword
        endif

      elseif (capitalline(1:4).EQ.'*THR') then ! Thresholds
        N=0 ! number of neurons
9200    read(keyfile,*) line
        if (line(1:1).NE.'*') then
          N=N+1
          backspace(keyfile)
          read(keyfile,*)  NNB(N)
          write(*,*) NNB(N)
          go to 9200  ! To read out all dats under the current keyword
        else
          backspace(keyfile)
          go to 9000 !  Go to read the next keyword
        endif

      elseif(capitalline(1:4).EQ.'*OUT') then ! Outer weights
        N=0 ! number of neurons
9300    read(keyfile,*) line
        if (line(1:1).NE.'*') then
          N=N+1
          backspace(keyfile)
          read(keyfile,*)  NNC(N)
          write(*,*) NNC(N)
          go to 9300  ! To read out all dats under the current keyword
        else
          backspace(keyfile)
          go to 9000 !  Go to read the next keyword
        endif

      elseif(capitalline(1:4).EQ.'*DSQ') then ! change of square root of Gth
        i=0
9400    read(keyfile,*) line
        if (line(1:1).NE.'*') then
          i=i+1
          backspace(keyfile)
          if (i.eq.1) then
            read(keyfile,*)  NNdSRGthP(i,1),NNdSRGthP(i,2), &
                             NNdSRGthP(i,3),NNdSRGthP(i,4),NNdSRGthP(i,5)
          else
            read(keyfile,*)  NNdSRGthP(i,1),NNdSRGthP(i,2), &
                                     NNdSRGthP(i,3),NNdSRGthP(i,4)
          endif
          write(*,*) NNdSRGthP(i,1),NNdSRGthP(i,2), &
                          NNdSRGthP(i,3),NNdSRGthP(i,4)
          go to 9400  ! To read out all dats under the current keyword
        else
          backspace(keyfile)
          go to 9000 !  Go to read the next keyword
        endif

      elseif(capitalline(1:4).EQ.'*INP') then ! Input limits
        i=0
9500    read(keyfile,*) line
        if (line(1:1).NE.'*') then
          i=i+1
          backspace(keyfile)
          read(keyfile,*)  NNMinMax(i,1), NNMinMax(i,2)
          write(*,*) NNMinMax(i,1), NNMinMax(i,2)
          go to 9500  ! To read out all dats under the current keyword
        else
          backspace(keyfile)
          go to 9000 !  Go to read the next keyword
        endif

      endif
    endif
    goto 9000  ! Go to read next line
9600 continue
    close(keyfile)

    NNN=N

  end subroutine ReadNeuralNetworksParam


  !> Non-Abaqus subroutines of Neural Networks for da/dN
  subroutine NeuralNetworksDaDN(dadN,E11,E33,G13,v13, &
            FILimit,FIILimit,GIC,GIIC_n,R,Gc,Gmax,Gmin,NumIPs,Phi, &
            IfElemTip,chif,chis,Z,NNW,NNB,NNC,NNMinMax,NNdSRGthP,NNN)

    real*8, parameter :: zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, &
            third = one/three, four = 4.d0,  half = 0.5d0, sf=1.e-20
    real*8 dadN,E11,E33,G13,v13,FILimit,FIILimit, &
           GIC,GIIC_n,R,Gc,Gmax,Gmin,Phi(4),IfElemTip
    real*8 alphaI,alphaII,chif,chis,Z,LIC,LIIC,LC, &
           AvgdSRGth, dSRG,lamda,rho,sigmaC,Y
    integer NumIPs,i,j

    real*8 NNW(500,3),NNB(500),NNC(500),NNMinMax(4,2),NNdSRGthP(5,5)
    integer NNN

    ! Compute Chif and scaling
    dSRG=dsqrt(Gmax)*(one-R) ! by Eq. 9 in Allegri, M&D paper
    ! Compute the averaged delta_square root of Gth
    AvgdSRGth=zero
    if (IfElemTip.ne.zero)then
      j=1
    else
      j=NumIPs
    endif
    do i=1,j
      if (Phi(i).le.NNdSRGthP(1,2)) then
        AvgdSRGth=AvgdSRGth+ &
                     NNdSRGthP(2,1)*(Phi(i)-NNdSRGthP(1,1))**3+ &
                     NNdSRGthP(2,2)*(Phi(i)-NNdSRGthP(1,1))**2+ &
                     NNdSRGthP(2,3)*(Phi(i)-NNdSRGthP(1,1))+ &
                     NNdSRGthP(2,4)
      elseif (Phi(i).le.NNdSRGthP(1,3)) then
        AvgdSRGth=AvgdSRGth+ &
                     NNdSRGthP(3,1)*(Phi(i)-NNdSRGthP(1,2))**3+ &
                     NNdSRGthP(3,2)*(Phi(i)-NNdSRGthP(1,2))**2+ &
                     NNdSRGthP(3,3)*(Phi(i)-NNdSRGthP(1,2))+ &
                     NNdSRGthP(3,4)
      elseif (Phi(i).le.NNdSRGthP(1,4)) then
        AvgdSRGth=AvgdSRGth+ &
                     NNdSRGthP(4,1)*(Phi(i)-NNdSRGthP(1,3))**3+ &
                     NNdSRGthP(4,2)*(Phi(i)-NNdSRGthP(1,3))**2+ &
                     NNdSRGthP(4,3)*(Phi(i)-NNdSRGthP(1,3))+ &
                     NNdSRGthP(4,4)
      else
        AvgdSRGth=AvgdSRGth+ &
                     NNdSRGthP(5,1)*(Phi(i)-NNdSRGthP(1,4))**3+ &
                     NNdSRGthP(5,2)*(Phi(i)-NNdSRGthP(1,4))**2+ &
                     NNdSRGthP(5,3)*(Phi(i)-NNdSRGthP(1,4))+ &
                     NNdSRGthP(5,4)
      endif
    enddo
    AvgdSRGth=AvgdSRGth/j

    ! Compute Gc for each IP and average
    Gc=zero
    do i=1,j
      Gc=Gc+GIC*GIIC_n/(GIIC_n*(1-Phi(i))+GIC*Phi(i))
    enddo
    Gc=Gc/j

    chif=(dSRG-AvgdSRGth)/(dsqrt(Gc)+sf)   ! chif
    if (chif.lt.zero) then
      chif=zero
    endif
    chif=(chif-NNMinMax(1,1))/(NNMinMax(1,2)-NNMinMax(1,1)+sf)! scaled to [0,1]

    ! Compute Chis and scaling
    chis=one-dSRG/(one-R)/(dsqrt(Gc)+sf)  ! chis
    if (chis.lt.zero) then
      chis=zero
    endif
    chis=(chis-NNMinMax(2,1))/(NNMinMax(2,2)-NNMinMax(2,1)+sf)! scaled to [0,1]

    ! Compute Z and scaling
    lamda=E33/E11; rho=dsqrt(E33/E11)*(E11/two/G13-v13)
    alphaI=(lamda**(-one/four))*dsqrt((1+rho)/E11/E33/two)
    alphaII=(lamda**(one/four))*dsqrt((1+rho)/E11/E33/two)
    LIC=GIC/two/3.1415926/alphaI/FILimit/FILimit
    LIIC=GIIC_n/two/3.1415926/alphaII/FIILimit/FIILimit
    LC=zero; sigmaC=zero
    do i=1,j
      LC = LC + (LIC*(one-Phi(i))*GIIC_n+LIIC*Phi(i)*GIC) &
                  /(GIC*Phi(i)+GIIC_n*(one-Phi(i)))
      sigmaC = sigmaC+FILimit*FIILimit/ &
             dsqrt(FILimit*FILimit*Phi(i)+FIILimit*FIILimit*(one-Phi(i)))
    enddo
    LC=LC/j;  sigmaC=sigmaC/j

    Z=sigmaC*LC/Gc  !Z(Phi)
    Z=(Z-NNMinMax(3,1))/(NNMinMax(3,2)-NNMinMax(3,1)+sf)! scaled to [0,1]

    ! Substiute into Neural networks formula to get da/dN
    Y=zero
    do j=1,NNN
      Y=Y+NNC(j)* &
             1/(1+dexp(-(NNW(j,1)*chif+NNW(j,2)*chis+NNW(j,3)*Z+NNB(j))))
    enddo
    !scaled to between -1 and 1
    Y=Y*(NNMinMax(4,2)-NNMinMax(4,1))+NNMinMax(4,1)
    if (Y.lt.-1.d0) then
      Y=-0.99
    endif
    if (Y.gt.1.d0) then
      Y=0.99
    endif
    !Scaled back to real da/dN, by Eq. 39
    dadN=LC*10**(tan(3.1415926*Y/2));

  end subroutine NeuralNetworksDaDN

end module Fatigue_ANN
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132


