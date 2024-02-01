!DIR$ FREEFORM
!> Subroutines for Weibull's fibre failure criterion
module Weibull_Fibre_Failure
  use Abaqus_Definitions, only: wp=>abaqus_real_kind
  implicit none

contains

  !> VUMAT interface
  !>
  !> State variables:
  !>   hsvNew(i,1)  element ID
  !>   hsvNew(i,2)  element volume
  !>   hsvNew(i,3)  element ID with the maximum S11 for the last step
  !>   hsvNew(i,4)  Weibull integral of the last step
  !>   hsvNew(i,5)  element
  !>   hsvNew(i,6)  element deletion flag, assigned 0 for deletion
  !>   hsvNew(i,7)  element damage variable; 0, healthy; 3, failed
  !>
  subroutine vumat_mlaminate ( &
         nblock, ndir, nshr, nstatev, nfieldv, nprops, lanneal, &
         stepTime, totalTime, dt, cmname, coordMp, charLength, &
         props, density, strainInc, relSpinInc, &
         tempOld, stretchOld, defgradOld, fieldOld, &
         stressOld, hsvOld, enerInternOld, enerInelasOld, &
         tempNew, stretchNew, defgradNew, fieldNew, &
         stressNew, hsvNew, enerInternNew, enerInelasNew, &
         nMatPoint, nLayer, nSecPoint, nElement)

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

    real*8, parameter :: zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0, &
              three = 3.d0, four = 4.d0

    double precision E11,E22,E33,v12,v21,v13,v31,v23,v32,G12,G13,G23,G31, &
     C11,C12,C13,C21,C22,C23,C31,C32,C33,UDS1,tempVar,dTemper, &
     dEpsilon1,dEpsilon2,dEpsilon3,Alpha1,Alpha2,Alpha3,A,B,TauTrial, &
     TempTEps,Eep,dlambda,lambda0, ifNL, FFCID, S0

    integer :: i, j, m, jrcd
    integer locnum, LastMaxS11ElemID,CurrentMaxS11ElemID

    character*256 outdir,filepath,jobname
    integer lenoutdir,lenjobname,keyfile

    real*8 weibsum, lastweibsum, comTime, maxS11
    common /comm_Real/ weibsum, lastweibsum, comTime, maxS11
    common /comm_Integral/ LastMaxS11ElemID, CurrentMaxS11ElemID


    E11=props(1)
    E22=props(2)
    E33=props(3)
    v12=props(4)
    v13=props(5)
    v23=props(6)
    G12=props(7)
    G13=props(8)
    G23=props(9)
    Alpha1=props(10)
    Alpha2=props(11)
    Alpha3=props(12)
    S0=props(13)
    m=props(14)
    
    v21=v12*E22/E11
    v31=v13*E33/E11
    v32=v23*E33/E22
    
    tempVar=(1-v23*v32-v21*(v12+v13*v32)-v31*(v12*v23+v13))
    
    C11=(1-v23*v32)*E11/tempVar
    C22=(1-v13*v31)*E22/tempVar
    C33=(1-v12*v21)*E33/tempVar
    
    C12=(v21+v31*v23)*E11/tempVar
    C23=(v32+v12*v31)*E22/tempVar
    C13=(v13+v12*v23)*E33/tempVar
    
    C21=C12
    C32=C23
    C31=C13

    do i = 1, nblock

      ! Make history variable storage method as in Ls-Dyna
      do j=1,nstatev
        hsvNew(i,j)=hsvOld(i,j)
      end do

      ! save element ID into state variable 1 at the first time increment
      if(totalTime.le.dt) then
        call vgetpartinfo(nElement(i),1,cpname,locnum,jrcd)
        hsvNew(i,1)=locnum
        hsvNew(i,7)=zero
        hsvNew(i,6)=one
      end if

      ! Zero the weibull integral and the maximum S11 and save the maximum S11 element ID of the last step,
      !  at the first time when calling VUMAT in the current time step
      if (totalTime.gt.comTime) then
        lastweibsum=weibsum
        weibsum=zero
        maxS11=zero
        LastMaxS11ElemID=CurrentMaxS11ElemID
      end if
      ! assign the totalTime into a common block
      comTime=totalTime
      ! assign the last step element ID with the maximum S11
      hsvNew(i,3)=LastMaxS11ElemID
      ! assign the last step weibull sum to a history variable
      hsvNew(i,4)=lastweibsum

      !  save element volume into a history variable
      hsvNew(i,2)=charLength(i)**three

      ! If Weibull's criterion is satisfied, then set the damage variable
      if (lastweibsum.gt.one) then
        if (hsvNew(i,1).eq.LastMaxS11ElemID) then
          hsvNew(i,7)=three
        end if
      end if

      ! Update stresses
      if (hsvNew(i,7).eq.three) then

        stressNew(i,1)=zero
        stressNew(i,2)=zero
        stressNew(i,3)=zero
        stressNew(i,4)=zero
        stressNew(i,5)=zero
        stressNew(i,6)=zero
        hsvNew(i,6)=zero

      else

        dTemper=tempNew(i)-tempOld(i)
        dEpsilon1=strainInc(i,1)-Alpha1*dTemper
        dEpsilon2=strainInc(i,2)-Alpha2*dTemper
        dEpsilon3=strainInc(i,3)-Alpha3*dTemper

        stressNew(i,1)=stressOld(i,1)+C11*dEpsilon1+ &
                                 C12*dEpsilon2+C13*dEpsilon3
        stressNew(i,2)=stressOld(i,2)+C21*dEpsilon1+ &
                                 C22*dEpsilon2+C23*dEpsilon3
        stressNew(i,3)=stressOld(i,3)+C31*dEpsilon1+ &
                                 C32*dEpsilon2+C33*dEpsilon3
        stressNew(i,4)=stressOld(i,4)+two*G12*strainInc(i,4)
        stressNew(i,5)=stressOld(i,5)+two*G23*strainInc(i,5)
        stressNew(i,6)=stressOld(i,6)+two*G13*strainInc(i,6)

      end if

      ! Weibull failure criterion
      if (stressNew(i,1).gt.zero) then
        weibsum=weibsum+((stressNew(i,1)/S0)**m)*hsvNew(i,2)
      end if

      ! find out the element ID with the maximum S11 value
      if (stressNew(i,1).gt.maxS11) then
        maxS11=stressNew(i,1)
        CurrentMaxS11ElemID=hsvNew(i,1)
      end if

    end do

  end subroutine vumat_mlaminate

end module Weibull_Fibre_Failure
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132