!DIR$ FREEFORM

!>  Subroutines to compute the characteristic length by ERR/Ds gradient
!>
!>  Calculate the characteristic length of the element
!>
!>  Depends on the type of characteristic length required as specified by FCLType
!>
!>   FCLType:
!>     1 - Square root of the mid plane area
!>     2 - the length of the element
!>     3 - the width of the element
!>     otherwise - the 'dynamic' characteristic length.
!>
module Fatigue_CharacLength
  use Abaqus_Definitions, only: wp=>abaqus_real_kind
  use Mesh_Utils, only: cross_ab, dot_ab
  implicit none

contains

  !> Calculate the characteristic length for an element
  subroutine calcCharacteristicLength(hsvNew, FCLType, ElemType,FCL, &
    nblock, nstatev, &
    locnum, i)

    use Fatigue_Globals
    implicit none
    real(wp) hsvNew(nblock, nstatev)
    real*8 ElemType, IPCoord(4, 3), &
      ERR4IPs(4), CERR4IPs(4), FCL, FCLType
    integer i, nblock, nstatev, locnum

    if (FCLType.eq.1.d0) then
      FCL=sqrt(hsvNew(i,47)) !mid-plane area square root
    elseif (FCLType.eq.2.d0) then
      FCL=hsvNew(i,45) !length
    elseif (FCLType.eq.3.d0) then
      FCL=hsvNew(i,46) !width
    else
      ! If COH3D8 used, compute dynamic cha. L, only once for efficiency
      if (ElemType.ne.1.d0) then
        if (hsvNew(i,48).eq.0.d0) then  ! comment it if need update each dt
          IPCoord(1,1:3)=hsvNew(i,58:60)
          IPCoord(2,1:3)=hsvNew(i,61:63)
          IPCoord(3,1:3)=hsvNew(i,64:66)
          IPCoord(4,1:3)=hsvNew(i,67:69)
          ERR4IPs(1:4)=(/GV4IP(locnum,1),GV4IP(locnum,6), &
            GV4IP(locnum,11),GV4IP(locnum,16)/)
          CERR4IPs(1:4)=(/GV4IP(locnum,2),GV4IP(locnum,7), &
            GV4IP(locnum,12),GV4IP(locnum,17)/)
          call CharacLength(IPCoord,ERR4IPs,CERR4IPs,hsvNew(i,48))
        end if ! if hsvNew(i,48).eq.zero
      end if
      FCL=hsvNew(i,48)
    end if

  end subroutine calcCharacteristicLength

  
  subroutine CharacLength(IPCoord,ERRs,CERRs,FCL)

    double precision, parameter :: zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0, &
      third = one/three, four = 4.d0,  half = 0.5d0, sf=1.e-20

    real(wp) :: FCL
    real*8 IPCoord(4,3),IPLocalCoord(4,2),ERRs(4),CERRs(4),Xi,Eta, &
      dN1dXi,dN2dXi,dN3dXi,dN4dXi,dN1dEta,dN2dEta,dN3dEta,dN4dEta, &
      dN1dX,dN2dX,dN3dX,dN4dX,dN1dY,dN2dY,dN3dY,dN4dY, &
      A(2,4),B(4,2),e1(3),e2(3),e3(3),V12(3),V13(3),V14(3), &
      temp,CGV(2),JM(2,2),JMDet,JMR(2,2),EndNode(2,2),sgn(4), &
      ISPoint(4,2),LocalEC(2)
    integer i,j,k,ISNum,ifEdgeCut(4)

    ! Sort out the axis vectors of the coodinate system that has it x axis
    ! pointing from IP 1 to IP 2 and X-Y plane determined by the IPs 1, 2 and 4,
    !  below is sort out of x- and y- axis unit vecors e1 and e2.
    V12=(/IPCoord(2,1)-IPCoord(1,1),IPCoord(2,2)-IPCoord(1,2), &
      IPCoord(2,3)-IPCoord(1,3)/)
    V14=(/IPCoord(4,1)-IPCoord(1,1),IPCoord(4,2)-IPCoord(1,2), &
      IPCoord(4,3)-IPCoord(1,3)/)
    call cross_ab(V12,V14,e3);  call cross_ab(e3,V12,e2)
    e1=V12
    temp=dsqrt(e1(1)*e1(1)+e1(2)*e1(2)+e1(3)*e1(3))
    e1(1)=e1(1)/temp; e1(2)=e1(2)/temp; e1(3)=e1(3)/temp;
    temp=dsqrt(e2(1)*e2(1)+e2(2)*e2(2)+e2(3)*e2(3))
    e2(1)=e2(1)/temp; e2(2)=e2(2)/temp; e2(3)=e2(3)/temp;

    ! Coorindates of four IPs in the local coodinate system
    V13=(/IPCoord(3,1)-IPCoord(1,1),IPCoord(3,2)-IPCoord(1,2), &
      IPCoord(3,3)-IPCoord(1,3)/)
    IPLocalCoord(1,1:2)=(/zero,zero/)
    call dot_ab(e1,V12,temp); IPLocalCoord(2,1:2)=(/temp,0.d0/)
    call dot_ab(e1,V13,temp); IPLocalCoord(3,1)=temp
    call dot_ab(e2,V13,temp); IPLocalCoord(3,2)=temp
    call dot_ab(e1,V14,temp); IPLocalCoord(4,1)=temp
    call dot_ab(e2,V14,temp); IPLocalCoord(4,2)=temp

    ! Shape function partial derivative relative to natural axes at elem centroid
    Xi=0; Eta=0;
    dN1dXi=0.25*(Eta-1);  dN2dXi=0.25*(1-Eta);
    dN3dXi=0.25*(1+Eta);  dN4dXi=0.25*(-1-Eta);
    dN1dEta=0.25*(Xi-1);  dN2dEta=0.25*(-Xi-1);
    dN3dEta=0.25*(1+Xi);  dN4dEta=0.25*(1-Xi);

    ! Jacobian matrix JM at [Xi, Eta] and its inverse matrix
    A(1,1:4)=(/dN1dXi,dN2dXi,dN3dXi,dN4dXi/)
    A(2,1:4)=(/dN1dEta,dN2dEta,dN3dEta,dN4dEta/)
    B=IPLocalCoord
    JM(1,1)=A(1,1)*B(1,1)+A(1,2)*B(2,1)+A(1,3)*B(3,1)+A(1,4)*B(4,1)
    JM(1,2)=A(1,1)*B(1,2)+A(1,2)*B(2,2)+A(1,3)*B(3,2)+A(1,4)*B(4,2)
    JM(2,1)=A(2,1)*B(1,1)+A(2,2)*B(2,1)+A(2,3)*B(3,1)+A(2,4)*B(4,1)
    JM(2,2)=A(2,1)*B(1,2)+A(2,2)*B(2,2)+A(2,3)*B(3,2)+A(2,4)*B(4,2)
    JMDet=JM(1,1)*JM(2,2)-JM(2,1)*JM(1,2)+sf
    JMR(1,1)=JM(2,2)/JMDet; JMR(2,2)=JM(1,1)/JMDet;
    JMR(1,2)=-JM(1,2)/JMDet; JMR(2,1)=-JM(2,1)/JMDet;

    ! Shape function partial derivative relative to local coord. system e1 and e2
    dN1dX=JMR(1,1)*dN1dXi+JMR(1,2)*dN1dEta
    dN1dY=JMR(2,1)*dN1dXi+JMR(2,2)*dN1dEta
    dN2dX=JMR(1,1)*dN2dXi+JMR(1,2)*dN2dEta
    dN2dY=JMR(2,1)*dN2dXi+JMR(2,2)*dN2dEta
    dN3dX=JMR(1,1)*dN3dXi+JMR(1,2)*dN3dEta
    dN3dY=JMR(2,1)*dN3dXi+JMR(2,2)*dN3dEta
    dN4dX=JMR(1,1)*dN4dXi+JMR(1,2)*dN4dEta
    dN4dY=JMR(2,1)*dN4dXi+JMR(2,2)*dN4dEta

    ! Crack growth direction vector CGV(2)
    CGV(1)=dN1dX*ERRs(1)/CERRs(1)+dN2dX*ERRs(2)/CERRs(2)+ &
      dN3dX*ERRs(3)/CERRs(3)+dN4dX*ERRs(4)/CERRs(4)
    CGV(2)=dN1dY*ERRs(1)/CERRs(1)+dN2dY*ERRs(2)/CERRs(2)+ &
      dN3dY*ERRs(3)/CERRs(3)+dN4dY*ERRs(4)/CERRs(4)


    ! Characteristic length by the crack vector and element centroid
    ! first find the very far end points of the cutting line
    LocalEC(1)=(IPLocalCoord(1,1)+IPLocalCoord(2,1)+ &
      IPLocalCoord(3,1)+IPLocalCoord(4,1))/four
    LocalEC(2)=(IPLocalCoord(1,2)+IPLocalCoord(2,2)+ &
      IPLocalCoord(3,2)+IPLocalCoord(4,2))/four
    EndNode(1,1)=LocalEC(1)-2.E10*CGV(1)
    EndNode(1,2)=LocalEC(2)-2.E10*CGV(2)
    EndNode(2,1)=LocalEC(1)+2.E10*CGV(1)
    EndNode(2,2)=LocalEC(2)+2.E10*CGV(2)
    call SubPoin2Line(EndNode(1,1:2),EndNode(2,1:2), &
      IPLocalCoord(1,1:2), sgn(1))
    call SubPoin2Line(EndNode(1,1:2),EndNode(2,1:2), &
      IPLocalCoord(2,1:2), sgn(2))
    call SubPoin2Line(EndNode(1,1:2),EndNode(2,1:2), &
      IPLocalCoord(3,1:2), sgn(3))
    call SubPoin2Line(EndNode(1,1:2),EndNode(2,1:2), &
      IPLocalCoord(4,1:2), sgn(4))

    !  Check if the line intesects with the edge of nodes 1 and 2.
    do j=1,4; ISPoint(j,1:2)=(/1.E10,1.E10/); ifEdgeCut(j)=0; end do
    if ((sgn(1)*sgn(2)).le.zero) then ! intersect with edge of nodes 1 and 2
      ifEdgeCut(1)=1
      call IntersectPoint(IPLocalCoord(1,1:2),IPLocalCoord(2,1:2), &
        sgn(1),sgn(2),ISPoint(1,1:2))
    end if
    if ((sgn(2)*sgn(3)).le.zero) then ! intersect with edge of nodes 2 and 3
      ifEdgeCut(2)=1
      call IntersectPoint(IPLocalCoord(2,1:2),IPLocalCoord(3,1:2), &
        sgn(2),sgn(3),ISPoint(2,1:2))
    end if
    if ((sgn(3)*sgn(4)).le.zero) then ! intersect with edge of nodes 3 and 4
      ifEdgeCut(3)=1
      call IntersectPoint(IPLocalCoord(3,1:2),IPLocalCoord(4,1:2), &
        sgn(3),sgn(4),ISPoint(3,1:2))
    end if
    if ((sgn(1)*sgn(4)).le.zero) then ! intersect with edge of nodes 1 and 4
      ifEdgeCut(4)=1
      call IntersectPoint(IPLocalCoord(1,1:2),IPLocalCoord(4,1:2), &
        sgn(1),sgn(4),ISPoint(4,1:2))
    end if

    ! The characteristic length is assigned the maximum of the distances between
    ! all two inserction points, to consider there are more than two insersection
    ! points achived above in the case the cutting line exactly passes hrough one
    ! or two of the element corners
    FCL=zero
    do j=1,4
      do k=1,4
        if (j.ne.k) then
          if ((ifEdgeCut(j).eq.one).and.(ifEdgeCut(k).eq.one)) then
            temp=dsqrt((ISPoint(j,1)-ISPoint(k,1))**2+ &
              (ISPoint(j,2)-ISPoint(k,2))**2)
            FCL=max(FCL,temp)
          end if
        end if
      end do
    end do

  end subroutine CharacLength


  !>  Subroutines for computing the value by a linear equation
  subroutine SubPoin2Line(LineP1,LineP2,PointInQues,SubValue)

    real*8 LineP1(2),LineP2(2),PointInQues(2),SubValue

    SubValue=(PointInQues(2)-LineP2(2))*(LineP2(1)-LineP1(1))- &
      (PointInQues(1)-LineP2(1))*(LineP2(2)-LineP1(2))

  end subroutine SubPoin2Line


  !> Subroutine for computing insersection point on an element edge
  subroutine IntersectPoint(P1,P2,WF1,WF2,ISP)

    real*8  P1(2),P2(2),WF1,WF2,ISP(2),ratio1,ratio2

    ratio1=abs(WF1)/(abs(WF1)+abs(WF2))
    ratio2=abs(WF2)/(abs(WF1)+abs(WF2))

    ISP(1)=ratio2*P1(1)+ratio1*P2(1)
    ISP(2)=ratio2*P1(2)+ratio1*P2(2)

  end subroutine IntersectPoint

end module Fatigue_CharacLength
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132

