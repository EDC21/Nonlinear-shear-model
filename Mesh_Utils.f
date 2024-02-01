!DIR$ FREEFORM

!> 
module Mesh_Utils
  use iso_fortran_env, only: dp=>real64
  implicit none

contains

  !> Neighbour search and Computing Dimensions
  subroutine SNBCD(filepath,outdir,lenoutdir,maxCohElemNum,&
            maxNodeNum,Alpha_0,CohElem,IPEdgNbrElem,&
            IP1NbrElem,IP2NbrElem,IP3NbrElem,IP4NbrElem,&
            CohElemNum,ElemDimens,NodeCoords)

    character*256, intent(in) :: filepath
    character*256, intent(in) :: outdir
    integer, intent(in) :: lenoutdir
    integer, intent(in) :: maxCohElemNum
    integer, intent(in) :: maxNodeNum
    real*8, intent(in) :: Alpha_0
    integer :: CohElem(maxCohElemNum,9), &
                             IPEdgNbrElem(maxCohElemNum,9), &
                             IP1NbrElem(maxCohElemNum,12), &
                             IP2NbrElem(maxCohElemNum,12), &
                             IP3NbrElem(maxCohElemNum,12), &
                             IP4NbrElem(maxCohElemNum,12)
    integer, intent(out) :: CohElemNum
    real*8, intent(out) :: ElemDimens(maxCohElemNum,4), NodeCoords(maxNodeNum,3)

    real*8, parameter :: zero = 0.d0, half = 0.5d0, one = 1.d0, two = 2.d0,&
     three = 3.d0, four = 4.d0, &
     tol=1.e-6, scantol=2.e-2, GapToler=0.012
    integer, parameter :: NEPBX=2, NEPBY=2, NEPBZ=2,NumOfBands=1000
    ! NEPBX: Number of Elements Per Block in the global X axis.

    integer cohnodenum,tempID,tempIPEdgeNBNum,NSBX,NSBY,NSBZ,& !No. of Slave Block
     SBIDX,SBIDY,SBIDZ,&  ! Slave Block ID in X/Y/Z
     MBIDX2,MBIDY2,MBIDZ2,&  ! Master Block ID in X/Y/Z
     tempMBID_X(2),tempMBID_Y(2),tempMBID_Z(2),& !temporary master block ID
     MaxNumElemSB, MaxNumElemMB,& ! Max Number of Elements in one block
     tempIPEdgNB(8),BandID,ScanStartEID, CoinNodNum(8),&
     tempIP1NbrNum,tempIP2NbrNum,tempIP3NbrNum,tempIP4NbrNum
    real*8 temp(3),temp1,&
     Va(3),Vb(3),Vn(3),EBC(3),ETC(3),EBC2TC(3),&
     S15C(3),S26C(3),S37C(3),S48C(3),L,W,T,area,&
     VnDotVn,Alpha,&
     XLimits(2),YLimits(2),ZLimits(2),& !Min and max X/Y/Z of cohsive elems
     MaxESizeX,MaxESizeY,MaxESizeZ,& ! Maximum Element Size in X/Y/Z
     SBSizeX,SBSizeY,SBSizeZ,&  ! Slave Block Size in X/Y/Z
     tempEC(3),BandWid,&
     ScanLine,SndScanLine,temp2,temp3,IPEdgeNBFound,&
     IP1NBFound,IP2NBFound,IP3NBFound,IP4NBFound,IfIPNbr

    integer*8, allocatable :: SlaveBlock(:,:,:,:),& !Slave Blocks
     NumElemSB(:,:,:),&   !Number of Elements in each Slave Block
     MasterBlock(:,:,:,:),& !Master Blocks
     NumElemMB(:,:,:),&  !Number of Elements in each Master Block
     XBands(:,:),& !Save chopped bands in X direction
     YBands(:,:),& !Save chopped bands in Y direction
     ZBands(:,:),& !Save chopped bands in Z direction
     NumElem1Band(:,:) !Number of elements in each band and in X/Y/Z
    character*256 line,capitalline,&
     COH3D8Keyword,C3D8RKeyword,NodeKeyword

    integer keyfile,e,f,g,i,j,k,m,n,p,q,r

    integer ElemNumSorted
    integer, allocatable, dimension(:,:) :: IP1NbrElemLoID,&
    &IP2NbrElemLoID, IP3NbrElemLoID, IP4NbrElemLoID, IP1CoinElemIP,&
    &IP2CoinElemIP, IP3CoinElemIP, IP4CoinElemIP, ElemBlkIDs
    real*8, allocatable, dimension(:,:) :: EMinAxXYZ,&
    &SortedEMin, EIDAscMinXYZ, AllVn, EC

    ! ********************** Allocate temporary arrays **********************
    allocate(IP1NbrElemLoID(maxCohElemNum,11))
    allocate(IP2NbrElemLoID(maxCohElemNum,11))
    allocate(IP3NbrElemLoID(maxCohElemNum,11))
    allocate(IP4NbrElemLoID(maxCohElemNum,11))
    allocate(IP1CoinElemIP(maxCohElemNum,11))
    allocate(IP2CoinElemIP(maxCohElemNum,11))
    allocate(IP3CoinElemIP(maxCohElemNum,11))
    allocate(IP4CoinElemIP(maxCohElemNum,11))
    allocate(ElemBlkIDs(maxCohElemNum,6))

    allocate(EMinAxXYZ(maxCohElemNum,6))
    allocate(SortedEMin(maxCohElemNum,3))
    allocate(EIDAscMinXYZ(maxCohElemNum,3))
    allocate(AllVn(maxCohElemNum,3))
    allocate(EC(maxCohElemNum,3))

    !***************************************************************************
    !***********Open a file to record the running status of this subroutine*****
    !***************************************************************************
    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*)&
     'Start the neighbour searching & dimension calculation subroutine'
    write(keyfile,*) 'Start reading cohesive elements'
    close(keyfile)

    !Read out cohesive nodes and elements from the model .inp file
    block

      real(dp), allocatable :: NodeCoords_dyn(:,:)
      integer, allocatable :: CohElem_dyn(:,:)
      integer :: nnode

      call ReadCohesiveMesh(filepath,NodeCoords_dyn,CohElem_dyn)
      
      nnode = size(NodeCoords_dyn,1)
      CohElemNum = size(CohElem_dyn,1)

      NodeCoords(1:nnode,:) = NodeCoords_dyn
      CohElem(1:CohElemNum,:) = CohElem_dyn

    end block

!
    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'Finished reading nodes'
    write(keyfile,*)  'The total number of cohesive elements is: ',&
     CohElemNum
    write(keyfile,*) 'Finished reading all nodes & cohesive elements.'
    write(keyfile,*)&
     'Start computing element dimensions and the global limits etc'
    close(keyfile)
!
    !Sort out the centres, local axes and dimensions of elements
    call ElemCentreAxesDimens(NodeCoords,CohElemNum,CohElem,&
     MaxESizeX,MaxESizeY,MaxESizeZ,EC,AllVn,ElemDimens)
!
    ! Sort out the X/Y/Z limits for the whole cohesive zone,
    ! and the X/Y/Z limits for each cohesive element
    call CohZoneElemXYZLimits(NodeCoords,CohElemNum,CohElem,&
     XLimits,YLimits,ZLimits,EMinAxXYZ)
!


    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*)&
     'Finished computing element dimensions and the global limits etc'
    write(keyfile,*) 'The max element sizes in the global X, Y, Z are:',&
     MaxESizeX, MaxESizeY, MaxESizeZ
    write(keyfile,*) 'The cohesive-zone lower and upper X limits are:',&
     XLimits(1), XLimits(2)
    write(keyfile,*) 'The cohesive-zone lower and upper Y limits are:',&
     YLimits(1), YLimits(2)
    write(keyfile,*) 'The cohesive-zone lower and upper Z limits are:',&
     ZLimits(1), ZLimits(2)
    write(keyfile,*) 'Start model chopping and element allocation'
    write(keyfile,*) 'into a slave block and at least one master block'
    close(keyfile)

    !****************************************************************************
    !To speed up neighbour search, the whole cohesive region is chopped to small
    !blocks by three scanning lines in global X, Y and Z directions. Scan starts
    !from the minimum coordinate to the maximum coordinate of the whole cohesive
    !zone. Each scanning loop follows these steps: 1)slightly move the scan
    !line by a very small value 1e-6 towards the positive X/Y/Z direction,
    !2)find all the elements crossed by the scan line, 3)move the scanning line
    !to the maximum cooridates of all the crossed elements, 4)find all the elems
    !on the nagnetive-axis side of the scan line, but not in the previous bands.
    !
    !Neighbour search is done for each element within its affiliation block and
    ! adjacent blocks(3*9 blocks at maximum)
    !**************************************************************************

    !First sort the cohesive elements in the scending order by their 'minimum'
    !coordinates in global X,Y and Z directions. To speed up, the sorting
    !is also done by chopping the whole cohesive zone into multiple bands.
    ! The sorted element IDs are saved in the variable 'EIDAscMinXYZ'


    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*)&
     'Start sorting elements in the acending X/Y/Z coorinates'
    close(keyfile)
    !Sort elements in the ascending order of the minimum X coordinates of elements
    call SortElementAscendingMinCoord(CohElemNum,XLimits,&
     EMinAxXYZ(1:maxCohElemNum,1),EIDAscMinXYZ(1:maxCohElemNum,1),&
     ElemNumSorted)

    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'Total sorted element number in X:', ElemNumSorted
    close(keyfile)

    !Sort elements in the ascending order of the minimum Y coordinates of elements
    call SortElementAscendingMinCoord(CohElemNum,YLimits,&
     EMinAxXYZ(1:maxCohElemNum,3),EIDAscMinXYZ(1:maxCohElemNum,2),&
     ElemNumSorted)

    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'Total sorted element number in Y:', ElemNumSorted
    close(keyfile)


    call SortElementAscendingMinCoord(CohElemNum,ZLimits,&
     EMinAxXYZ(1:maxCohElemNum,5),EIDAscMinXYZ(1:maxCohElemNum,3),&
     ElemNumSorted)

    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'Total sorted element number in Z:', ElemNumSorted
    write(keyfile,*)&
     'Finished sorting elements in the acending X/Y/Z coorinates'
    close(keyfile)

    !Scan and save element affiliation block IDs in the X direction
    NSBX=1
    ScanLine=XLimits(1)
    ScanStartEID=1
4000 SndScanLine=ScanLine+scantol
    !Check out all the elements with their minimum coordinate points within
    !the band between the scan line and the secondary scan line
    do j=ScanStartEID,CohElemNum ! scan from left to right based on the min X
      if (EMinAxXYZ(EIDAscMinXYZ(j,1),1).le.SndScanLine) then
        ElemBlkIDs(EIDAscMinXYZ(j,1),1)=NSBX! save X block ID
        ScanLine=max(ScanLine,EMinAxXYZ(EIDAscMinXYZ(j,1),2))
      else
        go to 4020
      end if
    end do
4020 ScanStartEID=j
    !Check out all the elements with their minimum coordinate points within
    !the band between the updated scan line and the secondary scan line
    temp1=ScanLine ! update the scanline for next searching loop
    do j=ScanStartEID,CohElemNum ! scan from left to right based on the min X
      if (EMinAxXYZ(EIDAscMinXYZ(j,1),1).le.(ScanLine-tol)) then
        ElemBlkIDs(EIDAscMinXYZ(j,1),1)=NSBX! save X block ID
        temp1=max(temp1,EMinAxXYZ(EIDAscMinXYZ(j,1),2))
      else
        go to 4040
      end if
    end do
4040 ScanStartEID=j
    ScanLine=temp1
    if (ScanLine.ge.(XLimits(2)-tol)) then
      go to 4050
    else
      if (ScanLine.gt.SndScanLine) then
        NSBX=NSBX+1
      else
        ScanLine=SndScanLine
      end if
      go to 4000
    end if
4050 continue

    !Scan and save element affiliation block IDs in the Y direction
    NSBY=1
    ScanLine=YLimits(1)
    ScanStartEID=1
5000 SndScanLine=ScanLine+scantol
    !Check out all the elements that overlap with the band between
    !the secondary scan line and the previous scan line
    do j=ScanStartEID,CohElemNum ! scan from left to right based on the min Y
      if (EMinAxXYZ(EIDAscMinXYZ(j,2),3).le.SndScanLine) then
        ElemBlkIDs(EIDAscMinXYZ(j,2),2)=NSBY! save Y block ID
        ScanLine=max(ScanLine,EMinAxXYZ(EIDAscMinXYZ(j,2),4))
      else
        go to 5020
      end if
    end do
5020 ScanStartEID=j
    !Check out all the elements that overlap with the band between the secondary
    !and the current scan line
    temp1=ScanLine ! update the scanline for next searching loop
    do j=ScanStartEID,CohElemNum ! scan from left to right based on the min Y
      if (EMinAxXYZ(EIDAscMinXYZ(j,2),3).le.(ScanLine-tol)) then
        ElemBlkIDs(EIDAscMinXYZ(j,2),2)=NSBY! save Y block ID
        temp1=max(temp1,EMinAxXYZ(EIDAscMinXYZ(j,2),4))
      else
        go to 5040
      end if
    end do
5040 ScanStartEID=j
    ScanLine=temp1
    if (ScanLine.ge.(YLimits(2)-tol)) then
      go to 5050
    else
      if (ScanLine.gt.SndScanLine) then
        NSBY=NSBY+1
      else
        ScanLine=SndScanLine
      end if
      go to 5000
    end if
5050 continue

    !Scan and save element affiliation block IDs in the Z direction
    NSBZ=1
    ScanLine=ZLimits(1)
    ScanStartEID=1
6000 SndScanLine=ScanLine+scantol
    !Check out all the elements that overlap with the band between
    !the secondary scan line and the previous scan line
    do j=ScanStartEID,CohElemNum ! scan from left to right based on the min Z
      if (EMinAxXYZ(EIDAscMinXYZ(j,3),5).le.SndScanLine) then
        ElemBlkIDs(EIDAscMinXYZ(j,3),3)=NSBZ! save Z block ID
        ScanLine=max(ScanLine,EMinAxXYZ(EIDAscMinXYZ(j,3),6))
      else
        go to 6020
      end if
    end do
6020 ScanStartEID=j
    !Check out all the elements that overlap with the band between the secondary
    !and the current scan line
    temp1=ScanLine ! update the scanline for next searching loop
    do j=ScanStartEID,CohElemNum ! scan from left to right based on the min Z
      if (EMinAxXYZ(EIDAscMinXYZ(j,3),5).le.(ScanLine-tol)) then
        ElemBlkIDs(EIDAscMinXYZ(j,3),3)=NSBZ! save Z block ID into
        temp1=max(temp1,EMinAxXYZ(EIDAscMinXYZ(j,3),6))
      else
        go to 6040
      end if
    end do
6040 ScanStartEID=j
    ScanLine=temp1
    if (ScanLine.ge.(ZLimits(2)-tol)) then
      go to 6050
    else
      if (ScanLine.gt.SndScanLine) then
        NSBZ=NSBZ+1
      else
        ScanLine=SndScanLine
      end if
      go to 6000
    end if
6050 continue

    ! allocate and initiliase the variables used to count the numbers of elements
    ! for each slave/master blocks
    allocate(NumElemSB(NSBX,NSBY,NSBZ))
    do i=1,NSBX
      do j=1,NSBY
        do k=1,NSBZ
          NumElemSB(i,j,k)=0
        end do
      end do
    end do

    !Loop within all cohesive elements to sort out the sizes of all blocks
    !and the maximum number of elements over all blocks
    MaxNumElemSB=0
    do i=1,CohElemNum
      SBIDX=ElemBlkIDs(i,1);SBIDY=ElemBlkIDs(i,2);SBIDZ=ElemBlkIDs(i,3)
      NumElemSB(SBIDX,SBIDY,SBIDZ)=NumElemSB(SBIDX,SBIDY,SBIDZ)+1
      MaxNumElemSB=max(MaxNumElemSB,NumElemSB(SBIDX,SBIDY,SBIDZ))
    end do

    ! allocate element into a slave block
    allocate(SlaveBlock(NSBX,NSBY,NSBZ,MaxNumElemSB))
    do i=1,NSBX
      do j=1,NSBY
        do k=1,NSBZ
          NumElemSB(i,j,k)=0
        end do
      end do
    end do

    do i=1,CohElemNum
      NumElemSB(ElemBlkIDs(i,1),ElemBlkIDs(i,2),ElemBlkIDs(i,3))=&
      &NumElemSB(ElemBlkIDs(i,1),ElemBlkIDs(i,2),ElemBlkIDs(i,3))+1
      SlaveBlock(ElemBlkIDs(i,1),ElemBlkIDs(i,2),ElemBlkIDs(i,3),&
      &NumElemSB(ElemBlkIDs(i,1),ElemBlkIDs(i,2),ElemBlkIDs(i,3)))=i
    end do

    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'Finished model chopping and element allocation'
    write(keyfile,*) 'The numbers of blocks in X,Y and Z are:',&
     NSBX, NSBY, NSBZ
    write(keyfile,*) 'The max number of elements in 1 block is:',&
     MaxNumElemSB
    write(keyfile,*) 'Start neighbour seraching'
    close(keyfile)

    ! Start search neighbour for each elemnt and IP
    do i=1,NSBX
      do j=1,NSBY
        do k=1,NSBZ
          do m=1,NumElemSB(i,j,k)
            p=SlaveBlock(i,j,k,m) ! local ID of the current element in [1,CohElemNum]
            tempIPEdgNB(1:8)=(/0,0,0,0,0,0,0,0/)
            IP1NbrElemLoID(p,1:11)=(/0,0,0,0,0,0,0,0,0,0,0/)
            IP2NbrElemLoID(p,1:11)=(/0,0,0,0,0,0,0,0,0,0,0/)
            IP3NbrElemLoID(p,1:11)=(/0,0,0,0,0,0,0,0,0,0,0/)
            IP4NbrElemLoID(p,1:11)=(/0,0,0,0,0,0,0,0,0,0,0/)
            tempIPEdgeNBNum=0
            tempIP1NbrNum=0
            tempIP2NbrNum=0
            tempIP3NbrNum=0
            tempIP4NbrNum=0
            IPEdgeNBFound=zero
            IP1NBFound=zero
            IP2NBFound=zero
            IP3NBFound=zero
            IP4NBFound=zero
            do e=i-1,i+1 ! loop in all adjacent blocks if any
              do f=j-1,j+1
                do g=k-1,k+1
                  if ((e.ge.one).and.(e.le.NSBX).and.(f.ge.one).and.(f.le.NSBY)&
                      .and.(g.ge.one).and.(g.le.NSBZ)) then
                    do n=1,NumElemSB(e,f,g)
                      q=SlaveBlock(e,f,g,n) ! local element ID in [1,CohElemNum]
                      if (p.ne.q) then  ! to avoid detecting itself as a neighbour
                        call dot_ab(AllVn(p,1:3),AllVn(q,1:3),VnDotVn)
                        VnDotVn=abs(VnDotVn)
                        if (VnDotVn.gt.one) then !acos gives 'NaN' if entry larger than 1
                          Alpha=90.d0
                        else
                          Alpha=180.0*acos(VnDotVn)/3.1415926
                        end if ! Alpha is actually the 'minimum' angle between the normals
                        ! distinguish between inter or intra laminar element
                        if (Alpha.lt.Alpha_0) then
                          ! search in-plane edge neighbours on all possbile connections
                          if (IPEdgeNBFound.ne.one) then
                            call WhetherConnectInPlaneEdge(IPEdgeNBFound,&
                             CoinNodNum,CohElem(p,1:9),CohElem(q,1:9),&
                             tempIPEdgeNBNum,tempIPEdgNB)
                          end if

                          ! search IP 1 host vertical edge adjacent neighbours
                          IfIPNbr=zero
                          if (IP1NBFound.ne.one) then
                            call IPConnectedElement(IP1NBFound,IfIPNbr,&
                             CohElem(p,2),CohElem(p,6),q,CohElem(q,1:9),&
                             NodeCoords,tempIP1NbrNum,&
                             IP1NbrElemLoID(p,1:11),IP1CoinElemIP(p,1:11))
                          end if

                          ! search IP 2 host vertical edge adjacent neighbours
                          if (IP2NBFound.ne.one) then
                            call IPConnectedElement(IP2NBFound,IfIPNbr,&
                             CohElem(p,3),CohElem(p,7),q,CohElem(q,1:9),&
                             NodeCoords,tempIP2NbrNum,&
                             IP2NbrElemLoID(p,1:11),IP2CoinElemIP(p,1:11))
                          end if
!
                          ! search IP 3 host vertical edge adjacent neighbours
                          if (IP3NBFound.ne.one) then
                            call IPConnectedElement(IP3NBFound,IfIPNbr,&
                             CohElem(p,4),CohElem(p,8),q,CohElem(q,1:9),&
                             NodeCoords,tempIP3NbrNum,&
                             IP3NbrElemLoID(p,1:11),IP3CoinElemIP(p,1:11))
                          end if

                          ! search IP 4 host vertical edge adjacent neighbours
                          if (IP4NBFound.ne.one) then
                            call IPConnectedElement(IP4NBFound,IfIPNbr,&
                             CohElem(p,5),CohElem(p,9),q,CohElem(q,1:9),&
                             NodeCoords,tempIP4NbrNum,&
                             IP4NbrElemLoID(p,1:11),IP4CoinElemIP(p,1:11))
                          end if

                        end if ! end of 'if (Alpha.lt.Alpha_0) then'
                      end if ! end of if (p.ne.q) then
                    end do ! end of do n=1,ElemNum1MC

                  end if ! end of 'if((p.ge.1).and.(p.le.NSBX)...'
                end do !end of 'do e=k-1,k+1'
              end do !end of 'do f=j-1,j+1'
            end do !end of 'do g=i-1,i+1'
8000        continue
            ! Save current element in-plane edge neighbours into whole matrix
            IPEdgNbrElem(p,2:9)=tempIPEdgNB(1:8)
            IPEdgNbrElem(p,1)=CohElem(p,1)
            ! Save current element integration point neighbours into whole matrix
            do e=1,11
              if (IP1NbrElemLoID(p,e).ne.zero) then
                IP1NbrElem(p,e)=CohElem(IP1NbrElemLoID(p,e),1)
              else
                IP1NbrElem(p,e)=zero
              end if
              if (IP2NbrElemLoID(p,e).ne.zero) then
                IP2NbrElem(p,e)=CohElem(IP2NbrElemLoID(p,e),1)
              else
                IP2NbrElem(p,e)=zero
              end if
              if (IP3NbrElemLoID(p,e).ne.zero) then
                IP3NbrElem(p,e)=CohElem(IP3NbrElemLoID(p,e),1)
              else
                IP3NbrElem(p,e)=zero
              end if
              if (IP4NbrElemLoID(p,e).ne.zero) then
                IP4NbrElem(p,e)=CohElem(IP4NbrElemLoID(p,e),1)
              else
                IP4NbrElem(p,e)=zero
              end if
            end do
          end do ! end of 'do m=1,NumElemSB(i,j,k)'
        end do ! end of 'do k=1,NSBZ'
      end do ! end of 'do j=1,NSBY'
    end do	! end of 'do i=1,NSBX'

    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'Finished neighbour seraching'
    write(keyfile,*) 'Start checkout IPs at discontinuous locations'
    close(keyfile)


    !Judge if elements' IPs are at an edge or discontinuous location of the model
    !by checking if the corresponding in-plane edges of its host elements and
    !neighbouring elements are connected to an element.
    do i=1,NSBX
      do j=1,NSBY
        do k=1,NSBZ
          do m=1,NumElemSB(i,j,k)
            p=SlaveBlock(i,j,k,m) ! local ID of the current element in [1,CohElemNum]
            IP1NbrElem(p,12)=zero
            IP2NbrElem(p,12)=zero
            IP3NbrElem(p,12)=zero
            IP4NbrElem(p,12)=zero

            !Check whether IP 1 is an IP at edge or discontinuous location
            call WhetherEdgeIP(p,one,IPEdgNbrElem,IP1NbrElemLoID(p,1:11),&
             IP1CoinElemIP(p,1:11),IP1NbrElem(p,12))
!
            !Check whether IP 2 is an IP at edge or discontinuous location
            call WhetherEdgeIP(p,two,IPEdgNbrElem,IP2NbrElemLoID(p,1:11),&
             IP2CoinElemIP(p,1:11),IP2NbrElem(p,12))
!
            !Check whether IP 3 is an IP at edge or discontinuous location
            call WhetherEdgeIP(p,three,IPEdgNbrElem,IP3NbrElemLoID(p,1:11),&
             IP3CoinElemIP(p,1:11),IP3NbrElem(p,12))
!
            !Check whether IP 4 is an IP at edge or discontinuous location
            call WhetherEdgeIP(p,four,IPEdgNbrElem,IP4NbrElemLoID(p,1:11),&
             IP4CoinElemIP(p,1:11),IP4NbrElem(p,12))

          end do ! end of 'do m=1,NumElemSB(i,j,k)'
        end do ! end of 'do k=1,NSBZ'
      end do ! end of 'do j=1,NSBY'
    end do	! end of 'do i=1,NSBX'

    open (newunit = keyfile, file=outdir(1:lenoutdir)//'/'&
     //'Status_SNBCD.txt',position='append',status='unknown')
    write(keyfile,*) 'All DONE!'
    close(keyfile)

  end subroutine SNBCD

  
  !> Subroutine to read cohesive nodes and elements from input deck (.inp)
  subroutine ReadCohesiveMesh(inp_file,NodeCoords,CohElem)

    !> File path to input deck file (.inp)
    character(*), intent(in) :: inp_file

    !> Matrix of node coordinates (nnode x 3)
    real(dp), intent(out), allocatable :: NodeCoords(:,:)

    !> Matrix of element connectivities (nelem x 9)
    integer, intent(out), allocatable :: CohElem(:,:)

    !> Local vars
    integer :: fh, stat, pass, line_len
    integer :: nnode, nelem
    integer, parameter :: MAX_LINE_LENGTH = 1024
    character(len=MAX_LINE_LENGTH) :: line
    integer :: temp_int(9)
    real(dp) :: temp_real(3)

    !> Keyword identifiers
    character(*), parameter :: NODEKEYWORD='*NODE, NSET=SEARCHREQUIRED'
    character(*), parameter :: COH3D8Keyword='*ELEMENT, TYPE=COH3D8, ELSET=INTERFACE'
    character(*), parameter :: C3D8RKeyword ='*ELEMENT, TYPE=C3D8R, ELSET=INTERFACE'

    do pass=1,2

      open(newunit=fh,file=trim(inp_file),status='old')

      nnode = 0
      nelem = 0

      do

        read(fh,'(A)',iostat=stat) line

        if (stat /= 0) then
          exit
        end if

        line_len = len_trim(line)
        line(1:line_len) = to_upper(line(1:line_len))

        if (line(1:26) == NODEKEYWORD) then

          ! Read cohesive nodes
          do

            read(fh,*,iostat=stat) temp_int(1), temp_real(1), temp_real(2), temp_real(3)

            if (stat /= 0) then
              backspace(fh)
              exit
            end if

            nnode = nnode + 1

            if (pass == 2) then

              NodeCoords(nnode,1:3) = temp_real(1:3)

            end if

          end do ! read cohesive nodes

        end if ! node keyword

        ! Read cohesive elements
        if (line(1:38) == COH3D8Keyword .or. &
            line(1:37) == C3D8RKeyword) then

          do

            read(fh,*,iostat=stat) temp_int(1:9)

            if (stat /= 0) then
              backspace(fh)
              exit
            end if

            nelem = nelem + 1

            if (pass == 2) then

              CohElem(nelem,1:9) = temp_int(1:9)

            end if

          end do ! read cohesive elements

        end if ! elements keyword

      end do ! read line loop

      close(fh)

      if (pass == 1) then
        allocate(NodeCoords(nnode,3))
        allocate(CohElem(nelem,9))
      end if

    end do ! two pass loop

  end subroutine ReadCohesiveMesh


  !> Compute centres/axes/dimensions of cohesive elements
  subroutine ElemCentreAxesDimens(NodeCoords,CohElemNum,CohElem,&
                  MaxESizeX,MaxESizeY,MaxESizeZ,EC,AllVn,ElemDimens)
    use Fatigue_Globals, only: maxCohElemNum, maxNodeNum
    real*8, parameter :: half = 0.5d0, two = 2.d0, four = 4.d0
    integer CohElem(maxCohElemNum,9),CohElemNum,n,j,k
    real*8  temp(3),temp1,AllVn(maxCohElemNum,3),&
     Va(3),Vb(3),Vn(3),EBC(3),ETC(3),EBC2TC(3),&
     S15C(3),S26C(3),S37C(3),S48C(3),L,W,T,area,&
     NodeCoords(maxNodeNum,3),EC(maxCohElemNum,3),&
     ElemDimens(maxCohElemNum,4),tempEC(3),&
     MaxESizeX,MaxESizeY,MaxESizeZ


    MaxESizeX=0; MaxESizeY=0; MaxESizeZ=0;
    do j=1,CohElemNum  ! loop in all cohesive elements
      do n=1,3
        S15C(n)=half*(NodeCoords(CohElem(j,2),n)+&
         NodeCoords(CohElem(j,6),n))
        S26C(n)=half*(NodeCoords(CohElem(j,3),n)+&
         NodeCoords(CohElem(j,7),n))
        S37C(n)=half*(NodeCoords(CohElem(j,4),n)+&
         NodeCoords(CohElem(j,8),n))
        S48C(n)=half*(NodeCoords(CohElem(j,5),n)+&
         NodeCoords(CohElem(j,9),n))
        Va(n)=(S37C(n)+S48C(n))/two-(S15C(n)+S26C(n))/two
        Vb(n)=(S26C(n)+S37C(n))/two-(S15C(n)+S48C(n))/two
        EBC(n)=(NodeCoords(CohElem(j,2),n)+NodeCoords(CohElem(j,3),n)&
         +NodeCoords(CohElem(j,4),n)+NodeCoords(CohElem(j,5),n))/four
        ETC(n)=(NodeCoords(CohElem(j,6),n)+NodeCoords(CohElem(j,7),n)&
         +NodeCoords(CohElem(j,8),n)+NodeCoords(CohElem(j,9),n))/four
        EBC2TC(n)=ETC(n)-EBC(n)
        tempEC(n)=(ETC(n)+EBC(n))/two
      end do
      ! Update maximum element dimensions in global X, Y and Z axes.
      MaxESizeX=max(MaxESizeX,abs(Va(1)),abs(Vb(1)),abs(EBC2TC(1)))
      MaxESizeY=max(MaxESizeY,abs(Va(2)),abs(Vb(2)),abs(EBC2TC(2)))
      MaxESizeZ=max(MaxESizeZ,abs(Va(3)),abs(Vb(3)),abs(EBC2TC(3)))
      ! Element length, width and mid cross-section area
      call cross_ab(Va,Vb,Vn)
      area=(Vn(1)*Vn(1)+Vn(2)*Vn(2)+Vn(3)*Vn(3))**half
      L=   (Va(1)*Va(1)+Va(2)*Va(2)+Va(3)*Va(3))**half
      W=   (Vb(1)*Vb(1)+Vb(2)*Vb(2)+Vb(3)*Vb(3))**half
      if (L.lt.W) then
        temp1=W
        W=L
        L=temp1
        temp=Va
        Va=Vb
        Vb=temp
      end if
      ! Normalised local vectors of Va and Vn
      do k=1,3
        Va(k)=Va(k)/L
        Vn(k)=Vn(k)/area
      end do
      ! recompute unit Vb using Va and Vn
      call cross_ab(Vn,Va,Vb)
      ! Save the unit normal vectors
      AllVn(j,1:3)=Vn(1:3)
      ! element thickness
      T=abs(Vn(1)*EBC2TC(1)+Vn(2)*EBC2TC(2)+Vn(3)*EBC2TC(3))
      ! save unit local vectors, element centres and element dimensions
      EC(j,1:3)=tempEC
      ElemDimens(j,1)=T
      ElemDimens(j,2)=L
      ElemDimens(j,3)=W
      ElemDimens(j,4)=area

    end do ! end of do j=1,CohElemNum

  end subroutine ElemCentreAxesDimens


  !> Compute the X/Y/Z limits for the whole cohesive zone,
  !>   and the X/Y/Z limits for each cohesive element
  subroutine CohZoneElemXYZLimits(NodeCoords,CohElemNum,CohElem,&
                    XLimits,YLimits,ZLimits,EMinAxXYZ)
    use Fatigue_Globals, only: maxCohElemNum, maxNodeNum

    integer CohElem(maxCohElemNum,9),CohElemNum,j,k,tempID
    real*8  XLimits(2),YLimits(2),ZLimits(2),&
     NodeCoords(maxNodeNum,3),EMinAxXYZ(maxCohElemNum,6)

    ! Initialise variables
    XLimits(1)=NodeCoords(CohElem(1,2),1); XLimits(2)=XLimits(1);
    YLimits(1)=NodeCoords(CohElem(1,2),2); YLimits(2)=YLimits(1);
    ZLimits(1)=NodeCoords(CohElem(1,2),3); ZLimits(2)=ZLimits(1);

    do j=1, CohElemNum

      ! Update global min and max X/Y/Z of the cohesive element domain
      do k =2,9
        tempID=CohElem(j,k)
        XLimits(1)=min(XLimits(1),NodeCoords(tempID,1))
        XLimits(2)=max(XLimits(2),NodeCoords(tempID,1))
        YLimits(1)=min(YLimits(1),NodeCoords(tempID,2))
        YLimits(2)=max(YLimits(2),NodeCoords(tempID,2))
        ZLimits(1)=min(ZLimits(1),NodeCoords(tempID,3))
        ZLimits(2)=max(ZLimits(2),NodeCoords(tempID,3))
      end do

      ! Update the local min and max X/Y/Z of each cohesive element
      tempID=CohElem(j,2)
      EMinAxXYZ(j,1)=NodeCoords(tempID,1) ! X min
      EMinAxXYZ(j,2)=NodeCoords(tempID,1) ! X max
      EMinAxXYZ(j,3)=NodeCoords(tempID,2) ! Y min
      EMinAxXYZ(j,4)=NodeCoords(tempID,2) ! Y max
      EMinAxXYZ(j,5)=NodeCoords(tempID,3) ! Z min
      EMinAxXYZ(j,6)=NodeCoords(tempID,3) ! Z max
      do k=3,9
        tempID=CohElem(j,k)
        EMinAxXYZ(j,1)=min(EMinAxXYZ(j,1),NodeCoords(tempID,1))
        EMinAxXYZ(j,2)=max(EMinAxXYZ(j,2),NodeCoords(tempID,1))
        EMinAxXYZ(j,3)=min(EMinAxXYZ(j,3),NodeCoords(tempID,2))
        EMinAxXYZ(j,4)=max(EMinAxXYZ(j,4),NodeCoords(tempID,2))
        EMinAxXYZ(j,5)=min(EMinAxXYZ(j,5),NodeCoords(tempID,3))
        EMinAxXYZ(j,6)=max(EMinAxXYZ(j,6),NodeCoords(tempID,3))
      end do

    end do ! end of do j=1,CohElemNum

  end subroutine CohZoneElemXYZLimits


  !> Sort elements in the ascending order of
  !>   the minimum X coordinates of elements
  subroutine SortElementAscendingMinCoord(CohElemNum,Limits,&
               EMinCoord,EIDAscMinCoord,ElemNumSorted)
    use Fatigue_Globals, only: maxCohElemNum
    real*8, parameter :: zero=0.0d0, one = 1.d0, tol=1.e-6
    integer, parameter :: NumOfBands=1000
    integer CohElemNum,j,k,m,n,p,BandID,ElemNumSorted
    real*8  Limits(2),EMinCoord(maxCohElemNum),&
     EIDAscMinCoord(maxCohElemNum),temp2,temp3, BandWid

    integer*8, allocatable :: Bands(:,:),& !Save chopped bands in X direction
     NumElem1Band(:) !Number of elements in each band and in X/Y/Z

    !Create storage to save the No. of elements for each band
    allocate(NumElem1Band(NumOfBands))
    do j=1,NumOfBands; NumElem1Band(j)=0; end do

    !Chop into bands and allocate elements into bands
    allocate(Bands(NumOfBands,maxCohElemNum/10))
    BandWid=(Limits(2)-Limits(1))/NumOfBands
    do j=1,CohElemNum
      BandID=ceiling((EMinCoord(j)+tol-Limits(1))/BandWid)
      NumElem1Band(BandID)=NumElem1Band(BandID)+1
      Bands(BandID,NumElem1Band(BandID))=j
    end do
    !Loop over bands, sort elements in the ascending order of their mini coord.
    ElemNumSorted=0
    do j=1,NumOfBands
      if (NumElem1Band(j).gt.one) then !check if an empty band
        do k=1,(NumElem1Band(j)-1) !
          temp2=EMinCoord(Bands(j,k))
          do m=k+1,NumElem1Band(j)
            temp3=EMinCoord(Bands(j,m))
            if (temp3.lt.temp2) then
              n=Bands(j,k)
              Bands(j,k)=Bands(j,m)
              Bands(j,m)=n
            end if
          end do
        end do
        do p=1,NumElem1Band(j)
          EIDAscMinCoord(ElemNumSorted+p)=Bands(j,p)
        end do
        ElemNumSorted=ElemNumSorted+NumElem1Band(j)
      elseif (NumElem1Band(j).gt.zero) then !only 1 elem. in current band
        EIDAscMinCoord(ElemNumSorted+1)=Bands(j,1)
        ElemNumSorted=ElemNumSorted+1
      end if
    end do

  end subroutine SortElementAscendingMinCoord


  !> Check whether element q connects at least with one
  !>   of the 8 in-plane edges of element p.
  subroutine WhetherConnectInPlaneEdge(IPEdgeNBFound,CoinNodNum,&
            ElemP,ElemQ,tempIPEdgeNBNum,tempIPEdgNB)

    integer ElemP(9),ElemQ(9),r,CoinNodNum(8),tempIPEdgeNBNum,&
     tempIPEdgNB(8)
    real*8 IPEdgeNBFound

    CoinNodNum=(/0,0,0,0,0,0,0,0/)
    do r=2,9
      if ((ElemP(2).eq.ElemQ(r)).or.&
       (ElemP(3).eq.ElemQ(r)))then
        CoinNodNum(1)=CoinNodNum(1)+1
      end if
      if ((ElemP(4).eq.ElemQ(r)).or.&
       (ElemP(5).eq.ElemQ(r))) then
        CoinNodNum(2)=CoinNodNum(2)+1
      end if
      if ((ElemP(3).eq.ElemQ(r)).or.&
       (ElemP(4).eq.ElemQ(r))) then
        CoinNodNum(3)=CoinNodNum(3)+1
      end if
      if ((ElemP(2).eq.ElemQ(r)).or.&
       (ElemP(5).eq.ElemQ(r))) then
        CoinNodNum(4)=CoinNodNum(4)+1
      end if
      if ((ElemP(6).eq.ElemQ(r)).or.&
       (ElemP(7).eq.ElemQ(r))) then
        CoinNodNum(5)=CoinNodNum(5)+1
      end if
      if ((ElemP(8).eq.ElemQ(r)).or.&
       (ElemP(9).eq.ElemQ(r))) then
        CoinNodNum(6)=CoinNodNum(6)+1
      end if
      if ((ElemP(7).eq.ElemQ(r)).or.&
       (ElemP(8).eq.ElemQ(r))) then
        CoinNodNum(7)=CoinNodNum(7)+1
      end if
      if ((ElemP(6).eq.ElemQ(r)).or.&
       (ElemP(9).eq.ElemQ(r))) then
        CoinNodNum(8)=CoinNodNum(8)+1
      end if
    end do ! end of do r=2,9

    do r=1,8
      if (CoinNodNum(r).ge.2) then
        tempIPEdgNB(r)=ElemQ(1)
        tempIPEdgeNBNum=tempIPEdgeNBNum+1
      end if
    end do
    if (tempIPEdgeNBNum.eq.8) then
      IPEdgeNBFound=1.0d0
    end if

  end subroutine WhetherConnectInPlaneEdge


  !>  
  subroutine IPConnectedElement(IPNBFound,IfIPNbr,&
             IPNode1,IPNode2,q,ElemQ,NodeCoords,&
             tempIPNbrNum,IPNbrElemLoID,IPCoinElemIP)
    use Fatigue_Globals, only: maxNodeNum

    real*8, parameter :: zero=0.d0, one = 1.d0, GapToler=0.012

    integer IPNode1,IPNode2,ElemQ(9),tempIPNbrNum,q,&
     IPNbrElemLoID(11),IPCoinElemIP(11), r
    real*8 IfIPNbr,IPNBFound,NodeCoords(maxNodeNum,3),&
     temp1

    IfIPNbr=zero
    do r=2,9
      if ((IPNode1.eq.ElemQ(r)).or.&
       (IPNode2.eq.ElemQ(r))) then
        tempIPNbrNum=tempIPNbrNum+1;
        IPNbrElemLoID(tempIPNbrNum)=q
        IfIPNbr=one
        go to 7000
      end if
      ! by small radian (GapToler) sphere, considering that a gap
      ! may exist at an interlamianr interface location where
      ! top and bottom splitting cohesive intersect.
      ! by local IP associated node 1
      temp1=dsqrt((NodeCoords(IPNode1,1)-&
       NodeCoords(ElemQ(r),1))**2+&
       (NodeCoords(IPNode1,2)-&
       NodeCoords(ElemQ(r),2))**2+&
       (NodeCoords(IPNode1,3)-&
       NodeCoords(ElemQ(r),3))**2)
      if (temp1.lt.GapToler) then
        tempIPNbrNum=tempIPNbrNum+1;
        IPNbrElemLoID(tempIPNbrNum)=q
        IfIPNbr=one
        go to 7000
      end if
      ! by local IP associated node 2
      temp1=dsqrt((NodeCoords(IPNode2,1)-&
       NodeCoords(ElemQ(r),1))**2+&
       (NodeCoords(IPNode2,2)-&
       NodeCoords(ElemQ(r),2))**2+&
       (NodeCoords(IPNode2,3)-&
       NodeCoords(ElemQ(r),3))**2)
      if (temp1.lt.GapToler) then
        tempIPNbrNum=tempIPNbrNum+1;
        IPNbrElemLoID(tempIPNbrNum)=q
        IfIPNbr=one
        go to 7000
      end if
    end do
    ! check out which IP the current element's IP is next to
7000 if (IfIPNbr.eq.one) then
      if ((r.eq.2).or.(r.eq.6)) then
        IPCoinElemIP(tempIPNbrNum)=1
      elseif ((r.eq.3).or.(r.eq.7)) then
        IPCoinElemIP(tempIPNbrNum)=2
      elseif ((r.eq.4).or.(r.eq.8)) then
        IPCoinElemIP(tempIPNbrNum)=3
      else
        IPCoinElemIP(tempIPNbrNum)=4
      end if
    end if
    if (tempIPNbrNum.eq.11) then
      IPNBFound=one
    end if

  end subroutine IPConnectedElement


  !> Check whether an integration point is at edge or discontinuous location
  subroutine WhetherEdgeIP(p,IPNum,IPEdgNbrElem,&
                 IPNbrElemLoID,IPCoinElemIP,IfEdge)
    use Fatigue_Globals, only: maxCohElemNum

    real*8, parameter :: zero = 0.d0, one = 1.d0, two = 2.d0, three = 3.d0,&
                    &four = 4.d0

    integer p,IPEdgNbrElem(maxCohElemNum,9),NBNums(4),e,&
     IfEdge,IPNbrElemLoID(11),IPCoinElemIP(11)
    real(kind(zero)) IPNum



    ! find in-plane edge neighbour store numbers for current IP in IPEdgNbrElem
    if (IPNum.eq.one) then
      NBNums(1:4)=(/2,5,6,9/)
    elseif (IPNum.eq.two) then
      NBNums(1:4)=(/2,4,6,8/)
    elseif (IPNum.eq.three) then
      NBNums(1:4)=(/4,3,8,7/)
    elseif (IPNum.eq.four) then
      NBNums(1:4)=(/3,5,7,9/)
    end if

    ! Edge judgement for current IP by cheking whether its host element's
    ! corresponding in-plane edges are connected to elements or not.
    if((IPEdgNbrElem(p,NBNums(1)).eq.zero).or.&
     (IPEdgNbrElem(p,NBNums(2)).eq.zero).or.&
     (IPEdgNbrElem(p,NBNums(3)).eq.zero).or.&
     (IPEdgNbrElem(p,NBNums(4)).eq.zero))then
      IfEdge=one
      goto 8100
    end if
    ! check if the edge-sharing IPs' in-plane edges connect to elements
    do e=1,11
      if (IPNbrElemLoID(e).ne.zero) then
        if (IPCoinElemIP(e).eq.1) then
          if((IPEdgNbrElem(IPNbrElemLoID(e),2).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),5).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),6).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),9).eq.zero))then
            IfEdge=one
            goto 8100
          end if
        elseif (IPCoinElemIP(e).eq.2) then
          if((IPEdgNbrElem(IPNbrElemLoID(e),2).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),4).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),6).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),8).eq.zero))then
            IfEdge=one
            goto 8100
          end if
        elseif (IPCoinElemIP(e).eq.3) then
          if((IPEdgNbrElem(IPNbrElemLoID(e),4).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),3).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),8).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),7).eq.zero))then
            IfEdge=one
            goto 8100
          end if
        elseif (IPCoinElemIP(e).eq.4) then
          if((IPEdgNbrElem(IPNbrElemLoID(e),3).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),5).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),7).eq.zero).or.&
           (IPEdgNbrElem(IPNbrElemLoID(e),9).eq.zero))then
            IfEdge=one
            goto 8100
          end if
        end if
      end if
    end do
8100 continue

  end subroutine WhetherEdgeIP


  !> Compute the CROSS product of two vectors
  subroutine cross_ab(va,vb,axb)
    real*8 va(3),vb(3),axb(3)

    axb(1)=va(2)*vb(3)-va(3)*vb(2)
    axb(2)=va(3)*vb(1)-va(1)*vb(3)
    axb(3)=va(1)*vb(2)-va(2)*vb(1)

  end subroutine cross_ab


  !> Compute the DOT product of two vectors
  subroutine dot_ab(va,vb,ab)
    real*8 va(3),vb(3)
    real*8 ab

    ab=va(1)*vb(1)+va(2)*vb(2)+va(3)*vb(3)

  end subroutine dot_ab


  !> Convert string to upper case
  function to_upper(line) result(capitalChar)
    character(*), intent(in) :: line
    character(len=len(line)) :: capitalChar
    
    integer i

    do i = 1,len(line)
      if(line(i:i) .ge. "a" .and. line(i:i) .le. "z") then
           capitalChar(i:i) = char(ichar(line(i:i)) - 32)
      else
           capitalChar(i:i) = line(i:i)
      end if
    end do

  end function to_upper

end module Mesh_Utils
!DIR$ NOFREEFORM
!DIR$ FIXEDFORMLINESIZE:132