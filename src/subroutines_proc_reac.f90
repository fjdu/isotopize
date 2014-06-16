! Created on 2011-02-01 Tue 18:17:19

MODULE CMReac

use CMDT

implicit none

integer, parameter :: nItemInASpecies = 32, nMaxGrp = 32, &
  nDeutSideMax = 32, lenStrSideMax = 64, &
  idxSep = 20, eleCountMax = 2, &
  idxH=6, idxD=7

character(len=constLenNameSpecies) elementOld, elementNew

! idxOld = 6, idxNew = 7
integer idxOld, idxNew

integer nOtherDeutMax, noDMaxMetal
integer nDeutDegree
double precision noDEleAbundance

TYPE :: EleSpecies
  integer nItem
  integer, dimension(nItemInASpecies) :: vecIdEle, vecCtEle
END TYPE EleSpecies

TYPE :: EleGroup
  integer nItem, nEle, nWeight
  integer, dimension(:), allocatable :: EleList
  integer, dimension(:,:), allocatable :: ArrCt
END TYPE EleGroup

type (EleSpecies), dimension(:), allocatable :: SpeciesEleAll

CONTAINS

subroutine Name2Ele (SpeciesName, SpeciesEle)
  use CMDT
  implicit none
  character(len=*) SpeciesName
  type (EleSpecies) SpeciesEle
  integer i, i1, i2, nlen
  nlen = len(SpeciesName)
  SpeciesEle%nItem=0
  i=1
  do
    ! Start from position i, try to match the longest possible element
    ! name, return the index of the matched element in i1, and the end
    ! position of the matched element in i2.
    if (i .GT. nlen) exit
    i1 = MatchEle(SpeciesName, nlen, i, i2)
    if (i1 .GT. 0) then
      SpeciesEle%nItem = SpeciesEle%nItem + 1
      SpeciesEle%vecIdEle(SpeciesEle%nItem) = i1
      SpeciesEle%vecCtEle(SpeciesEle%nItem) = &
        getThisEleCount(SpeciesName, nlen, i2+1, i1)
      i = i1 + 1
    else
      i = i + 1
    end if
  end do
  SpeciesEle%vecIdEle(SpeciesEle%nItem+1 : nItemInASpecies) = 0
  SpeciesEle%vecCtEle(SpeciesEle%nItem+1 : nItemInASpecies) = 0
end subroutine Name2Ele


function MatchEle(SpeciesName, nlen, i1, i2)
  use CMDT
  implicit none
  character(len=*) SpeciesName
  integer i, i1, i2, nlen, MatchEle, lenMatch, lenMatchTmp
  MatchEle = 0
  lenMatch = 0
  do i=1, nElement
    if (MatchStr(SpeciesName(i1:nlen), trim(nameElements(i)), &
      lenMatchTmp)) then
      if (lenMatchTmp .GT. lenMatch) then
        MatchEle = i
        lenMatch = lenMatchTmp
      end if
    end if
  end do
  if (MatchEle .GT. 0) then
    i2 = i1 + lenMatch - 1
  else
    i2 = i1
  end if
end function MatchEle


function getThisEleCount(SpeciesName, nlen, i1, i2)
  use CMDT
  implicit none
  character(len=*) SpeciesName
  integer i, i1, i2, nlen, getThisEleCount
  logical flagFound
  flagFound = .FALSE.
  do i=i1, nlen
    if (IsDigit(SpeciesName(i:i))) then
      i2 = i
      flagFound = .TRUE.
    else
      exit
    end if
  end do
  if (flagFound) then
    read (SpeciesName(i1:i2), '(I16)') getThisEleCount
  else
    getThisEleCount = 1
    i2 = i1 - 1
  end if
end function getThisEleCount


function MatchStr(Str1, Str2, len2)
! Determine if Str2 is a substring of Str1.
! The match point must be the leading position.
  implicit none
  character(len=*) Str1, Str2
  integer len1, len2, i
  logical MatchStr
  len1 = len_trim(Str1)
  len2 = len_trim(Str2)
  MatchStr = .TRUE.
  if (len1 .GE. len2) then
    do i=1, len2
      if (Str1(i:i) .NE. Str2(i:i)) then
        MatchStr = .FALSE.
        exit
      end if
    end do
  else
    MatchStr = .FALSE.
  end if
end function MatchStr


function IsDigit(ch)
  logical IsDigit
  character ch
  if (LGE(ch, '0') .AND. LLE(ch, '9')) then
    IsDigit = .TRUE.
  else
    IsDigit = .FALSE.
  end if
end function IsDigit


subroutine Ele2Elements(SpeciesEle, SpeElements)
  use CMDT
  implicit none
  type (EleSpecies) SpeciesEle
  integer, dimension(nElement) :: SpeElements
  integer i
  SpeElements = 0
  do i=1, SpeciesEle%nItem
    SpeElements(SpeciesEle%vecIdEle(i)) = &
      SpeElements(SpeciesEle%vecIdEle(i)) + &
      SpeciesEle%vecCtEle(i)
  end do
end subroutine Ele2Elements


function IsEquiv(SpeciesEle1, SpeciesEle2)
  use CMDT
  implicit none
  integer i, nLenTmp
  logical IsEquiv
  type (EleSpecies) SpeciesEle1, SpeciesEle2
  IsEquiv = .TRUE.
  if (SpeciesEle1%nItem .EQ. SpeciesEle2%nItem) then
    if (SpeciesEle1%nItem .EQ. 1) then
      if ((SpeciesEle1%vecIdEle(1) .NE. SpeciesEle2%vecIdEle(1)) .OR. &
          (SpeciesEle1%vecCtEle(1) .NE. SpeciesEle2%vecCtEle(1))) then
        IsEquiv = .FALSE.
      end if
    else
      ! H3COH will be considered to be equivalent to CH3OH.
      if (((SpeciesEle1%vecIdEle(1) .NE. SpeciesEle2%vecIdEle(1)) .OR. &
           (SpeciesEle1%vecCtEle(1) .NE. SpeciesEle2%vecCtEle(1)) .OR. &
           (SpeciesEle1%vecIdEle(2) .NE. SpeciesEle2%vecIdEle(2)) .OR. &
           (SpeciesEle1%vecCtEle(2) .NE. SpeciesEle2%vecCtEle(2))) &
           .AND. &
          ((SpeciesEle1%vecIdEle(1) .NE. SpeciesEle2%vecIdEle(2)) .OR. &
           (SpeciesEle1%vecCtEle(1) .NE. SpeciesEle2%vecCtEle(2)) .OR. &
           (SpeciesEle1%vecIdEle(2) .NE. SpeciesEle2%vecIdEle(1)) .OR. &
           (SpeciesEle1%vecCtEle(2) .NE. SpeciesEle2%vecCtEle(1)))) then
        IsEquiv = .FALSE.
      else
        do i=3, SpeciesEle1%nItem
          if ((SpeciesEle1%vecIdEle(i) .NE. SpeciesEle2%vecIdEle(i)) .OR. &
              (SpeciesEle1%vecCtEle(i) .NE. SpeciesEle2%vecCtEle(i))) then
            IsEquiv = .FALSE.
            exit
          end if
        end do
      end if
      ! With the following part, HCN and NCH will be considered equivalent.
      ! Use this part with caution!
      if (.NOT. IsEquiv) then
        IsEquiv = .TRUE.
        if ((SpeciesEle1%vecIdEle(SpeciesEle1%nItem) .EQ. 1) .OR. &
            (SpeciesEle1%vecIdEle(SpeciesEle1%nItem) .EQ. 2)) then
          nLenTmp = SpeciesEle1%nItem - 1
        else
          nLenTmp = SpeciesEle1%nItem
        end if
        ! Assume the charge symbol is always put at the last position.
        if ((SpeciesEle1%vecIdEle(SpeciesEle1%nItem) .NE. SpeciesEle2%vecIdEle(SpeciesEle1%nItem)) .OR. &
            (SpeciesEle1%vecCtEle(SpeciesEle1%nItem) .NE. SpeciesEle2%vecCtEle(SpeciesEle1%nItem))) then
          IsEquiv = .FALSE.
        else
          do i=1, nLenTmp
            if ((SpeciesEle1%vecIdEle(i) .NE. SpeciesEle2%vecIdEle(nLenTmp+1-i)) .OR. &
                (SpeciesEle1%vecCtEle(i) .NE. SpeciesEle2%vecCtEle(nLenTmp+1-i))) then
              IsEquiv = .FALSE.
              exit
            end if
          end do
        end if
      end if
      !!!!!!!!!
    end if
  else
    IsEquiv = .FALSE.
  end if
end function IsEquiv


subroutine DeutGroups (HydrGroup, nDeu, vecDeutGroup, nDeutGroup)
  use CMDT
  implicit none
  integer i, j, k
  integer nDeu, nDeutGroup, nTotalH, nGamb, binRepre, binRepreBak
  integer, dimension(:), allocatable :: ArrCt
  type (EleGroup), dimension(*) :: vecDeutGroup
  type (EleGroup) HydrGroup
  logical flagMatch
  nDeutGroup = 0
  nTotalH = sum(HydrGroup%ArrCt(1:HydrGroup%nItem, 1))
  if (nTotalH .LT. nDeu) then
    !write (*,*) 'nTotalH .LT. nDeu'
    return
  end if
  if (nTotalH .GT. 32) then
    write (*,*) 'nTotalH .GT. 32'
    return
  end if
  nGamb = binomcoeff(nTotalH, nDeu)
  allocate(ArrCt(HydrGroup%nItem))
  binRepre = ISHFT(1, nDeu) - 1
  do i=1, nGamb
    binRepreBak = binRepre
    ArrCt = 0
    do j=1, HydrGroup%nItem
      do k=1, HydrGroup%ArrCt(j,1)
        ArrCt(j) = ArrCt(j) + IAND(binRepre, 1)
        binRepre = ISHFT(binRepre, -1)
      end do
    end do
    binRepre = getNextColex(binRepreBak)
    flagMatch = .FALSE.
    do j=1, nDeutGroup
      flagMatch = .TRUE.
      do k=1, HydrGroup%nItem
        if (vecDeutGroup(j)%ArrCt(k,2) .NE. ArrCt(k)) then
          flagMatch = .FALSE.
          exit
        end if
      end do
      if (flagMatch) then
        vecDeutGroup(j)%nWeight = vecDeutGroup(j)%nWeight + 1
        exit
      end if
    end do
    if (.NOT. flagMatch) then
      nDeutGroup = nDeutGroup + 1
      vecDeutGroup(nDeutGroup)%nItem = HydrGroup%nItem
      vecDeutGroup(nDeutGroup)%nEle = 2
      vecDeutGroup(nDeutGroup)%nWeight = 1
      vecDeutGroup(nDeutGroup)%EleList(1) = idxOld
      vecDeutGroup(nDeutGroup)%EleList(2) = idxNew
      do j=1, HydrGroup%nItem
        vecDeutGroup(nDeutGroup)%ArrCt(j, 2) = ArrCt(j)
        vecDeutGroup(nDeutGroup)%ArrCt(j, 1) = &
          HydrGroup%ArrCt(j, 1) - ArrCt(j)
      end do
    end if
  end do
  deallocate(ArrCt)
end subroutine DeutGroups


subroutine DeutIns (DeutGroup, SpeciesEle, SpeciesEleD)
  implicit none
  type (EleSpecies) SpeciesEle, SpeciesEleD
  type (EleGroup) DeutGroup
  integer i, i1, i2
  i1 = 0
  i2 = 0
  do i=1, SpeciesEle%nItem
    if (SpeciesEle%vecIdEle(i) .NE. idxOld) then
      i1 = i1 + 1
      SpeciesEleD%vecIdEle(i1) = SpeciesEle%vecIdEle(i)
      SpeciesEleD%vecCtEle(i1) = SpeciesEle%vecCtEle(i)
    else
      i2 = i2 + 1
      i1 = i1 + 1
      SpeciesEleD%vecIdEle(i1) = idxOld
      SpeciesEleD%vecCtEle(i1) = DeutGroup%ArrCt(i2,1)
      i1 = i1 + 1
      SpeciesEleD%vecIdEle(i1) = idxNew
      SpeciesEleD%vecCtEle(i1) = DeutGroup%ArrCt(i2,2)
    end if
  end do
  SpeciesEleD%nItem = i1
end subroutine DeutIns



function IsEquivSide(Side1, Side2, n1, n2)
  use CMDT
  implicit none
  logical IsEquivSide
  integer i, j, n1, n2
  character(len=constLenNameSpecies), dimension(:) :: &
    Side1, Side2
  character(len=constLenNameSpecies), dimension(:), allocatable :: &
    Side1_sorted, Side2_sorted
  type (EleSpecies) Sp1, Sp2
  !
  IsEquivSide = .TRUE.
  if (n1 .NE. n2) then
    IsEquivSide = .FALSE.
  else
    allocate(Side1_sorted(n1), Side2_sorted(n1))
    call SortStrList(Side1, Side1_sorted, n1)
    call SortStrList(Side2, Side2_sorted, n1)
    do i=1, n1
      if (Side1_sorted(i) .NE. Side2_sorted(i)) then
        IsEquivSide = .FALSE.
        exit
      end if
    end do
    if (.not. IsEquivSide) then
      do i=1, n1
        call Name2Ele(Side1(i), Sp1)
        call Name2Ele(Side2(i), Sp2)
        IsEquivSide = isEquivMirror(Sp1, Sp2)
        if (.not. IsEquivSide) then
          return
        end if
      end do
    end if
  end if
end function IsEquivSide



subroutine SortStrList(Str, Str_sorted, Len)
  use CMDT
  implicit none
  integer i, j, Len
  character(len=constLenNameSpecies), dimension(Len) :: Str, Str_sorted
  character(len=constLenNameSpecies) str_tmp
  Str_sorted = Str
  do i=1, Len
    do j=1, i-1
      if (LGT(Str_sorted(j), Str_sorted(i))) then
        str_tmp = Str_sorted(j)
        Str_sorted(j) = Str_sorted(i)
        Str_sorted(i) = str_tmp
      end if
    end do
  end do
end subroutine SortStrList


subroutine DeutReac (iReac, nDeut, fU)
  use CMDT
  implicit none
  character(len=9) strtmp
  integer i, j, i1, j1, k, iReac, nDeut, &
    nDeutThis, nDeutedLeft, nDeutedRight, fU
  integer TotalWeight, nSplittedLeft, nSplittedRight, nSplitted_tmp
  character(len=lenStrSideMax), dimension(nDeutSideMax) :: &
    StrDeutedLeft, StrDeutedRight
  integer, dimension(nDeutSideMax) :: WeightsLeft, WeightsRight
  character(len=constLenNameSpecies), dimension(nReactants) :: &
    StrSplittedLeft
  character(len=constLenNameSpecies), dimension(nProducts) :: &
    StrSplittedRight
  character(len=constLenNameSpecies), dimension(nProducts) :: &
    StrSplitted_tmp
  logical flag_same

  do nDeutThis=1, nDeut
    call DeutSide(Reactions(1:nRealReactants(iReac), iReac), &
           nRealReactants(iReac), nDeutThis, &
           StrDeutedLeft, WeightsLeft, nDeutedLeft)
    call DeutSide( &
           Reactions(nReactants+1:nReactants+nRealProducts(iReac), iReac), &
           nRealProducts(iReac), nDeutThis, &
           StrDeutedRight, WeightsRight, nDeutedRight)
    TotalWeight = sum(WeightsRight(1:nDeutedRight))
    do i=1, nDeutedLeft
      StrSplittedLeft=""
      call SplitSideStr(StrDeutedLeft(i), StrSplittedLeft, nSplittedLeft)
      flag_same = .FALSE.
      do i1=1, i-1
        StrSplitted_tmp = ""
        call SplitSideStr(StrDeutedLeft(i1), StrSplitted_tmp, nSplitted_tmp)
        if (IsEquivSide(StrSplittedLeft, StrSplitted_tmp, &
                        nSplittedLeft, nSplitted_tmp)) then
          flag_same = .TRUE.
          exit
        end if
      end do
      if (flag_same) cycle
      do j=1, nDeutedRight
        !
        if (WeightsRight(j) .LE. 0) cycle
        !
        StrSplittedRight=""
        call SplitSideStr(StrDeutedRight(j), StrSplittedRight, nSplittedRight)
        do j1=j+1, nDeutedRight
          StrSplitted_tmp = ""
          call SplitSideStr(StrDeutedRight(j1), StrSplitted_tmp, nSplitted_tmp)
          if (IsEquivSide(StrSplittedRight, StrSplitted_tmp, &
              nSplittedRight, nSplitted_tmp)) then
            WeightsRight(j) = WeightsRight(j) + WeightsRight(j1)
            WeightsRight(j1) = 0
          end if
        end do
        do k=nRealReactants(iReac)+1, nReactants
          if (len_trim(strReactants(k, iReac)) .gt. 0) then
            StrSplittedLeft(k) = strReactants(k, iReac)
          end if
        end do
        do k=nRealProducts(iReac)+1, nProducts
          if (len_trim(strProducts(k, iReac)) .gt. 0) then
            StrSplittedRight(k) = strProducts(k, iReac)
          end if
        end do
        call double2str(strtmp, dblABC(3, iReac), 9, 1)
        write (fU, &
          '(7A12, ES9.2, F9.2, A9, 2I6, I3, X,A1,X,A2,X,A1, " !", 4I3)') &
          StrSplittedLeft, StrSplittedRight, &
          dblABC(1, iReac) * dble(WeightsRight(j))/dble(TotalWeight), &
          dblABC(2, iReac), &
          strtmp, &
          int(T_min(iReac)), int(T_max(iReac)), typeReac(iReac), &
          cquality(iReac), ctype(iReac), stype(iReac), &
          WeightsRight(j), j, i, nDeutThis
      end do
    end do
  end do
end subroutine DeutReac


subroutine DeutSide (SpeciesGroup, &
           nGroupSize, nDeutThis, &
           StrDeutedSide, WeightsSide, nDeutedSide)
  implicit none
  type (EleSpecies) SpeciesEleD, SpeciesEleBig
  type (EleGroup) HydrGroup
  type (EleGroup), dimension(nMaxGrp) :: vecDeutGroup
  integer i, nDeutedSide, nGroupSize, nDeutThis
  integer, dimension(:) :: WeightsSide, SpeciesGroup
  character(len=*), dimension(:) :: StrDeutedSide
  call makeBigSpecies (SpeciesGroup, nGroupSize, SpeciesEleBig)
  HydrGroup%nItem = getHAprCt(SpeciesEleBig)
  HydrGroup%nEle = 1
  if (allocated(HydrGroup%EleList)) deallocate(HydrGroup%EleList)
  if (allocated(HydrGroup%ArrCt)) deallocate(HydrGroup%ArrCt)
  allocate(HydrGroup%EleList(1))
  allocate(HydrGroup%ArrCt(HydrGroup%nItem,1))
  do i=1, nMaxGrp
    if (.NOT. allocated(vecDeutGroup(i)%EleList)) &
      allocate(vecDeutGroup(i)%EleList(2))
    if (allocated(vecDeutGroup(i)%ArrCt)) &
      deallocate(vecDeutGroup(i)%ArrCt)
    allocate(vecDeutGroup(i)%ArrCt(HydrGroup%nItem, 2))
  end do
  HydrGroup%EleList(1) = idxOld
  call getHGroup (SpeciesEleBig, HydrGroup)
  call DeutGroups (HydrGroup, nDeutThis, vecDeutGroup, nDeutedSide)
  do i=1, nDeutedSide
    call DeutIns(vecDeutGroup(i), SpeciesEleBig, SpeciesEleD)
    StrDeutedSide(i) = ele2str(SpeciesEleD)
    WeightsSide(i) = vecDeutGroup(i)%nWeight
  end do
end subroutine DeutSide


subroutine makeBigSpecies (SpeciesGroup, nGroupSize, SpeciesEleBig)
  implicit none
  integer, dimension(:) :: SpeciesGroup
  integer nGroupSize, i
  type (EleSpecies) :: SpeciesEleBig
  SpeciesEleBig = SpeciesEleAll(SpeciesGroup(1))
  !Combine the EleSpecies structures of all the species in a group together.
  do i=2, nGroupSize
    SpeciesEleBig%nItem = SpeciesEleBig%nItem + 1
    SpeciesEleBig%vecIdEle(SpeciesEleBig%nItem) = idxSep
    SpeciesEleBig%vecCtEle(SpeciesEleBig%nItem) = 1
    SpeciesEleBig%vecIdEle(SpeciesEleBig%nItem+1 : &
      SpeciesEleBig%nItem+SpeciesEleAll(SpeciesGroup(i))%nItem) = &
      SpeciesEleAll(SpeciesGroup(i))%vecIdEle(1:&
        SpeciesEleAll(SpeciesGroup(i))%nItem)
    SpeciesEleBig%vecCtEle(SpeciesEleBig%nItem+1 : &
      SpeciesEleBig%nItem+SpeciesEleAll(SpeciesGroup(i))%nItem) = &
      SpeciesEleAll(SpeciesGroup(i))%vecCtEle(1:&
        SpeciesEleAll(SpeciesGroup(i))%nItem)
    SpeciesEleBig%nItem = SpeciesEleBig%nItem + &
      SpeciesEleAll(SpeciesGroup(i))%nItem
  end do
end subroutine makeBigSpecies


subroutine SplitSideStr (SideStr, StrSplitted, nSpeSplitted)
  use CMDT
  implicit none
  character(len=lenStrSideMax) SideStr
  character(len=*), dimension(:) :: StrSplitted
  integer nSpeSplitted, i, j
  nSpeSplitted = 1
  j = 1
  do i = 1, len_trim(SideStr)
    if (SideStr(i:i) .NE. nameElements(idxSep)) then
      StrSplitted(nSpeSplitted)(j:j) = SideStr(i:i)
      j = j + 1
    else
      nSpeSplitted = nSpeSplitted + 1
      j = 1
    end if
  end do
end subroutine SplitSideStr


function ele2str (EleSpe)
  use CMDT
  implicit none
  character(len=lenStrSideMax) ele2str
  character(len=eleCountMax) ctmp
  TYPE (EleSpecies) EleSpe
  integer i
  ele2str=""
  do i=1, EleSpe%nItem
    if (EleSpe%vecCtEle(i) .GT. 1) then
      write(ctmp, '(I2)') , EleSpe%vecCtEle(i)
      ele2str = trim(ele2str)//trim(nameElements(EleSpe%vecIdEle(i)))//&
        trim(ADJUSTL(ctmp))
    end if
    if (EleSpe%vecCtEle(i) .EQ. 1) then
      ele2str = trim(ele2str)//trim(nameElements(EleSpe%vecIdEle(i)))
    end if
  end do
end function ele2str


function getHAprCt (EleSpeciesA)
  implicit none
  integer i, getHAprCt
  TYPE(EleSpecies) EleSpeciesA
  getHAprCt = 0
  do i=1, EleSpeciesA%nItem
    if (EleSpeciesA%vecIdEle(i) .EQ. idxOld) then
      getHAprCt = getHAprCt + 1
    end if
  end do
end function getHAprCt


subroutine getHGroup (EleSpeciesA, HydrGroup)
  use CMDT
  implicit none
  integer i, i1
  TYPE(EleSpecies) EleSpeciesA
  type (EleGroup) HydrGroup
  i1 = 0
  do i=1, EleSpeciesA%nItem
    if (EleSpeciesA%vecIdEle(i) .EQ. idxOld) then
      i1 = i1 + 1
      HydrGroup%ArrCt(i1,1)  =  EleSpeciesA%vecCtEle(i)
    end if
  end do
end subroutine getHGroup


function binomcoeff (n, m)
  implicit none
  integer binomcoeff, n, m, i
  binomcoeff = 1
  do i=max(m,n-m)+1, n
    binomcoeff = binomcoeff * i
  end do
  do i=2, min(m,n-m)
    binomcoeff = binomcoeff / i
  end do
end function binomcoeff


function getNextColex (x)
! Gosper's hack.
  implicit none
  integer getNextColex, x, s, r
  s = IAND(x, -x)
  r = s + x
  getNextColex = IOR(r, ISHFT(IEOR(x,r), -2)/s)
end function getNextColex


subroutine intExchange(x, y)
! Or a even easier one: x=x+y; y=x-y; x=x-y
  implicit none
  integer x, y
  x = ieor(x, y)
  y = ieor(x, y)
  x = ieor(x, y)
end subroutine intExchange


function getAltAB(x, A, B)
  implicit none
  integer getAltAB, x, A, B
  getAltAB = A - x + B
  !getAltAB = IEOR(A, IEOR(B, x))
end function getAltAB


function isEquivMirror(Sp1, Sp2) result(eq)
  logical eq
  type(EleSpecies), intent(in) :: Sp1, Sp2
  type(EleSpecies) St
  integer i, j
  !
  eq = isEquivSimple(Sp1, Sp2)
  if (eq) then
    return
  end if
  !
  do j=1, 3
    St = mirror_eleSpecies(Sp1, j)
    eq = isEquivSimple(St, Sp2)
    if (eq) then
      return
    end if
  end do
end function isEquivMirror



function isEquivSimple(Sp1, Sp2) result(eq)
  logical eq
  type(EleSpecies), intent(in) :: Sp1, Sp2
  integer i
  eq = .true.
  if (Sp1%nItem .ne. Sp2%nItem) then
    eq = .false.
    return
  end if
  !
  do i=1, Sp1%nItem
    if ((Sp1%vecIdEle(i) .ne. Sp2%vecIdEle(i)) .or. &
        (Sp1%vecCtEle(i) .ne. Sp2%vecCtEle(i))) then
      eq = .false.
      return
    end if
  end do
end function isEquivSimple


function mirror_eleSpecies(eS, cases) result(m)
  type(EleSpecies) m
  type(EleSpecies), intent(in) :: eS
  integer, intent(in) :: cases
  integer i, n1, n2
  !
  select case(cases)
    case(1)
      ! Transform CH3OCH3 into H3COH3C
      m = simple_mirror(eS)
    case(2)
      ! Transform CH3OCH2D into CH2DOCH3
      !   First transform into DH2COH3C
      !
      m = simple_mirror(eS)
      !
      if (m%nItem .eq. 1) then
        return
      end if
      !
      call identify_head_H_group(m, n1, n2)
      if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
        call sort_sub_group(m, n1, n2)
      end if
      !
      call identify_tail_H_group(m, n1, n2)
      if ((n1 .gt. 0) .or. (n2 .gt. 0)) then
        call sort_sub_group(m, n2, n1)
      end if
      !
    case(3)
      ! Transform H3COCH2D into CH3OCH2D
      !
      m = eS
      !
      if (m%nItem .eq. 1) then
        return
      end if
      !
      call identify_head_H_group(m, n1, n2)
      if ((n1 .gt. 0) .and. (n2 .gt. 0)) then
        call sort_sub_group(m, n1, n2)
      end if
      !
      call identify_tail_H_group(m, n1, n2)
      if ((n1 .gt. 0) .or. (n2 .gt. 0)) then
        call sort_sub_group(m, n2, n1)
      end if
      !
  end select
end function mirror_eleSpecies


function simple_mirror(eS) result(m)
  type(EleSpecies) m
  type(EleSpecies), intent(in) :: eS
  integer i, i1, i2
  if (eS%vecIdEle(1) .lt. idxH) then
    i1 = 2
  else
    i1 = 1
  end if
  if (eS%vecIdEle(eS%nItem) .lt. idxH) then
    i2 = eS%nItem - 1
  else
    i2 = eS%nItem
  end if
  !
  m%nItem = eS%nItem
  !
  do i=1, i1-1
    m%vecIdEle(i) = eS%vecIdEle(i)
    m%vecCtEle(i) = eS%vecCtEle(i)
  end do
  do i=i2+1, eS%nItem
    m%vecIdEle(i) = eS%vecIdEle(i)
    m%vecCtEle(i) = eS%vecCtEle(i)
  end do
  do i=i1, i2
    m%vecIdEle(i) = eS%vecIdEle(i2-i+1)
    m%vecCtEle(i) = eS%vecCtEle(i2-i+1)
  end do
end function simple_mirror


subroutine identify_head_H_group(eS, n1, n2)
  type(EleSpecies), intent(in) :: eS
  integer, intent(out) :: n1, n2
  integer i
  n1 = 0; n2 = 0
  do i=1, eS%nItem
    if (eS%vecIdEle(i) .ge. idxH) then
      n1 = i
      exit
    end if
  end do
  do i=n1+1, eS%nItem
    if ((eS%vecIdEle(i) .ne. idxH) .and. (es%vecIdEle(i) .ne. idxD)) then
      n2 = i
      exit
    end if
  end do
  if ((es%vecIdEle(n1) .ne. idxH) .and. (es%vecIdEle(n1) .ne. idxD)) then
    n2 = n2 - 1
  end if
end subroutine identify_head_H_group


subroutine identify_tail_H_group(eS, n1, n2)
  type(EleSpecies), intent(in) :: eS
  integer, intent(out) :: n1, n2
  integer i
  n1 = 0; n2 = 0
  do i=eS%nItem, 1, -1
    if (eS%vecIdEle(i) .ge. idxH) then
      n1 = i
      exit
    end if
  end do
  do i=n1-1, 1, -1
    if ((eS%vecIdEle(i) .ne. idxH) .and. (es%vecIdEle(i) .ne. idxD)) then
      n2 = i
      exit
    end if
  end do
  if ((es%vecIdEle(n1) .ne. idxH) .and. (es%vecIdEle(n1) .ne. idxD)) then
    n2 = n2 + 1
  end if
end subroutine identify_tail_H_group


subroutine sort_sub_group(eS, i1, i2)
  type(EleSpecies), intent(inout) :: eS
  integer, intent(in) :: i1, i2
  integer i, j, ia, ib
  do i=i1, i2
    if (es%vecIdEle(i) .eq. idxD) then
      cycle
    end if
    do j=i1, i-1
      if ((eS%vecIdEle(j) .eq. idxH) .or. (eS%vecIdEle(j) .eq. idxD)) then
        ia = eS%vecIdEle(j)
        ib = eS%vecCtEle(j)
        eS%vecIdEle(j) = es%vecIdEle(i)
        eS%vecCtEle(j) = es%vecCtEle(i)
        eS%vecIdEle(i) = ia
        eS%vecCtEle(i) = ib
      end if
    end do
  end do
end subroutine sort_sub_group


pure subroutine double2str(str, x, nW, nP)
  double precision, intent(in) :: x
  integer, intent(in) :: nW, nP
  character(len=*), intent(out) :: str
  character(len=16) :: fmtstr
  write(fmtstr, '("(F", I2, ".", I2, ")")') nW, nP
  write(str, fmtstr) x
  if (str(1:1) .eq. '*') then
    write(fmtstr, '("(ES", I2, ".", I2, ")")') nW, nP
    write(str, fmtstr) x
  end if
end subroutine double2str

END MODULE
