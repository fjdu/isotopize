      program makedeut
      use CMDT
      use CMReac
      implicit none

      external f, jac, getPhyPar

      integer fU, fUyHistory_ascii, fUyHistory_bin, fUPhyPar
      !
      !
      integer ios, statALLOC
      integer nLineAll, nLineData
      character strTMP*128, strTMP1*128
      character FMTstr*128, FMTstryHistory*128, FMTstrPhyPar*128
      character(LEN=constLenNameSpecies) nameSpecies_tmp
      character(LEN=9) strtmp2
      logical IsWordChar

      integer i, j, k, h, i1, i2, i3, nNumMMTmp
      double precision, dimension(3) :: dblTmpVec
      double precision initialAbundance_tmp
      integer, dimension(8) :: intTmpvec
      logical flag, flag1, flag2, flagChange

      real ProgStartTime, ProgCurrentTime

      integer, dimension(:), allocatable :: &
        nElementReac, nElementProd

      integer nIterationBase, nTargetDecades, nerr, nChange
      double precision ratioTout, tStep
      double precision t, tout, tFinal
      double precision, dimension(:), allocatable :: y, yBak, RWORK
      double precision, dimension(:), allocatable :: ydotTMP
      integer IOPT, ISTATE, ITASK, ITOL, LIW, LRW, NEQ, MF, NNZ
      integer, dimension(:), allocatable :: IWORK

      logical getPhyPar, ExternalAction, getFileUnit, FileUnitOpened
      integer, parameter :: nPhyPar = 8

      integer getIdxSpecies, idxSpeciesAna
      
      character(len=128) getFilePreName
      character(len=256) combineStrArr

      type (EleSpecies), dimension(:), allocatable :: SpeciesEleD
      type (EleGroup) HydrGroup
      type (EleGroup), dimension(nMaxGrp) :: vecDeutGroup
      integer nDeutGroup
      !
      character(len=128) fDeuterated
      !
      character(len=160) :: fmtHeaderHerbst = &
        '#23456123456789ABCD123456789ABCD123456789A&
        &BCD123456789ABCD123456789ABCD123456789ABCD&
        E123456789A123456789A1234567891234'
      character(len=160) :: fmtHeaderMine = &
        '!23456789ABC123456789ABC123456789ABC123456789ABC&
        &123456789ABC123456789ABC123456789ABC123456789123456789&
        &12345678912345612345612312123      w  j  i  t'

      CALL CPU_TIME(ProgStartTime)
      call date_and_time(DATE=strTMP, TIME=strTMP1)

      write (*,'(A, A60, A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program start', '', char(27)//'[0m'
      write (*,'("Current time: ", A8, 2X, A10)') strTMP, strTMP1

      CALL GET_COMMAND_ARGUMENT(1, strTMP, i, ios)
      if (i .EQ. 0) strTMP = 'config.dat'

      write (*, '(/A)') 'Initializing...'
      CALL initialize (strTMP)

      idxOld = -1
      idxNew = -1
      do i=1, nElement
        if (nameElements(i) .eq. adjustl(elementOld)) then
          idxOld = i
        end if
        if (nameElements(i) .eq. adjustl(elementNew)) then
          idxNew = i
        end if
      end do
      if (idxOld .eq. -1) then
        write(*, '(A, A12, A)') 'Error: ', elementOld, ' is unknown!'
        stop
      end if
      if (idxNew .eq. -1) then
        idxNew = idxSep + 1
        if (idxNew .gt. nElement) then
          write(*, '(A)') 'The length of the elements list is too small!'
          stop
        end if
        nameElements(idxNew) = adjustl(elementNew)
      end if
      !
      nameElements(4) = inputGrainEleName
      nameElements(5) = inputGrainName
      !
      if (inputFormat .eq. 'Herbst') then
        nReactants = 2
      else
        nReactants = 3
      end if
      nProducts = 4

! Import all the reactions

      write (*, '(/A)') 'Importing reactions...'

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      CALL openFileSequentialRead &
        (fU, trim(path)//trim(fReactions), 999)
      CALL GetNLineF_alt (fU, nLineAll, nLineData, commentChar)
      write(*,*) nLineData, nLineAll

      nReactions = nLineData

      if (nSpecies_Est .LE. 0) & ! This estimation would be surely
        nSpecies_Est = nReactions * (nReactants+nProducts) ! enough.

      allocate &
        (strReactants(nReactants, nReactions), &
         strProducts(nProducts, nReactions), &
         nameSpecies(nSpecies_Est), &
         nRealReactants(nReactions), &
         nRealProducts(nReactions), &
         reactions(nReactants+nProducts, nReactions), &
         dblABC(3, nReactions), &
         T_min(nReactions), &
         T_max(nReactions), &
         typeReac(nReactions), &
         ctype(nReactions), &
         stype(nReactions), &
         cquality(nReactions), &
         rates(nReactions), STAT=statALLOC)

      CALL ReadReactions (fU, nLineAll, commentChar, inputFormat)
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
      !
      CALL MakeReactions

      allocate &
        (SpeciesElements(nElement, nSpecies), &
         massSpecies(nSpecies), &
         nElementReac(nElement), &
         nElementProd(nElement), &
         nReac_AsReactants(nSpecies), &
         nReac_AsProducts(nSpecies), &
         idxReac_AsReactants(nReactions, nSpecies), &
         idxReac_AsProducts(nReactions, nSpecies), &
         idxGrReactions(nReactions), &
         SpeciesEleAll(nSpecies), &
         SpeciesEleD(nSpecies), STAT=statALLOC)

      if (statALLOC .NE. 0) then
        write (*,'(/A32/)') "Error in allocating!"
        stop
      end if
      !
      nReacCounter = 1
      !
      do i=1, nSpecies
        !
        call Name2Ele (nameSpecies(i), SpeciesEleAll(i))
        call Ele2Elements(SpeciesEleAll(i), SpeciesElements(:, i))
        !
        massSpecies(i) = sum(SpeciesElements(:, i) * ElementMassNumber)
        !
        do j=1, i-1
          if (IsEquiv(SpeciesEleAll(i), SpeciesEleAll(j))) then
            write (*,*) nameSpecies(i), ' and ', nameSpecies(j), &
              'might be equivalent.'
            !write (*, '(48I3)') SpeciesEleAll(i)%vecIdEle, SpeciesEleAll(i)%vecCtEle
            !write (*, '(48I3)') SpeciesEleAll(j)%vecIdEle, SpeciesEleAll(j)%vecCtEle
          end if
        end do
      end do

      write (*,*) 'Working on the reactions...'
      !
      nameElements(4) = outputGrainEleName
      !
      fDeuterated = &
        trim(getFilePreName(fReactions))//'_isotopized.dat'
      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      !
      CALL openFileSequentialWrite &
        (fU, trim(path)//trim(fDeuterated), 9999)
      if (outputFormat .eq. 'Herbst') then
        write (fU, '(A)') trim(fmtHeaderHerbst)
      else
        write (fU, '(A)') trim(fmtHeaderMine)
      end if
      nReac_AsReactants = 0
      nReac_AsProducts = 0
      do i=1, nReactions
        !
        do j=1, nRealReactants(i)
          flag = .TRUE.
          do k=1, j-1
            if (reactions(k, i) .EQ. reactions(j, i)) then
              flag = .FALSE.
              exit
            end if
          end do
          if (flag) then
            nReac_AsReactants(reactions(j, i)) = &
              nReac_AsReactants(reactions(j, i)) + 1
            idxReac_AsReactants(nReac_AsReactants(reactions(j, i)), &
              reactions(j, i)) = i
          end if
        end do
        !
        do j=1+nReactants, nRealProducts(i)+nReactants
          flag = .TRUE.
          do k=1+nReactants, j-1
            if (reactions(k, i) .EQ. reactions(j, i)) then
              flag = .FALSE.
              exit
            end if
          end do
          if (flag) then
            nReac_AsProducts(reactions(j, i)) = &
              nReac_AsProducts(reactions(j, i)) + 1
            idxReac_AsProducts(nReac_AsProducts(reactions(j, i)), &
              reactions(j, i)) = i
          end if
        end do
        !
        nElementReac = 0
        nElementProd = 0
        do j=1, nRealReactants(i)
          nElementReac = nElementReac + &
            SpeciesElements(:, reactions(j, i))
        end do
        do j=1, nRealProducts(i)
          nElementProd = nElementProd + &
            SpeciesElements(:, reactions(nReactants+j, i))
        end do
        if ((abs(nElementReac(1) - nElementReac(2) &
               - nElementProd(1) + nElementProd(2)) + &
          sum(abs(nElementReac(6:nElement) - &
                  nElementProd(6:nElement)))) .NE. 0) then
          ! Start from H, hence Grain is not taken into account for the elemental conservation.
          write (*, '(2A, I6, A, 2X, 2A12, " -> ", 5A12)') &
            'Elements not conserved: ', &
            char(27)//'[41m', i, char(27)//'[0m', & !]]
            strReactants(1:nReactants, i), strProducts(1:nProducts, i)
          write(*, '(24A3)') nameElements(1:24)
          write(*, '(24I3)') nElementReac(1:24)
          write(*, '(24I3)') nElementProd(1:24)
          write(*, '(A, 2I6)') 'nRealReac, nRealProd:', nRealReactants(i), nRealProducts(i)
        end if
        !
        do j=1, i-1
          if ((sum(abs(reactions(:, j)-reactions(:, i))) .EQ. 0) &
              .AND. (typeReac(j) .EQ. typeReac(i))) then
            write (*,'(2A, 2I6, A)') &
              'Duplicate reaction pair: ', &
              char(27)//'[45m', i, j, char(27)//'[0m'
          end if
        end do
        !
        if (outputFormat .eq. 'Herbst') then
          !
          ! Format back.
          do j=1, nReactants
            call formatBack(strReactants(j, i))
          end do
          do j=1, nProducts
            call formatBack(strProducts(j, i))
          end do
          !
          write (fU, '(I5, X, 6(A12,X), X, ES9.2, X, ES9.2, X, ES9.2, X, I3)') &
            nReacCounter, strReactants(1:2, i), strProducts(1:4, i), &
            dblABC(1:2, i), &
            dblABC(3, i), &
            typeReac(i)
        else
          write (fU, '(7A12, ES9.2, ES9.2, ES9.2, 2I6, I3, X,A1,X,A2,X,A1)') &
            strReactants(:, i), strProducts(:, i), &
            dblABC(1:2, i), &
            dblABC(3, i), &
            int(T_min(i)), int(T_max(i)), typeReac(i), &
            cquality(i), ctype(i), stype(i)
        end if
        !
        nReacCounter = nReacCounter + 1
        !
        if (sum(nElementReac(idxD+1:nElement)) .GE. noDMaxMetal) then
          cycle
        end if
        !
        flag = .FALSE.
        do j=idxD+1, nElement
          if (nElementReac(j) .GT. 0) then
            if (ElementTypicalAbundance(j) .LT. noDEleAbundance) then
              flag = .TRUE.
              exit
            end if
          end if
        end do
        if (flag) cycle
        !
        if (sum(nElementReac(idxD+1:nElement)) .LE. nOtherDeutMax) then
          call DeutReac(i, nDeutDegree, fU)
        else
          call DeutReac(i, 1, fU)
        end if
        !call DeutReac (i, nDeutDegree, fU)
        !
      end do
      close (UNIT=fU, IOSTAT=ios, ERR=999, STATUS='KEEP')
      !
      write (*, FMT='(7(/, A32, I16))') &
        'Number of species: ', nSpecies, &
        'Number of reactions: ', nReactions, &
        'Max number of reactants: ', maxval(nRealReactants), &
        'Max number of products: ', maxval(nRealProducts), &
        'Number of isotopized reactions: ', nReacCounter
      deallocate (nElementReac, nElementProd, STAT=statALLOC)


998   CALL CPU_TIME(ProgCurrentTime)
      write (*, '(/A, F10.3)') 'Seconds elapsed: ', &
        ProgCurrentTime - ProgStartTime

999   write (*,'(A)') &
        char(27)//'[47m'//char(27)//'[34m'//&
        'Program end'//char(27)//'[0m'

      end program


      subroutine initialize (fileInitialize)
      ! Read all the configuration parameters
      use CMDT
      use CMReac
      implicit none
      character(LEN=128) fileInitialize
      integer fU, ios
      logical getFileUnit

      namelist /PhysicalParameters/ &
        Temperature, n_H, &
        Av, omega_Albedo, rateCosIon, MeanMolWeight, &
        BarrierWidth, ratioGrainToH, ratioDustGasMass, &
        GrainRadius, GrainDensity, SitesDensity, &
        nSpecies_Est, nIteration, timeUnit, sto_threshold, &
        rateThreshold, nOrderLim, nDeutDegree, &
        nOtherDeutMax, noDMaxMetal, noDEleAbundance, &
        elementOld, elementNew, &
        inputGrainEleName, inputGrainName, &
        outputFormat, outputGrainEleName, &
        commentChar, inputFormat, grain_special, &
        copy_rates
        !idxNew, idxOld

      namelist /Paths/ &
        path, fReactions, fMobility, fDissRecom, &
        fInitialCondition, &
        fReactionsSave, fSpeciesSave, fRatesSave, &
        fSaveFinalResult, fSaveElementalAb, fSaveALLHistory, &
        fAnalyseResult, fSaveALLHistory_ASCII, fFinalJac, &
        fNameMoments, fEqMoments, fOutLog, fPhyParHistory

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if

      CALL openFileSequentialRead (fU, fileInitialize, 999)
      read (UNIT=fU, IOSTAT=ios, NML=PhysicalParameters)
      read (UNIT=fU, IOSTAT=ios, NML=Paths)
      close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')

      end subroutine initialize
