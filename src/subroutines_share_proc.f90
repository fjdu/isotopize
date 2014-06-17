      subroutine CalcRates
      use CMDT
      implicit none
      integer i
      double precision dblTmp
!     All the formulae used here can be found in Hasegawa et al. 1992, 1993.
!       Notes added on 2011-03-30 Wed 15:06:30.

      T300 = Temperature / 300.0D0
!      TemperatureReduced = kBoltzmann * Temperature &
!        / (elementaryCharge**2 * coulombConst / (GrainRadius*1D-6))
!      JNegaPosi = (1D0 + 1D0/TemperatureReduced) &
!                * (1D0 + SQRT(2D0/(2D0+TemperatureReduced)))
!      JChargeNeut = (1D0 + SQRT(DPI/2D0/TemperatureReduced))

      if (ratioGrainToH .LT. 1D-50) &
        ratioGrainToH = & ! ratioGrainToH := N_Grain / N_H
          ratioDustGasMass * (mProton * MeanMolWeight) &
          / (4.0D0*DPI/3.0D0 * (GrainRadius)**3 * GrainDensity)
      SitesPerGrain = 4D0*DPI * GrainRadius**2 * (SitesDensity*1D4)
      VolumnPerGrain = 1D0 / (ratioGrainToH * n_H) ! in cm3

      do i=1, nReactions
        if ((nRealReactants(i) .LE. 0) .OR. &
            (nRealProducts(i) .LE. 0)) then
          rates(i) = 0.0D0
          cycle
        end if

        select case (typeReac(i))
        case (5, 53) ! Two body
          rates(i) = dblABC(1, i) * (T300**dblABC(2, i)) &
              * exp(-dblABC(3, i)/Temperature) &
              / VolumnPerGrain
        case (1) ! One body
          rates(i) = dblABC(1, i)
        case (2) ! One body
          rates(i) = dblABC(1, i) * (T300**dblABC(2, i)) &
              * dblABC(3, i) / (1-omega_Albedo)
        case (3) ! One body
          rates(i) = dblABC(1, i) * exp(-dblABC(3, i) * Av)
        case (61)
          ! rates(i) * Population =  number of i molecules 
          !   accreted per grain per unit time
          ! Pi * r**2 * V * n
          ! Also take into account the possible effect of
          !   stick coefficient and other temperature dependence 
          !   (e.g., coloumb focus).
          rates(i) = dblABC(1,i) * Temperature**dblABC(2,i) &
            * exp(-dblABC(3,i)/Temperature) &
            * DPI * GrainRadius**2 &
            * sqrt(8D0/DPI*kBoltzmann*Temperature &
            / (massSpecies(reactions(1, i)) * mProton)) &
            / (VolumnPerGrain*1D-6) ! VolumnPerGrain is in cm3.
        case (62)
          rates(i) = dblABC(1,i) * Temperature**dblABC(2,i) &
            * exp(-dblABC(3,i)/Temperature) &
            + 3.16D-19 * & ! Cosmic ray desorption rate
            dblABC(1,i) * (70D0)**dblABC(2,i) & !From Hasegawa1993
            * exp(-dblABC(3,i)/70D0)
             !exp(-2D0 * BarrierWidth / hbarPlanck &
             ! * sqrt(2D0 * massSpecies(reactions(1,i)) * mProton &
             ! * kBoltzmann * dblABC(3,i)))
             ! Quantum tunneling not included.
             ! The value of BarrierWidth is uncertain.
             ! Notice the difference bwtween hbarPlanck and hPlanck.
        case (63)
          ! Use the reduced mass in two-body collision reaction.
          dblTmp = massSpecies(reactions(1,i)) &
            * massSpecies(reactions(2,i)) &
            / (massSpecies(reactions(1,i)) &
            + massSpecies(reactions(2,i))) * mProton
          rates(i) = mobilitySpecies(reactions(1,i)) &
            * dblABC(1,i) * Temperature**dblABC(2,i) &
            * exp(-2D0 * BarrierWidth / hbarPlanck &
               * sqrt(2D0 * dblTmp * kBoltzmann * dblABC(3,i))) &
            / SitesPerGrain
        case (64)
          dblTmp = massSpecies(reactions(1,i)) &
            * massSpecies(reactions(2,i)) &
            / (massSpecies(reactions(1,i)) &
            + massSpecies(reactions(2,i))) * mProton
          rates(i) = &
            (mobilitySpecies(reactions(1,i)) &
            + mobilitySpecies(reactions(2,i))) &
            * dblABC(1,i) * Temperature**dblABC(2,i) &
            * exp(-2D0 * BarrierWidth / hbarPlanck &
               * sqrt(2D0 * dblTmp * kBoltzmann * dblABC(3,i))) &
            / SitesPerGrain
        case default
          rates(i) = 0D0
        end select
      end do

      end subroutine CalcRates




      subroutine CalcMobilities
      use CMDT
      implicit none
      integer i, j, i1, nLineAll, fU, ios
      double precision, dimension(3) :: dblTmpVec
      character(Len=128) FMTstr, strTMP
      character commentChar
      character(LEN=constLenNameSpecies) nameSpecies_tmp
      logical IsWordChar, getFileUnit

      mobilitySpecies = 0D0
      i1 = 0
      commentChar = '!'

      if (.NOT. getFileUnit(fU)) then
        write (*,*) 'Cannot allocate an output file unit!'
        stop
      end if
      if (IsWordChar(fMobility(1:1))) then
        CALL openFileSequentialRead &
            (fU, trim(path)//trim(fMobility), 99999)
        CALL GetNLineF (fU, nLineAll, nMobilityList, commentChar)
        !write (*, '(/A32, I16)') 'Number of mobile species: ', &
        !  nMobilityList

        rewind (UNIT=fU, IOSTAT=ios)
        write (FMTstr, FMT=&
          '("(", "A", I2, ", 6X, ", I1, "E", I1, ".", I1, ")")') &
          constLenNameSpecies, 3, 9, 2

        do i=1, nLineAll
          read (UNIT=fU, FMT='(A128)', IOSTAT=ios) strTMP
          if ((strTMP(1:1) .EQ. commentChar) .OR. &
              (strTMP(1:1) .EQ. ' ')) cycle
          read (strTMP, FMT=FMTstr, IOSTAT=ios) &
            nameSpecies_tmp, dblTmpVec
          if (ios .NE. 0) then
            write (*, *) 'Error in reading mobility, ios = ', ios
            stop
          end if
          do j=1, nSpecies
            if (trim(nameSpecies(j)) .EQ. &
              trim(nameSpecies_tmp)) then  ! The standard formula
              if (massSpecies(j) .GE. 2D0) then
                mobilitySpecies(j) &
                  = dblTmpVec(1) * Temperature**dblTmpVec(2) &
                    * exp(-dblTmpVec(3)/Temperature)
              else ! Quantum tunelling of H (and only H) on the grain
                mobilitySpecies(j) &
                  = dblTmpVec(1) * Temperature**dblTmpVec(2) &
                    * exp(-2D0 * BarrierWidth / hbarPlanck &
                    * sqrt(2D0 * massSpecies(j) * mProton &
                    * kBoltzmann * dblTmpVec(3)))
              end if
              i1 = i1 + 1
              exit
            end if
          end do
        end do
        close (UNIT=fU, IOSTAT=ios, STATUS='KEEP')
        !write (*, '(A32, I16)') 'Really used: ', i1
      end if
      end subroutine CalcMobilities





      subroutine ReadReactions (fU, nLineAll, commentChar, inputFormat)
      ! Read the reaction file.
      use CMDT
      implicit none
      integer i, j, k, fU, ios, nLineAll
      character(LEN=256) strtmp
      character(LEN=64) FMTstring
      character, intent(in) :: commentChar
      character(len=*), intent(in), optional :: inputFormat
      character(len=16) inputF
      !
      if (present(inputFormat)) then
        inputF = inputFormat
      else
        inputF = 'nonHerbst'
      end if
      !
      nRealReactants = 0
      nRealProducts = 0
      !
      if (inputF .eq. 'Herbst') then
        FMTstring = '(6X, 6(A12,X), X, 3(F9.0, X), I3)'
      else
        FMTstring = '(7A12, 3F9.0, 2F6.0, I3, X,A1, X,A2, X,A1)'
      end if
      !
      write(*,*) 'Reading with format: ', FMTstring
      !
      rewind (UNIT=fU, IOSTAT=ios)
      !
      i = 1
      !
      do
        read (UNIT=fU, FMT='(A256)', IOSTAT=ios) strtmp
        if (ios .LT. 0) exit
        !
        if (inputF .eq. 'Herbst') then
          ! Commented or empty lines are treated not the same.
          if (strtmp(1:1) .EQ. commentChar) then
            cycle
          end if
          read (UNIT=strtmp, FMT=FMTstring, IOSTAT=ios) &
            strReactants(:,i), strProducts(:,i), dblABC(:,i), &
            typeReac(i)
           T_min(i)=0D0; T_max(i)=0D0; cquality(i)=''; ctype(i)=''; stype(i)=''
          !write (*, FMTstring, IOSTAT=ios) &
          !  strReactants(:,i), strProducts(:,i), dblABC(:,i), &
          !  typeReac(i)
          !write(*,'(7A12)') strReactants(:,i), ' -> ', strProducts(:,i)
        else
          ! Commented or empty lines are treated the same.
          if ((strtmp(1:1) .EQ. commentChar) .OR. &
              (strtmp(1:1) .EQ. ' ')) then
            cycle
          end if
          read (UNIT=strtmp, FMT=FMTstring, IOSTAT=ios) &
            strReactants(:,i), strProducts(:,i), dblABC(:,i), &
            T_min(i), T_max(i), typeReac(i), cquality(i), ctype(i), stype(i)
        end if
        !
        do j=1, nReactants
          do k=1, constLenNameSpecies
            if (strReactants(j, i)(k:k) .NE. ' ') then
              nRealReactants(i) = nRealReactants(i) + 1
              exit
            end if
          end do
          if (inputF .eq. 'Herbst') then
            if (trim(strReactants(j, i)) .EQ. 'E') then
              strReactants(j, i) = 'E-'
            end if
          end if
          if (trim(strReactants(j, i)) .EQ. 'PHOTON') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
          if (trim(strReactants(j, i)) .EQ. 'CRPHOT') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
          if (trim(strReactants(j, i)) .EQ. 'CRP') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
          if (trim(strReactants(j, i)) .EQ. 'C-RAY') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
          if (trim(strReactants(j, i)) .EQ. 'TEMP') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
          if (trim(strReactants(j, i)) .EQ. 'LYAPHOTON') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
          if (trim(strReactants(j, i)) .EQ. 'XRAY') then
            nRealReactants(i) = nRealReactants(i) - 1
          end if
        end do
        do j=1, nProducts
          do k=1, constLenNameSpecies
            if (strProducts(j, i)(k:k) .NE. ' ') then
              nRealProducts(i) = nRealProducts(i) + 1
              exit
            end if
          end do
          if (inputF .eq. 'Herbst') then
            if (trim(strProducts(j, i)) .EQ. 'E') then
              strProducts(j, i) = 'E-'
            end if
          end if
          if (trim(strProducts(j, i)) .EQ. 'PHOTON') then
            nRealProducts(i) = nRealProducts(i) - 1
          end if
        end do
        i = i + 1
      end do
      end subroutine ReadReactions




      subroutine MakeReactions
      use CMDT
      implicit none

      integer i, j, k
      logical flag

      reactions = 0

      nameSpecies(1) = strReactants(1, 1)
      nSpecies = 1
      do i=1, nReactions
        do k=1, nRealReactants(i)
          flag = .TRUE.
          do j=1, nSpecies
            if (trim(nameSpecies(j)) .EQ. trim(strReactants(k, i))) then
              flag = .FALSE.
              reactions(k, i) = j
              exit
            end if
          end do
          if (flag) then
            nSpecies = nSpecies + 1
            nameSpecies(nSpecies) = strReactants(k, i)
            reactions(k, i) = nSpecies
          end if
        end do
        do k=1, nRealProducts(i)
          flag = .TRUE.
          do j=1, nSpecies
            if (trim(nameSpecies(j)) .EQ. trim(strProducts(k, i))) then
              flag = .FALSE.
              reactions(k+nReactants, i) = j
              exit
            end if
          end do
          if (flag) then
            nSpecies = nSpecies + 1
            nameSpecies(nSpecies) = strProducts(k, i)
            reactions(k+nReactants, i) = nSpecies
          end if
        end do
      end do
      end subroutine MakeReactions



      function getPhyPar (t)
      ! Input
      !         t:  time. The unit should be in year.
      ! Output
      !         getPhyPar:
      !              .TRUE.: at least one physical parameter has changed
      !             .FALSE.: no physical parameter has changed
      ! If there are some parameters don't change with time, then
      !   just ignore them.
      use CMDT
      implicit none
      logical getPhyPar
      double precision t, Temperature_Tmp
      double precision, parameter :: t0 = 3.0D5, t1 = 3.001D5
      double precision, parameter :: &
        Temperature0 = 1.0D1, Temperature1 = 6.0D1, ratioChange=1D-1
      Temperature_Tmp = Temperature
      !
      if (t .LT. t0) then
        Temperature = Temperature0
      else if (t .LT. t1) then
        Temperature = Temperature0 + &
          (t-t0)/(t1-t0) * (Temperature1 - Temperature0)
      else
        Temperature = Temperature1
      end if
      !
      if (abs(Temperature_Tmp-Temperature) .GE. &
        Temperature_Tmp*ratioChange) then
        getPhyPar = .TRUE.
      else
        Temperature = Temperature_Tmp
        getPhyPar = .FALSE.
      end if
      return
      end function getPhyPar


      function ExternalAction (t, y, NEQ)
      ! Some external actions which modify the system status
      ! Input
      !         t:  time. The unit should be in year.
      ! Output
      !         ExternalAction:
      !              .TRUE.: something important has changed
      !             .FALSE.: nothing important has changed
      use CMDT
      implicit none
      logical ExternalAction
      double precision, parameter :: t0 = 3.0D5, t1 = 3.03D5
      integer i, i1, NEQ
      double precision t, y(NEQ)
      if ((t .GE. t0) .AND. (t .LT. t1)) then
        do i=1, nGrReactions
          i1 = idxGrReactions(i)
          if (typeReac(i1) .EQ. 62) then
            y(reactions(1+nReactants, i1)) = &
              y(reactions(1+nReactants, i1)) + y(reactions(1, i1))
            y(reactions(1, i1)) = 0D0
          end if
        end do
        do i=nSpecies+1, NEQ
          y(i) = 0D0
        end do
        ExternalAction = .TRUE.
      else
        ExternalAction = .FALSE.
      end if
      return
      end function ExternalAction
