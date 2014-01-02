!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, either version 3 of the License, or
!!    (at your option) any later version.
!!
!!    Galacticus is distributed in the hope that it will be useful,
!!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!!    GNU General Public License for more details.
!!
!!    You should have received a copy of the GNU General Public License
!!    along with Galacticus.  If not, see <http://www.gnu.org/licenses/>.

!% Contains a module which implements an arbitrary piecewise power-law stellar initial mass function.

module Star_Formation_IMF_PiecewisePowerLaw
  !% Implements an arbitrary piecewise power-law stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_PiecewisePowerLaw, Star_Formation_IMF_Register_Name_PiecewisePowerLaw,&
       & Star_Formation_IMF_Recycled_Instantaneous_PiecewisePowerLaw, Star_Formation_IMF_Yield_Instantaneous_PiecewisePowerLaw,&
       & Star_Formation_IMF_Tabulate_PiecewisePowerLaw, Star_Formation_IMF_Minimum_Mass_PiecewisePowerLaw,&
       & Star_Formation_IMF_Maximum_Mass_PiecewisePowerLaw , Star_Formation_IMF_Phi_PiecewisePowerLaw

  ! Index assigned to this IMF.
  integer                                     :: imfIndex                                 =-1

  ! Flag indicating if the module has been initialized.
  logical                                     :: imfPiecewisePowerLawInitialized          =.false.

  ! Parameters of the IMF.
  double precision                            :: imfPiecewisePowerLawRecycledInstantaneous        , imfPiecewisePowerLawYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer                                     :: imfPieceCount
  double precision, allocatable, dimension(:) :: imfNormalization                                 , massExponent                           , &
       &                                         massLower                                        , massUpper

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_PiecewisePowerLaw</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_PiecewisePowerLaw(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex         =imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_PiecewisePowerLaw

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_PiecewisePowerLaw</unitName>
  !#  <name>PiecewisePowerLaw</name>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_PiecewisePowerLaw(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type     (varying_string), intent(inout) :: imfDescriptors(:), imfNames(:)
    type     (varying_string)                :: imfDescriptor
    character(len=7         )                :: label
    integer                                  :: iPiece

    call Star_Formation_IMF_Initialize_PiecewisePowerLaw
    imfDescriptor="PiecewisePowerLaw"
    do iPiece=1,imfPieceCount
       write (label,'(f7.2)') massLower   (iPiece)
       imfDescriptor=imfDescriptor//":m"//trim(adjustl(label))
       write (label,'(f7.3)') massExponent(iPiece)
       imfDescriptor=imfDescriptor//":e"//trim(adjustl(label))
    end do
    write (label,'(f7.2)') massUpper(imfPieceCount)
    imfDescriptor=imfDescriptor//":m"//trim(adjustl(label))
    imfNames      (imfIndex)="PiecewisePowerLaw"
    imfDescriptors(imfIndex)=imfDescriptor
    return
  end subroutine Star_Formation_IMF_Register_Name_PiecewisePowerLaw

  subroutine Star_Formation_IMF_Initialize_PiecewisePowerLaw
    !% Initialize the PiecewisePowerLaw IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    use Memory_Management
    use Galacticus_Error
    implicit none
    double precision, allocatable, dimension(:) :: massPoints
    logical                                     :: pieceWiseImfIsDefined

    if (.not.imfPiecewisePowerLawInitialized) then
       !$omp critical (IMF_PiecewisePowerLaw_Initialize)
       if (.not.imfPiecewisePowerLawInitialized) then
          !@ <inputParameter>
          !@   <name>imfPiecewisePowerLawRecycledInstantaneous</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>0.39</defaultValue>
          !@   <description>
          !@     The recycled fraction for piecewise power-law stellar initial mass functions in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfPiecewisePowerLawRecycledInstantaneous',imfPiecewisePowerLawRecycledInstantaneous,defaultValue=0.39d0)
          !@ <inputParameter>
          !@   <name>imfPiecewisePowerLawYieldInstantaneous</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>0.02</defaultValue>
          !@   <description>
          !@     The yield for piecewise power-law stellar initial mass functions in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfPiecewisePowerLawYieldInstantaneous'   ,imfPiecewisePowerLawYieldInstantaneous   ,defaultValue=0.02d0)

          ! Get the number of intervals and allocate arrays appropriately.
          imfPieceCount=Get_Input_Parameter_Array_Size('imfPiecewisePowerLawMassPoints')-1
          if (imfPieceCount == -1 ) then
             pieceWiseImfIsDefined=.false.
             imfPieceCount        =1
          else
             pieceWiseImfIsDefined=.true.
             if (imfPieceCount < 1) call Galacticus_Error_Report('Star_Formation_IMF_Initialize_PiecewisePowerLaw','at least 2 mass points are required to define the IMF')
          end if
          call Alloc_Array(massLower       ,[imfPieceCount])
          call Alloc_Array(massUpper       ,[imfPieceCount])
          call Alloc_Array(massExponent    ,[imfPieceCount])
          call Alloc_Array(imfNormalization,[imfPieceCount])
          allocate(massPoints(imfPieceCount+1))

          if (pieceWiseImfIsDefined) then
             ! Read the mass intervals.
             !@ <inputParameter>
             !@   <name>imfPiecewisePowerLawMassPoints</name>
             !@   <attachedTo>module</attachedTo>
             !@   <defaultValue>0.1, 125</defaultValue>
             !@   <description>
             !@     The mass points used to define a piecewise power-law initial mass function.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>initialMassFunction</group>
             !@ </inputParameter>
             call Get_Input_Parameter('imfPiecewisePowerLawMassPoints',massPoints )
             !@ <inputParameter>
             !@   <name>imfPiecewisePowerLawExponents</name>
             !@   <attachedTo>module</attachedTo>
             !@   <defaultValue>-2.35</defaultValue>
             !@   <description>
             !@     The exponents used to define a piecewise power-law initial mass function.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>1..*</cardinality>
             !@   <group>initialMassFunction</group>
             !@ </inputParameter>
             call Get_Input_Parameter('imfPiecewisePowerLawExponents' ,massExponent)
          else
             ! Set defaults (a Salpeter IMF).
             massPoints  =[0.1d0,125.0d0]
             massExponent=[-2.35d0]
          end if

          ! Extract lower and upper limits of the mass ranges.
          massLower=massPoints(1:imfPieceCount  )
          massUpper=massPoints(2:imfPieceCount+1)
          deallocate(massPoints)

          ! Get the normalization for this IMF.
          call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)

          imfPiecewisePowerLawInitialized=.true.
       end if
       !$omp end critical (IMF_PiecewisePowerLaw_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_PiecewisePowerLaw

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_PiecewisePowerLaw</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_PiecewisePowerLaw(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_PiecewisePowerLaw
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_PiecewisePowerLaw

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_PiecewisePowerLaw</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_PiecewisePowerLaw(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_PiecewisePowerLaw
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_PiecewisePowerLaw

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_PiecewisePowerLaw</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_PiecewisePowerLaw(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(in   ) :: initialMass
    double precision, intent(  out) :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_PiecewisePowerLaw
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_PiecewisePowerLaw

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_PiecewisePowerLaw</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_PiecewisePowerLaw(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_PiecewisePowerLaw
       recycledFraction=imfPiecewisePowerLawRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_PiecewisePowerLaw

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_PiecewisePowerLaw</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_PiecewisePowerLaw(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_PiecewisePowerLaw
       yield=imfPiecewisePowerLawYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_PiecewisePowerLaw

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_PiecewisePowerLaw</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_PiecewisePowerLaw(imfSelected,imfMatched,imf)
    !% Register the name of this IMF.
    use Tables
    use Star_Formation_IMF_PPL
    implicit none
    integer                      , intent(in   ) :: imfSelected
    logical                      , intent(inout) :: imfMatched
    class  (table1D), allocatable, intent(inout) :: imf
    integer         , parameter                  :: nPoints    =100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_PiecewisePowerLaw()
       allocate(table1DLogarithmicLinear :: imf)
       select type (imf)
       type is (table1DLogarithmicLinear)
          call imf%create  (                                              &
               &            massLower(1            )                    , &
               &            massUpper(imfPieceCount)                    , &
               &            nPoints                                       &
               &           )
          call imf%populate(                                              &
               &            Piecewise_Power_Law_IMF_Phi(                  &
               &                                        massLower       , &
               &                                        massUpper       , &
               &                                        massExponent    , &
               &                                        imfNormalization, &
               &                                        imf%xs()          &
               &                                       )                  &
               &           )
          imfMatched=.true.
       end select
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_PiecewisePowerLaw

end module Star_Formation_IMF_PiecewisePowerLaw
