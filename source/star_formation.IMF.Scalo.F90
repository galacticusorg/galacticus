!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the Scalo stellar initial mass function \citep{scalo_stellar_1986}.

module Star_Formation_IMF_Scalo
  !% Implements the Scalo stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Scalo, Star_Formation_IMF_Register_Name_Scalo,&
       & Star_Formation_IMF_Recycled_Instantaneous_Scalo, Star_Formation_IMF_Yield_Instantaneous_Scalo,&
       & Star_Formation_IMF_Tabulate_Scalo, Star_Formation_IMF_Minimum_Mass_Scalo, Star_Formation_IMF_Maximum_Mass_Scalo&
       &, Star_Formation_IMF_Phi_Scalo

  ! Index assigned to this IMF.
  integer                                    :: imfIndex                     =-1

  ! Flag indicating if the module has been initialized.
  logical                                    :: imfScaloInitialized          =.false.

  ! Parameters of the IMF.
  double precision                           :: imfScaloRecycledInstantaneous                                           , imfScaloYieldInstantaneous

  ! Fixed parameters of the IMF.
  integer         , parameter                :: imfPieceCount                =6
  double precision, dimension(imfPieceCount) :: imfNormalization                                                        , massExponent              =[1.60d0,-1.01d0,-2.75d0,-2.08d0,-3.50d0,-2.63d0], &
       &                                        massLower                    =[0.1d0,0.18d0,0.42d0,0.62d0,1.18d0,3.50d0], massUpper                 =[0.18d0,0.42d0,0.62d0,1.18d0,3.50d0,125.0d0]

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Scalo</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Scalo(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Scalo

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Scalo</unitName>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Scalo(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type(varying_string), intent(inout) :: imfDescriptors(:), imfNames(:)

    imfNames      (imfIndex)="Scalo"
    imfDescriptors(imfIndex)="Scalo"
    return
  end subroutine Star_Formation_IMF_Register_Name_Scalo

  subroutine Star_Formation_IMF_Initialize_Scalo
    !% Initialize the Scalo IMF module.
    use Input_Parameters
    use Star_Formation_IMF_PPL
    implicit none

    if (.not.imfScaloInitialized) then
       !$omp critical (IMF_Scalo_Initialize)
       if (.not.imfScaloInitialized) then
          !@ <inputParameter>
          !@   <name>imfScaloRecycledInstantaneous</name>
          !@   <defaultValue>0.24 (computed internally)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The recycled fraction for the Scalo \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfScaloRecycledInstantaneous',imfScaloRecycledInstantaneous,defaultValue=0.24d0)
          !@ <inputParameter>
          !@   <name>imfScaloYieldInstantaneous</name>
          !@   <defaultValue>0.086 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The yield for the Scalo \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfScaloYieldInstantaneous'   ,imfScaloYieldInstantaneous   ,defaultValue=0.0086d0)

          ! Get the normalization for this IMF.
          call Piecewise_Power_Law_IMF_Normalize(massLower,massUpper,massExponent,imfNormalization)

          imfScaloInitialized=.true.
       end if
       !$omp end critical (IMF_Scalo_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_Scalo

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Scalo</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Scalo(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Scalo
       minimumMass=massLower(1)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Scalo

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Scalo</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Scalo(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Scalo
       maximumMass=massUpper(imfPieceCount)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Scalo

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Scalo</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Scalo(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(in   ) :: initialMass
    double precision, intent(  out) :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Scalo
       imfPhi=Piecewise_Power_Law_IMF_Phi(massLower,massUpper,massExponent,imfNormalization,initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Scalo

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Scalo</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Scalo(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Scalo
       recycledFraction=imfScaloRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Scalo

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Scalo</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Scalo(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Scalo
       yield=imfScaloYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Scalo

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Scalo</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Scalo(imfSelected,imfMatched,imf)
    !% Register the name of this IMF.
    use Star_Formation_IMF_PPL
    use Tables
    implicit none
    integer                      , intent(in   ) :: imfSelected
    logical                      , intent(inout) :: imfMatched
    class  (table1D), allocatable, intent(inout) :: imf
    integer         , parameter                  :: nPoints    =100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Scalo()
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
  end subroutine Star_Formation_IMF_Tabulate_Scalo

end module Star_Formation_IMF_Scalo
