!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements the Chabrier stellar initial mass function \citep{chabrier_galactic_2001}.

module Star_Formation_IMF_Chabrier
  !% Implements the Chabrier stellar initial mass function.
  implicit none
  private
  public :: Star_Formation_IMF_Register_Chabrier, Star_Formation_IMF_Register_Name_Chabrier,&
       & Star_Formation_IMF_Recycled_Instantaneous_Chabrier, Star_Formation_IMF_Yield_Instantaneous_Chabrier,&
       & Star_Formation_IMF_Tabulate_Chabrier, Star_Formation_IMF_Minimum_Mass_Chabrier, Star_Formation_IMF_Maximum_Mass_Chabrier&
       &, Star_Formation_IMF_Phi_Chabrier

  ! Index assigned to this IMF.
  integer           :: imfIndex                        =-1

  ! Flag indicating if the module has been initialized.
  logical          :: imfChabrierInitialized          =.false.

  ! Parameters of the IMF.
  double precision :: imfChabrierRecycledInstantaneous        , imfChabrierYieldInstantaneous
  double precision :: imfChabrierMassLower                    , imfChabrierMassTransition    , &
       &              imfChabrierMassUpper               
  double precision :: imfChabrierExponent                     , imfChabrierMass              , &
       &              imfChabrierSigma                   
  double precision :: normalizationExponential                , normalizationLogNormal       

contains

  !# <imfRegister>
  !#  <unitName>Star_Formation_IMF_Register_Chabrier</unitName>
  !# </imfRegister>
  subroutine Star_Formation_IMF_Register_Chabrier(imfAvailableCount)
    !% Register this IMF by incrementing the count and keeping a record of the assigned index.
    implicit none
    integer, intent(inout) :: imfAvailableCount

    imfAvailableCount=imfAvailableCount+1
    imfIndex=imfAvailableCount
    return
  end subroutine Star_Formation_IMF_Register_Chabrier

  !# <imfRegisterName>
  !#  <unitName>Star_Formation_IMF_Register_Name_Chabrier</unitName>
  !#  <name>Chabrier</name>
  !# </imfRegisterName>
  subroutine Star_Formation_IMF_Register_Name_Chabrier(imfNames,imfDescriptors)
    !% Register the name of this IMF.
    use ISO_Varying_String
    implicit none
    type     (varying_string), intent(inout) :: imfDescriptors(:), imfNames(:)
    character(len=7         )                :: label
    
    call Star_Formation_IMF_Initialize_Chabrier()
    imfNames      (imfIndex)="Chabrier"
    imfDescriptors(imfIndex)="Chabrier"
    write (label,'(f7.3)') imfChabrierMassLower
    imfDescriptors(imfIndex)=imfDescriptors(imfIndex)//":mLow"  //trim(adjustl(label))
    write (label,'(f7.3)') imfChabrierMassUpper
    imfDescriptors(imfIndex)=imfDescriptors(imfIndex)//":mUp"   //trim(adjustl(label))
    write (label,'(f7.3)') imfChabrierMassTransition
    imfDescriptors(imfIndex)=imfDescriptors(imfIndex)//":mTrans"//trim(adjustl(label))
    write (label,'(f7.3)') imfChabrierMass
    imfDescriptors(imfIndex)=imfDescriptors(imfIndex)//":m"     //trim(adjustl(label))
    write (label,'(f7.3)') imfChabrierSigma
    imfDescriptors(imfIndex)=imfDescriptors(imfIndex)//":sigma" //trim(adjustl(label))
    write (label,'(f7.3)') imfChabrierExponent
    imfDescriptors(imfIndex)=imfDescriptors(imfIndex)//":alpha" //trim(adjustl(label))
    return
  end subroutine Star_Formation_IMF_Register_Name_Chabrier

  subroutine Star_Formation_IMF_Initialize_Chabrier
    !% Initialize the Chabrier IMF module.
    use Input_Parameters
    use Error_Functions
    use Numerical_Constants_Math
    implicit none
    double precision :: normalization

    if (.not.imfChabrierInitialized) then
       !$omp critical (IMF_Chabrier_Initialize)
       if (.not.imfChabrierInitialized) then
          !@ <inputParameter>
          !@   <name>imfChabrierRecycledInstantaneous</name>
          !@   <defaultValue>0.46 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The recycled fraction for the Chabrier \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierRecycledInstantaneous',imfChabrierRecycledInstantaneous,defaultValue=0.46d0)
          !@ <inputParameter>
          !@   <name>imfChabrierYieldInstantaneous</name>
          !@   <defaultValue>0.035 (internally computed)</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The yield for the Chabrier \gls{imf} in the instantaneous recycling approximation.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierYieldInstantaneous'   ,imfChabrierYieldInstantaneous   ,defaultValue=0.035d0)
          !@ <inputParameter>
          !@   <name>imfChabrierMassUpper</name>
          !@   <defaultValue>125</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The upper mass limit for the Chabrier \gls{imf}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierMassUpper'            ,imfChabrierMassUpper            ,defaultValue=125.0d0)
          !@ <inputParameter>
          !@   <name>imfChabrierMassLower</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The lower mass limit for the Chabrier \gls{imf}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierMassLower'            ,imfChabrierMassLower            ,defaultValue=  0.1d0)
          !@ <inputParameter>
          !@   <name>imfChabrierMassTransition</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The transition limit for the Chabrier \gls{imf}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierMassTransition'       ,imfChabrierMassTransition       ,defaultValue=  1.0d0)
          !@ <inputParameter>
          !@   <name>imfChabrierSigma</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The width of the lognormal part of the Chabrier \gls{imf}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierSigma'                ,imfChabrierSigma                ,defaultValue=  0.69d0)
          !@ <inputParameter>
          !@   <name>imfChabrierExponent</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The exponent of the power law part of the Chabrier \gls{imf}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierExponent'            ,imfChabrierExponent              ,defaultValue=  -2.3d0)
          !@ <inputParameter>
          !@   <name>imfChabrierMass</name>
          !@   <defaultValue>0.08</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Characteristic mass of the lognormal part of the Chabrier \gls{imf}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfChabrierMass'                ,imfChabrierMass                  ,defaultValue= 0.08d0)
          ! Compute normalizations.
          normalizationLogNormal  =+sqrt(Pi/2.0d0)                                     &
                &                  *imfChabrierSigma                                   &
                &                  *imfChabrierMass                                    &
                &                  *log(10.0d0)                                        &
                &                  *exp(                                               &
                &                       +0.5d0                                         &
                &                       *imfChabrierSigma**2                           &
                &                       *log(10.0d0)     **2                           &
                &                      )                                               &
                &                  *(                                                  &
                &                    -Error_Function(                                  &
                &                                    +imfChabrierSigma                 &
                &                                    *log (10.0d0)                     &
                &                                    /sqrt( 2.0d0)                     &
                &                                    -log10(                           &
                &                                           +imfChabrierMassTransition &
                &                                           /imfChabrierMass           &
                &                                          )                           &
                &                                    /sqrt( 2.0d0)                     &
                &                                    /imfChabrierSigma                 &
                &                                   )                                  &
                &                    +Error_Function(                                  &
                &                                    +imfChabrierSigma                 &
                &                                    *log (10.0d0)                     &
                &                                    /sqrt( 2.0d0)                     &
                &                                    -log10(                           &
                &                                           +imfChabrierMassLower      &
                &                                           /imfChabrierMass           &
                &                                          )                           &
                &                                    /sqrt( 2.0d0)                     &
                &                                    /imfChabrierSigma                 &
                &                                   )                                  &
                &                   )
          normalizationExponential=+exp(                                                       &
               &                        -0.50d0                                                &
               &                        *log10(                                                &
               &                               +imfChabrierMassTransition                      &
               &                               /imfChabrierMass                                &
               &                              )                          **2                   &
               &                        /imfChabrierSigma                **2                   &
               &                       )                                                       &
               &                   *imfChabrierMassTransition                                  &
               &                   /                               (2.0d0+imfChabrierExponent) &
               &                   *(                                                          &
               &                     +(                                                        &
               &                       +imfChabrierMassUpper                                   &
               &                       /imfChabrierMassTransition                              &
               &                      )                          **(2.0d0+imfChabrierExponent) &
               &                     -1.0d0                                                    &
               &                    )
          normalization           =normalizationLogNormal+normalizationExponential
          normalizationLogNormal  =1.0d0/normalization
          normalizationExponential=+exp(                                     &
               &                        -0.50d0                              &
               &                        *log10(                              &
               &                               +imfChabrierMassTransition    &
               &                               /imfChabrierMass              &
               &                              )                          **2 &
               &                        /imfChabrierSigma                **2 &
               &                       )                                     &
               &                   /imfChabrierMassTransition                &
               &                   /normalization
          ! Record that initialization is done.
          imfChabrierInitialized=.true.
       end if
       !$omp end critical (IMF_Chabrier_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize_Chabrier

  !# <imfMinimumMass>
  !#  <unitName>Star_Formation_IMF_Minimum_Mass_Chabrier</unitName>
  !# </imfMinimumMass>
  subroutine Star_Formation_IMF_Minimum_Mass_Chabrier(imfSelected,imfMatched,minimumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: minimumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       minimumMass=imfChabrierMassLower
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Minimum_Mass_Chabrier

  !# <imfMaximumMass>
  !#  <unitName>Star_Formation_IMF_Maximum_Mass_Chabrier</unitName>
  !# </imfMaximumMass>
  subroutine Star_Formation_IMF_Maximum_Mass_Chabrier(imfSelected,imfMatched,maximumMass)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: maximumMass

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       maximumMass=imfChabrierMassUpper
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Maximum_Mass_Chabrier

  !# <imfPhi>
  !#  <unitName>Star_Formation_IMF_Phi_Chabrier</unitName>
  !# </imfPhi>
  subroutine Star_Formation_IMF_Phi_Chabrier(imfSelected,imfMatched,initialMass,imfPhi)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(in   ) :: initialMass
    double precision, intent(  out) :: imfPhi

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       imfPhi=Chabrier_Phi(initialMass)
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Phi_Chabrier

  !# <imfRecycledInstantaneous>
  !#  <unitName>Star_Formation_IMF_Recycled_Instantaneous_Chabrier</unitName>
  !# </imfRecycledInstantaneous>
  subroutine Star_Formation_IMF_Recycled_Instantaneous_Chabrier(imfSelected,imfMatched,recycledFraction)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: recycledFraction

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       recycledFraction=imfChabrierRecycledInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Recycled_Instantaneous_Chabrier

  !# <imfYieldInstantaneous>
  !#  <unitName>Star_Formation_IMF_Yield_Instantaneous_Chabrier</unitName>
  !# </imfYieldInstantaneous>
  subroutine Star_Formation_IMF_Yield_Instantaneous_Chabrier(imfSelected,imfMatched,yield)
    !% Register the name of this IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    logical         , intent(inout) :: imfMatched
    double precision, intent(  out) :: yield

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier
       yield=imfChabrierYieldInstantaneous
       imfMatched=.true.
    end if
    return
  end subroutine Star_Formation_IMF_Yield_Instantaneous_Chabrier

  !# <imfTabulate>
  !#  <unitName>Star_Formation_IMF_Tabulate_Chabrier</unitName>
  !# </imfTabulate>
  subroutine Star_Formation_IMF_Tabulate_Chabrier(imfSelected,imfMatched,imf)
    !% Register the name of this IMF.
    use Tables
    implicit none
    integer                      , intent(in   ) :: imfSelected
    logical                      , intent(inout) :: imfMatched
    class  (table1D), allocatable, intent(inout) :: imf
    integer         , parameter                  :: nPoints    =100

    if (imfSelected == imfIndex) then
       call Star_Formation_IMF_Initialize_Chabrier()
       allocate(table1DLogarithmicLinear :: imf)
       select type (imf)
       type is (table1DLogarithmicLinear)
          call imf%create  (                      &
               &            imfChabrierMassLower, &
               &            imfChabrierMassUpper, &
               &            nPoints               &
               &           )
          call imf%populate(Chabrier_Phi(imf%xs()))
          imfMatched=.true.
       end select
    end if
    return
  end subroutine Star_Formation_IMF_Tabulate_Chabrier

  elemental double precision function Chabrier_Phi(mass)
    !% Evaluates the Chabrier initial mass function.
    implicit none
    double precision, intent(in   ) :: mass

    if (mass >= imfChabrierMassLower .and. mass < imfChabrierMassTransition) then
       Chabrier_Phi=normalizationLogNormal*exp(-0.5d0*(log10(mass/imfChabrierMass)/imfChabrierSigma)**2)/mass
    else if (mass >= imfChabrierMassTransition .and. mass < imfChabrierMassUpper) then
       Chabrier_Phi=normalizationExponential*(mass**imfChabrierExponent)
    else
       Chabrier_Phi=0.0d0
    end if
    return
  end function Chabrier_Phi

end module Star_Formation_IMF_Chabrier
