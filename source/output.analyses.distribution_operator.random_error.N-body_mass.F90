!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
!!    Andrew Benson <abenson@carnegiescience.edu>
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

  !!{
  Implements a random error output analysis distribution operator class providing errors in $\log_{10}$
  of N-body halo mass.
  !!}

  use :: Statistics_NBody_Halo_Mass_Errors, only : nbodyHaloMassErrorClass

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRndmErrNbodyMass">
   <description>A random error output analysis distribution operator class providing errors in $\log_{10}$ of N-body halo mass.</description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRndmErrNbodyMass
     !!{
     A random error output distribution operator class providing errors in $\log_{10}$ of N-body halo mass.
     !!}
     private
     class(nbodyHaloMassErrorClass), pointer :: nbodyHaloMassError_ => null()
   contains
     final     ::                 randomErrorNbodyMassDestructor
     procedure :: rootVariance => randomErrorNbodyMassRootVariance
  end type outputAnalysisDistributionOperatorRndmErrNbodyMass

  interface outputAnalysisDistributionOperatorRndmErrNbodyMass
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorRndmErrNbodyMass} output analysis distribution operator class.
     !!}
     module procedure randomErrorNbodyMassConstructorParameters
     module procedure randomErrorNbodyMassConstructorInternal
  end interface outputAnalysisDistributionOperatorRndmErrNbodyMass

contains

  function randomErrorNbodyMassConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorRndmErrNbodyMass} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorRndmErrNbodyMass)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (nbodyHaloMassErrorClass                           ), pointer       :: nbodyHaloMassError_

    !![
    <objectBuilder class="nbodyHaloMassError" name="nbodyHaloMassError_" source="parameters"/>
    !!]
    self=outputAnalysisDistributionOperatorRndmErrNbodyMass(nbodyHaloMassError_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="nbodyHaloMassError_"/>
    !!]
    return
  end function randomErrorNbodyMassConstructorParameters

  function randomErrorNbodyMassConstructorInternal(nbodyHaloMassError_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisDistributionOperatorRndmErrNbodyMass} output analysis distribution operator class.
    !!}
    implicit none
    type (outputAnalysisDistributionOperatorRndmErrNbodyMass)                        :: self
    class(nbodyHaloMassErrorClass                           ), intent(in   ), target :: nbodyHaloMassError_
    !![
    <constructorAssign variables="*nbodyHaloMassError_"/>
    !!]

    return
  end function randomErrorNbodyMassConstructorInternal

  subroutine randomErrorNbodyMassDestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisDistributionOperatorRndmErrNbodyMass} output analysis distribution operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorRndmErrNbodyMass), intent(inout) :: self

    !![
    <objectDestructor name="self%nbodyHaloMassError_" />
    !!]
    return
  end subroutine randomErrorNbodyMassDestructor

  double precision function randomErrorNbodyMassRootVariance(self,propertyValue,node)
    !!{
    Computes errors on $\log_{10}($halo masses$)$ for N-body halos.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorRndmErrNbodyMass), intent(inout) :: self
    double precision                                                    , intent(in   ) :: propertyValue
    type            (treeNode                                          ), intent(inout) :: node
    !$GLC attributes unused :: propertyValue

    ! Return the fractional error in halo mass, divided by log(10) to convert from natural to base-10 logarithm.
    randomErrorNbodyMassRootVariance=+self%nbodyHaloMassError_%errorFractional(node  ) &
         &                           /                         log            (10.0d0)
    return
  end function randomErrorNbodyMassRootVariance
