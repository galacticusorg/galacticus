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

!+    Contributions to this file made by: Xiaolong Du, Andrew Benson.

  !!{
  An implementation of exponentially truncated dark matter halo profiles \cite{kazantzidis_2006}.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Mass_Distributions     , only : enumerationNonAnalyticSolversType

  !![
  <darkMatterProfileDMO name="darkMatterProfileDMOTruncatedExponential">
    <description>
      Exponentially truncated dark matter halo profiles \cite{kazantzidis_2006} are constructed via the
      \refClass{massDistributionSphericalTruncatedExponential} class.
    </description>
  </darkMatterProfileDMO>
  !!]
  type, extends(darkMatterProfileDMOClass) :: darkMatterProfileDMOTruncatedExponential
     !!{
     A dark matter halo profile class implementing exponentially truncated dark matter halos.
     !!}
     private
     class           (darkMatterProfileDMOClass        ), pointer :: darkMatterProfileDMO_ => null()
     class           (darkMatterHaloScaleClass         ), pointer :: darkMatterHaloScale_  => null()
     double precision                                             :: radiusFractionalDecay
     type            (enumerationNonAnalyticSolversType)          :: nonAnalyticSolver
   contains
     final     ::        truncatedExponentialDestructor
     procedure :: get => truncatedExponentialGet     
  end type darkMatterProfileDMOTruncatedExponential

  interface darkMatterProfileDMOTruncatedExponential
     !!{
     Constructors for the \refClass{darkMatterProfileDMOTruncatedExponential} dark matter halo profile class.
     !!}
     module procedure truncatedExponentialConstructorParameters
     module procedure truncatedExponentialConstructorInternal
  end interface darkMatterProfileDMOTruncatedExponential

contains

  function truncatedExponentialConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{darkMatterProfileDMOTruncatedExponential} dark matter halo profile class which takes a parameter set as input.
    !!}
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversEncode
    use :: Input_Parameters  , only : inputParameter                     , inputParameters
    implicit none
    type            (darkMatterProfileDMOTruncatedExponential)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass               ), pointer       :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    type            (varying_string                          )                :: nonAnalyticSolver
    double precision                                                          :: radiusFractionalDecay

    !![
    <inputParameter>
      <name>nonAnalyticSolver</name>
      <defaultValue>var_str('fallThrough')</defaultValue>
      <source>parameters</source>
      <description>Selects how solutions are computed when no analytic solution is available. If set to ``{\normalfont \ttfamily fallThrough}'' then the solution ignoring heating is used, while if set to ``{\normalfont \ttfamily numerical}'' then numerical solvers are used to find solutions.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusFractionalDecay</name>
      <defaultValue>1.0d0</defaultValue>
      <source>parameters</source>
      <description>The truncation scale (in units of the virial radius).</description>
    </inputParameter>    
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    !!]
    self=darkMatterProfileDMOTruncatedExponential(radiusFractionalDecay,enumerationNonAnalyticSolversEncode(char(nonAnalyticSolver),includesPrefix=.false.),darkMatterProfileDMO_,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    !!]
    return
  end function truncatedExponentialConstructorParameters

  function truncatedExponentialConstructorInternal(radiusFractionalDecay,nonAnalyticSolver,darkMatterProfileDMO_,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{darkMatterProfileDMOTruncatedExponential} dark matter profile class.
    !!}
    use :: Error             , only : Error_Report
    use :: Mass_Distributions, only : enumerationNonAnalyticSolversIsValid
    implicit none
    type            (darkMatterProfileDMOTruncatedExponential)                        :: self
    class           (darkMatterProfileDMOClass               ), intent(in   ), target :: darkMatterProfileDMO_
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                          , intent(in   )         :: radiusFractionalDecay
    type            (enumerationNonAnalyticSolversType       ), intent(in   )         :: nonAnalyticSolver
    !![
    <constructorAssign variables="radiusFractionalDecay, nonAnalyticSolver, *darkMatterProfileDMO_, *darkMatterHaloScale_"/>
    !!]

    ! Validate.
    if (.not.enumerationNonAnalyticSolversIsValid(nonAnalyticSolver)) call Error_Report('invalid non-analytic solver type'//{introspection:location})  
    return
  end function truncatedExponentialConstructorInternal

  subroutine truncatedExponentialDestructor(self)
    !!{
    Destructor for the \refClass{darkMatterProfileDMOTruncatedExponential} dark matter halo profile class.
    !!}
    implicit none
    type(darkMatterProfileDMOTruncatedExponential), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    <objectDestructor name="self%darkMatterHaloScale_" />
    !!]
    return
  end subroutine truncatedExponentialDestructor

  function truncatedExponentialGet(self,node,weightBy,weightIndex) result(massDistribution_)
    !!{
    Return the dark matter mass distribution for the given {\normalfont \ttfamily node}.
    !!}
    use :: Galactic_Structure_Options, only : componentTypeDarkHalo                        , massTypeDark                   , weightByMass
    use :: Mass_Distributions        , only : massDistributionSphericalTruncatedExponential, kinematicsDistributionTruncated, massDistributionSpherical
    implicit none
    class           (massDistributionClass                   ), pointer                 :: massDistribution_
    type            (kinematicsDistributionTruncated         ), pointer                 :: kinematicsDistribution_
    class           (darkMatterProfileDMOTruncatedExponential), intent(inout)           :: self
    type            (treeNode                                ), intent(inout)           :: node
    type            (enumerationWeightByType                 ), intent(in   ), optional :: weightBy
    integer                                                   , intent(in   ), optional :: weightIndex
    double precision                                                                    :: radiusVirial
    class           (massDistributionClass                   ), pointer                 :: massDistributionDecorated
    !![
    <optionalArgument name="weightBy" defaultsTo="weightByMass" />
    !!]

    ! Assume a null distribution by default.
    massDistribution_ => null()
    ! If weighting is not by mass, return a null profile.
    if (weightBy_ /= weightByMass) return
    ! Create the mass distribution.
    allocate(massDistributionSphericalTruncatedExponential :: massDistribution_)
    select type(massDistribution_)
    type is (massDistributionSphericalTruncatedExponential)
       radiusVirial              =  self%darkMatterHaloScale_ %radiusVirial(node                     )
       massDistributionDecorated => self%darkMatterProfileDMO_%get         (node,weightBy,weightIndex)
       select type (massDistributionDecorated)
       class is (massDistributionSpherical)
          !![
	  <referenceConstruct object="massDistribution_">
	    <constructor>
              massDistributionSphericalTruncatedExponential(                                                               &amp;
	      &amp;                                         radiusTruncateMinimum=                           radiusVirial, &amp;
	      &amp;                                         radiusTruncateDecay  =self%radiusFractionalDecay*radiusVirial, &amp;
	      &amp;                                         nonAnalyticSolver    =self%nonAnalyticSolver                 , &amp;
	      &amp;                                         massDistribution_    =     massDistributionDecorated         , &amp;
              &amp;                                         componentType        =     componentTypeDarkHalo             , &amp;
              &amp;                                         massType             =     massTypeDark                        &amp;
              &amp;                                        )
	    </constructor>
	  </referenceConstruct>
          !!]
       class default
          call Error_Report('expected a spherical mass distribution'//{introspection:location})
       end select
       !![
       <objectDestructor name="massDistributionDecorated"/>
       !!]       
    end select    
    allocate(kinematicsDistribution_)
    !![
    <referenceConstruct object="kinematicsDistribution_">
      <constructor>
        kinematicsDistributionTruncated( &amp;
	 &amp;                         )
      </constructor>
    </referenceConstruct>
    !!]
    call massDistribution_%setKinematicsDistribution(kinematicsDistribution_)
    !![
    <objectDestructor name="kinematicsDistribution_"/>
    !!]
    return
  end function truncatedExponentialGet
