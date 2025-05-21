!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements an N-body dark matter halo mass error class using the model of \cite{trenti_how_2010}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorTrenti2010">
   <description>An N-body dark matter halo mass error class using the model of \cite{trenti_how_2010}.</description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorPowerLaw) :: nbodyHaloMassErrorTrenti2010
     !!{
     An N-body halo mass error class using the model of \cite{trenti_how_2010}.
     !!}
     private
     ! Parameters of the correlation model.
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_         => null()
     double precision                                   :: correlationNormalization             , correlationMassExponent, &
          &                                                correlationRedshiftExponent          , massParticle
   contains
     final     ::                trenti2010Destructor
     procedure :: correlation => trenti2010Correlation
  end type nbodyHaloMassErrorTrenti2010

  interface nbodyHaloMassErrorTrenti2010
     !!{
     Constructors for the \refClass{nbodyHaloMassErrorTrenti2010} N-body halo mass error class.
     !!}
     module procedure nbodyHaloMassErrorTrenti2010Parameters
     module procedure nbodyHaloMassErrorTrenti2010Internal
  end interface nbodyHaloMassErrorTrenti2010

contains

  function nbodyHaloMassErrorTrenti2010Parameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nbodyHaloMassErrorTrenti2010} N-body halo mass error class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyHaloMassErrorTrenti2010)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctions_
    double precision                                              :: massParticle           , correlationNormalization   , &
         &                                                           correlationMassExponent, correlationRedshiftExponent

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>massParticle</name>
      <source>parameters</source>
      <variable>massParticle</variable>
      <description>The mass of the particle in the N-body simulation in which halos were found.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationNormalization</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <variable>correlationNormalization</variable>
      <description>Variable $C_0$ in the model for the correlation between halo mass errors: $C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta$.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationMassExponent</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <variable>correlationMassExponent</variable>
      <description>Variable $\alpha$ in the model for the correlation between halo mass errors: $C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta$.</description>
    </inputParameter>
    <inputParameter>
      <name>correlationRedshiftExponent</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <variable>correlationRedshiftExponent</variable>
      <description>Variable $\beta$ in the model for the correlation between halo mass errors: $C_{12} = C_0 [M_2/M_1]^\alpha [(1+z_2)/(1+z_1)]^\beta$.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nbodyHaloMassErrorTrenti2010(massParticle,correlationNormalization,correlationMassExponent,correlationRedshiftExponent,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function nbodyHaloMassErrorTrenti2010Parameters

  function nbodyHaloMassErrorTrenti2010Internal(massParticle,correlationNormalization,correlationMassExponent,correlationRedshiftExponent,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{nbodyHaloMassErrorTrenti2010} N-body halo mass error class. \cite{trenti_how_2010} report
    a normalization of the fractional error in particle number of 0.15 at $N=1000$ particles. Since this is based on
    comparisons of halos in simulations differing in number of particles by a factor $8$ this actually overestimates the
    normalization by a factor $\sqrt{5/4}$. Therefore, we use a normalization of $0.135$ here.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (nbodyHaloMassErrorTrenti2010)                                  :: self
    double precision                              , intent(in   )                   :: massParticle
    double precision                              , intent(in   ), optional         :: correlationNormalization               , correlationMassExponent, &
         &                                                                             correlationRedshiftExponent
    class           (cosmologyFunctionsClass     ), intent(in   ), optional, target :: cosmologyFunctions_
    double precision                              , parameter                       :: exponent                =-1.000d0/3.0d0
    double precision                              , parameter                       :: normalization           =+0.135d0
    double precision                              , parameter                       :: particleNumberReference =+1.000d3

    self%massParticle                  =massParticle
    self%normalizationSquared          =(normalization*(massReference/particleNumberReference/massParticle)**exponent)**2
    self%exponent                      =                                                                     exponent
    self%fractionalErrorHighMassSquared=+0.0d0
    ! Set correlation properties.
    if (present(correlationNormalization).or.present(correlationMassExponent).or.present(correlationRedshiftExponent)) then
       if (.not.(present(correlationNormalization).and.present(correlationMassExponent).and.present(correlationRedshiftExponent))) &
            & call Error_Report('all parameters of correlation model must be provided'//{introspection:location})
       if (.not.present(cosmologyFunctions_)) &
            & call Error_Report('cosmology functions must be provided for correlation model'//{introspection:location})
       self%correlationNormalization    =  correlationNormalization
       self%correlationMassExponent     =  correlationMassExponent
       self%correlationRedshiftExponent =  correlationRedshiftExponent
       self%cosmologyFunctions_         => cosmologyFunctions_
       !![
       <referenceCountIncrement owner="self" isResult="yes" object="cosmologyFunctions_"/>
       !!]
    else
       self%correlationNormalization    =  0.0d0
       self%correlationMassExponent     =  0.0d0
       self%correlationRedshiftExponent =  0.0d0
       self%cosmologyFunctions_         => null()
    end if
    return
  end function nbodyHaloMassErrorTrenti2010Internal

  subroutine trenti2010Destructor(self)
    !!{
    Destructor for the \refClass{nbodyHaloMassErrorTrenti2010} N-body statistics class.
    !!}
    implicit none
    type(nbodyHaloMassErrorTrenti2010), intent(inout) :: self

    if (associated(self%cosmologyFunctions_)) then
       !![
       <objectDestructor name="self%cosmologyFunctions_"/>
       !!]
    end if
    return
  end subroutine trenti2010Destructor

  double precision function trenti2010Correlation(self,node1,node2)
    !!{
    Return the correlation of the masses of a pair of N-body halos.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (nbodyHaloMassErrorTrenti2010), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node1    , node2
    class           (nodeComponentBasic          ), pointer       :: basic1   , basic2
    double precision                                              :: massRatio, expansionFactorRatio

    ! Extract mass and expansion factor ratios.
    basic1               =>                                           node1 %basic()
    basic2               =>                                           node2 %basic()
    massRatio            =  +                                         basic2%mass ()  &
         &                  /                                         basic1%mass ()
    expansionFactorRatio =  +self%cosmologyFunctions_%expansionFactor(basic2%time ()) &
         &                  /self%cosmologyFunctions_%expansionFactor(basic1%time ())
    ! Ensure ratios are below unity, invert otherwise.
    if (           massRatio > 1.0d0)            massRatio=1.0d0/           massRatio
    if (expansionFactorRatio > 1.0d0) expansionFactorRatio=1.0d0/expansionFactorRatio
    ! Evaluate the correlation.
    trenti2010Correlation=+self%correlationNormalization                                   &
         &                *                    massRatio**self%correlationMassExponent     &
         &                *         expansionFactorRatio**self%correlationRedshiftExponent
    return
  end function trenti2010Correlation
