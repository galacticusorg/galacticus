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
Implements an N-body dark matter halo mass error class using the model of \cite{trenti_how_2010}.
!!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <nbodyHaloMassError name="nbodyHaloMassErrorTrenti2010">
   <description>An N-body dark matter halo mass error class that models correlated halo mass errors between different halos and epochs using the fitting formula of \cite{trenti_how_2010}. The simulation particle mass and correlation model parameters $C_0$, $\alpha$, and $\beta$ are set via the corresponding input parameters. This class is implemented as a specialization of the \refClass{nbodyHaloMassErrorPowerLaw} class, with the parent class parameters set to reproduce the \cite{trenti_how_2010} fractional error formula and to enable the power-law correlation model.</description>
  </nbodyHaloMassError>
  !!]
  type, extends(nbodyHaloMassErrorPowerLaw) :: nbodyHaloMassErrorTrenti2010
     !!{
     An N-body halo mass error class using the model of \cite{trenti_how_2010}.
     !!}
     private
  end type nbodyHaloMassErrorTrenti2010

  interface nbodyHaloMassErrorTrenti2010
     !!{
     Constructors for the \refClass{nbodyHaloMassErrorTrenti2010} N-body halo mass error class.
     !!}
     module procedure nbodyHaloMassErrorTrenti2010ConstructorParameters
     module procedure nbodyHaloMassErrorTrenti2010ConstructorInternal
  end interface nbodyHaloMassErrorTrenti2010

contains

  function nbodyHaloMassErrorTrenti2010ConstructorParameters(parameters) result(self)
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
  end function nbodyHaloMassErrorTrenti2010ConstructorParameters

  function nbodyHaloMassErrorTrenti2010ConstructorInternal(massParticle,correlationNormalization,correlationMassExponent,correlationRedshiftExponent,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{nbodyHaloMassErrorTrenti2010} N-body halo mass error class. \cite{trenti_how_2010} report
    a normalization of the fractional error in particle number of 0.15 at $N=1000$ particles. Since this is based on
    comparisons of halos in simulations differing in number of particles by a factor $8$ this actually overestimates the
    normalization by a factor $\sqrt{5/4}$. Therefore, we use a normalization of $0.135$ here.

    The \cite{trenti_how_2010} fractional error model, $\sigma(M) = N_0 / M^{1/3}$ with
    $N_0 = 0.135 \, (1000 \, m_\mathrm{p})^{1/3}$, is reproduced by the parent
    \refClass{nbodyHaloMassErrorPowerLaw} class by setting the power-law exponent to $-1/3$, the
    high-mass error to zero, and the normalization $\sigma_{12} = N_0 / M_\mathrm{ref}^{1/3}$ where
    $M_\mathrm{ref} = 10^{12} \mathrm{M}_\odot$ is the reference mass used by the parent class.
    !!}
    use :: Math_Exponentiation, only : cubeRoot
    implicit none
    type            (nbodyHaloMassErrorTrenti2010)                        :: self
    double precision                              , intent(in   )         :: massParticle
    double precision                              , intent(in   )         :: correlationNormalization             , correlationMassExponent, &
         &                                                                   correlationRedshiftExponent
    class           (cosmologyFunctionsClass     ), intent(in   ), target :: cosmologyFunctions_
    double precision                              , parameter             :: normalizationTrenti         =+0.135d0
    double precision                              , parameter             :: particleNumberReference     =+1.000d3
    double precision                                                      :: normalization

    ! Compute the power-law normalization, sigma_12, that reproduces the Trenti et al. (2010) fractional error formula when
    ! combined with an exponent of -1/3 and zero high-mass error.
    normalization=+         normalizationTrenti                  &
         &        *cubeRoot(                                     &
         &                  +particleNumberReference             &
         &                  *massParticle                        &
         &                  /                   massReference    &
         &                 )
    self%nbodyHaloMassErrorPowerLaw=nbodyHaloMassErrorPowerLaw(                                             &
         &                                                     normalization              =normalization              , &
         &                                                     exponent                   =-1.0d0/3.0d0               , &
         &                                                     fractionalErrorHighMass    =+0.0d0                     , &
         &                                                     correlationModelTrivial    =.false.                    , &
         &                                                     correlationNormalization   =correlationNormalization   , &
         &                                                     correlationMassExponent    =correlationMassExponent    , &
         &                                                     correlationRedshiftExponent=correlationRedshiftExponent, &
         &                                                     cosmologyFunctions_        =cosmologyFunctions_          &
         &                                                    )
    return
  end function nbodyHaloMassErrorTrenti2010ConstructorInternal
