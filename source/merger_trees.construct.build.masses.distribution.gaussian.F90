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
  Implementation of a merger tree halo mass function sampling class in which the sampling rate is given by a Gaussian distribution in halo mass.
  !!}

  !![
  <mergerTreeBuildMassDistribution name="mergerTreeBuildMassDistributionGaussian">
   <description>A merger tree halo mass function sampling class in which the sampling rate is given by a Gaussian distribution in halo mass.</description>
  </mergerTreeBuildMassDistribution>
  !!]
  type, extends(mergerTreeBuildMassDistributionClass) :: mergerTreeBuildMassDistributionGaussian
     !!{
     Implementation of merger tree halo mass function sampling class in which the sampling rate is given by a Gaussian distribution in halo mass.
     !!}
     private
     double precision :: mean, rootVariance
   contains
     procedure :: sample => gaussianSample
  end type mergerTreeBuildMassDistributionGaussian

  interface mergerTreeBuildMassDistributionGaussian
     !!{
     Constructors for the \refClass{mergerTreeBuildMassDistributionGaussian} merger tree halo mass function sampling class.
     !!}
     module procedure gaussianConstructorParameters
     module procedure gaussianConstructorInternal
  end interface mergerTreeBuildMassDistributionGaussian

contains

  function gaussianConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeBuildMassDistributionGaussian} merger tree halo mass function sampling class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeBuildMassDistributionGaussian)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    double precision                                                         :: mean      , rootVariance

    !![
    <inputParameter>
      <name>mean</name>
      <description>The mean mass of halo to simulate when using a Gaussian sampling of the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>rootVariance</name>
      <description>The dispersion in mass of halo to simulate when using a Gaussian sampling of the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=mergerTreeBuildMassDistributionGaussian(mean,rootVariance)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function gaussianConstructorParameters

  function gaussianConstructorInternal(mean,rootVariance) result(self)
    !!{
    Internal constructor for the \refClass{mergerTreeBuildMassDistributionGaussian} merger tree halo mass function sampling class.
    !!}
    implicit none
    type            (mergerTreeBuildMassDistributionGaussian)                :: self
    double precision                                         , intent(in   ) :: mean, rootVariance
    !![
    <constructorAssign variables="mean, rootVariance"/>
    !!]

    return
  end function gaussianConstructorInternal

  double precision function gaussianSample(self,mass,time,massMinimum,massMaximum)
    !!{
    Computes the halo mass function sampling rate using a volume-limited sampling.
    !!}
    implicit none
    class           (mergerTreeBuildMassDistributionGaussian), intent(inout) :: self
    double precision                                         , intent(in   ) :: mass       , massMaximum, &
         &                                                                      massMinimum, time
    !$GLC attributes unused :: time

    if (mass <= massMinimum .or. mass > massMaximum) then
       gaussianSample=0.0d0
    else
       gaussianSample=+exp(                     &
            &              -0.5d0               &
            &              *(                   &
            &                +(                 &
            &                  +     mass       &
            &                  -self%mean       &
            &                 )                 &
            &                /self%rootVariance &
            &               )**2                &
            &             )
    end if
    return
  end function gaussianSample
