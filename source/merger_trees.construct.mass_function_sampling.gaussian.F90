!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implementation of a merger tree halo mass function sampling class in which the sampling rate is given by a Gaussian distribution in halo mass.

  !# <mergerTreeHaloMassFunctionSampling name="mergerTreeHaloMassFunctionSamplingGaussian" defaultThreadPrivate="yes">
  !#  <description>A merger tree halo mass function sampling class in which the sampling rate is given by a Gaussian distribution in halo mass.</description>
  !# </mergerTreeHaloMassFunctionSampling>
  type, extends(mergerTreeHaloMassFunctionSamplingClass) :: mergerTreeHaloMassFunctionSamplingGaussian
     !% Implementation of merger tree halo mass function sampling class in which the sampling rate is given by a Gaussian distribution in halo mass.
     private
     double precision :: mean, rootVariance
   contains
     procedure :: sample => gaussianSample
  end type mergerTreeHaloMassFunctionSamplingGaussian

  interface mergerTreeHaloMassFunctionSamplingGaussian
     !% Constructors for the {\normalfont \ttfamily gaussian} merger tree halo mass function sampling class.
     module procedure gaussianConstructorParameters
     module procedure gaussianConstructorInternal
  end interface mergerTreeHaloMassFunctionSamplingGaussian

contains

  function gaussianConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily gaussian} merger tree halo mass function sampling class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (mergerTreeHaloMassFunctionSamplingGaussian)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    double precision                                                            :: mean      , rootVariance
    
    !# <inputParameter>
    !#   <name>mean</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The mean mass of halo to simulate when using a Gaussian sampling of the halo mass function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>rootVariance</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The dispersion in mass of halo to simulate when using a Gaussian sampling of the halo mass function.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    self=mergerTreeHaloMassFunctionSamplingGaussian(mean,rootVariance)
    !# <inputParametersValidate source="parameters"/>
    return
  end function gaussianConstructorParameters

  function gaussianConstructorInternal(mean,rootVariance) result(self)
    !% Internal constructor for the {\normalfont \ttfamily gaussian} merger tree halo mass function sampling class.
    implicit none
    type            (mergerTreeHaloMassFunctionSamplingGaussian)                :: self
    double precision                                            , intent(in   ) :: mean, rootVariance
    !# <constructorAssign variables="mean, rootVariance"/>
    
    return
  end function gaussianConstructorInternal

  double precision function gaussianSample(self,mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate using a volume-limited sampling.
    implicit none
    class           (mergerTreeHaloMassFunctionSamplingGaussian), intent(inout) :: self
    double precision                                            , intent(in   ) :: mass       , massMaximum, &
         &                                                                         massMinimum, time
    !GCC$ attributes unused :: time

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
