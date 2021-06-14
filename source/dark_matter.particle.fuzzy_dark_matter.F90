!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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

!% Contains a module which implements a fuzzy dark matter particle class.

  !# <darkMatterParticle name="darkMatterParticleFuzzyDarkMatter">
  !#  <description>Provides a fuzzy dark matter particle.</description>
  !# </darkMatterParticle>
  type, extends(darkMatterParticleClass) :: darkMatterParticleFuzzyDarkMatter
     !% A fuzzy dark matter particle class.
     private
     double precision :: massValue
   contains
     procedure :: mass => fuzzyDMMass
  end type darkMatterParticleFuzzyDarkMatter

  interface darkMatterParticleFuzzyDarkMatter
     !% Constructors for the ``{\normalfont \ttfamily fuzzyDarkMatter}'' dark matter particle class.
     module procedure fuzzyDMConstructorParameters
     module procedure fuzzyDMConstructorInternal
  end interface darkMatterParticleFuzzyDarkMatter

contains

  function fuzzyDMConstructorParameters(parameters) result(self)
    !% Constructor for the ``{\normalfont \ttfamily fuzzyDarkMatter}'' dark matter particle class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterParticleFuzzyDarkMatter)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: massValue

    !# <inputParameter>
    !#   <name>mass</name>
    !#   <source>parameters</source>
    !#   <variable>massValue</variable>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The mass (in units of $10^{-22}$~eV) of the fuzzy dark matter particle.</description>
    !# </inputParameter>
    self=darkMatterParticleFuzzyDarkMatter(massValue)
    !# <inputParametersValidate source="parameters"/>
    return
  end function fuzzyDMConstructorParameters

  function fuzzyDMConstructorInternal(mass) result(self)
    !% Internal constructor for the ``{\normalfont \ttfamily fuzzyDarkMatter}'' dark matter particle class.
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    type            (darkMatterParticleFuzzyDarkMatter)                :: self
    double precision                                   , intent(in   ) :: mass
    double precision                                   , parameter     :: massUnitsFuzzyDM=1.0d-22/kilo

    self%massValue=+mass             &
         &         /massUnitsFuzzyDM
    return
  end function fuzzyDMConstructorInternal

  double precision function fuzzyDMMass(self)
    !% Return the mass, in units of keV, of a fuzzy dark matter particle.
    implicit none
    class(darkMatterParticleFuzzyDarkMatter), intent(inout) :: self

    fuzzyDMMass=self%massValue
    return
  end function fuzzyDMMass
