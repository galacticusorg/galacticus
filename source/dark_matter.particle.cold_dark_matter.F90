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
Implements a cold dark matter particle class.
!!}


  !![
  <darkMatterParticle name="darkMatterParticleCDM">
   <description>Provides a cold dark matter particle.</description>
  </darkMatterParticle>
  !!]
  type, extends(darkMatterParticleClass) :: darkMatterParticleCDM
     !!{
     A cold dark matter particle class.
     !!}
     private
   contains
     procedure :: mass => CDMMass
  end type darkMatterParticleCDM

  interface darkMatterParticleCDM
     !!{
     Constructors for the ``{\normalfont \ttfamily CDM}'' dark matter particle class.
     !!}
     module procedure CDMConstructorParameters
  end interface darkMatterParticleCDM

contains

  function CDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``{\normalfont \ttfamily CDM}'' dark matter particle class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(darkMatterParticleCDM)                :: self
    type(inputParameters      ), intent(inout) :: parameters

    self=darkMatterParticleCDM()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function CDMConstructorParameters

  double precision function CDMMass(self)
    !!{
    Return the mass, in units of keV, of a cold dark matter particle. An infinite mass is assumed.
    !!}
    implicit none
    class(darkMatterParticleCDM), intent(inout) :: self

    CDMMass=huge(0.0d0)
    return
  end function CDMMass
