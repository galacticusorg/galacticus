!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Contains a module which implements a cold dark matter particle class.

  use Cosmology_Parameters

  !# <darkMatterParticle name="darkMatterParticleCDM">
  !#  <description>Provides a cold dark matter particle.</description>
  !# </darkMatterParticle>
  type, extends(darkMatterParticleClass) :: darkMatterParticleCDM
     !% A cold dark matter particle class.
     private
   contains
  end type darkMatterParticleCDM
  
  interface darkMatterParticleCDM
     !% Constructors for the ``{\normalfont \ttfamily CDM}'' dark matter particle class.
     module procedure CDMConstructorParameters
  end interface darkMatterParticleCDM

contains

  function CDMConstructorParameters(parameters)
    !% Constructor for the ``{\normalfont \ttfamily CDM}'' dark matter particle class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(darkMatterParticleCDM)                :: CDMConstructorParameters
    type(inputParameters      ), intent(inout) :: parameters
    !GCC$ attributes unused :: parameters

    CDMConstructorParameters=darkMatterParticleCDM()
    return
  end function CDMConstructorParameters

  subroutine CDMDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(darkMatterParticleCDM), intent(inout) :: self
    type (inputParameters      ), intent(inout) :: descriptor
    !GCC$ attributes unused :: self
    
    call descriptor%addParameter("darkMatterParticleMethod","CDM")
    return
  end subroutine CDMDescriptor
