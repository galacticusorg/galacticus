!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a nonlinear power spectrum class in which the nonlinear power spectrum is just the linear
!% power spectrum. Intended primarily for testing purposes.

  use Linear_Growth
  use Power_Spectra

  !# <powerSpectrumNonlinear name="powerSpectrumNonlinearLinear">
  !#  <description>Provides a nonlinear power spectrum class in which the power spectrum equals the linear theory power spectrum. Intended primarily for testing purposes.</description>
  !# </powerSpectrumNonlinear>
  type, extends(powerSpectrumNonlinearClass) :: powerSpectrumNonlinearLinear
     !% A linear transfer function class.
     private
     class(linearGrowthClass ), pointer :: linearGrowth_
     class(powerSpectrumClass), pointer :: powerSpectrum_
   contains
     final     ::          linearDestructor
     procedure :: value => linearValue
  end type powerSpectrumNonlinearLinear

  interface powerSpectrumNonlinearLinear
     !% Constructors for the linear nonlinear power spectrum class.
     module procedure linearConstructorParameters
     module procedure linearConstructorInternal
  end interface powerSpectrumNonlinearLinear

contains

  function linearConstructorParameters(parameters)
    !% Constructor for the linear nonlinear power spectrum class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(powerSpectrumNonlinearLinear)                :: linearConstructorParameters
    type(inputParameters             ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    

    !# <objectBuilder class="powerSpectrum" name="linearConstructorParameters%powerSpectrum_" source="parameters"/>
    !# <objectBuilder class="linearGrowth"  name="linearConstructorParameters%linearGrowth_"  source="parameters"/>
    return
  end function linearConstructorParameters

  function linearConstructorInternal(powerSpectrum_,linearGrowth_)
    !% Internal constructor for the linear nonlinear power spectrum class.
    implicit none
    type (powerSpectrumNonlinearLinear)                        :: linearConstructorInternal
    class(powerSpectrumClass          ), intent(in   ), target :: powerSpectrum_
    class(linearGrowthClass           ), intent(in   ), target :: linearGrowth_

    linearConstructorInternal%powerSpectrum_ => powerSpectrum_
    linearConstructorInternal%linearGrowth_  => linearGrowth_
    return
  end function linearConstructorInternal

  subroutine linearDestructor(self)
    !% Destructor for the linear nonlinear power spectrum class.
    implicit none
    type(powerSpectrumNonlinearLinear), intent(inout) :: self

    !# <objectDestructor name="self%powerSpectrum_"/>
    !# <objectDestructor name="self%linearGrowth_" />
    return
  end subroutine linearDestructor

  double precision function linearValue(self,wavenumber,time)
    !% Return the nonlinear power spectrum at the given wavenumber.
    implicit none
    class           (powerSpectrumNonlinearLinear), intent(inout) :: self
    double precision                              , intent(in   ) :: wavenumber, time

    linearValue=self%powerSpectrum_%power(wavenumber)*self%linearGrowth_%value(time)**2
    return
  end function linearValue
