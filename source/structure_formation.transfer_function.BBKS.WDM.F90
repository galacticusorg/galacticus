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

!% Contains a module which implements a transfer function class based on the \gls{wdm} modifier of \cite{bardeen_statistics_1986}.

  use Cosmology_Parameters
  use Dark_Matter_Particles

  !# <transferFunction name="transferFunctionBBKSWDM">
  !#  <description>Provides a transfer function based on the \gls{wdm} modifier of \cite{bardeen_statistics_1986}.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionBBKSWDM
     !% A transfer function class which modifies another transfer function using the \gls{wdm} modifier of \cite{bardeen_statistics_1986}.
     private
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
     class           (darkMatterParticleClass ), pointer :: darkMatterParticle_
     double precision                                    :: lengthFreeStreaming
   contains
     final     ::                          bbksWDMDestructor
     procedure :: value                 => bbksWDMValue
     procedure :: logarithmicDerivative => bbksWDMLogarithmicDerivative
     procedure :: halfModeMass          => bbksWDMHalfModeMass
     procedure :: descriptor            => bbksWDMDescriptor
  end type transferFunctionBBKSWDM

  interface transferFunctionBBKSWDM
     !% Constructors for the ``{\normalfont \ttfamily bbksWDM}'' transfer function class.
     module procedure bbksWDMConstructorParameters
     module procedure bbksWDMConstructorInternal
  end interface transferFunctionBBKSWDM

contains

  function bbksWDMConstructorParameters(parameters)
    !% Constructor for the ``{\normalfont \ttfamily bbksWDM}'' transfer function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (transferFunctionBBKSWDM )                :: bbksWDMConstructorParameters
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_    
    class           (darkMatterParticleClass ), pointer       :: darkMatterParticle_    
    !# <inputParameterList label="allowedParameterNames" />
    
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !# <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    ! Call the internal constructor
    bbksWDMConstructorParameters =bbksWDMConstructorInternal(transferFunctionCDM,cosmologyParameters_,darkMatterParticle_)
    return
  end function bbksWDMConstructorParameters

  function bbksWDMConstructorInternal(transferFunctionCDM,cosmologyParameters_,darkMatterParticle_)
    !% Internal constructor for the ``{\normalfont \ttfamily bbksWDM}'' transfer function class.
    use Galacticus_Error
    implicit none
    type            (transferFunctionBBKSWDM )                                  :: bbksWDMConstructorInternal
    class           (transferFunctionClass   ), target, intent(in   )           :: transferFunctionCDM
    class           (cosmologyParametersClass), target, intent(in   ), optional :: cosmologyParameters_    
    class           (darkMatterParticleClass ), target, intent(in   ), optional :: darkMatterParticle_    
    double precision                                                            :: degreesOfFreedomEffectiveDecoupling

    bbksWDMConstructorInternal%transferFunctionCDM => transferFunctionCDM
    ! Determine the cosmological parameters to use.
    if (present(cosmologyParameters_)) then
       bbksWDMConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    else
       bbksWDMConstructorInternal%cosmologyParameters_ => cosmologyParameters()
    end if
    ! Determine the dark matter particle to use.
    if (present(darkMatterParticle_)) then
       bbksWDMConstructorInternal%darkMatterParticle_ => darkMatterParticle_
    else
       bbksWDMConstructorInternal%darkMatterParticle_ => darkMatterParticle()
    end if
    ! Get degrees of freedom at the time at which the dark matter particle decoupled.
    select type (particle => bbksWDMConstructorInternal%darkMatterParticle_)
    class is (darkMatterParticleWDMThermal)
       degreesOfFreedomEffectiveDecoupling=particle%degreesOfFreedomEffectiveDecoupling()
    class default
       degreesOfFreedomEffectiveDecoupling=0.0d0
       call Galacticus_Error_Report('bbksWDMConstructorInternal','transfer function expects a thermal warm dark matter particle')
    end select
    ! Compute the free-streaming length-like parameter (equation G6 of BBKS).
    bbksWDMConstructorInternal%lengthFreeStreaming=+0.2d0                                                                                     &
         &                                         /(                                                                                         &
         &                                           +degreesOfFreedomEffectiveDecoupling                                                     &
         &                                           /100.d0                                                                                  &
         &                                          )**(4.0d0/3.0d0)                                                                          &
         &                                         /(                                                                                         &
         &                                           +(                                                                                       &
         &                                             +bbksWDMConstructorInternal%cosmologyParameters_%OmegaMatter   (                  )    &
         &                                             -bbksWDMConstructorInternal%cosmologyParameters_%OmegaBaryon   (                  )    &
         &                                            )                                                                                       &
         &                                           /  bbksWDMConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)**2 &
         &                                          )
    return
  end function bbksWDMConstructorInternal
  
  subroutine bbksWDMDestructor(self)
    !% Destructor for the bbksWDM transfer function class.
    implicit none
    type(transferFunctionBBKSWDM), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%darkMatterParticle_" />
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine bbksWDMDestructor

  double precision function bbksWDMValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBBKSWDM), intent(inout) :: self
    double precision                         , intent(in   ) :: wavenumber
    double precision                                         :: wavenumberScaleFree

    wavenumberScaleFree=+wavenumber               &
         &              *self%lengthFreeStreaming
    bbksWDMValue       =+self%transferFunctionCDM%value(wavenumber) &
         &              *exp(                                       &
         &                   -0.5d0                                 &
         &                   *(                                     &
         &                     +  wavenumberScaleFree               &
         &                     *(                                   &
         &                       +1.0d0                             &
         &                       +wavenumberScaleFree               &
         &                      )                                   &
         &                    )                                     &
         &                  )
    return
  end function bbksWDMValue

  double precision function bbksWDMLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBBKSWDM), intent(inout) :: self
    double precision                         , intent(in   ) :: wavenumber
    double precision                                         :: wavenumberScaleFree
    
    wavenumberScaleFree=+wavenumber               &
         &              *self%lengthFreeStreaming
    bbksWDMLogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber) &
         &                       -  wavenumberScaleFree                                      &
         &                       *(                                                          &
         &                         +0.5d0                                                    &
         &                         +wavenumberScaleFree                                      &
         &                        )    
    return
  end function bbksWDMLogarithmicDerivative
  
  double precision function bbksWDMHalfModeMass(self)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use Numerical_Constants_Math
    implicit none
    class           (transferFunctionBBKSWDM), intent(inout) :: self
    double precision                         , parameter     :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                         :: matterDensity                                                  , wavenumberHalfMode
    
    wavenumberHalfMode =+wavenumberHalfModeScaleFree &
         &              /self%lengthFreeStreaming
    matterDensity      =+self%cosmologyParameters_%OmegaMatter    () &
         &              *self%cosmologyParameters_%densityCritical()
    bbksWDMHalfModeMass=+4.0d0                &
         &              *Pi                   &
         &              /3.0d0                &
         &              *matterDensity        &
         &              *(                    &
         &                +Pi                 &
         &                /wavenumberHalfMode &
         &              )**3
    return
  end function bbksWDMHalfModeMass

  subroutine bbksWDMDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (transferFunctionBBKSWDM), intent(inout) :: self
    type     (inputParameters        ), intent(inout) :: descriptor
    type     (inputParameters        )                :: subParameters

    call descriptor%addParameter("transferFunctionMethod","BBKSWDM")
    subParameters=descriptor%subparameters("transferFunctionMethod")
    call self%transferFunctionCDM %descriptor(subParameters)
    call self%cosmologyParameters_%descriptor(subParameters)
    call self%darkMatterParticle_ %descriptor(subParameters)
   return
  end subroutine bbksWDMDescriptor
