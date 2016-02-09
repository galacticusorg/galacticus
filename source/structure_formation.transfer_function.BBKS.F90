!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
!!    Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Contains a module which implements the transfer function fitting function of
  !% \cite{bardeen_statistics_1986}.

  use Cosmology_Parameters

  !# <transferFunction name="transferFunctionBBKS">
  !#  <description>Provides the \cite{bardeen_statistics_1986} fitting function for the transfer function.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionBBKS
     !% A bbks transfer function class.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
     double precision                                    :: Gamma
   contains
     final     ::                          bbksDestructor
     procedure :: value                 => bbksValue
     procedure :: logarithmicDerivative => bbksLogarithmicDerivative
     procedure :: halfModeMass          => bbksHalfModeMass
     procedure :: descriptor            => bbksDescriptor
  end type transferFunctionBBKS

  interface transferFunctionBBKS
     !% Constructors for the ``BBKS'' transfer function class.
     module procedure bbksConstructorParameters
     module procedure bbksConstructorInternal
  end interface transferFunctionBBKS

  ! Fitting function parameters.
  double precision, parameter :: bbksP=2.34d0, bbksA=3.89d0, bbksB=16.10d0, bbksC=5.46d0, bbksD=6.71d0
  
contains

  function bbksConstructorParameters(parameters)
    !% Constructor for the ``BBKS'' transfer function class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type (transferFunctionBBKS    )                :: bbksConstructorParameters
    type (inputParameters         ), intent(in   ) :: parameters
    class(cosmologyParametersClass), pointer       :: cosmologyParameters_
    
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    bbksConstructorParameters=bbksConstructorInternal(cosmologyParameters_)
    return
  end function bbksConstructorParameters

  function bbksConstructorInternal(cosmologyParameters_)
    !% Internal constructor for the ``BBKS'' transfer function class.
    implicit none
    type (transferFunctionBBKS    )                                  :: bbksConstructorInternal
    class(cosmologyParametersClass), intent(in   ), target, optional :: cosmologyParameters_    

    ! Determine the cosmological parameters to use.
    if (present(cosmologyParameters_)) then
       bbksConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    else
       bbksConstructorInternal%cosmologyParameters_ => cosmologyParameters()
    end if
    ! Compute the Gamma parameter.
    bbksConstructorInternal%Gamma=+             bbksConstructorInternal%cosmologyParameters_%OmegaMatter   (                  ) &
         &                        *             bbksConstructorInternal%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH) &
         &                        *exp(                                                                                         &
         &                             -        bbksConstructorInternal%cosmologyParameters_%OmegaBaryon   (                  ) &
         &                             *(                                                                                       &
         &                               +1.0d0                                                                                 &
         &                               +sqrt(                                                                                 &
         &                                     +2.0d0                                                                           &
         &                                     *bbksConstructorInternal%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH) &
         &                                    )                                                                                 &
         &                               /      bbksConstructorInternal%cosmologyParameters_%OmegaMatter   (                  ) &
         &                              )                                                                                       &
         &                            )                                                                                         &
         &                        /(                                                                                            &
         &                          +           bbksConstructorInternal%cosmologyParameters_%temperatureCMB(                  ) &
         &                          /2.7d0                                                                                      &
         &                         )**2
    return
  end function bbksConstructorInternal

  subroutine bbksDestructor(self)
    !% Destructor for the ``BBKS'' transfer function class.
    implicit none
    type(transferFunctionBBKS), intent(inout) :: self
    
    if     (                                                       &
         &   associated(self%cosmologyParameters_                ) &
         &  .and.                                                  &
         &              self%cosmologyParameters_%isFinalizable()  &
         & ) deallocate(self%cosmologyParameters_                )
    return
  end subroutine bbksDestructor

  double precision function bbksValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBBKS), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber
    double precision                                      :: wavenumberHUnits, q

    ! Get wavenumber in "little-h" units.
    wavenumberHUnits=+wavenumber                                                   &
         &           /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    q               =+wavenumberHUnits &
         &           /self%Gamma
    bbksValue       =+(              &
         &             +log(         &
         &                  +1.00d0  &
         &                  +bbksP   &
         &                  *q       &
         &                 )         &
         &             /bbksP        &
         &             /q            &
         &            )              &
         &           /(              &
         &             +1.0d0        &
         &             +(bbksA*q)    &
         &             +(bbksB*q)**2 &
         &             +(bbksC*q)**3 &
         &             +(bbksD*q)**4 &
         &            )**0.25d0
    return
  end function bbksValue

  double precision function bbksLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBBKS), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber
    double precision                                      :: wavenumberHUnits, q

    ! Get wavenumber in "little-h" units.
    wavenumberHUnits=+wavenumber                                                   &
         &           /self%cosmologyParameters_%HubbleConstant(hubbleUnitsLittleH)
    q               =+wavenumberHUnits &
         &           /self%Gamma
    bbksLogarithmicDerivative=-2.0d0                &
         &                    +(                    &
         &                      +4.0d0              &
         &                      +3.0d0*(bbksA*q)    &
         &                      +2.0d0*(bbksB*q)**2 &
         &                      +      (bbksC*q)**3 &
         &                    )                     &
         &                    /(                    &
         &                      +4.0d0              &
         &                      *(                  &
         &                        +1.0d0            &
         &                        +      q          &
         &                        *(                &
         &                          +    bbksA      &
         &                          +    q          &
         &                          *(              &
         &                            +  bbksB**2   &
         &                            +  q          &
         &                            *(            &
         &                              +bbksC**3   &
         &                              +bbksD**4   &
         &                              *q          &
         &                             )            &
         &                           )              &
         &                         )                &
         &                       )                  &
         &                     )                    &
         &                    +(                    &
         &                      +bbksP              &
         &                      *q                  &
         &                     )                    &
         &                    /(                    &
         &                      +(                  &
         &                        +1.0d0            &
         &                        +bbksP            &
         &                        *q                &
         &                       )                  &
         &                      *log(               &
         &                           +1.0d0         &
         &                           +bbksP         &
         &                           *q             &
         &                          )               &
         &                     )
    return
  end function bbksLogarithmicDerivative

  double precision function bbksHalfModeMass(self)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function. Not supported in this implementation.
    use Galacticus_Error
    implicit none
    class(transferFunctionBBKS), intent(inout) :: self

    call Galacticus_Error_Report('bbksHalfModeMass','not supported by this implementation')
    return
  end function bbksHalfModeMass

  subroutine bbksDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (transferFunctionBBKS), intent(inout) :: self
    type     (inputParameters     ), intent(inout) :: descriptor
    type     (node                ), pointer       :: parameterNode
    type     (inputParameters     )                :: subParameters

    call descriptor%addParameter("transferFunctionMethod","BBKS")
    parameterNode => descriptor%node("transferFunctionMethod")
    subParameters=inputParameters(parameterNode)
    call self%cosmologyParameters_%descriptor(subParameters)
    return
  end subroutine bbksDescriptor
