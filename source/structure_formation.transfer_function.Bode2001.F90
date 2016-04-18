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

!% Contains a module which implements a transfer function class based on the \gls{wdm} modifier of \cite{bode_halo_2001}.

  use Cosmology_Parameters

  !# <transferFunction name="transferFunctionBode2001">
  !#  <description>Provides a transfer function based on the \gls{wdm} modifier of \cite{bode_halo_2001}.</description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionBode2001
     !% A transfer function class which modifies another transfer function using the \gls{wdm} modifier of \cite{bode_halo_2001}.
     private
     double precision                                    :: scaleCutOff         , epsilon, &
          &                                                 eta                 , nu
     class           (transferFunctionClass   ), pointer :: transferFunctionCDM
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
   contains
     final     ::                          bode2001Destructor
     procedure :: value                 => bode2001Value
     procedure :: logarithmicDerivative => bode2001LogarithmicDerivative
     procedure :: halfModeMass          => bode2001HalfModeMass
     procedure :: descriptor            => bode2001Descriptor
  end type transferFunctionBode2001

  interface transferFunctionBode2001
     !% Constructors for the ``{\normalfont \ttfamily bode2001}'' transfer function class.
     module procedure bode2001ConstructorParameters
     module procedure bode2001ConstructorInternal
  end interface transferFunctionBode2001

contains

  function bode2001ConstructorParameters(parameters)
    !% Constructor for the ``{\normalfont \ttfamily bode2001}'' transfer function class which takes a parameter set as input.
    use Input_Parameters2
    use Galacticus_Error
    implicit none
    type            (transferFunctionBode2001)                :: bode2001ConstructorParameters
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_    
    double precision                                          :: scaleCutOff                  , epsilon, &
         &                                                       eta                          , nu
    !# <inputParameterList label="allowedParameterNames" />

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunctionMethod')) call Galacticus_Error_Report("bode2001ConstructorParameters","an explicit 'transferFunctionMethod' must be given"//{introspection:location})
    ! Read parameters.
    !# <inputParameter>
    !#   <name>scaleCutOff</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The cut-off scale in the transfer function due to warm dark matter.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>epsilon</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.361d0</defaultValue>
    !#   <defaultSource>\citep{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta</name>
    !#   <source>parameters</source>
    !#   <defaultValue>5.0d0</defaultValue>
    !#   <defaultSource>\citep{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>nu</name>
    !#   <source>parameters</source>
    !#   <defaultValue>1.2d0</defaultValue>
    !#   <defaultSource>\citep{barkana_constraints_2001}</defaultSource>
    !#   <description>The parameter $\epsilon$ appearing in the warm dark matter transfer function \citep{barkana_constraints_2001}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    ! Call the internal constructor
    bode2001ConstructorParameters=bode2001ConstructorInternal(transferFunctionCDM,scaleCutOff,epsilon,eta,nu,cosmologyParameters_)
    return
  end function bode2001ConstructorParameters

  function bode2001ConstructorInternal(transferFunctionCDM,scaleCutOff,epsilon,eta,nu,cosmologyParameters_)
    !% Internal constructor for the ``{\normalfont \ttfamily bode2001}'' transfer function class.
    implicit none
    type            (transferFunctionBode2001)                                  :: bode2001ConstructorInternal
    class           (transferFunctionClass   ), target, intent(in   )           :: transferFunctionCDM
    double precision                                  , intent(in   )           :: scaleCutOff                , epsilon, &
         &                                                                         eta                        , nu
    class           (cosmologyParametersClass), target, intent(in   ), optional :: cosmologyParameters_    

    bode2001ConstructorInternal%transferFunctionCDM => transferFunctionCDM
    bode2001ConstructorInternal%scaleCutOff         =  scaleCutOff
    bode2001ConstructorInternal%epsilon             =  epsilon
    bode2001ConstructorInternal%eta                 =  eta
    bode2001ConstructorInternal%nu                  =  nu
    ! Determine the cosmological parameters to use.
    if (present(cosmologyParameters_)) then
       bode2001ConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    else
       bode2001ConstructorInternal%cosmologyParameters_ => cosmologyParameters()
    end if
    return
  end function bode2001ConstructorInternal

  subroutine bode2001Destructor(self)
    !% Destructor for the bode2001 transfer function class.
    implicit none
    type(transferFunctionBode2001), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%transferFunctionCDM" />
    return
  end subroutine bode2001Destructor

  double precision function bode2001Value(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    bode2001Value       =+self%transferFunctionCDM%value(wavenumber)
    if (self%scaleCutOff > 0.0d0)                 &
         & bode2001Value=+bode2001Value           &
         &               /(                       &
         &                  1.0d0                 &
         &                 +                      &
         &                  (                     &
         &                   +self%epsilon        &
         &                   *wavenumber          &
         &                   *self%scaleCutOff    &
         &                  )**(2.0d0   *self%nu) &
         &                )  **(self%eta/self%nu)
    return
  end function bode2001Value

  double precision function bode2001LogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , intent(in   ) :: wavenumber

    bode2001LogarithmicDerivative=+self%transferFunctionCDM%logarithmicDerivative(wavenumber)    
    if (self%scaleCutOff > 0.0d0)                                       &
         & bode2001LogarithmicDerivative=+bode2001LogarithmicDerivative &
         &                               +2.0d0                         &
         &                               *self%eta                      &
         &                               *(                             &
         &                                 +1.0d0                       &
         &                                 /(                           &
         &                                   +1.0d0                     &
         &                                   +(                         &
         &                                     +self%epsilon            &
         &                                     *wavenumber              &
         &                                     *self%scaleCutOff        &
         &                                    )**(2.0d0*self%nu)        &
         &                                  )                           &
         &                                 -1.0d0                       &
         &                                )
    return
  end function bode2001LogarithmicDerivative
  
  double precision function bode2001HalfModeMass(self)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is suppressed by a factor of two relative
    !% to a \gls{cdm} transfer function.
    use Numerical_Constants_Math
    implicit none
    class           (transferFunctionBode2001), intent(inout) :: self
    double precision                          , parameter     :: wavenumberHalfModeScaleFree=sqrt(0.25d0+2.0d0*log(2.0d0))-0.5d0
    double precision                                          :: matterDensity                                                  , wavenumberHalfMode
    
    matterDensity       =+self%cosmologyParameters_%OmegaMatter    () &
         &               *self%cosmologyParameters_%densityCritical()
    wavenumberHalfMode  =+(                            &
         &                 +2.0d0**(+self%nu/self%eta) &
         &                 -1.0d0                      &
         &                )      **(+0.5d0  /self%nu ) &
         &                /self%epsilon                &
         &                /self%scaleCutOff
    bode2001HalfModeMass=+4.0d0                &
         &               *Pi                   &
         &               /3.0d0                &
         &               *matterDensity        &
         &               *(                    &
         &                 +Pi                 &
         &                 /wavenumberHalfMode &
         &               )**3
    return
  end function bode2001HalfModeMass

  subroutine bode2001Descriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class    (transferFunctionBode2001), intent(inout) :: self
    type     (inputParameters         ), intent(inout) :: descriptor
    type     (inputParameters         )                :: subParameters
    character(len=10                  )                :: parameterLabel

    call descriptor%addParameter("transferFunctionMethod","bode2001")
    subParameters=descriptor%subparameters("transferFunctionMethod")
    write (parameterLabel,'(f10.6)') self%scaleCutOff
    call subParameters%addParameter("scaleCutOff",trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%epsilon
    call subParameters%addParameter("epsilon"    ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%eta
    call subParameters%addParameter("eta"        ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(f10.6)') self%nu
    call subParameters%addParameter("nu"         ,trim(adjustl(parameterLabel)))
    call self%transferFunctionCDM% descriptor(subParameters)
    call self%cosmologyParameters_%descriptor(subParameters)
    return
  end subroutine bode2001Descriptor
