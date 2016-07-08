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

!% Contains a module which implements an N-body dark matter halo mass error class in which
!% errors are a power-law in halo mass.

  !# <nbodyHaloMassError name="nbodyHaloMassErrorPowerLaw">
  !#  <description>An N-body dark matter halo mass error class in which errors are a power-law in halo mass.</description>
  !# </nbodyHaloMassError>
  type, extends(nbodyHaloMassErrorClass) :: nbodyHaloMassErrorPowerLaw
     !% An N-body halo mass error class in which errors are a power-law in halo mass.
     private
     double precision :: normalization, exponent
   contains
     procedure :: errorFractional => powerLawErrorFractional
  end type nbodyHaloMassErrorPowerLaw

  interface nbodyHaloMassErrorPowerLaw
     !% Constructors for the {\normalfont \ttfamily powerLaw} N-body halo mass error class.
     module procedure nbodyHaloMassErrorPowerLawParameters
     module procedure nbodyHaloMassErrorPowerLawInternal
  end interface nbodyHaloMassErrorPowerLaw

contains

  function nbodyHaloMassErrorPowerLawParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily powerLaw} N-body halo mass error class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(nbodyHaloMassErrorPowerLaw)                :: nbodyHaloMassErrorPowerLawParameters
    type(inputParameters           ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>normalization</name>
    !#   <source>parameters</source>
    !#   <variable>nbodyHaloMassErrorPowerLawParameters%normalization</variable>
    !#   <description>Parameter $\sigma_{12}$ appearing in model for random errors in the halo mass function.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>exponent</name>
    !#   <source>parameters</source>
    !#   <variable>nbodyHaloMassErrorPowerLawParameters%exponent</variable>
    !#   <description>
    !#    Parameter $\gamma$ appearing in model for random errors in the halo mass
    !#    function. Specifically, the fractional error is given by
    !#    \begin{equation}
    !#    \sigma(M) = \sigma_{12} \left({M_{\rm halo} \over 10^{12}M_\odot}\right)^\gamma,
    !#    \end{equation}
    !#    where $\sigma_{12}=${\normalfont \ttfamily [normalization]}, and $\gamma=${\normalfont
    !#    \ttfamily [exponent]}.
    !#   </description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    return
  end function nbodyHaloMassErrorPowerLawParameters

  function nbodyHaloMassErrorPowerLawInternal(normalization,exponent)
    !% Internal constructor for the {\normalfont \ttfamily powerLaw} N-body halo mass error class.
    implicit none
    type            (nbodyHaloMassErrorPowerLaw)                :: nbodyHaloMassErrorPowerLawInternal
    double precision                            , intent(in   ) :: normalization                     , exponent

    nbodyHaloMassErrorPowerLawInternal%normalization=normalization
    nbodyHaloMassErrorPowerLawInternal%exponent     =exponent
    return
  end function nbodyHaloMassErrorPowerLawInternal

  double precision function powerLawErrorFractional(self,node)
    !% Return the fractional error on the mass of an N-body halo in the power-law error model.
    implicit none
    class           (nbodyHaloMassErrorPowerLaw), intent(inout)            :: self
    type            (treeNode                  ), intent(inout), pointer   :: node
    double precision                                           , parameter :: massNormalization=1.0d12
    class           (nodeComponentBasic        )               , pointer   :: basic

    basic                   =>  node%basic        ()
    powerLawErrorFractional =  +self%normalization   &
         &                     *(                    &
         &                       +basic%mass()       &
         &                       /massNormalization  &
         &                      )**self%exponent
    return
  end function powerLawErrorFractional
  
