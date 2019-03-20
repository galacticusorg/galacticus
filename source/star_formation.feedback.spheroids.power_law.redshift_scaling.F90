!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Implementation of a power-law outflow rate due to star formation feedback in galactic spheroids in which the characteristic velocity scales as a power of $(1+z)$.

  use Cosmology_Functions, only : cosmologyFunctionsClass, cosmologyFunctions
  
  !# <starFormationFeedbackSpheroids name="starFormationFeedbackSpheroidsPowerLawRedshiftScaling">
  !#  <description>A power-law outflow rate due to star formation feedback in galactic spheroids in which the characteristic velocity scales as a power of $(1+z)$.</description>
  !# </starFormationFeedbackSpheroids>
  type, extends(starFormationFeedbackSpheroidsPowerLaw) :: starFormationFeedbackSpheroidsPowerLawRedshiftScaling
     !% Implementation of a power-law outflow rate due to star formation feedback in galactic spheroids in which the characteristic velocity scales as a power of $(1+z)$.
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: exponentRedshift
   contains
     final     ::                           powerLawRedshiftScalingDestructor
     procedure :: velocityCharacteristic => powerLawRedshiftScalingVelocityCharacteristic
  end type starFormationFeedbackSpheroidsPowerLawRedshiftScaling

  interface starFormationFeedbackSpheroidsPowerLawRedshiftScaling
     !% Constructors for the power-law redshift-scaling star formation feedback in spheroids class.
     module procedure powerLawRedshiftScalingConstructorParameters
     module procedure powerLawRedshiftScalingConstructorInternal
  end interface starFormationFeedbackSpheroidsPowerLawRedshiftScaling

contains

  function powerLawRedshiftScalingConstructorParameters(parameters) result(self)
    !% Constructor for the power-law redshift-scaling star formation feedback in spheroids class which takes a parameter set as input.
    implicit none
    type(starFormationFeedbackSpheroidsPowerLawRedshiftScaling)                :: self
    type(inputParameters                                      ), intent(inout) :: parameters

    !# <inputParameter>
    !#   <name>exponentRedshift</name>
    !#   <source>parameters</source>
    !#   <variable>self%exponentRedshift</variable>
    !#   <defaultValue>0.0d0</defaultValue>
    !#   <description>The redshift scaling of characteristic velocity for \gls{sne}-driven outflow rates in spheroids.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="self%cosmologyFunctions_" source="parameters"/>
    self%starFormationFeedbackSpheroidsPowerLaw=starFormationFeedbackSpheroidsPowerLaw(parameters)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end function powerLawRedshiftScalingConstructorParameters

  function powerLawRedshiftScalingConstructorInternal(velocityCharacteristic_,exponent,exponentRedshift,cosmologyFunctions_) result(self)
    !% Internal constructor for the power-law redshift-scaling star formation feedback from spheroids class.
    implicit none
    type            (starFormationFeedbackSpheroidsPowerLawRedshiftScaling)                        :: self
    double precision                                                       , intent(in   )         :: velocityCharacteristic_, exponent, &
         &                                                                                            exponentRedshift
    class           (cosmologyFunctionsClass                              ), intent(in   ), target :: cosmologyFunctions_
    !# <constructorAssign variables="exponentRedshift, *cosmologyFunctions_"/>    

    self%starFormationFeedbackSpheroidsPowerLaw=starFormationFeedbackSpheroidsPowerLaw(velocityCharacteristic_,exponent)
    return
  end function powerLawRedshiftScalingConstructorInternal

  subroutine powerLawRedshiftScalingDestructor(self)
    !% Destructor for the power-law redshift-scaling feedback from star formation in spheroids class.
    implicit none
    type(starFormationFeedbackSpheroidsPowerLawRedshiftScaling), intent(inout) :: self
  
    !# <objectDestructor name="self%cosmologyFunctions_"/>
    return
  end subroutine powerLawRedshiftScalingDestructor
  
  double precision function powerLawRedshiftScalingVelocityCharacteristic(self,node)
    !% Return the characteristic velocity for power-law feedback models in spheroids. In this case the characteristic velocity is a
    !% constant.
    use Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(starFormationFeedbackSpheroidsPowerLawRedshiftScaling), intent(inout) :: self
    type (treeNode                                             ), intent(inout) :: node
    class(nodeComponentBasic                                   ), pointer       :: basic

    basic                                         =>  node%basic                                      (            )
    powerLawRedshiftScalingVelocityCharacteristic =  +self                    %velocityCharacteristic_                                      &
         &                                           /self%cosmologyFunctions_%expansionFactor        (basic%time())**self%exponentRedshift
    return
  end function powerLawRedshiftScalingVelocityCharacteristic
