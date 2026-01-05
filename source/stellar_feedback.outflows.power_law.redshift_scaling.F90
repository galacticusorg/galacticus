!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
  Implementation of a power-law stellar feedback model in which the characteristic velocity scales as a power of $(1+z)$.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass

  !![
  <stellarFeedbackOutflows name="stellarFeedbackOutflowsPowerLawRedshiftScaling">
   <description>A power-law stellar feedback model in which the characteristic velocity scales as a power of $(1+z)$.</description>
  </stellarFeedbackOutflows>
  !!]
  type, extends(stellarFeedbackOutflowsPowerLaw) :: stellarFeedbackOutflowsPowerLawRedshiftScaling
     !!{
     Implementation of a power-law stellar feedback model in which the characteristic velocity scales as a power of $(1+z)$.
     !!}
     private
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_ => null()
     double precision                                   :: exponentRedshift
   contains
     final     ::                           powerLawRedshiftScalingDestructor
     procedure :: velocityCharacteristic => powerLawRedshiftScalingVelocityCharacteristic
  end type stellarFeedbackOutflowsPowerLawRedshiftScaling

  interface stellarFeedbackOutflowsPowerLawRedshiftScaling
     !!{
     Constructors for the power-law redshift-scaling stellar feedback class.
     !!}
     module procedure powerLawRedshiftScalingConstructorParameters
     module procedure powerLawRedshiftScalingConstructorInternal
  end interface stellarFeedbackOutflowsPowerLawRedshiftScaling

contains

  function powerLawRedshiftScalingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the power-law redshift-scaling stellar feedback class which takes a parameter set as input.
    !!}
    implicit none
    type(stellarFeedbackOutflowsPowerLawRedshiftScaling)                :: self
    type(inputParameters                               ), intent(inout) :: parameters

    !![
    <inputParameter>
      <name>exponentRedshift</name>
      <source>parameters</source>
      <variable>self%exponentRedshift</variable>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift scaling of characteristic velocity for \gls{sne}-driven outflow rates in disks.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="self%cosmologyFunctions_" source="parameters"/>
    !!]
    self%stellarFeedbackOutflowsPowerLaw=stellarFeedbackOutflowsPowerLaw(parameters)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function powerLawRedshiftScalingConstructorParameters

  function powerLawRedshiftScalingConstructorInternal(velocityCharacteristic_,exponent,exponentRedshift,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the power-law redshift-scaling stellar feedback class.
    !!}
    implicit none
    type            (stellarFeedbackOutflowsPowerLawRedshiftScaling)                        :: self
    double precision                                                , intent(in   )         :: velocityCharacteristic_, exponent, &
         &                                                                                     exponentRedshift
    class           (cosmologyFunctionsClass                       ), intent(in   ), target :: cosmologyFunctions_
    !![
    <constructorAssign variables="exponentRedshift, *cosmologyFunctions_"/>
    !!]

    self%stellarFeedbackOutflowsPowerLaw=stellarFeedbackOutflowsPowerLaw(velocityCharacteristic_,exponent)
    return
  end function powerLawRedshiftScalingConstructorInternal

  subroutine powerLawRedshiftScalingDestructor(self)
    !!{
    Destructor for the power-law redshift-scaling feedback from star formation in disks class.
    !!}
    implicit none
    type(stellarFeedbackOutflowsPowerLawRedshiftScaling), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine powerLawRedshiftScalingDestructor

  double precision function powerLawRedshiftScalingVelocityCharacteristic(self,component)
    !!{
    Return the characteristic velocity for power-law feedback models in disks. In this case the characteristic velocity is a
    constant.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(stellarFeedbackOutflowsPowerLawRedshiftScaling), intent(inout) :: self
    class(nodeComponent                                 ), intent(inout) :: component
    class(nodeComponentBasic                            ), pointer       :: basic

    basic                                         =>  component%hostNode           %basic                  (            )
    powerLawRedshiftScalingVelocityCharacteristic =  +self                         %velocityCharacteristic_                                      &
         &                                           /self     %cosmologyFunctions_%expansionFactor        (basic%time())**self%exponentRedshift
    return
  end function powerLawRedshiftScalingVelocityCharacteristic
