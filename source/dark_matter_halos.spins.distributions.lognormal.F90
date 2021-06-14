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

  !% An implementation of the dark matter halo spin distribution which assumes a
  !% log-normal distribution.

  !# <haloSpinDistribution name="haloSpinDistributionLogNormal">
  !#  <description>
  !#   A halo spin distribution class in which the spin is drawn from a lognormal distribution with median {\normalfont \ttfamily
  !#   [median]} and width {\normalfont \ttfamily [sigma]}.
  !#  </description>
  !# </haloSpinDistribution>
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionLogNormal
     !% A dark matter halo spin distribution concentration class which assumes a
     !% log-normal distribution.
     private
     double precision :: median, sigma
   contains
     procedure :: sample        => logNormalSample
     procedure :: distribution  => logNormalDistribution
  end type haloSpinDistributionLogNormal

  interface haloSpinDistributionLogNormal
     !% Constructors for the {\normalfont \ttfamily logNormal} dark matter halo spin
     !% distribution class.
     module procedure logNormalConstructorParameters
     module procedure logNormalConstructorInternal
  end interface haloSpinDistributionLogNormal

contains

  function logNormalConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily logNormal} dark matter halo spin
    !% distribution class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(haloSpinDistributionLogNormal)                :: logNormalConstructorParameters
    type(inputParameters              ), intent(inout) :: parameters

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>median</name>
    !#   <source>parameters</source>
    !#   <variable>logNormalConstructorParameters%median</variable>
    !#   <defaultValue>0.03687d0</defaultValue>
    !#   <defaultSource>\citep{bett_spin_2007}</defaultSource>
    !#   <description>The median spin in a log-normal spin distribution.</description>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>sigma</name>
    !#   <source>parameters</source>
    !#   <variable>logNormalConstructorParameters%sigma</variable>
    !#   <defaultValue>0.5102d0</defaultValue>
    !#   <defaultSource>(\citealt{bett_spin_2007}; note that in this reference the value of $\sigma$ quoted is for $\log_{10}\lambda$, while here we use $\log\lambda$)</defaultSource>
    !#   <description>The width of a log-normal spin distribution.</description>
    !# </inputParameter>
    logNormalConstructorParameters%median=log(logNormalConstructorParameters%median)
    !# <inputParametersValidate source="parameters"/>
    return
  end function logNormalConstructorParameters

  function logNormalConstructorInternal(median,sigma)
    !% Internal constructor for the {\normalfont \ttfamily logNormal} dark matter halo spin
    !% distribution class.
    implicit none
    type            (haloSpinDistributionLogNormal)                :: logNormalConstructorInternal
    double precision                               , intent(in   ) :: median                      , sigma

    logNormalConstructorInternal%median=log(median)
    logNormalConstructorInternal%sigma =    sigma
    return
  end function logNormalConstructorInternal

  double precision function logNormalSample(self,node)
    !% Sample from a log-normal spin parameter distribution for the given {\normalfont
    !% \ttfamily node}.
    implicit none
    class(haloSpinDistributionLogNormal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    !$GLC attributes unused :: node

    logNormalSample=exp(                                                             &
         &              +self%median                                                 &
         &              +self%sigma                                                  &
         &              *node%hostTree%randomNumberGenerator_%standardNormalSample() &
         &             )
    return
  end function logNormalSample

  double precision function logNormalDistribution(self,node)
    !% Return the spin parameter distribution for the given {\normalfont \ttfamily node}
    !% assuming a log-normal distribution.
    use :: Galacticus_Nodes        , only : nodeComponentSpin, treeNode
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(haloSpinDistributionLogNormal), intent(inout) :: self
    type (treeNode                     ), intent(inout) :: node
    class(nodeComponentSpin            ), pointer       :: spin

    spin                  => node%spin()
    logNormalDistribution =  +exp(                    &
         &                        -(                  &
         &                          +log(spin%spin()) &
         &                          -self%median      &
         &                         )**2               &
         &                        /2.0d0              &
         &                        /self%sigma**2      &
         &                       )                    &
         &                   /sqrt(                   &
         &                         +2.0d0             &
         &                         *Pi                &
         &                        )                   &
         &                   /self%sigma              &
         &                   /spin%spin()
    return
  end function logNormalDistribution
