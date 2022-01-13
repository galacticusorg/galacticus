!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An implementation of the dark matter halo spin distribution which uses the fitting function proposed by
  \cite{bett_spin_2007}.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMO, darkMatterProfileDMOClass
  use :: Tables                  , only : table1D             , table1DLogarithmicLinear

  !![
  <haloSpinDistribution name="haloSpinDistributionBett2007">
   <description>
    A halo spin distribution in which the spin is drawn from the distribution found by \cite{bett_spin_2007}. The $\lambda_0$
    and $\alpha$ parameter of Bett et al.'s distribution are set by the {\normalfont \ttfamily [lambda0]} and {\normalfont
    \ttfamily [alpha]} input parameters.
   </description>
  </haloSpinDistribution>
  !!]
  type, extends(haloSpinDistributionClass) :: haloSpinDistributionBett2007
     !!{
     A dark matter halo spin distribution class which assumes a \cite{bett_spin_2007} distribution.
     !!}
     private
     class           (darkMatterProfileDMOClass), pointer     :: darkMatterProfileDMO_ => null()
     double precision                                         :: alpha                          , lambda0, &
          &                                                      normalization
     type            (table1DLogarithmicLinear )              :: distributionTable
     class           (table1D                  ), allocatable :: distributionInverse
     logical                                                  :: isInvertible
   contains
     final     ::                 bett2007Destructor
     procedure :: sample       => bett2007Sample
     procedure :: distribution => bett2007Distribution
  end type haloSpinDistributionBett2007

  interface haloSpinDistributionBett2007
     !!{
     Constructors for the {\normalfont \ttfamily bett2007} dark matter halo spin
     distribution class.
     !!}
     module procedure bett2007ConstructorParameters
     module procedure bett2007ConstructorInternal
  end interface haloSpinDistributionBett2007

  ! Tabulation parameters.
  integer         , parameter :: bett2007TabulationPointsCount=1000
  double precision, parameter :: bett2007SpinMaximum          =   0.2d+0
  double precision, parameter :: bett2007SpinMinimum          =   1.0d-6

contains

  function bett2007ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily bett2007} dark matter halo spin
    distribution class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloSpinDistributionBett2007)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterProfileDMOClass   ), pointer       :: darkMatterProfileDMO_
    double precision                                              :: lambda0              , alpha

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>lambda0</name>
      <source>parameters</source>
      <defaultValue>0.04326d0</defaultValue>
      <defaultSource>\citep{bett_spin_2007}</defaultSource>
      <description>The parameter $\lambda_0$ in the halo spin distribution of \cite{bett_spin_2007}.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha</name>
      <source>parameters</source>
      <defaultValue>2.509d0</defaultValue>
      <defaultSource>\citep{bett_spin_2007}</defaultSource>
      <description>The parameter $\alpha$ in the halo spin distribution of \cite{bett_spin_2007}.</description>
    </inputParameter>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=haloSpinDistributionBett2007(lambda0,alpha,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function bett2007ConstructorParameters

  function bett2007ConstructorInternal(lambda0,alpha,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily bett2007} dark matter halo spin
    distribution class.
    !!}
    use :: Gamma_Functions, only : Gamma_Function      , Gamma_Function_Incomplete_Complementary
    use :: Table_Labels   , only : extrapolationTypeFix
    implicit none
    type            (haloSpinDistributionBett2007)                        :: self
    double precision                              , intent(in   )         :: lambda0              , alpha
    class           (darkMatterProfileDMOClass   ), intent(in   ), target :: darkMatterProfileDMO_
    double precision                                                      :: spinDimensionless    , tableMaximum
    integer                                                               :: iSpin
    !![
    <constructorAssign variables="alpha, lambda0, *darkMatterProfileDMO_"/>
    !!]
    
    ! Compute the normalization.
    self%normalization=+3.0d0                          &
         &             *self%alpha**(self%alpha-1.0d0) &
         &             /Gamma_Function(self%alpha)
    ! Tabulate the cumulative distribution.
    tableMaximum=(bett2007SpinMaximum/lambda0)**(3.0d0/alpha)
    call self%distributionTable%destroy()
    call self%distributionTable%create (                                                                                        &
         &                              lambda0*bett2007SpinMinimum**(alpha/3.0d0),                                             &
         &                              lambda0*tableMaximum       **(alpha/3.0d0)                                            , &
         &                              bett2007TabulationPointsCount                                                         , &
         &                              extrapolationType                         =[extrapolationTypeFix,extrapolationTypeFix]  &
         &                             )
    ! Compute the cumulative probability distribution.
    self%isInvertible=.true.
    do iSpin=1,bett2007TabulationPointsCount
       spinDimensionless=(                                 &
            &             +self%distributionTable%x(iSpin) &
            &             /lambda0                         &
            &            )**(3.0d0/alpha)
       call self%distributionTable%populate(                                                            &
            &                               Gamma_Function_Incomplete_Complementary(                    &
            &                                                                       +alpha            , &
            &                                                                       +alpha              &
            &                                                                       *spinDimensionless  &
            &                                                                      )                  , &
            &                               iSpin                                                       &
            &                              )
       if (iSpin > 1) then
          if (self%distributionTable%y(iSpin) <= self%distributionTable%y(iSpin-1)) self%isInvertible=.false.
       end if
    end do
    if (self%isInvertible) call self%distributionTable%reverse(self%distributionInverse)
    return
  end function bett2007ConstructorInternal

  subroutine bett2007Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily bett2007} dark matter halo spin
    distribution class.
    !!}
    implicit none
    type(haloSpinDistributionBett2007), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    call                                          self%distributionTable  %destroy()
    if (allocated(self%distributionInverse)) call self%distributionInverse%destroy()
    return
  end subroutine bett2007Destructor

  double precision function bett2007Sample(self,node)
    !!{
    Sample from a \cite{bett_spin_2007} spin parameter distribution for the given {\normalfont
    \ttfamily node}.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(haloSpinDistributionBett2007), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    if (self%isInvertible) then
       bett2007Sample=self%distributionInverse%interpolate(node%hostTree%randomNumberGenerator_%uniformSample())
    else
       bett2007Sample=0.0d0
       call Galacticus_Error_Report('can not sample - cumulative distribution table was not monotonic'//{introspection:location})
    end if
    return
  end function bett2007Sample

  double precision function bett2007Distribution(self,node)
    !!{
    Compute the spin parameter distribution for the given {\normalfont \ttfamily node} assuming the fitting function of
    \cite{bett_spin_2007}.
    !!}
    use :: Dark_Matter_Halo_Spins, only : Dark_Matter_Halo_Angular_Momentum_Scale
    use :: Galacticus_Nodes      , only : nodeComponentSpin                      , treeNode
    implicit none
    class           (haloSpinDistributionBett2007), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    class           (nodeComponentSpin           ), pointer       :: spin
    double precision                                              :: spin_

    spin                 =>  node%spin           ()
    spin_                =  +spin%angularMomentum()                                                   &
         &                  /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterProfileDMO_)
    bett2007Distribution =  +self%normalization         &
         &                  *(                          &
         &                    +spin_                    &
         &                    /self%lambda0             &
         &                   )**3                       &
         &                  *exp(                       &
         &                       -self%alpha            &
         &                       *(                     &
         &                         +spin_               &
         &                         /self%lambda0        &
         &                        )**(3.0d0/self%alpha) &
         &                      )                       &
         &                  /spin_
    return
  end function bett2007Distribution
