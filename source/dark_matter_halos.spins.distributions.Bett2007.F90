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
  An implementation of the dark matter halo spin distribution which uses the fitting function proposed by
  \cite{bett_spin_2007}.
  !!}

  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Tables                 , only : table1D                 , table1DLogarithmicLinear

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
     class           (darkMatterHaloScaleClass ), pointer     :: darkMatterHaloScale_  => null()
     double precision                                         :: alpha                          , lambda0, &
          &                                                      normalization
     type            (table1DLogarithmicLinear )              :: distributionTable
     class           (table1D                  ), allocatable :: distributionInverse
   contains
     final     ::                 bett2007Destructor
     procedure :: sample       => bett2007Sample
     procedure :: distribution => bett2007Distribution
  end type haloSpinDistributionBett2007

  interface haloSpinDistributionBett2007
     !!{
     Constructors for the \refClass{haloSpinDistributionBett2007} dark matter halo spin
     distribution class.
     !!}
     module procedure bett2007ConstructorParameters
     module procedure bett2007ConstructorInternal
  end interface haloSpinDistributionBett2007

  ! Tabulation parameters.
  double precision, parameter :: countPointsPerDecade=200.0d+0
  double precision, parameter :: spinMaximum         =  1.0d+0
  double precision, parameter :: spinMinimum         =  1.0d-6

contains

  function bett2007ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloSpinDistributionBett2007} dark matter halo spin
    distribution class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloSpinDistributionBett2007)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    double precision                                              :: lambda0             , alpha

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
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=haloSpinDistributionBett2007(lambda0,alpha,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function bett2007ConstructorParameters

  function bett2007ConstructorInternal(lambda0,alpha,darkMatterHaloScale_) result(self)
    !!{
    Internal constructor for the \refClass{haloSpinDistributionBett2007} dark matter halo spin
    distribution class.
    !!}
    use :: Gamma_Functions, only : Gamma_Function      , Gamma_Function_Incomplete_Complementary
    use :: Table_Labels   , only : extrapolationTypeFix
    implicit none
    type            (haloSpinDistributionBett2007)                        :: self
    double precision                              , intent(in   )         :: lambda0             , alpha
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                                      :: spinDimensionless   , spinMaximum_
    integer                                                               :: iSpin               , countPointsTabulation
    logical                                                               :: success
    !![
    <constructorAssign variables="alpha, lambda0, *darkMatterHaloScale_"/>
    !!]
    
    ! Compute the normalization.
    self%normalization=+3.0d0                          &
         &             *self%alpha**(self%alpha-1.0d0) &
         &             /Gamma_Function(self%alpha)
    ! Tabulate the cumulative distribution.
    success     =.false.
    spinMaximum_=spinMaximum
    do while (.not.success)
       success=.true.
       countPointsTabulation=+int(                      &
            &                     +log10(               &
            &                            +spinMaximum_  &
            &                            /spinMinimum   &
            &                           )               &
            &                     *countPointsPerDecade &
            &                    )                      &
            &                +1
       call self%distributionTable%destroy()
       call self%distributionTable%create (                                                                   &
            &                              spinMinimum                                                      , &
            &                              spinMaximum_                                                     , &
            &                              countPointsTabulation                                            , &
            &                              extrapolationType    =[extrapolationTypeFix,extrapolationTypeFix]  &
            &                             )
       ! Compute the cumulative probability distribution.
       do iSpin=1,countPointsTabulation
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
             if (self%distributionTable%y(iSpin) <= self%distributionTable%y(iSpin-1)) then
                success     =.false.
                spinMaximum_=self%distributionTable%x(iSpin-1)
                exit
             end if
          end if
       end do
    end do
    call self%distributionTable%reverse(self%distributionInverse)
    return
  end function bett2007ConstructorInternal

  subroutine bett2007Destructor(self)
    !!{
    Destructor for the \refClass{haloSpinDistributionBett2007} dark matter halo spin
    distribution class.
    !!}
    implicit none
    type(haloSpinDistributionBett2007), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
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
    implicit none
    class(haloSpinDistributionBett2007), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    bett2007Sample=self%distributionInverse%interpolate(node%hostTree%randomNumberGenerator_%uniformSample())
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
    spin_                =  +spin%angularMomentum()                                                  &
         &                  /Dark_Matter_Halo_Angular_Momentum_Scale(node,self%darkMatterHaloScale_)
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
