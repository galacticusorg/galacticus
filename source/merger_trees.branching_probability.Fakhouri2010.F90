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

!% Implements a merger tree branching probability class using the algorithm of \cite{fakhouri_merger_2010}.

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass

  !# <mergerTreeBranchingProbability name="mergerTreeBranchingProbabilityFakhouri2010">
  !#  <description>Merger tree branching probabilities using the algorithm of \cite{fakhouri_merger_2010}.</description>
  !# </mergerTreeBranchingProbability>
  type, extends(mergerTreeBranchingProbabilityClass) :: mergerTreeBranchingProbabilityFakhouri2010
     !% A merger tree branching probability class using the algorithm of \cite{fakhouri_merger_2010}.
     private
     double precision                                         :: alpha                              , beta , &
          &                                                      gamma                              , eta  , &
          &                                                      A                                  , xiBar
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
   contains
     final     ::         fakhouri2010Destructor
     procedure :: rate => fakhouri2010Rate
  end type mergerTreeBranchingProbabilityFakhouri2010

  interface mergerTreeBranchingProbabilityFakhouri2010
     !% Constructors for the {\normalfont \ttfamily fakhouri2010} merger tree builder class.
     module procedure fakhouri2010ConstructorParameters
     module procedure fakhouri2010ConstructorInternal
  end interface mergerTreeBranchingProbabilityFakhouri2010

contains

  function fakhouri2010ConstructorParameters(parameters) result(self)
    !% Constructor for the ``fakhouri2010'' merger tree branching probability class which reads parameters from a provided
    !% parameter list.
    implicit none
    type            (mergerTreeBranchingProbabilityFakhouri2010)                :: self
    type            (inputParameters                           ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                   ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass             ), pointer       :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                  ), pointer       :: criticalOverdensity_
    double precision                                                            :: alpha                    , beta , &
         &                                                                         gamma                    , eta  , &
         &                                                                         A                        , xiBar

    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <defaultValue>0.133d0</defaultValue>
    !#   <defaultSource>\citep{fakhouri_merger_2010}</defaultSource>
    !#   <description>The parameter $\alpha$ appearing in equation~(1) of \cite{fakhouri_merger_2010}.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <defaultValue>-1.995d0</defaultValue>
    !#   <defaultSource>\citep{fakhouri_merger_2010}</defaultSource>
    !#   <description>The parameter $\beta$ appearing in equation~(1) of \cite{fakhouri_merger_2010}.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>gamma</name>
    !#   <defaultValue>0.263d0</defaultValue>
    !#   <defaultSource>\citep{fakhouri_merger_2010}</defaultSource>
    !#   <description>The parameter $\gamma$ appearing in equation~(1) of \cite{fakhouri_merger_2010}.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>eta</name>
    !#   <defaultValue>0.0993d0</defaultValue>
    !#   <defaultSource>\citep{fakhouri_merger_2010}</defaultSource>
    !#   <description>The parameter $\eta$ appearing in equation~(1) of \cite{fakhouri_merger_2010}.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>A</name>
    !#   <defaultValue>0.0104d0</defaultValue>
    !#   <defaultSource>\citep{fakhouri_merger_2010}</defaultSource>
    !#   <description>The parameter $A$ appearing in equation~(1) of \cite{fakhouri_merger_2010}.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>xiBar</name>
    !#   <defaultValue>9.72d-3</defaultValue>
    !#   <defaultSource>\citep{fakhouri_merger_2010}</defaultSource>
    !#   <description>The parameter $\bar{\xi}$ appearing in equation~(1) of \cite{fakhouri_merger_2010}.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    self=mergerTreeBranchingProbabilityFakhouri2010(alpha,beta,gamma,eta,A,xiBar,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyFunctions_"      />
    !# <objectDestructor name="cosmologicalMassVariance_"/>
    !# <objectDestructor name="criticalOverdensity_"     />
    return
  end function fakhouri2010ConstructorParameters

  function fakhouri2010ConstructorInternal(alpha,beta,gamma,eta,A,xiBar,cosmologyFunctions_,cosmologicalMassVariance_,criticalOverdensity_) result(self)
    !% Internal constructor for the ``fakhouri2010'' merger tree branching probability class.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    type            (mergerTreeBranchingProbabilityFakhouri2010)                        :: self
    double precision                                            , intent(in   )         :: alpha                    , beta , &
         &                                                                                 gamma                    , eta  , &
         &                                                                                 A                        , xiBar
    class           (cosmologyFunctionsClass                   ), intent(in   ), target :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass             ), intent(in   ), target :: cosmologicalMassVariance_
    class           (criticalOverdensityClass                  ), intent(in   ), target :: criticalOverdensity_
    !# <constructorAssign variables="alpha, beta, gamma, eta, A, xiBar, *cosmologyFunctions_, *cosmologicalMassVariance_, *criticalOverdensity_"/>
    
    return
  end function fakhouri2010ConstructorInternal

  subroutine fakhouri2010Destructor(self)
    implicit none
    type(mergerTreeBranchingProbabilityFakhouri2010), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyFunctions_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_"/>
    !# <objectDestructor name="self%criticalOverdensity_"     />
    return
  end subroutine fakhouri2010Destructor

  double precision function fakhouri2010Rate(self,mass,deltaCritical,time,massBranch,node)
    !% Return the rate per unit mass and per unit change in $\delta_\mathrm{crit}$ that a halo of mass {\normalfont \ttfamily haloMass} at time
    !% {\normalfont \ttfamily deltaCritical} will undergo a branching to progenitors with mass {\normalfont \ttfamily massBranch}.
    implicit none
    class           (mergerTreeBranchingProbabilityFakhouri2010), intent(inout), target :: self
    double precision                                            , intent(in   )         :: deltaCritical                      , mass                             , &
         &                                                                                 massBranch                         , time
    type            (treeNode                                  ), intent(inout), target :: node
    double precision                                            , parameter             :: massReference               =1.0d12
    double precision                                                                    :: massRatio                          , redshift                         , &
         &                                                                                 expansionFactor                    , rootVarianceLogarithmicGrowthRate, &
         &                                                                                 deltaCritEffectiveGrowthRate       , deltaCriticalGrowthRate

    ! Evaluate the fitting function of Fakhouri, Ma, & Boylan-Kolchin (2010): d²N/dξ/dz.
    massRatio       =      +massBranch  &
         &           /(mass-massBranch)
    expansionFactor =+self%cosmologyFunctions_%expansionFactor            (time           )
    redshift        =+self%cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)
    fakhouri2010Rate=+self%A                                          &
         &           *    (mass     /     massReference)**self%alpha  &
         &           *     massRatio                    **self%beta   &
         &           *exp((massRatio/self%xiBar        )**self%gamma) &
         &           *    (1.0d0    +     redshift     )**self%eta
    ! Evaluate the rate of change of the effective barrier.
    deltaCriticalGrowthRate          =self%criticalOverdensity_     %gradientTime                       (time=time,mass=mass,node=node)
    rootVarianceLogarithmicGrowthRate=self%cosmologicalMassVariance_%rootVarianceLogarithmicGradientTime(time=time,mass=mass          )
    deltaCritEffectiveGrowthRate     =abs(                                   &
         &                                +deltaCriticalGrowthRate           &
         &                                -deltaCritical                     &
         &                                *rootVarianceLogarithmicGrowthRate &
         &                                /time                              &
         &                             )
    ! Convert the Fakhouri, Ma, & Boylan-Kolchin (2010) merger rate, d²N/dξ/dz, to the form we want, d²N/dB/dm, where B=δ/σ is the barrier.
    fakhouri2010Rate=+fakhouri2010Rate                                        &
         &           *self%cosmologyFunctions_%expansionRate(expansionFactor) &
         &           /                                       expansionFactor  &
         &           /deltaCritEffectiveGrowthRate                            &
         &           *  mass                                                  &
         &           /(                                                       &
         &             +mass                                                  &
         &             -massBranch                                            &
         &            )**2
    return
  end function fakhouri2010Rate
