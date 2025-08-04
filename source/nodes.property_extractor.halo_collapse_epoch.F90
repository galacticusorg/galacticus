!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
  Implements a node property extractor which extracts the redshift of collapse for a halo using the definition of
  \cite{schneider_structure_2015}, which is based on the conditional first crossing distribution from excursion set theory.
  !!}
  
  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass, criticalOverdensityClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  
  !![
  <nodePropertyExtractor name="nodePropertyExtractorHaloCollapseEpoch">
   <description>
     A node property extractor which extracts the redshift of collapse, $z_\mathrm{c}$, for a halo using the definition of
     \cite{schneider_structure_2015}, which is based on the conditional first crossing distribution from excursion set theory:     
     \begin{equation}
      \delta_\mathrm{c}(z_\mathrm{c}) = \left( {\pi \over 2} \left[ \sigma^2(f M) - \sigma^2(M) \right]
      \right)^{1/2}+\delta_\mathrm{c})(z_0),
    \end{equation}
    where $\delta_\mathrm{c}(z)$ is the critical overdensity for collapse at redshift $z$, and $f$ is the fraction of a halo's
    mass assembled at formation time (given by the {\normalfont \ttfamily [massFractionFormation]} parameter.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorHaloCollapseEpoch
     !!{     
     A node property extractor which extracts the redshift of collapse for a halo using the definition of
     \cite{schneider_structure_2015}, which is based on the conditional first crossing distribution from excursion set theory.
     !!}
     private
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_      => null()
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()
     double precision                                         :: massFractionFormation
   contains
     final     ::                haloCollapseEpochDestructor
     procedure :: extract     => haloCollapseEpochExtract
     procedure :: name        => haloCollapseEpochName
     procedure :: description => haloCollapseEpochDescription
     procedure :: unitsInSI   => haloCollapseEpochUnitsInSI
  end type nodePropertyExtractorHaloCollapseEpoch

  interface nodePropertyExtractorHaloCollapseEpoch
     !!{
     Constructors for the \refClass{nodePropertyExtractorHaloCollapseEpoch} output extractor class.
     !!}
     module procedure haloCollapseEpochConstructorParameters
     module procedure haloCollapseEpochConstructorInternal
  end interface nodePropertyExtractorHaloCollapseEpoch

contains

  function haloCollapseEpochConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorHaloCollapseEpoch} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nodePropertyExtractorHaloCollapseEpoch)                :: self
    type            (inputParameters                       ), intent(inout) :: parameters
    class           (criticalOverdensityClass              ), pointer       :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass         ), pointer       :: cosmologicalMassvariance_
    class           (cosmologyFunctionsClass               ), pointer       :: cosmologyFunctions_
    double precision                                                        :: massFractionFormation

    !![
    <inputParameter>
      <name>massFractionFormation</name>
      <source>parameters</source>
      <defaultValue>0.05d0</defaultValue>
      <description>The fraction of a halo's mass assembled at ``formation'' in the halo concentration algorithm of \cite{schneider_structure_2015}.</description>
    </inputParameter>
    <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !!]
    self=nodePropertyExtractorHaloCollapseEpoch(massFractionFormation,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_)
    !![
    <objectDestructor name="criticalOverdensity_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function haloCollapseEpochConstructorParameters

  function haloCollapseEpochConstructorInternal(massFractionFormation,criticalOverdensity_,cosmologicalMassvariance_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorHaloCollapseEpoch} output extractor property extractor class.
    !!}
    implicit none
    type            (nodePropertyExtractorHaloCollapseEpoch)                        :: self
    class           (criticalOverdensityClass              ), intent(in   ), target :: criticalOverdensity_
    class           (cosmologicalMassVarianceClass         ), intent(in   ), target :: cosmologicalMassvariance_
    class           (cosmologyFunctionsClass               ), intent(in   ), target :: cosmologyFunctions_
    double precision                                        , intent(in   )         :: massFractionFormation
    !![
    <constructorAssign variables="massFractionFormation, *criticalOverdensity_, *cosmologicalMassvariance_, *cosmologyFunctions_"/>
    !!]

    return
  end function haloCollapseEpochConstructorInternal

  subroutine haloCollapseEpochDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorHaloCollapseEpoch} output extractor property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorHaloCollapseEpoch), intent(inout) :: self

    !![
    <objectDestructor name="self%criticalOverdensity_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    <objectDestructor name="self%cosmologyFunctions_"      />
    !!]
    return
  end subroutine haloCollapseEpochDestructor

  double precision function haloCollapseEpochExtract(self,node,instance) result(redshiftCollapse)
    !!{
    Extract the redshift of collapse.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    implicit none
    class           (nodePropertyExtractorHaloCollapseEpoch), intent(inout), target   :: self
    type            (treeNode                              ), intent(inout), target   :: node
    type            (multiCounter                          ), intent(inout), optional :: instance
    class           (nodeComponentBasic                    )               , pointer  :: basic
    double precision                                                                  :: mass                       , variance        , &
         &                                                                               time_                      , timeCollapse    , &
         &                                                                               collapseCriticalOverdensity, timeNow
    !$GLC attributes unused :: instance

    ! Get the basic component and the halo mass and time.
    basic => node %basic           ()
    mass  =  basic%mass            ()
    time_ =  basic%timeLastIsolated()
    ! Find critical overdensity at collapse for this node. The critical overdensity at collapse is scaled by a factor
    ! σ(M,t₀)/σ(M,t), as the "timeOfCollapse" function expects a critical overdensity divided by the linear growth factor. (This
    ! is because the definition of critical overdensity in Galacticus does not include the 1/D(t) factor that is included in the
    ! definition used by Schneider et al. (2015).)
    timeNow                    =     self%cosmologyFunctions_      %cosmicTime  (1.0d0                                )
    variance                   =max(                                                                                        &
         &                          +0.0d0                                                                                , &
         &                          +self%cosmologicalMassVariance_%rootVariance(mass*self%massFractionFormation,time_)**2  &
         &                          -self%cosmologicalMassVariance_%rootVariance(mass                           ,time_)**2  &
         &                         )
    collapseCriticalOverdensity=+(                                                                                 &
         &                        +sqrt(                                                                           &
         &                              +Pi                                                                        &
         &                              /2.0d0                                                                     &
         &                              *variance                                                                  &
         &                             )                                                                           &
         &                        +  self%criticalOverdensity_     %value       (time=time_  ,mass=mass,node=node) &
         &                       )                                                                                 &
         &                      *    self%cosmologicalMassVariance_%rootVariance(time=timeNow,mass=mass          ) &
         &                      /    self%cosmologicalMassVariance_%rootVariance(time=time_  ,mass=mass          )

    ! Compute the corresponding epoch of collapse.
    timeCollapse         =self%criticalOverdensity_%timeOfCollapse              (collapseCriticalOverdensity,mass,node)
    redshiftCollapse     =self%cosmologyFunctions_ %redshiftFromExpansionFactor(                                        &
         &                self%cosmologyFunctions_ %expansionFactor             (                                       &
         &                                                                       timeCollapse                           &
         &                                                                      )                                       &
         &                                                                     )
    return
  end function haloCollapseEpochExtract
  
  function haloCollapseEpochName(self)
    !!{
    Return the names of the {\normalfont \ttfamily haloCollapseEpoch} properties.
    !!}
    implicit none
    type (varying_string                        )                :: haloCollapseEpochName
    class(nodePropertyExtractorHaloCollapseEpoch), intent(inout) :: self
    !$GLC attributes unused :: self

    haloCollapseEpochName=var_str('haloCollapseRedshift')
    return
  end function haloCollapseEpochName

  function haloCollapseEpochDescription(self)
    !!{
    Return the descriptions of the {\normalfont \ttfamily haloCollapseEpoch} properties.
    !!}
    implicit none
    type (varying_string                        )                :: haloCollapseEpochDescription
    class(nodePropertyExtractorHaloCollapseEpoch), intent(inout) :: self
    !$GLC attributes unused :: self

    haloCollapseEpochDescription=var_str('The redshift of halo collapse..')
    return
  end function haloCollapseEpochDescription

  double precision function haloCollapseEpochUnitsInSI(self)
    !!{
    Return the units of the {\normalfont \ttfamily haloCollapseEpoch} properties in the SI system.
    !!}
    implicit none
    class(nodePropertyExtractorHaloCollapseEpoch), intent(inout) :: self
    !$GLC attributes unused :: self

    haloCollapseEpochUnitsInSI=1.0d0
    return
  end function haloCollapseEpochUnitsInSI
