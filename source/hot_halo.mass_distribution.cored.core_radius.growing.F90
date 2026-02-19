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
  An implementation of the hot halo mass distribution core radius class in which the core grows as the hot halo content is depleted.
  !!}

  use :: Cosmology_Parameters   , only : cosmologyParametersClass
  use :: Dark_Matter_Halo_Scales, only : darkMatterHaloScaleClass
  use :: Tables                 , only : table1D                 , table1DLogarithmicLinear

  !![
  <hotHaloMassDistributionCoreRadius name="hotHaloMassDistributionCoreRadiusGrowing">
   <description>
    A hot halo mass distribution core radius class which implements a core radius equal to a fraction {\normalfont \ttfamily
    [coreRadiusOverScaleRadius]} of the node's dark matter profile scale radius for nodes containing a mass of hot gas equal to
    the universal baryon fraction times their total mass. For nodes containing less hot gas mass, the core radius is expanded
    to maintain the same gas density at the virial radius, with a maximum core radius of {\normalfont \ttfamily
    [coreRadiusOverVirialRadiusMaximum]} times the node's virial radius.
   </description>
  </hotHaloMassDistributionCoreRadius>
  !!]
  type, extends(hotHaloMassDistributionCoreRadiusClass) :: hotHaloMassDistributionCoreRadiusGrowing
     !!{
     An implementation of the hot halo mass distribution core radius class in which the core grows as the hot halo content is depleted.
     !!}
     private
     class           (cosmologyParametersClass), pointer     :: cosmologyParameters_            => null()
     class           (darkMatterHaloScaleClass), pointer     :: darkMatterHaloScale_            => null()
     double precision                                        :: coreRadiusOverScaleRadius                 , coreRadiusOverVirialRadiusMaximum
     double precision                                        :: coreRadiusMaximum                         , coreRadiusMinimum
     double precision                                        :: hotGasFractionSaved                       , coreRadiusOverVirialRadiusInitialSaved, &
          &                                                     coreRadiusOverVirialRadiusSaved
     integer                                                 :: coreRadiusTableCount
     logical                                                 :: coreRadiusTableInitialized      =  .false.
     type            (table1DLogarithmicLinear)              :: coreRadiusTable
     class           (table1D                 ), allocatable :: coreRadiusTableInverse
   contains
     final     ::           growingDestructor
     procedure :: radius => growingRadius
  end type hotHaloMassDistributionCoreRadiusGrowing

  interface hotHaloMassDistributionCoreRadiusGrowing
     !!{
     Constructors for the \refClass{hotHaloMassDistributionCoreRadiusGrowing} hot halo mass distribution core radius class.
     !!}
     module procedure growingConstructorParameters
     module procedure growingConstructorInternal
  end interface hotHaloMassDistributionCoreRadiusGrowing

  integer, parameter :: tablePointsPerDecade=100

contains

  function growingConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{hotHaloMassDistributionCoreRadiusGrowing} hot halo mass distribution core radius class which builds the object
    from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloMassDistributionCoreRadiusGrowing)                :: self
    type            (inputParameters                         ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass                ), pointer       :: darkMatterHaloScale_
    class           (cosmologyParametersClass                ), pointer       :: cosmologyParameters_
    double precision                                                          :: coreRadiusOverScaleRadius, coreRadiusOverVirialRadiusMaximum

    !![
    <inputParameter>
      <name>coreRadiusOverScaleRadius</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The core radius in the hot halo density profile in units of the dark matter profile scale radius.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>coreRadiusOverVirialRadiusMaximum</name>
      <defaultValue>10.0d0</defaultValue>
      <description>The maximum core radius in the ``growing'' hot halo density profile in units of the virial radius.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !!]
    self=hotHaloMassDistributionCoreRadiusGrowing(coreRadiusOverScaleRadius,coreRadiusOverVirialRadiusMaximum,darkMatterHaloScale_,cosmologyParameters_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function growingConstructorParameters

  function growingConstructorInternal(coreRadiusOverScaleRadius,coreRadiusOverVirialRadiusMaximum,darkMatterHaloScale_,cosmologyParameters_) result(self)
    !!{
    Default constructor for the {\normalfont \ttfamily growing} hot halo mass distribution core radius class.
    !!}
    use :: Error           , only : Component_List                   , Error_Report
    use :: Galacticus_Nodes, only : defaultDarkMatterProfileComponent
    implicit none
    type            (hotHaloMassDistributionCoreRadiusGrowing)                        :: self
    double precision                                          , intent(in   )         :: coreRadiusOverScaleRadius, coreRadiusOverVirialRadiusMaximum
    class           (darkMatterHaloScaleClass                ), intent(in   ), target :: darkMatterHaloScale_
    class           (cosmologyParametersClass                ), intent(in   ), target :: cosmologyParameters_
    !![
    <constructorAssign variables="coreRadiusOverScaleRadius, coreRadiusOverVirialRadiusMaximum, *darkMatterHaloScale_, *cosmologyParameters_"/>
    !!]

    ! Ensure that the dark matter profile supports the scale property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                &
         & call Error_Report                                                                                     &
         &      (                                                                                                &
         &       'method requires a dark matter profile component that provides a gettable "scale" property.'//  &
         &       Component_List(                                                                                 &
         &                      'darkMatterProfile'                                                           ,  &
         &                       defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)   &
         &                     )                                                                              // &
         &       {introspection:location}                                                                        &
         &      )
    ! Initialize memoized values and table status.
    self%hotGasFractionSaved                   =-huge(0.0d0)
    self%coreRadiusOverVirialRadiusInitialSaved=-huge(0.0d0)
    self%coreRadiusOverVirialRadiusSaved       =-huge(0.0d0)
    self%coreRadiusTableInitialized            =.false.
    return
  end function growingConstructorInternal

  subroutine growingDestructor(self)
    !!{
    Destructor for the \refClass{hotHaloMassDistributionCoreRadiusGrowing} hot halo mass distribution class.
    !!}
    implicit none
    type(hotHaloMassDistributionCoreRadiusGrowing), intent(inout) :: self

    if (self%coreRadiusTableInitialized) call self%coreRadiusTable%destroy()
    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%cosmologyParameters_"/>
    !!]
    return
  end subroutine growingDestructor

  double precision function growingRadius(self,node)
    !!{
    Return the core radius of the hot halo mass distribution.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, nodeComponentDarkMatterProfile, nodeComponentHotHalo, treeNode
    implicit none
    class           (hotHaloMassDistributionCoreRadiusGrowing), intent(inout) :: self
    type            (treeNode                                ), intent(inout) :: node
    class           (nodeComponentBasic                      ), pointer       :: basic
    class           (nodeComponentHotHalo                    ), pointer       :: hotHalo
    class           (nodeComponentDarkMatterProfile          ), pointer       :: darkMatterProfile
    double precision                                                          :: hotGasFraction                   , coreRadiusOverVirialRadius, &
         &                                                                       coreRadiusOverVirialRadiusInitial, targetValue
    logical                                                                   :: makeTable

    ! Get components.
    basic             => node%basic            ()
    hotHalo           => node%hotHalo          ()
    darkMatterProfile => node%darkMatterProfile()
    ! Find the fraction of gas in the hot halo relative to that expected from the universal baryon fraction.
    hotGasFraction=+     hotHalo             %mass       () &
         &         /     basic               %mass       () &
         &         *self%cosmologyParameters_%OmegaMatter() &
         &         /self%cosmologyParameters_%OmegaBaryon()
    ! Return an arbitrary value for empty halos.
    if (hotGasFraction <= 0.0d0) then
       growingRadius=self%coreRadiusOverScaleRadius
       return
    end if
    ! Compute the desired core radius (in units of the virial radius) for a fully populated halo.
    coreRadiusOverVirialRadiusInitial=+self                     %coreRadiusOverScaleRadius       &
         &                            *     darkMatterProfile   %scale                    (    ) &
         &                            /self%darkMatterHaloScale_%radiusVirial             (node)
    ! Check if the initial core radius and hot gas fraction equal the previously stored values.
    if     (                                                                                        &
         &  .not.                                                                                   &
         &       (                                                                                  &
         &         coreRadiusOverVirialRadiusInitial == self%coreRadiusOverVirialRadiusInitialSaved &
         &        .and.                                                                             &
         &         hotGasFraction                    == self%hotGasFractionSaved                    &
         &       )                                                                                  &
         & ) then
       ! Create a tabulation of core radius vs. virial density factor if necessary.
       makeTable=.not.self%coreRadiusTableInitialized
       if (.not.makeTable) makeTable=(coreRadiusOverVirialRadiusInitial < self%coreRadiusTable%x(-1))
       if (makeTable) then
          self%coreRadiusMinimum   =min(self%coreRadiusOverScaleRadius,coreRadiusOverVirialRadiusInitial)
          self%coreRadiusMaximum   =self%coreRadiusOverVirialRadiusMaximum
          self%coreRadiusTableCount=int(log10(self%coreRadiusMaximum/self%coreRadiusMinimum)*dble(tablePointsPerDecade))+1
          call self%coreRadiusTable%destroy (                                                                       )
          call self%coreRadiusTable%create  (self%coreRadiusMaximum,self%coreRadiusMinimum,self%coreRadiusTableCount)
          call self%coreRadiusTable%populate(growingCoreVirialDensityFunction(self%coreRadiusTable%xs())            )
          call self%coreRadiusTable%reverse (self%coreRadiusTableInverse                                            )
          self%coreRadiusTableInitialized=.true.
       end if
       ! Compute the target value of the function giving the density at the virial radius per unit gas mass.
       targetValue=growingCoreVirialDensityFunction(coreRadiusOverVirialRadiusInitial)*hotGasFraction
       ! Interpolate to get the required core radius.
       if      (hotGasFraction >= 1.0d0                    ) then
          coreRadiusOverVirialRadius=coreRadiusOverVirialRadiusInitial
       else if (targetValue    <= self%coreRadiusTable%y(1)) then
          coreRadiusOverVirialRadius=self%coreRadiusTable%x(1)
       else
          coreRadiusOverVirialRadius=self%coreRadiusTableInverse%interpolate(targetValue)
       end if
       self%coreRadiusOverVirialRadiusInitialSaved=coreRadiusOverVirialRadiusInitial
       self%hotGasFractionSaved                   =hotGasFraction
       self%coreRadiusOverVirialRadiusSaved       =coreRadiusOverVirialRadius
    end if
    ! Compute the resulting core radius.
    growingRadius=+self                     %coreRadiusOverVirialRadiusSaved       &
         &        *self%darkMatterHaloScale_%radiusVirial                   (node)
    return
  end function growingRadius

  elemental double precision function growingCoreVirialDensityFunction(radiusOverVirialRadius)
    !!{
    Returns the function $(1+r_\mathrm{c}^2)[1-r_\mathrm{c} \tan^{-1}(1/r_\mathrm{c}]$ which is proportional to the density at the
    virial radius of a cored isothermal profile with core radius $r_\mathrm{c}$ (in units of the virial radius) per unit mass.
    !!}
    implicit none
    double precision, intent(in   ) :: radiusOverVirialRadius

    growingCoreVirialDensityFunction=+(                                 &
         &                             +      1.0d0                     &
         &                             +      radiusOverVirialRadius**2 &
         &                            )                                 &
         &                           *(                                 &
         &                             +      1.0d0                     &
         &                             -      radiusOverVirialRadius    &
         &                             *atan(                           &
         &                                   +1.0d0                     &
         &                                   /radiusOverVirialRadius    &
         &                                  )                           &
         &                            )
    return
  end function growingCoreVirialDensityFunction
