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

  !% An implementation of dark matter halo virial density contrasts based on spherical collapse in a matter plus cosmological constant universe.
  
  use Tables                               , only : table1D
  use Cosmology_Parameters                 , only : cosmologyParameters             , cosmologyParametersClass
  use Cosmology_Functions                  , only : cosmologyFunctions              , cosmologyFunctionsClass
  use Intergalactic_Medium_Filtering_Masses, only : intergalacticMediumFilteringMass, intergalacticMediumFilteringMassClass

  !# <virialDensityContrast name="virialDensityContrastSphericalCollapseBaryonsDM">
  !#  <description>Dark matter halo virial density contrasts based on the spherical collapse in a matter plus cosmological constant universe.</description>
  !# </virialDensityContrast>
  type, extends(virialDensityContrastClass) :: virialDensityContrastSphericalCollapseBaryonsDM
     !% A dark matter halo virial density contrast class based on spherical collapse in a matter plus cosmological constant universe.
     private
     logical                                                              :: tableInitialized                  =  .false., turnaroundInitialized           =.false.
     double precision                                                     :: tableClusteredTimeMinimum                   , tableClusteredTimeMaximum               , &
          &                                                                  tableUnclusteredTimeMinimum                 , tableUnclusteredTimeMaximum             , &
          &                                                                  turnaroundClusteredTimeMinimum              , turnaroundClusteredTimeMaximum          , &
          &                                                                  turnaroundUnclusteredTimeMinimum            , turnaroundUnclusteredTimeMaximum
     logical                                                              :: tableStore
     integer                                                              :: energyFixedAt
     class           (table1D                              ), allocatable :: deltaVirialClustered                        , deltaVirialUnclustered                  , &
          &                                                                  turnaroundClustered                         , turnaroundUnclustered
     class           (cosmologyParametersClass             ), pointer     :: cosmologyParameters_              => null()
     class           (cosmologyFunctionsClass              ), pointer     :: cosmologyFunctions_               => null()
     class           (intergalacticMediumFilteringMassClass), pointer     :: intergalacticMediumFilteringMass_ => null()
   contains
     !@ <objectMethods>
     !@   <object>virialDensityContrastSphericalCollapseBaryonsDM</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse virial density contrast.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>retabulateTurnaround</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate spherical collapse turnaround radius.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                sphericalCollapseBaryonsDMDestructor
     procedure :: densityContrast             => sphericalCollapseBaryonsDMDensityContrast
     procedure :: densityContrastRateOfChange => sphericalCollapseBaryonsDMDensityContrastRateOfChange
     procedure :: turnAroundOverVirialRadii   => sphericalCollapseBaryonsDMTurnAroundOverVirialRadii
     procedure :: retabulate                  => sphericalCollapseBaryonsDMRetabulate
     procedure :: retabulateTurnaround        => sphericalCollapseBaryonsDMRetabulateTurnaround
  end type virialDensityContrastSphericalCollapseBaryonsDM

  interface virialDensityContrastSphericalCollapseBaryonsDM
     !% Constructors for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} dark matter halo virial density contrast class.
     module procedure sphericalCollapseBaryonsDMConstructorParameters
     module procedure sphericalCollapseBaryonsDMConstructorInternal
  end interface virialDensityContrastSphericalCollapseBaryonsDM

contains
  
  function sphericalCollapseBaryonsDMConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} dark matter halo virial density contrast class that takes a parameter set as input.
    use Input_Parameters                     , only : inputParameter                                           , inputParameters
    use Spherical_Collapse_Matter_Dark_Energy, only : enumerationDarkEnergySphericalCollapseEnergyFixedAtEncode
    implicit none
    type   (virialDensityContrastSphericalCollapseBaryonsDM)                :: self
    type   (inputParameters                                ), intent(inout) :: parameters
    class  (cosmologyParametersClass                       ), pointer       :: cosmologyParameters_
    class  (cosmologyFunctionsClass                        ), pointer       :: cosmologyFunctions_
    class  (intergalacticMediumFilteringMassClass          ), pointer       :: intergalacticMediumFilteringMass_
    logical                                                                 :: tableStore
    type   (varying_string                                 )                :: energyFixedAt

    !# <inputParameter>
    !#   <name>tableStore</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <description>If true, store/restore the tabulated solution to/from file when possible.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>energyFixedAt</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>var_str('turnaround')</defaultValue>
    !#   <description>Selects the epoch at which the energy of a spherical top hat perturbation in a dark energy cosmology should be
    !#     ``fixed'' for the purposes of computing virial density contrasts. (See the discussion in
    !#     \citealt{percival_cosmological_2005}; \S8.)</description>
    !#   <source>parameters</source>
    !#   <type>string</type>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters"              name="cosmologyParameters_"              source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"               name="cosmologyFunctions_"               source="parameters"/>
    !# <objectBuilder class="intergalacticMediumFilteringMass" name="intergalacticMediumFilteringMass_" source="parameters"/>
    self=virialDensityContrastSphericalCollapseBaryonsDM(tableStore,enumerationDarkEnergySphericalCollapseEnergyFixedAtEncode(char(energyFixedAt),includesPrefix=.false.),cosmologyParameters_,cosmologyFunctions_,intergalacticMediumFilteringMass_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="cosmologyParameters_"             />
    !# <objectDestructor name="cosmologyFunctions_"              />
    !# <objectDestructor name="intergalacticMediumFilteringMass_"/>
    return
  end function sphericalCollapseBaryonsDMConstructorParameters

  function sphericalCollapseBaryonsDMConstructorInternal(tableStore,energyFixedAt,cosmologyParameters_,cosmologyFunctions_,intergalacticMediumFilteringMass_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} dark matter halo virial density contrast class.
    implicit none
    type   (virialDensityContrastSphericalCollapseBaryonsDM)                        :: self
    class  (cosmologyParametersClass                       ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologyFunctionsClass                        ), intent(in   ), target :: cosmologyFunctions_
    class  (intergalacticMediumFilteringMassClass          ), intent(in   ), target :: intergalacticMediumFilteringMass_
    integer                                                 , intent(in   )         :: energyFixedAt
    logical                                                 , intent(in   )         :: tableStore
    !# <constructorAssign variables="tableStore, energyFixedAt, *cosmologyParameters_, *cosmologyFunctions_, *intergalacticMediumFilteringMass_"/>

    self%tableInitialized=.false.
    return
  end function sphericalCollapseBaryonsDMConstructorInternal
  
  subroutine sphericalCollapseBaryonsDMDestructor(self)
    !% Destructor for the {\normalfont \ttfamily sphericalCollapseBaryonsDM} dark matter halo virial density contrast class.
    implicit none
    type (virialDensityContrastSphericalCollapseBaryonsDM), intent(inout) :: self

    if (self%tableInitialized) then
       call self%deltaVirialClustered  %destroy()
       call self%deltaVirialUnclustered%destroy()
       deallocate(self%deltaVirialClustered  )
       deallocate(self%deltaVirialUnclustered)
    end if
    if (self%turnaroundInitialized) then
       call self%turnaroundClustered  %destroy()
       call self%turnaroundUnclustered%destroy()
       deallocate(self%turnaroundClustered  )
       deallocate(self%turnaroundUnclustered)
    end if
    !# <objectDestructor name="self%cosmologyParameters_"             />
    !# <objectDestructor name="self%cosmologyFunctions_"              />
    !# <objectDestructor name="self%intergalacticMediumFilteringMass_"/>
    return
  end subroutine sphericalCollapseBaryonsDMDestructor

  subroutine sphericalCollapseBaryonsDMRetabulate(self,time)
    !% Recompute the look-up tables for virial density contrast.
    use Spherical_Collapse_BDM, only : Spherical_Collapse_BDM_Virial_Density_Contrast_Tabulate
    implicit none
    class           (virialDensityContrastSphericalCollapseBaryonsDM), intent(inout) :: self
    double precision                                                 , intent(in   ) :: time
    logical                                                                          :: remakeTable

    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(time < self%tableClusteredTimeMinimum .or. time > self%tableClusteredTimeMaximum .or. time < self%tableUnclusteredTimeMinimum .or. time > self%tableUnclusteredTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_BDM_Virial_Density_Contrast_Tabulate(time,.false.,self%energyFixedAt,self%deltaVirialUnclustered,self%cosmologyParameters_,self%cosmologyFunctions_)
       call Spherical_Collapse_BDM_Virial_Density_Contrast_Tabulate(time,.true. ,self%energyFixedAt,self%deltaVirialClustered  ,self%cosmologyParameters_,self%cosmologyFunctions_)
       self%tableInitialized=.true.
       self%tableClusteredTimeMinimum  =self%deltaVirialClustered  %x(+1)
       self%tableClusteredTimeMaximum  =self%deltaVirialClustered  %x(-1)
       self%tableUnclusteredTimeMinimum=self%deltaVirialUnclustered%x(+1)
       self%tableUnclusteredTimeMaximum=self%deltaVirialUnclustered%x(-1)
    end if
    return
  end subroutine sphericalCollapseBaryonsDMRetabulate

  double precision function sphericalCollapseBaryonsDMDensityContrast(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseBaryonsDM), intent(inout)           :: self
    double precision                                                 , intent(in   )           :: mass
    double precision                                                 , intent(in   ), optional :: time      , expansionFactor
    logical                                                          , intent(in   ), optional :: collapsing
    double precision                                                                           :: time_     , interpolator
    
    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator=self%intergalacticMediumFilteringMass_%fractionBaryons(mass,time_)
    ! Interpolate the virial density contrast between the clustered and unclustered baryons case.
    sphericalCollapseBaryonsDMDensityContrast=+self%deltaVirialClustered  %interpolate(time_)*       interpolator  &
         &                                    +self%deltaVirialUnclustered%interpolate(time_)*(1.0d0-interpolator) 
    return
  end function sphericalCollapseBaryonsDMDensityContrast

  double precision function sphericalCollapseBaryonsDMDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based spherical collapse in a matter plus cosmological constant universe.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseBaryonsDM), intent(inout)           :: self
    double precision                                                 , intent(in   )           :: mass
    double precision                                                 , intent(in   ), optional :: time                    , expansionFactor
    logical                                                          , intent(in   ), optional :: collapsing
    double precision                                                                           :: time_                   , interpolator   , &
         &                                                                                        interpolatorRateOfChange

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator            =self%intergalacticMediumFilteringMass_%fractionBaryons            (mass,time_)
    interpolatorRateOfChange=self%intergalacticMediumFilteringMass_%fractionBaryonsRateOfChange(mass,time_)
    ! Interpolate the virial density contrast rate of change between the clustered and unclustered baryons case.
    sphericalCollapseBaryonsDMDensityContrastRateOfChange=+  self%deltaVirialClustered  %interpolateGradient(time_)*       interpolator              &
         &                                                +  self%deltaVirialUnclustered%interpolateGradient(time_)*(1.0d0-interpolator            ) &
         &                                                +(                                                                                         &
         &                                                  +self%deltaVirialClustered  %interpolate        (time_)                                  &
         &                                                  -self%deltaVirialUnclustered%interpolate        (time_)                                  &
         &                                                 )                                                                                         &
         &                                                *                                                                interpolatorRateOfChange
   return
  end function sphericalCollapseBaryonsDMDensityContrastRateOfChange

  subroutine sphericalCollapseBaryonsDMRetabulateTurnaround(self,time)
    !% Recompute the look-up tables for virial density contrast.
    use Spherical_Collapse_BDM, only : Spherical_Collapse_BDM_Turnaround_Radius_Tabulate
    implicit none
    class           (virialDensityContrastSphericalCollapseBaryonsDM), intent(inout) :: self
    double precision                                                 , intent(in   ) :: time
    logical                                                                          :: remakeTable

    ! Check if we need to recompute our table.
    if (self%turnaroundInitialized) then
       remakeTable=(time < self%turnaroundClusteredTimeMinimum .or. time > self%turnaroundClusteredTimeMaximum .or. time < self%turnaroundUnclusteredTimeMinimum .or. time > self%turnaroundUnclusteredTimeMaximum)
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       call Spherical_Collapse_BDM_Turnaround_Radius_Tabulate(time,.false.,self%energyFixedAt,self%turnaroundUnclustered,self%cosmologyParameters_,self%cosmologyFunctions_)
       call Spherical_Collapse_BDM_Turnaround_Radius_Tabulate(time,.true. ,self%energyFixedAt,self%turnaroundClustered  ,self%cosmologyParameters_,self%cosmologyFunctions_)
       self%turnaroundInitialized=.true.
       self%turnaroundClusteredTimeMinimum  =self%turnaroundClustered  %x(+1)
       self%turnaroundClusteredTimeMaximum  =self%turnaroundClustered  %x(-1)
       self%turnaroundUnclusteredTimeMinimum=self%turnaroundUnclustered%x(+1)
       self%turnaroundUnclusteredTimeMaximum=self%turnaroundUnclustered%x(-1)
    end if
    return
  end subroutine sphericalCollapseBaryonsDMRetabulateTurnaround

  double precision function sphericalCollapseBaryonsDMTurnAroundOverVirialRadii(self,mass,time,expansionFactor,collapsing)
    !% Return the ratio of turnaround and virial radii at the given epoch, based spherical collapse in a matter plus cosmological
    !% constant universe.
    use Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (virialDensityContrastSphericalCollapseBaryonsDM), intent(inout)           :: self
    double precision                                                 , intent(in   )           :: mass
    double precision                                                 , intent(in   ), optional :: time      , expansionFactor
    logical                                                          , intent(in   ), optional :: collapsing
    double precision                                                                           :: time_      , interpolator
    
    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulateTurnaround(time_)
    ! Construct an interpolation between the cases where baryons are clustered and unclustered. We use the same factor that
    ! appears in the suppression of baryonic accretion as a reasonable measure.
    interpolator=self%intergalacticMediumFilteringMass_%fractionBaryons(mass,time_)
    ! Interpolate the virial density contrast between the clustered and unclustered baryons case.
    sphericalCollapseBaryonsDMTurnAroundOverVirialRadii=+self%turnaroundClustered  %interpolate(time_)*       interpolator  &
         &                                              +self%turnaroundUnclustered%interpolate(time_)*(1.0d0-interpolator) 
    return
  end function sphericalCollapseBaryonsDMTurnAroundOverVirialRadii
