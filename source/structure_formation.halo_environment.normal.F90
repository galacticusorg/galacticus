!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

!% Contains a module which implements a normally-distributed halo environment.

  use Cosmology_Parameters
  use Cosmology_Functions
  use Linear_Growth
  use Tables

  !# <haloEnvironment name="haloEnvironmentNormal">
  !#  <description>Implements a normally-distributed halo environment.</description>
  !# </haloEnvironment>
  type, extends(haloEnvironmentClass) :: haloEnvironmentNormal
     !% A normal halo environment class.
     private
     class           (cosmologyParametersClass     ), pointer :: cosmologyParameters_
     class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_
     class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_
     class           (linearGrowthClass            ), pointer :: linearGrowth_
     class           (criticalOverdensityClass     ), pointer :: criticalOverdensity_
     type            (table2DLinLinLin             )          :: linearToNonLinear
     double precision                                         :: radiusEnvironment        , variance
   contains
     final     ::                         normalDestructor
     procedure :: overdensityLinear    => normalOverdensityLinear
     procedure :: overdensityNonLinear => normalOverdensityNonLinear
     procedure :: environmentRadius    => normalEnvironmentRadius
     procedure :: environmentMass      => normalEnvironmentMass
  end type haloEnvironmentNormal

  interface haloEnvironmentNormal
     !% Constructors for the {\normalfont \ttfamily normal} halo environment class.
     module procedure normalConstructorParameters
     module procedure normalConstructorInternal
  end interface haloEnvironmentNormal

contains
  
  function normalConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily normal} halo environment class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (haloEnvironmentNormal        )                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), pointer       :: linearGrowth_
    class           (criticalOverdensityClass     ), pointer       :: criticalOverdensity_
    double precision                                               :: radiusEnvironment

    !# <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    !# <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !# <objectBuilder class="linearGrowth"             name="linearGrowth_"             source="parameters"/>
    !# <objectBuilder class="criticalOverdensity"      name="criticalOverdensity_"      source="parameters"/>
    !# <inputParameter>
    !#   <name>radiusEnvironment</name>
    !#   <source>parameters</source>
    !#   <variable>radiusEnvironment</variable>
    !#   <defaultValue>7.0d0</defaultValue>
    !#   <description>The radius of the sphere used to determine the variance in the environmental density.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=haloEnvironmentNormal(radiusEnvironment,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function normalConstructorParameters

  function normalConstructorInternal(radiusEnvironment,cosmologyParameters_,cosmologyFunctions_,cosmologicalMassVariance_,linearGrowth_,criticalOverdensity_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily normal} halo mass function class.
    use Spherical_Collapse_Matter_Lambda
    use Numerical_Constants_Math
    implicit none
    type            (haloEnvironmentNormal        )                        :: self
    class           (cosmologyParametersClass     ), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), target, intent(in   ) :: cosmologyFunctions_
    class           (cosmologicalMassVarianceClass), target, intent(in   ) :: cosmologicalMassVariance_
    class           (linearGrowthClass            ), target, intent(in   ) :: linearGrowth_
    class           (criticalOverdensityClass     ), target, intent(in   ) :: criticalOverdensity_
    double precision                                       , intent(in   ) :: radiusEnvironment
    !# <constructorAssign variables="radiusEnvironment, *cosmologyParameters_, *cosmologyFunctions_, *cosmologicalMassVariance_, *linearGrowth_, *criticalOverdensity_" />

    ! Find the root-variance in the linear density field on the given scale.
    self%variance=self%cosmologicalMassVariance_%rootVariance(                                                  &
         &                                                    +4.0d0                                            &
         &                                                    /3.0d0                                            &
         &                                                    *Pi                                               &
         &                                                    *self%cosmologyParameters_%OmegaMatter      ()    &
         &                                                    *self%cosmologyParameters_%densityCritical  ()    &
         &                                                    *self                     %radiusEnvironment  **3 &
         &                                                   )                                              **2
    ! Get a table of linear vs. nonlinear density.
    call Spherical_Collapse_Matter_Lambda_Nonlinear_Mapping(self%cosmologyFunctions_%cosmicTime(1.0d0),self%linearToNonLinear,self%linearGrowth_,self%cosmologyFunctions_)
    return
  end function normalConstructorInternal

  subroutine normalDestructor(self)
    !% Destructor for the {\normalfont \ttfamily normal} halo mass function class.
    implicit none
    type(haloEnvironmentNormal), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"      />
    !# <objectDestructor name="self%cosmologicalMassVariance_" />
    !# <objectDestructor name="self%linearGrowth_"             />
    return
  end subroutine normalDestructor

  double precision function normalOverdensityLinear(self,node,presentDay)
    !% Return the environment of the given {\normalfont \ttfamily node}.
    use Kind_Numbers
    use Statistics_Distributions
    implicit none
    class           (haloEnvironmentNormal     ), intent(inout)           :: self
    type            (treeNode                  ), intent(inout)           :: node
    logical                                     , intent(in   ), optional :: presentDay
    type            (treeNode                  ), pointer                 :: nodeRoot
    class           (nodeComponentBasic        ), pointer                 :: basic
    integer         (kind_int8                 ), save                    :: uniqueIDPrevious       =-1_kind_int8
    double precision                            , save                    :: overdensityPrevious
    !$omp threadprivate(uniqueIDPrevious,overdensityPrevious)
    double precision                            , parameter               :: overdensityMean        =0.0d+0
    double precision                            , parameter               :: limitUpperBuffer       =1.0d-3
    double precision                                                      :: overdensityVariance
    type            (distributionPeakBackground)                          :: distributionOverdensity
    !# <optionalArgument name="presentDay" defaultsTo=".false." />

    if (node%hostTree%baseNode%uniqueID() /= uniqueIDPrevious) then
       uniqueIDPrevious=node%hostTree%baseNode%uniqueID()
       if (node%hostTree%properties%exists('haloEnvironmentOverdensity')) then
          overdensityPrevious=node%hostTree%properties%value('haloEnvironmentOverdensity')
       else
          ! Find the root node of the tree.
          nodeRoot                =>  node%hostTree     %baseNode
          ! Find variance for the root node.
          overdensityVariance     =  +self%variance                                            &
               &                     *self%linearGrowth_%value      (expansionFactor=1.0d0)**2
          ! Construct the distribution for δ. This assumes a normal distribution for the densities, but conditioned on the fact
          ! that the region has not collapsed on any larger scale. The resulting distribution is given by eqn. (9) of Mo & White
          ! (1996; MNRAS; 282; 347). We include some small buffer to the collapse threshold to avoid rounding errors.
          distributionOverdensity =   distributionPeakBackground       (                                                           &
               &                                                        +     overdensityVariance                                , &
               &                                                        +self%criticalOverdensity_%value(expansionFactor=1.0d0)    &
               &                                                        *(                                                         &
               &                                                          +1.0d0                                                   &
               &                                                          -limitUpperBuffer                                        &
               &                                                         )                                                         &
               &                                                       )
          ! Choose an overdensity.
          overdensityPrevious     =  +distributionOverdensity   %sample(                                                           &
               &                                                        randomNumberGenerator= node%hostTree%randomNumberGenerator &
               &                                                       )
          call node%hostTree%properties%set('haloEnvironmentOverdensity',overdensityPrevious)
       end if
    end if
    normalOverdensityLinear=overdensityPrevious
    if (.not.presentDay_) then
       basic                   =>  node                                 %basic(                 )
       normalOverdensityLinear =  +normalOverdensityLinear                                        &
            &                     *self                   %linearGrowth_%value(time=basic%time())
    end if
    return
  end function normalOverdensityLinear

  double precision function normalOverdensityNonLinear(self,node)
    !% Return the environment of the given {\normalfont \ttfamily node}.
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self
    type (treeNode             ), intent(inout) :: node
    class(nodeComponentBasic   ), pointer       :: basic

    basic                      => node                  %basic      (                                         )
    normalOverdensityNonLinear =  self%linearToNonLinear%interpolate(self%overdensityLinear(node),basic%time())
    return
  end function normalOverdensityNonLinear

  double precision function normalEnvironmentRadius(self)
    !% Return the radius of the environment.
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    normalEnvironmentRadius=self%radiusEnvironment
    return
  end function normalEnvironmentRadius

  double precision function normalEnvironmentMass(self)
    !% Return the mass of the environment.
    use Numerical_Constants_Math
    implicit none
    class(haloEnvironmentNormal), intent(inout) :: self

    normalEnvironmentMass=+4.0d0                                                                   &
         &                *Pi                                                                      &
         &                *self%radiusEnvironment                                              **3 &
         &                *self%cosmologyFunctions_%matterDensityEpochal(expansionFactor=1.0d0)    &
         &                /3.0d0
    return
  end function normalEnvironmentMass

