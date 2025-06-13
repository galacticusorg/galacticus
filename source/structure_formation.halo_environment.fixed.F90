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
Implements a fixed halo environment.
!!}

  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Linear_Growth             , only : linearGrowthClass
  use :: Spherical_Collapse_Solvers, only : sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt
  use :: Tables                    , only : table2DLinLinLin

  !![
  <haloEnvironment name="haloEnvironmentFixed">
   <description>Implements a fixed halo environment.</description>
   <deepCopy>
    <functionClass variables="sphericalCollapseSolver_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="sphericalCollapseSolver_"/>
   </stateStorable>
  </haloEnvironment>
  !!]
  type, extends(haloEnvironmentClass) :: haloEnvironmentFixed
     !!{
     A fixed halo environment class.
     !!}
     private
     class           (cosmologyFunctionsClass                          ), pointer :: cosmologyFunctions_          => null()
     class           (linearGrowthClass                                ), pointer :: linearGrowth_                => null()
     type            (sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt ), pointer :: sphericalCollapseSolver_     => null()
     type            (table2DLinLinLin                                 )          :: linearToNonLinear
     double precision                                                             :: radiusEnvironment                     , massEnvironment, &
          &                                                                          overdensity
     logical                                                                      :: linearToNonLinearInitialized
   contains
     final     ::                                  fixedHEDestructor
     procedure :: overdensityLinear             => fixedHEOverdensityLinear
     procedure :: overdensityLinearGradientTime => fixedHEOverdensityLinearGradientTime
     procedure :: overdensityNonLinear          => fixedHEOverdensityNonLinear
     procedure :: environmentRadius             => fixedHEEnvironmentRadius
     procedure :: environmentMass               => fixedHEEnvironmentMass
     procedure :: overdensityLinearMaximum      => fixedHEOverdensityLinearMaximum
     procedure :: pdf                           => fixedHEPDF
     procedure :: cdf                           => fixedHECDF
     procedure :: overdensityLinearSet          => fixedHEOverdensityLinearSet
     procedure :: overdensityIsSettable         => fixedHEOverdensityIsSettable
     procedure :: isNodeDependent               => fixedHEIsNodeDependent
     procedure :: isTreeDependent               => fixedHEIsTreeDependent
  end type haloEnvironmentFixed

  interface haloEnvironmentFixed
     !!{
     Constructors for the \refClass{haloEnvironmentFixed} halo environment class.
     !!}
     module procedure fixedHEConstructorParameters
     module procedure fixedHEConstructorInternal
  end interface haloEnvironmentFixed

contains

  function fixedHEConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{haloEnvironmentFixed} halo environment class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (haloEnvironmentFixed   )                :: self
    type            (inputParameters        ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass), pointer       :: cosmologyFunctions_
    class           (linearGrowthClass      ), pointer       :: linearGrowth_
    double precision                                         :: radiusEnvironment  , massEnvironment, &
         &                                                      overdensity

    !![
    <inputParameter>
      <name>overdensity</name>
      <source>parameters</source>
      <description>The overdensity of the environment.</description>
    </inputParameter>
    <inputParameter>
      <name>radiusEnvironment</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The radius of the sphere used to determine the variance in the environmental density.</description>
    </inputParameter>
    <inputParameter>
      <name>massEnvironment</name>
      <source>parameters</source>
      <defaultValue>1.0d15</defaultValue>
      <description>The mass within the sphere sphere used to determine the variance in the environmental density.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    <objectBuilder class="linearGrowth"       name="linearGrowth_"       source="parameters"/>
    <conditionalCall>
     <call>self=haloEnvironmentFixed(cosmologyFunctions_,linearGrowth_,overdensity{conditions})</call>
     <argument name="massEnvironment"   value="massEnvironment"   parameterPresent="parameters"/>
     <argument name="radiusEnvironment" value="radiusEnvironment" parameterPresent="parameters"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="linearGrowth_"            />
    !!]
    return
  end function fixedHEConstructorParameters

  function fixedHEConstructorInternal(cosmologyFunctions_,linearGrowth_,overdensity,radiusEnvironment,massEnvironment) result(self)
    !!{
    Internal constructor for the \refClass{haloEnvironmentFixed} halo mass function class.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (haloEnvironmentFixed   )                           :: self
    class           (cosmologyFunctionsClass) , intent(in   ), target   :: cosmologyFunctions_
    class           (linearGrowthClass      ) , intent(in   ), target   :: linearGrowth_
    double precision                          , intent(in   )           :: overdensity
    double precision                          , intent(in   ), optional :: radiusEnvironment  , massEnvironment
    !![
    <constructorAssign variables="*cosmologyFunctions_, *linearGrowth_, overdensity, radiusEnvironment, massEnvironment"/>
    !!]

    ! Construct a spherical collapse solver.
    allocate(self%sphericalCollapseSolver_)
    !![
    <referenceConstruct owner="self" isResult="yes" object="sphericalCollapseSolver_" constructor="sphericalCollapseSolverCllsnlssMttrCsmlgclCnstnt(self%cosmologyFunctions_,self%linearGrowth_)"/>
    !!]

    
    if (present(radiusEnvironment).and.present(massEnvironment)) call Error_Report('only one of radiusEnvironment and massEnvironment may be specified'//{introspection:location})
    if (present(radiusEnvironment)) then
       self%massEnvironment=+4.0d0                                                                   &
            &               *Pi                                                                      &
            &               *self%radiusEnvironment                                              **3 &
            &               *self%cosmologyFunctions_%matterDensityEpochal(expansionFactor=1.0d0)    &
            &               /3.0d0
    else if (present(massEnvironment)) then
       self%radiusEnvironment=+(                                                                      &
            &                   +3.0d0                                                                &
            &                   *self%massEnvironment                                                 &
            &                   /self%cosmologyFunctions_%matterDensityEpochal(expansionFactor=1.0d0) &
            &                   /4.0d0                                                                &
            &                   /Pi                                                                   &
            &                  )**(1.0d0/3.0d0)
    else
       call Error_Report('one of radiusEnvironment and massEnvironment must be specified'//{introspection:location})
    end if
    return
  end function fixedHEConstructorInternal

  subroutine fixedHEDestructor(self)
    !!{
    Destructor for the \refClass{haloEnvironmentFixed} halo mass function class.
    !!}
    implicit none
    type(haloEnvironmentFixed), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"     />
    <objectDestructor name="self%linearGrowth_"           />
    <objectDestructor name="self%sphericalCollapseSolver_"/>
    !!]
    return
  end subroutine fixedHEDestructor

  double precision function fixedHEOverdensityLinear(self,node,presentDay)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    implicit none
    class           (haloEnvironmentFixed              ), intent(inout)           :: self
    type            (treeNode                          ), intent(inout)           :: node
    logical                                             , intent(in   ), optional :: presentDay
    class           (nodeComponentBasic                ), pointer                 :: basic
    !![
    <optionalArgument name="presentDay" defaultsTo=".false." />
    !!]

    fixedHEOverdensityLinear=self%overdensity
    if (.not.presentDay_) then
       basic                    =>  node                                 %basic(                 )
       fixedHEOverdensityLinear = +fixedHEOverdensityLinear                                        &
            &                      *self                   %linearGrowth_%value(time=basic%time())
    end if
    return
  end function fixedHEOverdensityLinear

  double precision function fixedHEOverdensityLinearGradientTime(self,node)
    !!{
    Return the time gradient of the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentBasic  ), pointer       :: basic

    basic                                =>  node%basic()
    fixedHEOverdensityLinearGradientTime =  +self%overdensityLinear(node)                                                      &
         &                                  *self%linearGrowth_      %logarithmicDerivativeExpansionFactor( time=basic%time()) &
         &                                  *self%cosmologyFunctions_%expansionRate                       (                    &
         &                                   self%cosmologyFunctions_%expansionFactor                      (     basic%time()) &
         &                                                                                              )
    return
  end function fixedHEOverdensityLinearGradientTime

  double precision function fixedHEOverdensityNonLinear(self,node)
    !!{
    Return the environment of the given {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self
    type (treeNode            ), intent(inout) :: node
    class(nodeComponentBasic  ), pointer       :: basic

    ! Get a table of linear vs. nonlinear density.
    if (.not.self%linearToNonLinearInitialized) then
       call self%sphericalCollapseSolver_%linearNonlinearMap(self%cosmologyFunctions_%cosmicTime(1.0d0),self%linearToNonLinear)
       self%linearToNonLinearInitialized=.true.
    end if
    ! Find the nonlinear overdensity.
    basic                       => node                  %basic      (                                         )
    fixedHEOverdensityNonLinear =  self%linearToNonLinear%interpolate(self%overdensityLinear(node),basic%time())
    return
  end function fixedHEOverdensityNonLinear

  double precision function fixedHEEnvironmentRadius(self)
    !!{
    Return the radius of the environment.
    !!}
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self

    fixedHEEnvironmentRadius=self%radiusEnvironment
    return
  end function fixedHEEnvironmentRadius

  double precision function fixedHEEnvironmentMass(self)
    !!{
    Return the mass of the environment.
    !!}
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self

    fixedHEEnvironmentMass=self%massEnvironment
    return
  end function fixedHEEnvironmentMass

  double precision function fixedHEOverdensityLinearMaximum(self)
    !!{
    Return the maximum overdensity for which the \gls{pdf} is non-zero.
    !!}
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self

    fixedHEOverdensityLinearMaximum=self%overdensity
    return
  end function fixedHEOverdensityLinearMaximum

  double precision function fixedHEPDF(self,overdensity)
    !!{
    Return the PDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentFixed), intent(inout) :: self
    double precision                      , intent(in   ) :: overdensity
    !$GLC attributes unused :: self, overdensity
    
    fixedHEPDF=0.0d0
    call Error_Report('PDF is a delta-function'//{introspection:location})
    return
  end function fixedHEPDF

  double precision function fixedHECDF(self,overdensity)
    !!{
    Return the CDF of the environmental overdensity.
    !!}
    implicit none
    class           (haloEnvironmentFixed), intent(inout) :: self
    double precision                      , intent(in   ) :: overdensity
    !$GLC attributes unused :: self

    if (overdensity > self%overdensity) then
       fixedHECDF=1.0d0
    else
       fixedHECDF=0.0d0
    end if
    return
  end function fixedHECDF

  subroutine fixedHEOverdensityLinearSet(self,node,overdensity)
    !!{
    Set the overdensity of the environmental overdensity.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (haloEnvironmentFixed), intent(inout) :: self
    type            (treeNode            ), intent(inout) :: node
    double precision                      , intent(in   ) :: overdensity
    !$GLC attributes unused :: self, node, overdensity

    call Error_Report('can not set overdensity'//{introspection:location})
    return
  end subroutine fixedHEOverdensityLinearSet

  logical function fixedHEOverdensityIsSettable(self)
    !!{
    Return false as the overdensity is not settable.
    !!}
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedHEOverdensityIsSettable=.false.
    return
  end function fixedHEOverdensityIsSettable

  logical function fixedHEIsNodeDependent(self)
    !!{
    Return false as the environment is not dependent on the node.
    !!}
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedHEIsNodeDependent=.false.
    return
  end function fixedHEIsNodeDependent

  logical function fixedHEIsTreeDependent(self)
    !!{
    Return false as the environment is not dependent on the tree.
    !!}
    implicit none
    class(haloEnvironmentFixed), intent(inout) :: self
    !$GLC attributes unused :: self

    fixedHEIsTreeDependent=.false.
    return
  end function fixedHEIsTreeDependent
