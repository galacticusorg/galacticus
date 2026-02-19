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
Contains a module which provides a class that implements critical overdensity.
!!}

module Cosmological_Density_Field
  !!{
  Provides an object that implements critical overdensities and halo environments.
  !!}
  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Galacticus_Nodes   , only : treeNode
  use :: Linear_Growth      , only : linearGrowthClass
  use :: Tables             , only : table1DLinearLinear
  use :: Root_Finder        , only : rootFinder

  private

  !![
  <functionClass>
   <name>criticalOverdensity</name>
   <descriptiveName>Critical Overdensity</descriptiveName>
   <description>Object providing critical overdensities.</description>
   <default>sphericalCollapseClsnlssMttrCsmlgclCnstnt</default>
   <data>integer         (kind_int8                    )              :: lastUniqueID                =  -1_kind_int8, lastTreeID                 =-1_kind_int8</data>
   <data>double precision                                             :: criticalOverdensityTarget                  , mass                                    </data>
   <data>double precision                                             :: time                                       , timeNow                    =-huge(0.0d0)</data>
   <data>double precision                                             :: timeOfCollapsePrevious      =  -huge(0.0d0), criticalOverdensityPrevious=-huge(0.0d0)</data>
   <data>double precision                                             :: massPrevious                =  -huge(0.0d0)                                          </data>
   <data>type            (table1DLinearLinear          )              :: collapseThreshold                                                                    </data>
   <data>double precision                                             :: collapseThresholdMinimum                   , collapseThresholdMaximum                </data>
   <data>logical                                                      :: collapseThresholdInitialized=.false.                                                 </data>
   <data>type            (treeNode                     ), pointer     :: node                                                                                 </data>
   <data>logical                                                      :: massPresent                                , nodePresent                             </data>
   <data>logical                                                      :: treePresent                                                                          </data>
   <data>logical                                                      :: dependenciesInitialized     =  .false.     , isMassDependent_                        </data>
   <data>logical                                                      :: isNodeDependent_                           , isTreeDependent_                        </data>
   <data>class           (cosmologyFunctionsClass      ), pointer     :: cosmologyFunctions_         => null()                                                </data>
   <data>class           (linearGrowthClass            ), pointer     :: linearGrowth_               => null()                                                </data>
   <data>class           (cosmologicalMassVarianceClass), pointer     :: cosmologicalMassVariance_   => null()                                                </data>
   <data>type            (rootFinder                   ), allocatable :: finderTimeOfCollapse                                                                 </data>
   <data>
    <scope>module</scope>
    <threadprivate>yes</threadprivate>
    <content>class(criticalOverdensityClass), pointer :: globalSelf => null()</content>
   </data>
   <method name="value" >
    <description>Return the critical overdensity at the given time and mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical                   , intent(in   ), optional :: collapsing                 </argument>
    <argument>double precision          , intent(in   ), optional :: mass                       </argument>
    <argument>type            (treeNode), intent(inout), optional :: node                       </argument>
   </method>
   <method name="timeOfCollapse" >
    <description>Returns the time of collapse for a perturbation of linear theory overdensity {\normalfont \ttfamily criticalOverdensity}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   )                   :: criticalOverdensity</argument>
    <argument>double precision          , intent(in   ), optional         :: mass               </argument>
    <argument>type            (treeNode), intent(inout), optional, target :: node               </argument>
    <function>criticalOverdensityTimeOfCollapse</function>
   </method>
   <method name="collapsingMass" >
    <description>Return the mass scale just collapsing at the given cosmic time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <selfTarget>yes</selfTarget>
    <argument>double precision          , intent(in   ), optional         :: time      , expansionFactor</argument>
    <argument>logical                   , intent(in   ), optional         :: collapsing                 </argument>
    <argument>type            (treeNode), intent(inout), optional, target :: node                       </argument>
    <modules>Root_Finder Error</modules>
    <code>
     double precision            , parameter :: massGuess        =1.0d+13, toleranceAbsolute=0.0d+00, &amp;
          &amp;                                 toleranceRelative=1.0d-06, massTiny         =1.0d-30
     type            (rootFinder), save      :: finder
     logical                     , save      :: finderInitialized=.false.
     !$omp threadprivate(finder,finderInitialized)
     integer                                                    :: status
     double precision                                           :: collapseTime
     call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=collapseTime)
     if (.not.finderInitialized) then
        finder=rootFinder(                                                             &amp;
             &amp;        rootFunction                 =collapsingMassRoot           , &amp;
             &amp;        toleranceAbsolute            =toleranceAbsolute            , &amp;
             &amp;        toleranceRelative            =toleranceRelative            , &amp;
             &amp;        rangeExpandUpward            =2.0d0                        , &amp;
             &amp;        rangeExpandDownward          =0.5d0                        , &amp;
             &amp;        rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &amp;
             &amp;        rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &amp;
             &amp;        rangeDownwardLimit           =massTiny                     , &amp;
             &amp;        rangeExpandType              =rangeExpandMultiplicative      &amp;
             &amp;       )
        finderInitialized=.true.
     end if
     globalSelf                           => self
     self      %time                      =  collapseTime
     self      %nodePresent               =  present(node)
     if (self%nodePresent) self%node      =>         node
     criticalOverdensityCollapsingMass=finder%find(rootGuess=massGuess,status=status)
     if (status == errorStatusOutOfRange) criticalOverdensityCollapsingMass=0.0d0
     return
    </code>
   </method>
   <method name="gradientTime" >
    <description>Return the derivative with respect to time of the linear theory critical overdensity for collapse at the given cosmic time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical                   , intent(in   ), optional :: collapsing                 </argument>
    <argument>double precision          , intent(in   ), optional :: mass                       </argument>
    <argument>type            (treeNode), intent(inout), optional :: node                       </argument>
   </method>
   <method name="gradientMass" >
    <description>Return the derivative with respect to mass of the linear theory critical overdensity for collapse at the given cosmic time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision          , intent(in   ), optional :: time      , expansionFactor</argument>
    <argument>logical                   , intent(in   ), optional :: collapsing                 </argument>
    <argument>double precision          , intent(in   ), optional :: mass                       </argument>
    <argument>type            (treeNode), intent(inout), optional :: node                       </argument>
   </method>
   <method name="isMassDependent" >
    <description>Return true if the critical overdensity is dependent on the mass of the halo.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="isNodeDependent" >
    <description>Return true if the critical overdensity is dependent on the {\normalfont \ttfamily node} object (not just the {\normalfont \ttfamily node\%hostTree} object).</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <method name="isTreeDependent" >
    <description>Return true if the critical overdensity is dependent on the {\normalfont \ttfamily node\%hostTree} object.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
   <autoHook>
    <modules>
     <name>Events_Hooks</name>
     <only>calculationResetEvent, openMPThreadBindingAllLevels</only>
    </modules>
    <code>
     call calculationResetEvent%attach(self,criticalOverdensityCalculationReset,openMPThreadBindingAllLevels,label='criticalOverdensity')
     return
    </code>
   </autoHook>
   <destructor>
    <modules>
     <name>Events_Hooks</name>
     <only>calculationResetEvent, openMPThreadBindingAllLevels</only>
    </modules>
    <code>
     if (calculationResetEvent%isAttached(self,criticalOverdensityCalculationReset)) call calculationResetEvent%detach(self,criticalOverdensityCalculationReset)
     return
    </code>
   </destructor>
  </functionClass>
  !!]

  !![
  <functionClass>
   <name>haloEnvironment</name>
   <descriptiveName>Halo Environment</descriptiveName>
   <description>Class providing halo environment.</description>
   <default>uniform</default>
   <method name="overdensityLinear" >
    <description>Return the environmental linear overdensity for the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)           :: node</argument>
    <argument>logical          , intent(in   ), optional :: presentDay</argument>
   </method>
   <method name="overdensityLinearGradientTime" >
    <description>Return the gradient with time of the environmental linear overdensity for the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type   (treeNode), intent(inout)           :: node</argument>
   </method>
   <method name="overdensityNonLinear" >
    <description>Return the environmental non-linear overdensity for the given {\normalfont \ttfamily node}.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>type(treeNode), intent(inout) :: node</argument>
   </method>
   <method name="environmentRadius" >
    <description>Return the radius of the region used to defined the environment.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="environmentMass" >
    <description>Return the mean mass contained in the region used to defined the environment.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>  
   <method name="overdensityLinearMinimum" >
    <description>Return the minimum linear overdensity for which the environmental overdensity \gls{pdf} is non-zero.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      haloEnvironmentOverdensityLinearMinimum=-huge(0.0d0)
    </code>
   </method>
   <method name="overdensityLinearMaximum" >
    <description>Return the maximum linear overdensity for which the environmental overdensity \gls{pdf} is non-zero.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      haloEnvironmentOverdensityLinearMaximum=+huge(0.0d0)
    </code>
   </method>
   <method name="pdf" >
    <description>Return the \gls{pdf} of the environmental overdensity for the given overdensity.</description>
    <type>double precision</type>
    <argument>double precision, intent(in   ) :: overdensity</argument>
    <pass>yes</pass>
   </method>
   <method name="cdf" >
    <description>Return the \gls{cdf} of the environmental overdensity for the given overdensity.</description>
    <type>double precision</type>
    <argument>double precision, intent(in   ) :: overdensity</argument>
    <pass>yes</pass>
   </method>
   <method name="overdensityLinearSet" >
    <description>Set the environmental overdensity for the give node.</description>
    <type>void</type>
    <argument>type            (treeNode), intent(inout) :: node</argument>
    <argument>double precision          , intent(in   ) :: overdensity</argument>
    <pass>yes</pass>
   </method>
   <method name="overdensityIsSettable" >
    <description>Return true if the overdensity is settable.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      haloEnvironmentOverdensityIsSettable=.true.
    </code>
   </method>
   <method name="volumeFractionOccupied" >
    <description>Return the fraction of the volume occupied by the regions described by this environment.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      haloEnvironmentVolumeFractionOccupied=1.0d0
    </code>
   </method>
   <method name="isNodeDependent" >
    <description>Return true if the environment is node dependent (but false if the only dependency on the node is via its host tree).</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      haloEnvironmentIsNodeDependent=.true.
    </code>
   </method>
   <method name="isTreeDependent" >
    <description>Return true if the environment is tree dependent.</description>
    <type>logical</type>
    <pass>yes</pass>
    <code>
      !$GLC attributes unused :: self
      haloEnvironmentIsTreeDependent=.true.
    </code>
   </method>
  </functionClass>
  !!]

  !![
  <functionClass>
   <name>cosmologicalMassVariance</name>
   <descriptiveName>Mass Variance of Cosmological Density Field</descriptiveName>
   <description>
    A class providing the mass variance, $\sigma(M)$, of the cosmological density field.
   </description>
   <default>filteredPower</default>
   <method name="powerNormalization" >
    <description>Return the normalization of the power spectrum.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="sigma8" >
    <description>Return the value of $\sigma_8$.</description>
    <type>double precision</type>
    <pass>yes</pass>
   </method>
   <method name="rootVariance" >
    <description>Return the root-variance of the cosmological density field.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
   <method name="rootVarianceLogarithmicGradient" >
    <description>Return the logarithmic gradient of the root-variance of the cosmological density field with respect to mass.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
   <method name="rootVarianceLogarithmicGradientTime" >
    <description>Return the logarithmic gradient of the root-variance of the cosmological density field with respect to time.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass, time</argument>
   </method>
   <method name="rootVarianceAndLogarithmicGradient" >
    <description>Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field.</description>
    <type>void</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: mass        , time</argument>
    <argument>double precision, intent(  out) :: rootVariance, rootVarianceLogarithmicGradient</argument>
   </method>
   <method name="mass" >
    <description>Return the mass corresponding to the given {\normalfont \ttfamily rootVariance} of the cosmological density field.</description>
    <type>double precision</type>
    <pass>yes</pass>
    <argument>double precision, intent(in   ) :: rootVariance, time</argument>
   </method>
   <method name="growthIsMassDependent" >
    <description>Return true if the growth of the variance with time is mass-dependent.</description>
    <type>logical</type>
    <pass>yes</pass>
   </method>
  </functionClass>
  !!]

contains

  double precision function criticalOverdensityTimeOfCollapse(self,criticalOverdensity,mass,node,status)
    !!{
    Returns the time of collapse for a perturbation of linear theory overdensity {\normalfont \ttfamily criticalOverdensity}.
    !!}
    use :: Cosmology_Functions, only : timeToleranceRelativeBigCrunch
    use :: Root_Finder        , only : rangeExpandMultiplicative     , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive
    use :: Error              , only : Error_Report                  , errorStatusSuccess           , errorStatusOutOfRange
    implicit none
    class           (criticalOverdensityClass), intent(inout)              , target :: self
    double precision                          , intent(in   )                       :: criticalOverdensity
    double precision                          , intent(in   ), optional             :: mass
    type            (treeNode                ), intent(inout), optional    , target :: node
    integer                                   , intent(  out), optional             :: status
    double precision                          , parameter                           :: toleranceRelative                =1.0d-12, toleranceAbsolute       =0.0d0, &
         &                                                                             fractionTimeCollapseGrowthMinimum=1.0d-03
    integer                                   , parameter                           :: countPerUnit                     =10000
    double precision                          , allocatable  , dimension(:)         :: threshold
    double precision                                                                :: timeBigCrunch                            , timeGuess                     , &
         &                                                                             collapseThresholdMinimum                 , collapseThresholdMaximum      , &
         &                                                                             collapseTimePrevious                     , collapseTimeUpperLimit        , &
         &                                                                             timeUpperLimit
    logical                                                                         :: updateResult                             , remakeTable
    integer                                                                         :: i                                        , countThresholds               , &
         &                                                                             countNewLower                            , countNewUpper

    ! Assume a successful calculation by default.
    if (present(status)) status=errorStatusSuccess    
    ! Determine dependencies.
    if (.not.self%dependenciesInitialized) then
       self%isMassDependent_       =self%isMassDependent() .or. self%cosmologicalMassVariance_%growthIsMassDependent()
       self%isNodeDependent_       =self%isNodeDependent()
       self%isTreeDependent_       =self%isTreeDependent()
       self%dependenciesInitialized=.true.
    end if
    ! Determine which arguments are present and not ignorable.
    self%massPresent=present(mass).and.self%isMassDependent_
    self%nodePresent=present(node).and.self%isNodeDependent_
    self%treePresent=present(node).and.self%isTreeDependent_
    ! Determine if the memoized value must be updated.
    updateResult=.false.
    if (criticalOverdensity /= self%criticalOverdensityPrevious) updateResult=.true.
    if (self%massPresent) then
       if (mass                              /= self%massPrevious) updateResult=.true.
    else
       if (-huge(0.0d0)                      /= self%massPrevious) updateResult=.true.
    end if
    if (self%nodePresent) then
       if (node                  %uniqueID() /= self%lastUniqueID) updateResult=.true.
    else
       if (-1_kind_int8                      /= self%lastUniqueID) updateResult=.true.
    end if
    if (self%treePresent) then
       if (node%hostTree%nodeBase%uniqueID() /= self%lastTreeID  ) updateResult=.true.
    else
       if (-1_kind_int8                      /= self%lastTreeID  ) updateResult=.true.
    end if
     ! Recompute memoized value if necessary.
    if (updateResult) then
       if (self%massPresent) then
          self%massPrevious=mass
       else
          self%massPrevious=-huge(0.0d0)
       end if
       if (self%nodePresent) then
          self%lastUniqueID=node         %uniqueID()
       else
          self%lastUniqueID=-1_kind_int8
       end if
       self%criticalOverdensityPrevious=criticalOverdensity
       if (self%timeNow < 0.0d0) self%timeNow=self%cosmologyFunctions_%cosmicTime(1.0d0)
       timeBigCrunch=self%cosmologyFunctions_%timeBigCrunch()
       if (timeBigCrunch < 0.0d0) then
          timeBigCrunch=huge(0.0d0)
       else
          timeBigCrunch=timeBigCrunch*(1.0d0-timeToleranceRelativeBigCrunch)
          ! Check for critical overdensity for collapse exceeding that at the Big Crunch.
          if     (                                                                                           &
               &   +criticalOverdensity                                                                      &
               &  <=                                                                                         &
               &   +self                          %value       (time=     timeBigCrunch,mass=mass,node=node) &
               &   *self%cosmologicalMassVariance_%rootVariance(time=self%timeNow      ,mass=mass          ) &
               &   /self%cosmologicalMassVariance_%rootVariance(time=     timeBigCrunch,mass=mass          ) &
               & ) then
             ! Return a collapse time close to the Big Crunch in this case.
             criticalOverdensityTimeOfCollapse=timeBigCrunch
             return
          end if
       end if
       if (.not.allocated(self%finderTimeOfCollapse)) then
          allocate(self%finderTimeOfCollapse)
          self%finderTimeOfCollapse=rootFinder(rootFunction=collapseTimeRoot,toleranceAbsolute=toleranceAbsolute,toleranceRelative=toleranceRelative)
          call self%finderTimeOfCollapse%rangeExpand(                                                             &
               &                                     rangeExpandUpward            =2.0d0                        , &
               &                                     rangeExpandDownward          =0.5d0                        , &
               &                                     rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                                     rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                                     rangeUpwardLimit             =timeBigCrunch                , &
               &                                     rangeExpandType              =rangeExpandMultiplicative      &
               &                                    )
       end if
       globalSelf => self
       if     (                   &
            &   self%massPresent  &
            & ) self%mass =  mass
       if     (                   &
            &   self%nodePresent  &
            &  .or.               &
            &   self%treePresent  &
            & ) self%node => node
       if (self%massPresent) then
          self%criticalOverdensityTarget=criticalOverdensity/self%cosmologicalMassVariance_%rootVariance(mass,self%timeNow)
       else
          self%criticalOverdensityTarget=criticalOverdensity
       end if
       if (self%timeOfCollapsePrevious < 0.0d0) self%timeOfCollapsePrevious=self%cosmologyFunctions_%cosmicTime(1.0d0)
       if (.not.self%massPresent.and..not.self%nodePresent) then
          ! There is no dependency on mass, or the node. If there is dependency on the tree, then, if the tree differs from the
          ! previously seen case we must destroy the previously made table, as our tree has changed.
          if (self%treePresent) then
             if (node%hostTree%nodeBase%uniqueID() /= self%lastTreeID) then
                self%collapseThresholdInitialized=.false.
                call self%collapseThreshold%destroy()
             end if
             self%lastTreeID  =node%hostTree%nodeBase%uniqueID()
          else
             self%lastTreeID  =-1_kind_int8
          end if
          ! Neither the mass or the node are provided, so we can use a simple tabulation of collapse thresholds for rapid
          ! inversion. Note that we do not tabulate lower than the requested threshold as in some cosmologies (e.g. with dark
          ! energy or a cosmological constant) this can require tabulating to extremely large cosmic times).
          remakeTable=.false.
          if (.not.self%collapseThresholdInitialized) then
             remakeTable                  =.true.
             self%collapseThresholdMinimum=      criticalOverdensity
             self%collapseThresholdMaximum=2.0d0*criticalOverdensity
             countThresholds              =int(dble(countPerUnit)*(self%collapseThresholdMaximum-self%collapseThresholdMinimum))+2
             ! Ensure the maximum of the table is precisely an integer number of steps above the minimum.
             self%collapseThresholdMaximum=self%collapseThresholdMinimum+dble(countThresholds-1)/dble(countPerUnit)
             allocate(threshold(countThresholds))
             threshold=-huge(0.0d0)
          else if (criticalOverdensity < self%collapseThresholdMinimum .or. criticalOverdensity > self%collapseThresholdMaximum) then
             remakeTable                      =.true.
             collapseThresholdMinimum=min(      criticalOverdensity,self%collapseThresholdMinimum)
             collapseThresholdMaximum=max(2.0d0*criticalOverdensity,self%collapseThresholdMaximum)
             ! Determine how many points the table must be extended by in each direction to span the new required range.
             countNewLower=0
             countNewUpper=0
             if (self%collapseThresholdMinimum > collapseThresholdMinimum) countNewLower=int((+self%collapseThresholdMinimum-collapseThresholdMinimum)*dble(countPerUnit)+1.0d0)
             if (self%collapseThresholdMaximum < collapseThresholdMaximum) countNewUpper=int((-self%collapseThresholdMaximum+collapseThresholdMaximum)*dble(countPerUnit)+1.0d0)
             countThresholds=self%collapseThreshold%size()+countNewLower+countNewUpper
             ! Adjust the limits of the table by an integer number of steps.
             self%collapseThresholdMinimum=self%collapseThresholdMinimum-dble(countNewLower)/dble(countPerUnit)
             self%collapseThresholdMaximum=self%collapseThresholdMaximum+dble(countNewUpper)/dble(countPerUnit)
             allocate(threshold(countThresholds))
             threshold=-huge(0.0d0)
             ! Populate the table with pre-existing results.
             threshold(countNewLower+1:countNewLower+self%collapseThreshold%size())=reshape(self%collapseThreshold%ys(),[self%collapseThreshold%size()])
          end if
          if (remakeTable) then
             if (self%collapseThresholdInitialized) call self%collapseThreshold%destroy()
             call self%collapseThreshold%create(self%collapseThresholdMinimum,self%collapseThresholdMaximum,countThresholds)
             ! Populate the table in regions where it was not previously populated.
             do i=1,countThresholds
                self%criticalOverdensityTarget=self%collapseThreshold%x(i)
                if (threshold(i) < 0.0d0) then
                   if (i == 1) then
                      timeGuess=self%cosmologyFunctions_%cosmicTime(expansionFactor=1.0d0)
                   else
                      timeGuess=threshold(i-1)
                   end if
                   threshold(i)=self%finderTimeOfCollapse%find(rootGuess=timeGuess)
                end if
             end do
             call self%collapseThreshold%populate(threshold)
             deallocate(threshold)
             self%criticalOverdensityTarget   =criticalOverdensity
             self%collapseThresholdInitialized=.true.
          end if
          self%timeOfCollapsePrevious=self%collapseThreshold%interpolate(criticalOverdensity)
       else
          ! Check that we can find an upper bound to the collapse time. In some cosmologies (e.g. those with a cosmological
          ! constant) there is a limit to how much the linear growth factor can increase - and therefore some critical
          ! overdensities will never be reached.
          timeUpperLimit        =+self%cosmologyFunctions_%cosmicTime               (expansionFactor=1.0d0         )
          collapseTimeUpperLimit=+self                    %criticalOverdensityTarget                                 &
               &                 -                         collapseTimeRoot         (                timeUpperLimit)
          do while (collapseTimeRoot(timeUpperLimit) < 0.0d0)
             timeUpperLimit        =+2.0d0                                          &
                  &                 *                               timeUpperLimit
             collapseTimePrevious  =+     collapseTimeUpperLimit
             collapseTimeUpperLimit=+self%criticalOverdensityTarget                 &
                  &                 -     collapseTimeRoot         (timeUpperLimit)
             ! Check if the collapsing density has increased by some minimal amount as we doubled the time of collapse.
             if     (                                                                                                      &
                  &   collapseTimePrevious-collapseTimeUpperLimit < fractionTimeCollapseGrowthMinimum*collapseTimePrevious &
                  &  .and.                                                                                                 &
                  &   collapseTimeRoot(timeUpperLimit)            < 0.0d0                                                  &
                  & ) then
                ! It did not, suggesting that there is no solution for the collapse time.
                if (present(status)) then
                   criticalOverdensityTimeOfCollapse=-huge(0.0d0)
                   status                           =errorStatusOutOfRange
                   return
                else
                   call Error_Report('unable to bracket collapse time - this can happen in cosmologies with a upper bound to the linear growth (e.g. cosmological constant models)'//{introspection:location})
                end if
             end if
          end do
          self%timeOfCollapsePrevious=self%finderTimeOfCollapse%find(rootGuess=self%timeOfCollapsePrevious)
       end if
    end if
    ! Return the memoized value.
    criticalOverdensityTimeOfCollapse=self%timeOfCollapsePrevious
    return
  end function criticalOverdensityTimeOfCollapse

  double precision function collapseTimeRoot(time)
    !!{
    Function used in root finding for the collapse time at a given critical overdensity. We have some target linear theory
    overdensity, $\delta_0$, extrapolated to the present epoch, $t_0$, and want to know at what epoch is would collapse. In a
    pure dark matter universe, in which modes on all scales grow at the same rate (i.e. the growth factor $D(k,t)$ is
    independent of wavenumber $k$) we would simply have $\delta(t) = \delta_0 D(t)$, and so would solve for
    $\delta_\mathrm{c}(t)/D(t) = \delta_0$. We want to generalize this to cases where the growth factor is scale dependent. In
    these cases we actually want to consider the growth rate for a region of fixed mass. Since modes of a range of different
    wavenumber contribute to the growth of such a region, we write $\delta(t) = \delta_0 \sigma(M,t) / \sigma(M,t_0)$, and
    therefore solve for $\delta_\mathrm{c}(t) \sigma(M,t_0) / \sigma(M,t) = \delta_0$. If no mass argument was provided then we
    revert to assuming that $D(t,k)$ is independent of wavenumber and simply solve for $\delta_\mathrm{c}(t)/D(t) = \delta_0$.
    !!}
    implicit none
    double precision, intent(in   ) :: time

    if (globalSelf%massPresent) then
       if (globalSelf%nodePresent) then
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf                          %value       (time=time,mass=globalSelf%mass,node=globalSelf%node) &
               &                                               /globalSelf%cosmologicalMassVariance_%rootVariance(time=time,mass=globalSelf%mass                     )
       else
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf                          %value       (time=time,mass=globalSelf%mass                     ) &
               &                                               /globalSelf%cosmologicalMassVariance_%rootVariance(time=time,mass=globalSelf%mass                     )
       end if
    else
       if (globalSelf%nodePresent) then
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf                          %value       (time=time                     ,node=globalSelf%node) &
               &                                               /globalSelf%linearGrowth_            %value       (time=time                                          )
       else
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf                          %value       (time=time                     ,node=globalSelf%node) &
               &                                               /globalSelf%linearGrowth_            %value       (time=time                                          )
       end if
    end if
    return
  end function collapseTimeRoot

  double precision function collapsingMassRoot(mass)
    !!{
    Function used in root finding for the collapsing mass at a given time.
    !!}
    implicit none
    double precision, intent(in   ) :: mass

    if (globalSelf%nodePresent) then
       collapsingMassRoot=globalSelf%cosmologicalMassVariance_%rootVariance(mass,globalSelf%time)-globalSelf%value(time=globalSelf%time,mass=mass,node=globalSelf%node)
    else
       collapsingMassRoot=globalSelf%cosmologicalMassVariance_%rootVariance(mass,globalSelf%time)-globalSelf%value(time=globalSelf%time,mass=mass                     )
    end if
    return
  end function collapsingMassRoot

  subroutine criticalOverdensityCalculationReset(self,node,uniqueID)
    !!{
    Reset the critical overdensity calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class(  criticalOverdensityClass), intent(inout) :: self
    type   (treeNode                ), intent(inout) :: node
    integer(kind_int8               ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%lastUniqueID=uniqueID
    self%lastTreeID  =-1_kind_int8
    if (associated(node%hostTree)) then
       if (associated(node%hostTree%nodeBase)) self%lastTreeID=node%hostTree%nodeBase%uniqueID()
    end if
    return
  end subroutine criticalOverdensityCalculationReset

end module Cosmological_Density_Field
