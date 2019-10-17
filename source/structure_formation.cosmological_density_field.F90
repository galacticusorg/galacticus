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

!% Contains a module which provides a class that implements critical overdensity.

module Cosmological_Density_Field
  !% Provides an object that implements critical overdensities and halo environments.
  use :: Cosmology_Functions, only : cosmologyFunctionsClass, timeToleranceRelativeBigCrunch
  use :: Galacticus_Nodes   , only : treeNode
  private

  !# <functionClass>
  !#  <name>criticalOverdensity</name>
  !#  <descriptiveName>Critical Overdensity</descriptiveName>
  !#  <description>Object providing critical overdensities.</description>
  !#  <default>sphericalCollapseMatterLambda</default>
  !#  <data>double precision                                         :: criticalOverdensityTarget          , mass       , time</data>
  !#  <data>type            (treeNode                     ), pointer :: node                                                  </data>
  !#  <data>logical                                                  :: massPresent                        , nodePresent      </data>
  !#  <data>class           (cosmologyFunctionsClass      ), pointer :: cosmologyFunctions_       => null()                   </data>
  !#  <data>class           (cosmologicalMassVarianceClass), pointer :: cosmologicalMassVariance_ => null()                   </data>
  !#  <data>
  !#   <scope>module</scope>
  !#   <threadprivate>yes</threadprivate>
  !#   <content>class(criticalOverdensityClass), pointer :: globalSelf</content>
  !#  </data>
  !#  <method name="value" >
  !#   <description>Return the critical overdensity at the given time and mass.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision          , intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical                   , intent(in   ), optional :: collapsing                 </argument>
  !#   <argument>double precision          , intent(in   ), optional :: mass                       </argument>
  !#   <argument>type            (treeNode), intent(inout), optional :: node                       </argument>
  !#  </method>
  !#  <method name="timeOfCollapse" >
  !#   <description>Returns the time of collapse for a perturbation of linear theory overdensity {\normalfont \ttfamily criticalOverdensity}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>double precision          , intent(in   )                   :: criticalOverdensity</argument>
  !#   <argument>double precision          , intent(in   ), optional         :: mass               </argument>
  !#   <argument>type            (treeNode), intent(inout), optional, target :: node               </argument>
  !#   <modules>Root_Finder</modules>
  !#   <code>
  !#    double precision            , parameter :: toleranceRelative=1.0d-12, toleranceAbsolute=0.0d0
  !#    double precision                        :: timeBigCrunch
  !#    type            (rootFinder)            :: finder
  !#    call finder%rootFunction(collapseTimeRoot                   )
  !#    call finder%tolerance   (toleranceAbsolute,toleranceRelative)
  !#    timeBigCrunch=self%cosmologyFunctions_%timeBigCrunch()
  !#    if (timeBigCrunch &lt; 0.0d0) then
  !#       timeBigCrunch=huge(0.0d0)
  !#    else
  !#       timeBigCrunch=timeBigCrunch*(1.0d0-timeToleranceRelativeBigCrunch)
  !#       ! Check for critical overdensity for collapse exceeding that at the Big Crunch.
  !#       if (criticalOverdensity &lt;= self%value(time=timeBigCrunch,mass=mass,node=node)) then
  !#          ! Return a collapse time close to the Big Crunch in this case.
  !#          criticalOverdensityTimeOfCollapse=timeBigCrunch
  !#          return
  !#       end if
  !#    end if
  !#    call finder%rangeExpand(                                                             &amp;
  !#       &amp;                rangeExpandUpward            =2.0d0                        , &amp;
  !#       &amp;                rangeExpandDownward          =0.5d0                        , &amp;
  !#       &amp;                rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
  !#       &amp;                rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &amp;
  !#       &amp;                rangeUpwardLimit             =timeBigCrunch                , &amp;
  !#       &amp;                rangeExpandType              =rangeExpandMultiplicative      &amp;
  !#       &amp;               )
  !#    globalSelf                           => self
  !#    self      %criticalOverdensityTarget =  criticalOverdensity
  !#    self      %massPresent               =  present(mass)
  !#    self      %nodePresent               =  present(node)
  !#    if (self%massPresent) self%mass      =          mass
  !#    if (self%nodePresent) self%node      =>         node
  !#    criticalOverdensityTimeOfCollapse=finder%find(rootGuess=self%cosmologyFunctions_%cosmicTime(1.0d0))
  !#    return
  !#   </code>
  !#  </method>
  !#  <method name="collapsingMass" >
  !#   <description>Return the mass scale just collapsing at the given cosmic time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <selfTarget>yes</selfTarget>
  !#   <argument>double precision          , intent(in   ), optional         :: time      , expansionFactor</argument>
  !#   <argument>logical                   , intent(in   ), optional         :: collapsing                 </argument>
  !#   <argument>type            (treeNode), intent(inout), optional, target :: node                       </argument>
  !#   <modules>Root_Finder Galacticus_Error</modules>
  !#   <code>
  !#    double precision                               , parameter :: massGuess           =1.0d+13, toleranceAbsolute=0.0d+00, &amp;
  !#         &amp;                                                    toleranceRelative   =1.0d-06, massTiny         =1.0d-30
  !#    type            (rootFinder                   ), save      :: finder
  !#    !$omp threadprivate(finder)
  !#    integer                                                    :: status
  !#    double precision                                           :: collapseTime
  !#    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=collapseTime)
  !#    if (.not.finder%isInitialized()) then
  !#       call finder%rootFunction(collapsingMassRoot                 )
  !#       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
  !#       call finder%rangeExpand (                                                                 &amp;
  !#            &amp;                   rangeExpandUpward            =2.0d0                        , &amp;
  !#            &amp;                   rangeExpandDownward          =0.5d0                        , &amp;
  !#            &amp;                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &amp;
  !#            &amp;                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &amp;
  !#            &amp;                   rangeDownwardLimit           =massTiny                     , &amp;
  !#            &amp;                   rangeExpandType              =rangeExpandMultiplicative      &amp;
  !#            &amp;                  )
  !#    end if
  !#    globalSelf                           => self
  !#    self      %time                      =  collapseTime
  !#    self      %nodePresent               =  present(node)
  !#    if (self%nodePresent) self%node      =>         node
  !#    criticalOverdensityCollapsingMass=finder%find(rootGuess=massGuess,status=status)
  !#    if (status == errorStatusOutOfRange) criticalOverdensityCollapsingMass=0.0d0
  !#    return
  !#   </code>
  !#  </method>
  !#  <method name="gradientTime" >
  !#   <description>Return the derivative with respect to time of the linear theory critical overdensity for collapse at the given cosmic time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision          , intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical                   , intent(in   ), optional :: collapsing                 </argument>
  !#   <argument>double precision          , intent(in   ), optional :: mass                       </argument>
  !#   <argument>type            (treeNode), intent(inout), optional :: node                       </argument>
  !#  </method>
  !#  <method name="gradientMass" >
  !#   <description>Return the derivative with respect to mass of the linear theory critical overdensity for collapse at the given cosmic time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision          , intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical                   , intent(in   ), optional :: collapsing                 </argument>
  !#   <argument>double precision          , intent(in   ), optional :: mass                       </argument>
  !#   <argument>type            (treeNode), intent(inout), optional :: node                       </argument>
  !#  </method>
  !#  <method name="isMassDependent" >
  !#   <description>Return true if the critical overdensity is dependent on the mass of the halo.</description>
  !#   <type>logical</type>
  !#   <pass>yes</pass>
  !#  </method>
  !# </functionClass>

  !# <functionClass>
  !#  <name>haloEnvironment</name>
  !#  <descriptiveName>Halo Environment</descriptiveName>
  !#  <description>Class providing halo environment.</description>
  !#  <default>uniform</default>
  !#  <method name="overdensityLinear" >
  !#   <description>Return the environmental linear overdensity for the given {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type   (treeNode), intent(inout)           :: node</argument>
  !#   <argument>logical          , intent(in   ), optional :: presentDay</argument>
  !#  </method>
  !#  <method name="overdensityNonLinear" >
  !#   <description>Return the environmental non-linear overdensity for the given {\normalfont \ttfamily node}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>type(treeNode), intent(inout) :: node</argument>
  !#  </method>
  !#  <method name="environmentRadius" >
  !#   <description>Return the radius of the region used to defined the environmental.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="environmentMass" >
  !#   <description>Return the mean mass contained in the region used to defined the environmental.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="overdensityLinearMinimum" >
  !#   <description>Return the minimum linear overdensity for which the environmental overdensity \gls{pdf} is non-zero.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <code>
  !#     !GCC$ attributes unused :: self
  !#     haloEnvironmentOverdensityLinearMinimum=-huge(0.0d0)
  !#   </code>
  !#  </method>
  !#  <method name="overdensityLinearMaximum" >
  !#   <description>Return the maximum linear overdensity for which the environmental overdensity \gls{pdf} is non-zero.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <code>
  !#     !GCC$ attributes unused :: self
  !#     haloEnvironmentOverdensityLinearMaximum=+huge(0.0d0)
  !#   </code>
  !#  </method>
  !#  <method name="pdf" >
  !#   <description>Return the \gls{pdf} of the environmental overdensity for the given overdensity.</description>
  !#   <type>double precision</type>
  !#   <argument>double precision, intent(in   ) :: overdensity</argument>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="cdf" >
  !#   <description>Return the \gls{cdf} of the environmental overdensity for the given overdensity.</description>
  !#   <type>double precision</type>
  !#   <argument>double precision, intent(in   ) :: overdensity</argument>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="overdensityLinearSet" >
  !#   <description>Set the environmental overdensity for the give node.</description>
  !#   <type>void</type>
  !#   <argument>type            (treeNode), intent(inout) :: node</argument>
  !#   <argument>double precision          , intent(in   ) :: overdensity</argument>
  !#   <pass>yes</pass>
  !#  </method>
  !# </functionClass>

  !# <functionClass>
  !#  <name>cosmologicalMassVariance</name>
  !#  <descriptiveName>Mass Variance of Cosmological Density Field</descriptiveName>
  !#  <description>Object providing mass variance of the cosmological density field.</description>
  !#  <default>filteredPower</default>
  !#  <method name="powerNormalization" >
  !#   <description>Return the normalization of the power spectrum.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="sigma8" >
  !#   <description>Return the value of $\sigma_8$.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#  </method>
  !#  <method name="rootVariance" >
  !#   <description>Return the root-variance of the cosmological density field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#  </method>
  !#  <method name="rootVarianceLogarithmicGradient" >
  !#   <description>Return the logarithmic gradient of the root-variance of the cosmological density field with respect to mass.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#  </method>
  !#  <method name="rootVarianceAndLogarithmicGradient" >
  !#   <description>Return the value and logarithmic gradient with respect to mass of the root-variance of the cosmological density field.</description>
  !#   <type>void</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: mass</argument>
  !#   <argument>double precision, intent(  out) :: rootVariance, rootVarianceLogarithmicGradient</argument>
  !#  </method>
  !#  <method name="mass" >
  !#   <description>Return the mass corresponding to the given {\normalfont \ttfamily rootVariance} of the cosmological density field.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ) :: rootVariance</argument>
  !#  </method>
  !# </functionClass>

contains

  double precision function collapseTimeRoot(time)
    !% Function used in root finding for the collapse time at a given critical overdensity.
    implicit none
    double precision, intent(in   ) :: time

    if (globalSelf%massPresent) then
       if (globalSelf%nodePresent) then
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf%value(time=time,mass=globalSelf%mass,node=globalSelf%node)
       else
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf%value(time=time,mass=globalSelf%mass                     )
       end if
    else
       if (globalSelf%nodePresent) then
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf%value(time=time                     ,node=globalSelf%node)
       else
          collapseTimeRoot=globalSelf%criticalOverdensityTarget-globalSelf%value(time=time                                          )
       end if
    end if
    return
  end function collapseTimeRoot

  double precision function collapsingMassRoot(mass)
    !% Function used in root finding for the collapsing mass at a given time.
    implicit none
    double precision, intent(in   ) :: mass

    if (globalSelf%nodePresent) then
       collapsingMassRoot=globalSelf%cosmologicalMassVariance_%rootVariance(mass)-globalSelf%value(time=globalSelf%time,mass=mass,node=globalSelf%node)
    else
       collapsingMassRoot=globalSelf%cosmologicalMassVariance_%rootVariance(mass)-globalSelf%value(time=globalSelf%time,mass=mass                     )
    end if
    return
  end function collapsingMassRoot

end module Cosmological_Density_Field
