!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

module Critical_Overdensities
  !% Provides an object that implements critical overdensities.
  use FGSL
  
  !# <functionClass>
  !#  <name>criticalOverdensity</name>
  !#  <descriptiveName>Critical Overdensity</descriptiveName>
  !#  <description>Object providing critical overdensities.</description>
  !#  <default>sphericalCollapseMatterLambda</default>
  !#  <stateful>yes</stateful>
  !#  <method name="value" >
  !#   <description>Return the critical overdensity at the given time and mass.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#   <argument>double precision, intent(in   ), optional :: mass</argument>
  !#  </method>
  !#  <method name="timeOfCollapse" >
  !#   <description>Returns the time of collapse for a perturbation of linear theory overdensity {\normalfont \ttfamily criticalOverdensity}.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   )           :: criticalOverdensity</argument>
  !#   <argument>double precision, intent(in   ), optional :: mass</argument>
  !#   <modules>Cosmology_Functions Root_Finder</modules>
  !#   <code>
  !#    double precision                         , parameter :: toleranceRelative   =1.0d-6, toleranceAbsolute=0.0d0
  !#    type            (rootFinder             ), save      :: finder
  !#    !$omp threadprivate(finder)
  !#    class           (cosmologyFunctionsClass), pointer   :: cosmologyFunctions_
  !#    if (.not.finder%isInitialized()) then
  !#       call finder%rootFunction(collapseTimeRoot                   )
  !#       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
  !#       call finder%rangeExpand (                                                                 &amp;
  !#            &amp;                   rangeExpandUpward            =2.0d0                        , &amp;
  !#            &amp;                   rangeExpandDownward          =0.5d0                        , &amp;
  !#            &amp;                   rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &amp;
  !#            &amp;                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &amp;
  !#            &amp;                   rangeExpandType              =rangeExpandMultiplicative      &amp;
  !#            &amp;                  )
  !#    end if
  !#    cosmologyFunctions_ => cosmologyFunctions()
  !#    criticalOverdensityTimeOfCollapse=finder%find(rootGuess=cosmologyFunctions_%cosmicTime(1.0d0))
  !#    return
  !#    contains
  !#      double precision function collapseTimeRoot(time)
  !#        double precision, intent(in   ) :: time
  !#        collapseTimeRoot=criticalOverdensity-self%value(time=time,mass=mass)
  !#        return
  !#      end function collapseTimeRoot
  !#   </code>
  !#  </method>
  !#  <method name="collapsingMass" >
  !#   <description>Return the mass scale just collapsing at the given cosmic time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#   <modules>Cosmology_Functions Root_Finder Power_Spectra</modules>
  !#   <code>
  !#    double precision                         , parameter :: massGuess           =1.0d13, toleranceAbsolute=0.0d0, &amp;
  !#         &amp;                                              toleranceRelative   =1.0d-6
  !#    type            (rootFinder             ), save      :: finder
  !#    !$omp threadprivate(finder)
  !#    class           (cosmologyFunctionsClass), pointer   :: cosmologyFunctions_
  !#    double precision                                     :: collapseTime
  !#    cosmologyFunctions_ => cosmologyFunctions()
  !#    call cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=collapseTime)
  !#    if (.not.finder%isInitialized()) then
  !#       call finder%rootFunction(collapsingMassRoot                 )
  !#       call finder%tolerance   (toleranceAbsolute,toleranceRelative)
  !#       call finder%rangeExpand (                                                                 &amp;
  !#            &amp;                   rangeExpandUpward            =2.0d0                        , &amp;
  !#            &amp;                   rangeExpandDownward          =0.5d0                        , &amp;
  !#            &amp;                   rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &amp;
  !#            &amp;                   rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &amp;
  !#            &amp;                   rangeExpandType              =rangeExpandMultiplicative      &amp;
  !#            &amp;                  )
  !#    end if
  !#    criticalOverdensityCollapsingMass=finder%find(rootGuess=massGuess)
  !#    return
  !#    contains
  !#      double precision function collapsingMassRoot(mass)
  !#        use Power_Spectra
  !#        implicit none
  !#        double precision, intent(in   ) :: mass        
  !#        collapsingMassRoot=Cosmological_Mass_Root_Variance(mass)-self%value(time=collapseTime,mass=mass)
  !#        return
  !#      end function collapsingMassRoot
  !#   </code>
  !#  </method>
  !#  <method name="gradientTime" >
  !#   <description>Return the derivative with respect to time of the linear theory critical overdensity for collapse at the given cosmic time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#   <argument>double precision, intent(in   ), optional :: mass</argument>
  !#  </method>
  !#  <method name="gradientMass" >
  !#   <description>Return the derivative with respect to mass of the linear theory critical overdensity for collapse at the given cosmic time.</description>
  !#   <type>double precision</type>
  !#   <pass>yes</pass>
  !#   <argument>double precision, intent(in   ), optional :: time      , expansionFactor</argument>
  !#   <argument>logical         , intent(in   ), optional :: collapsing</argument>
  !#   <argument>double precision, intent(in   ), optional :: mass</argument>
  !#  </method>
  !# </functionClass>

end module Critical_Overdensities
