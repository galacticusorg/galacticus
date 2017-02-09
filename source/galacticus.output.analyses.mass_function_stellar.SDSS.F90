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

!% Contains a module which implements a stellar mass function output analysis class.

  !# <outputAnalysis name="outputAnalysisMassFunctionStellarSDSS" defaultThreadPrivate="yes">
  !#  <description>An SDSS stellar mass function output analysis class.</description>
  !# </outputAnalysis>
  type, extends(outputAnalysisMassFunctionStellar) :: outputAnalysisMassFunctionStellarSDSS
     !% An SDSS stellar mass function output analysis class.
     private
  end type outputAnalysisMassFunctionStellarSDSS

  interface outputAnalysisMassFunctionStellarSDSS
     !% Constructors for the ``massFunctionStellarSDSS'' output analysis class.
     module procedure massFunctionStellarSDSSConstructorParameters
     module procedure massFunctionStellarSDSSConstructorInternal
  end interface outputAnalysisMassFunctionStellarSDSS

contains

  function massFunctionStellarSDSSConstructorParameters(parameters)
    !% Constructor for the ``massFunctionStellarSDSS'' output analysis class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type            (outputAnalysisMassFunctionStellarSDSS)                :: massFunctionStellarSDSSConstructorParameters
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass              ), pointer       :: cosmologyFunctions_
    double precision                                                       :: randomError
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>randomError</name>
    !#   <source>parameters</source>
    !#   <variable>randomError</variable>
    !#   <defaultValue>0.07d0</defaultValue>
    !#   <description>The random error for stellar masses [in dex].</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    ! Build the object.
    massFunctionStellarSDSSConstructorParameters=outputAnalysisMassFunctionStellarSDSS(cosmologyFunctions_,randomError)
    return
  end function massFunctionStellarSDSSConstructorParameters

  function massFunctionStellarSDSSConstructorInternal(cosmologyFunctions_,randomError)
    !% Constructor for the ``massFunctionStellarSDSS'' output analysis class for internal use.
    use Input_Parameters2
    use Galacticus_Input_Paths
    implicit none
    type            (outputAnalysisMassFunctionStellarSDSS        )                         :: massFunctionStellarSDSSConstructorInternal
    class           (cosmologyFunctionsClass                      ), intent(in   ), target  :: cosmologyFunctions_
    double precision                                               , intent(in   )          :: randomError
    type            (galacticFilterStellarMass                    )               , pointer :: galacticFilter_
    type            (surveyGeometryLiWhite2009SDSS                )               , pointer :: surveyGeometry_
    type            (outputAnalysisDistributionOperatorRandomError)               , pointer :: outputAnalysisDistributionOperator_

    ! Build the SDSS survey geometry of Li & White (2008) with their imposed redshift limits.
    allocate(surveyGeometry_)
    surveyGeometry_=surveyGeometryLiWhite2009SDSS(redshiftMinimum=1.0d-3,redshiftMaximum=0.5d0)
    ! Build a filter which select galaxies with stellar mass 10⁶M☉ or greater.
    allocate(galacticFilter_)
    galacticFilter_=galacticFilterStellarMass(massThreshold=1.0d6)
    ! Build a random error distribution operator.
    allocate(outputAnalysisDistributionOperator_)
    outputAnalysisDistributionOperator_=outputAnalysisDistributionOperatorRandomError(randomError)
    ! Build the object.
    massFunctionStellarSDSSConstructorInternal%outputAnalysisMassFunctionStellar=                                                                               &
         & outputAnalysisMassFunctionStellar(                                                                                                                   &
         &                                   char(Galacticus_Input_Path()//'/data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.hdf5'), &
         &                                   galacticFilter_                                                                                                  , &
         &                                   surveyGeometry_                                                                                                  , &
         &                                   cosmologyFunctions_                                                                                              , &
         &                                   outputAnalysisDistributionOperator_                                                                                &
         &                                  )
    ! Clean up.
    nullify(surveyGeometry_                    )
    nullify(galacticFilter_                    )
    nullify(outputAnalysisDistributionOperator_)
    return
  end function massFunctionStellarSDSSConstructorInternal
