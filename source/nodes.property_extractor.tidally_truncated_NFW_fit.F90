!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
Implements a node property extractor that fits for a tidal truncation radius for an NFW profile.
!!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass, darkMatterProfileDMONFW
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Galactic_Structure      , only : galacticStructureClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTidallyTruncatedNFWFit">
   <description>A node property extractor that fits for a tidal truncation radius for an NFW profile.</description>
  <deepCopy>
    <functionClass variables="darkMatterProfileDMONFW_"/>
   </deepCopy>
   <stateStorable>
    <functionClass variables="darkMatterProfileDMONFW_"/>
   </stateStorable>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorTuple) :: nodePropertyExtractorTidallyTruncatedNFWFit
     !!{
     A node property extractor that fits for a tidal truncation radius for an NFW profile.
     !!}
     private
     type (darkMatterProfileDMONFW  ), pointer :: darkMatterProfileDMONFW_ => null()
     class(darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_    => null()
     class(darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_     => null()
     class(galacticStructureClass   ), pointer :: galacticStructure_       => null()
   contains
     final     ::                 tidallyTruncatedNFWFitDestructor
     procedure :: elementCount => tidallyTruncatedNFWFitElementCount
     procedure :: extract      => tidallyTruncatedNFWFitExtract
     procedure :: names        => tidallyTruncatedNFWFitNames
     procedure :: descriptions => tidallyTruncatedNFWFitDescriptions
     procedure :: unitsInSI    => tidallyTruncatedNFWFitUnitsInSI
  end type nodePropertyExtractorTidallyTruncatedNFWFit

  interface nodePropertyExtractorTidallyTruncatedNFWFit
     !!{
     Constructors for the ``tidallyTruncatedNFWFit'' output analysis class.
     !!}
     module procedure tidallyTruncatedNFWFitConstructorParameters
     module procedure tidallyTruncatedNFWFitConstructorInternal
  end interface nodePropertyExtractorTidallyTruncatedNFWFit

contains

  function tidallyTruncatedNFWFitConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``tidallyTruncatedNFWFit'' output analysis property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodePropertyExtractorTidallyTruncatedNFWFit)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(darkMatterProfileDMOClass                  ), pointer       :: darkMatterProfileDMO_
    class(darkMatterHaloScaleClass                   ), pointer       :: darkMatterHaloScale_
    class(galacticStructureClass                     ), pointer       :: galacticStructure_

    !![
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="galacticStructure"    name="galacticStructure_"    source="parameters"/>
    !!]
    self=nodePropertyExtractorTidallyTruncatedNFWFit(darkMatterHaloScale_,galacticStructure_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="galacticStructure_"   />
    !!]
    return
  end function tidallyTruncatedNFWFitConstructorParameters

  function tidallyTruncatedNFWFitConstructorInternal(darkMatterHaloScale_,galacticStructure_,darkMatterProfileDMO_) result(self)
    !!{
    Internal constructor for the ``tidallyTruncatedNFWFit'' output analysis property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorTidallyTruncatedNFWFit)                        :: self
    class(darkMatterProfileDMOClass                  ), intent(in   ), target :: darkMatterProfileDMO_
    class(darkMatterHaloScaleClass                   ), intent(in   ), target :: darkMatterHaloScale_
    class(galacticStructureClass                     ), intent(in   ), target :: galacticStructure_
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *galacticStructure_, *darkMatterProfileDMO_"/>
    !!]
    
    allocate(self%darkMatterProfileDMONFW_)
    !![
    <referenceConstruct isResult="yes" owner="self" object="darkMatterProfileDMONFW_">
      <constructor>
	darkMatterProfileDMONFW(                                                           &amp;
	 &amp;                  velocityDispersionUseSeriesExpansion=.false.             , &amp;
	 &amp;                  darkMatterHaloScale_                =darkMatterHaloScale_  &amp;
	 &amp;                 )
      </constructor>
    </referenceConstruct>
    !!]
    return
  end function tidallyTruncatedNFWFitConstructorInternal

  subroutine tidallyTruncatedNFWFitDestructor(self)
    !!{
    Destructor for the ``tidallyTruncatedNFWFit'' output analysis property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorTidallyTruncatedNFWFit), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterProfileDMONFW_"/>
    <objectDestructor name="self%darkMatterProfileDMO_"   />
    <objectDestructor name="self%darkMatterHaloScale_"    />
    <objectDestructor name="self%galacticStructure_"      />
    !!]
    return
  end subroutine tidallyTruncatedNFWFitDestructor

  integer function tidallyTruncatedNFWFitElementCount(self,time)
    !!{
    Return the number of elements in the {\normalfont \ttfamily tidallyTruncatedNFWFit} property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorTidallyTruncatedNFWFit), intent(inout) :: self
    double precision                                             , intent(in   ) :: time
    !$GLC attributes unused :: self, time
    
    tidallyTruncatedNFWFitElementCount=3   
    return
  end function tidallyTruncatedNFWFitElementCount

   function tidallyTruncatedNFWFitExtract(self,node,time,instance)
    !!{
    Implement a tidallyTruncatedNFWFit output analysis.
    !!}
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: Galacticus_Nodes          , only : nodeComponentDarkMatterProfile, nodeComponentSatellite
    use            :: Numerical_Ranges          , only : Make_Range                    , rangeTypeLogarithmic
    use            :: Multidimensional_Minimizer, only : multiDMinimizer
    implicit none
    double precision                                             , allocatable  , dimension(:) :: tidallyTruncatedNFWFitExtract
    class           (nodePropertyExtractorTidallyTruncatedNFWFit), intent(inout), target       :: self
    type            (treeNode                                   ), intent(inout), target       :: node
    double precision                                             , intent(in   )               :: time
    type            (multiCounter                               ), intent(inout), optional     :: instance
    class           (nodeComponentDarkMatterProfile             ), pointer                     :: darkMatterProfile
    class           (nodeComponentSatellite                     ), pointer                     :: satellite
    double precision                                             , allocatable  , dimension(:) :: radii                     , fractionDensity
    type            (multiDMinimizer                            ), allocatable                 :: minimizer_
    double precision                                                            , dimension(1) :: locationMinimum
    double precision                                             , parameter                   :: fractionRadiusOuter=0.95d0, fractionRadiusScale=0.1d0, &
         &                                                                                        fractionMaximum    =0.10d0, fractionStep       =0.1d0
    integer                                                      , parameter                   :: countRadiiPerDecade=10
    integer                                                                                    :: countRadii                , i                        , &
         &                                                                                        iteration
    logical                                                                                    :: converged
    double precision                                                                           :: radiusOuter               , massTotal                , &
         &                                                                                        radiusMinimum             , radiusMaximum            , &
         &                                                                                        radiusScale               , radiusVirial
    !$GLC attributes unused :: instance

    allocate(tidallyTruncatedNFWFitExtract(3))
    darkMatterProfile => node                                  %darkMatterProfile  ()
    tidallyTruncatedNFWFitExtract(3)=self%darkMatterProfileDMONFW_%density(node,darkMatterProfile%scale())
    if (node%isSatellite()) then
       ! Extract required properties.
       satellite         => node                                  %satellite          (                                              )
       massTotal         =  self             %galacticStructure_  %massEnclosed       (node                                          )
       radiusOuter       =  self             %galacticStructure_  %radiusEnclosingMass(node,mass=min(satellite%boundMass(),massTotal))
       radiusScale       =  darkMatterProfile                     %scale              (                                              )
       radiusVirial      =  self             %darkMatterHaloScale_%radiusVirial       (node                                          )
       ! Choose radii for fitting.
       radiusMaximum=    fractionRadiusOuter*radiusOuter
       radiusMinimum=min(fractionRadiusScale*radiusScale,fractionMaximum*radiusMaximum)
       countRadii   =int(log10(radiusMaximum/radiusMinimum)*dble(countRadiiPerDecade)+1.0d0)
       radii        =Make_Range(radiusMinimum,radiusMaximum,countRadii,rangeTypeLogarithmic)
       ! Tabulate the density ratio relative to an NFW profile.
       allocate(fractionDensity(countRadii))
       do i=1,countRadii
          fractionDensity(i)=+self%darkMatterProfileDMO_   %density(node,radii(i)) &
               &             /self%darkMatterProfileDMONFW_%density(node,radii(i))
       end do
       ! Optimize the fit.
       allocate(minimizer_)
       minimizer_=multiDMinimizer(1_c_size_t,fitMetric)
       call minimizer_%set(x=[log(radiusVirial)],stepSize=[fractionStep])
       iteration=0
       converged=.false.
       do while (.not.converged .and. iteration < 100)
          call minimizer_%iterate()
          iteration=iteration+1
          converged=minimizer_%testSize(toleranceAbsolute=1.0d-3*radiusScale)
       end do
       locationMinimum                 =minimizer_%x()
       tidallyTruncatedNFWFitExtract(1)=exp      (locationMinimum(1))
       tidallyTruncatedNFWFitExtract(2)=fitMetric(locationMinimum   )
    else
       tidallyTruncatedNFWFitExtract(1)=huge(0.0d0)
       tidallyTruncatedNFWFitExtract(2)=     0.0d0
    end if
    return
    
  contains
    
    double precision function fitMetric(properties)
      !!{
      Evaluate the fit metric.
      !!}
      implicit none
      double precision, intent(in   ), dimension(         :) :: properties
      double precision               , dimension(countRadii) :: fractionTruncation
      double precision                                       :: radiusTruncation
      
      radiusTruncation  =+exp(properties(1))
      fractionTruncation=+  1.0d0              &
           &             /(                    &
           &               +1.0d0              &
           &               +(                  &
           &                 +radii            &
           &                 /radiusTruncation &
           &                )**2               &
           &              )
      fitMetric         =+sum(                           &
           &                  +log10(                    &
           &                         +fractionDensity    &
           &                         /fractionTruncation &
           &                        )**2                 &
           &                 )                           &
           &             /dble(countRadii)
      return
    end function fitMetric
    
  end function tidallyTruncatedNFWFitExtract
   
  subroutine tidallyTruncatedNFWFitNames(self,time,names)
    !!{
    Return the name of the best-fit radius of a tidally-truncated NFW profile.
    !!}
    implicit none
    class           (nodePropertyExtractorTidallyTruncatedNFWFit), intent(inout)                             :: self
    double precision                                             , intent(in   )                             :: time
    type            (varying_string                             ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: self, time
    
    allocate(names(3))
    names(1)=var_str('radiusTidalTruncationNFW')
    names(2)=var_str('metricTidalTruncationNFW')
    names(3)=var_str('densityNormalizationTidalTruncationNFW')
    return
  end subroutine tidallyTruncatedNFWFitNames

  subroutine tidallyTruncatedNFWFitDescriptions(self,time,descriptions)
    !!{
    Return a description of a tidally-truncated NFW profile.
    !!}
    implicit none
    class           (nodePropertyExtractorTidallyTruncatedNFWFit), intent(inout)                             :: self
    double precision                                             , intent(in   )                             :: time
    type            (varying_string                             ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: self, time

    allocate(descriptions(3))
    descriptions(1)=var_str('The best-fit tidal truncation radius assuming an underlying NFW profile.'    )
    descriptions(2)=var_str('The best-fit tidal truncation fit metric assuming an underlying NFW profile.')
    descriptions(3)=var_str('The best-fit tidal truncation density normalization, assuming and underling NFW Profile.')
    return
  end subroutine tidallyTruncatedNFWFitDescriptions

  function tidallyTruncatedNFWFitUnitsInSI(self,time)
    !!{
    Return the units of a tidally-truncated NFW profile.
    !!}
    use :: Numerical_Constants_Astronomical, only : megaParsec, massSolar
    implicit none
    double precision                                             , dimension(:) , allocatable :: tidallyTruncatedNFWFitUnitsInSI
    class           (nodePropertyExtractorTidallyTruncatedNFWFit), intent(inout)              :: self
    double precision                                             , intent(in   )              :: time
    !$GLC attributes unused :: self, time

    allocate(tidallyTruncatedNFWFitUnitsInSI(3))
    tidallyTruncatedNFWFitUnitsInSI(1)=megaParsec
    tidallyTruncatedNFWFitUnitsInSI(2)=1.0d0
    tidallyTruncatedNFWFitUnitsInSI(3)=massSolar / megaParsec**3
    return
  end function tidallyTruncatedNFWFitUnitsInSI
