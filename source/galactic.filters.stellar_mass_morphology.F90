!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025, 2026
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
Implements a galactic high-pass filter for stellar mass-weighted morphology (i.e. spheroid-to-total ratio).
!!}

  !![
  <galacticFilter name="galacticFilterStellarMassMorphology">
   <description>
   A galactic high-pass filter for stellar mass-weighted morphology (i.e. spheroid-to-total ratio). Galaxies with a
   spheroid-to-total ratio (by stellar mass) greater than or equal to a fixed threshold, $R_{\star,0}=${\normalfont \ttfamily
   [spheroidToTotalThreshold]}.
   </description>
  </galacticFilter>
  !!]
  type, extends(galacticFilterClass) :: galacticFilterStellarMassMorphology
     !!{
     A galactic high-pass filter class for stellar mass-weighted morphology (i.e. spheroid-to-total ratio).
     !!}
     private
     double precision :: spheroidToTotalRatioThreshold
   contains
     procedure :: passes => stellarMassMorphologyPasses
  end type galacticFilterStellarMassMorphology

  interface galacticFilterStellarMassMorphology
     !!{
     Constructors for the \refClass{galacticFilterStellarMassMorphology} galactic filter class.
     !!}
     module procedure stellarMassMorphologyConstructorParameters
     module procedure stellarMassMorphologyConstructorInternal
  end interface galacticFilterStellarMassMorphology

contains

  function stellarMassMorphologyConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{galacticFilterStellarMassMorphology} galactic filter class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (galacticFilterStellarMassMorphology)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    double precision                                                     :: spheroidToTotalRatioThreshold
    !![
    <inputParameter>
      <name>spheroidToTotalRatioThreshold</name>
      <source>parameters</source>
      <description>The parameter $R_0$ appearing in the stellar mass-weight morphology threshold for the stellar mass-weighted morphology galactic filter class.</description>
    </inputParameter>
    !!]
    self=galacticFilterStellarMassMorphology(spheroidToTotalRatioThreshold)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function stellarMassMorphologyConstructorParameters

  function stellarMassMorphologyConstructorInternal(spheroidToTotalRatioThreshold) result(self)
    !!{
    Internal constructor for the ``stellarMassMorphology'' galactic filter class.
    !!}
    implicit none
    type            (galacticFilterStellarMassMorphology)                :: self
    double precision                                     , intent(in   ) :: spheroidToTotalRatioThreshold
    !![
    <constructorAssign variables="spheroidToTotalRatioThreshold"/>
    !!]
    
    return
  end function stellarMassMorphologyConstructorInternal

  logical function stellarMassMorphologyPasses(self,node)
    !!{
    Implement a stellar mass-weighted morphology high-pass galactic filter.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentDisk, nodeComponentSpheroid, treeNode
    implicit none
    class           (galacticFilterStellarMassMorphology), intent(inout)         :: self
    type            (treeNode                           ), intent(inout), target :: node
    class           (nodeComponentDisk                  ), pointer               :: disk
    class           (nodeComponentSpheroid              ), pointer               :: spheroid
    double precision                                                             :: stellarMass

    disk              => node    %disk       ()
    spheroid          => node    %spheroid   ()
    stellarMass       = +disk    %massStellar() &
         &              +spheroid%massStellar()
    if (stellarMass > 0.0d0) then
       stellarMassMorphologyPasses= +spheroid%massStellar()            &
            &                       /stellarMass                       &
            &                      >=                                  &
            &                       self%spheroidToTotalRatioThreshold
    else
       stellarMassMorphologyPasses=.false.
    end if
    return
  end function stellarMassMorphologyPasses
