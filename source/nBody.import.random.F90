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
Implements an N-body data importer which generates random points.
!!}

  use :: Numerical_Random_Numbers, only : randomNumberGeneratorClass

  !![
  <nbodyImporter name="nbodyImporterRandom">
   <description>An importer which generates random points.</description>
  </nbodyImporter>
  !!]
  type, extends(nbodyImporterClass) :: nbodyImporterRandom
     !!{
     An importer which generates random points.
     !!}
     private
     class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_ => null()
     double precision                            , dimension(2)  :: xRange                          , yRange, &
          &                                                         zRange
     integer         (c_size_t                  )                :: countPoints
   contains
     final     ::           randomDestructor
     procedure :: import => randomImport
     procedure :: isHDF5 => randomIsHDF5
  end type nbodyImporterRandom

  interface nbodyImporterRandom
     !!{
     Constructors for the {\normalfont \ttfamily random} N-body importer class.
     !!}
     module procedure randomConstructorParameters
     module procedure randomConstructorInternal
  end interface nbodyImporterRandom

contains
  
  function randomConstructorParameters(parameters) result (self)
    !!{
    Constructor for the {\normalfont \ttfamily random} N-body importer class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nbodyImporterRandom       )                :: self
    type            (inputParameters           ), intent(inout) :: parameters
    class           (randomNumberGeneratorClass), pointer       :: randomNumberGenerator_ => null()
    double precision                            , dimension(2)  :: xRange                          , yRange, &
         &                                                         zRange
    integer         (c_size_t                  )                :: countPoints

    !![
    <inputParameter>
      <name>countPoints</name>
      <source>parameters</source>
      <description>The number of random points to generate.</description>
    </inputParameter>
    <inputParameter>
      <name>xRange</name>
      <source>parameters</source>
      <description>The range within which to generate points in the $x$-direction.</description>
    </inputParameter>
    <inputParameter>
      <name>yRange</name>
      <source>parameters</source>
      <description>The range within which to generate points in the $y$-direction.</description>
    </inputParameter>
    <inputParameter>
      <name>zRange</name>
      <source>parameters</source>
      <description>The range within which to generate points in the $z$-direction.</description>
    </inputParameter>
    <objectBuilder class="randomNumberGenerator" name="randomNumberGenerator_" source="parameters"/>
    !!]
    self=nbodyImporterRandom(countPoints,xRange,yRange,zRange,randomNumberGenerator_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="randomNumberGenerator_"/>
    !!]
    return
  end function randomConstructorParameters

  function randomConstructorInternal(countPoints,xRange,yRange,zRange,randomNumberGenerator_) result (self)
    !!{
    Internal constructor for the {\normalfont \ttfamily random} N-body importer class.
    !!}
    implicit none
    type            (nbodyImporterRandom       )                              :: self
    integer         (c_size_t                  ), intent(in   )               :: countPoints
    double precision                            , intent(in   ), dimension(2) :: xRange                , yRange, &
         &                                                                       zRange
    class           (randomNumberGeneratorClass), intent(in   ), target       :: randomNumberGenerator_
    !![
    <constructorAssign variables="countPoints, xRange, yRange, zRange, *randomNumberGenerator_"/>
    !!]

    return
  end function randomConstructorInternal

  subroutine randomDestructor(self)
    !!{
    Destructor for {\normalfont \ttfamily random} importer class.
    !!}
    implicit none
    type(nbodyImporterRandom), intent(inout) :: self
    
    !![
    <objectDestructor name="self%randomNumberGenerator_"/>
    !!]
    return
  end subroutine randomDestructor

  subroutine randomImport(self,simulations)
    !!{
    Import data from a Random file.
    !!}
    use :: Display, only : displayIndent     , displayUnindent         , verbosityLevelStandard
    use :: Hashes , only : rank1DoublePtrHash, rank1IntegerSizeTPtrHash, rank2DoublePtrHash    , rank2IntegerSizeTPtrHash, &
         &                 doubleHash        , integerSizeTHash        , varyingStringHash     , genericHash
    implicit none
    class           (nbodyImporterRandom), intent(inout)                              :: self
    type            (nBodyData          ), intent(  out), dimension(  :), allocatable :: simulations
    double precision                                    , dimension(:,:), pointer     :: position   , velocity
    integer         (c_size_t           )               , dimension(  :), pointer     :: particleID
    integer         (c_size_t           )                                             :: i

    call displayIndent('import simulation from Random file',verbosityLevelStandard)
    allocate(simulations(1                 ))
    allocate(particleID (  self%countPoints))
    allocate(position   (3,self%countPoints))
    allocate(velocity   (3,self%countPoints))
    simulations(1)%label='random'
    simulations(1)%propertiesInteger     =rank1IntegerSizeTPtrHash()
    simulations(1)%propertiesReal        =rank1DoublePtrHash      ()
    simulations(1)%propertiesIntegerRank1=rank2IntegerSizeTPtrHash()
    simulations(1)%propertiesRealRank1   =rank2DoublePtrHash      ()
    simulations(1)%attributesInteger     =integerSizeTHash        ()
    simulations(1)%attributesReal        =doubleHash              ()
    simulations(1)%attributesText        =varyingStringHash       ()
    simulations(1)%attributesGeneric     =genericHash             ()
    do i=1_c_size_t,self%countPoints
       particleID(  i)=i
       position  (:,i)=[self%randomNumberGenerator_%uniformSample(),self%randomNumberGenerator_%uniformSample(),self%randomNumberGenerator_%uniformSample()]
    end do
    velocity(:,:)=0.0d0
    position(1,:)=position(1,:)*(self%xRange(2)-self%xRange(1))+self%xRange(1)
    position(2,:)=position(2,:)*(self%yRange(2)-self%yRange(1))+self%yRange(1)
    position(3,:)=position(3,:)*(self%zRange(2)-self%zRange(1))+self%zRange(1)
    call simulations(1)%propertiesInteger  %set('particleID',particleID)
    call simulations(1)%propertiesRealRank1%set('position'  ,position  )
    call simulations(1)%propertiesRealRank1%set('velocity'  ,velocity  )
    call displayUnindent('done',verbosityLevelStandard)
    return
  end subroutine randomImport

  logical function randomIsHDF5(self)
    !!{
    Return whether or not the imported data is from an HDF5 file.
    !!}
    implicit none
    class(nbodyImporterRandom), intent(inout) :: self

    randomIsHDF5=.false.
    return
  end function randomIsHDF5
