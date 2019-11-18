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
  
  !% Implements a merger tree evolve profiler that collects simple data.

  use :: Hashes, only : integerScalarHash
  
  !# <mergerTreeEvolveProfiler name="mergerTreeEvolveProfilerSimple">
  !#  <description>A merger tree evolve profiler that collects simple data.</description>
  !# </mergerTreeEvolveProfiler>
  type, extends(mergerTreeEvolveProfilerClass) :: mergerTreeEvolveProfilerSimple
     !% A merger tree evolve profiler that collects simple data.
     private
     double precision                                               :: timeStepMaximum        , timeStepMinimum
     integer                                                        :: timeStepPointsPerDecade
     double precision                   , allocatable, dimension(:) :: timeStep
     integer                            , allocatable, dimension(:) :: timeStepCount
     type            (integerScalarHash)                            :: propertyHits
   contains
     final     ::            simpleDestructor
     procedure :: profile => simpleProfile
  end type mergerTreeEvolveProfilerSimple

  interface mergerTreeEvolveProfilerSimple
     !% Constructors for the {\normalfont \ttfamily simple} merger tree evolve profiler class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface mergerTreeEvolveProfilerSimple

contains
  
  function simpleConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily simple} merger tree evolve profiler class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameters
    implicit none
    type            (mergerTreeEvolveProfilerSimple)                :: self
    type            (inputParameters               ), intent(inout) :: parameters
    double precision                                                :: timeStepMinimum        , timeStepMaximum
    integer                                                         :: timeStepPointsPerDecade

    !# <inputParameter>
    !#   <name>timeStepMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d-6</defaultValue>
    !#   <description>The smallest timestep to use in profiling ODE solver steps.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeStepMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d+1</defaultValue>
    !#   <description>The largest timestep to use in profiling ODE solver steps.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeStepPointsPerDecade</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>3</defaultValue>
    !#   <description>The number of bins per decade of timestep to use when profiling ODE solver steps.</description>
    !#   <source>parameters</source>
    !#   <type>integer</type>
    !# </inputParameter>
    self=mergerTreeEvolveProfilerSimple(timeStepMinimum,timeStepMaximum,timeStepPointsPerDecade)
    !# <inputParametersValidate source="parameters"/>
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(timeStepMinimum,timeStepMaximum,timeStepPointsPerDecade) result(self)
    use :: Numerical_Ranges, only : Make_Range, rangeTypeLogarithmic
    implicit none
    type            (mergerTreeEvolveProfilerSimple)                :: self
    double precision                                , intent(in   ) :: timeStepMinimum        , timeStepMaximum
    integer                                         , intent(in   ) :: timeStepPointsPerDecade
    integer                                                         :: timeStepPoints
    !# <constructorAssign variables="timeStepMinimum, timeStepMaximum, timeStepPointsPerDecade"/>
    
    timeStepPoints=int(log10(timeStepMaximum/timeStepMinimum)*dble(timeStepPointsPerDecade))+1
    allocate(self%timeStep     (timeStepPoints))
    allocate(self%timeStepCount(timeStepPoints))
    self%timeStep     =Make_Range(timeStepMinimum,timeStepMaximum,timeStepPoints,rangeTypeLogarithmic)
    self%timeStepCount=0
    call self%propertyHits%initialize()
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !% Output collected meta-data on tree evolution.
    use :: Galacticus_HDF5                 , only : galacticusOutputFile
    use :: IO_HDF5                         , only : hdf5Object          , hdf5Access
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    type            (mergerTreeEvolveProfilerSimple), intent(inout)               :: self
    type            (varying_string                ), allocatable  , dimension(:) :: propertyNames
    integer                                         , allocatable  , dimension(:) :: propertyHitCount , timeStepCountPrevious
    type            (hdf5Object                    )                              :: metaDataDataset  , metaDataGroup, &
         &                                                                           profilerDataGroup
    integer                                                                       :: i                , hitCount
    
    if (all(self%timestepCount == 0)) return
    !$ call hdf5Access%set  ()
    metaDataGroup    =galacticusOutputFile%openGroup('metaData'       ,'Galacticus meta data.'     )
    profilerDataGroup=metaDataGroup       %openGroup('evolverProfiler','Meta-data on tree evolver.')
    if (profilerDataGroup%hasDataset('timeStepCount')) then
       call profilerDataGroup%readDataset('timeStepCount',timeStepCountPrevious)
       self%timeStepCount=self%timeStepCount+timeStepCountPrevious
    end if
    if (profilerDataGroup%hasDataset('propertyNames')) then
       call profilerDataGroup%readDataset("propertyNames"   ,propertyNames   )
       call profilerDataGroup%readDataset("propertyHitCount",propertyHitCount)
       do i=1,size(propertyNames)
          if (self%propertyHits%exists(propertyNames(i))) then
             hitCount=self%propertyHits%value(propertyNames(i))+propertyHitCount(i)
          else
             hitCount=                                          propertyHitCount(i)
          end if
          call self%propertyHits%set(propertyNames(i),hitCount)
       end do
       deallocate(propertyNames   )
       deallocate(propertyHitCount)
       call profilerDataGroup%remove("propertyNames"   )
       call profilerDataGroup%remove("propertyHitCount")
    end if
    call self%propertyHits%keys  (propertyNames   )
    call self%propertyHits%values(propertyHitCount)
    call profilerDataGroup%writeDataset  (self%timeStep        ,"timeStep"        ,"Timestep [Gyr]"        ,datasetReturned=metaDataDataset)
    call metaDataDataset  %writeAttribute(     gigaYear        ,"unitsInSI"                                                                )
    call metaDataDataset  %close         (                                                                                                 )
    call profilerDataGroup%writeDataset  (self%timeStepCount   ,"timeStepCount"   ,"Timestep histogram []"                                 )
    call profilerDataGroup%writeDataset  (     propertyNames   ,"propertyNames"   ,"Property names"                                        )
    call profilerDataGroup%writeDataset  (     propertyHitCount,"propertyHitCount","Property hit count"                                    )
    call profilerDataGroup%close         (                                                                                                 )
    call metaDataGroup    %close         (                                                                                                 )
    !$ call hdf5Access%unset()
    return
  end subroutine simpleDestructor

  subroutine simpleProfile(self,timestep,propertyName)
    !% Profile the differential evolution step.
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    use            :: Arrays_Search, only : Search_Array
    implicit none
    class           (mergerTreeEvolveProfilerSimple), intent(inout) :: self
    double precision                                , intent(in   ) :: timeStep
    type            (varying_string                ), intent(in   ) :: propertyName
    integer                                                         :: hitCount
    integer         (c_size_t                      )                :: i

    i                    =Search_Array(self%timeStep,timeStep)
    self%timeStepCount(i)=self%timeStepCount(i)+1
    if (self%propertyHits%exists(propertyName)) then
       hitCount=self%propertyHits%value(propertyName)+1
    else
       hitCount=1
    end if
    call self%propertyHits%set(propertyName,hitCount)
    return
  end subroutine simpleProfile
