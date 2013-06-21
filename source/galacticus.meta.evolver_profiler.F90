!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which constructs a profile of \glc\ ODE evolver statistics.

module Galacticus_Meta_Evolver_Profiler
  !% Constructs a profile of \glc\ ODE evolver statistics.
  use Hashes
  implicit none
  private
  public :: Galacticus_Meta_Evolver_Profile, Galacticus_Meta_Evolver_Profiler_Output

  ! Arrays used to store profiling data.
  logical                                                        :: metaProfileInitialized    =.false.
  double precision                                               :: metaProfileTimeStepMaximum        , metaProfileTimeStepMinimum
  integer                                                        :: metaProfileTimeStepPoints         , metaProfileTimeStepPointsPerDecade
  double precision                   , allocatable, dimension(:) :: metaProfileTimeStep
  integer                            , allocatable, dimension(:) :: metaProfileTimeStepCount

  ! Hash for storing property hits.
  type            (integerScalarHash)                            :: propertyHits

contains

  subroutine Galacticus_Meta_Evolver_Profile(timeStep,propertyName)
    !% Record profiling information on the ODE evolver.
    use Input_Parameters
    use Memory_Management
    use Numerical_Ranges
    use ISO_Varying_String
    use Arrays_Search
    implicit none
    double precision                , intent(in   ) :: timeStep
    type            (varying_string), intent(in   ) :: propertyName
    integer                                         :: hitCount    , iStep

    !$omp critical (Meta_Profile_Record)
    if (.not.metaProfileInitialized) then
       ! Get parameters controlling profiling.
       !@ <inputParameter>
       !@   <name>metaProfileTimeStepMinimum</name>
       !@   <defaultValue>$10^{-6}$ Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The smallest timestep to use in profiling ODE solver steps.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('metaProfileTimeStepMinimum',metaProfileTimeStepMinimum,defaultValue=1.0d-6)
       !@ <inputParameter>
       !@   <name>metaProfileTimeStepMaximum</name>
       !@   <defaultValue>$10$ Gyr</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The largest timestep to use in profiling ODE solver steps.
       !@   </description>
       !@   <type>real</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('metaProfileTimeStepMaximum',metaProfileTimeStepMaximum,defaultValue=1.0d+1)
       !@ <inputParameter>
       !@   <name>metaProfileTimeStepPointsPerDecade</name>
       !@   <defaultValue>3</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The number of bins per decade of timestep to use when profiling ODE solver steps.
       !@   </description>
       !@   <type>integer</type>
       !@   <cardinality>1</cardinality>
       !@ </inputParameter>
       call Get_Input_Parameter('metaProfileTimeStepPointsPerDecade',metaProfileTimeStepPointsPerDecade,defaultValue=3)
       ! Create an array of timesteps.
       metaProfileTimeStepPoints=int(log10(metaProfileTimeStepMaximum/metaProfileTimeStepMinimum)*dble(metaProfileTimeStepPointsPerDecade))+1
       call Alloc_Array(metaProfileTimeStep     ,[metaProfileTimeStepPoints])
       call Alloc_Array(metaProfileTimeStepCount,[metaProfileTimeStepPoints])
       metaProfileTimeStep=Make_Range(metaProfileTimeStepMinimum,metaProfileTimeStepMaximum,metaProfileTimeStepPoints &
            &,rangeTypeLogarithmic)
       metaProfileTimeStepCount=0
       ! Initialize the property hits hash.
       call propertyHits%initialize()
       ! Flag that profiling is now initialized.
       metaProfileInitialized=.true.
    end if

    ! Accumulate the count of step sizes.
    iStep=Search_Array(metaProfileTimeStep,timeStep)
    metaProfileTimeStepCount(iStep)=metaProfileTimeStepCount(iStep)+1

    ! Record the limiting property.
    if (propertyHits%exists(propertyName)) then
       hitCount=propertyHits%value(propertyName)+1
    else
       hitCount=1
    end if
    call propertyHits%set(propertyName,hitCount)
    !$omp end critical (Meta_Profile_Record)
    return
  end subroutine Galacticus_Meta_Evolver_Profile

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Meta_Evolver_Profiler_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Meta_Evolver_Profiler_Output
    !% Outputs collected meta-data on tree evolution.
    use ISO_Varying_String
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    type   (varying_string), allocatable, dimension(:) :: metaProfilePropertyNames
    integer                , allocatable, dimension(:) :: metaProfilePropertyHitCount
    type   (hdf5Object    )                            :: metaDataDataset            , metaDataGroup, &
         &                                                profilerDataGroup

    ! Output tree evolution meta-data if any was collected.
    if (metaProfileInitialized) then

       ! Open output groups.
       metaDataGroup    =galacticusOutputFile%openGroup('metaData'       ,'Galacticus meta data.'     )
       profilerDataGroup=metaDataGroup       %openGroup('evolverProfiler','Meta-data on tree evolver.')

       ! Write timestep histogram.
       call profilerDataGroup%writeDataset(metaProfileTimeStep     ,"metaProfileTimeStep"     ,"Timestep [Gyr]"        ,datasetReturned=metaDataDataset)
       call metaDataDataset%writeAttribute(gigaYear,"unitsInSI")
       call metaDataDataset%close()
       call profilerDataGroup%writeDataset(metaProfileTimeStepCount,"metaProfileTimeStepCount","Timestep historgram []"                                )

       ! Write property histogram.
       call propertyHits%keys  (metaProfilePropertyNames   )
       call propertyHits%values(metaProfilePropertyHitCount)
       call profilerDataGroup%writeDataset(metaProfilePropertyNames   ,"metaProfilePropertyNames"   ,"Property names"    )
       call profilerDataGroup%writeDataset(metaProfilePropertyHitCount,"metaProfilePropertyHitCount","Property hit count")
       deallocate(metaProfilePropertyNames   )
       deallocate(metaProfilePropertyHitCount)

       ! Close output groups.
       call profilerDataGroup%close()
       call metaDataGroup    %close()

    end if

    return
  end subroutine Galacticus_Meta_Evolver_Profiler_Output

end module Galacticus_Meta_Evolver_Profiler
