!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

!% Contains a module which implements writing of \gls{igm} state to the \glc\ output file.

module Intergalactic_Medium_Outputs
  !% Implements writing of \gls{igm} state to the \glc\ output file.
  implicit none
  private
  public :: Intergalactic_Medium_Output

contains

  !# <universePostEvolveTask>
  !#  <unitName>Intergalactic_Medium_Output</unitName>
  !# </universePostEvolveTask>
  subroutine Intergalactic_Medium_Output()
    !% Output \gls{igm} state to the \glc\ output file.
    use IO_HDF5
    use Galacticus_HDF5
    use Numerical_Ranges
    use Input_Parameters
    use Cosmology_Functions
    use Intergalactic_Medium_State
    use Numerical_Constants_Astronomical
    use Memory_Management
    implicit none
    class           (cosmologyFunctionsClass      ), pointer                   :: cosmologyFunctions_
    class           (intergalacticMediumStateClass), pointer                   :: intergalacticMediumState_
    double precision                               , allocatable, dimension(:) :: time                     , redshift       , &
         &                                                                        temperature              , massFiltering
    double precision                                                           :: redshiftMinimum          , redshiftMaximum, &
         &                                                                        timeMinimum              , timeMaximum
    type            (hdf5Object                   )                            :: intergalacticMediumGroup , dataSet
    logical                                                                    :: intergalacticMediumOutput
    integer                                                                    :: stepsPerDecade           , stepCount      , &
         &                                                                        i
    type            (inputParameters              )                            :: subParameters
    
    ! Test for output.
    !# <inputParameter>
    !#   <name>intergalacticMediumOutput</name>
    !#   <source>globalParameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>intergalacticMediumOutput</variable>
    !#   <description>If true, intergalactic medium state will be output.</description>
    !#   <type>boolean</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    if (intergalacticMediumOutput) then
       ! Get output epochs.
       subParameters=globalParameters%subParameters('intergalacticMediumOutput')
       !# <inputParameter>
       !#   <name>redshiftMinimum</name>
       !#   <source>subParameters</source>
       !#   <defaultValue>0.0d0</defaultValue>
       !#   <variable>redshiftMinimum</variable>
       !#   <description>Minimum redshift for which \gls{igm} properties should be output.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>redshiftMaximum</name>
       !#   <source>subParameters</source>
       !#   <defaultValue>1000.0d0</defaultValue>
       !#   <variable>redshiftMaximum</variable>
       !#   <description>Maximum redshift for which \gls{igm} properties should be output.</description>
       !#   <type>real</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       !# <inputParameter>
       !#   <name>stepsPerDecade</name>
       !#   <source>subParameters</source>
       !#   <defaultValue>10</defaultValue>
       !#   <variable>stepsPerDecade</variable>
       !#   <description>Number of steps per decade of time for at which \gls{igm} properties should be output.</description>
       !#   <type>integer</type>
       !#   <cardinality>1</cardinality>
       !# </inputParameter>
       ! Get required objects.
       cosmologyFunctions_       => cosmologyFunctions      ()
       intergalacticMediumState_ => intergalacticMediumState()
       ! Determine number of output epochs.
       timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
       timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
       stepCount  =int(log10(timeMaximum/timeMinimum)*dble(stepsPerDecade))+1
       ! Allocate IGM state arrays.
       call allocateArray(time         ,[stepCount])
       call allocateArray(redshift     ,[stepCount])
       call allocateArray(temperature  ,[stepCount])
       call allocateArray(massFiltering,[stepCount])
       ! Build array of times.
       time=Make_Range(timeMinimum,timeMaximum,stepCount,rangeType=rangeTypeLogarithmic)
       ! Populate arrays with IGM state.
       do i=1,stepCount
          redshift     (i)=cosmologyFunctions_      %redshiftFromExpansionFactor(cosmologyFunctions_%expansionFactor(time(i)))
          temperature  (i)=intergalacticMediumState_%temperature                (                                    time(i) )
          massFiltering(i)=intergalacticMediumState_%filteringMass              (                                    time(i) )
       end do
       call hdf5Access%set()
       intergalacticMediumGroup=galacticusOutputFile%openGroup('intergalacticMedium','Intergalactic medium state.')
       call intergalacticMediumGroup%writeDataset  (time           ,'time'         ,'Cosmic time'    ,datasetReturned=dataset)
       call dataset                 %writeAttribute('Gyr'          ,'units'                                                  )
       call dataset                 %writeAttribute(gigayear       ,'unitsInSI'                                              )
       call dataset                 %close()
       call intergalacticMediumGroup%writeDataset  (redshift       ,'redshift'     ,'Redshift'       ,datasetReturned=dataset)
       call dataset                 %writeAttribute('dimensionless','units'                                                  )
       call dataset                 %close()
       call intergalacticMediumGroup%writeDataset  (temperature    ,'temperature'  ,'IGM temperature',datasetReturned=dataset)
       call dataset                 %writeAttribute('K'            ,'units'                                                  )
       call dataset                 %writeAttribute(1.0d0          ,'unitsInSI'                                              )
       call dataset                 %close()
       call intergalacticMediumGroup%writeDataset  (massFiltering  ,'massFiltering','Filtering mass' ,datasetReturned=dataset)
       call dataset                 %writeAttribute('M☉'          ,'units'                                                  )
       call dataset                 %writeAttribute(massSolar      ,'unitsInSI'                                              )
       call dataset                 %close()
       call intergalacticMediumGroup%close()
       call hdf5Access%unset()
       ! Deallocate arrays.
       call deallocateArray(time         )       
       call deallocateArray(redshift     )       
       call deallocateArray(temperature  )       
       call deallocateArray(massFiltering)       
    end if
    return
  end subroutine Intergalactic_Medium_Output

end module Intergalactic_Medium_Outputs
