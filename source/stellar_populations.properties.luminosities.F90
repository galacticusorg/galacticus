!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

module Stellar_Population_Properties_Luminosities
  use ISO_Varying_String
  implicit none
  private
  public :: Stellar_Population_Luminosities_Count, Stellar_Population_Luminosities_Get, Stellar_Population_Luminosities_Name,&
       & Stellar_Population_Luminosities_Output_Count, Stellar_Population_Luminosities_Output

  ! Flag specifying if module has been initialized.
  logical                                         :: luminositiesInitialized=.false.

  ! Arrays which hold the luminosity specifications.
  integer                                         :: luminosityCount
  integer,              allocatable, dimension(:) :: luminosityFilterIndex,luminosityIndex
  double precision,     allocatable, dimension(:) :: luminosityRedshift,luminosityCosmicTime
  type(varying_string), allocatable, dimension(:) :: luminosityFilter,luminosityType,luminosityName
  
  ! Luminosity output options.
  integer            :: luminosityOutputOption
  integer, parameter :: luminosityOutputOptionAll=0, luminosityOutputOptionFuture=1, luminosityOutputOptionPresent=2 

contains

  integer function Stellar_Population_Luminosities_Output_Count(time)
    !% Return a count of the number of luminosities to be output at {\tt time}.
    implicit none
    double precision, intent(in) :: time
    integer                      :: luminosityIndex

    Stellar_Population_Luminosities_Output_Count=0
    do luminosityIndex=1,luminosityCount
       if (Stellar_Population_Luminosities_Output(luminosityIndex,time)) Stellar_Population_Luminosities_Output_Count&
            &=Stellar_Population_Luminosities_Output_Count+1
    end do
    return
  end function Stellar_Population_Luminosities_Output_Count

  logical function Stellar_Population_Luminosities_Output(luminosityIndex,time)
    !% Return true or false depending on whether {\tt luminosityIndex} should be output at {\tt time}.
    use Galacticus_Error
    implicit none
    integer,          intent(in) :: luminosityIndex
    double precision, intent(in) :: time
    double precision, parameter  :: timeTolerance=1.0d-3

    select case (luminosityOutputOption)
    case (luminosityOutputOptionAll)
       Stellar_Population_Luminosities_Output=.true.
    case (luminosityOutputOptionFuture)
       Stellar_Population_Luminosities_Output=(     luminosityCosmicTime(luminosityIndex)       >= time*(1.0d0-timeTolerance))
    case (luminosityOutputOptionPresent)
       Stellar_Population_Luminosities_Output=(dabs(luminosityCosmicTime(luminosityIndex)-time) <= time*       timeTolerance )
    case default
       call Galacticus_Error_Report('Stellar_Population_Luminosities_Output','unknown luminosity output option')
    end select
    return
  end function Stellar_Population_Luminosities_Output

  function Stellar_Population_Luminosities_Name(index)
    !% Return a name for the specified entry in the luminosity list.
    use ISO_Varying_String
    implicit none
    type(varying_string)             :: Stellar_Population_Luminosities_Name
    integer,              intent(in) :: index
    
    Stellar_Population_Luminosities_Name=luminosityName(index)
    return
  end function Stellar_Population_Luminosities_Name
  
  function Stellar_Population_Luminosities_Get(imfSelected,time,abundancesStellar)
    !% Return the luminosity in each requested band for a stellar population of $1M_\odot$ with the specified {\tt abundancesStellar} and
    !% which formed at cosmological {\tt time} with IMF specified by {\tt imfSelected}.
    use Abundances_Structure
    use Stellar_Population_Luminosities
    implicit none
    double precision,          dimension(luminosityCount) :: Stellar_Population_Luminosities_Get
    integer,                   intent(in)                 :: imfSelected
    double precision,          intent(in)                 :: time
    type(abundances),          intent(in)                 :: abundancesStellar
    double precision,          dimension(luminosityCount) :: ages

    ! Ensure module is initialized.
    call Stellar_Population_Properties_Luminosities_Initialize

    ! Get the ages that this stellar population will have at the various output times.
    ages=luminosityCosmicTime-time

    ! Get the luminosities for each requested band.
    Stellar_Population_Luminosities_Get=Stellar_Population_Luminosity(luminosityIndex,luminosityFilterIndex,imfSelected,abundancesStellar&
         &,ages,luminosityRedshift)
    return
  end function Stellar_Population_Luminosities_Get

  integer function Stellar_Population_Luminosities_Count()
    !% Return the number of stellar luminosities being tracked.
    implicit none

    ! Ensure module is initialized.
    call Stellar_Population_Properties_Luminosities_Initialize

    ! Return the number of stellar luminosities being tracked.
    Stellar_Population_Luminosities_Count=luminosityCount
    return
  end function Stellar_Population_Luminosities_Count

  subroutine Stellar_Population_Properties_Luminosities_Initialize
    use Input_Parameters
    use Galacticus_Error
    use Memory_Management
    use Instruments_Filters
    use Cosmology_Functions
    implicit none
    integer              :: iLuminosity
    double precision     :: expansionFactor
    character(len=10)    :: redshiftLabel
    type(varying_string) :: luminosityOutputOptionText

    ! Initialize the module if necessary.
    if (.not.luminositiesInitialized) then
       !$omp critical (Stellar_Population_Properties_Luminosities_Initialize)
       if (.not.luminositiesInitialized) then
          ! Get luminosity output option.
          !@ <inputParameter>
          !@   <name>luminosityOutputOption</name>
          !@   <defaultValue>present</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Selects which luminosities will be output at each output time:
          !@     \begin{description}
          !@     \item [all] Output all luminosities;
          !@     \item [future] Output only those luminosities computed for the present output or future times;
          !@     \item [present] Output only those luminosities computed for the present output time.
          !@     \end{description}
          !@   </description>
          !@   <type>integer</type>
          !@   <cardinality>1</cardinality>
          !@   <group>output</group>
          !@ </inputParameter>
          call Get_Input_Parameter('luminosityOutputOption',luminosityOutputOptionText,defaultValue="present")
          select case (char(luminosityOutputOptionText))
          case ("all")
             luminosityOutputOption=luminosityOutputOptionAll
          case ("future")
             luminosityOutputOption=luminosityOutputOptionFuture
          case ("present")
             luminosityOutputOption=luminosityOutputOptionPresent
          case default
             call Galacticus_Error_Report("Stellar_Population_Properties_Luminosities_Initialize","unrecognized luminosityOutputOption")
          end select
          
          ! Read in the parameters which specify the luminosities to be computed.
          luminosityCount=Get_Input_Parameter_Array_Size('luminosityRedshift')
          if (Get_Input_Parameter_Array_Size('luminosityFilter') /= luminosityCount) call&
               & Galacticus_Error_Report('Stellar_Population_Properties_Luminosities_Initialize','luminosityFilter and luminosityCount&
               & input arrays must have same dimension')
          if (Get_Input_Parameter_Array_Size('luminosityType') /= luminosityCount) call&
               & Galacticus_Error_Report('Stellar_Population_Properties_Luminosities_Initialize','luminosityType and luminosityCount&
               & input arrays must have same dimension')
          
          if (luminosityCount > 0) then
             call Alloc_Array(luminosityRedshift   ,[luminosityCount])
             call Alloc_Array(luminosityCosmicTime ,[luminosityCount])
             allocate(luminosityFilter(luminosityCount))
             allocate(luminosityType  (luminosityCount))
             allocate(luminosityName  (luminosityCount))
             call Memory_Usage_Record(sizeof(luminosityFilter)+sizeof(luminosityType)+sizeof(luminosityName),blockCount=3)
             call Alloc_Array(luminosityFilterIndex,[luminosityCount])
             call Alloc_Array(luminosityIndex      ,[luminosityCount])
             !@ <inputParameter>
             !@   <name>luminosityRedshift</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The redshift for which to compute each specified stellar luminosity.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>0..*</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('luminosityRedshift',luminosityRedshift)
             !@ <inputParameter>
             !@   <name>luminosityFilter</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The filter name for each stellar luminosity to be computed.
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>0..*</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('luminosityFilter'  ,luminosityFilter  )
             !@ <inputParameter>
             !@   <name>luminosityType</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The luminosity type for each stellar luminosity to be computed:
             !@     \begin{description}
             !@      \item[rest] Compute luminosity in the galaxy rest frame;
             !@      \item[observed] Compute luminosity in the observer frame\footnote{The luminosity computed in this way is that in the galaxy rest
             !@                      frame using a filter blueshifted to the galaxy's redshift. This means that to compute an apparent magnitude you
             !@                      must add not only the distance modulus, but a factor of $-2.5\log_{10}(1+z)$ to account for compression of photon
             !@                      frequencies.}.
             !@     \end{description}
             !@   </description>
             !@   <type>string</type>
             !@   <cardinality>0..*</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('luminosityType'    ,luminosityType    )
             
             ! Process the list of luminosities.
             do iLuminosity=1,luminosityCount
                ! Assign a name to this luminosity.
                write (redshiftLabel,'(f7.4)') luminosityRedshift(iLuminosity)
                luminosityName(iLuminosity)=luminosityFilter(iLuminosity)//":"//luminosityType(iLuminosity)//":z"&
                     &//trim(adjustl(redshiftLabel))
                ! Assign an index for this luminosity.
                luminosityIndex(iLuminosity)=iLuminosity
                ! Get the index of the specified filter.
                luminosityFilterIndex(iLuminosity)=Filter_Get_Index(luminosityFilter(iLuminosity))
                ! Set the reference time (i.e. cosmological time corresponding to the specified redshift) for this filter.
                expansionFactor=Expansion_Factor_from_Redshift(luminosityRedshift(iLuminosity))
                luminosityCosmicTime(iLuminosity)=Cosmology_Age(expansionFactor)
                ! Set the filter redshifting factor. This is equal to the requested redshift if an observed frame was specified, otherwise
                ! it is set to zero to indicate a rest-frame filter.
                select case(char(luminosityType(iLuminosity)))
                case ("rest")
                   luminosityRedshift(iLuminosity)=0.0d0
                case ("observed")
                   ! Do nothing, we already have the correct redshift.
                case default
                   call Galacticus_Error_Report('Stellar_Population_Properties_Luminosities_Initialize','unrecognized filter type')
                end select
             end do
          end if
          luminositiesInitialized=.true.
       end if
       !$omp end critical (Stellar_Population_Properties_Luminosities_Initialize)
    end if
    return
  end subroutine Stellar_Population_Properties_Luminosities_Initialize
  
end module Stellar_Population_Properties_Luminosities
