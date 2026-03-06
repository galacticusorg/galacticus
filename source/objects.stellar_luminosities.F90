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

!+    Contributions to this file made by:  Alex Merson.

!!{
Contains a module which defines the stellar luminosities object.
!!}

module Stellar_Luminosities_Structure
  !!{
  Defines the stellar luminosities object.
  !!}
  use :: ISO_Varying_String                    , only : varying_string
  use :: Stellar_Population_Spectra_Postprocess, only : stellarPopulationSpectraPostprocessorBuilderClass, stellarPopulationSpectraPostprocessorList
  implicit none
  private
  public :: stellarLuminosities                    , max                                      , &
       &    abs                                    , operator(*)                              , &
       &    Stellar_Luminosities_Parameter_Map     , Stellar_Luminosities_State_Store         , &
       &    Stellar_Luminosities_State_Restore     , Stellar_Luminosities_Initializor         , &
       &    Stellar_Luminosities_Thread_Initializor, Stellar_Luminosities_Thread_Uninitializor
  
  ! Interface to parameter mapping functions.
  interface Stellar_Luminosities_Parameter_Map
     module procedure Stellar_Luminosities_Parameter_Map_Double
  end interface Stellar_Luminosities_Parameter_Map

  ! Interface to max() function for stellar luminosities objects.
  interface max
     module procedure stellarLuminositiesMax
  end interface max

  ! Interface to abs() function for stellar luminosities objects.
  interface abs
     module procedure stellarLuminositiesAbs
  end interface abs

  ! Interface to multiplication operators with stellar luminosities objects as their second argument.
  interface operator(*)
     module procedure Stellar_Luminosities_Multiply_Switched
  end interface operator(*)
 
  !![
  <enumeration>
    <name>frame</name>
    <description>Frame for luminosity calculations.</description>
    <encodeFunction>yes</encodeFunction>
    <decodeFunction>yes</decodeFunction>
    <entry label="rest"    />
    <entry label="observed"/>
  </enumeration>
  !!]

  type stellarLuminosities
     !!{
     The stellar luminosities structure.
     !!}
     private
     double precision                      , allocatable, dimension(:) :: luminosityValue
     type            (enumerationFrameType)                            :: frame
   contains
     !![
     <methods>
       <method description="Multiply stellar luminosities by a scalar." method="operator(*)" />
       <method description="Divide stellar luminosities by a scalar." method="operator(/)" />
       <method description="Add two stellarLuminosities." method="operator(+)" />
       <method description="Subtract one abundance from another." method="operator(-)" />
       <method description="Increment a stellar luminosities object." method="increment" />
       <method description="Return a count of the number of properties in a serialized stellar luminosities object." method="serializeCount" />
       <method description="Serialize a stellar luminosities object to an array." method="serialize" />
       <method description="Deserialize a stellar luminosities object from an array." method="deserialize" />
       <method description="Return true if a stellar luminosities object is zero." method="isZero" />
       <method description="Destroy a stellar luminosities object." method="destroy" />
       <method description="Reset a stellar luminosities object." method="reset" />
       <method description="Build a stellar luminosities object from a provided XML description." method="builder" />
       <method description="Dump a stellar luminosities object." method="dump" />
       <method description="Dump a stellar luminosities object to binary." method="dumpRaw" />
       <method description="Read a stellar luminosities object from binary." method="readRaw" />
       <method description="Set a stellar luminosities object to unity." method="setToUnity" />
       <method description="Return the $i^\mathrm{th}$ luminosity." method="luminosity" />
       <method description="Store a stellar luminosities object in the output buffers." method="output" />
       <method description="Store a stellar luminosities object in the output buffers." method="postOutput" />
       <method description="Return the number of luminosities to be output at the given time." method="luminosityOutputCount" />
       <method description="Specify the count of a stellar luminosities object for output." method="outputCount" />
       <method description="Specify the names of stellar luminosities object properties for output." method="outputNames" />
       <method description="Return the total number of luminosities tracked. If {\normalfont \ttfamily unmapped} is true, then the number of luminosities prior to mapping is returned." method="luminosityCount" />
       <method description="Set the luminosities using a single stellar population." method="setLuminosities" />
       <method description="Return true if the indexed luminosity is to be output at the given time." method="isOutput" />
       <method description="Return the index to a luminosity specified by name or properties." method="index" />
       <method description="Return the name of a luminosity specified by index." method="name" />
       <method description="Truncate the number of stellar luminosities stored to match that in the given {\normalfont \ttfamily templateLuminosities}." method="truncate" />
       <method description="Returns the size of any non-static components of the type." method="nonStaticSizeOf" />
     </methods>
     !!]
     procedure         ::                          Stellar_Luminosities_Add
     procedure         ::                          Stellar_Luminosities_Subtract
     procedure         ::                          Stellar_Luminosities_Multiply
     procedure         ::                          Stellar_Luminosities_Divide
     generic           :: operator(+)           => Stellar_Luminosities_Add
     generic           :: operator(-)           => Stellar_Luminosities_Subtract
     generic           :: operator(*)           => Stellar_Luminosities_Multiply
     generic           :: operator(/)           => Stellar_Luminosities_Divide
     procedure         :: nonStaticSizeOf       => Stellar_Luminosities_Non_Static_Size_Of
     procedure         :: isZero                => Stellar_Luminosities_Is_Zero
     procedure         :: destroy               => Stellar_Luminosities_Destroy
     procedure         :: reset                 => Stellar_Luminosities_Reset
     procedure         :: builder               => Stellar_Luminosities_Builder
     procedure         :: dump                  => Stellar_Luminosities_Dump
     procedure         :: dumpRaw               => Stellar_Luminosities_Dump_Raw
     procedure         :: readRaw               => Stellar_Luminosities_Read_Raw
     procedure         :: setToUnity            => Stellar_Luminosities_Set_To_Unity
     procedure         :: luminosity            => Stellar_Luminosities_Luminosity
     procedure         :: setLuminosities       => Stellar_Luminosities_Set
     procedure, nopass :: luminosityCount       => Stellar_Luminosities_Property_Count
     procedure         :: serializeCount        => Stellar_Luminosities_Serialize_Count
     procedure         :: serialize             => Stellar_Luminosities_Serialize
     procedure         :: deserialize           => Stellar_Luminosities_Deserialize
     procedure         :: increment             => Stellar_Luminosities_Increment
     procedure         :: output                => Stellar_Luminosities_Output
     procedure         :: postOutput            => Stellar_Luminosities_Post_Output
     procedure, nopass :: luminosityOutputCount => Stellar_Luminosities_Output_Count_Get
     procedure         :: outputCount           => Stellar_Luminosities_Output_Count
     procedure         :: outputNames           => Stellar_Luminosities_Output_Names
     procedure, nopass :: isOutput              => Stellar_Luminosities_Is_Output
     procedure, nopass ::                          Stellar_Luminosities_Index_From_Name
     procedure, nopass ::                          Stellar_Luminosities_Index_From_Properties
     generic           :: index                 => Stellar_Luminosities_Index_From_Name      , &
          &                                        Stellar_Luminosities_Index_From_Properties
     procedure, nopass :: name                  => Stellar_Luminosities_Name
     procedure         :: truncate              => Stellar_Luminosities_Truncate
  end type stellarLuminosities
  
  ! Arrays which hold the luminosity specifications.
  integer                                                                                        :: luminosityCount                                      , luminosityCountUnmapped
  integer                                                            , allocatable, dimension(:) :: luminosityFilterIndex                                , luminosityIndex               , &
       &                                                                                            luminosityMap
  double precision                                                   , allocatable, dimension(:) :: luminosityBandRedshift                               , luminosityCosmicTime          , &
       &                                                                                            luminosityRedshift                                   , luminosityWavelengthEffective , &
       &                                                                                            luminosityVegaOffset
  type            (stellarPopulationSpectraPostprocessorList        ), allocatable, dimension(:) :: luminosityPostprocessor
  type            (varying_string                                   ), allocatable, dimension(:) :: luminosityFilter                                     , luminosityName                , &
       &                                                                                            luminosityPostprocessSet                             , luminosityType                , &
       &                                                                                            luminosityRedshiftText                               , luminosityBandRedshiftText

  ! Luminosity output options.
  integer                                                                                        :: luminosityOutputOption
  integer                                                            , parameter                 :: luminosityOutputOptionAll                    =0      , luminosityOutputOptionFuture=1, &
       &                                                                                            luminosityOutputOptionPresent                =2

  ! Unit and zero stellarLuminosities objects.
  type            (stellarLuminosities                              ), public                    :: unitStellarLuminosities                              , zeroStellarLuminosities

  ! Stellar population postprocessor builder used during initialization and state restoration.
  class           (stellarPopulationSpectraPostprocessorBuilderClass), pointer                   :: stellarPopulationSpectraPostprocessorBuilder__
  !$omp threadprivate(stellarPopulationSpectraPostprocessorBuilder__)

contains

  !![
  <nodeComponentInitializationTask>
   <unitName>Stellar_Luminosities_Initializor</unitName>
  </nodeComponentInitializationTask>
  !!]
  subroutine Stellar_Luminosities_Initializor(parameters)
    !!{
    Initialize the {\normalfont \ttfamily stellarLuminositiesStructure} object module. Determines which stellar luminosities are to be tracked.
    !!}
    use            :: Array_Utilities    , only : Array_Reverse
    use            :: Cosmology_Functions, only : cosmologyFunctions, cosmologyFunctionsClass
    use            :: Error              , only : Error_Report
    use, intrinsic :: ISO_C_Binding      , only : c_size_t
    use            :: ISO_Varying_String , only : assignment(=)     , char                       , operator(//)      , operator(/=), &
          &                                       operator(==)      , var_str
    use            :: Input_Parameters   , only : inputParameters
    use            :: Instruments_Filters, only : Filter_Get_Index  , Filter_Wavelength_Effective, Filter_Vega_Offset
    use            :: Sorting            , only : sortByIndex       , sortIndex
    use            :: String_Handling    , only : operator(//)
    implicit none
    type            (inputParameters                                  ), intent(inout)             :: parameters
    class           (cosmologyFunctionsClass                          ), pointer                   :: cosmologyFunctions_
    class           (stellarPopulationSpectraPostprocessorBuilderClass), pointer                   :: stellarPopulationSpectraPostprocessorBuilder_
    integer                                                                                        :: iLuminosity                                  , jLuminosity
    double precision                                                                               :: expansionFactor
    character       (len=10                                           )                            :: redshiftLabel
    type            (varying_string                                   )                            :: luminosityOutputOptionText
    integer         (c_size_t                                         ), allocatable, dimension(:) :: luminosityTimeIndex

    ! Get luminosity output option.
    !![
    <inputParameter>
      <name>luminosityOutputOption</name>
      <defaultValue>var_str('present')</defaultValue>
      <description>
         Selects which luminosities will be output at each output time:
         \begin{description}
         \item [all] Output all luminosities;
         \item [future] Output only those luminosities computed for the present output or future times;
         \item [present] Output only those luminosities computed for the present output time.
         \end{description}
      </description>
      <source>parameters</source>
      <variable>luminosityOutputOptionText</variable>
    </inputParameter>
    !!]
    select case (char(luminosityOutputOptionText))
    case ("all"    )
       luminosityOutputOption=luminosityOutputOptionAll
    case ("future" )
       luminosityOutputOption=luminosityOutputOptionFuture
    case ("present")
       luminosityOutputOption=luminosityOutputOptionPresent
    case default
       call Error_Report("unrecognized luminosityOutputOption"//{introspection:location})
    end select
    
    ! Get required objects.
    !![
    <objectBuilder class="cosmologyFunctions"                           name="cosmologyFunctions_"                           source="parameters"/>
    <objectBuilder class="stellarPopulationSpectraPostprocessorBuilder" name="stellarPopulationSpectraPostprocessorBuilder_" source="parameters"/>
    !!]
 
    ! Read in the parameters which specify the luminosities to be computed.
    luminosityCount=parameters%count('luminosityRedshift',zeroIfNotPresent=.true.)
    luminosityCountUnmapped=luminosityCount
    if (parameters%count('luminosityFilter',zeroIfNotPresent=.true.) /= luminosityCount) &
         & call Error_Report(var_str('luminosityFilter [')//parameters%count('luminosityFilter',zeroIfNotPresent=.true.)//'] and luminosityRedshift ['//luminosityCount//'] input arrays must have same dimension'//{introspection:location})
    if (parameters%count('luminosityType',zeroIfNotPresent=.true.) /= luminosityCount) &
         & call Error_Report(var_str('luminosityType [')//parameters%count('luminosityType',zeroIfNotPresent=.true.)//'] and luminosityRedshift ['//luminosityCount//'] input arrays must have same dimension'//{introspection:location})
    if (parameters%isPresent('luminosityBandRedshift')) then
       if (parameters%count('luminosityBandRedshift',zeroIfNotPresent=.true.) /= luminosityCount) &
            & call Error_Report('luminosityBandRedshift and luminosityRedshift input arrays must have same dimension'//{introspection:location})
    end if

    if (luminosityCount > 0) then
       if (allocated(luminosityMap             )) deallocate(luminosityMap             )
       if (allocated(luminosityRedshift        )) deallocate(luminosityRedshift        )
       if (allocated(luminosityBandRedshift    )) deallocate(luminosityBandRedshift    )
       if (allocated(luminosityFilter          )) deallocate(luminosityFilter          )
       if (allocated(luminosityType            )) deallocate(luminosityType            )
       if (allocated(luminosityPostprocessSet  )) deallocate(luminosityPostprocessSet  )
       if (allocated(luminosityRedshiftText    )) deallocate(luminosityRedshiftText    )       
       if (allocated(luminosityBandRedshiftText)) deallocate(luminosityBandRedshiftText)       
       allocate(luminosityMap             (luminosityCount))
       allocate(luminosityRedshift        (luminosityCount))
       allocate(luminosityBandRedshift    (luminosityCount))
       allocate(luminosityFilter          (luminosityCount))
       allocate(luminosityType            (luminosityCount))
       allocate(luminosityPostprocessSet  (luminosityCount))
       allocate(luminosityRedshiftText    (luminosityCount))
       allocate(luminosityBandRedshiftText(luminosityCount))
       !![
       <inputParameter>
         <name>luminosityRedshift</name>
         <description>The redshift for which to compute each specified stellar luminosity.</description>
         <source>parameters</source>
         <variable>luminosityRedshiftText</variable>
       </inputParameter>
       !!]
       do iLuminosity=1,size(luminosityRedshiftText)
          if (luminosityRedshiftText(iLuminosity) /= "all") then
             redshiftLabel=char(luminosityRedshiftText(iLuminosity))
             read (redshiftLabel,*) luminosityRedshift(iLuminosity)
          else
             luminosityRedshift(iLuminosity)=-2.0d0
          end if
          ! Assign a mapping from initial to final array of luminosities (this is initially an identity mapping).
          luminosityMap(iLuminosity)=iLuminosity
       end do
       if (parameters%isPresent('luminosityBandRedshift')) then
          !![
          <inputParameter>
            <name>luminosityBandRedshift</name>
	    <description>If present, force filters to be shifted to this redshift rather than that specified by {\normalfont \ttfamily [luminosityRedshift]}. Allows sampling of the SED at wavelengths corresponding to other redshifts.</description>
	    <source>parameters</source>
            <variable>luminosityBandRedshiftText</variable>
          </inputParameter>
          !!]
       do iLuminosity=1,size(luminosityBandRedshiftText)
          if (luminosityBandRedshiftText(iLuminosity) /= "all") then
             if (luminosityRedshiftText(iLuminosity) == "all") call Error_Report('luminosityBandRedshift=="all" only allowed where luminosityRedshift=="all"'//{introspection:location})
             redshiftLabel=char(luminosityBandRedshiftText(iLuminosity))
             read (redshiftLabel,*) luminosityBandRedshift(iLuminosity)
          else
             if (luminosityRedshiftText(iLuminosity) /= "all") call Error_Report('luminosityBandRedshift=="all" required where luminosityRedshift=="all"'    //{introspection:location})
             luminosityBandRedshift(iLuminosity)=-2.0d0
          end if
       end do
       else
          luminosityBandRedshift=luminosityRedshift
       end if
       !![
       <inputParameter>
         <name>luminosityFilter</name>
         <description>The filter name for each stellar luminosity to be computed.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>luminosityType</name>
         <description>
            The luminosity type for each stellar luminosity to be computed:
            \begin{description}
             \item[rest] Compute luminosity in the galaxy rest frame;
             \item[observed] Compute luminosity in the observer frame\footnote{The luminosity computed in this way is that in the galaxy rest
                             frame using a filter blueshifted to the galaxy's redshift. This means that to compute an apparent magnitude you
                             must add not only the distance modulus, but a factor of $-2.5\log_{10}(1+z)$ to account for compression of photon
                             frequencies.}.
            \end{description}
         </description>
         <source>parameters</source>
       </inputParameter>
       !!]
       ! Read postprocessing set information.
       if (parameters%count('luminosityPostprocessSet',zeroIfNotPresent=.true.) > 0) then
          if (parameters%count('luminosityPostprocessSet') /= luminosityCount) &
               & call Error_Report('luminosityPostprocessSet and luminosityFilter input arrays must have same dimension'//{introspection:location})
          !![
          <inputParameter>
            <name>luminosityPostprocessSet</name>
	    <description>The name of the set of postprocessing algorithms to apply to this filter.</description>
	    <source>parameters</source>
          </inputParameter>
          !!]
       else
          luminosityPostprocessSet="default"
       end if
       ! Handle luminosity definition special cases.
       call Stellar_Luminosities_Special_Cases(luminosityMap,luminosityRedshiftText,luminosityRedshift,luminosityBandRedshift,luminosityFilter,luminosityType,luminosityPostprocessSet,parameters)
       luminosityCount=size(luminosityRedshift)
       ! Allocate remaining required arrays.
       if (allocated(luminosityName               )) deallocate(luminosityName               )
       if (allocated(luminosityPostprocessor      )) deallocate(luminosityPostprocessor      )
       if (allocated(luminosityFilterIndex        )) deallocate(luminosityFilterIndex        )
       if (allocated(luminosityIndex              )) deallocate(luminosityIndex              )
       if (allocated(luminosityCosmicTime         )) deallocate(luminosityCosmicTime         )
       if (allocated(luminosityTimeIndex          )) deallocate(luminosityTimeIndex          )
       if (allocated(luminosityWavelengthEffective)) deallocate(luminosityWavelengthEffective)
       if (allocated(luminosityVegaOffset         )) deallocate(luminosityVegaOffset         )
       allocate(luminosityName               (luminosityCount))
       allocate(luminosityPostprocessor      (luminosityCount))
       allocate(luminosityFilterIndex        (luminosityCount))
       allocate(luminosityIndex              (luminosityCount))
       allocate(luminosityCosmicTime         (luminosityCount))
       allocate(luminosityTimeIndex          (luminosityCount))
       allocate(luminosityWavelengthEffective(luminosityCount))
       allocate(luminosityVegaOffset         (luminosityCount))
       ! Process the list of luminosities.
       do iLuminosity=1,luminosityCount
          ! Assign a name to this luminosity.
          write (redshiftLabel,'(f7.4)')     luminosityBandRedshift  (iLuminosity)
          luminosityName       (iLuminosity)=luminosityFilter        (iLuminosity)//":"// &
               &                             luminosityType          (iLuminosity)//":"// &
               &                             "z"   //trim(adjustl(redshiftLabel))
          if (parameters%isPresent('luminosityBandRedshift')) then
             write (redshiftLabel,'(f7.4)')  luminosityRedshift      (iLuminosity)
             luminosityName    (iLuminosity)=luminosityName          (iLuminosity)//":"// &
                  &                          "zOut"//trim(adjustl(redshiftLabel))
          end if
          if (luminosityPostprocessSet(iLuminosity) /= "default") &
               & luminosityName(iLuminosity)=luminosityName          (iLuminosity)//":"// &
               &                             luminosityPostprocessSet(iLuminosity)
          ! Check for duplicated luminosities.
          if (iLuminosity > 1) then
             do jLuminosity=1,iLuminosity-1
                if (luminosityName(iLuminosity) == luminosityName(jLuminosity)) &
                     & call Error_Report('luminosity '//luminosityName(iLuminosity)//' appears more than once in the input parameter file'//{introspection:location})
             end do
          end if
          ! Assign an index for this luminosity.
          luminosityIndex              (iLuminosity)=                                                                      iLuminosity
          ! Get the index of the specified filter.
          luminosityFilterIndex        (iLuminosity)=Filter_Get_Index                               (luminosityFilter     (iLuminosity))
          ! Get effective wavelength and Vega offset.
          luminosityWavelengthEffective(iLuminosity)=Filter_Wavelength_Effective                    (luminosityFilterIndex(iLuminosity))
          luminosityVegaOffset         (iLuminosity)=Filter_Vega_Offset                             (luminosityFilterIndex(iLuminosity))
          ! Set the reference time (i.e. cosmological time corresponding to the specified redshift) for this filter.
          expansionFactor                           =cosmologyFunctions_%expansionFactorFromRedshift(luminosityRedshift   (iLuminosity))
          luminosityCosmicTime         (iLuminosity)=cosmologyFunctions_%cosmicTime                 (expansionFactor                   )
          ! Set the filter redshifting factor. This is equal to the requested redshift if an observed frame was specified, otherwise
          ! it is set to zero to indicate a rest-frame filter.
          select case(char(luminosityType(iLuminosity)))
          case ("rest")
             luminosityBandRedshift(iLuminosity)=0.0d0
          case ("observed")
             ! Do nothing, we already have the correct redshift.
          case default
             call Error_Report('unrecognized filter type - must be "rest" or "observed"'//{introspection:location})
          end select
          ! Find the index for the postprocessing chain to be applied to this filter.
          luminosityPostprocessor(iLuminosity)%stellarPopulationSpectraPostprocessor_ => stellarPopulationSpectraPostprocessorBuilder_%build(luminosityPostprocessSet(iLuminosity))
       end do
       ! Sort the luminosities such that the latest luminosities are stored first.
       luminosityTimeIndex=Array_Reverse(sortIndex(luminosityCosmicTime))
       call sortByIndex             (luminosityFilterIndex        ,luminosityTimeIndex)
       call sortByIndexPostprocessor(luminosityPostprocessor      ,luminosityTimeIndex)
       call sortByIndex             (luminosityCosmicTime         ,luminosityTimeIndex)
       call sortByIndex             (luminosityName               ,luminosityTimeIndex)
       call sortByIndex             (luminosityRedshift           ,luminosityTimeIndex)
       call sortByIndex             (luminosityBandRedshift       ,luminosityTimeIndex)
       call sortByIndex             (luminosityFilter             ,luminosityTimeIndex)
       call sortByIndex             (luminosityType               ,luminosityTimeIndex)
       call sortByIndex             (luminosityPostprocessSet     ,luminosityTimeIndex)
       call sortByIndex             (luminosityWavelengthEffective,luminosityTimeIndex)
       call sortByIndex             (luminosityVegaOffset         ,luminosityTimeIndex)
       ! Allocate unit and zero stellar abundance objects.
       if (allocated(unitStellarLuminosities%luminosityValue)) deallocate(unitStellarLuminosities%luminosityValue)
       if (allocated(zeroStellarLuminosities%luminosityValue)) deallocate(zeroStellarLuminosities%luminosityValue)
       allocate(unitStellarLuminosities%luminosityValue(luminosityCount))
       allocate(zeroStellarLuminosities%luminosityValue(luminosityCount))
       unitStellarLuminosities%luminosityValue=1.0d0
       zeroStellarLuminosities%luminosityValue=0.0d0
    end if
    !![
    <objectDestructor name="cosmologyFunctions_"                          />
    <objectDestructor name="stellarPopulationSpectraPostprocessorBuilder_"/>
    !!]
    return
  end subroutine Stellar_Luminosities_Initializor
  
  !![
  <nodeComponentThreadInitializationTask>
   <unitName>Stellar_Luminosities_Thread_Initializor</unitName>
  </nodeComponentThreadInitializationTask>
  !!]
  subroutine Stellar_Luminosities_Thread_Initializor(parameters)
    use :: Input_Parameters, only : inputParameters
    implicit none
    type(inputParameters), intent(inout) :: parameters

    !![
    <objectBuilder class="stellarPopulationSpectraPostprocessorBuilder" name="stellarPopulationSpectraPostprocessorBuilder__" source="parameters"/>
    !!]
    return
  end subroutine Stellar_Luminosities_Thread_Initializor
  
  !![
  <nodeComponentThreadUninitializationTask>
   <unitName>Stellar_Luminosities_Thread_Uninitializor</unitName>
  </nodeComponentThreadUninitializationTask>
  !!]
  subroutine Stellar_Luminosities_Thread_Uninitializor()
    implicit none
    
    !![
    <objectDestructor name="stellarPopulationSpectraPostprocessorBuilder__"/>
    !!]
    return
  end subroutine Stellar_Luminosities_Thread_Uninitializor
  
  subroutine Stellar_Luminosities_Destroy(self)
    !!{
    Destroy an stellarLuminosities object.
    !!}
    implicit none
    class(stellarLuminosities), intent(inout) :: self

    if (allocated(self%luminosityValue)) deallocate(self%luminosityValue)
    return
  end subroutine Stellar_Luminosities_Destroy

  subroutine Stellar_Luminosities_Builder(self,stellarLuminositiesDefinition)
    !!{
    Build a {\normalfont \ttfamily stellarLuminosities} object from the given XML {\normalfont \ttfamily stellarLuminositiesDefinition}.
    !!}
    use :: FoX_DOM, only : node                        , extractDataContent
    use :: Error  , only : Error_Report
    use :: IO_XML , only : XML_Get_Elements_By_Tag_Name, xmlNodeList
    implicit none
    class  (stellarLuminosities), intent(inout)              :: self
    type   (node               ), intent(in   ), pointer     :: stellarLuminositiesDefinition
    type   (node               )               , pointer     :: luminosity
    type   (xmlNodeList        ), dimension(:) , allocatable :: luminosityList
    integer                                                  :: i

    ! Get the luminosities.
    !$omp critical (FoX_DOM_Access)
    call XML_Get_Elements_By_Tag_Name(stellarLuminositiesDefinition,'luminosity',luminosityList)
    !$omp end critical (FoX_DOM_Access)
    if (luminosityCount > 0) then
       do i=0,luminosityCount-1
          !$omp critical (FoX_DOM_Access)
          luminosity => luminosityList(i)%element
          call extractDataContent(luminosity,self%luminosityValue(i+1))
          !$omp end critical (FoX_DOM_Access)
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Builder

  subroutine Stellar_Luminosities_Dump(self,verbosityLevel)
    !!{
    Dump a stellar luminosities object.
    !!}
    use :: Display           , only : displayMessage, enumerationVerbosityLevelType
    use :: ISO_Varying_String, only : operator(//)
    implicit none
    class    (stellarLuminosities          ), intent(in   ) :: self
    type     (enumerationVerbosityLevelType), intent(in   ) :: verbosityLevel
    integer                                                 :: i
    character(len=22                       )                :: label
    type     (varying_string               )                :: message

    ! Dump the contents.
    if (luminosityCount > 0) then
       do i=1,luminosityCount
          if (i <= size(self%luminosityValue)) then
             write (label,'(e22.16)') self%luminosityValue(i)
          else
             label="pruned"
          end if
          message=luminosityName(i)//':          '//label
          call displayMessage(message,verbosityLevel)
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Dump

  subroutine Stellar_Luminosities_Dump_Raw(self,fileHandle)
    !!{
    Dump an stellarLuminosities object to binary.
    !!}
    implicit none
    class  (stellarLuminosities), intent(in   ) :: self
    integer                     , intent(in   ) :: fileHandle

    ! Dump the content.
    if (luminosityCount > 0) then
       write (fileHandle) size(self%luminosityValue)
       write (fileHandle) self%luminosityValue
    end if
    return
  end subroutine Stellar_Luminosities_Dump_Raw

  subroutine Stellar_Luminosities_Read_Raw(self,fileHandle)
    !!{
    Read an stellarLuminosities object from binary.
    !!}
    implicit none
    class  (stellarLuminosities), intent(inout) :: self
    integer                     , intent(in   ) :: fileHandle
    integer                                     :: luminosityActiveCount

    ! Read the content.
    if (luminosityCount > 0) then
       call Stellar_Luminosities_Create(self)
       read (fileHandle) luminosityActiveCount
       deallocate(self%luminosityValue                       )
       allocate  (self%luminosityValue(luminosityActiveCount))
       read (fileHandle) self%luminosityValue
    end if
    return
  end subroutine Stellar_Luminosities_Read_Raw

  subroutine Stellar_Luminosities_Reset(self)
    !!{
    Reset an stellarLuminosities object.
    !!}
    implicit none
    class(stellarLuminosities), intent(inout) :: self

    ! Ensure object is initialized.
    call Stellar_Luminosities_Create(self)
    ! Zero all properties.
    if (luminosityCount > 0) self%luminosityValue=0.0d0
    return
  end subroutine Stellar_Luminosities_Reset

  subroutine Stellar_Luminosities_Set_To_Unity(self)
    !!{
    Set an stellarLuminosities object to unity.
    !!}
    implicit none
    class(stellarLuminosities), intent(inout) :: self

    ! Ensure object is initialized.
    call Stellar_Luminosities_Create(self)
    ! Set values to unity.
    if (luminosityCount > 0) self%luminosityValue=1.0d0
    return
  end subroutine Stellar_Luminosities_Set_To_Unity

  logical function Stellar_Luminosities_Is_Zero(self)
    !!{
    Test whether an stellarLuminosities object is zero.
    !!}
    implicit none
    class(stellarLuminosities), intent(in   ) :: self

    ! Detect if all stellar luminosities are zero.
    Stellar_Luminosities_Is_Zero=.true.
    if (luminosityCount > 0 .and. allocated(self%luminosityValue)) then
       if (any(self%luminosityValue /= 0.0d0)) Stellar_Luminosities_Is_Zero=.false.
    end if
    return
  end function Stellar_Luminosities_Is_Zero

  double precision function Stellar_Luminosities_Luminosity(self,index)
    !!{
    Return the requested luminosity from a {\normalfont \ttfamily stellarLuminosities} object.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class  (stellarLuminosities), intent(inout) :: self
    integer                     , intent(in   ) :: index

    ! Return the requested luminosity.
    if (allocated(self%luminosityValue)) then
       if (index > 0 .and. index <= size(self%luminosityValue)) then
          Stellar_Luminosities_Luminosity=self%luminosityValue(index)
       else
          Stellar_Luminosities_Luminosity=0.0d0
          call Error_Report('index out of range'//{introspection:location})
       end if
    else
       Stellar_Luminosities_Luminosity=0.0d0
    end if
    return
  end function Stellar_Luminosities_Luminosity

  integer function stellarLuminositiesCountMaximum(luminosities1,luminosities2)
    implicit none
    type   (stellarLuminosities), intent(in   ) :: luminosities1, luminosities2

    if (allocated(luminosities1%luminosityValue).or.allocated(luminosities2%luminosityValue)) then
       stellarLuminositiesCountMaximum=0
       if (allocated(luminosities1%luminosityValue)) stellarLuminositiesCountMaximum=max(stellarLuminositiesCountMaximum,size(luminosities1%luminosityValue))
       if (allocated(luminosities2%luminosityValue)) stellarLuminositiesCountMaximum=max(stellarLuminositiesCountMaximum,size(luminosities2%luminosityValue))
    else
       stellarLuminositiesCountMaximum=luminosityCount
    end if
    return
  end function stellarLuminositiesCountMaximum

  function stellarLuminositiesMax(luminosities1,luminosities2)
    !!{
    Return an element-by-element {\normalfont \ttfamily max()} on two stellar luminosity objects.
    !!}
    implicit none
    type   (stellarLuminosities)                :: stellarLuminositiesMax
    type   (stellarLuminosities), intent(in   ) :: luminosities1         , luminosities2
    integer                                     :: luminosityCountActual

    if (luminosityCount > 0) then
       luminosityCountActual=stellarLuminositiesCountMaximum(luminosities1,luminosities2)
       stellarLuminositiesMax%luminosityValue=spread(-huge(0.0d0),1,luminosityCountActual)
       if (allocated(luminosities1%luminosityValue))                                              &
            & stellarLuminositiesMax     %luminosityValue(1:size(luminosities1%luminosityValue))= &
            &  max(                                                                               &
            &      stellarLuminositiesMax%luminosityValue(1:size(luminosities1%luminosityValue)), &
            &      luminosities1         %luminosityValue(1:size(luminosities1%luminosityValue))  &
            &     )
       if (allocated(luminosities2%luminosityValue))                                              &
            & stellarLuminositiesMax     %luminosityValue(1:size(luminosities2%luminosityValue))= &
            &  max(                                                                               &
            &      stellarLuminositiesMax%luminosityValue(1:size(luminosities2%luminosityValue)), &
            &      luminosities2         %luminosityValue(1:size(luminosities2%luminosityValue))  &
            &     )
    end if
    return
  end function stellarLuminositiesMax

  function stellarLuminositiesAbs(luminosities)
    !!{
    Return an element-by-element {\normalfont \ttfamily abs()} on a stellar luminosity object.
    !!}
    implicit none
    type(stellarLuminosities)                :: stellarLuminositiesAbs
    type(stellarLuminosities), intent(in   ) :: luminosities

    if (luminosityCount > 0) stellarLuminositiesAbs%luminosityValue=abs(luminosities%luminosityValue)
    return
  end function stellarLuminositiesAbs

  function Stellar_Luminosities_Add(luminosities1,luminosities2)
    !!{
    Add two stellar luminosities objects.
    !!}
    implicit none
    type   (stellarLuminosities)                          :: Stellar_Luminosities_Add
    class  (stellarLuminosities), intent(in   )           :: luminosities1
    class  (stellarLuminosities), intent(in   ), optional :: luminosities2
    integer                                               :: luminosityCountActual

    if (luminosityCount > 0) then
       if (present(luminosities2)) then
          luminosityCountActual=stellarLuminositiesCountMaximum(luminosities1,luminosities2)
          Stellar_Luminosities_Add%luminosityValue=spread(0.0d0,1,luminosityCountActual)
          if (allocated(luminosities1%luminosityValue))                                                 &
               & Stellar_Luminosities_Add      %luminosityValue(1:size(luminosities1%luminosityValue))= &
               &      +Stellar_Luminosities_Add%luminosityValue(1:size(luminosities1%luminosityValue))  &
               &      +luminosities1           %luminosityValue(1:size(luminosities1%luminosityValue))
          if (allocated(luminosities2%luminosityValue))                                                 &
               & Stellar_Luminosities_Add      %luminosityValue(1:size(luminosities2%luminosityValue))= &
               &      +Stellar_Luminosities_Add%luminosityValue(1:size(luminosities2%luminosityValue))  &
               &      +luminosities2           %luminosityValue(1:size(luminosities2%luminosityValue))
       else
          Stellar_Luminosities_Add%luminosityValue=+luminosities1%luminosityValue
       end if
    end if
    return
  end function Stellar_Luminosities_Add

  subroutine Stellar_Luminosities_Increment(self,increment)
    !!{
    Increment a stellar luminosities object.
    !!}
    implicit none
    class  (stellarLuminosities), intent(inout) :: self
    class  (stellarLuminosities), intent(in   ) :: increment
    integer                                     :: luminosityCountActual

    ! Increment.
    if (luminosityCount > 0) then
       luminosityCountActual=luminosityCount
       if (allocated(self     %luminosityValue)) luminosityCountActual=min(luminosityCountActual,size(self     %luminosityValue))
       if (allocated(increment%luminosityValue)) luminosityCountActual=min(luminosityCountActual,size(increment%luminosityValue))
       self%luminosityValue(1:luminosityCountActual)=+self     %luminosityValue(1:luminosityCountActual) &
            &                                        +increment%luminosityValue(1:luminosityCountActual)
    end if
    return
  end subroutine Stellar_Luminosities_Increment

  function Stellar_Luminosities_Subtract(luminosities1,luminosities2)
    !!{
    Subtract two stellar luminosities objects.
    !!}
    implicit none
    type   (stellarLuminosities)                          :: Stellar_Luminosities_Subtract
    class  (stellarLuminosities), intent(in   )           :: luminosities1
    class  (stellarLuminosities), intent(in   ), optional :: luminosities2
    integer                                               :: luminosityCountActual

    if (luminosityCount > 0) then
       if (present(luminosities2)) then
          luminosityCountActual=stellarLuminositiesCountMaximum(luminosities1,luminosities2)
          Stellar_Luminosities_Subtract%luminosityValue=spread(0.0d0,1,luminosityCountActual)
          if (allocated(luminosities1%luminosityValue))                                                      &
               & Stellar_Luminosities_Subtract      %luminosityValue(1:size(luminosities1%luminosityValue))= &
               &      +Stellar_Luminosities_Subtract%luminosityValue(1:size(luminosities1%luminosityValue))  &
               &      +luminosities1                %luminosityValue(1:size(luminosities1%luminosityValue))
          if (allocated(luminosities2%luminosityValue))                                                      &
               & Stellar_Luminosities_Subtract      %luminosityValue(1:size(luminosities2%luminosityValue))= &
               &      +Stellar_Luminosities_Subtract%luminosityValue(1:size(luminosities2%luminosityValue))  &
               &      -luminosities2                %luminosityValue(1:size(luminosities2%luminosityValue))
       else
          Stellar_Luminosities_Subtract%luminosityValue=-luminosities1%luminosityValue
       end if
    end if
    return
  end function Stellar_Luminosities_Subtract

  function Stellar_Luminosities_Multiply(stellarLuminosities1,multiplier)
    !!{
    Multiply a stellar luminosities object by a scalar.
    !!}
    implicit none
    type            (stellarLuminosities)                :: Stellar_Luminosities_Multiply
    class           (stellarLuminosities), intent(in   ) :: stellarLuminosities1
    double precision                     , intent(in   ) :: multiplier

    if (luminosityCount > 0) Stellar_Luminosities_Multiply%luminosityValue=stellarLuminosities1%luminosityValue*multiplier
    return
  end function Stellar_Luminosities_Multiply

  function Stellar_Luminosities_Multiply_Switched(multiplier,stellarLuminosities1)
    !!{
    Multiply a stellar luminosities object by a scalar.
    !!}
    implicit none
    type            (stellarLuminosities)                :: Stellar_Luminosities_Multiply_Switched
    class           (stellarLuminosities), intent(in   ) :: stellarLuminosities1
    double precision                     , intent(in   ) :: multiplier

    if (luminosityCount > 0) Stellar_Luminosities_Multiply_Switched%luminosityValue=stellarLuminosities1%luminosityValue*multiplier
    return
  end function Stellar_Luminosities_Multiply_Switched

  function Stellar_Luminosities_Divide(stellarLuminosities1,divisor)
    !!{
    Divide a stellar luminosities object by a scalar.
    !!}
    implicit none
    type            (stellarLuminosities)                :: Stellar_Luminosities_Divide
    class           (stellarLuminosities), intent(in   ) :: stellarLuminosities1
    double precision                     , intent(in   ) :: divisor

    if (luminosityCount > 0) Stellar_Luminosities_Divide%luminosityValue=stellarLuminosities1%luminosityValue/divisor
    return
  end function Stellar_Luminosities_Divide

  integer function Stellar_Luminosities_Property_Count(unmapped)
    !!{
    Return the number of properties required to track stellar luminosities.
    !!}
    implicit none
    logical, intent(in   ), optional :: unmapped

    ! Return the relevant count.
    if (present(unmapped).and.unmapped) then
       Stellar_Luminosities_Property_Count=luminosityCountUnmapped
    else
       Stellar_Luminosities_Property_Count=luminosityCount
    end if
    return
  end function Stellar_Luminosities_Property_Count

  integer function Stellar_Luminosities_Serialize_Count(self)
    !!{
    Return the number of properties required to track stellar luminosities.
    !!}
    implicit none
    class(stellarLuminosities), intent(in   ) :: self

    if (allocated(self%luminosityValue)) then
       Stellar_Luminosities_Serialize_Count=size(self%luminosityValue)
    else
       Stellar_Luminosities_Serialize_Count=luminosityCount
    end if
    return
  end function Stellar_Luminosities_Serialize_Count

  function Stellar_Luminosities_Name(index)
    !!{
    Return a name for the specified entry in the stellar luminosities structure.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : trim
    implicit none
    type   (varying_string)                :: Stellar_Luminosities_Name
    integer                , intent(in   ) :: index

    ! Check for index in range.
    if (index > 0 .and. index <= luminosityCount) then
       Stellar_Luminosities_Name=trim(luminosityName(index))
    else
       call Error_Report('index out of range'//{introspection:location})
    end if
    return
  end function Stellar_Luminosities_Name

  subroutine Stellar_Luminosities_Create(self)
    !!{
    Ensure that the {\normalfont \ttfamily luminosity} array in a {\normalfont \ttfamily stellarLuminosities} is allocated.
    !!}
    implicit none
    type(stellarLuminosities), intent(inout) :: self

    if (.not.allocated(self%luminosityValue)) allocate(self%luminosityValue(luminosityCount))
    return
  end subroutine Stellar_Luminosities_Create

  subroutine Stellar_Luminosities_Deserialize(self,stellarLuminositiesArray)
    !!{
    Pack stellar luminosities from an array into a {\normalfont \ttfamily stellarLuminosities} structure.
    !!}
    implicit none
    class           (stellarLuminosities)              , intent(inout) :: self
    double precision                     , dimension(:), intent(in   ) :: stellarLuminositiesArray

    select type (self)
    type is (stellarLuminosities)
       ! Ensure luminosities array exists.
       call Stellar_Luminosities_Create(self)
       ! Extract luminosity values from array.
       self%luminosityValue=stellarLuminositiesArray(1:size(self%luminosityValue))
    end select
    return
  end subroutine Stellar_Luminosities_Deserialize

  subroutine Stellar_Luminosities_Serialize(self,stellarLuminositiesArray)
    !!{
    Unpack stellar luminosities from a {\normalfont \ttfamily stellarLuminosities} structure into an array.
    !!}
    implicit none
    double precision                     , dimension(:), intent(  out) :: stellarLuminositiesArray(:)
    class           (stellarLuminosities)              , intent(in   ) :: self

    ! Place luminosities into array.
    if (allocated(self%luminosityValue)) then
       stellarLuminositiesArray(1:size(self%luminosityValue))=self%luminosityValue
    else
       stellarLuminositiesArray(1:luminosityCount           )=0.0d0
    end if
    return
  end subroutine Stellar_Luminosities_Serialize

  subroutine Stellar_Luminosities_Output(self,integerProperty,integerBufferCount,integerProperties,doubleProperty,doubleBufferCount,doubleProperties,time,outputInstance)
    !!{
    Store a {\normalfont \ttfamily stellarLuminosities} object in the output buffers.
    !!}
    use :: Kind_Numbers                      , only : kind_int8
    use :: Multi_Counters                    , only : multiCounter
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger , outputPropertyDouble
    implicit none
    class           (stellarLuminosities  ), intent(inout)               :: self
    double precision                       , intent(in   )               :: time
    integer                                , intent(inout)               :: doubleBufferCount , doubleProperty , &
         &                                                                  integerBufferCount, integerProperty
    type            (outputPropertyInteger), intent(inout), dimension(:) :: integerProperties
    type            (outputPropertyDouble ), intent(inout), dimension(:) :: doubleProperties
    type            (multiCounter         ), intent(in   )               :: outputInstance
    integer                                                              :: i
    !$GLC attributes unused :: integerProperty, integerBufferCount, integerProperties, outputInstance

    if (luminosityCount > 0) then
       do i=1,luminosityCount
          if (Stellar_Luminosities_Is_Output(i,time)) then
             doubleProperties(doubleProperty+1)%scalar(doubleBufferCount)=self%luminosityValue(i)
             doubleProperty=doubleProperty+1
          end if
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Output

  subroutine Stellar_Luminosities_Post_Output(self,time)
    !!{
    Clean up a {\normalfont \ttfamily stellarLuminosities} object after output.
    !!}
    implicit none
    class           (stellarLuminosities)                , intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    double precision                     , dimension(:  ), allocatable   :: luminosityTmp
    integer                                                              :: i            , luminosityRemainingCount

    if (luminosityCount > 0) then
       select case (luminosityOutputOption)
       case (luminosityOutputOptionFuture,luminosityOutputOptionPresent)
          ! Luminosities from this and earlier outputs no longer needed, so prune them.
          call Move_Alloc(self%luminosityValue,luminosityTmp)
          luminosityRemainingCount=luminosityCount
          do i=1,luminosityCount
             if (.not.Stellar_Luminosities_Is_Output(i,time,luminosityOutputOptionFuture)) &
                  & luminosityRemainingCount=luminosityRemainingCount-1
          end do
          allocate(self%luminosityValue(luminosityRemainingCount))
          self%luminosityValue=luminosityTmp(1:luminosityRemainingCount)
          deallocate(luminosityTmp)
       end select
    end if
    return
  end subroutine Stellar_Luminosities_Post_Output

  subroutine Stellar_Luminosities_Output_Count(self,integerPropertyCount,doublePropertyCount,time)
    !!{
    Increment the output count to account for a {\normalfont \ttfamily stellarLuminosities} object.
    !!}
    implicit none
    class           (stellarLuminosities), intent(in   ) :: self
    integer                              , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision                     , intent(in   ) :: time
    !$GLC attributes unused :: self, integerPropertyCount

    doublePropertyCount=doublePropertyCount+Stellar_Luminosities_Output_Count_Get(time)
    return
  end subroutine Stellar_Luminosities_Output_Count

  integer function Stellar_Luminosities_Output_Count_Get(time)
    !!{
    Compute the number of luminosities to be output at a given time.
    !!}
    implicit none
    double precision, intent(in   ) :: time
    integer                         :: i

    Stellar_Luminosities_Output_Count_Get=0
    do i=1,luminosityCount
       if (Stellar_Luminosities_Is_Output(i,time)) Stellar_Luminosities_Output_Count_Get=Stellar_Luminosities_Output_Count_Get+1
    end do
    return
  end function Stellar_Luminosities_Output_Count_Get

  subroutine Stellar_Luminosities_Output_Names(self,integerProperty,integerProperties,doubleProperty,doubleProperties,time,prefix,comment,unitsInSI)
    !!{
    Assign names to output buffers for a {\normalfont \ttfamily stellarLuminosities} object.
    !!}
    use :: ISO_Varying_String                , only : assignment(=)        , operator(//)        , trim
    use :: Merger_Tree_Outputter_Buffer_Types, only : outputPropertyInteger, outputPropertyDouble
    implicit none
    class           (stellarLuminosities  )              , intent(in   ) :: self
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleProperty   , integerProperty
    type            (outputPropertyInteger), dimension(:), intent(inout) :: integerProperties
    type            (outputPropertyDouble ), dimension(:), intent(inout) :: doubleProperties
    character       (len=*                )              , intent(in   ) :: comment          , prefix
    double precision                                     , intent(in   ) :: unitsInSI
    integer                                                              :: i
    !$GLC attributes unused :: self, integerProperty, integerProperties

    if (luminosityCount > 0) then
       do i=1,luminosityCount
          if (Stellar_Luminosities_Is_Output(i,time)) then
             doubleProperty=doubleProperty+1
             doubleProperties(doubleProperty)%name     =trim(prefix )// ':'//trim(luminosityName(i))
             doubleProperties(doubleProperty)%comment  =trim(comment)//' ['//trim(luminosityName(i))//']'
             doubleProperties(doubleProperty)%unitsInSI=unitsInSI
             call doubleProperties(doubleProperty)%metaDataRank0%set('wavelengthEffective',luminosityWavelengthEffective(i))
             call doubleProperties(doubleProperty)%metaDataRank0%set('vegaOffset'         ,luminosityVegaOffset         (i))
          end if
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Output_Names

  logical function Stellar_Luminosities_Is_Output(luminosityIndex,time,outputOption)
    !!{
    Return true or false depending on whether {\normalfont \ttfamily luminosityIndex} should be output at {\normalfont \ttfamily time}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer         , intent(in   )           :: luminosityIndex
    double precision, intent(in   )           :: time
    integer         , intent(in   ), optional :: outputOption
    double precision, parameter               :: timeTolerance     =1.0d-3
    integer                                   :: outputOptionActual

    ! Determine output option to use.
    if (present(outputOption)) then
       outputOptionActual=outputOption
    else
       outputOptionActual=luminosityOutputOption
    end if
    select case (outputOptionActual)
    case (luminosityOutputOptionAll)
       Stellar_Luminosities_Is_Output=.true.
    case (luminosityOutputOptionFuture)
       Stellar_Luminosities_Is_Output=(    luminosityCosmicTime(luminosityIndex)       >= time*(1.0d0-timeTolerance))
    case (luminosityOutputOptionPresent)
       Stellar_Luminosities_Is_Output=(abs(luminosityCosmicTime(luminosityIndex)-time) <= time*       timeTolerance )
    case default
       Stellar_Luminosities_Is_Output=.false.
       call Error_Report('unknown luminosity output option'//{introspection:location})
    end select
    return
  end function Stellar_Luminosities_Is_Output

  subroutine Stellar_Luminosities_Set(self,mass,stellarPopulation_,stellarPopulationBroadBandLuminosities_,time,abundancesStellar)
    !!{
    Set the luminosity in each band for a single {\normalfont \ttfamily stellarPopulation\_} of given {\normalfont \ttfamily
    mass} with the specified {\normalfont \ttfamily abundancesStellar} and which formed at cosmological {\normalfont \ttfamily
    time}.
    !!}
    use :: Abundances_Structure                      , only : abundances
    use :: Stellar_Population_Broad_Band_Luminosities, only : stellarPopulationBroadBandLuminositiesClass
    use :: Stellar_Populations                       , only : stellarPopulationClass
    implicit none
    class           (stellarLuminosities                        )                             :: self
    class           (stellarPopulationClass                     ), intent(inout)              :: stellarPopulation_
    class           (stellarPopulationBroadBandLuminositiesClass), intent(inout)              :: stellarPopulationBroadBandLuminosities_
    double precision                                             , intent(in   )              :: mass                                   , time
    type            (abundances                                 ), intent(in   )              :: abundancesStellar
    double precision                                             , dimension(:) , allocatable :: ages                                   , massToLightRatio

    ! Return if no luminosities are tracked.
    if (luminosityCount == 0) return

    ! Allocate workspace.
    allocate(ages            (luminosityCount))
    allocate(massToLightRatio(luminosityCount))

    ! Get the ages that this stellar population will have at the various output times.
    ages=luminosityCosmicTime-time

    ! Get the luminosities for each requested band.
    massToLightRatio=stellarPopulationBroadBandLuminosities_%luminosities(                         &
         &                                                                luminosityIndex        , &
         &                                                                luminosityFilterIndex  , &
         &                                                                luminosityPostprocessor, &
         &                                                                stellarPopulation_     , &
         &                                                                abundancesStellar      , &
         &                                                                ages                   , &
         &                                                                luminosityBandRedshift   &
         &                                                               )
    call Stellar_Luminosities_Create(self)
    self%luminosityValue=mass*massToLightRatio(1:size(self%luminosityValue))
    return
  end subroutine Stellar_Luminosities_Set

  integer function Stellar_Luminosities_Index_From_Name(name)
    !!{
    Return the index of and specified entry in the luminosity list given its name.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(==)
    implicit none
    type   (varying_string), intent(in   ) :: name
    integer                                :: i

    Stellar_Luminosities_Index_From_Name=-1
    do i=1,luminosityCount
       if (name == luminosityName(i)) then
          Stellar_Luminosities_Index_From_Name=i
          return
       end if
    end do
    call Error_Report('unmatched name'//{introspection:location})
    return
  end function Stellar_Luminosities_Index_From_Name

  integer function Stellar_Luminosities_Index_From_Properties(filterName,filterType,redshift,redshiftBand,postprocessChain)
    !!{
    Return the index of and specified entry in the luminosity list given its properties.
    !!}
    use :: Display             , only : displayReset , displayGreen
    use :: Error               , only : Error_Report
    use :: ISO_Varying_String  , only : assignment(=), operator(//)   , operator(==), varying_string, &
         &                              len
    use :: Numerical_Comparison, only : Values_Agree
    use :: String_Handling     , only : operator(//) , stringXMLFormat
    implicit none
    character       (len=*         ), intent(in   )           :: filterName             , filterType
    double precision                , intent(in   )           :: redshift
    double precision                , intent(in   ), optional :: redshiftBand
    character       (len=*         ), intent(in   ), optional :: postprocessChain
    integer                                                   :: i                      , lengthNameMaximum      , &
         &                                                       countDigitsMaximum     , lengthDescriptorMaximum, &
         &                                                       countCountDigitsMaximum
    character       (len=64        )                          :: label                  , digitFormat            , &
         &                                                       labelCount             , digitDigitFormat
    type            (varying_string)                          :: message

    Stellar_Luminosities_Index_From_Properties=-1
    do i=1,luminosityCount
       if     (                                                                                             &
            &                  filterName       == luminosityFilter        (i)                              &
            &  .and.                                                                                        &
            &                  filterType       == luminosityType          (i)                              &
            &  .and.                                                                                        &
            &     Values_Agree(redshift    ,       luminosityRedshift      (i),relTol=1.0d-6,absTol=1.0d-6) &
            &  .and.                                                                                        &
            &   (                                                                                           &
            &     .not.present(redshiftBand    )                                                            &
            &    .or.                                                                                       &
            &     Values_Agree(redshiftBand,       luminosityBandRedshift  (i),relTol=1.0d-6,absTol=1.0d-6) &
            &   )                                                                                           &
            &  .and.                                                                                        &
            &   (                                                                                           &
            &     .not.present(postprocessChain)                                                            &
            &    .or.                                                                                       &
            &                  postprocessChain == luminosityPostprocessSet(i)                              &
            &   )                                                                                           &
            & ) then
          Stellar_Luminosities_Index_From_Properties=i
          return
       end if
    end do
    if (luminosityCount > 0) then
       lengthNameMaximum =max(4,max(len(filterName),maxval(len(luminosityFilter))))
       countDigitsMaximum=int(log10(dble(luminosityCount)))+1
    else
       lengthNameMaximum =4
       countDigitsMaximum=1
    end if
    countCountDigitsMaximum=int(log10(dble(countDigitsMaximum)))+1
    lengthDescriptorMaximum=max(9,5+2*countDigitsMaximum)
    write (digitDigitFormat,'(a,i1,a,i1,a)') '(a,i',countCountDigitsMaximum,',a,i',countCountDigitsMaximum,',a)'
    write (digitFormat,digitDigitFormat) "(i",countDigitsMaximum,".",countDigitsMaximum,")"
    write (label,'(f7.4)') redshift
    message='unmatched properties for requested stellar luminosity:'//char(10)
    message=message//char(10)//repeat(" ",lengthDescriptorMaximum+2)//'name'    //repeat(" ",lengthNameMaximum-4              )//' '//'type    '//' '//'redshift'
    if (present(redshiftBand)) then
       message=message//' '//'band-redshift'
    end if
    if (present(postprocessChain)) then
       message=message//' '//'postprocessor'
    end if
    message=message//char(10)//'requested:'//repeat(" ",lengthDescriptorMaximum-9)//' '//filterName//repeat(" ",lengthNameMaximum-len(filterName))//' '//filterType//repeat(" ",8-len(filterType))//' '//trim(adjustl(label))//' '
    if (present(redshiftBand)) then
       write (label,'(f7.4)') redshiftBand
       message=message//' '//trim(adjustl(label))//'       '
    end if
    if (present(postprocessChain)) then
       message=message//' '//postprocessChain
    end if
    if (luminosityCount > 0) then
       message=message//char(10)//'available:'
       write (labelCount,digitFormat) luminosityCount
       do i=1,luminosityCount
          write (label,digitFormat) i
          message=message//char(10)//" "//trim(label)//" of "//trim(labelCount)//repeat(" ",lengthDescriptorMaximum-4-2*countCountDigitsMaximum)
          message=message//" "//luminosityFilter(i)//repeat(" ",lengthNameMaximum-len(luminosityFilter(i)))
          message=message//" "//luminosityType  (i)//repeat(" ",8                -len(luminosityType  (i)))
          write    (label,'(f7.4)') luminosityRedshift    (i)
          message   =message//" "//trim(adjustl(label))//" "
          if (present(redshiftBand)) then
             write (label,'(f7.4)') luminosityBandRedshift(i)
             message=message//" "//trim(adjustl(label))//"      "
          end if
          if (present(postprocessChain)) then
             message=message//" "//luminosityPostprocessSet(i)
          end if
       end do
    else
       message=message//char(10)//'no luminosities available'
    end if
    message=message//char(10)//char(10)//displayGreen()//'HELP:'//displayReset()//' you can resolve this by adding the appropriate entries (highlighted in bold) to the following parameters:'//char(10)
    message   =message//char(10)//stringXMLFormat('<luminosityFilter         value="'//             filterName        //'"/>')
    message   =message//char(10)//stringXMLFormat('<luminosityType           value="'//             filterType        //'"/>')
    write (label,'(f7.4)') redshift
    message   =message//char(10)//stringXMLFormat('<luminosityRedshift       value="'//trim(adjustl(label           ))//'"/>')
    if (present(redshiftBand)) then
       write (label,'(f7.4)') redshiftBand
       message=message//char(10)//stringXMLFormat('<luminosityBandRedshift   value="'//trim(adjustl(label           ))//'"/>')
    end if
    if (present(postprocessChain)) then
       message=message//char(10)//stringXMLFormat('<luminosityPostprocessSet value="'//             postprocessChain  //'"/>')
    end if
    call Error_Report(message//{introspection:location})
    return
  end function Stellar_Luminosities_Index_From_Properties

  subroutine Stellar_Luminosities_SED_Top_Hat_Step(wavelengthCentral,filterWidth,wavelengthMinimum,wavelengthMaximum,observedWidth,redshift,stellarPopulationSpectra_)
    !!{
    Given a top hat filter central wavelength and filter width, determine the position and width of the next top hat filter in the array.
    !!}
    use :: Stellar_Population_Spectra, only : stellarPopulationSpectraClass
    implicit none
    double precision                               , intent(inout)          :: wavelengthCentral        , filterWidth
    double precision                               , intent(in   )          :: wavelengthMinimum        , wavelengthMaximum    , &
         &                                                                     observedWidth            , redshift
    class           (stellarPopulationSpectraClass), intent(in   ), pointer :: stellarPopulationSpectra_
    double precision                                                        :: restWavelengthMinimum    , restWavelengthMaximum, &
         &                                                                     restWidth                , wavelengthLowerEdge  , &
         &                                                                     tabulatedWidth

    ! Determine rest-frame wavelength extent.
    restWavelengthMinimum=wavelengthMinimum/(1.0d0+redshift)
    restWavelengthMaximum=wavelengthMaximum/(1.0d0+redshift)
    restWidth            =observedWidth    /(1.0d0+redshift)
    ! Move to lower edge of next filter.
    wavelengthLowerEdge=wavelengthCentral+filterWidth/2.0d0
    ! Get wavelength interval in stellar population spectra at lower edge wavelength
    tabulatedWidth=stellarPopulationSpectra_%wavelengthInterval(wavelengthLowerEdge)
    ! Determine where the new filter is: (i) inside the observed wavelength range, (ii) inside the rest wavelength range,
    ! or (iii) in between the observed and rest wavelength ranges.
    if      (wavelengthLowerEdge < restWavelengthMaximum) then
       ! Option (i): still inside rest-frame wavelength range.
       if (tabulatedWidth > restWidth) filterWidth=tabulatedWidth
       wavelengthCentral=wavelengthLowerEdge+filterWidth/2.0d0
    else if (wavelengthLowerEdge > wavelengthMinimum    ) then
       ! Option (ii): inside observed-frame wavelength range.
       filterWidth=observedWidth
       if(tabulatedWidth > observedWidth) filterWidth=tabulatedWidth
       wavelengthCentral=wavelengthLowerEdge+filterWidth/2.0d0
    else
       ! Option (iii): between rest-frame and observed-frame wavelength ranges.
       wavelengthCentral=wavelengthMinimum
       filterWidth      =observedWidth
       tabulatedWidth   =stellarPopulationSpectra_%wavelengthInterval(wavelengthCentral)
       if(tabulatedWidth > observedWidth) filterWidth=tabulatedWidth
       ! If the gap between the rest-frame and observed-frame ranges is small, check to avoid overlap of filters.
       if (wavelengthCentral-filterWidth/2.0d0 < wavelengthLowerEdge) then
          ! Overlap possible. Adjust filter position.
          tabulatedWidth=stellarPopulationSpectra_%wavelengthInterval(wavelengthLowerEdge)
          filterWidth   =observedWidth
          if(tabulatedWidth > observedWidth) filterWidth=tabulatedWidth
          wavelengthCentral=wavelengthLowerEdge+filterWidth/2.0d0
       end if
    end if
  end subroutine Stellar_Luminosities_SED_Top_Hat_Step

  subroutine Stellar_Luminosities_Special_Cases(luminosityMap,luminosityRedshiftText,luminosityRedshift,luminosityBandRedshift,luminosityFilter,luminosityType,luminosityPostprocessSet,parameters)
    !!{
    Modify the input list of luminosities for special cases.
    !!}
    use            :: Cosmology_Functions       , only : cosmologyFunctions      , cosmologyFunctionsClass
    use            :: HII_Region_Emission_Lines , only : emissionLineWavelength
    use, intrinsic :: ISO_C_Binding             , only : c_size_t
    use            :: ISO_Varying_String        , only : assignment(=)           , char                         , extract, operator(==), &
          &                                              var_str
    use            :: Input_Parameters          , only : inputParameters
    use            :: Stellar_Luminosities_Data , only : outputCount             , outputRedshifts
    use            :: Stellar_Population_Spectra, only : stellarPopulationSpectra, stellarPopulationSpectraClass
    use            :: String_Handling           , only : String_Split_Words      , char
    implicit none
    integer                                        , intent(inout), allocatable, dimension(:) :: luminosityMap
    type            (varying_string               ), intent(inout), allocatable, dimension(:) :: luminosityRedshiftText   , luminosityFilter           , &
         &                                                                                       luminosityType           , luminosityPostprocessSet
    double precision                               , intent(inout), allocatable, dimension(:) :: luminosityRedshift       , luminosityBandRedshift
    type            (inputParameters              ), intent(inout)                            :: parameters
    integer         (c_size_t                     )                                           :: i                        , j                          , &
         &                                                                                       k                        , newFilterCount,luminosityCount
    integer                                                       , allocatable, dimension(:) :: luminosityMapTmp
    type            (varying_string               )               , allocatable, dimension(:) :: luminosityRedshiftTextTmp, luminosityFilterTmp        , &
         &                                                                                       luminosityTypeTmp        , luminosityPostprocessSetTmp
    type            (varying_string               )                            , dimension(5) :: specialFilterWords
    double precision                                              , allocatable, dimension(:) :: luminosityRedshiftTmp    , luminosityBandRedshiftTmp
    class           (stellarPopulationSpectraClass), pointer                                  :: stellarPopulationSpectra_
    class           (cosmologyFunctionsClass      ), pointer                                  :: cosmologyFunctions_
    character       (len= 32                      )                                           :: redshiftLabel            , word                       , &
         &                                                                                       wavelengthCentralLabel   , resolutionLabel
    character       (len=256                      )                                           :: newFilterName            , lineName
    double precision                                                                          :: wavelengthMinimum        , wavelengthMaximum          , &
         &                                                                                       restWavelengthMinimum    , restWavelengthMaximum      , &
         &                                                                                       wavelengthRatio          , wavelengthCentral          , &
         &                                                                                       observedWidth            , restWidth                  , &
         &                                                                                       tabulatedWidth           , filterWidth                , &
         &                                                                                       resolution
    
    !![
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="stellarPopulationSpectra" name="stellarPopulationSpectra_" source="parameters"/>
    !!]
    ! Iterate over all luminosities.
    i=1
    do while (i <= size(luminosityRedshiftText))
       ! Check for special cases.
       luminosityCount=size(luminosityMap)
       if (luminosityRedshiftText(i) == "all") then
          ! Resize the arrays.
          call Move_Alloc (luminosityMap           ,luminosityMapTmp           )
          call Move_Alloc (luminosityRedshiftText  ,luminosityRedshiftTextTmp  )
          call Move_Alloc (luminosityRedshift      ,luminosityRedshiftTmp      )
          call Move_Alloc (luminosityBandRedshift  ,luminosityBandRedshiftTmp  )
          call Move_Alloc (luminosityFilter        ,luminosityFilterTmp        )
          call Move_Alloc (luminosityType          ,luminosityTypeTmp          )
          call Move_Alloc (luminosityPostprocessSet,luminosityPostprocessSetTmp)
          allocate(luminosityRedshiftText  (luminosityCount+outputCount-1))
          allocate(luminosityFilter        (luminosityCount+outputCount-1))
          allocate(luminosityType          (luminosityCount+outputCount-1))
          allocate(luminosityPostprocessSet(luminosityCount+outputCount-1))
          allocate(luminosityMap           (luminosityCount+outputCount-1))
          allocate(luminosityRedshift      (luminosityCount+outputCount-1))
          allocate(luminosityBandRedshift  (luminosityCount+outputCount-1))
          call Stellar_Luminosities_Expand_Filter_Set( &
               & i                          ,          &
               & outputCount                ,          &
               & luminosityMap              ,          &
               & luminosityRedshiftText     ,          &
               & luminosityFilter           ,          &
               & luminosityType             ,          &
               & luminosityPostprocessSet   ,          &
               & luminosityRedshift         ,          &
               & luminosityBandRedshift     ,          &
               & luminosityMapTmp           ,          &
               & luminosityRedshiftTextTmp  ,          &
               & luminosityFilterTmp        ,          &
               & luminosityTypeTmp          ,          &
               & luminosityPostprocessSetTmp,          &
               & luminosityRedshiftTmp      ,          &
               & luminosityBandRedshiftTmp             &
               &                                     )
          ! Modify new filters.
          do j=1,outputCount
             write (redshiftLabel,*) outputRedshifts(j)
             luminosityRedshiftText   (j+i-1)=redshiftLabel
             luminosityRedshift       (j+i-1)=outputRedshifts(j)
             if (luminosityBandRedshiftTmp  (i) <= -2.0d0) then
                luminosityBandRedshift(j+i-1)=outputRedshifts(j)
             else
                luminosityBandRedshift(j+i-1)=luminosityBandRedshiftTmp  (i)
             end if
          end do
          deallocate(luminosityRedshiftTextTmp  )
          deallocate(luminosityFilterTmp        )
          deallocate(luminosityTypeTmp          )
          deallocate(luminosityPostprocessSetTmp)
          deallocate(luminosityMapTmp           )
          deallocate(luminosityRedshiftTmp      )
          deallocate(luminosityBandRedshiftTmp  )
          ! Update count of luminosities.
          luminosityCount=size(luminosityMap)
       end if
       ! Alias for filters required for emission line calculations.
       if (luminosityFilter(i) == "emissionLineFilters") then
          call Move_Alloc (luminosityMap           ,luminosityMapTmp           )
          call Move_Alloc (luminosityRedshiftText  ,luminosityRedshiftTextTmp  )
          call Move_Alloc (luminosityRedshift      ,luminosityRedshiftTmp      )
          call Move_Alloc (luminosityBandRedshift  ,luminosityBandRedshiftTmp  )
          call Move_Alloc (luminosityFilter        ,luminosityFilterTmp        )
          call Move_Alloc (luminosityType          ,luminosityTypeTmp          )
          call Move_Alloc (luminosityPostprocessSet,luminosityPostprocessSetTmp)
          allocate(luminosityRedshiftText  (luminosityCount+3_c_size_t-1))
          allocate(luminosityFilter        (luminosityCount+3_c_size_t-1))
          allocate(luminosityType          (luminosityCount+3_c_size_t-1))
          allocate(luminosityPostprocessSet(luminosityCount+3_c_size_t-1))
          allocate(luminosityMap           (luminosityCount+3_c_size_t-1))
          allocate(luminosityRedshift      (luminosityCount+3_c_size_t-1))
          allocate(luminosityBandRedshift  (luminosityCount+3_c_size_t-1))
          call Stellar_Luminosities_Expand_Filter_Set( &
               & i                          ,          &
               & 3_c_size_t                 ,          &
               & luminosityMap              ,          &
               & luminosityRedshiftText     ,          &
               & luminosityFilter           ,          &
               & luminosityType             ,          &
               & luminosityPostprocessSet   ,          &
               & luminosityRedshift         ,          &
               & luminosityBandRedshift     ,          &
               & luminosityMapTmp           ,          &
               & luminosityRedshiftTextTmp  ,          &
               & luminosityFilterTmp        ,          &
               & luminosityTypeTmp          ,          &
               & luminosityPostprocessSetTmp,          &
               & luminosityRedshiftTmp      ,          &
               & luminosityBandRedshiftTmp             &
               &                                     )
             ! Create new filters.
             luminosityRedshiftText   (i:i+2)=luminosityRedshiftTextTmp  (i)
             luminosityRedshift       (i:i+2)=luminosityRedshiftTmp      (i)
             luminosityBandRedshift   (i:i+2)=luminosityBandRedshiftTmp  (i)
             luminosityPostprocessSet (i:i+2)=luminosityPostprocessSetTmp(i)
             luminosityFilter         (  i+0)=var_str("Lyc"            )
             luminosityFilter         (  i+1)=var_str("HeliumContinuum")
             luminosityFilter         (  i+2)=var_str("OxygenContinuum")
             luminosityType           (  i+0)=var_str("rest"           )
             luminosityType           (  i+1)=var_str("rest"           )
             luminosityType           (  i+2)=var_str("rest"           )
          deallocate(luminosityRedshiftTextTmp  )
          deallocate(luminosityFilterTmp        )
          deallocate(luminosityTypeTmp          )
          deallocate(luminosityPostprocessSetTmp)
          deallocate(luminosityMapTmp           )
          deallocate(luminosityRedshiftTmp      )
          deallocate(luminosityBandRedshiftTmp  )
          ! Update count of luminosities.
          luminosityCount=size(luminosityMap)
       end if
       ! Arrays of top-hat filters.
       if (extract(luminosityFilter(i),1,27) == "fixedResolutionTopHatArray_") then
          call String_Split_Words(specialFilterWords,char(luminosityFilter(i)),separator="_")
          word=char(specialFilterWords(2))
          read (word,*) wavelengthMinimum
          word=char(specialFilterWords(3))
          read (word,*) wavelengthMaximum
          word=char(specialFilterWords(4))
          read (word,*) resolution
          ! Determine the ratio of central wavelengths for successive filters.
          wavelengthRatio=                                &
               &  (sqrt(4.0d0*resolution**2+1.0d0)+1.0d0) &
               & /(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)
          newFilterCount=0
          wavelengthCentral=wavelengthMinimum/((sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution)
          do while (wavelengthCentral < wavelengthMaximum)
             newFilterCount=newFilterCount+1
             wavelengthCentral=wavelengthCentral*wavelengthRatio
          end do
          ! Resize the arrays.
          call Move_Alloc (luminosityMap           ,luminosityMapTmp           )
          call Move_Alloc (luminosityRedshiftText  ,luminosityRedshiftTextTmp  )
          call Move_Alloc (luminosityRedshift      ,luminosityRedshiftTmp      )
          call Move_Alloc (luminosityBandRedshift  ,luminosityBandRedshiftTmp  )
          call Move_Alloc (luminosityFilter        ,luminosityFilterTmp        )
          call Move_Alloc (luminosityType          ,luminosityTypeTmp          )
          call Move_Alloc (luminosityPostprocessSet,luminosityPostprocessSetTmp)
          allocate(luminosityRedshiftText  (luminosityCount+newFilterCount-1))
          allocate(luminosityFilter        (luminosityCount+newFilterCount-1))
          allocate(luminosityType          (luminosityCount+newFilterCount-1))
          allocate(luminosityPostprocessSet(luminosityCount+newFilterCount-1))
          allocate(luminosityMap           (luminosityCount+newFilterCount-1))
          allocate(luminosityRedshift      (luminosityCount+newFilterCount-1))
          allocate(luminosityBandRedshift  (luminosityCount+newFilterCount-1))
          call Stellar_Luminosities_Expand_Filter_Set( &
               & i                          ,          &
               & newFilterCount             ,          &
               & luminosityMap              ,          &
               & luminosityRedshiftText     ,          &
               & luminosityFilter           ,          &
               & luminosityType             ,          &
               & luminosityPostprocessSet   ,          &
               & luminosityRedshift         ,          &
               & luminosityBandRedshift     ,          &
               & luminosityMapTmp           ,          &
               & luminosityRedshiftTextTmp  ,          &
               & luminosityFilterTmp        ,          &
               & luminosityTypeTmp          ,          &
               & luminosityPostprocessSetTmp,          &
               & luminosityRedshiftTmp      ,          &
               & luminosityBandRedshiftTmp             &
               &                                     )
          ! Compute central wavelength of the initial filter.
          j=0
          wavelengthCentral=wavelengthMinimum/((sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution)
          do while (wavelengthCentral < wavelengthMaximum)
             j=j+1
             ! Compute the appropriate filter name.
             write (wavelengthCentralLabel,'(f11.3)') wavelengthCentral
             write (       resolutionLabel,'(f10.2)') resolution
             write (newFilterName,'(a,a,a,a)') "fixedResolutionTopHat_",trim(adjustl(wavelengthCentralLabel)),"_",trim(adjustl(resolutionLabel))
             ! Create new filter.
             luminosityRedshiftText   (j+i-1)=luminosityRedshiftTextTmp  (i)
             luminosityRedshift       (j+i-1)=luminosityRedshiftTmp      (i)
             luminosityBandRedshift   (j+i-1)=luminosityBandRedshiftTmp  (i)
             luminosityFilter         (j+i-1)=trim(newFilterName)
             luminosityType           (j+i-1)=luminosityTypeTmp          (i)
             luminosityPostprocessSet (j+i-1)=luminosityPostprocessSetTmp(i)
             ! Increase the central wavelength.
             wavelengthCentral=wavelengthCentral*wavelengthRatio
          end do
          deallocate(luminosityRedshiftTextTmp  )
          deallocate(luminosityFilterTmp        )
          deallocate(luminosityTypeTmp          )
          deallocate(luminosityPostprocessSetTmp)
          deallocate(luminosityMapTmp           )
          deallocate(luminosityRedshiftTmp      )
          deallocate(luminosityBandRedshiftTmp  )
          ! Update count of luminosities.
          luminosityCount=size(luminosityMap)
       end if
       ! Arrays of top-hat filters for SEDs
       if (extract(luminosityFilter(i),1,30) == "adaptiveResolutionTopHatArray_") then
          call String_Split_Words(specialFilterWords,char(luminosityFilter(i)),separator="_")
          word=char(specialFilterWords(3))
          read (word,*) wavelengthMinimum
          word=char(specialFilterWords(4))
          read (word,*) wavelengthMaximum
          word=char(specialFilterWords(5))
          read (word,*) observedWidth
          ! Set rest wavelength limits and rest wavelength width
          restWavelengthMinimum = wavelengthMinimum/(1.0d0+luminosityRedshift(i))
          restWavelengthMaximum = wavelengthMaximum/(1.0d0+luminosityRedshift(i))
          restWidth = observedWidth/(1.0d0+luminosityRedshift(i))
          ! Count number of filters that need to be added.
          newFilterCount   =0
          wavelengthCentral=restWavelengthMinimum
          tabulatedWidth   =stellarPopulationSpectra_%wavelengthInterval(wavelengthCentral)
          filterWidth      =restWidth
          if(tabulatedWidth > restWidth) filterWidth=tabulatedWidth
          do while (wavelengthCentral < wavelengthMaximum)
             if(wavelengthCentral < wavelengthMaximum) newFilterCount=newFilterCount+1
             call Stellar_Luminosities_SED_Top_Hat_Step(                              &
                  &                                     wavelengthCentral           , &
                  &                                     filterWidth                 , &
                  &                                     wavelengthMinimum           , &
                  &                                     wavelengthMaximum           , &
                  &                                     observedWidth               , &
                  &                                     luminosityRedshift       (i), &
                  &                                     stellarPopulationSpectra_     &
                  &                                    )
          end do
          ! Resize the arrays.
          call Move_Alloc (luminosityMap           ,luminosityMapTmp           )
          call Move_Alloc (luminosityRedshiftText  ,luminosityRedshiftTextTmp  )
          call Move_Alloc (luminosityRedshift      ,luminosityRedshiftTmp      )
          call Move_Alloc (luminosityBandRedshift  ,luminosityBandRedshiftTmp  )
          call Move_Alloc (luminosityFilter        ,luminosityFilterTmp        )
          call Move_Alloc (luminosityType          ,luminosityTypeTmp          )
          call Move_Alloc (luminosityPostprocessSet,luminosityPostprocessSetTmp)
          allocate(luminosityRedshiftText  (luminosityCount+newFilterCount-1))
          allocate(luminosityFilter        (luminosityCount+newFilterCount-1))
          allocate(luminosityType          (luminosityCount+newFilterCount-1))
          allocate(luminosityPostprocessSet(luminosityCount+newFilterCount-1))
          allocate(luminosityMap           (luminosityCount+newFilterCount-1))
          allocate(luminosityRedshift      (luminosityCount+newFilterCount-1))
          allocate(luminosityBandRedshift  (luminosityCount+newFilterCount-1))
          call Stellar_Luminosities_Expand_Filter_Set( &
               & i                          ,          &
               & newFilterCount             ,          &
               & luminosityMap              ,          &
               & luminosityRedshiftText     ,          &
               & luminosityFilter           ,          &
               & luminosityType             ,          &
               & luminosityPostprocessSet   ,          &
               & luminosityRedshift         ,          &
               & luminosityBandRedshift     ,          &
               & luminosityMapTmp           ,          &
               & luminosityRedshiftTextTmp  ,          &
               & luminosityFilterTmp        ,          &
               & luminosityTypeTmp          ,          &
               & luminosityPostprocessSetTmp,          &
               & luminosityRedshiftTmp      ,          &
               & luminosityBandRedshiftTmp             &
               &                                     )
          ! Compute central wavelength of the initial filter.
          j                =0
          wavelengthCentral=restWavelengthMinimum
          filterWidth      =restWidth
          tabulatedWidth   =stellarPopulationSpectra_%wavelengthInterval(wavelengthCentral)
          if(tabulatedWidth > restWidth) filterWidth = tabulatedWidth
          do while (wavelengthCentral < wavelengthMaximum)
             if(wavelengthCentral < wavelengthMaximum) then
                j=j+1
                ! Compute the appropriate filter name.
                write (wavelengthCentralLabel,'(f11.3)') wavelengthCentral
                write (       resolutionLabel,'(f10.2)') filterWidth
                write (newFilterName,'(a,a,a,a)') "adaptiveResolutionTopHat_",trim(adjustl(wavelengthCentralLabel)),"_",trim(adjustl(resolutionLabel))
                ! Create new filter.
                luminosityRedshiftText   (j+i-1)=luminosityRedshiftTextTmp  (i)
                luminosityRedshift       (j+i-1)=luminosityRedshiftTmp      (i)
                luminosityBandRedshift   (j+i-1)=luminosityBandRedshiftTmp  (i)
                luminosityFilter         (j+i-1)=trim(newFilterName)
                luminosityType           (j+i-1)=luminosityTypeTmp          (i)
                luminosityPostprocessSet (j+i-1)=luminosityPostprocessSetTmp(i)
             end if
             ! Compute central wavelength and width of next top hat filter
             call Stellar_Luminosities_SED_Top_Hat_Step(                              &
                  &                                     wavelengthCentral           , &
                  &                                     filterWidth                 , &
                  &                                     wavelengthMinimum           , &
                  &                                     wavelengthMaximum           , &
                  &                                     observedWidth               , &
                  &                                     luminosityRedshift       (i), &
                  &                                     stellarPopulationSpectra_     &
                  &                                    )
          end do
          deallocate(luminosityRedshiftTextTmp  )
          deallocate(luminosityFilterTmp        )
          deallocate(luminosityTypeTmp          )
          deallocate(luminosityPostprocessSetTmp)
          deallocate(luminosityMapTmp           )
          deallocate(luminosityRedshiftTmp      )
          deallocate(luminosityBandRedshiftTmp  )
          ! Update count of luminosities.
          luminosityCount=size(luminosityMap)
       end if
       ! Arrays of top-hat filters for equivalent width calculations
       if (extract(luminosityFilter(i),1,26) == "emissionLineContinuumPair_") then
          call String_Split_Words(specialFilterWords,char(luminosityFilter(i)),separator="_")
          lineName=char(specialFilterWords(2))
          ! Determine emission line wavelength
          wavelengthCentral=emissionLineWavelength(lineName)
          ! Read resolution
          word=char(specialFilterWords(3))
          read (word,*) resolution
          ! Determine minimum and maximum wavelengths to draw filters between
          wavelengthRatio=+(sqrt(4.0d0*resolution**2+1.0d0)+1.0d0) &
               &          /(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)
          wavelengthMinimum = wavelengthCentral*(sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution
          wavelengthMinimum = wavelengthMinimum/wavelengthRatio
          wavelengthMaximum = wavelengthCentral*(sqrt(4.0d0*resolution**2+1.0d0)+1.0d0)/2.0d0/resolution
          wavelengthMaximum = wavelengthMaximum*wavelengthRatio
          ! Determine the ratio of central wavelengths for successive filters.
          newFilterCount=0
          wavelengthCentral=wavelengthMinimum/((sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution)
          wavelengthCentral = wavelengthCentral/wavelengthRatio
          do k=1,3,1
             if (k == 1 .or. k == 3) newFilterCount=newFilterCount+1
             wavelengthCentral=wavelengthCentral*wavelengthRatio
          end do
          ! Resize the arrays.
          call Move_Alloc (luminosityMap           ,luminosityMapTmp           )
          call Move_Alloc (luminosityRedshiftText  ,luminosityRedshiftTextTmp  )
          call Move_Alloc (luminosityRedshift      ,luminosityRedshiftTmp      )
          call Move_Alloc (luminosityBandRedshift  ,luminosityBandRedshiftTmp  )
          call Move_Alloc (luminosityFilter        ,luminosityFilterTmp        )
          call Move_Alloc (luminosityType          ,luminosityTypeTmp          )
          call Move_Alloc (luminosityPostprocessSet,luminosityPostprocessSetTmp)
          allocate(luminosityRedshiftText  (luminosityCount+newFilterCount-1))
          allocate(luminosityFilter        (luminosityCount+newFilterCount-1))
          allocate(luminosityType          (luminosityCount+newFilterCount-1))
          allocate(luminosityPostprocessSet(luminosityCount+newFilterCount-1))
          allocate(luminosityMap           (luminosityCount+newFilterCount-1))
          allocate(luminosityRedshift      (luminosityCount+newFilterCount-1))
          allocate(luminosityBandRedshift  (luminosityCount+newFilterCount-1))
          call Stellar_Luminosities_Expand_Filter_Set( &
               & i                          ,          &
               & newFilterCount             ,          &
               & luminosityMap              ,          &
               & luminosityRedshiftText     ,          &
               & luminosityFilter           ,          &
               & luminosityType             ,          &
               & luminosityPostprocessSet   ,          &
               & luminosityRedshift         ,          &
               & luminosityBandRedshift     ,          &
               & luminosityMapTmp           ,          &
               & luminosityRedshiftTextTmp  ,          &
               & luminosityFilterTmp        ,          &
               & luminosityTypeTmp          ,          &
               & luminosityPostprocessSetTmp,          &
               & luminosityRedshiftTmp      ,          &
               & luminosityBandRedshiftTmp             &
               &                                     )
          ! Compute central wavelength of the initial filter.
          j=0
          wavelengthCentral=wavelengthMinimum/((sqrt(4.0d0*resolution**2+1.0d0)-1.0d0)/2.0d0/resolution)
          !wavelengthCentral = wavelengthCentral/wavelengthRatio
          do k = 1,3,1
             ! Compute the appropriate filter name.
             write (wavelengthCentralLabel,'(f11.3)') wavelengthCentral
             write (       resolutionLabel,'(f10.2)') resolution
             write (newFilterName,'(a,a,a,a,a,a)') "emissionLineContinuumBracketed_",trim(adjustl(lineName)),&
                  "_",trim(adjustl(wavelengthCentralLabel)),"_",trim(adjustl(resolutionLabel))
             ! Create new filter.
             if (k == 1 .or. k == 3) then
                j=j+1
                luminosityRedshiftText   (j+i-1)=luminosityRedshiftTextTmp  (i)
                luminosityRedshift       (j+i-1)=luminosityRedshiftTmp      (i)
                luminosityBandRedshift   (j+i-1)=luminosityBandRedshiftTmp  (i)
                luminosityFilter         (j+i-1)=trim(newFilterName)
                luminosityType           (j+i-1)=luminosityTypeTmp          (i)
                luminosityPostprocessSet (j+i-1)=luminosityPostprocessSetTmp(i)
             end if
             ! Increase the central wavelength.
             wavelengthCentral=wavelengthCentral*wavelengthRatio
          end do
          deallocate(luminosityRedshiftTextTmp  )
          deallocate(luminosityFilterTmp        )
          deallocate(luminosityTypeTmp          )
          deallocate(luminosityPostprocessSetTmp)
          deallocate(luminosityMapTmp           )
          deallocate(luminosityRedshiftTmp      )
          deallocate(luminosityBandRedshiftTmp  )
          ! Update count of luminosities.
          luminosityCount=size(luminosityMap)
       end if
       ! Next luminosity.
       i=i+1
    end do
    !![
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="stellarPopulationSpectra_"/>
    !!]
    return
  end subroutine Stellar_Luminosities_Special_Cases

  subroutine Stellar_Luminosities_Expand_Filter_Set( &
       & expandFrom                 ,                &
       & expandCount                ,                &
       & luminosityMap              ,                &
       & luminosityRedshiftText     ,                &
       & luminosityFilter           ,                &
       & luminosityType             ,                &
       & luminosityPostprocessSet   ,                &
       & luminosityRedshift         ,                &
       & luminosityBandRedshift     ,                &
       & luminosityMapTmp           ,                &
       & luminosityRedshiftTextTmp  ,                &
       & luminosityFilterTmp        ,                &
       & luminosityTypeTmp          ,                &
       & luminosityPostprocessSetTmp,                &
       & luminosityRedshiftTmp      ,                &
       & luminosityBandRedshiftTmp                   &
       &                                           )
    !!{
    Expand the filter set by removing the filter at index {\normalfont \ttfamily expandFrom} by adding {\normalfont \ttfamily expandCount} replicas of the filter at that point.
    !!}
    use, intrinsic :: ISO_C_Binding    , only : c_size_t
    implicit none
    integer         (c_size_t      ), intent(in   )               :: expandFrom               , expandCount
    integer                         , intent(inout), dimension(:) :: luminosityMap
    type            (varying_string), intent(inout), dimension(:) :: luminosityRedshiftText   , luminosityFilter           , &
         &                                                           luminosityType           , luminosityPostprocessSet
    double precision                , intent(inout), dimension(:) :: luminosityRedshift       , luminosityBandRedshift
    integer                         , intent(inout), dimension(:) :: luminosityMapTmp
    type            (varying_string), intent(inout), dimension(:) :: luminosityRedshiftTextTmp, luminosityFilterTmp        , &
         &                                                           luminosityTypeTmp        , luminosityPostprocessSetTmp
    double precision                , intent(inout), dimension(:) :: luminosityRedshiftTmp    , luminosityBandRedshiftTmp
    integer                                                       :: luminosityCount

    luminosityCount=size(luminosityMapTmp)
    if (expandFrom > 1              ) then
       luminosityMap           (1            :expandFrom                          -1)=luminosityMapTmp           (1  :expandFrom-1            )
       luminosityRedshiftText  (1            :expandFrom                          -1)=luminosityRedshiftTextTmp  (1  :expandFrom-1            )
       luminosityRedshift      (1            :expandFrom                          -1)=luminosityRedshiftTmp      (1  :expandFrom-1            )
       luminosityBandRedshift  (1            :expandFrom                          -1)=luminosityBandRedshiftTmp  (1  :expandFrom-1            )
       luminosityFilter        (1            :expandFrom                          -1)=luminosityFilterTmp        (1  :expandFrom-1            )
       luminosityType          (1            :expandFrom                          -1)=luminosityTypeTmp          (1  :expandFrom-1            )
       luminosityPostprocessSet(1            :expandFrom                          -1)=luminosityPostprocessSetTmp(1  :expandFrom-1            )
    end if
    if (expandFrom < luminosityCount) then
       luminosityMap           (expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityMapTmp           (expandFrom+1:luminosityCount)
       luminosityRedshiftText  (expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityRedshiftTextTmp  (expandFrom+1:luminosityCount)
       luminosityRedshift      (expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityRedshiftTmp      (expandFrom+1:luminosityCount)
       luminosityBandRedshift  (expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityBandRedshiftTmp  (expandFrom+1:luminosityCount)
       luminosityFilter        (expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityFilterTmp        (expandFrom+1:luminosityCount)
       luminosityType          (expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityTypeTmp          (expandFrom+1:luminosityCount)
       luminosityPostprocessSet(expandFrom+expandCount:luminosityCount+expandCount-1)=luminosityPostprocessSetTmp(expandFrom+1:luminosityCount)
    end if
    luminosityMap              (expandFrom            :expandFrom     +expandCount-1)=luminosityMapTmp           (expandFrom                  )
    luminosityRedshiftText     (expandFrom            :expandFrom     +expandCount-1)=luminosityRedshiftTextTmp  (expandFrom                  )
    luminosityRedshift         (expandFrom            :expandFrom     +expandCount-1)=luminosityRedshiftTmp      (expandFrom                  )
    luminosityBandRedshift     (expandFrom            :expandFrom     +expandCount-1)=luminosityBandRedshiftTmp  (expandFrom                  )
    luminosityFilter           (expandFrom            :expandFrom     +expandCount-1)=luminosityFilterTmp        (expandFrom                  )
    luminosityType             (expandFrom            :expandFrom     +expandCount-1)=luminosityTypeTmp          (expandFrom                  )
    luminosityPostprocessSet   (expandFrom            :expandFrom     +expandCount-1)=luminosityPostprocessSetTmp(expandFrom                  )
    return
  end subroutine Stellar_Luminosities_Expand_Filter_Set

  subroutine Stellar_Luminosities_Truncate(self,templateLuminosities)
    !!{
    Truncate (or pad) the stellar luminosities to match the number in the given {\normalfont \ttfamily templateLuminosities}.
    !!}
    implicit none
    class           (stellarLuminosities), intent(inout)               :: self
    type            (stellarLuminosities), intent(in   )               :: templateLuminosities
    double precision                     , allocatable  , dimension(:) :: luminositiesTmp
    integer                                                            :: templateCount       , selfCount, &
         &                                                                minCount

    if (allocated(templateLuminosities%luminosityValue)) then
       templateCount=size(templateLuminosities%luminosityValue)
       if (allocated(self%luminosityValue)) then
          ! Our luminosities are allocated. Check size.
          selfCount=size(self%luminosityValue)
          if (selfCount /= templateCount) then
             ! Size does not match template. Reallocate and pad as necessary.
             minCount=min(templateCount,selfCount)
             call Move_Alloc(self%luminosityValue,luminositiesTmp)
             allocate(self%luminosityValue(templateCount))
             self%luminosityValue(1:minCount)=luminositiesTmp(1:minCount)
             if (templateCount > selfCount) self%luminosityValue(selfCount+1:templateCount)=0.0d0
             deallocate(luminositiesTmp)
          end if
       else
          ! Our luminosities are not allocated. Allocate and set to zero.
          allocate(self%luminosityValue(templateCount))
          self%luminosityValue=0.0d0
       end if
    else if (allocated(self%luminosityValue)) then
       ! Template luminosities are not allocated, simply deallocate our luminosities to match.
       deallocate(self%luminosityValue)
    end if
    return
  end subroutine Stellar_Luminosities_Truncate

  subroutine Stellar_Luminosities_Parameter_Map_Double(parameters)
    !!{
    Map an array of luminosity-related input parameters into a new array accounting for special case processing.
    !!}
    implicit none
    double precision, intent(inout), allocatable, dimension(:) :: parameters
    double precision               , allocatable, dimension(:) :: parametersMapped
    integer                                                    :: i

    ! Allocate new array.
    allocate(parametersMapped(luminosityCount))
    ! Map from the old array.
    do i=1,luminosityCount
       parametersMapped(i)=parameters(luminosityMap(i))
    end do
    ! Copy the new array.
    deallocate(parameters)
    call Move_Alloc(parametersMapped,parameters)
    return
  end subroutine Stellar_Luminosities_Parameter_Map_Double

  !![
  <stateStoreTask>
   <unitName>Stellar_Luminosities_State_Store</unitName>
  </stateStoreTask>
  !!]
  subroutine Stellar_Luminosities_State_Store(stateFile,gslStateFile,stateOperationID)
    !!{
    Write the luminosities state to file.
    !!}
    use            :: Display           , only : displayIndent, displayMessage, displayUnindent, verbosityLevelWorking
    use, intrinsic :: ISO_C_Binding     , only : c_ptr        , c_size_t
    use            :: ISO_Varying_String, only : operator(//) , var_str
    use            :: String_Handling   , only : operator(//)
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer                          :: i
    !$GLC attributes unused :: gslStateFile, stateOperationID

    call displayIndent  (var_str('storing state for "stellar luminosities" [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    call displayMessage(var_str('storing "luminosityCount" [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    write (stateFile) luminosityCount
    call displayMessage(var_str('storing luminosities [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    write (stateFile) luminosityIndex,luminosityCosmicTime,luminosityRedshift,luminosityBandRedshift
    do i=1,luminosityCount
       call displayMessage(var_str('storing luminosity ')//i//' of '//luminosityCount//' [position: '//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
       call luminosityName          (i)%stateStore(stateFile)
       call luminosityType          (i)%stateStore(stateFile)
       call luminosityFilter        (i)%stateStore(stateFile)
       call luminosityPostprocessSet(i)%stateStore(stateFile)
    end do
    call displayUnindent(var_str('done [position: ')//FTell(stateFile)//']'                                    ,verbosity=verbosityLevelWorking)
   return
  end subroutine Stellar_Luminosities_State_Store

  !![
  <stateRetrieveTask>
   <unitName>Stellar_Luminosities_State_Restore</unitName>
  </stateRetrieveTask>
  !!]
  subroutine Stellar_Luminosities_State_Restore(stateFile,gslStateFile,stateOperationID)
    !!{
    Retrieve the luminosities state from the file.
    !!}
    use            :: Display            , only : displayIndent   , displayMessage, displayUnindent, verbosityLevelWorking
    use, intrinsic :: ISO_C_Binding      , only : c_ptr           , c_size_t
    use            :: ISO_Varying_String , only : operator(//)    , var_str
    use            :: Instruments_Filters, only : Filter_Get_Index
    use            :: String_Handling    , only : operator(//)
    implicit none
    integer          , intent(in   ) :: stateFile
    integer(c_size_t), intent(in   ) :: stateOperationID
    type   (c_ptr   ), intent(in   ) :: gslStateFile
    integer                          :: i
    !$GLC attributes unused :: gslStateFile, stateOperationID

    call displayIndent  (var_str('restoring state for "stellar luminosities" [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    if (allocated(luminosityFilterIndex                  )) deallocate(luminosityFilterIndex                  )
    if (allocated(luminosityIndex                        )) deallocate(luminosityIndex                        )
    if (allocated(luminosityPostprocessor                )) deallocate(luminosityPostprocessor                )
    if (allocated(luminosityCosmicTime                   )) deallocate(luminosityCosmicTime                   )
    if (allocated(luminosityRedshift                     )) deallocate(luminosityRedshift                     )
    if (allocated(luminosityBandRedshift                 )) deallocate(luminosityBandRedshift                 )
    if (allocated(luminosityName                         )) deallocate(luminosityName                         )
    if (allocated(luminosityType                         )) deallocate(luminosityType                         )
    if (allocated(luminosityFilter                       )) deallocate(luminosityFilter                       )
    if (allocated(luminosityPostprocessSet               )) deallocate(luminosityPostprocessSet               )
    if (allocated(unitStellarLuminosities%luminosityValue)) deallocate(unitStellarLuminosities%luminosityValue)
    if (allocated(zeroStellarLuminosities%luminosityValue)) deallocate(zeroStellarLuminosities%luminosityValue)
    call displayMessage(var_str('restoring "luminosityCount" [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    read (stateFile) luminosityCount
    allocate(luminosityFilterIndex                  (luminosityCount))
    allocate(luminosityIndex                        (luminosityCount))
    allocate(luminosityPostprocessor                (luminosityCount))
    allocate(luminosityCosmicTime                   (luminosityCount))
    allocate(luminosityRedshift                     (luminosityCount))
    allocate(luminosityBandRedshift                 (luminosityCount))
    allocate(luminosityName                         (luminosityCount))
    allocate(luminosityType                         (luminosityCount))
    allocate(luminosityFilter                       (luminosityCount))
    allocate(luminosityPostprocessSet               (luminosityCount))
    allocate(unitStellarLuminosities%luminosityValue(luminosityCount))
    allocate(zeroStellarLuminosities%luminosityValue(luminosityCount))
    unitStellarLuminosities%luminosityValue=1.0d0
    zeroStellarLuminosities%luminosityValue=0.0d0
    call displayMessage(var_str('restoring luminosities [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    read (stateFile) luminosityIndex,luminosityCosmicTime,luminosityRedshift,luminosityBandRedshift
    do i=1,luminosityCount
       call displayMessage(var_str('restoring luminosity ')//i//' of '//luminosityCount//' [position: '//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
       call luminosityName          (i)%stateRestore(stateFile)
       call luminosityType          (i)%stateRestore(stateFile)
       call luminosityFilter        (i)%stateRestore(stateFile)
       call luminosityPostprocessSet(i)%stateRestore(stateFile)
       luminosityFilterIndex  (i)                                        =  Filter_Get_Index                                    (luminosityFilter        (i))
       luminosityPostprocessor(i)%stellarPopulationSpectraPostprocessor_ => stellarPopulationSpectraPostprocessorBuilder__%build(luminosityPostprocessSet(i))
    end do
    call displayUnindent(var_str('done [position: ')//FTell(stateFile)//']',verbosity=verbosityLevelWorking)
    return
  end subroutine Stellar_Luminosities_State_Restore

  subroutine sortByIndexPostprocessor(array,index)
    !!{
    Given an {\normalfont \ttfamily array}, sort it in place using the supplied index.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    type   (stellarPopulationSpectraPostprocessorList), dimension(:          ), intent(inout) :: array
    integer(kind=c_size_t                            ), dimension(:          ), intent(in   ) :: index
    type   (stellarPopulationSpectraPostprocessorList), dimension(size(array))                :: arrayTmp
    integer(kind=c_size_t                            )                                        :: i

    do i=1,size(array)
       arrayTmp(i)=array(index(i))
    end do
    array=arrayTmp
    return
  end subroutine sortByIndexPostprocessor

  function Stellar_Luminosities_Non_Static_Size_Of(self)
    !!{
    Return the size of any non-static components of the object.
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    integer(c_size_t           )                :: Stellar_Luminosities_Non_Static_Size_Of
    class  (stellarLuminosities), intent(in   ) :: self

    if (allocated(self%luminosityValue)) then
       Stellar_Luminosities_Non_Static_Size_Of=sizeof(self%luminosityValue)
    else
       Stellar_Luminosities_Non_Static_Size_Of=0_c_size_t
    end if
    return
  end function Stellar_Luminosities_Non_Static_Size_Of

end module Stellar_Luminosities_Structure
