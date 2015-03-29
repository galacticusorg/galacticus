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

!% Contains a module which defines the stellar luminosities object.

module Stellar_Luminosities_Structure
  !% Defines the stellar luminosities object.
  use ISO_Varying_String
  implicit none
  private
  public :: stellarLuminosities, max, operator(*), Stellar_Luminosities_Parameter_Map

  ! Interface to parameter mapping functions.
  interface Stellar_Luminosities_Parameter_Map
     module procedure Stellar_Luminosities_Parameter_Map_Double
  end interface Stellar_Luminosities_Parameter_Map
  
  ! Interface to max() function for stellar luminosities objects.
  interface max
     module procedure stellarLuminositiesMax
  end interface max

  ! Interface to multiplication operators with stellar luminosities objects as their second argument.
  interface operator(*)
     module procedure Stellar_Luminosities_Multiply_Switched
  end interface operator(*)

  type stellarLuminosities
     !% The stellar luminosities structure.
     private
     double precision, allocatable, dimension(:) :: luminosityValue
   contains
     !@ <objectMethods>
     !@   <object>stellarLuminosities</object>
     !@   <objectMethod>
     !@     <method>multiply</method>
     !@     <type>\textcolor{red}{\textless type(stellarLuminosities)\textgreater}</type>
     !@     <arguments>\doublezero\ multiplier\argin</arguments>
     !@     <description>Multiply stellar luminosities by a scalar.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>divide</method>
     !@     <type>\textcolor{red}{\textless type(stellarLuminosities)\textgreater}</type>
     !@     <arguments>\doublezero\ divisor\argin</arguments>
     !@     <description>Divide stellar luminosities by a scalar.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>add</method>
     !@     <type>\textcolor{red}{\textless type(stellarLuminosities)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(stellarLuminosities)\textgreater} stellarLuminosities2\argin</arguments>
     !@     <description>Add two stellarLuminosities.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>subtract</method>
     !@     <type>\textcolor{red}{\textless type(stellarLuminosities)\textgreater}</type>
     !@     <arguments>\textcolor{red}{\textless type(stellarLuminosities)\textgreater} stellarLuminosities2\argin</arguments>
     !@     <description>Subtract one abundance from another.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>increment</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(stellarLuminosities)\textgreater} addStellarLuminosities\argin</arguments>
     !@     <description>Increment a stellar luminosities object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serializeCount</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return a count of the number of properties in a serialized stellar luminosities object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>serialize</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ array\argout</arguments>
     !@     <description>Serialize a stellar luminosities object to an array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>deserialize</method>
     !@     <type>\void</type>
     !@     <arguments>\doubleone\ array\argin</arguments>
     !@     <description>Deserialize a stellar luminosities object from an array.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isZero</method>
     !@     <type>\logicalzero</type>
     !@     <arguments></arguments>
     !@     <description>Return true if a stellar luminosities object is zero.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>destroy</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Destroy a stellar luminosities object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>reset</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Reset a stellar luminosities object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>builder</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *type(node)\textgreater} stellarLuminositiesDefinition\argin</arguments>
     !@     <description>Build a stellar luminosities object from a provided XML description.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>dump</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Dump a stellar luminosities object.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <type>\void</type>
     !@     <method>dumpRaw</method>
     !@     <arguments>\intzero\ fileHandle\argin</arguments>
     !@     <description>Dump a stellar luminosities object to binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readRaw</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ fileHandle\argin</arguments>
     !@     <description>Read a stellar luminosities object from binary.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setToUnity</method>
     !@     <type>\void</type>
     !@     <arguments></arguments>
     !@     <description>Set a stellar luminosities object to unity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>luminosity</method>
     !@     <type>\doublezero</type>
     !@     <arguments>\intzero\ index\argin</arguments>
     !@     <description>Return the $i^{\mathrm th}$ luminosity.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>output</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ integerProperty\arginout, \intzero\ integerBufferCount\arginout, \inttwo\ integerBuffer\arginout, \intzero doubleProperty\arginout, \intzero\ doubleBufferCount\arginout, \doubletwo\ doubleBuffer\arginout, \doublezero\ time\argin, \intzero\ instance\argin</arguments>
     !@     <description>Store a stellar luminosities object in the output buffers.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>luminosityOutputCount</method>
     !@     <type>\intzero</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Return the number of luminosities to be output at the given time.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>outputCount</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ integerPropertyCount\arginout, \intzero\ doublePropertyCount\arginout, \doublezero\ time\argin, \intzero\ instance\argin</arguments>
     !@     <description>Specify the count of a stellar luminosities object for output.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>outputNames</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ integerProperty\arginout, \textcolor{red}{\textless char[*](:)\textgreater} integerPropertyNames\arginout, \textcolor{red}{\textless char[*](:)\textgreater} integerPropertyComments\arginout, \doubleone\ integerPropertyUnitsSI\arginout, \intzero\ doubleProperty\arginout, \textcolor{red}{\textless char[*](:)\textgreater} doublePropertyNames\arginout, \textcolor{red}{\textless char[*](:)\textgreater} doublePropertyComments\arginout, \doubleone\ doublePropertyUnitsSI\arginout, \doublezero\ time\argin, \intzero\ instance\argin</arguments>
     !@     <description>Specify the names of stellar luminosities object properties for output.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>luminosityCount</method>
     !@     <type>\intzero</type>
     !@     <arguments>\logicalzero\ [unmapped]</arguments>
     !@     <description>Return the total number of luminosities tracked. If {\normalfont \ttfamily unmapped} is true, then the number of luminosities prior to mapping is returned.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>setLuminosities</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ mass\argin,\intzero\ imfSelected\argin,\doublezero currentTime\argin,\textcolor{red}{\textless type(abundances)\textgreater} fuelAbundances\argin</arguments>
     !@     <description>Set the luminosities using a single stellar population.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>isOutput</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\intzero\ index\argin, \doublezero\ time\argin</arguments>
     !@     <description>Return true if the indexed luminosity is to be output at the given time.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>index</method>
     !@     <type>\intzero</type>
     !@     <arguments>\textcolor{red}{\textless type(varying\_string)\textgreater} name\argin</arguments>
     !@     <description>Return the index to a luminosity specified by name.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>name</method>
     !@     <type>\textcolor{red}{\textless type(varying\_string)\textgreater}</type>
     !@     <arguments>\intzero\ index\argin</arguments>
     !@     <description>Return the name of a luminosity specified by index.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>truncate</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(stellarLuminosities)\textgreater} templateLuminosities\argin</arguments>
     !@     <description>Truncate the number of stellar luminosities stored to match that in the given {\normalfont \ttfamily templateLuminosities}.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final             ::                          Stellar_Luminosities_Destructor
     procedure         :: add                   => Stellar_Luminosities_Add
     procedure         :: subtract              => Stellar_Luminosities_Subtract
     procedure         :: multiply              => Stellar_Luminosities_Multiply
     procedure         :: divide                => Stellar_Luminosities_Divide
     generic           :: operator(+)           => add
     generic           :: operator(-)           => subtract
     generic           :: operator(*)           => multiply
     generic           :: operator(/)           => divide
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
     procedure, nopass :: luminosityOutputCount => Stellar_Luminosities_Output_Count_Get
     procedure         :: outputCount           => Stellar_Luminosities_Output_Count
     procedure         :: outputNames           => Stellar_Luminosities_Output_Names
     procedure, nopass :: isOutput              => Stellar_Luminosities_Is_Output
     procedure, nopass :: index                 => Stellar_Luminosities_Index
     procedure, nopass :: name                  => Stellar_Luminosities_Name
     procedure         :: truncate              => Stellar_Luminosities_Truncate
  end type stellarLuminosities

  ! Flag specifying if module has been initialized.
  logical                                                          :: luminositiesInitialized           =.false.

  ! Arrays which hold the luminosity specifications.
  integer                                                          :: luminosityCount                           , luminosityCountUnmapped
  integer                              , allocatable, dimension(:) :: luminosityFilterIndex                     , luminosityIndex               , &
       &                                                              luminosityPostprocessingChainIndex        , luminosityMap
  double precision                     , allocatable, dimension(:) :: luminosityBandRedshift                    , luminosityCosmicTime          , &
       &                                                              luminosityRedshift
  type            (varying_string     ), allocatable, dimension(:) :: luminosityFilter                          , luminosityName                , &
       &                                                              luminosityPostprocessSet                  , luminosityType                , &
       &                                                              luminosityRedshiftText

  ! Luminosity output options.
  integer                                                          :: luminosityOutputOption
  integer                              , parameter                 :: luminosityOutputOptionAll         =0      , luminosityOutputOptionFuture=1, &
       &                                                              luminosityOutputOptionPresent     =2

  ! Unit and zero stellarLuminosities objects.
  type            (stellarLuminosities), public                    :: unitStellarLuminosities                   , zeroStellarLuminosities

contains

  subroutine Stellar_Luminosities_Initialize
    !% Initialize the {\normalfont \ttfamily stellarLuminositiesStructure} object module. Determines which stellar luminosities are to be tracked.
    use, intrinsic :: ISO_C_Binding
    use Input_Parameters
    use Galacticus_Error
    use Memory_Management
    use Instruments_Filters
    use Cosmology_Functions
    use Stellar_Population_Spectra_Postprocess
    use Array_Utilities
    use Sort
    implicit none
    class           (cosmologyFunctionsClass), pointer                   :: cosmologyFunctionsDefault
    integer                                                              :: iLuminosity               , jLuminosity
    double precision                                                     :: expansionFactor
    character       (len=10                 )                            :: redshiftLabel
    type            (varying_string         )                            :: luminosityOutputOptionText
    integer         (c_size_t               ), allocatable, dimension(:) :: luminosityTimeIndex
    
    ! Initialize the module if necessary.
    if (.not.luminositiesInitialized) then
       !$omp critical (Stellar_Luminosities_Initialize)
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
             call Galacticus_Error_Report("Stellar_Luminosities_Initialize","unrecognized luminosityOutputOption")
          end select

          ! Read in the parameters which specify the luminosities to be computed.
          luminosityCount        =Get_Input_Parameter_Array_Size('luminosityRedshift')
          luminosityCountUnmapped=luminosityCount
          if (Get_Input_Parameter_Array_Size('luminosityFilter') /= luminosityCount) call&
               & Galacticus_Error_Report('Stellar_Luminosities_Initialize','luminosityFilter and luminosityRedshift&
               & input arrays must have same dimension')
          if (Get_Input_Parameter_Array_Size('luminosityType') /= luminosityCount) call&
               & Galacticus_Error_Report('Stellar_Luminosities_Initialize','luminosityType and luminosityRedshift&
               & input arrays must have same dimension')
          if (Input_Parameter_Is_Present('luminosityBandRedshift')) then
             if (Get_Input_Parameter_Array_Size('luminosityBandRedshift') /= luminosityCount) call&
                  & Galacticus_Error_Report('Stellar_Luminosities_Initialize','luminosityBandRedshift and luminosityRedshift&
                  & input arrays must have same dimension')
          end if

          if (luminosityCount > 0) then
             call Alloc_Array(luminosityMap                     ,[luminosityCount])
             call Alloc_Array(luminosityRedshift                ,[luminosityCount])
             call Alloc_Array(luminosityBandRedshift            ,[luminosityCount])
             allocate(luminosityFilter        (luminosityCount))
             allocate(luminosityType          (luminosityCount))
             allocate(luminosityPostprocessSet(luminosityCount))
             allocate(luminosityRedshiftText  (luminosityCount))
             call Memory_Usage_Record(sizeof(luminosityFilter)+sizeof(luminosityType)+sizeof(luminosityPostprocessSet),blockCount=4)
             !@ <inputParameter>
             !@   <name>luminosityRedshift</name>
             !@   <attachedTo>module</attachedTo>
             !@   <description>
             !@     The redshift for which to compute each specified stellar luminosity.
             !@   </description>
             !@   <type>real</type>
             !@   <cardinality>0..*</cardinality>
             !@ </inputParameter>
             call Get_Input_Parameter('luminosityRedshift',luminosityRedshiftText)
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
             if (Input_Parameter_Is_Present('luminosityBandRedshift')) then
                !@ <inputParameter>
                !@   <name>luminosityBandRedshift</name>
                !@   <attachedTo>module</attachedTo>
                !@   <description>
                !@     If present, force filters to be shifted to this redshift rather than that specified by {\normalfont \ttfamily [luminosityRedshift]}. Allows sampling of the SED at wavelengths corresponding to other redshifts.
                !@   </description>
                !@   <type>real</type>
                !@   <cardinality>0..*</cardinality>
                !@ </inputParameter>
                call Get_Input_Parameter('luminosityBandRedshift',luminosityBandRedshift)
             else                
                luminosityBandRedshift=luminosityRedshift
             end if
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
             ! Read postprocessing set information.
             if (Get_Input_Parameter_Array_Size('luminosityPostprocessSet') > 0) then
                if (Get_Input_Parameter_Array_Size('luminosityPostprocessSet') /= luminosityCount) call&
                     & Galacticus_Error_Report('Stellar_Luminosities_Initialize','luminosityPostprocessSet and luminosityFilter&
                     & input arrays must have same dimension')
                !@ <inputParameter>
                !@   <name>luminosityPostprocessSet</name>
                !@   <attachedTo>module</attachedTo>
                !@   <description>
                !@     The name of the set of postprocessing algorithms to apply to this filter.
                !@   </description>
                !@   <type>string</type>
                !@   <cardinality>0..*</cardinality>
                !@ </inputParameter>
                call Get_Input_Parameter('luminosityPostprocessSet',luminosityPostprocessSet)
             else
                luminosityPostprocessSet="default"
             end if
             ! Handle luminosity definition special cases.
             call Stellar_Luminosities_Special_Cases(luminosityMap,luminosityRedshiftText,luminosityRedshift,luminosityBandRedshift,luminosityFilter,luminosityType,luminosityPostprocessSet)
             luminosityCount=size(luminosityRedshift)
             ! Allocate remaining required arrays.
             allocate(luminosityName(luminosityCount))
             call Alloc_Array(luminosityFilterIndex             ,[luminosityCount])
             call Alloc_Array(luminosityIndex                   ,[luminosityCount])
             call Alloc_Array(luminosityPostprocessingChainIndex,[luminosityCount])
             call Alloc_Array(luminosityCosmicTime              ,[luminosityCount])
             call Alloc_Array(luminosityTimeIndex               ,[luminosityCount])
             ! Get the default cosmology functions object.
             cosmologyFunctionsDefault => cosmologyFunctions()
             ! Process the list of luminosities.
             do iLuminosity=1,luminosityCount
                ! Assign a name to this luminosity.
                write (redshiftLabel,'(f7.4)') luminosityBandRedshift(iLuminosity)
                luminosityName(iLuminosity)=luminosityFilter(iLuminosity)//":"//luminosityType(iLuminosity)//":z"&
                     &//trim(adjustl(redshiftLabel))
                if (luminosityPostprocessSet(iLuminosity) /= "default") luminosityName(iLuminosity)=luminosityName(iLuminosity)&
                     &//":"//char(luminosityPostprocessSet(iLuminosity))
                ! Check for duplicated luminosities.
                if (iLuminosity > 1) then
                   do jLuminosity=1,iLuminosity-1
                      if (luminosityName(iLuminosity) == luminosityName(jLuminosity)) then
                         call Galacticus_Error_Report('Stellar_Luminosities_Initialize','luminosity '//luminosityName(iLuminosity)//' appears more than once in the input parameter file')
                      end if
                   end do
                end if
                ! Assign an index for this luminosity.
                luminosityIndex(iLuminosity)=iLuminosity
                ! Get the index of the specified filter.
                luminosityFilterIndex(iLuminosity)=Filter_Get_Index(luminosityFilter(iLuminosity))
                ! Set the reference time (i.e. cosmological time corresponding to the specified redshift) for this filter.
                expansionFactor=cosmologyFunctionsDefault%expansionFactorFromRedshift(luminosityRedshift(iLuminosity))
                luminosityCosmicTime(iLuminosity)=cosmologyFunctionsDefault%cosmicTime(expansionFactor)
                ! Set the filter redshifting factor. This is equal to the requested redshift if an observed frame was specified, otherwise
                ! it is set to zero to indicate a rest-frame filter.
                select case(char(luminosityType(iLuminosity)))
                case ("rest")
                   luminosityBandRedshift(iLuminosity)=0.0d0
                case ("observed")
                   ! Do nothing, we already have the correct redshift.
                case default
                   call Galacticus_Error_Report('Stellar_Luminosities_Initialize','unrecognized filter type')
                end select
                ! Find the index for the postprocessing chain to be applied to this filter.
                luminosityPostprocessingChainIndex(iLuminosity)=Stellar_Population_Spectrum_Postprocess_Index(luminosityPostprocessSet(iLuminosity))
             end do
             ! Sort the luminosities such that the latest luminosities are stored first.
             luminosityTimeIndex=Array_Reverse(Sort_Index_Do(luminosityCosmicTime))
             call Sort_By_Index(luminosityFilterIndex             ,luminosityTimeIndex)
             call Sort_By_Index(luminosityPostprocessingChainIndex,luminosityTimeIndex)
             call Sort_By_Index(luminosityCosmicTime              ,luminosityTimeIndex)
             call Sort_By_Index(luminosityName                    ,luminosityTimeIndex)
             call Sort_By_Index(luminosityRedshift                ,luminosityTimeIndex)
             call Sort_By_Index(luminosityBandRedshift            ,luminosityTimeIndex)             
             ! Allocate unit and zero stellar abundance objects.
             call Alloc_Array(unitStellarLuminosities%luminosityValue,[luminosityCount])
             call Alloc_Array(zeroStellarLuminosities%luminosityValue,[luminosityCount])
             unitStellarLuminosities%luminosityValue=1.0d0
             zeroStellarLuminosities%luminosityValue=0.0d0
          end if
          luminositiesInitialized=.true.
       end if
       !$omp end critical (Stellar_Luminosities_Initialize)
    end if
    return
  end subroutine Stellar_Luminosities_Initialize

  subroutine Stellar_Luminosities_Destructor(self)
    !% Destructor for a {\normalfont \ttfamily stellarLuminosities} object.
    use Memory_Management
    implicit none
    type(stellarLuminosities), intent(inout) :: self

    if (allocated(self%luminosityValue)) call Dealloc_Array(self%luminosityValue)
    return
  end subroutine Stellar_Luminosities_Destructor

  subroutine Stellar_Luminosities_Destroy(self)
    !% Destroy an stellarLuminosities object.
    use Memory_Management
    implicit none
    class(stellarLuminosities), intent(inout) :: self

    if (allocated(self%luminosityValue)) call Dealloc_Array(self%luminosityValue)
    return
  end subroutine Stellar_Luminosities_Destroy

  subroutine Stellar_Luminosities_Builder(self,stellarLuminositiesDefinition)
    !% Build a {\normalfont \ttfamily stellarLuminosities} object from the given XML {\normalfont \ttfamily stellarLuminositiesDefinition}.
    use FoX_DOM
    use Galacticus_Error
    implicit none
    class  (stellarLuminosities), intent(inout)          :: self
    type   (node               ), intent(in   ), pointer :: stellarLuminositiesDefinition
    type   (node               )               , pointer :: luminosity
    type   (nodeList           )               , pointer :: luminosityList
    integer                                              :: i

    ! Get the metallicity.
    luminosityList => getElementsByTagName(stellarLuminositiesDefinition,'luminosity')
    if (luminosityCount > 0) then
       do i=0,luminosityCount-1
          luminosity => item(luminosityList,i)
          call extractDataContent(luminosity,self%luminosityValue(i))
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Builder

  subroutine Stellar_Luminosities_Dump(self)
    !% Dump a stellar luminosities object.
    use Galacticus_Display
    implicit none
    class    (stellarLuminosities), intent(in   ) :: self
    integer                                       :: i
    character(len=12             )                :: label
    type     (varying_string     )                :: message

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Dump the contents.
    if (luminosityCount > 0) then
       do i=1,luminosityCount
          if (i <= size(self%luminosityValue)) then
             write (label,'(e12.6)') self%luminosityValue(i)
          else
             label="pruned"
          end if
          message=luminosityName(i)//':          '//label
          call Galacticus_Display_Message(message)
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Dump

  subroutine Stellar_Luminosities_Dump_Raw(self,fileHandle)
    !% Dump an stellarLuminosities object to binary.
    implicit none
    class  (stellarLuminosities), intent(in   ) :: self
    integer                     , intent(in   ) :: fileHandle

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Dump the content.
    if (luminosityCount > 0) then
       write (fileHandle) size(self%luminosityValue)
       write (fileHandle) self%luminosityValue
    end if
    return
  end subroutine Stellar_Luminosities_Dump_Raw

  subroutine Stellar_Luminosities_Read_Raw(self,fileHandle)
    !% Read an stellarLuminosities object from binary.
    use Memory_Management
    implicit none
    class  (stellarLuminosities), intent(inout) :: self
    integer                     , intent(in   ) :: fileHandle
    integer                                     :: luminosityActiveCount
    
    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Read the content.
    if (luminosityCount > 0) then
       call Stellar_Luminosities_Create(self)
       read (fileHandle) luminosityActiveCount
       call Dealloc_Array(self%luminosityValue                        )
       call Alloc_Array  (self%luminosityValue,[luminosityActiveCount])
       read (fileHandle) self%luminosityValue
    end if
    return
  end subroutine Stellar_Luminosities_Read_Raw

  subroutine Stellar_Luminosities_Reset(self)
    !% Reset an stellarLuminosities object.
    implicit none
    class(stellarLuminosities), intent(inout) :: self

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Ensure object is initialized.
    call Stellar_Luminosities_Create(self)
    ! Zero all properties.
    if (luminosityCount > 0) self%luminosityValue=0.0d0
    return
  end subroutine Stellar_Luminosities_Reset

  subroutine Stellar_Luminosities_Set_To_Unity(self)
    !% Set an stellarLuminosities object to unity.
    implicit none
    class(stellarLuminosities), intent(inout) :: self

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Ensure object is initialized.
    call Stellar_Luminosities_Create(self)
    ! Set values to unity.
    if (luminosityCount > 0) self%luminosityValue=1.0d0
    return
  end subroutine Stellar_Luminosities_Set_To_Unity

  logical function Stellar_Luminosities_Is_Zero(self)
    !% Test whether an stellarLuminosities object is zero.
    implicit none
    class(stellarLuminosities), intent(in   ) :: self

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Detect if all stellar luminosities are zero.
    Stellar_Luminosities_Is_Zero=.true.
    if (luminosityCount > 0 .and. allocated(self%luminosityValue)) then
       if (any(self%luminosityValue /= 0.0d0)) Stellar_Luminosities_Is_Zero=.false.
    end if
    return
  end function Stellar_Luminosities_Is_Zero

  double precision function Stellar_Luminosities_Luminosity(self,index)
    !% Return the requested luminosity from a {\normalfont \ttfamily stellarLuminosities} object.
    use Galacticus_Error
    implicit none
    class  (stellarLuminosities), intent(inout) :: self
    integer                     , intent(in   ) :: index

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Return the requested luminosity.
    if (allocated(self%luminosityValue)) then
       if (index > 0 .and. index <= size(self%luminosityValue)) then
          Stellar_Luminosities_Luminosity=self%luminosityValue(index)
       else
          call Galacticus_Error_Report('Stellar_Luminosities_Luminosity','index out of range')
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
    !% Return an element-by-element {\normalfont \ttfamily max()} on two stellar luminosity objects.
    implicit none
    type   (stellarLuminosities)                :: stellarLuminositiesMax
    type   (stellarLuminosities), intent(in   ) :: luminosities1         , luminosities2
    integer                                     :: luminosityCountActual

    if (luminosityCount > 0) then
       luminosityCountActual=stellarLuminositiesCountMaximum(luminosities1,luminosities2)
       stellarLuminositiesMax%luminosityValue=spread(-HUGE(0.0d0),1,luminosityCountActual)
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

  function Stellar_Luminosities_Add(luminosities1,luminosities2)
    !% Add two stellar luminosities objects.
    implicit none
    type   (stellarLuminosities)                          :: Stellar_Luminosities_Add
    class  (stellarLuminosities), intent(in   )           :: luminosities1
    class  (stellarLuminosities), intent(in   ), optional :: luminosities2
    integer                                               :: luminosityCountActual

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
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
    !% Increment a stellar luminosities object.
    implicit none
    class  (stellarLuminosities), intent(inout) :: self
    class  (stellarLuminosities), intent(in   ) :: increment
    integer                                     :: luminosityCountActual

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
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
    !% Subtract two stellar luminosities objects.
    implicit none
    type   (stellarLuminosities)                          :: Stellar_Luminosities_Subtract
    class  (stellarLuminosities), intent(in   )           :: luminosities1
    class  (stellarLuminosities), intent(in   ), optional :: luminosities2
    integer                                               :: luminosityCountActual

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
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
    !% Multiply a stellar luminosities object by a scalar.
    implicit none
    type            (stellarLuminosities)                :: Stellar_Luminosities_Multiply
    class           (stellarLuminosities), intent(in   ) :: stellarLuminosities1
    double precision                     , intent(in   ) :: multiplier

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    if (luminosityCount > 0) Stellar_Luminosities_Multiply%luminosityValue=stellarLuminosities1%luminosityValue*multiplier
    return
  end function Stellar_Luminosities_Multiply

  function Stellar_Luminosities_Multiply_Switched(multiplier,stellarLuminosities1)
    !% Multiply a stellar luminosities object by a scalar.
    implicit none
    type            (stellarLuminosities)                :: Stellar_Luminosities_Multiply_Switched
    class           (stellarLuminosities), intent(in   ) :: stellarLuminosities1
    double precision                     , intent(in   ) :: multiplier

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    if (luminosityCount > 0) Stellar_Luminosities_Multiply_Switched%luminosityValue=stellarLuminosities1%luminosityValue*multiplier
    return
  end function Stellar_Luminosities_Multiply_Switched

  function Stellar_Luminosities_Divide(stellarLuminosities1,divisor)
    !% Divide a stellar luminosities object by a scalar.
    implicit none
    type            (stellarLuminosities)                :: Stellar_Luminosities_Divide
    class           (stellarLuminosities), intent(in   ) :: stellarLuminosities1
    double precision                     , intent(in   ) :: divisor

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    if (luminosityCount > 0) Stellar_Luminosities_Divide%luminosityValue=stellarLuminosities1%luminosityValue/divisor
    return
  end function Stellar_Luminosities_Divide

  integer function Stellar_Luminosities_Property_Count(unmapped)
    !% Return the number of properties required to track stellar luminosities.
    implicit none
    logical, intent(in   ), optional :: unmapped

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Return the relevant count.
    if (present(unmapped).and.unmapped) then
       Stellar_Luminosities_Property_Count=luminosityCountUnmapped
    else
       Stellar_Luminosities_Property_Count=luminosityCount
    end if
    return
  end function Stellar_Luminosities_Property_Count

  integer function Stellar_Luminosities_Serialize_Count(self)
    !% Return the number of properties required to track stellar luminosities.
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
    !% Return a name for the specified entry in the stellar luminosities structure.
    use Galacticus_Error
    implicit none
    type   (varying_string)                :: Stellar_Luminosities_Name
    integer                , intent(in   ) :: index

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Check for index in range.
    if (index > 0 .and. index <= luminosityCount) then
       Stellar_Luminosities_Name=trim(luminosityName(index))
    else
       call Galacticus_Error_Report('Stellar_Luminosities_Name','index out of range')
    end if
    return
  end function Stellar_Luminosities_Name

  subroutine Stellar_Luminosities_Create(self)
    !% Ensure that the {\normalfont \ttfamily luminosity} array in a {\normalfont \ttfamily stellarLuminosities} is allocated.
    use Memory_Management
    implicit none
    type(stellarLuminosities), intent(inout) :: self

    if (.not.allocated(self%luminosityValue)) call Alloc_Array(self%luminosityValue,[luminosityCount])
    return
  end subroutine Stellar_Luminosities_Create

  subroutine Stellar_Luminosities_Deserialize(self,stellarLuminositiesArray)
    !% Pack stellar luminosities from an array into a {\normalfont \ttfamily stellarLuminosities} structure.
    implicit none
    class           (stellarLuminosities)              , intent(inout) :: self
    double precision                     , dimension(:), intent(in   ) :: stellarLuminositiesArray

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
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
    !% Unpack stellar luminosities from a {\normalfont \ttfamily stellarLuminosities} structure into an array.
    implicit none
    double precision                     , dimension(:), intent(  out) :: stellarLuminositiesArray(:)
    class           (stellarLuminosities)              , intent(in   ) :: self

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Place luminosities into array.
    if (allocated(self%luminosityValue)) then
       stellarLuminositiesArray(1:size(self%luminosityValue))=self%luminosityValue
    else
       stellarLuminositiesArray(1:luminosityCount           )=0.0d0
    end if
    return
  end subroutine Stellar_Luminosities_Serialize

  subroutine Stellar_Luminosities_Output(self,integerProperty,integerBufferCount,integerBuffer,doubleProperty,doubleBufferCount&
       &,doubleBuffer,time)
    !% Store a {\normalfont \ttfamily stellarLuminosities} object in the output buffers.
    use Kind_Numbers
    use Memory_Management
    implicit none
    class           (stellarLuminosities)                , intent(inout) :: self
    double precision                                     , intent(in   ) :: time
    integer                                              , intent(inout) :: doubleBufferCount, doubleProperty, integerBufferCount, &
         &                                                                  integerProperty
    integer         (kind=kind_int8     ), dimension(:,:), intent(inout) :: integerBuffer
    double precision                     , dimension(:,:), intent(inout) :: doubleBuffer
    double precision                     , dimension(:  ), allocatable   :: luminosityTmp
    integer                                                              :: i                , luminosityRemainingCount
    
    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    if (luminosityCount > 0) then
       do i=1,luminosityCount
          if (Stellar_Luminosities_Is_Output(i,time)) then
             doubleBuffer(doubleBufferCount,doubleProperty+1:doubleProperty+1)=self%luminosityValue(i)
             doubleProperty=doubleProperty+1
          end if
       end do
       ! Test if we can prune luminosities.
       select case (luminosityOutputOption)
       case (luminosityOutputOptionFuture,luminosityOutputOptionPresent)
          ! Luminosities from this and earlier outputs no longer needed, so prune them.
          call Move_Alloc(self%luminosityValue,luminosityTmp)
          luminosityRemainingCount=size(luminosityTmp)
          do i=1,luminosityCount
             if (Stellar_Luminosities_Is_Output(i,time,luminosityOutputOptionPresent)) &
                  & luminosityRemainingCount=luminosityRemainingCount-1
          end do          
          call Alloc_Array(self%luminosityValue,[luminosityRemainingCount])
          self%luminosityValue=luminosityTmp(1:luminosityRemainingCount)
          call Dealloc_Array(luminosityTmp)
       end select
    end if
    return
  end subroutine Stellar_Luminosities_Output

  subroutine Stellar_Luminosities_Output_Count(self,integerPropertyCount,doublePropertyCount,time)
    !% Increment the output count to account for a {\normalfont \ttfamily stellarLuminosities} object.
    implicit none
    class           (stellarLuminosities), intent(in   ) :: self
    integer                              , intent(inout) :: doublePropertyCount, integerPropertyCount
    double precision                     , intent(in   ) :: time

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    doublePropertyCount=doublePropertyCount+Stellar_Luminosities_Output_Count_Get(time)
    return
  end subroutine Stellar_Luminosities_Output_Count

  integer function Stellar_Luminosities_Output_Count_Get(time)
    !% Compute the number of luminosities to be output at a given time.
    implicit none
    double precision, intent(in   ) :: time
    integer                         :: i

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    Stellar_Luminosities_Output_Count_Get=0
    do i=1,luminosityCount
       if (Stellar_Luminosities_Is_Output(i,time)) Stellar_Luminosities_Output_Count_Get=Stellar_Luminosities_Output_Count_Get+1
    end do
    return
  end function Stellar_Luminosities_Output_Count_Get

  subroutine Stellar_Luminosities_Output_Names(self,integerProperty,integerPropertyNames,integerPropertyComments,integerPropertyUnitsSI&
       &,doubleProperty,doublePropertyNames,doublePropertyComments,doublePropertyUnitsSI,time,prefix,comment,unitsInSI)
    !% Assign names to output buffers for a {\normalfont \ttfamily stellarLuminosities} object.
    implicit none
    class           (stellarLuminosities)              , intent(in   ) :: self
    double precision                                   , intent(in   ) :: time
    integer                                            , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*              ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                                integerPropertyComments, integerPropertyNames
    double precision                     , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    character       (len=*              )              , intent(in   ) :: comment                , prefix
    double precision                                   , intent(in   ) :: unitsInSI
    integer                                                            :: i

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    if (luminosityCount > 0) then
       do i=1,luminosityCount
          if (Stellar_Luminosities_Is_Output(i,time)) then
             doubleProperty=doubleProperty+1
             doublePropertyNames   (doubleProperty)=trim(prefix )// ':'//trim(luminosityName(i))
             doublePropertyComments(doubleProperty)=trim(comment)//' ['//trim(luminosityName(i))//']'
             doublePropertyUnitsSI (doubleProperty)=unitsInSI
          end if
       end do
    end if
    return
  end subroutine Stellar_Luminosities_Output_Names

  logical function Stellar_Luminosities_Is_Output(luminosityIndex,time,outputOption)
    !% Return true or false depending on whether {\normalfont \ttfamily luminosityIndex} should be output at {\normalfont \ttfamily time}.
    use Galacticus_Error
    implicit none
    integer         , intent(in   )           :: luminosityIndex
    double precision, intent(in   )           :: time
    integer         , intent(in   ), optional :: outputOption
    double precision, parameter               :: timeTolerance  =1.0d-3
    integer                                   :: outputOptionActual

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
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
       Stellar_Luminosities_Is_Output=(     luminosityCosmicTime(luminosityIndex)       >= time*(1.0d0-timeTolerance))
    case (luminosityOutputOptionPresent)
       Stellar_Luminosities_Is_Output=(abs(luminosityCosmicTime(luminosityIndex)-time) <= time*       timeTolerance )
    case default
       call Galacticus_Error_Report('Stellar_Luminosities_Is_Output','unknown luminosity output option')
    end select
    return
  end function Stellar_Luminosities_Is_Output

  subroutine Stellar_Luminosities_Set(self,mass,imfSelected,time,abundancesStellar)
    !% Set the luminosity in each band for a single stellar population of given {\normalfont \ttfamily mass} with the specified {\normalfont \ttfamily abundancesStellar} and
    !% which formed at cosmological {\normalfont \ttfamily time} with IMF specified by {\normalfont \ttfamily imfSelected}.
    use Abundances_Structure
    use Stellar_Population_Luminosities
    implicit none
    class           (stellarLuminosities)                             :: self
    integer                              , intent(in   )              :: imfSelected
    double precision                     , intent(in   )              :: mass             , time
    type            (abundances         ), intent(in   )              :: abundancesStellar
    double precision                     , dimension(luminosityCount) :: ages             , massToLightRatio

    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()

    ! Return if no luminosities are tracked.
    if (luminosityCount == 0) return

    ! Get the ages that this stellar population will have at the various output times.
    ages=luminosityCosmicTime-time

    ! Get the luminosities for each requested band.
    massToLightRatio=Stellar_Population_Luminosity(                                    &
         &                                         luminosityIndex                   , &
         &                                         luminosityFilterIndex             , &
         &                                         luminosityPostprocessingChainIndex, &
         &                                         imfSelected                       , &
         &                                         abundancesStellar                 , &
         &                                         ages                              , &
         &                                         luminosityBandRedshift              &
         &                                        )
    call Stellar_Luminosities_Create(self)
    self%luminosityValue=mass*massToLightRatio(1:size(self%luminosityValue))
    return
  end subroutine Stellar_Luminosities_Set

  integer function Stellar_Luminosities_Index(name)
    !% Return the index of and specified entry in the luminosity list given its name.
    use Galacticus_Error
    implicit none
    type   (varying_string), intent(in   ) :: name
    integer                                :: i

    do i=1,luminosityCount
       if (name == luminosityName(i)) then
          Stellar_Luminosities_Index=i
          return
       end if
    end do
    call Galacticus_Error_Report('Stellar_Population_Luminosities_Index','unmatched name')
    return
  end function Stellar_Luminosities_Index

  subroutine Stellar_Luminosities_Special_Cases(luminosityMap,luminosityRedshiftText,luminosityRedshift,luminosityBandRedshift,luminosityFilter,luminosityType,luminosityPostprocessSet)
    !% Modify the input list of luminosities for special cases.
    use, intrinsic :: ISO_C_Binding
    use Galacticus_Output_Times
    use Cosmology_Functions
    use Memory_Management
    use String_Handling
    implicit none
    integer                                  , intent(inout), allocatable, dimension(:) :: luminosityMap
    type            (varying_string         ), intent(inout), allocatable, dimension(:) :: luminosityRedshiftText   , luminosityFilter           , &
         &                                                                                 luminosityType           , luminosityPostprocessSet
    double precision                         , intent(inout), allocatable, dimension(:) :: luminosityRedshift       , luminosityBandRedshift
    integer         (c_size_t               )                                           :: j                        , outputCount                , &
         &                                                                                 i                        , newFilterCount
    integer                                                 , allocatable, dimension(:) :: luminosityMapTmp
    type            (varying_string         )               , allocatable, dimension(:) :: luminosityRedshiftTextTmp, luminosityFilterTmp        , &
         &                                                                                 luminosityTypeTmp        , luminosityPostprocessSetTmp
    type            (varying_string         )                            , dimension(4) :: specialFilterWords
    double precision                                        , allocatable, dimension(:) :: luminosityRedshiftTmp    , luminosityBandRedshiftTmp
    class           (cosmologyFunctionsClass), pointer                                  :: cosmologyFunctions_
    character       (len=32                 )                                           :: redshiftLabel            , word                     , &
         &                                                                                 wavelengthCentralLabel   , resolutionLabel
    character       (len=128                )                                           :: newFilterName
    double precision                                                                    :: outputRedshift           , resolution               , &
         &                                                                                 wavelengthMinimum        , wavelengthMaximum        , &
         &                                                                                 wavelengthRatio          , wavelengthCentral
    ! Get cosmology functions.
    cosmologyFunctions_ => cosmologyFunctions()
    ! Get number of output redshifts.
    outputCount=Galacticus_Output_Time_Count()
    ! Iterate over all luminosities.
    i=1
    do while (i <= size(luminosityRedshiftText))
       ! Check for special cases.
       if (luminosityRedshiftText(i) == "all") then
          ! Resize the arrays.
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
             outputRedshift=Galacticus_Output_Redshift(j)
             write (redshiftLabel,*) outputRedshift
             luminosityRedshiftText   (j+i-1)=redshiftLabel
             luminosityRedshift       (j+i-1)=outputRedshift
             if (luminosityBandRedshiftTmp  (i) <= -2.0d0) then
                luminosityBandRedshift(j+i-1)=outputRedshift
             else
                luminosityBandRedshift(j+i-1)=luminosityBandRedshiftTmp  (i)
             end if
          end do
          deallocate        (luminosityRedshiftTextTmp  )
          deallocate        (luminosityFilterTmp        )
          deallocate        (luminosityTypeTmp          )
          deallocate        (luminosityPostprocessSetTmp)     
          call Dealloc_Array(luminosityMapTmp           )
          call Dealloc_Array(luminosityRedshiftTmp      )
          call Dealloc_Array(luminosityBandRedshiftTmp  )
       end if
       ! Arrays of top-hat filters.
       if (extract(luminosityFilter(i),1,12) == "topHatArray_") then
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
             write (newFilterName,'(a,a,a,a)') "topHat_",trim(adjustl(wavelengthCentralLabel)),"_",trim(adjustl(resolutionLabel))
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
          deallocate        (luminosityRedshiftTextTmp  )
          deallocate        (luminosityFilterTmp        )
          deallocate        (luminosityTypeTmp          )
          deallocate        (luminosityPostprocessSetTmp)     
          call Dealloc_Array(luminosityMapTmp           )
          call Dealloc_Array(luminosityRedshiftTmp      )
          call Dealloc_Array(luminosityBandRedshiftTmp  )
       end if
       ! Next luminosity.
       i=i+1
    end do
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
    !% Expand the filter set by removing the filter at index {\normalfont \ttfamily expandFrom} by adding {\normalfont \ttfamily expandCount} replicas of the filter at that point.
    use, intrinsic :: ISO_C_Binding
    use Memory_Management
    implicit none
    integer         (c_size_t      ), intent(in   )                            :: expandFrom               , expandCount
    integer                         , intent(inout), allocatable, dimension(:) :: luminosityMap
    type            (varying_string), intent(inout), allocatable, dimension(:) :: luminosityRedshiftText   , luminosityFilter           , &
         &                                                                        luminosityType           , luminosityPostprocessSet
    double precision                , intent(inout), allocatable, dimension(:) :: luminosityRedshift       , luminosityBandRedshift
    integer                         , intent(inout), allocatable, dimension(:) :: luminosityMapTmp
    type            (varying_string), intent(inout), allocatable, dimension(:) :: luminosityRedshiftTextTmp, luminosityFilterTmp        , &
         &                                                                        luminosityTypeTmp        , luminosityPostprocessSetTmp
    double precision                , intent(inout), allocatable, dimension(:) :: luminosityRedshiftTmp    , luminosityBandRedshiftTmp
    integer                                                                    :: luminosityCount

    luminosityCount=size(luminosityMap)
    call Move_Alloc (luminosityMap           ,luminosityMapTmp           )
    call Move_Alloc (luminosityRedshiftText  ,luminosityRedshiftTextTmp  )
    call Move_Alloc (luminosityRedshift      ,luminosityRedshiftTmp      )
    call Move_Alloc (luminosityBandRedshift  ,luminosityBandRedshiftTmp  )
    call Move_Alloc (luminosityFilter        ,luminosityFilterTmp        )
    call Move_Alloc (luminosityType          ,luminosityTypeTmp          )
    call Move_Alloc (luminosityPostprocessSet,luminosityPostprocessSetTmp)
    allocate        (luminosityRedshiftText   (size(luminosityRedshiftTextTmp  )+expandCount-1))
    allocate        (luminosityFilter         (size(luminosityFilterTmp        )+expandCount-1))
    allocate        (luminosityType           (size(luminosityTypeTmp          )+expandCount-1))
    allocate        (luminosityPostprocessSet (size(luminosityPostprocessSetTmp)+expandCount-1))
    call Alloc_Array(luminosityMap           ,[size(luminosityRedshiftTmp      )+expandCount-1])
    call Alloc_Array(luminosityRedshift      ,[size(luminosityRedshiftTmp      )+expandCount-1])
    call Alloc_Array(luminosityBandRedshift  ,[size(luminosityBandRedshiftTmp  )+expandCount-1])
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
    !% Truncate (or pad) the stellar luminosities to match the number in the given {\normalfont \ttfamily templateLuminosities}.
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
    !% Map an array of luminosity-related input parameters into a new array accounting for special case processing.
    use Memory_Management
    implicit none
    double precision, intent(inout), allocatable, dimension(:) :: parameters
    double precision               , allocatable, dimension(:) :: parametersMapped
    integer                                                    :: i
    
    ! Ensure module is initialized.
    call Stellar_Luminosities_Initialize()
    ! Allocate new array.
    call Alloc_Array(parametersMapped,[luminosityCount])
    ! Map from the old array.
    do i=1,luminosityCount
       parametersMapped(i)=parameters(luminosityMap(i))
    end do
    ! Copy the new array.
    call Dealloc_Array(parameters)
    call Move_Alloc(parametersMapped,parameters)
    return
  end subroutine Stellar_Luminosities_Parameter_Map_Double
  
end module Stellar_Luminosities_Structure
