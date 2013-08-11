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

!% Contains a module which implements functionality related to the stellar initial mass function.

module Star_Formation_IMF
  !% Implements functionality related to the stellar initial mass function.
  use ISO_Varying_String
  use Abundances_Structure
  !# <include directive="imfRegister" type="moduleUse">
  include 'star_formation.IMF.register.modules.inc'
  !# </include>
  !# <include directive="imfSelectionMethod" type="moduleUse">
  !#  <functionArgs>imfSelectionMethod,IMF_Select,imfNames</functionArgs>
  include 'star_formation.IMF.select.modules.inc'
  !# </include>
  implicit none
  private
  public :: IMF_Select, IMF_Recycled_Fraction_Instantaneous, IMF_Recycling_Rate_NonInstantaneous, IMF_Yield_Instantaneous,&
       & IMF_Metal_Yield_Rate_NonInstantaneous, IMF_Energy_Input_Rate_NonInstantaneous, IMF_Name, IMF_Tabulate, IMF_Descriptor

  ! Flag to indicate if this module has been initialized.
  logical                                                                           :: imfInitialized                         =.false.

  ! Count of the number of available IMFs.
  integer                                                                           :: imfAvailableCount                      =0

  ! Array of IMF names.
  type            (varying_string     )           , allocatable, dimension(:      ) :: imfDescriptors                                 , imfNames

  ! Tables of recycled fractions.
  logical                                         , allocatable, dimension(:      ) :: recycledFractionTabulated
  integer                                         , allocatable, dimension(:      ) :: recycledFractionIndex
  double precision                                , allocatable, dimension(:      ) :: recycledFractionTableAge                       , recycledFractionTableMetallicity
  double precision                                , allocatable, dimension(:,:,:  ) :: recycledFractionTable
  integer                              , parameter                                  :: recycledFractionTableMetallicityCount  =10
  integer                              , parameter                                  :: recycledFractionTableAgeCount          =50
  double precision                     , parameter                                  :: recycledFractionTableMetallicityMinimum=1.0d-4
  double precision                     , parameter                                  :: recycledFractionTableMetallicityMaximum=0.6d-1
  double precision                     , parameter                                  :: recycledFractionTableAgeMinimum        =1.0d-3
  double precision                     , parameter                                  :: recycledFractionTableAgeMaximum        =1.0d+2

  ! Tables of metal yields fractions.
  logical                                         , allocatable, dimension(:      ) :: metalYieldTabulated
  integer                                         , allocatable, dimension(:      ) :: metalYieldIndex
  double precision                                , allocatable, dimension(:      ) :: metalYieldTableAge                             , metalYieldTableMetallicity
  double precision                                , allocatable, dimension(:,:,:,:) :: metalYieldTable
  integer                              , parameter                                  :: metalYieldTableMetallicityCount        =10
  integer                              , parameter                                  :: metalYieldTableAgeCount                =50
  double precision                     , parameter                                  :: metalYieldTableMetallicityMinimum      =1.0d-4
  double precision                     , parameter                                  :: metalYieldTableMetallicityMaximum      =0.6d-1
  double precision                     , parameter                                  :: metalYieldTableAgeMinimum              =1.0d-3
  double precision                     , parameter                                  :: metalYieldTableAgeMaximum              =1.0d+2

  ! Tables of cumulative energy inputs.
  logical                                         , allocatable, dimension(:    )   :: energyInputTabulated
  integer                                         , allocatable, dimension(:    )   :: energyInputIndex
  double precision                                , allocatable, dimension(:    )   :: energyInputTableAge                            , energyInputTableMetallicity
  double precision                                , allocatable, dimension(:,:,:)   :: energyInputTable
  integer                              , parameter                                  :: energyInputTableMetallicityCount       =10
  integer                              , parameter                                  :: energyInputTableAgeCount               =50
  double precision                     , parameter                                  :: energyInputTableMetallicityMinimum     =1.0d-4
  double precision                     , parameter                                  :: energyInputTableMetallicityMaximum     =0.6d-1
  double precision                     , parameter                                  :: energyInputTableAgeMinimum             =1.0d-3
  double precision                     , parameter                                  :: energyInputTableAgeMaximum             =1.0d+2

  ! Module global variables used in integration.
  integer                                                                           :: atomIndexGlobal                                , imfSelectedGlobal
  !$omp threadprivate(atomIndexGlobal,imfSelectedGlobal)
  double precision                                                                  :: lifetime                                       , metallicity
  !$omp threadprivate(lifetime,metallicity)

  ! Count of number of individual elements tracked.
  integer                                                                           :: elementCount

  ! Pointer to the function that selects which IMF to use.
  procedure       (IMF_Select_Template), pointer                                    :: IMF_Select_Do                          =>null()
  abstract interface
     integer function IMF_Select_Template(starFormationRate,fuelAbundances,component)
       import abundances
       double precision            , intent(in   ) :: starFormationRate
       type            (abundances), intent(in   ) :: fuelAbundances
       integer                     , intent(in   ) :: component
     end function IMF_Select_Template
  end interface

  ! Current file format used for stellar population properties.
  integer                         , parameter :: fileFormatCurrent                                      =1

  ! A unique label for output files.
  type            (varying_string)            :: imfUniqueLabel

  ! Options controlling instantaneous stellar evolution approximations.
  logical                                     :: starFormationImfInstantaneousApproximation
  double precision                            :: starFormationImfInstantaneousApproximationMassLongLived
  double precision                            :: starFormationImfInstantaneousApproximationEffectiveAge

contains

  integer function IMF_Select(starFormationRate,fuelAbundances,component)
    !% Selects an IMF give an input {\tt starFormationRate} and {\tt fuelAbundances}.
    implicit none
    double precision            , intent(in   ) :: starFormationRate
    type            (abundances), intent(in   ) :: fuelAbundances
    integer                     , intent(in   ) :: component

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Call the function that makes the selection.
    IMF_Select=IMF_Select_Do(starFormationRate,fuelAbundances,component)
    return
  end function IMF_Select

  function IMF_Name(imfIndex)
    !% Return the name of the IMF with the specified index.
    use Galacticus_Error
    implicit none
    type   (varying_string)                :: IMF_Name
    integer                , intent(in   ) :: imfIndex

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    if (imfIndex <= size (imfNames)) then
       IMF_Name=imfNames(imfIndex)
    else
       call Galacticus_Error_Report('IMF_Name','imfIndex is out of range')
    end if
    return
  end function IMF_Name

  function IMF_Descriptor(imfIndex)
    !% Return a full descriptor for the IMF with the specified index.
    use Galacticus_Error
    implicit none
    type   (varying_string)                :: IMF_Descriptor
    integer                , intent(in   ) :: imfIndex

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    if (imfIndex <= size (imfDescriptors)) then
       IMF_Descriptor=imfDescriptors(imfIndex)
    else
       call Galacticus_Error_Report('IMF_Descriptor','imfIndex is out of range')
    end if
    return
  end function IMF_Descriptor

  subroutine Star_Formation_IMF_Initialize
    !% Initialize the IMF subsystem.
    use Memory_Management
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string) :: imfSelectionMethod

    ! Initialize the IMF subsystem if necessary.
    if (.not.imfInitialized) then
       !$omp critical(IMF_Initialize)
       if (.not.imfInitialized) then
          ! Register all available IMFs.
          !# <include directive="imfRegister" type="functionCall" functionType="void">
          !#  <functionArgs>imfAvailableCount</functionArgs>
          include 'star_formation.IMF.register.inc'
          !# </include>

          ! Get a list of IMF names and descriptors.
          allocate(imfNames      (imfAvailableCount))
          allocate(imfDescriptors(imfAvailableCount))
          call Memory_Usage_Record(sizeof(imfNames)+sizeof(imfDescriptors),file=__FILE__,line=__LINE__)
          !# <include directive="imfRegisterName" type="functionCall" functionType="void">
          !#  <functionArgs>imfNames,imfDescriptors</functionArgs>
          include 'star_formation.IMF.register_names.inc'
          !# </include>

          ! Register the IMF selection method.
          !@ <inputParameter>
          !@   <name>imfSelectionMethod</name>
          !@   <defaultValue>fixed</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The name of the method to be used for selecting which \gls{imf} to use.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('imfSelectionMethod',imfSelectionMethod,defaultValue='fixed')
          ! Include file that makes calls to all available method initialization routines.
          !# <include directive="imfSelectionMethod" type="functionCall" functionType="void">
          !#  <functionArgs>imfSelectionMethod,IMF_Select_Do,imfNames</functionArgs>
          include 'star_formation.IMF.select.inc'
          !# </include>
          if (.not.associated(IMF_Select_Do)) call Galacticus_Error_Report('Star_Formation_IMF_Initialize'&
               &,'method '//char(imfSelectionMethod)//' is unrecognized')

          ! Get options controlling the instantaneous stellar evolution approximation.
          !@ <inputParameter>
          !@   <name>starFormationImfInstantaneousApproximation</name>
          !@   <defaultValue>false</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     Option controlling whether stellar evolution should follow the instantaneous approximation.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationImfInstantaneousApproximation',starFormationImfInstantaneousApproximation,defaultValue=.false.)
          !@ <inputParameter>
          !@   <name>starFormationImfInstantaneousApproximationMassLongLived</name>
          !@   <defaultValue>$1M_\odot$</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The mass below which stars are assumed to be infinitely long-lived in the instantaneous approximation for stellar evolution.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationImfInstantaneousApproximationMassLongLived',starFormationImfInstantaneousApproximationMassLongLived,defaultValue=1.0d0)
          !@ <inputParameter>
          !@   <name>starFormationImfInstantaneousApproximationEffectiveAge</name>
          !@   <defaultValue>$13.8$Gyr</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The effective age to use for computing SNeIa yield when using the instantaneous stellar evolution approximation.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@   <group>initialMassFunction</group>
          !@ </inputParameter>
          call Get_Input_Parameter('starFormationImfInstantaneousApproximationEffectiveAge',starFormationImfInstantaneousApproximationEffectiveAge,defaultValue=13.8d0)

          ! Get a count of the number of individual elements that must be tracked.
          elementCount=Abundances_Property_Count()

          ! Get a unique label for IMF files. Ignore all parameters that define the choice of IMF and element since we will write
          ! these to separate files anyway.
          !# <uniqueLabel>
          !#  <function>Star_Formation_IMF_Label</function>
          !#  <ignore>elementsToTrack</ignore>
          !#  <ignore>imfSelectionMethod</ignore>
          !#  <ignore>imfPiecewisePowerLawExponents</ignore>
          !#  <ignore>imfPiecewisePowerLawMassPoints</ignore>
          !#  <ignoreRegex>^imfSelection.*</ignoreRegex>
          !#  <ignoreRegex>^imf.*RecycledInstantaneous$</ignoreRegex>
          !#  <ignoreRegex>^imf.*YieldInstantaneous$</ignoreRegex>
          !# </uniqueLabel>
          imfUniqueLabel=Star_Formation_IMF_Label(includeVersion=.true.,asHash=.true.)
          ! Flag that the module is now initialized.
          imfInitialized=.true.
       end if
       !$omp end critical(IMF_Initialize)
    end if
    return
  end subroutine Star_Formation_IMF_Initialize

  double precision function IMF_Recycled_Fraction_Instantaneous(starFormationRate,fuelAbundances,component)
    !% Returns a recycled fraction for the IMF suitable for use in the instantaneous recycling approximation.
    implicit none
    double precision            , intent(in   ) :: starFormationRate
    type            (abundances), intent(in   ) :: fuelAbundances
    integer                     , intent(in   ) :: component
    integer                                     :: imfSelected
    logical                                     :: imfMatched

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    ! Get the recycled fraction from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfRecycledInstantaneous" type="functionCall" functionType="void">
    !#  <functionArgs>imfSelected,imfMatched,IMF_Recycled_Fraction_Instantaneous</functionArgs>
    !#  <onReturn>if (imfMatched) return</onReturn>
    include 'star_formation.IMF.recycled_instantaneous.inc'
    !# </include>
    return
  end function IMF_Recycled_Fraction_Instantaneous

  double precision function IMF_Yield_Instantaneous(starFormationRate,fuelAbundances,component)
    !% Returns a yield for the IMF suitable for use in the instantaneous recycling approximation.
    implicit none
    double precision            , intent(in   ) :: starFormationRate
    type            (abundances), intent(in   ) :: fuelAbundances
    integer                     , intent(in   ) :: component
    integer                                     :: imfSelected
    logical                                     :: imfMatched

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    ! Get the yield from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfYieldInstantaneous" type="functionCall" functionType="void">
    !#  <functionArgs>imfSelected,imfMatched,IMF_Yield_Instantaneous</functionArgs>
    !#  <onReturn>if (imfMatched) return</onReturn>
    include 'star_formation.IMF.yield_instantaneous.inc'
    !# </include>
    return
  end function IMF_Yield_Instantaneous

  subroutine IMF_Tabulate(imfIndex,imf)
    !% Returns a tabulation of the IMF with sufficient resolution to resolve all features.
    use Tables
    implicit none
    class  (table1D), allocatable, intent(inout) :: imf
    integer                      , intent(in   ) :: imfIndex
    logical                                      :: imfMatched

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Get the recycled fraction from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfTabulate" type="functionCall" functionType="void">
    !#  <functionArgs>imfIndex,imfMatched,imf</functionArgs>
    !#  <onReturn>if (imfMatched) return</onReturn>
    include 'star_formation.IMF.tabulate.inc'
    !# </include>
    return
  end subroutine IMF_Tabulate

  double precision function IMF_Minimum_Mass(imfSelected)
    !% Returns the minimum mass in the selected IMF.
    implicit none
    integer, intent(in   ) :: imfSelected
    logical                :: imfMatched

    ! Get the minimum mass from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfMinimumMass" type="functionCall" functionType="void">
    !#  <functionArgs>imfSelected,imfMatched,IMF_Minimum_Mass</functionArgs>
    !#  <onReturn>if (imfMatched) return</onReturn>
    include 'star_formation.IMF.minimum_mass.inc'
    !# </include>
    return
  end function IMF_Minimum_Mass

  double precision function IMF_Maximum_Mass(imfSelected)
    !% Returns the maximum mass in the selected IMF.
    implicit none
    integer, intent(in   ) :: imfSelected
    logical                :: imfMatched

    ! Get the maximum mass from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfMaximumMass" type="functionCall" functionType="void">
    !#  <functionArgs>imfSelected,imfMatched,IMF_Maximum_Mass</functionArgs>
    !#  <onReturn>if (imfMatched) return</onReturn>
    include 'star_formation.IMF.maximum_mass.inc'
    !# </include>
    return
  end function IMF_Maximum_Mass

  double precision function IMF_Phi(initialMass,imfSelected)
    !% Returns the IMF, $\Phi(M)$, at mass $M=${\tt initialMass} for the selected IMF.
    implicit none
    integer         , intent(in   ) :: imfSelected
    double precision, intent(in   ) :: initialMass
    logical                         :: imfMatched

    ! Get the IMF for the selected IMF and initial stellar mass.
    imfMatched=.false.
    !# <include directive="imfPhi" type="functionCall" functionType="void">
    !#  <functionArgs>imfSelected,imfMatched,initialMass,IMF_Phi</functionArgs>
    !#  <onReturn>if (imfMatched) return</onReturn>
    include 'star_formation.IMF.phi.inc'
    !# </include>
    return
  end function IMF_Phi

  double precision function IMF_Recycling_Rate_NonInstantaneous(starFormationRate,fuelAbundances,component,ageMinimum,ageMaximum)
    !% Returns the recycling rate for a simple stellar population. The \gls{imf} is determined from the given {\tt starFormationRate}
    !% and {\tt fuelAbundances}. The recycling rate (in the fraction of the population's mass returned to the \gls{ism} per Gyr) is
    !% computed for the given {\tt age} (in Gyr). The recycled fraction is computed on a grid of age and metallicity. This is
    !% stored to file and will be read back in on subsequent runs. This is useful as computation of the table is relatively slow.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Ranges
    use Numerical_Constants_Astronomical
    use Memory_Management
    use Galacticus_Display
    use File_Utilities
    use FoX_wxml
    use FoX_dom
    use IO_XML
    use Galacticus_Error
    use Dates_and_Times
    use Galacticus_Input_Paths
    implicit none
    double precision                                                                          , intent(in   ) ::        ageMinimum                                         , starFormationRate
    double precision                                                                          , intent(in   ) , optional::                                ageMaximum
    type            (abundances                )                                              , intent(in   ) ::        fuelAbundances
    integer                                                                                   , intent(in   ) ::        component
    logical                                                    , allocatable, dimension(:    )                ::        recycledFractionTabulatedTemporary
    integer                                                    , allocatable, dimension(:    )                ::        recycledFractionIndexTemporary
    double precision                                           , allocatable, dimension(:    )                ::        tableTemporary
    double precision                                           , allocatable, dimension(:,:,:)                ::        recycledFractionTableTemporary
    double precision                                                        , dimension(2    )                ::        metallicityFactors                                 , recycleRate
    type            (fgsl_interp               )         , save                                               ::        interpolationAgeObject
    type            (fgsl_interp_accel         )         , save                                               ::        interpolationAgeAccelerator                        , interpolationMetallicityAccelerator
    logical                                              , save                                               ::        interpolationAgeReset                       =.true., interpolationMetallicityReset      =.true.
    !$omp threadprivate(interpolationAgeObject,interpolationMetallicityAccelerator &
    !$omp ,interpolationAgeAccelerator,interpolationMetallicityReset,interpolationAgeReset)
    type            (Node                      ), pointer                                                     ::        doc                                                , thisItem
    type            (NodeList                  ), pointer                                                     ::        dataList
    type            (c_ptr                     )                                                              ::        parameterPointer
    type            (fgsl_function             )                                                              ::        integrandFunction
    type            (fgsl_integration_workspace)                                                              ::        integrationWorkspace
    integer                                                                                                   ::        fileFormat                                         , iAge                                      , &
         &                                                                                                              iMetallicity                                       , iRecycledFraction                         , &
         &                                                                                                              imfCount                                           , imfSelected                               , &
         &                                                                                                              ioErr                                              , loopCount                                 , &
         &                                                                                                              loopCountTotal                                     , metallicityIndex                          , &
         &                                                                                                              tableIndex
    double precision                                                                                          ::        maximumMass                                        , minimumMass                               , &
         &                                                                                                              recycledFractionMaximum                            , recycledFractionMinimum
    character       (len=20                    )                                                              ::        parameterValue                                     , progressMessage
    type            (xmlf_t                    )                                                              ::        recycledFractionDoc
    type            (varying_string            )                                                              ::        fileName
    logical                                                                                                   ::        makeFile

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    !$omp critical(IMF_Recycling_Rate_NonInstantaneous_Initialize)
    ! Check that flag and index arrays exist.
    if (.not.allocated(recycledFractionTabulated)) then
       call Alloc_Array(recycledFractionTabulated,[imfSelected],file=__FILE__,line=__LINE__)
       call Alloc_Array(recycledFractionIndex    ,[imfSelected],file=__FILE__,line=__LINE__)
       recycledFractionTabulated=.false.
       recycledFractionIndex    =0
    end if

    ! Check that flag and index arrays are large enough.
    if (size(recycledFractionTabulated) < imfSelected) then
       call Move_Alloc (recycledFractionTabulated,recycledFractionTabulatedTemporary)
       call Move_Alloc (recycledFractionIndex    ,recycledFractionIndexTemporary    )
       call Alloc_Array(recycledFractionTabulated,[imfSelected],file=__FILE__,line=__LINE__)
       call Alloc_Array(recycledFractionIndex    ,[imfSelected],file=__FILE__,line=__LINE__)
       recycledFractionTabulated(1:size(recycledFractionTabulatedTemporary))            =recycledFractionTabulatedTemporary
       recycledFractionIndex    (1:size(recycledFractionIndexTemporary    ))            =recycledFractionIndexTemporary
       recycledFractionTabulated(  size(recycledFractionTabulatedTemporary):imfSelected)=.false.
       recycledFractionIndex    (  size(recycledFractionTabulatedTemporary):imfSelected)=0
       call Dealloc_Array(recycledFractionTabulatedTemporary,file=__FILE__,line=__LINE__)
       call Dealloc_Array(recycledFractionIndexTemporary    ,file=__FILE__,line=__LINE__)
    end if

    ! Tabulate the recycled fraction for this IMF if it has not already been computed.
    if (.not.recycledFractionTabulated(imfSelected)) then

       ! Expand the tabulations array by enough to accomodate a new IMF.
       if (allocated(recycledFractionTable)) then
          imfCount=size(recycledFractionTable,dim=3)
          call Move_Alloc(recycledFractionTable,recycledFractionTableTemporary)
          call Alloc_Array(recycledFractionTable,[recycledFractionTableAgeCount,recycledFractionTableMetallicityCount,imfCount],file=__FILE__,line=__LINE__)
          recycledFractionTable(:,:,1:imfCount)=recycledFractionTableTemporary
          call Dealloc_Array(recycledFractionTableTemporary,file=__FILE__,line=__LINE__)
       else
          call Alloc_Array(recycledFractionTableAge        ,[recycledFractionTableAgeCount        ],file=__FILE__,line=__LINE__)
          call Alloc_Array(recycledFractionTableMetallicity,[recycledFractionTableMetallicityCount],file=__FILE__,line=__LINE__)
          recycledFractionTableAge                                                 =Make_Range(recycledFractionTableAgeMinimum&
               &,recycledFractionTableAgeMaximum,recycledFractionTableAgeCount,rangeType=rangeTypeLogarithmic)
          recycledFractionTableMetallicity(1)                                      =0.0d0
          recycledFractionTableMetallicity(2:recycledFractionTableMetallicityCount)&
               &=Make_Range(recycledFractionTableMetallicityMinimum,recycledFractionTableMetallicityMaximum&
               &,recycledFractionTableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          call Alloc_Array(recycledFractionTable,[recycledFractionTableAgeCount,recycledFractionTableMetallicityCount,1],file=__FILE__,line=__LINE__)
       end if

       ! Record the index in the array where this IMF will be stored.
       recycledFractionIndex(imfSelected)=size(recycledFractionTable,dim=3)

       ! Check if the table has been computed and stored previously.
       fileName=char(Galacticus_Input_Path())//'data/stellarPopulations/Stellar_Recycled_Fraction_'//imfNames(imfSelected)//'_'//imfUniqueLabel//'.xml'
       makeFile=.false.
       if (File_Exists(fileName)) then
          ! Open the XML file containing energy input.
          call Galacticus_Display_Indent('Parsing file: '//fileName,verbosityDebug)
          doc => parseFile(char(fileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','Unable to parse recycled fractions file "'//fileName//'"'//char(10)//'HELP: file may be corrupted - try deleting this file and rerunning Galacticus')
          ! Check the version number.
          if (XML_Path_Exists(doc,"fileFormat")) then
             thisItem => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
             call extractDataContent(thisItem,fileFormat)
             if (fileFormat /= fileFormatCurrent) makeFile=.true.
          else
             makeFile=.true.
          end if
          if (makeFile) then
             call destroy(doc)
             call Galacticus_Display_Unindent('done',verbosityDebug)
          end if
       else
          makeFile=.true.
       end if

       if (.not.makeFile) then

          ! Find the ages element and extract data.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"ages")
          if (XML_Array_Length(thisItem,"data") /= recycledFractionTableAgeCount) call&
               & Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous' ,'ages array in XML file does not match internal&
               & expectation')
          call XML_Array_Read(thisItem,"data",tableTemporary)
          if (recycledFractionIndex(imfSelected) == 1) then
             recycledFractionTableAge=tableTemporary
          else
             if (any(recycledFractionTableAge /= tableTemporary)) call&
                  & Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous' ,'mismatch in ages array in XML file')
          end if
          deallocate(tableTemporary)

          ! Find the metallicities element and extract data.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"metallicities")
          dataList => getElementsByTagname(thisItem,"data"         )
          if (XML_Array_Length(thisItem,"data") /= recycledFractionTableMetallicityCount) call&
               & Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous' ,'metallicities array in XML file does not match&
               & internal expectation')
          call XML_Array_Read(thisItem,"data",tableTemporary)
          if (recycledFractionIndex(imfSelected) == 1) then
             recycledFractionTableMetallicity=tableTemporary
          else
             if (any(recycledFractionTableMetallicity /= tableTemporary)) call&
                  & Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous' ,'mismatch in metallicities array in XML file')
          end if
          deallocate(tableTemporary)

          ! Find the recycledFraction element and extract data.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"recycledFraction")
          dataList => getElementsByTagname(thisItem,"data"            )
          if (getLength(dataList) /= recycledFractionTableAgeCount*recycledFractionTableMetallicityCount) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous'&
               & ,'recycled fractions array in XML file does not match internal expectation')
          iRecycledFraction=0
          do iAge=1,recycledFractionTableAgeCount
             do iMetallicity=1,recycledFractionTableMetallicityCount
                iRecycledFraction=iRecycledFraction+1
                thisItem => item(dataList,iRecycledFraction-1)
                call extractDataContent(thisItem,recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected)))
             end do
          end do

          ! Destroy the document.
          call destroy(doc)
          call Galacticus_Display_Unindent('done',verbosityDebug)

       else

          call Galacticus_Display_Indent('Tabulating mass recycling rate for '//char(imfNames(imfSelected))//' IMF',verbosityWorking)
          call Galacticus_Display_Counter(0,.true.,verbosityWorking)

          ! Open an XML file to output the data to.
          call xml_OpenFile(char(fileName),recycledFractionDoc)
          call xml_NewElement(recycledFractionDoc,"stellarPopulation")
          call xml_NewElement(recycledFractionDoc,"fileFormat")
          call xml_AddCharacters(recycledFractionDoc,fileFormatCurrent)
          call xml_EndElement(recycledFractionDoc,"fileFormat")
          call xml_NewElement(recycledFractionDoc,"description")
          call xml_AddCharacters(recycledFractionDoc,"Recycled fraction for a "//char(imfNames(imfSelected))//" IMF")
          call xml_EndElement(recycledFractionDoc,"description")
          call xml_NewElement(recycledFractionDoc,"source")
          call xml_AddCharacters(recycledFractionDoc,"Computed by Galacticus")
          call xml_EndElement(recycledFractionDoc,"source")
          call xml_NewElement(recycledFractionDoc,"date")
          call xml_AddCharacters(recycledFractionDoc,char(Formatted_Date_and_Time()))
          call xml_EndElement(recycledFractionDoc,"date")

          ! Write ages to the XML file.
          call xml_NewElement(recycledFractionDoc,"ages")
          call xml_NewElement(recycledFractionDoc,"description")
          call xml_AddCharacters(recycledFractionDoc,"Age of the stellar population in Gyr")
          call xml_EndElement(recycledFractionDoc,"description")
          do iAge=1,recycledFractionTableAgeCount
             call xml_NewElement(recycledFractionDoc,"data")
             write (parameterValue,'(e10.4)') recycledFractionTableAge(iAge)
             call xml_AddCharacters(recycledFractionDoc,trim(parameterValue))
             call xml_EndElement(recycledFractionDoc,"data")
          end do
          call xml_EndElement(recycledFractionDoc,"ages")

          ! Write metallicities to the XML file.
          call xml_NewElement(recycledFractionDoc,"metallicities")
          call xml_NewElement(recycledFractionDoc,"description")
          call xml_AddCharacters(recycledFractionDoc,"Metallicity (fractional mass of total metals) of the stellar population")
          call xml_EndElement(recycledFractionDoc,"description")
          do iMetallicity=1,recycledFractionTableMetallicityCount
             call xml_NewElement(recycledFractionDoc,"data")
             write (parameterValue,'(e10.4)') recycledFractionTableMetallicity(iMetallicity)
             call xml_AddCharacters(recycledFractionDoc,trim(parameterValue))
             call xml_EndElement(recycledFractionDoc,"data")
          end do
          call xml_EndElement(recycledFractionDoc,"metallicities")

          ! Loop over ages and metallicities and compute the recycled fraction.
          imfSelectedGlobal=imfSelected
          call xml_NewElement(recycledFractionDoc,"recycledFraction")
          loopCountTotal=recycledFractionTableMetallicityCount*recycledFractionTableAgeCount
          loopCount     =0
          !$omp parallel do private (iAge,iMetallicity,progressMessage,minimumMass,maximumMass,integrandFunction,integrationWorkspace) copyin(imfSelectedGlobal)
          do iAge=1,recycledFractionTableAgeCount
             lifetime=recycledFractionTableAge(iAge)
             write (progressMessage,'(a6,e8.2,a4)') 'age = ',lifetime,' Gyr'
             call Galacticus_Display_Message(progressMessage,verbosityDebug)
             do iMetallicity=1,recycledFractionTableMetallicityCount
                ! Set the metallicity. If using the instantaneous recycling approximation, assume Solar metallicity always.
                if (starFormationImfInstantaneousApproximation) then
                   metallicity=metallicitySolar
                else
                   metallicity=recycledFractionTableMetallicity(iMetallicity)
                end if
                ! Update the counter.
                !$omp atomic
                loopCount=loopCount+1
                call Galacticus_Display_Counter(                                                   &
                     &                           int(100.0d0*dble(loopCount)/dble(loopCountTotal)) &
                     &                          ,.false.                                           &
                     &                          ,verbosityWorking                                  &
                     &                         )
                ! Find the minimum and maximum masses to integrate over for this IMF.
                minimumMass=IMF_Minimum_Mass(imfSelected)
                maximumMass=IMF_Maximum_Mass(imfSelected)
                ! Integrate ejected mass over the IMF between these limits.
                recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected))=Integrate(minimumMass,maximumMass&
                     &,Recycled_Fraction_Integrand ,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute&
                     &=1.0d-3 ,toleranceRelative=1.0d-4)
                call Integrate_Done(integrandFunction,integrationWorkspace)
             end do
          end do
          !$omp end parallel do
          do iAge=1,recycledFractionTableAgeCount
             do iMetallicity=1,recycledFractionTableMetallicityCount
                ! Enforce monotonicity in the recycled fraction. Non-monotonicity can arise due to the vagaries of interpolating
                ! stellar lifetimes in an irregular grid of stellar models.
                if (iAge > 1 )                                                                               &
                     & recycledFractionTable       (iAge  ,iMetallicity,recycledFractionIndex(imfSelected))  &
                     &   =max(                                                                               &
                     &        recycledFractionTable(iAge  ,iMetallicity,recycledFractionIndex(imfSelected)), &
                     &        recycledFractionTable(iAge-1,iMetallicity,recycledFractionIndex(imfSelected))  &
                     &       )
                call xml_NewElement(recycledFractionDoc,"data")
                write (parameterValue,'(e10.4)') recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected))
                call xml_AddCharacters(recycledFractionDoc,trim(parameterValue))
                call xml_EndElement(recycledFractionDoc,"data")
             end do
          end do
          call xml_EndElement(recycledFractionDoc,"recycledFraction")
          call Galacticus_Display_Counter_Clear(           verbosityWorking)
          call Galacticus_Display_Unindent     ('finished',verbosityWorking)
          call xml_EndElement(recycledFractionDoc,"stellarPopulation")
          call xml_Close(recycledFractionDoc)
       end if

       ! Flag that this IMF has now been tabulated.
       recycledFractionTabulated(imfSelected)=.true.
    end if
    !$omp end critical(IMF_Recycling_Rate_NonInstantaneous_Initialize)

    ! Get the index where this IMF is stored in the table.
    tableIndex=recycledFractionIndex(imfSelected)

    ! Interpolate to get the derivative in the recycled rate at two adjacent metallicities.
    metallicity=max(Abundances_Get_Metallicity(fuelAbundances),0.0d0)
    if (metallicity > recycledFractionTableMetallicityMaximum) then
       metallicityIndex=recycledFractionTableMetallicityCount
       metallicityFactors=[1.0d0,0.0d0]
       if (present(ageMaximum)) then
          ! Get average recycling rate between ageMinimum and ageMaximum.
          if (ageMinimum > 0.0d0) then
             recycledFractionMinimum=Interpolate(recycledFractionTableAgeCount&
                  &,recycledFractionTableAge ,recycledFractionTable(: ,metallicityIndex,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator ,ageMinimum,reset=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
          else
             recycledFractionMinimum=0.0d0
          end if
          recycledFractionMaximum=Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge,recycledFractionTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMaximum,reset &
               &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
          recycleRate(1)=(recycledFractionMaximum-recycledFractionMinimum)/(ageMaximum-ageMinimum)
       else
          ! Get instantaneous recycling rate at ageMinimum.
          recycleRate(1)=Interpolate_Derivative(recycledFractionTableAgeCount,recycledFractionTableAge,recycledFractionTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMinimum,reset&
               &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
       end if
       recycleRate(2)=0.0d0
    else
       metallicityIndex=Interpolate_Locate(recycledFractionTableMetallicityCount,recycledFractionTableMetallicity&
            &,interpolationMetallicityAccelerator,metallicity,reset=interpolationMetallicityReset)
       metallicityFactors=Interpolate_Linear_Generate_Factors(recycledFractionTableMetallicityCount,recycledFractionTableMetallicity&
            &,metallicityIndex,metallicity)
       ! Interpolate in age at both metallicities.
       do iMetallicity=0,1
          if (present(ageMaximum)) then
             ! Get average recycling rate between ageMinimum and ageMaximum.
             if (ageMinimum > 0.0d0) then
                recycledFractionMinimum=Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge ,recycledFractionTable(:,metallicityIndex&
                     &+iMetallicity,tableIndex),interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset&
                     &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
             else
                recycledFractionMinimum=0.0d0
             end if
             recycledFractionMaximum=Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge&
                  &,recycledFractionTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator,ageMaximum,reset =interpolationAgeReset,extrapolationType=extrapolationTypeLinear)

             recycleRate(iMetallicity+1)=(recycledFractionMaximum-recycledFractionMinimum)/(ageMaximum-ageMinimum)
          else
             ! Get instantaneous recycling rate at ageMinimum.
             recycleRate(iMetallicity+1)=Interpolate_Derivative(recycledFractionTableAgeCount,recycledFractionTableAge&
                  &,recycledFractionTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
          end if
       end do
    end if

    ! Interpolate in metallicity to get the actual rate.
    IMF_Recycling_Rate_NonInstantaneous=sum(metallicityFactors*recycleRate)

    return
  end function IMF_Recycling_Rate_NonInstantaneous

  logical function Star_Is_Evolved(initialMass,metallicity,age)
    !% Returns true if the specified star is evolved by the given {\tt age}.
    use Stellar_Astrophysics
    implicit none
    double precision, intent(in   ) :: age, initialMass, metallicity

    if (starFormationImfInstantaneousApproximation) then
       ! Instantaneous calculation - star is evolved if it is more massive that the specified mass of long-lived stars.
       Star_Is_Evolved=(initialMass > starFormationImfInstantaneousApproximationMassLongLived)
    else
       ! Standard calculation - star is evolved if its lifeltime is less than the supplied age.
       Star_Is_Evolved=(Star_Lifetime(initialMass,metallicity) < age)
    end if
    return
  end function Star_Is_Evolved

  function Recycled_Fraction_Integrand(initialMass,parameterPointer) bind(c)
    !% Integrand used in evaluating recycled fractions.
    use, intrinsic :: ISO_C_Binding
    use Stellar_Astrophysics
    implicit none
    real(kind=c_double)        :: Recycled_Fraction_Integrand
    real(kind=c_double), value :: initialMass
    type(c_ptr        ), value :: parameterPointer

    if (Star_Is_Evolved(initialMass,metallicity,lifetime)) then
       Recycled_Fraction_Integrand=IMF_Phi(initialMass,imfSelectedGlobal)*Star_Ejected_Mass(initialMass,metallicity)
    else
       Recycled_Fraction_Integrand=0.0d0
    end if
    return
  end function Recycled_Fraction_Integrand

  double precision function IMF_Metal_Yield_Rate_NonInstantaneous(starFormationRate,fuelAbundances,component,ageMinimum,ageMaximum,abundanceIndex)
    !% Returns the metal yield rate for a simple stellar population, either for the total metallicity or, if {\tt atomIndex} is
    !% given, for the specified element. The \gls{imf} is determined from the given {\tt starFormationRate} and {\tt fuelAbundances}.
    !% The metal yield rate (in fraction of the population's mass returned to the \gls{ism} as new metals per Gyr) is computed for the
    !% given {\tt age} (in Gyr). The metal yield is computed on a grid of age and metallicity. This is stored to file and will be
    !% read back in on subsequent runs. This is useful as computation of the table is relatively slow.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Ranges
    use Numerical_Constants_Astronomical
    use Memory_Management
    use Galacticus_Display
    use File_Utilities
    use FoX_wxml
    use FoX_dom
    use IO_XML
    use Galacticus_Error
    use Dates_and_Times
    use Galacticus_Input_Paths
     implicit none
     double precision                                                                            , intent(in   )           :: ageMinimum                         , starFormationRate
     double precision                                                                            , intent(in   ), optional :: ageMaximum
     integer                                                                                     , intent(in   ), optional :: abundanceIndex
     type            (abundances                )                                                , intent(in   )           :: fuelAbundances
     integer                                                                                     , intent(in   )           :: component
     logical                                                    , allocatable, dimension(:      )                          :: metalYieldTabulatedTemporary
     integer                                                    , allocatable, dimension(:      )                          :: metalYieldIndexTemporary
     double precision                                           , allocatable, dimension(:    )                            :: tableTemporary
     double precision                                           , allocatable, dimension(:,:,:,:)                          :: metalYieldTableTemporary
     double precision                                                        , dimension(2      )                          :: metalYieldRate                     , metallicityFactors
     type            (fgsl_interp               )         , save                                                           :: interpolationAgeObject
     type            (fgsl_interp_accel         )         , save                                                           :: interpolationAgeAccelerator        , interpolationMetallicityAccelerator
     logical                                              , save                                                           :: interpolationAgeReset       =.true., interpolationMetallicityReset      =.true.
     !$omp threadprivate(interpolationAgeObject,interpolationMetallicityAccelerator &
     !$omp ,interpolationAgeAccelerator,interpolationMetallicityReset,interpolationAgeReset)
     type            (Node                      ), pointer                                                                 :: doc                                , thisItem
     type            (NodeList                  ), pointer                                                                 :: dataList
     type            (c_ptr                     )                                                                          :: parameterPointer
     type            (fgsl_function             )                                                                          :: integrandFunction
     type            (fgsl_integration_workspace)                                                                          :: integrationWorkspace
     integer                                                                                                               :: abundanceIndexActual               , fileFormat                                , &
          &                                                                                                                   iAge                               , iElement                                  , &
          &                                                                                                                   iMetalYield                        , iMetallicity                              , &
          &                                                                                                                   imfCount                           , imfSelected                               , &
          &                                                                                                                   ioErr                              , loopCount                                 , &
          &                                                                                                                   loopCountTotal                     , metallicityIndex                          , &
          &                                                                                                                   tableIndex
     double precision                                                                                                      :: maximumMass                        , minimumMass                               , &
          &                                                                                                                   yieldMaximum                       , yieldMinimum
     character       (len=20                    )                                                                          :: parameterValue                     , progressMessage
     type            (xmlf_t                    )                                                                          :: metalYieldDoc
     type            (varying_string            )                                                                          :: fileName
     logical                                                                                                               :: makeFile

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    !$omp critical(IMF_Metal_Yield_Rate_NonInstantaneous_Initialize)
    ! Check that flag and index arrays exist.
    if (.not.allocated(metalYieldTabulated)) then
       call Alloc_Array(metalYieldTabulated,[imfSelected],file=__FILE__,line=__LINE__)
       call Alloc_Array(metalYieldIndex    ,[imfSelected],file=__FILE__,line=__LINE__)
       metalYieldTabulated=.false.
       metalYieldIndex    =0
    end if

    ! Check that flag and index arrays are large enough.
    if (size(metalYieldTabulated) < imfSelected) then
       call Move_Alloc (metalYieldTabulated,metalYieldTabulatedTemporary)
       call Move_Alloc (metalYieldIndex    ,metalYieldIndexTemporary    )
       call Alloc_Array(metalYieldTabulated,[imfSelected],file=__FILE__,line=__LINE__)
       call Alloc_Array(metalYieldIndex    ,[imfSelected],file=__FILE__,line=__LINE__)
       metalYieldTabulated(1:size(metalYieldTabulatedTemporary))            =metalYieldTabulatedTemporary
       metalYieldIndex    (1:size(metalYieldIndexTemporary    ))            =metalYieldIndexTemporary
       metalYieldTabulated(  size(metalYieldTabulatedTemporary):imfSelected)=.false.
       metalYieldIndex    (  size(metalYieldTabulatedTemporary):imfSelected)=0
       call Dealloc_Array(metalYieldTabulatedTemporary,file=__FILE__,line=__LINE__)
       call Dealloc_Array(metalYieldIndexTemporary    ,file=__FILE__,line=__LINE__)
    end if

    ! Tabulate the metal yield for this IMF if it has not already been computed.
    if (.not.metalYieldTabulated(imfSelected)) then

       ! Expand the tabulations array by enough to accomodate a new IMF.
       if (allocated(metalYieldTable)) then
          imfCount=size(metalYieldTable,dim=4)
          call Move_Alloc(metalYieldTable,metalYieldTableTemporary)
          call Alloc_Array(metalYieldTable,[metalYieldTableAgeCount,metalYieldTableMetallicityCount,elementCount+1,imfCount],file=__FILE__,line=__LINE__)
          metalYieldTable(:,:,:,1:imfCount)=metalYieldTableTemporary
          call Dealloc_Array(metalYieldTableTemporary)
       else
          call Alloc_Array(metalYieldTableAge        ,[metalYieldTableAgeCount        ],file=__FILE__,line=__LINE__)
          call Alloc_Array(metalYieldTableMetallicity,[metalYieldTableMetallicityCount],file=__FILE__,line=__LINE__)
          metalYieldTableAge                                           =Make_Range(metalYieldTableAgeMinimum&
               &,metalYieldTableAgeMaximum,metalYieldTableAgeCount,rangeType=rangeTypeLogarithmic)
          metalYieldTableMetallicity(1)                                =0.0d0
          metalYieldTableMetallicity(2:metalYieldTableMetallicityCount)&
               &=Make_Range(metalYieldTableMetallicityMinimum,metalYieldTableMetallicityMaximum&
               &,metalYieldTableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          call Alloc_Array(metalYieldTable,[metalYieldTableAgeCount,metalYieldTableMetallicityCount,elementCount+1,1],file=__FILE__,line=__LINE__)
       end if

       ! Record the index in the array where this IMF will be stored.
       metalYieldIndex(imfSelected)=size(metalYieldTable,dim=4)

       ! Loop through all elements (and total metallicity) for which yield must be tabulated.
       elementsLoop : do iElement=1,elementCount ! iElement=1 will correspond to total metallicity.

          ! Check if the table has been computed and stored previously.
          select case (iElement)
          case (1)
             ! Total metallicity.
             fileName=char(Galacticus_Input_Path())//'data/stellarPopulations/Stellar_Metal_Yield_'//imfNames(imfSelected)//'_'//imfUniqueLabel//'.xml'
          case (2:)
             ! Individual element.
             fileName=char(Galacticus_Input_Path())//'data/stellarPopulations/Stellar_'//Abundances_Names(iElement)//'_Yield_'//imfNames(imfSelected)//'_'//imfUniqueLabel//'.xml'
          end select
          makeFile=.false.
          if (File_Exists(fileName)) then
             ! Open the XML file containing energy input.
             call Galacticus_Display_Indent('Parsing file: '//fileName,verbosityDebug)
             doc => parseFile(char(fileName),iostat=ioErr)
             if (ioErr /= 0) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','Unable to parse metal yields file "'//fileName//'"'//char(10)//'HELP: file may be corrupted - try deleting this file and rerunning Galacticus')
             ! Check the version number.
             if (XML_Path_Exists(doc,"fileFormat")) then
                thisItem => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
                call extractDataContent(thisItem,fileFormat)
                if (fileFormat /= fileFormatCurrent) makeFile=.true.
             else
                makeFile=.true.
             end if
             if (makeFile) then
                call destroy(doc)
                call Galacticus_Display_Unindent('done',verbosityDebug)
             end if
          else
             makeFile=.true.
          end if

          fileBuildCheck : if (.not.makeFile) then

             ! Find the ages element and extract data.
             thisItem => XML_Get_First_Element_By_Tag_Name(doc,"ages")
             if (XML_Array_Length(thisItem,"data") /= metalYieldTableAgeCount) call&
                  & Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous' ,'ages array in XML file does not match&
                  & internal expectation')
             call XML_Array_Read(thisItem,"data",tableTemporary)
             if (iElement == 1 .and. metalYieldIndex(imfSelected) == 1) then
                metalYieldTableAge=tableTemporary
             else
                if (any(metalYieldTableAge /= tableTemporary)) call&
                     & Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous' ,'mismatch in ages array in XML file')
             end if
             deallocate(tableTemporary)

             ! Find the metallicities element and extract data.
             thisItem => XML_Get_First_Element_By_Tag_Name(doc,"metallicities")
             if (XML_Array_Length(thisItem,"data") /= metalYieldTableMetallicityCount) call&
                  & Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous' ,'metallicities array in XML file does not&
                  & match internal expectation')
             call XML_Array_Read(thisItem,"data",tableTemporary)
             if (iElement == 1 .and. metalYieldIndex(imfSelected) == 1) then
                metalYieldTableMetallicity=tableTemporary
             else
                if (any(metalYieldTableMetallicity /= tableTemporary)) call&
                     & Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous' ,'mismatch in metallicities array in XML&
                     & file')
             end if
             deallocate(tableTemporary)

             ! Find the metalYield element and extract data.
             select case (iElement)
             case(1)
                ! Total metallicity.
                thisItem => XML_Get_First_Element_By_Tag_Name(doc,"metalYield"  )
             case (2:)
                ! Individual element.
                thisItem => XML_Get_First_Element_By_Tag_Name(doc,"elementYield")
             end select
             dataList => getElementsByTagname(thisItem,"data")
             if (getLength(dataList) /= metalYieldTableAgeCount*metalYieldTableMetallicityCount) call&
                  & Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous' ,'metal yield array in XML file does not&
                  & match internal expectation')
             iMetalYield=0
             do iAge=1,metalYieldTableAgeCount
                do iMetallicity=1,metalYieldTableMetallicityCount
                   iMetalYield=iMetalYield+1
                   thisItem => item(dataList,iMetalYield-1)
                   call extractDataContent(thisItem,metalYieldTable(iAge,iMetallicity,iElement &
                        &,metalYieldIndex(imfSelected)))
                end do
             end do

             ! Destroy the document.
             call destroy(doc)
             call Galacticus_Display_Unindent('done',verbosityDebug)

          else

             ! Display a message since this calculation will take a long time.
             select case (iElement)
             case (1)
                call Galacticus_Display_Indent('Tabulating metal yield rate for '//char(imfNames(imfSelected))//' IMF',verbosityWorking)
             case (2:)
                call Galacticus_Display_Indent('Tabulating '//char(Abundances_Names(iElement))//' yield rate for '&
                     &//char(imfNames(imfSelected))//' IMF',verbosityWorking)
             end select
             call Galacticus_Display_Counter(0,.true.,verbosityWorking)

             ! Open an XML file to output the data to.
             call xml_OpenFile(char(fileName),metalYieldDoc)
             call xml_NewElement(metalYieldDoc,"stellarPopulation")
             call xml_NewElement(metalYieldDoc,"fileFormat")
             call xml_AddCharacters(metalYieldDoc,fileFormatCurrent)
             call xml_EndElement(metalYieldDoc,"fileFormat")
             call xml_NewElement(metalYieldDoc,"description")
             select case (iElement)
             case (1)
                call xml_AddCharacters(metalYieldDoc,"Metal yield for a "//char(imfNames(imfSelected))//" IMF")
             case (2:)
                call xml_AddCharacters(metalYieldDoc,char(Abundances_Names(iElement))//" yield for a "//char(imfNames(imfSelected))//" IMF")
             end select
             call xml_EndElement(metalYieldDoc,"description")
             call xml_NewElement(metalYieldDoc,"source")
             call xml_AddCharacters(metalYieldDoc,"Computed by Galacticus")
             call xml_EndElement(metalYieldDoc,"source")
             call xml_NewElement(metalYieldDoc,"date")
             call xml_AddCharacters(metalYieldDoc,char(Formatted_Date_and_Time()))
             call xml_EndElement(metalYieldDoc,"date")

             ! Write ages to the XML file.
             call xml_NewElement(metalYieldDoc,"ages")
             call xml_NewElement(metalYieldDoc,"description")
             call xml_AddCharacters(metalYieldDoc,"Age of the stellar population in Gyr")
             call xml_EndElement(metalYieldDoc,"description")
             do iAge=1,metalYieldTableAgeCount
                call xml_NewElement(metalYieldDoc,"data")
                write (parameterValue,'(e10.4)') metalYieldTableAge(iAge)
                call xml_AddCharacters(metalYieldDoc,trim(parameterValue))
                call xml_EndElement(metalYieldDoc,"data")
             end do
             call xml_EndElement(metalYieldDoc,"ages")

             ! Write metallicities to the XML file.
             call xml_NewElement(metalYieldDoc,"metallicities")
             call xml_NewElement(metalYieldDoc,"description")
             call xml_AddCharacters(metalYieldDoc,"Metallicity (fractional mass of total metals) of the stellar population")
             call xml_EndElement(metalYieldDoc,"description")
             do iMetallicity=1,metalYieldTableMetallicityCount
                call xml_NewElement(metalYieldDoc,"data")
                write (parameterValue,'(e10.4)') metalYieldTableMetallicity(iMetallicity)
                call xml_AddCharacters(metalYieldDoc,trim(parameterValue))
                call xml_EndElement(metalYieldDoc,"data")
             end do
             call xml_EndElement(metalYieldDoc,"metallicities")

             ! Loop over ages and metallicities and compute the metal yield.
             imfSelectedGlobal=imfSelected
             atomIndexGlobal=Abundances_Atomic_Index(iElement)
             select case (iElement)
             case (1)
                ! Total metallicity.
                call xml_NewElement(metalYieldDoc,"metalYield")
             case default
                ! Individual element.
                call xml_NewElement(metalYieldDoc,"elementYield")
             end select
             loopCountTotal=metalYieldTableMetallicityCount*metalYieldTableAgeCount
             loopCount     =0
             !$omp parallel do private (iAge,iMetallicity,progressMessage,minimumMass,maximumMass,integrandFunction,integrationWorkspace) copyin(imfSelectedGlobal,atomIndexGlobal)
             do iAge=1,metalYieldTableAgeCount
                lifetime=metalYieldTableAge(iAge)
                write (progressMessage,'(a6,e8.2,a4)') 'age = ',lifetime,' Gyr'
                call Galacticus_Display_Message(progressMessage,verbosityDebug)
                do iMetallicity=1,metalYieldTableMetallicityCount
                   ! Set the metallicity. If using the instantaneous recycling approximation, assume Solar metallicity always.
                   if (starFormationImfInstantaneousApproximation) then
                      metallicity=metallicitySolar
                   else
                      metallicity=metalYieldTableMetallicity(iMetallicity)
                   end if
                   ! Update the counter.
                   !$omp atomic
                   loopCount=loopCount+1
                   call Galacticus_Display_Counter(                                                   &
                        &                           int(100.0d0*dble(loopCount)/dble(loopCountTotal)) &
                        &                          ,.false.                                           &
                        &                          ,verbosityWorking                                  &
                        &                         )
                   ! Find the minimum and maximum masses to integrate over for this IMF.
                   minimumMass=IMF_Minimum_Mass(imfSelected)
                   maximumMass=IMF_Maximum_Mass(imfSelected)
                   ! Integrate ejected mass over the IMF between these limits.
                   metalYieldTable(iAge,iMetallicity,iElement,metalYieldIndex(imfSelected))=Integrate(minimumMass,maximumMass&
                        &,Metal_Yield_Integrand,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=1.0d-4&
                        &,toleranceRelative=1.0d-5)
                   call Integrate_Done(integrandFunction,integrationWorkspace)
                end do
             end do
             !$omp end parallel do
             do iAge=1,metalYieldTableAgeCount
                do iMetallicity=1,metalYieldTableMetallicityCount
                   ! Enforce monotonicity in the metal yield. Non-monotonicity can arise due to the vagaries of interpolating
                   ! stellar lifetimes in an irregular grid of stellar models.
                   if (iAge > 1) metalYieldTable(iAge,iMetallicity,iElement,metalYieldIndex(imfSelected))&
                        &=max(metalYieldTable(iAge,iMetallicity,iElement,metalYieldIndex(imfSelected)),metalYieldTable(iAge-1&
                        &,iMetallicity,iElement,metalYieldIndex(imfSelected)))
                   call xml_NewElement(metalYieldDoc,"data")
                   write (parameterValue,'(e10.4)') metalYieldTable(iAge,iMetallicity,iElement,metalYieldIndex(imfSelected))
                   call xml_AddCharacters(metalYieldDoc,trim(parameterValue))
                   call xml_EndElement(metalYieldDoc,"data")
                end do
             end do
             select case (iElement)
             case (1)
                ! Total metallicity.
                call xml_EndElement(metalYieldDoc,"metalYield")
             case default
                ! Individual element.
                call xml_EndElement(metalYieldDoc,"elementYield")
             end select
             call Galacticus_Display_Counter_Clear(           verbosityWorking)
             call Galacticus_Display_Unindent     ('finished',verbosityWorking)
             call xml_EndElement(metalYieldDoc,"stellarPopulation")
             call xml_Close(metalYieldDoc)
          end if fileBuildCheck

       end do elementsLoop

       ! Flag that this IMF has now been tabulated.
       metalYieldTabulated(imfSelected)=.true.
    end if
    !$omp end critical(IMF_Metal_Yield_Rate_NonInstantaneous_Initialize)

    ! Get the index where this IMF is stored in the table.
    tableIndex=metalYieldIndex(imfSelected)

    ! Determine which element (or total metals) we've been asked for.
    if (present(abundanceIndex)) then
       ! Use the given abundance index.
       abundanceIndexActual=abundanceIndex
    else
       ! Assume that total metallicity is required.
       abundanceIndexActual=1
    end if

     ! Interpolate to get the derivative in the metal yield at two adjacent metallicities.
     metallicity=max(Abundances_Get_Metallicity(fuelAbundances),0.0d0)
     if (metallicity > metalYieldTableMetallicityMaximum) then
        metallicityIndex=metalYieldTableMetallicityCount
        metallicityFactors=[1.0d0,0.0d0]
        if (present(ageMaximum)) then
           ! Get average recycling rate between ageMinimum and ageMaximum.
           if (ageMinimum > 0.0d0) then
              yieldMinimum=Interpolate(metalYieldTableAgeCount,metalYieldTableAge &
                   &,metalYieldTable(: ,metallicityIndex,abundanceIndexActual,tableIndex) ,interpolationAgeObject&
                   &,interpolationAgeAccelerator ,ageMinimum,reset=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
           else
              yieldMinimum=0.0d0
           end if
           yieldMaximum=Interpolate(metalYieldTableAgeCount,metalYieldTableAge,metalYieldTable(: ,metallicityIndex &
                &,abundanceIndexActual,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMaximum,reset &
                &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
           metalYieldRate(1)=(yieldMaximum-yieldMinimum)/(ageMaximum-ageMinimum)
        else
           ! Get instantaneous recycling rate at ageMinimum.
           metalYieldRate(1)=Interpolate_Derivative(metalYieldTableAgeCount,metalYieldTableAge,metalYieldTable(: &
                &,metallicityIndex,abundanceIndexActual,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMinimum&
                & ,reset=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
        end if
        metalYieldRate(2)=0.0d0
     else
        metallicityIndex=Interpolate_Locate(metalYieldTableMetallicityCount,metalYieldTableMetallicity &
             &,interpolationMetallicityAccelerator,metallicity,reset=interpolationMetallicityReset)
        metallicityFactors=Interpolate_Linear_Generate_Factors(metalYieldTableMetallicityCount,metalYieldTableMetallicity&
             &,metallicityIndex,metallicity)
        ! Interpolate in age at both metallicities.
        do iMetallicity=0,1
          if (present(ageMaximum)) then
             ! Get average recycling rate between ageMinimum and ageMaximum.
             if (ageMinimum > 0.0d0) then
                yieldMinimum=Interpolate(metalYieldTableAgeCount ,metalYieldTableAge,metalYieldTable(:,metallicityIndex+iMetallicity&
                     &,abundanceIndexActual,tableIndex) ,interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset &
                     &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
             else
                yieldMinimum=0.0d0
             end if
             yieldMaximum=Interpolate(metalYieldTableAgeCount,metalYieldTableAge ,metalYieldTable(: &
                  &,metallicityIndex+iMetallicity,abundanceIndexActual,tableIndex),interpolationAgeObject &
                  &,interpolationAgeAccelerator,ageMaximum,reset =interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
             metalYieldRate(iMetallicity+1)=(yieldMaximum-yieldMinimum)/(ageMaximum-ageMinimum)
          else
             ! Get instantaneous recycling rate at ageMinimum.
             metalYieldRate(iMetallicity+1)=Interpolate_Derivative(metalYieldTableAgeCount,metalYieldTableAge ,metalYieldTable(: &
                  &,metallicityIndex+iMetallicity,abundanceIndexActual,tableIndex),interpolationAgeObject &
                  &,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
          end if
       end do
     end if

     ! Interpolate in metallicity to get the actual rate.
     IMF_Metal_Yield_Rate_NonInstantaneous=sum(metallicityFactors*metalYieldRate)

     return
   end function IMF_Metal_Yield_Rate_NonInstantaneous

  function Metal_Yield_Integrand(initialMass,parameterPointer) bind(c)
    !% Integrand used in evaluating metal yields.
    use, intrinsic :: ISO_C_Binding
    use Stellar_Astrophysics
    use Supernovae_Type_Ia
    implicit none
    real(kind=c_double)        :: Metal_Yield_Integrand
    real(kind=c_double), value :: initialMass
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: sneiaLifetime        , yieldMass

    ! Include yields from isolated stars.
    if (Star_Is_Evolved(initialMass,metallicity,lifetime)) then
       select case (atomIndexGlobal)
       case (0)
          ! Total metallicity required.
          yieldMass=Star_Metal_Yield_Mass(initialMass,metallicity                )
       case default
          ! Inidividual element required.
          yieldMass=Star_Metal_Yield_Mass(initialMass,metallicity,atomIndexGlobal)
       end select
       Metal_Yield_Integrand=IMF_Phi(initialMass,imfSelectedGlobal)*yieldMass
    else
       Metal_Yield_Integrand=0.0d0
    end if

    ! Include yield from Type Ia supernovae.
    if (starFormationImfInstantaneousApproximation) then
       ! In the instantaneous stellar evolution approximation use the effective age to compute the SNeIa yield.
       sneiaLifetime=starFormationImfInstantaneousApproximationEffectiveAge
    else
       ! In the standard calculation simply use the current age.
       sneiaLifetime=lifetime
    end if
    select case (atomIndexGlobal)
    case (0)
       ! Total metallicity required.
       yieldMass=SNeIa_Cumulative_Yield(initialMass,sneiaLifetime,metallicity                )
    case default
       yieldMass=SNeIa_Cumulative_Yield(initialMass,sneiaLifetime,metallicity,atomIndexGlobal)
    end select
    Metal_Yield_Integrand=Metal_Yield_Integrand+IMF_Phi(initialMass,imfSelectedGlobal)*yieldMass
    return
  end function Metal_Yield_Integrand

  double precision function IMF_Energy_Input_Rate_NonInstantaneous(starFormationRate,fuelAbundances,component,ageMinimum,ageMaximum)
    !% Returns the energy input rate for a simple stellar population in (km/s)$^2$ Gyr$^{-1}$. The \gls{imf} is determined from the
    !% given {\tt starFormationRate} and {\tt fuelAbundances}. The energy input rate is computed for the given {\tt age} (in
    !% Gyr). The cumulative energy input is computed on a grid of age and metallicity. This is stored to file and will be read
    !% back in on subsequent runs. This is useful as computation of the table is relatively slow.
    use, intrinsic :: ISO_C_Binding
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Ranges
    use Numerical_Constants_Astronomical
    use Memory_Management
    use Galacticus_Display
    use File_Utilities
    use FoX_wxml
    use FoX_dom
    use IO_XML
    use Galacticus_Error
    use Dates_and_Times
    use Galacticus_Input_Paths
    implicit none
    double precision                                                                          , intent(in   ) ::        ageMinimum                                    , starFormationRate
    double precision                                                                          , intent(in   ) , optional::                           ageMaximum
    type            (abundances                )                                              , intent(in   ) ::        fuelAbundances
    integer                                                                                   , intent(in   ) ::        component
    logical                                                    , allocatable, dimension(:    )                ::        energyInputTabulatedTemporary
    integer                                                    , allocatable, dimension(:    )                ::        energyInputIndexTemporary
    double precision                                           , allocatable, dimension(:    )                ::        tableTemporary
    double precision                                           , allocatable, dimension(:,:,:)                ::        energyInputTableTemporary
    double precision                                                        , dimension(2    )                ::        energyInputRate                               , metallicityFactors
    type            (fgsl_interp               )         , save                                               ::        interpolationAgeObject
    type            (fgsl_interp_accel         )         , save                                               ::        interpolationAgeAccelerator                   , interpolationMetallicityAccelerator
    logical                                              , save                                               ::        interpolationAgeReset                  =.true., interpolationMetallicityReset      =.true.
    !$omp threadprivate(interpolationAgeObject,interpolationMetallicityAccelerator &
    !$omp ,interpolationAgeAccelerator,interpolationMetallicityReset,interpolationAgeReset)
    type            (Node                      ), pointer                                                     ::        doc                                           , thisItem
    type            (NodeList                  ), pointer                                                     ::        dataList
    type            (c_ptr                     )                                                              ::        parameterPointer
    type            (fgsl_function             )                                                              ::        integrandFunction
    type            (fgsl_integration_workspace)                                                              ::        integrationWorkspace
    integer                                                                                                   ::        fileFormat                                    , iAge                                      , &
         &                                                                                                              iEnergyInput                                  , iMetallicity                              , &
         &                                                                                                              imfCount                                      , imfSelected                               , &
         &                                                                                                              ioErr                                         , loopCount                                 , &
         &                                                                                                              loopCountTotal                                , metallicityIndex                          , &
         &                                                                                                              tableIndex
    double precision                                                                                          ::        energyInputMaximum                            , energyInputMinimum                        , &
         &                                                                                                              maximumMass                                   , minimumMass
    character       (len=20                    )                                                              ::        parameterValue                                , progressMessage
    type            (xmlf_t                    )                                                              ::        energyInputDoc
    type            (varying_string            )                                                              ::        fileName
    logical                                                                                                   ::        makeFile

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances,component)

    !$omp critical(IMF_Energy_Input_Rate_NonInstantaneous_Initialize)
    ! Check that flag and index arrays exist.
    if (.not.allocated(energyInputTabulated)) then
       call Alloc_Array(energyInputTabulated,[imfSelected],file=__FILE__,line=__LINE__)
       call Alloc_Array(energyInputIndex    ,[imfSelected],file=__FILE__,line=__LINE__)
       energyInputTabulated=.false.
       energyInputIndex    =0
    end if

    ! Check that flag and index arrays are large enough.
    if (size(energyInputTabulated) < imfSelected) then
       call Move_Alloc (energyInputTabulated,energyInputTabulatedTemporary)
       call Move_Alloc (energyInputIndex    ,energyInputIndexTemporary    )
       call Alloc_Array(energyInputTabulated,[imfSelected],file=__FILE__,line=__LINE__)
       call Alloc_Array(energyInputIndex    ,[imfSelected],file=__FILE__,line=__LINE__)
       energyInputTabulated(1:size(energyInputTabulatedTemporary))=energyInputTabulatedTemporary
       energyInputIndex    (1:size(energyInputIndexTemporary    ))=energyInputIndexTemporary
       energyInputTabulated(  size(energyInputTabulatedTemporary):imfSelected)=.false.
       energyInputIndex    (  size(energyInputTabulatedTemporary):imfSelected)=0
       call Dealloc_Array(energyInputTabulatedTemporary,file=__FILE__,line=__LINE__)
       call Dealloc_Array(energyInputIndexTemporary    ,file=__FILE__,line=__LINE__)
    end if

    ! Tabulate the cumulative energy input for this IMF if it has not already been computed.
    if (.not.energyInputTabulated(imfSelected)) then

       ! Expand the tabulations array by enough to accomodate a new IMF.
       if (allocated(energyInputTable)) then
          imfCount=size(energyInputTable,dim=3)
          call Move_Alloc(energyInputTable,energyInputTableTemporary)
          call Alloc_Array(energyInputTable,[energyInputTableAgeCount,energyInputTableMetallicityCount,imfCount],file=__FILE__,line=__LINE__)
          energyInputTable(:,:,1:imfCount)=energyInputTableTemporary
          call Dealloc_Array(energyInputTableTemporary,file=__FILE__,line=__LINE__)
       else
          call Alloc_Array(energyInputTableAge        ,[energyInputTableAgeCount        ],file=__FILE__,line=__LINE__)
          call Alloc_Array(energyInputTableMetallicity,[energyInputTableMetallicityCount],file=__FILE__,line=__LINE__)
          energyInputTableAge=Make_Range(energyInputTableAgeMinimum&
               &,energyInputTableAgeMaximum,energyInputTableAgeCount,rangeType=rangeTypeLogarithmic)
          energyInputTableMetallicity(1)=0.0d0
          energyInputTableMetallicity(2:energyInputTableMetallicityCount)&
               &=Make_Range(energyInputTableMetallicityMinimum,energyInputTableMetallicityMaximum&
               &,energyInputTableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          call Alloc_Array(energyInputTable,[energyInputTableAgeCount,energyInputTableMetallicityCount,1],file=__FILE__,line=__LINE__)
       end if

       ! Record the index in the array where this IMF will be stored.
       energyInputIndex(imfSelected)=size(energyInputTable,dim=3)

       ! Check if the table has been computed and stored previously.
       fileName=char(Galacticus_Input_Path())//'data/stellarPopulations/Stellar_Energy_Input_'//imfNames(imfSelected)//'_'//imfUniqueLabel//'.xml'
       makeFile=.false.
       if (File_Exists(fileName)) then
           ! Open the XML file containing energy input.
          call Galacticus_Display_Indent('Parsing file: '//fileName,verbosityDebug)
          doc => parseFile(char(fileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','Unable to parse energy input filee "'//fileName//'"'//char(10)//'HELP: file may be corrupted - try deleting this file and rerunning Galacticus')
          ! Check the version number.
          if (XML_Path_Exists(doc,"fileFormat")) then
             thisItem => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
             call extractDataContent(thisItem,fileFormat)
             if (fileFormat /= fileFormatCurrent) makeFile=.true.
          else
             makeFile=.true.
          end if
          if (makeFile) then
             call destroy(doc)
             call Galacticus_Display_Unindent('done',verbosityDebug)
          end if
       else
          makeFile=.true.
       end if

       if (.not.makeFile) then

          ! Find the ages element and extract data.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"ages")
          if (XML_Array_Length(thisItem,"data") /= energyInputTableAgeCount) call&
               & Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','ages array in XML file does not match internal&
               & expectation')
          call XML_Array_Read(thisItem,"data",tableTemporary)
          if (energyInputIndex(imfSelected) == 1) then
             energyInputTableAge=tableTemporary
          else
             if (any(energyInputTableAge /= tableTemporary)) call&
                  & Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','mismatch in ages array in XML file')
          end if
          deallocate(tableTemporary)

          ! Find the metallicities element and extract data.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"metallicities")
          if (XML_Array_Length(thisItem,"data") /= energyInputTableMetallicityCount) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','metallicities array in XML file does not match internal expectation')
          call XML_Array_Read(thisItem,"data",tableTemporary)
          if (energyInputIndex(imfSelected) == 1) then
             energyInputTableMetallicity=tableTemporary
          else
             if (any(energyInputTableMetallicity /= tableTemporary)) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','mismatch in metallicities array in XML file')
          end if
          deallocate(tableTemporary)

          ! Find the energyInput element and extract data.
          thisItem => XML_Get_First_Element_By_Tag_Name(doc,"energyInput")
          dataList => getElementsByTagname(thisItem,"data"       )
          if (getLength(dataList) /= energyInputTableAgeCount*energyInputTableMetallicityCount) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','energy input array in XML file does not match internal expectation')
          iEnergyInput=0
          do iAge=1,energyInputTableAgeCount
             do iMetallicity=1,energyInputTableMetallicityCount
                iEnergyInput=iEnergyInput+1
                thisItem => item(dataList,iEnergyInput-1)
                call extractDataContent(thisItem,energyInputTable(iAge,iMetallicity,energyInputIndex(imfSelected)))
             end do
          end do

          ! Destroy the document.
          call destroy(doc)
          call Galacticus_Display_Unindent('done',verbosityDebug)

       else

          call Galacticus_Display_Indent('Tabulating cumulative energy input for '//char(imfNames(imfSelected))//' IMF',verbosityWorking)
          call Galacticus_Display_Counter(0,.true.,verbosityWorking)

          ! Open an XML file to output the data to.
          call xml_OpenFile(char(fileName),energyInputDoc)
          call xml_NewElement(energyInputDoc,"stellarPopulation")
          call xml_NewElement(energyInputDoc,"fileFormat")
          call xml_AddCharacters(energyInputDoc,fileFormatCurrent)
          call xml_EndElement(energyInputDoc,"fileFormat")
          call xml_NewElement(energyInputDoc,"description")
          call xml_AddCharacters(energyInputDoc,"Cumulative energy input for a "//char(imfNames(imfSelected))//" IMF")
          call xml_EndElement(energyInputDoc,"description")
          call xml_NewElement(energyInputDoc,"source")
          call xml_AddCharacters(energyInputDoc,"Computed by Galacticus")
          call xml_EndElement(energyInputDoc,"source")
          call xml_NewElement(energyInputDoc,"date")
          call xml_AddCharacters(energyInputDoc,char(Formatted_Date_and_Time()))
          call xml_EndElement(energyInputDoc,"date")

          ! Write ages to the XML file.
          call xml_NewElement(energyInputDoc,"ages")
          call xml_NewElement(energyInputDoc,"description")
          call xml_AddCharacters(energyInputDoc,"Age of the stellar population in Gyr")
          call xml_EndElement(energyInputDoc,"description")
          do iAge=1,energyInputTableAgeCount
             call xml_NewElement(energyInputDoc,"data")
             write (parameterValue,'(e10.4)') energyInputTableAge(iAge)
             call xml_AddCharacters(energyInputDoc,trim(parameterValue))
             call xml_EndElement(energyInputDoc,"data")
          end do
          call xml_EndElement(energyInputDoc,"ages")

          ! Write metallicities to the XML file.
          call xml_NewElement(energyInputDoc,"metallicities")
          call xml_NewElement(energyInputDoc,"description")
          call xml_AddCharacters(energyInputDoc,"Metallicity (fractional mass of total metals) of the stellar population")
          call xml_EndElement(energyInputDoc,"description")
          do iMetallicity=1,energyInputTableMetallicityCount
             call xml_NewElement(energyInputDoc,"data")
             write (parameterValue,'(e10.4)') energyInputTableMetallicity(iMetallicity)
             call xml_AddCharacters(energyInputDoc,trim(parameterValue))
             call xml_EndElement(energyInputDoc,"data")
          end do
          call xml_EndElement(energyInputDoc,"metallicities")
          ! Loop over ages and metallicities and compute the recycled fraction.
          imfSelectedGlobal=imfSelected
          call xml_NewElement(energyInputDoc,"energyInput")
          loopCountTotal=energyInputTableMetallicityCount*energyInputTableAgeCount
          loopCount     =0
          !$omp parallel do private (iAge,iMetallicity,progressMessage,minimumMass,maximumMass,integrandFunction,integrationWorkspace) copyin (imfSelectedGlobal)
          do iAge=1,energyInputTableAgeCount
             lifetime=energyInputTableAge(iAge)
             write (progressMessage,'(a6,e8.2,a4)') 'age = ',lifetime,' Gyr'
             call Galacticus_Display_Message(progressMessage,verbosityDebug)
             do iMetallicity=1,energyInputTableMetallicityCount
                ! Set the metallicity. If using the instantaneous recycling approximation, assume Solar metallicity always.
                if (starFormationImfInstantaneousApproximation) then
                   metallicity=metallicitySolar
                else
                   metallicity=energyInputTableMetallicity(iMetallicity)
                end if
                ! Update the counter.
                !$omp atomic
                loopCount=loopCount+1
                call Galacticus_Display_Counter(                                                   &
                     &                           int(100.0d0*dble(loopCount)/dble(loopCountTotal)) &
                     &                          ,.false.                                           &
                     &                          ,verbosityWorking                                  &
                     &                         )
                ! Find the minimum and maximum masses to integrate over for this IMF.
                minimumMass=IMF_Minimum_Mass(imfSelected)
                maximumMass=IMF_Maximum_Mass(imfSelected)
                ! Integrate cumulative energy input over the IMF between these limits.
                energyInputTable(iAge,iMetallicity,energyInputIndex(imfSelected))=Integrate(minimumMass,maximumMass&
                     &,Cumulative_Energy_Integrand,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0&
                     &,toleranceRelative=1.0d-3)
                call Integrate_Done(integrandFunction,integrationWorkspace)
             end do
          end do
          !$omp end parallel do
          do iAge=1,energyInputTableAgeCount
             do iMetallicity=1,energyInputTableMetallicityCount
                ! Enforce monotonicity in the cumulative energy input. Non-monotonicity can arise due to the vagaries of
                ! interpolating stellar lifetimes in an irregular grid of stellar models.
                if (iAge > 1) energyInputTable(iAge,iMetallicity,energyInputIndex(imfSelected))&
                     &=max(energyInputTable(iAge,iMetallicity,energyInputIndex(imfSelected)),energyInputTable(iAge-1&
                     &,iMetallicity,energyInputIndex(imfSelected)))
                call xml_NewElement(energyInputDoc,"data")
                write (parameterValue,'(e10.4)') energyInputTable(iAge,iMetallicity,energyInputIndex(imfSelected))
                call xml_AddCharacters(energyInputDoc,trim(parameterValue))
                call xml_EndElement(energyInputDoc,"data")
             end do
          end do
          call xml_EndElement(energyInputDoc,"energyInput")
          call Galacticus_Display_Counter_Clear(           verbosityWorking)
          call Galacticus_Display_Unindent     ('finished',verbosityWorking)
          call xml_EndElement(energyInputDoc,"stellarPopulation")
          call xml_Close(energyInputDoc)
       end if

       ! Flag that this IMF has now been tabulated.
       energyInputTabulated(imfSelected)=.true.
    end if
    !$omp end critical(IMF_Energy_Input_Rate_NonInstantaneous_Initialize)

    ! Get the index where this IMF is stored in the table.
    tableIndex=energyInputIndex(imfSelected)

    ! Interpolate to get the derivative in the recycled rate at two adjacent metallicities.
    metallicity=max(Abundances_Get_Metallicity(fuelAbundances),0.0d0)
    if (metallicity > energyInputTableMetallicityMaximum) then
       metallicityIndex=energyInputTableMetallicityCount
       metallicityFactors=[1.0d0,0.0d0]
       if (present(ageMaximum)) then
          ! Get average recycling rate between ageMinimum and ageMaximum.
          energyInputRate(1)=(Interpolate(energyInputTableAgeCount,energyInputTableAge,energyInputTable(: ,metallicityIndex&
               &,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMaximum,reset =interpolationAgeReset&
               &,extrapolationType=extrapolationTypeLinear)-Interpolate(energyInputTableAgeCount,energyInputTableAge ,energyInputTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset &
               &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear))/(ageMaximum-ageMinimum)
       else
          ! Get instantaneous energy input rate at ageMinimum.
          energyInputRate(1)=Interpolate_Derivative(energyInputTableAgeCount,energyInputTableAge,energyInputTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMinimum,reset &
               &=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
       end if
       energyInputRate(2)=0.0d0
    else
       metallicityIndex=Interpolate_Locate(energyInputTableMetallicityCount,energyInputTableMetallicity &
            &,interpolationMetallicityAccelerator,metallicity,reset=interpolationMetallicityReset)
       metallicityFactors=Interpolate_Linear_Generate_Factors(energyInputTableMetallicityCount,energyInputTableMetallicity &
            &,metallicityIndex,metallicity)
       ! Interpolate in age at both metallicities.
       do iMetallicity=0,1
          if (present(ageMaximum)) then
             ! Get average recycling rate between ageMinimum and ageMaximum.
             if (ageMinimum > 0.0d0) then
                energyInputMinimum=Interpolate(energyInputTableAgeCount&
                     &,energyInputTableAge ,energyInputTable(: ,metallicityIndex +iMetallicity,tableIndex),interpolationAgeObject&
                     &,interpolationAgeAccelerator ,ageMinimum ,reset =interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
             else
                energyInputMinimum=0.0d0
             end if
             energyInputMaximum=Interpolate(energyInputTableAgeCount,energyInputTableAge ,energyInputTable(: &
                  &,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject ,interpolationAgeAccelerator,ageMaximum &
                  &,reset =interpolationAgeReset,extrapolationType=extrapolationTypeLinear)

             energyInputRate(iMetallicity+1)=(energyInputMaximum-energyInputMinimum)/(ageMaximum-ageMinimum)
          else
             ! Get instantaneous recycling rate at ageMinimum.
             energyInputRate(iMetallicity+1)=Interpolate_Derivative(energyInputTableAgeCount,energyInputTableAge &
                  &,energyInputTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject &
                  &,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset,extrapolationType=extrapolationTypeLinear)
          end if
       end do
    end if

    ! Interpolate in metallicity to get the actual rate.
    IMF_Energy_Input_Rate_NonInstantaneous=sum(metallicityFactors*energyInputRate)

    return
  end function IMF_Energy_Input_Rate_NonInstantaneous

  function Cumulative_Energy_Integrand(initialMass,parameterPointer) bind(c)
    !% Integrand used in evaluating cumulative energy input.
    use, intrinsic :: ISO_C_Binding
    use Stellar_Feedback
    implicit none
    real(kind=c_double)        :: Cumulative_Energy_Integrand
    real(kind=c_double), value :: initialMass
    type(c_ptr        ), value :: parameterPointer
    real(kind=c_double)        :: energyLifetime

    if (starFormationImfInstantaneousApproximation) then
       ! In the instantaneous stellar evolution approximation, assume stars more massive than the long-lived star cut off
       ! contribute to the energy input (with an age equal to the specified effective age), while less massive stars contribute
       ! nothing.
       if (initialMass > starFormationImfInstantaneousApproximationMassLongLived) then
          energyLifetime=starFormationImfInstantaneousApproximationEffectiveAge
       else
          energyLifetime=0.0d0
       end if
    else
       ! In the standard calculation, simply use the current lifetime.
       energyLifetime=lifetime
    end if
    Cumulative_Energy_Integrand=IMF_Phi(initialMass,imfSelectedGlobal)*Stellar_Feedback_Cumulative_Energy_Input(initialMass&
         &,energyLifetime,metallicity)
    return
  end function Cumulative_Energy_Integrand

end module Star_Formation_IMF
