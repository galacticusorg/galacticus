!% Contains a module which implements functionality related to the stellar initial mass function.

module Star_Formation_IMF
  !% Implements functionality related to the stellar initial mass function.
  use ISO_Varying_String
  use Abundances_Structure
  use FGSL
  !# <include directive="imfRegister" type="moduleUse">
  include 'star_formation.IMF.register.modules.inc'
  !# </include>
  !# <include directive="imfSelectionMethod" type="moduleUse">
  !#  <subroutineArgs>imfSelectionMethod,IMF_Select,imfNames</subroutineArgs>
  include 'star_formation.IMF.select.modules.inc'
  !# </include>
  private
  public :: IMF_Select, IMF_Recycled_Fraction_Instantaneous, IMF_Recycling_Rate_NonInstantaneous, IMF_Yield_Instantaneous,&
       & IMF_Metal_Yield_Rate_NonInstantaneous, IMF_Energy_Input_Rate_NonInstantaneous, IMF_Name, IMF_Tabulate

  ! Flag to indicate if this module has been initialized.
  logical :: imfInitialized=.false.

  ! Count of the number of available IMFs.
  integer :: imfAvailableCount=0

  ! Array of IMF names.
  type(varying_string), allocatable, dimension(:    ) :: imfNames

  ! Tables of recycled fractions.
  logical,              allocatable, dimension(:    ) :: recycledFractionTabulated
  integer,              allocatable, dimension(:    ) :: recycledFractionIndex
  double precision,     allocatable, dimension(:    ) :: recycledFractionTableAge,recycledFractionTableMetallicity
  double precision,     allocatable, dimension(:,:,:) :: recycledFractionTable
  integer,              parameter                     :: recycledFractionTableMetallicityCount  =10
  integer,              parameter                     :: recycledFractionTableAgeCount          =50
  double precision,     parameter                     :: recycledFractionTableMetallicityMinimum=1.0d-4
  double precision,     parameter                     :: recycledFractionTableMetallicityMaximum=0.6d-1
  double precision,     parameter                     :: recycledFractionTableAgeMinimum        =1.0d-3
  double precision,     parameter                     :: recycledFractionTableAgeMaximum        =1.0d+2

  ! Tables of metal yields fractions.
  logical,              allocatable, dimension(:    ) :: metalYieldTabulated
  integer,              allocatable, dimension(:    ) :: metalYieldIndex
  double precision,     allocatable, dimension(:    ) :: metalYieldTableAge,metalYieldTableMetallicity
  double precision,     allocatable, dimension(:,:,:) :: metalYieldTable
  integer,              parameter                     :: metalYieldTableMetallicityCount  =10
  integer,              parameter                     :: metalYieldTableAgeCount          =50
  double precision,     parameter                     :: metalYieldTableMetallicityMinimum=1.0d-4
  double precision,     parameter                     :: metalYieldTableMetallicityMaximum=0.6d-1
  double precision,     parameter                     :: metalYieldTableAgeMinimum        =1.0d-3
  double precision,     parameter                     :: metalYieldTableAgeMaximum        =1.0d+2

  ! Tables of cumulative energy inputs.
  logical,              allocatable, dimension(:    ) :: energyInputTabulated
  integer,              allocatable, dimension(:    ) :: energyInputIndex
  double precision,     allocatable, dimension(:    ) :: energyInputTableAge,energyInputTableMetallicity
  double precision,     allocatable, dimension(:,:,:) :: energyInputTable
  integer,              parameter                     :: energyInputTableMetallicityCount  =10
  integer,              parameter                     :: energyInputTableAgeCount          =50
  double precision,     parameter                     :: energyInputTableMetallicityMinimum=1.0d-4
  double precision,     parameter                     :: energyInputTableMetallicityMaximum=0.6d-1
  double precision,     parameter                     :: energyInputTableAgeMinimum        =1.0d-3
  double precision,     parameter                     :: energyInputTableAgeMaximum        =1.0d+2

  ! Module global variables used in integration.
  integer                                             :: imfSelectedGlobal
  double precision                                    :: metallicity,lifetime

  ! Pointer to the function that selects which IMF to use.
  procedure(IMF_Select_Template), pointer :: IMF_Select => null()
  interface IMF_Select_Template
     integer function IMF_Select_Template(starFormationRate,fuelAbundances)
       import abundancesStructure
       double precision,          intent(in) :: starFormationRate
       type(abundancesStructure), intent(in) :: fuelAbundances
     end function IMF_Select_Template
  end interface

contains

  function IMF_Name(imfIndex)
    !% Return the name of the IMF with the specified index.
    use Galacticus_Error
    implicit none
    type(varying_string) :: IMF_Name
    integer, intent(in) :: imfIndex
    
    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    if (imfIndex <= size (imfNames)) then
       IMF_Name=imfNames(imfIndex)
    else
       call Galacticus_Error_Report('IMF_Name','imfIndex is out of range')
    end if
    return
  end function IMF_Name

  subroutine Star_Formation_IMF_Initialize
    !% Initialize the IMF subsystem.
    use Memory_Management
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type(varying_string) :: imfSelectionMethod

    ! Initialize the IMF subsystem if necessary.
    !$omp critical(IMF_Initialize)
    if (.not.imfInitialized) then
       ! Register all available IMFs.
       !# <include directive="imfRegister" type="code" action="subroutine">
       !#  <subroutineArgs>imfAvailableCount</subroutineArgs>
       include 'star_formation.IMF.register.inc'
       !# </include>

       ! Get a list of IMF names.
       call Alloc_Array(imfNames,imfAvailableCount,'imfNames')
       !# <include directive="imfRegisterName" type="code" action="subroutine">
       !#  <subroutineArgs>imfNames</subroutineArgs>
       include 'star_formation.IMF.register_names.inc'
       !# </include>

       ! Register the IMF selection method.
       !@ <inputParameter>
       !@   <name>imfSelectionMethod</name>
       !@   <defaultValue>fixed</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the method to be used for selecting which \IMF\ to use.
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('imfSelectionMethod',imfSelectionMethod,defaultValue='fixed')
       ! Include file that makes calls to all available method initialization routines.
       !# <include directive="imfSelectionMethod" type="code" action="subroutine">
       !#  <subroutineArgs>imfSelectionMethod,IMF_Select,imfNames</subroutineArgs>
       include 'star_formation.IMF.select.inc'
       !# </include>
       if (.not.associated(IMF_Select)) call Galacticus_Error_Report('Star_Formation_IMF_Initialize'&
            &,'method '//char(imfSelectionMethod)//' is unrecognized')

       ! Flag that the module is now initialized.
       imfInitialized=.true.
    end if
    !$omp end critical(IMF_Initialize)
    return
  end subroutine Star_Formation_IMF_Initialize

  double precision function IMF_Recycled_Fraction_Instantaneous(starFormationRate,fuelAbundances)
    !% Returns a recycled fraction for the IMF suitable for use in the instantaneous recycling approximation.
    use Abundances_Structure
    implicit none
    double precision,          intent(in) :: starFormationRate
    type(abundancesStructure), intent(in) :: fuelAbundances
    integer                               :: imfSelected
    logical                               :: imfMatched

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Get the recycled fraction from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfRecycledInstantaneous" type="code" action="subroutine">
    !#  <subroutineArgs>imfSelected,imfMatched,IMF_Recycled_Fraction_Instantaneous</subroutineArgs>
    !#  <subroutineAction>if (imfMatched) return</subroutineAction>
    include 'star_formation.IMF.recycled_instantaneous.inc'
    !# </include>
    return
  end function IMF_Recycled_Fraction_Instantaneous

  double precision function IMF_Yield_Instantaneous(starFormationRate,fuelAbundances)
    !% Returns a yield for the IMF suitable for use in the instantaneous recycling approximation.
    use Abundances_Structure
    implicit none
    double precision,          intent(in) :: starFormationRate
    type(abundancesStructure), intent(in) :: fuelAbundances
    integer                               :: imfSelected
    logical                               :: imfMatched

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Get the yield from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfYieldInstantaneous" type="code" action="subroutine">
    !#  <subroutineArgs>imfSelected,imfMatched,IMF_Yield_Instantaneous</subroutineArgs>
    !#  <subroutineAction>if (imfMatched) return</subroutineAction>
    include 'star_formation.IMF.yield_instantaneous.inc'
    !# </include>
    return
  end function IMF_Yield_Instantaneous

  subroutine IMF_Tabulate(imfIndex,imfMass,imfPhi)
    !% Returns a tabulation of the IMF with sufficient resolution to resolve all features.
    implicit none
    integer,          intent(in)                               :: imfIndex
    double precision, intent(inout), allocatable, dimension(:) :: imfMass,imfPhi
    logical                                                    :: imfMatched

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Get the recycled fraction from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfTabulate" type="code" action="subroutine">
    !#  <subroutineArgs>imfIndex,imfMatched,imfMass,imfPhi</subroutineArgs>
    !#  <subroutineAction>if (imfMatched) return</subroutineAction>
    include 'star_formation.IMF.tabulate.inc'
    !# </include>
    return
  end subroutine IMF_Tabulate

  double precision function IMF_Minimum_Mass(imfSelected)
    !% Returns the minimum mass in the selected IMF.
    implicit none
    integer, intent(in) :: imfSelected
    logical             :: imfMatched

    ! Get the minimum mass from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfMinimumMass" type="code" action="subroutine">
    !#  <subroutineArgs>imfSelected,imfMatched,IMF_Minimum_Mass</subroutineArgs>
    !#  <subroutineAction>if (imfMatched) return</subroutineAction>
    include 'star_formation.IMF.minimum_mass.inc'
    !# </include>
    return
  end function IMF_Minimum_Mass

  double precision function IMF_Maximum_Mass(imfSelected)
    !% Returns the maximum mass in the selected IMF.
    implicit none
    integer, intent(in) :: imfSelected
    logical             :: imfMatched

    ! Get the maximum mass from the appropriate IMF.
    imfMatched=.false.
    !# <include directive="imfMaximumMass" type="code" action="subroutine">
    !#  <subroutineArgs>imfSelected,imfMatched,IMF_Maximum_Mass</subroutineArgs>
    !#  <subroutineAction>if (imfMatched) return</subroutineAction>
    include 'star_formation.IMF.maximum_mass.inc'
    !# </include>
    return
  end function IMF_Maximum_Mass

  double precision function IMF_Phi(initialMass,imfSelected)
    !% Returns the IMF, $\Phi(M)$, at mass $M=${\tt initialMass} for the selected IMF.
    implicit none
    integer,          intent(in) :: imfSelected
    double precision, intent(in) :: initialMass
    logical                      :: imfMatched

    ! Get the IMF for the selected IMF and initial stellar mass.
    imfMatched=.false.
    !# <include directive="imfPhi" type="code" action="subroutine">
    !#  <subroutineArgs>imfSelected,imfMatched,initialMass,IMF_Phi</subroutineArgs>
    !#  <subroutineAction>if (imfMatched) return</subroutineAction>
    include 'star_formation.IMF.phi.inc'
    !# </include>
    return
  end function IMF_Phi

  double precision function IMF_Recycling_Rate_NonInstantaneous(starFormationRate,fuelAbundances,ageMinimum,ageMaximum)
    !% Returns the recycling rate for a simple stellar population. The \IMF\ is determined from the given {\tt starFormationRate}
    !% and {\tt fuelAbundances}. The recycling rate (in the fraction of the population's mass returned to the \ISM\ per Gyr) is
    !% computed for the given {\tt age} (in Gyr). The recycled fraction is computed on a grid of age and metallicity. This is
    !% stored to file and will be read back in on subsequent runs. This is useful as computation of the table is relatively slow.
    use, intrinsic :: ISO_C_Binding
    use Abundances_Structure
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Ranges
    use Stellar_Astrophysics
    use Memory_Management
    use Galacticus_Display
    use File_Utilities
    use FoX_wxml
    use FoX_dom
    use ISO_Varying_String
    use Galacticus_Error
    use Dates_and_Times
    implicit none
    double precision,          intent(in)                    :: starFormationRate,ageMinimum
    double precision,          intent(in),  optional         :: ageMaximum
    type(abundancesStructure), intent(in)                    :: fuelAbundances
    logical,                   allocatable, dimension(:    ) :: recycledFractionTabulatedTemporary
    integer,                   allocatable, dimension(:    ) :: recycledFractionIndexTemporary
    double precision,          allocatable, dimension(:,:,:) :: recycledFractionTableTemporary
    double precision,                       dimension(2    ) :: recycleRate,metallicityFactors
    type(fgsl_interp),         save                          ::                                      interpolationAgeObject
    type(fgsl_interp_accel),   save                          :: interpolationMetallicityAccelerator ,interpolationAgeAccelerator
    logical,                   save                          :: interpolationMetallicityReset=.true.,interpolationAgeReset=.true.
    !$omp threadprivate(interpolationAgeObject,interpolationMetallicityAccelerator &
    !$omp ,interpolationAgeAccelerator,interpolationMetallicityReset,interpolationAgeReset)
    type(Node),                pointer                       :: doc
    type(NodeList),            pointer                       :: columnList,dataList
    type(c_ptr)                                              :: parameterPointer
    type(fgsl_function)                                      :: integrandFunction
    type(fgsl_integration_workspace)                         :: integrationWorkspace
    integer                                                  :: imfSelected,iAge,iMetallicity,imfCount,tableIndex&
         &,metallicityIndex,iRecycledFraction,ioErr
    double precision                                         :: minimumMass,maximumMass,initialMass
    character(len=20)                                        :: progressMessage,parameterValue
    type(xmlf_t)                                             :: recycledFractionDoc
    type(varying_string)                                     :: fileName

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Check that flag and index arrays exist.
    if (.not.allocated(recycledFractionTabulated)) then
       call Alloc_Array(recycledFractionTabulated,imfSelected,'recycledFractionTabulated')
       call Alloc_Array(recycledFractionIndex    ,imfSelected,'recycledFractionIndex'    )
       recycledFractionTabulated=.false.
       recycledFractionIndex    =0
    end if

    ! Check that flag and index arrays are large enough.
    if (size(recycledFractionTabulated) < imfSelected) then
       call Move_Alloc (recycledFractionTabulated,recycledFractionTabulatedTemporary)
       call Move_Alloc (recycledFractionIndex    ,recycledFractionIndexTemporary    )
       call Alloc_Array(recycledFractionTabulated,imfSelected,'recycledFractionTabulated')
       call Alloc_Array(recycledFractionIndex    ,imfSelected,'recycledFractionIndex'    )
       recycledFractionTabulated(1:size(recycledFractionTabulatedTemporary))            =recycledFractionTabulatedTemporary
       recycledFractionIndex    (1:size(recycledFractionIndexTemporary    ))            =recycledFractionIndexTemporary
       recycledFractionTabulated(  size(recycledFractionTabulatedTemporary):imfSelected)=.false.
       recycledFractionIndex    (  size(recycledFractionTabulatedTemporary):imfSelected)=0
       call Dealloc_Array(recycledFractionTabulatedTemporary)
       call Dealloc_Array(recycledFractionIndexTemporary    )
    end if

    ! Tabulate the recycled fraction for this IMF if it has not already been computed.
    if (.not.recycledFractionTabulated(imfSelected)) then
       
       ! Expand the tabulations array by enough to accomodate a new IMF.
       if (allocated(recycledFractionTable)) then
          imfCount=size(recycledFractionTable,dim=3)
          call Move_Alloc(recycledFractionTable,recycledFractionTableTemporary)
          call Alloc_Array(recycledFractionTable,recycledFractionTableAgeCount,recycledFractionTableMetallicityCount,imfCount,'recycledFractionTable')
          recycledFractionTable(:,:,1:imfCount)=recycledFractionTableTemporary
          call Dealloc_Array(recycledFractionTableTemporary)
       else
          call Alloc_Array(recycledFractionTableAge        ,recycledFractionTableAgeCount        ,'recycledFractionTableAge'        )
          call Alloc_Array(recycledFractionTableMetallicity,recycledFractionTableMetallicityCount,'recycledFractionTableMetallicity')
          recycledFractionTableAge                                                 =Make_Range(recycledFractionTableAgeMinimum&
               &,recycledFractionTableAgeMaximum,recycledFractionTableAgeCount,rangeType=rangeTypeLogarithmic)
          recycledFractionTableMetallicity(1)                                      =0.0d0
          recycledFractionTableMetallicity(2:recycledFractionTableMetallicityCount)&
               &=Make_Range(recycledFractionTableMetallicityMinimum,recycledFractionTableMetallicityMaximum&
               &,recycledFractionTableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          call Alloc_Array(recycledFractionTable,recycledFractionTableAgeCount,recycledFractionTableMetallicityCount,1,'recycledFractionTable')
       end if
       
       ! Record the index in the array where this IMF will be stored.
       recycledFractionIndex(imfSelected)=size(recycledFractionTable,dim=3)

       ! Check if the table has been computed and stored previously.
       fileName='./data/Stellar_Recycled_Fraction_'//imfNames(imfSelected)//'.xml'
       if (File_Exists(fileName)) then
          
          ! Open the XML file containing recycled fractions.
          doc => parseFile(char(fileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','Unable to parse recycled fractions file')
          
          ! Find the ages element and extract data.
          columnList => getElementsByTagname(doc               ,"ages")
          dataList   => getElementsByTagname(item(columnList,0),"data")
          if (getLength(dataList) /= recycledFractionTableAgeCount) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','ages array in XML file does not match internal expectation')
          do iAge=1,recycledFractionTableAgeCount
             call extractDataContent(item(dataList,iAge-1),lifetime)
             if (recycledFractionIndex(imfSelected) == 1) then
                recycledFractionTableAge(iAge)=lifetime
             else
                if (recycledFractionTableAge(iAge) /= lifetime) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','mismatch in ages array in XML file')
             end if
          end do

          ! Find the metallicities element and extract data.
          columnList => getElementsByTagname(doc               ,"metallicities")
          dataList   => getElementsByTagname(item(columnList,0),"data"         )
          if (getLength(dataList) /= recycledFractionTableMetallicityCount) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','metallicities array in XML file does not match internal expectation')
          do iMetallicity=1,recycledFractionTableMetallicityCount
             call extractDataContent(item(dataList,iMetallicity-1),metallicity)
             if (recycledFractionIndex(imfSelected) == 1) then
                recycledFractionTableMetallicity(iMetallicity)=metallicity
             else
                if (recycledFractionTableMetallicity(iMetallicity) /= metallicity) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','mismatch in metallicities array in XML file')
             end if
          end do

          ! Find the recycledFraction element and extract data.
          columnList => getElementsByTagname(doc               ,"recycledFraction")
          dataList   => getElementsByTagname(item(columnList,0),"data"            )
          if (getLength(dataList) /= recycledFractionTableAgeCount*recycledFractionTableMetallicityCount) call Galacticus_Error_Report('IMF_Recycling_Rate_NonInstantaneous','recycled fractions array in XML file does not match internal expectation')
          iRecycledFraction=0
          do iAge=1,recycledFractionTableAgeCount
             do iMetallicity=1,recycledFractionTableMetallicityCount
                iRecycledFraction=iRecycledFraction+1
                call extractDataContent(item(dataList,iRecycledFraction-1),recycledFractionTable(iAge,iMetallicity&
                     &,recycledFractionIndex(imfSelected)))
             end do
          end do

       else

          call Galacticus_Display_Indent('Tabulating mass recycling rate for '//char(imfNames(imfSelected))//' IMF',2)
          
          ! Open an XML file to output the data to.
          call xml_OpenFile(char(fileName),recycledFractionDoc)
          call xml_NewElement(recycledFractionDoc,"stellarPopulation")
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
          do iAge=1,recycledFractionTableAgeCount
             lifetime=recycledFractionTableAge(iAge)
             write (progressMessage,'(a6,e8.2,a4)') 'age = ',lifetime,' Gyr'
             call Galacticus_Display_Message(progressMessage,3)
             do iMetallicity=1,recycledFractionTableMetallicityCount
                metallicity=recycledFractionTableMetallicity(iMetallicity)
                ! Find the minimum and maximum masses to integrate over for this IMF.
                minimumMass=IMF_Minimum_Mass(imfSelected)
                maximumMass=IMF_Maximum_Mass(imfSelected)
                ! Integrate ejected mass over the IMF between these limits.
                recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected))=Integrate(minimumMass,maximumMass&
                     &,Recycled_Fraction_Integrand ,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=1.0d-3&
                     &,toleranceRelative=1.0d-4)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                ! Enforce monotonicity in the recycled fraction. Non-monotonicity can arise due to the vagaries of interpolating stellar
                ! lifetimes in an irregular grid of stellar models.
                if (iAge > 1) recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected))&
                     &=max(recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected)),recycledFractionTable(iAge-1&
                     &,iMetallicity,recycledFractionIndex(imfSelected)))
                call xml_NewElement(recycledFractionDoc,"data")
                write (parameterValue,'(e10.4)') recycledFractionTable(iAge,iMetallicity,recycledFractionIndex(imfSelected))
                call xml_AddCharacters(recycledFractionDoc,trim(parameterValue))
                call xml_EndElement(recycledFractionDoc,"data")
             end do
          end do
          call xml_EndElement(recycledFractionDoc,"recycledFraction")
          call Galacticus_Display_Unindent('finished',2)
          call xml_EndElement(recycledFractionDoc,"stellarPopulation")
          call xml_Close(recycledFractionDoc)
       end if
       
       ! Flag that this IMF has now been tabulated.
       recycledFractionTabulated(imfSelected)=.true.
    end if

    ! Get the index where this IMF is stored in the table.
    tableIndex=recycledFractionIndex(imfSelected)
    
    ! Interpolate to get the derivative in the recycled rate at two adjacent metallicities.
    metallicity=Abundances_Get_Metallicity(fuelAbundances)
    if (metallicity > recycledFractionTableMetallicityMaximum) then
       metallicityIndex=recycledFractionTableMetallicityCount
       metallicityFactors=[1.0d0,0.0d0]
       if (present(ageMaximum)) then
          ! Get average recycling rate between ageMinimum and ageMaximum.
          recycleRate(1)=(Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge,recycledFractionTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMaximum,reset&
               &=interpolationAgeReset)-Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge&
               &,recycledFractionTable(: ,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator&
               &,ageMinimum,reset=interpolationAgeReset))/(ageMaximum-ageMinimum)
       else
          ! Get instantaneous recycling rate at ageMinimum.
          recycleRate(1)=Interpolate_Derivative(recycledFractionTableAgeCount,recycledFractionTableAge,recycledFractionTable(:&
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset)
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
             recycleRate(iMetallicity+1)=(Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge&
                  &,recycledFractionTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator,ageMaximum,reset =interpolationAgeReset)&
                  &-Interpolate(recycledFractionTableAgeCount,recycledFractionTableAge ,recycledFractionTable(:,metallicityIndex&
                  &+iMetallicity,tableIndex),interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset&
                  &=interpolationAgeReset))/(ageMaximum-ageMinimum)
          else
             ! Get instantaneous recycling rate at ageMinimum.
             recycleRate(iMetallicity+1)=Interpolate_Derivative(recycledFractionTableAgeCount,recycledFractionTableAge&
                  &,recycledFractionTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset)
          end if
       end do
    end if
    
    ! Interpolate in metallicity to get the actual rate.
    IMF_Recycling_Rate_NonInstantaneous=sum(metallicityFactors*recycleRate)

    return
  end function IMF_Recycling_Rate_NonInstantaneous

  function Recycled_Fraction_Integrand(initialMass,parameterPointer) bind(c)
    !% Integrand used in evaluating recycled fractions.
    use, intrinsic :: ISO_C_Binding
    use Stellar_Astrophysics
    implicit none
    real(c_double)        :: Recycled_Fraction_Integrand
    real(c_double), value :: initialMass
    type(c_ptr),    value :: parameterPointer

    if (Star_Lifetime(initialMass,metallicity) <= lifetime) then
       Recycled_Fraction_Integrand=IMF_Phi(initialMass,imfSelectedGlobal)*Star_Ejected_Mass(initialMass,metallicity)
    else
       Recycled_Fraction_Integrand=0.0d0
    end if
    return
  end function Recycled_Fraction_Integrand
  
  double precision function IMF_Metal_Yield_Rate_NonInstantaneous(starFormationRate,fuelAbundances,ageMinimum,ageMaximum)
     !% Returns the metal yield rate for a simple stellar population. The \IMF\ is determined from the given {\tt starFormationRate}
     !% and {\tt fuelAbundances}. The metal yield rate (in fraction of the population's mass returned to the \ISM\ as new metals per Gyr) is
     !% computed for the given {\tt age} (in Gyr). The metal yield is computed on a grid of age and metallicity. This is
     !% stored to file and will be read back in on subsequent runs. This is useful as computation of the table is relatively slow.
     use, intrinsic :: ISO_C_Binding
     use Abundances_Structure
     use Numerical_Integration
     use Numerical_Interpolation
     use Numerical_Ranges
     use Stellar_Astrophysics
     use Memory_Management
     use Galacticus_Display
     use File_Utilities
     use FoX_wxml
     use FoX_dom
     use ISO_Varying_String
     use Galacticus_Error
     use Dates_and_Times
     implicit none
     double precision,          intent(in)                    :: starFormationRate,ageMinimum
     double precision,          intent(in),  optional         :: ageMaximum
     type(abundancesStructure), intent(in)                    :: fuelAbundances
     logical,                   allocatable, dimension(:    ) :: metalYieldTabulatedTemporary
     integer,                   allocatable, dimension(:    ) :: metalYieldIndexTemporary
     double precision,          allocatable, dimension(:,:,:) :: metalYieldTableTemporary
     double precision,                       dimension(2    ) :: metalYieldRate,metallicityFactors
     type(fgsl_interp),         save                          ::                                      interpolationAgeObject
     type(fgsl_interp_accel),   save                          :: interpolationMetallicityAccelerator ,interpolationAgeAccelerator
     logical,                   save                          :: interpolationMetallicityReset=.true.,interpolationAgeReset=.true.
     !$omp threadprivate(interpolationAgeObject,interpolationMetallicityAccelerator &
     !$omp ,interpolationAgeAccelerator,interpolationMetallicityReset,interpolationAgeReset)
     type(Node),                pointer                       :: doc
     type(NodeList),            pointer                       :: columnList,dataList
     type(c_ptr)                                              :: parameterPointer
     type(fgsl_function)                                      :: integrandFunction
     type(fgsl_integration_workspace)                         :: integrationWorkspace
     integer                                                  :: imfSelected,iAge,iMetallicity,imfCount,tableIndex&
          &,metallicityIndex,iMetalYield,ioErr
     double precision                                         :: minimumMass,maximumMass,initialMass
     character(len=20)                                        :: progressMessage,parameterValue
     character(len= 8)                                        :: date
     character(len=10)                                        :: time
     character(len= 5)                                        :: zone
     type(xmlf_t)                                             :: metalYieldDoc
     type(varying_string)                                     :: fileName

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Check that flag and index arrays exist.
    if (.not.allocated(metalYieldTabulated)) then
       call Alloc_Array(metalYieldTabulated,imfSelected,'metalYieldTabulated')
       call Alloc_Array(metalYieldIndex    ,imfSelected,'metalYieldIndex'    )
       metalYieldTabulated=.false.
       metalYieldIndex    =0
    end if

    ! Check that flag and index arrays are large enough.
    if (size(metalYieldTabulated) < imfSelected) then
       call Move_Alloc (metalYieldTabulated,metalYieldTabulatedTemporary)
       call Move_Alloc (metalYieldIndex    ,metalYieldIndexTemporary    )
       call Alloc_Array(metalYieldTabulated,imfSelected,'metalYieldTabulated')
       call Alloc_Array(metalYieldIndex    ,imfSelected,'metalYieldIndex'    )
       metalYieldTabulated(1:size(metalYieldTabulatedTemporary))            =metalYieldTabulatedTemporary
       metalYieldIndex    (1:size(metalYieldIndexTemporary    ))            =metalYieldIndexTemporary
       metalYieldTabulated(  size(metalYieldTabulatedTemporary):imfSelected)=.false.
       metalYieldIndex    (  size(metalYieldTabulatedTemporary):imfSelected)=0
       call Dealloc_Array(metalYieldTabulatedTemporary)
       call Dealloc_Array(metalYieldIndexTemporary    )
    end if

    ! Tabulate the metal yield for this IMF if it has not already been computed.
    if (.not.metalYieldTabulated(imfSelected)) then
     
       ! Expand the tabulations array by enough to accomodate a new IMF.
       if (allocated(metalYieldTable)) then
          imfCount=size(metalYieldTable,dim=3)
          call Move_Alloc(metalYieldTable,metalYieldTableTemporary)
          call Alloc_Array(metalYieldTable,metalYieldTableAgeCount,metalYieldTableMetallicityCount,imfCount,'metalYieldTable')
          metalYieldTable(:,:,1:imfCount)=metalYieldTableTemporary
          call Dealloc_Array(metalYieldTableTemporary)
       else
          call Alloc_Array(metalYieldTableAge        ,metalYieldTableAgeCount        ,'metalYieldTableAge'        )
          call Alloc_Array(metalYieldTableMetallicity,metalYieldTableMetallicityCount,'metalYieldTableMetallicity')
          metalYieldTableAge                                           =Make_Range(metalYieldTableAgeMinimum&
               &,metalYieldTableAgeMaximum,metalYieldTableAgeCount,rangeType=rangeTypeLogarithmic)
          metalYieldTableMetallicity(1)                                =0.0d0
          metalYieldTableMetallicity(2:metalYieldTableMetallicityCount)&
               &=Make_Range(metalYieldTableMetallicityMinimum,metalYieldTableMetallicityMaximum&
               &,metalYieldTableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          call Alloc_Array(metalYieldTable,metalYieldTableAgeCount,metalYieldTableMetallicityCount,1,'metalYieldTable')
       end if
     
       ! Record the index in the array where this IMF will be stored.
       metalYieldIndex(imfSelected)=size(metalYieldTable,dim=3)

       ! Check if the table has been computed and stored previously.
       fileName='./data/Stellar_Metal_Yield_'//imfNames(imfSelected)//'.xml'
       if (File_Exists(fileName)) then
        
          ! Open the XML file containing metal yields.
          doc => parseFile(char(fileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','Unable to parse metal yields file')
        
          ! Find the ages element and extract data.
          columnList => getElementsByTagname(doc               ,"ages")
          dataList   => getElementsByTagname(item(columnList,0),"data")
          if (getLength(dataList) /= metalYieldTableAgeCount) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','ages array in XML file does not match internal expectation')
          do iAge=1,metalYieldTableAgeCount
             call extractDataContent(item(dataList,iAge-1),lifetime)
             if (metalYieldIndex(imfSelected) == 1) then
                metalYieldTableAge(iAge)=lifetime
             else
                if (metalYieldTableAge(iAge) /= lifetime) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','mismatch in ages array in XML file')
             end if
          end do

          ! Find the metallicities element and extract data.
          columnList => getElementsByTagname(doc               ,"metallicities")
          dataList   => getElementsByTagname(item(columnList,0),"data"         )
          if (getLength(dataList) /= metalYieldTableMetallicityCount) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','metallicities array in XML file does not match internal expectation')
          do iMetallicity=1,metalYieldTableMetallicityCount
             call extractDataContent(item(dataList,iMetallicity-1),metallicity)
             if (metalYieldIndex(imfSelected) == 1) then
                metalYieldTableMetallicity(iMetallicity)=metallicity
             else
                if (metalYieldTableMetallicity(iMetallicity) /= metallicity) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','mismatch in metallicities array in XML file')
             end if
          end do

          ! Find the metalYield element and extract data.
          columnList => getElementsByTagname(doc               ,"metalYield")
          dataList   => getElementsByTagname(item(columnList,0),"data"      )
          if (getLength(dataList) /= metalYieldTableAgeCount*metalYieldTableMetallicityCount) call Galacticus_Error_Report('IMF_Metal_Yield_Rate_NonInstantaneous','metal yield array in XML file does not match internal expectation')
          iMetalYield=0
          do iAge=1,metalYieldTableAgeCount
             do iMetallicity=1,metalYieldTableMetallicityCount
                iMetalYield=iMetalYield+1
                call extractDataContent(item(dataList,iMetalYield-1),metalYieldTable(iAge,iMetallicity &
                     &,metalYieldIndex(imfSelected)))
             end do
          end do

       else

          call Galacticus_Display_Indent('Tabulating metal yield rate for '//char(imfNames(imfSelected))//' IMF',2)
        
          ! Open an XML file to output the data to.
          call xml_OpenFile(char(fileName),metalYieldDoc)
          call xml_NewElement(metalYieldDoc,"stellarPopulation")
          call xml_NewElement(metalYieldDoc,"description")
          call xml_AddCharacters(metalYieldDoc,"Metal yield for a "//char(imfNames(imfSelected))//" IMF")
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
          call xml_NewElement(metalYieldDoc,"metalYield")
          do iAge=1,metalYieldTableAgeCount
             lifetime=metalYieldTableAge(iAge)
             write (progressMessage,'(a6,e8.2,a4)') 'age = ',lifetime,' Gyr'
             call Galacticus_Display_Message(progressMessage,3)
             do iMetallicity=1,metalYieldTableMetallicityCount
                metallicity=metalYieldTableMetallicity(iMetallicity)
                ! Find the minimum and maximum masses to integrate over for this IMF.
                minimumMass=IMF_Minimum_Mass(imfSelected)
                maximumMass=IMF_Maximum_Mass(imfSelected)
                ! Integrate ejected mass over the IMF between these limits.                
                metalYieldTable(iAge,iMetallicity,metalYieldIndex(imfSelected))=Integrate(minimumMass,maximumMass&
                     &,Metal_Yield_Integrand ,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=1.0d-3&
                     &,toleranceRelative=1.0d-4)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                ! Enforce monotonicity in the metal yield. Non-monotonicity can arise due to the vagaries of interpolating stellar
                ! lifetimes in an irregular grid of stellar models.
                if (iAge > 1) metalYieldTable(iAge,iMetallicity,metalYieldIndex(imfSelected))&
                     &=max(metalYieldTable(iAge,iMetallicity,metalYieldIndex(imfSelected)),metalYieldTable(iAge-1&
                     &,iMetallicity,metalYieldIndex(imfSelected)))
                call xml_NewElement(metalYieldDoc,"data")
                write (parameterValue,'(e10.4)') metalYieldTable(iAge,iMetallicity,metalYieldIndex(imfSelected))
                call xml_AddCharacters(metalYieldDoc,trim(parameterValue))
                call xml_EndElement(metalYieldDoc,"data")
             end do
          end do
          call xml_EndElement(metalYieldDoc,"metalYield")
          call Galacticus_Display_Unindent('finished',2)
          call xml_EndElement(metalYieldDoc,"stellarPopulation")
          call xml_Close(metalYieldDoc)
       end if
     
       ! Flag that this IMF has now been tabulated.
       metalYieldTabulated(imfSelected)=.true.
    end if

    ! Get the index where this IMF is stored in the table.
    tableIndex=metalYieldIndex(imfSelected)

     ! Interpolate to get the derivative in the metal yield at two adjacent metallicities.
     metallicity=Abundances_Get_Metallicity(fuelAbundances)
     if (metallicity > metalYieldTableMetallicityMaximum) then
        metallicityIndex=metalYieldTableMetallicityCount
        metallicityFactors=[1.0d0,0.0d0]
        if (present(ageMaximum)) then
           ! Get average recycling rate between ageMinimum and ageMaximum.
           metalYieldRate(1)=(Interpolate(metalYieldTableAgeCount,metalYieldTableAge,metalYieldTable(: ,metallicityIndex&
                &,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMaximum,reset =interpolationAgeReset)&
                &-Interpolate(metalYieldTableAgeCount,metalYieldTableAge ,metalYieldTable(: ,metallicityIndex,tableIndex)&
                &,interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset=interpolationAgeReset))/(ageMaximum&
                &-ageMinimum)
        else
           ! Get instantaneous recycling rate at ageMinimum.
           metalYieldRate(1)=Interpolate_Derivative(metalYieldTableAgeCount,metalYieldTableAge,metalYieldTable(: &
                &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset)
        end if
        metalYieldRate(2)=0.0d0
     else
        metallicityIndex=Interpolate_Locate(metalYieldTableMetallicityCount,metalYieldTableMetallicity&
             &,interpolationMetallicityAccelerator,metallicity,reset=interpolationMetallicityReset)
        metallicityFactors=Interpolate_Linear_Generate_Factors(metalYieldTableMetallicityCount,metalYieldTableMetallicity&
             &,metallicityIndex,metallicity)
        ! Interpolate in age at both metallicities.
        do iMetallicity=0,1
          if (present(ageMaximum)) then
             ! Get average recycling rate between ageMinimum and ageMaximum.
             metalYieldRate(iMetallicity+1)=(Interpolate(metalYieldTableAgeCount,metalYieldTableAge&
                  &,metalYieldTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator,ageMaximum,reset =interpolationAgeReset)&
                  &-Interpolate(metalYieldTableAgeCount,metalYieldTableAge ,metalYieldTable(:,metallicityIndex&
                  &+iMetallicity,tableIndex),interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset&
                  &=interpolationAgeReset))/(ageMaximum-ageMinimum)
          else
             ! Get instantaneous recycling rate at ageMinimum.
             metalYieldRate(iMetallicity+1)=Interpolate_Derivative(metalYieldTableAgeCount,metalYieldTableAge&
                  &,metalYieldTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject&
                  &,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset)
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
    real(c_double)        :: Metal_Yield_Integrand
    real(c_double), value :: initialMass
    type(c_ptr),    value :: parameterPointer

    ! Include yields from isolated stars.
    if (Star_Lifetime(initialMass,metallicity) <= lifetime) then
       Metal_Yield_Integrand=IMF_Phi(initialMass,imfSelectedGlobal)*Star_Metal_Yield_Mass(initialMass,metallicity)
    else
       Metal_Yield_Integrand=0.0d0
    end if

    ! Include yield from Type Ia supernovae.
    Metal_Yield_Integrand=Metal_Yield_Integrand+IMF_Phi(initialMass,imfSelectedGlobal)*SNeIa_Cumulative_Yield(initialMass&
         &,lifetime,metallicity)
    return
  end function Metal_Yield_Integrand

  double precision function IMF_Energy_Input_Rate_NonInstantaneous(starFormationRate,fuelAbundances,ageMinimum,ageMaximum)
    !% Returns the energy input rate for a simple stellar population in (km/s)$^2$ Gyr$^{-1}$. The \IMF\ is determined from the
    !% given {\tt starFormationRate} and {\tt fuelAbundances}. The energy input rate is computed for the given {\tt age} (in
    !% Gyr). The cumulative energy input is computed on a grid of age and metallicity. This is stored to file and will be read
    !% back in on subsequent runs. This is useful as computation of the table is relatively slow.
    use, intrinsic :: ISO_C_Binding
    use Abundances_Structure
    use Numerical_Integration
    use Numerical_Interpolation
    use Numerical_Ranges
    use Stellar_Astrophysics
    use Memory_Management
    use Galacticus_Display
    use File_Utilities
    use FoX_wxml
    use FoX_dom
    use ISO_Varying_String
    use Galacticus_Error
    use Dates_and_Times
    implicit none
    double precision,          intent(in)                    :: starFormationRate,ageMinimum
    double precision,          intent(in),  optional         :: ageMaximum
    type(abundancesStructure), intent(in)                    :: fuelAbundances
    logical,                   allocatable, dimension(:    ) :: energyInputTabulatedTemporary
    integer,                   allocatable, dimension(:    ) :: energyInputIndexTemporary
    double precision,          allocatable, dimension(:,:,:) :: energyInputTableTemporary
    double precision,                       dimension(2    ) :: energyInputRate,metallicityFactors
    type(fgsl_interp),         save                          ::                                      interpolationAgeObject
    type(fgsl_interp_accel),   save                          :: interpolationMetallicityAccelerator ,interpolationAgeAccelerator
    logical,                   save                          :: interpolationMetallicityReset=.true.,interpolationAgeReset=.true.
    !$omp threadprivate(interpolationAgeObject,interpolationMetallicityAccelerator &
    !$omp ,interpolationAgeAccelerator,interpolationMetallicityReset,interpolationAgeReset)
    type(Node),                pointer                       :: doc
    type(NodeList),            pointer                       :: columnList,dataList
    type(c_ptr)                                              :: parameterPointer
    type(fgsl_function)                                      :: integrandFunction
    type(fgsl_integration_workspace)                         :: integrationWorkspace
    integer                                                  :: imfSelected,iAge,iMetallicity,imfCount,tableIndex&
         &,metallicityIndex,iEnergyInput,ioErr
    double precision                                         :: minimumMass,maximumMass,initialMass,totalEnergyInput
    character(len=20)                                        :: progressMessage,parameterValue
    type(xmlf_t)                                             :: energyInputDoc
    type(varying_string)                                     :: fileName

    ! Initialize the IMF subsystem.
    call Star_Formation_IMF_Initialize

    ! Select which IMF to use.
    imfSelected=IMF_Select(starFormationRate,fuelAbundances)

    ! Check that flag and index arrays exist.
    if (.not.allocated(energyInputTabulated)) then
       call Alloc_Array(energyInputTabulated,imfSelected,'energyInputTabulated')
       call Alloc_Array(energyInputIndex    ,imfSelected,'energyInputIndex'    )
       energyInputTabulated=.false.
       energyInputIndex    =0
    end if

    ! Check that flag and index arrays are large enough.
    if (size(energyInputTabulated) < imfSelected) then
       call Move_Alloc (energyInputTabulated,energyInputTabulatedTemporary)
       call Move_Alloc (energyInputIndex    ,energyInputIndexTemporary    )
       call Alloc_Array(energyInputTabulated,imfSelected,'energyInputTabulated')
       call Alloc_Array(energyInputIndex    ,imfSelected,'energyInputIndex'    )
       energyInputTabulated(1:size(energyInputTabulatedTemporary))=energyInputTabulatedTemporary
       energyInputIndex    (1:size(energyInputIndexTemporary    ))=energyInputIndexTemporary
       energyInputTabulated(  size(energyInputTabulatedTemporary):imfSelected)=.false.
       energyInputIndex    (  size(energyInputTabulatedTemporary):imfSelected)=0
       call Dealloc_Array(energyInputTabulatedTemporary)
       call Dealloc_Array(energyInputIndexTemporary    )
    end if

    ! Tabulate the cumulative energy input for this IMF if it has not already been computed.
    if (.not.energyInputTabulated(imfSelected)) then
       
       ! Expand the tabulations array by enough to accomodate a new IMF.
       if (allocated(energyInputTable)) then
          imfCount=size(energyInputTable,dim=3)
          call Move_Alloc(energyInputTable,energyInputTableTemporary)
          call Alloc_Array(energyInputTable,energyInputTableAgeCount,energyInputTableMetallicityCount,imfCount,'energyInputTable')
          energyInputTable(:,:,1:imfCount)=energyInputTableTemporary
          call Dealloc_Array(energyInputTableTemporary)
       else
          call Alloc_Array(energyInputTableAge        ,energyInputTableAgeCount        ,'energyInputTableAge'        )
          call Alloc_Array(energyInputTableMetallicity,energyInputTableMetallicityCount,'energyInputTableMetallicity')
          energyInputTableAge=Make_Range(energyInputTableAgeMinimum&
               &,energyInputTableAgeMaximum,energyInputTableAgeCount,rangeType=rangeTypeLogarithmic)
          energyInputTableMetallicity(1)=0.0d0
          energyInputTableMetallicity(2:energyInputTableMetallicityCount)&
               &=Make_Range(energyInputTableMetallicityMinimum,energyInputTableMetallicityMaximum&
               &,energyInputTableMetallicityCount-1,rangeType=rangeTypeLogarithmic)
          call Alloc_Array(energyInputTable,energyInputTableAgeCount,energyInputTableMetallicityCount,1,'energyInputTable')
       end if
       
       ! Record the index in the array where this IMF will be stored.
       energyInputIndex(imfSelected)=size(energyInputTable,dim=3)

       ! Check if the table has been computed and stored previously.
       fileName='./data/Stellar_Energy_Input_'//imfNames(imfSelected)//'.xml'
       if (File_Exists(fileName)) then
          
          ! Open the XML file containing energy input.
          doc => parseFile(char(fileName),iostat=ioErr)
          if (ioErr /= 0) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','Unable to parse energy input file')
          
          ! Find the ages element and extract data.
          columnList => getElementsByTagname(doc               ,"ages")
          dataList   => getElementsByTagname(item(columnList,0),"data")
          if (getLength(dataList) /= energyInputTableAgeCount) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','ages array in XML file does not match internal expectation')
          do iAge=1,energyInputTableAgeCount
             call extractDataContent(item(dataList,iAge-1),lifetime)
             if (energyInputIndex(imfSelected) == 1) then
                energyInputTableAge(iAge)=lifetime
             else
                if (energyInputTableAge(iAge) /= lifetime) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','mismatch in ages array in XML file')
             end if
          end do

          ! Find the metallicities element and extract data.
          columnList => getElementsByTagname(doc               ,"metallicities")
          dataList   => getElementsByTagname(item(columnList,0),"data"         )
          if (getLength(dataList) /= energyInputTableMetallicityCount) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','metallicities array in XML file does not match internal expectation')
          do iMetallicity=1,energyInputTableMetallicityCount
             call extractDataContent(item(dataList,iMetallicity-1),metallicity)
             if (energyInputIndex(imfSelected) == 1) then
                energyInputTableMetallicity(iMetallicity)=metallicity
             else
                if (energyInputTableMetallicity(iMetallicity) /= metallicity) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','mismatch in metallicities array in XML file')
             end if
          end do

          ! Find the energyInput element and extract data.
          columnList => getElementsByTagname(doc               ,"energyInput")
          dataList   => getElementsByTagname(item(columnList,0),"data"            )
          if (getLength(dataList) /= energyInputTableAgeCount*energyInputTableMetallicityCount) call Galacticus_Error_Report('IMF_Energy_Input_Rate_NonInstantaneous','energy input array in XML file does not match internal expectation')
          iEnergyInput=0
          do iAge=1,energyInputTableAgeCount
             do iMetallicity=1,energyInputTableMetallicityCount
                iEnergyInput=iEnergyInput+1
                call extractDataContent(item(dataList,iEnergyInput-1),energyInputTable(iAge,iMetallicity &
                     &,energyInputIndex(imfSelected)))
             end do
          end do

       else

          call Galacticus_Display_Indent('Tabulating cumulative energy input for '//char(imfNames(imfSelected))//' IMF',2)
          
          ! Open an XML file to output the data to.
          call xml_OpenFile(char(fileName),energyInputDoc)
          call xml_NewElement(energyInputDoc,"stellarPopulation")
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
          do iAge=1,energyInputTableAgeCount
             lifetime=energyInputTableAge(iAge)
             write (progressMessage,'(a6,e8.2,a4)') 'age = ',lifetime,' Gyr'
             call Galacticus_Display_Message(progressMessage,3)
             do iMetallicity=1,energyInputTableMetallicityCount
                metallicity=energyInputTableMetallicity(iMetallicity)
                ! Find the minimum and maximum masses to integrate over for this IMF.
                minimumMass=IMF_Minimum_Mass(imfSelected)
                maximumMass=IMF_Maximum_Mass(imfSelected)
                ! Integrate cumulative energy input over the IMF between these limits.
                energyInputTable(iAge,iMetallicity,energyInputIndex(imfSelected))=Integrate(minimumMass,maximumMass&
                     &,Cumulative_Energy_Integrand,parameterPointer,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0&
                     &,toleranceRelative=1.0d-3)
                call Integrate_Done(integrandFunction,integrationWorkspace)
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
          call Galacticus_Display_Unindent('finished',2)
          call xml_EndElement(energyInputDoc,"stellarPopulation")
          call xml_Close(energyInputDoc)
       end if
       
       ! Flag that this IMF has now been tabulated.
       energyInputTabulated(imfSelected)=.true.
    end if

    ! Get the index where this IMF is stored in the table.
    tableIndex=energyInputIndex(imfSelected)
    
    ! Interpolate to get the derivative in the recycled rate at two adjacent metallicities.
    metallicity=Abundances_Get_Metallicity(fuelAbundances)
    if (metallicity > energyInputTableMetallicityMaximum) then
       metallicityIndex=energyInputTableMetallicityCount
       metallicityFactors=[1.0d0,0.0d0]
       if (present(ageMaximum)) then
          ! Get average recycling rate between ageMinimum and ageMaximum.
          energyInputRate(1)=(Interpolate(energyInputTableAgeCount,energyInputTableAge,energyInputTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMaximum,reset&
               &=interpolationAgeReset)-Interpolate(energyInputTableAgeCount,energyInputTableAge ,energyInputTable(: &
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum,reset&
               &=interpolationAgeReset))/(ageMaximum-ageMinimum)
       else
          ! Get instantaneous energy input rate at ageMinimum.
          energyInputRate(1)=Interpolate_Derivative(energyInputTableAgeCount,energyInputTableAge,energyInputTable(:&
               &,metallicityIndex,tableIndex),interpolationAgeObject,interpolationAgeAccelerator,ageMinimum,reset&
               &=interpolationAgeReset)
       end if
       energyInputRate(2)=0.0d0
    else
       metallicityIndex=Interpolate_Locate(energyInputTableMetallicityCount,energyInputTableMetallicity&
            &,interpolationMetallicityAccelerator,metallicity,reset=interpolationMetallicityReset)
       metallicityFactors=Interpolate_Linear_Generate_Factors(energyInputTableMetallicityCount,energyInputTableMetallicity&
            &,metallicityIndex,metallicity)
       ! Interpolate in age at both metallicities.
       do iMetallicity=0,1
          if (present(ageMaximum)) then
             ! Get average recycling rate between ageMinimum and ageMaximum.
             energyInputRate(iMetallicity+1)=(Interpolate(energyInputTableAgeCount,energyInputTableAge ,energyInputTable(: &
                  &,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject ,interpolationAgeAccelerator,ageMaximum&
                  &,reset =interpolationAgeReset) -Interpolate(energyInputTableAgeCount,energyInputTableAge ,energyInputTable(:&
                  &,metallicityIndex +iMetallicity,tableIndex),interpolationAgeObject,interpolationAgeAccelerator ,ageMinimum&
                  &,reset =interpolationAgeReset))/(ageMaximum-ageMinimum)
          else
             ! Get instantaneous recycling rate at ageMinimum.
             energyInputRate(iMetallicity+1)=Interpolate_Derivative(energyInputTableAgeCount,energyInputTableAge &
                  &,energyInputTable(: ,metallicityIndex+iMetallicity,tableIndex),interpolationAgeObject &
                  &,interpolationAgeAccelerator,ageMinimum,reset=interpolationAgeReset)
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
    real(c_double)        :: Cumulative_Energy_Integrand
    real(c_double), value :: initialMass
    type(c_ptr),    value :: parameterPointer

    Cumulative_Energy_Integrand=IMF_Phi(initialMass,imfSelectedGlobal)*Stellar_Feedback_Cumulative_Energy_Input(initialMass&
         &,lifetime,metallicity)
    return
  end function Cumulative_Energy_Integrand
  
end module Star_Formation_IMF
