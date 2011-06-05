!% Contains a module which implements calculation related to stellar astrophyics.

module Stellar_Astrophysics_File
  !% Implements calculation related to stellar astrophyics.
  private
  public :: Stellar_Astrophysics_File_Initialize

  ! Arrays to store stellar properties.
  double precision, allocatable, dimension(:) :: stellarLifetime,stellarLifetimeMass,stellarLifetimeMetallicity
  double precision, allocatable, dimension(:) :: ejectedMass    ,ejectedMassMass    ,ejectedMassMetallicity
  double precision, allocatable, dimension(:) :: metalYield     ,metalYieldMass     ,metalYieldMetallicity

contains

  !# <stellarAstrophysicsMethod>
  !#  <unitName>Stellar_Astrophysics_File_Initialize</unitName>
  !# </stellarAstrophysicsMethod>
  subroutine Stellar_Astrophysics_File_Initialize(stellarAstrophysicsMethod,Star_Ejected_Mass_Get,Star_Initial_Mass_Get,Star_Metal_Yield_Mass_Get,Star_Lifetime_Get)
    !% Initialize the stellar astrophysics module.
    use FoX_dom
    use Galacticus_Error
    use Memory_Management
    use ISO_Varying_String
    use Input_Parameters
    implicit none
    type(varying_string),          intent(in)    :: stellarAstrophysicsMethod
    procedure(),          pointer, intent(inout) :: Star_Ejected_Mass_Get,Star_Initial_Mass_Get,Star_Metal_Yield_Mass_Get,Star_Lifetime_Get
    type(Node),           pointer                :: doc,thisStar,thisDatum
    type(NodeList),       pointer                :: starList,propertyList
    type(varying_string)                         :: stellarPropertiesFile
    integer                                      :: ioErr,iStar,lifetimeCount,ejectedMassCount,metalYieldCount
    double precision                             :: initialMass,metallicity

    ! Check if our method is selected.
    if (stellarAstrophysicsMethod == 'file') then
       ! Set up procedure pointers.
       Star_Ejected_Mass_Get     => Star_Ejected_Mass_File
       Star_Initial_Mass_Get     => Star_Initial_Mass_File
       Star_Metal_Yield_Mass_Get => Star_Metal_Yield_Mass_File
       Star_Lifetime_Get         => Star_Lifetime_File
       
       ! Get the name of the file containing stellar data.
       !@ <inputParameter>
       !@   <name>stellarPropertiesFile</name>
       !@   <defaultValue>data/Stellar\_Properties\_Compilation.xml</defaultValue>
       !@   <attachedTo>module</attachedTo>
       !@   <description>
       !@     The name of the XML file from which to read stellar properties (ejected masses, yields, etc.).
       !@   </description>
       !@ </inputParameter>
       call Get_Input_Parameter('stellarPropertiesFile',stellarPropertiesFile,defaultValue='data/Stellar_Properties_Compilation.xml')

       !$omp critical (FoX_DOM_Access)
       ! Open the XML file containing stellar properties.
       doc => parseFile(char(stellarPropertiesFile),iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','Unable to parse stellar properties file')

       ! Get a list of all stars.
       starList => getElementsByTagname(doc,"star")

       ! Count up number of stars with given properties.
       lifetimeCount   =0
       ejectedMassCount=0
       metalYieldCount =0
       do iStar=0,getLength(starList)-1
          thisStar     => item(starList,iStar)
          propertyList => getElementsByTagname(thisStar,"initialMass")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','star must have precisely one initial mass')
          propertyList => getElementsByTagname(thisStar,"metallicity")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','star must have precisely one metallicity')
          propertyList => getElementsByTagname(thisStar,"lifetime")
          if (getLength(propertyList) == 1) lifetimeCount=lifetimeCount+1
          if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','star has multiple lifetimes')
          propertyList => getElementsByTagname(thisStar,"ejectedMass")
          if (getLength(propertyList) == 1) ejectedMassCount=ejectedMassCount+1
          if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','star has multiple ejected masses')
          propertyList => getElementsByTagname(thisStar,"metalYieldMass")
          if (getLength(propertyList) == 1) metalYieldCount=metalYieldCount+1
          if (getLength(propertyList) >  1) call Galacticus_Error_Report('Stellar_Astrophysics_Initialize','star has multiple metal yield masses')
       end do

       ! Allocate arrays to store stellar properties.
       call Alloc_Array(stellarLifetime           ,lifetimeCount   ,'stellarLifetime'           )
       call Alloc_Array(stellarLifetimeMass       ,lifetimeCount   ,'stellarLifetimeMass'       )
       call Alloc_Array(stellarLifetimeMetallicity,lifetimeCount   ,'stellarLifetimeMetallicity')
       call Alloc_Array(ejectedMass               ,ejectedMassCount,'ejectedMass'               )
       call Alloc_Array(ejectedMassMass           ,ejectedMassCount,'ejectedMassMass'           )
       call Alloc_Array(ejectedMassMetallicity    ,ejectedMassCount,'ejectedMassMetallicity'    )
       call Alloc_Array(metalYield                ,metalYieldCount ,'metalYield'                )
       call Alloc_Array(metalYieldMass            ,metalYieldCount ,'metalYieldMass'            )
       call Alloc_Array(metalYieldMetallicity     ,metalYieldCount ,'metalYieldMetallicity'     )

       ! Loop over stars to process their properties.
       lifetimeCount   =0
       ejectedMassCount=0
       metalYieldCount =0
       do iStar=0,getLength(starList)-1
          thisStar     => item(starList,iStar)
          propertyList => getElementsByTagname(thisStar,"initialMass")
          thisDatum    => item(propertyList,0)
          call extractDataContent(thisDatum,initialMass)
          propertyList => getElementsByTagname(thisStar,"metallicity")
          thisDatum    => item(propertyList,0)
          call extractDataContent(thisDatum,metallicity)

          ! Process stellar lifetimes.
          propertyList => getElementsByTagname(thisStar,"lifetime")
          if (getLength(propertyList) == 1) then
             thisDatum => item(propertyList,0)
             lifetimeCount=lifetimeCount+1
             call extractDataContent(thisDatum,stellarLifetime(lifetimeCount))
             stellarLifetimeMass       (lifetimeCount)=initialMass
             stellarLifetimeMetallicity(lifetimeCount)=metallicity
          end if

          ! Process ejected masses.
          propertyList => getElementsByTagname(thisStar,"ejectedMass")
          if (getLength(propertyList) == 1) then
             thisDatum => item(propertyList,0)
             ejectedMassCount=ejectedMassCount+1
             call extractDataContent(thisDatum,ejectedMass(ejectedMassCount))
             ejectedMassMass       (ejectedMassCount)=initialMass
             ejectedMassMetallicity(ejectedMassCount)=metallicity
          end if

          ! Process metal yields.
          propertyList => getElementsByTagname(thisStar,"metalYieldMass")
          if (getLength(propertyList) == 1) then
             thisDatum => item(propertyList,0)
             metalYieldCount=metalYieldCount+1
             call extractDataContent(thisDatum,metalYield(metalYieldCount))
             metalYieldMass       (metalYieldCount)=initialMass
             metalYieldMetallicity(metalYieldCount)=metallicity
          end if

       end do

       ! Destroy the document.
       call destroy(doc)
       !$omp end critical (FoX_DOM_Access)

    end if

    return
  end subroutine Stellar_Astrophysics_File_Initialize

  double precision function Star_Initial_Mass_File(lifetime,metallicity)
    !% Return the initial mass of a star of given {\tt lifetime} and {\tt metallicity}.
    use Numerical_Interpolation_2D_Irregular
    implicit none
    double precision,              intent(in) :: lifetime,metallicity
    type(interp2dIrregularObject), save       :: interpolationWorkspace
    !$omp threadprivate(interpolationWorkspace)

    Star_Initial_Mass_File=Interpolate_2D_Irregular(stellarLifetime,stellarLifetimeMetallicity,stellarLifetimeMass,lifetime,metallicity&
         &,interpolationWorkspace,reset=.true.)

    return
  end function Star_Initial_Mass_File

  double precision function Star_Lifetime_File(initialMass,metallicity)
    !% Return the lifetime of a star (in Gyr) given an {\tt initialMass} and {\tt metallicity}.
    use Numerical_Interpolation_2D_Irregular
    implicit none
    double precision,              intent(in) :: initialMass,metallicity
    type(interp2dIrregularObject), save       :: interpolationWorkspace
    logical,                       save       :: resetInterpolation=.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)

    Star_Lifetime_File=Interpolate_2D_Irregular(stellarLifetimeMass,stellarLifetimeMetallicity,stellarLifetime,initialMass,metallicity&
         & ,interpolationWorkspace,reset=resetInterpolation)

    return
  end function Star_Lifetime_File

  double precision function Star_Ejected_Mass_File(initialMass,metallicity)
    !% Return the mass ejected during the lifetime of a star of given {\tt initialMass} and {\tt metallicity}.
    use Numerical_Interpolation_2D_Irregular
    implicit none
    double precision,              intent(in) :: initialMass,metallicity
    type(interp2dIrregularObject), save       :: interpolationWorkspace
    logical,                       save       :: resetInterpolation=.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)

    ! Compute the ejected mass.
    Star_Ejected_Mass_File=max(Interpolate_2D_Irregular(ejectedMassMass,ejectedMassMetallicity,ejectedMass,initialMass,metallicity&
         &,interpolationWorkspace,reset=resetInterpolation),0.0d0)

    return
  end function Star_Ejected_Mass_File

  double precision function Star_Metal_Yield_Mass_File(initialMass,metallicity)
    !% Return the mass of metals yielded by a star of given {\tt initialMass} and {\tt metallicity}.
    use Numerical_Interpolation_2D_Irregular
    implicit none
    double precision,              intent(in) :: initialMass,metallicity
    type(interp2dIrregularObject), save       :: interpolationWorkspace
    logical,                       save       :: resetInterpolation=.true.
    !$omp threadprivate(interpolationWorkspace,resetInterpolation)

    ! Compute the metal mass yield.
    Star_Metal_Yield_Mass_File=max(Interpolate_2D_Irregular(metalYieldMass,metalYieldMetallicity,metalYield,initialMass,metallicity&
         &,interpolationWorkspace,reset=resetInterpolation),0.0d0)
    return
  end function Star_Metal_Yield_Mass_File

end module Stellar_Astrophysics_File
