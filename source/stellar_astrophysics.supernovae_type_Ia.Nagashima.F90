!% Contains a module which implements calculations related to Type Ia supernovae.

module Supernovae_Type_Ia_Nagashima
  !% Implements calculations related to Type Ia supernovae.
  private
  public :: Supernovae_Type_Ia_Nagashima_Initialize
  
  ! Parameters of the distribution of binaries from Nagashima et al. (2005; MNRAS; 358; 1427; eqn. 17).
  double precision :: binaryMassMinimum  =3.00d0, binaryMassMaximum=12.0d0
  double precision :: typeIaNormalization=0.07d0, gamma            = 2.0d0

  ! Total yield of metals from Type Ia supernova.
  double precision :: totalYield

contains

  !# <supernovaeIaMethod>
  !#  <unitName>Supernovae_Type_Ia_Nagashima_Initialize</unitName>
  !# </supernovaeIaMethod>
  subroutine Supernovae_Type_Ia_Nagashima_Initialize(supernovaeIaMethod,SNeIa_Cumulative_Number_Get,SNeIa_Cumulative_Yield_Get)
    !% Initialize the ``Nagashima'' Type Ia supernovae module.
    use Numerical_Constants_Units
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use ISO_Varying_String
    use Galacticus_Error
    use FoX_dom
    implicit none
    type(varying_string), intent(in) :: supernovaeIaMethod
    procedure(),          pointer    :: SNeIa_Cumulative_Number_Get,SNeIa_Cumulative_Yield_Get
    type(Node),           pointer    :: doc,thisIsotope,thisYield
    type(NodeList),       pointer    :: isotopesList,propertyList
    integer                          :: iIsotope,ioErr
    double precision                 :: isotopeYield

    if (supernovaeIaMethod == 'Nagashima') then
       ! Set up pointers to our procedures.
       SNeIa_Cumulative_Number_Get => SNeIa_Cumulative_Number_Nagashima
       SNeIa_Cumulative_Yield_Get  => SNeIa_Cumulative_Yield_Nagashima
       ! Read in Type Ia yields.
       !$omp critical (FoX_DOM_Access)
       ! Open the XML file containing yields.
       doc => parseFile('data/Supernovae_Type_Ia_Yields.xml',iostat=ioErr)
       if (ioErr /= 0) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize','Unable to parse yields file')

       ! Get a list of all isotopes.
       isotopesList => getElementsByTagname(doc,"isotope")

       ! Loop through isotopes and compute the net metal yield.
       do iIsotope=0,getLength(isotopesList)-1
          thisIsotope  => item(isotopesList,iIsotope)
          propertyList => getElementsByTagname(thisIsotope,"yield")
          if (getLength(propertyList) /= 1) call Galacticus_Error_Report('Supernovae_Type_Ia_Nagashima_Initialize','isotope must have precisely one yield')
          thisYield => item(propertyList,0)
          call extractDataContent(thisYield,isotopeYield)
          totalYield=totalYield+isotopeYield
       end do

       ! Destroy the document.
       call destroy(doc)

       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end subroutine Supernovae_Type_Ia_Nagashima_Initialize

  double precision function SNeIa_Cumulative_Number_Nagashima(initialMass,age,metallicity)
    !% Compute the cumulative number of Type Ia supernovae originating per unit mass of stars that form with given {\tt
    !% initialMass} and {\tt metallicity} after a time {\tt age}. The calculation is based on that of \cite{nagashima_metal_2005}. The
    !% number returned here assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an initial
    !% mass function.
    use Stellar_Astrophysics
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity
    double precision             :: dyingStarMass,muMinimum

    ! Check if initial mass is within the range of binary masses that lead to Type Ia supernovae.
    if (initialMass > binaryMassMinimum .and. initialMass < binaryMassMaximum) then
       
       ! Get the initial mass of a star which is just dying at this age.
       dyingStarMass=Star_Initial_Mass(age,metallicity)
       
       ! Compute the cumulative number of Type Ia supernovae originating from stars of this mass.
       muMinimum=max(dyingStarMass/initialMass,(1.0d0-binaryMassMaximum/2.0d0/initialMass))
       if (muMinimum < 0.5d0) then
          SNeIa_Cumulative_Number_Nagashima=typeIaNormalization*(1.0d0-(2.0d0*muMinimum)**(1.0d0+gamma))
       else
          SNeIa_Cumulative_Number_Nagashima=0.0d0      
       end if
      
    else
       ! Mass is not in range - assume that no Type Ia SNe are produced.
       SNeIa_Cumulative_Number_Nagashima=0.0d0
    end if
    return
  end function SNeIa_Cumulative_Number_Nagashima

  double precision function SNeIa_Cumulative_Yield_Nagashima(initialMass,age,metallicity)
    !% Compute the cumulative yield from Type Ia supernovae originating per unit mass of stars that form with given {\tt
    !% initialMass} and {\tt metallicity} after a time {\tt age}. The calculation is based on the Type Ia rate calculation of
    !% \cite{nagashima_metal_2005} and the Type Ia yields from \cite{nomoto_nucleosynthesis_1997}. The number returned here
    !% assumes a distribution of binary mass ratios and so only makes sense once it is integrated over an initial mass function.
    use Stellar_Astrophysics
    implicit none
    double precision, intent(in) :: initialMass,age,metallicity

    SNeIa_Cumulative_Yield_Nagashima=SNeIa_Cumulative_Number_Nagashima(initialMass,age,metallicity)*totalYield
    return
  end function SNeIa_Cumulative_Yield_Nagashima

end module Supernovae_Type_Ia_Nagashima
