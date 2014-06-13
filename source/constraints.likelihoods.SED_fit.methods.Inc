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

!% Defines methods for a SED fitting likelihood function.

function likelihoodSEDFitConstructor(definition)
  !% Constructor for the SED fitting likelihood class.
  use FoX_DOM
  use IO_XML
  use Instruments_Filters
  use ISO_Varying_String
  use Memory_Management
  use Stellar_Population_Spectra_Postprocess
  use Cosmology_Functions
  implicit none
  type     (likelihoodSEDFit       )                         :: likelihoodSEDFitConstructor
  type     (node                   ), pointer, intent(in   ) :: definition
  type     (node                   ), pointer                :: thisPhotometry             , thisElement
  type     (nodeList               ), pointer                :: photometry
  class    (cosmologyFunctionsClass), pointer                :: cosmologyFunctions_
  character(len=128                )                         :: filterName                 , systemName
  integer                                                    :: i
  double precision                                           :: redshift

  ! Find the photometry for the galaxy to fit.
  photometry => getElementsByTagName(definition,"photometry")
  ! Allocate storage for photometry.
  likelihoodSEDFitConstructor%photometryCount=getLength(photometry)
  call Alloc_Array(likelihoodSEDFitConstructor%magnitude               ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%error                   ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%filterIndex             ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%luminosityIndex         ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%redshift                ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%age                     ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%massToLightRatio        ,[likelihoodSEDFitConstructor%photometryCount])
  call Alloc_Array(likelihoodSEDFitConstructor%postprocessingChainIndex,[likelihoodSEDFitConstructor%photometryCount])
  ! Extract photometry information.
  do i=1,likelihoodSEDFitConstructor%photometryCount
     thisPhotometry => item(photometry,i-1)
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"magnitude")
     call extractDataContent(thisElement,likelihoodSEDFitConstructor%magnitude(i))
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"error"    )
     call extractDataContent(thisElement,likelihoodSEDFitConstructor%error    (i))
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"filter"   )
     call extractDataContent(thisElement,filterName                              )
     thisElement    => XML_Get_First_Element_By_Tag_Name(thisPhotometry,"system"   )
     call extractDataContent(thisElement,systemName                              )
     ! Set a luminosity index.
     likelihoodSEDFitConstructor%luminosityIndex(i)=i
     ! Find the index for this filter.
     likelihoodSEDFitConstructor%filterIndex    (i)=Filter_Get_Index(var_str(trim(filterName)))
     ! Convert to AB system.
     if (trim(systemName) == "vega")                                         &
          & likelihoodSEDFitConstructor                     %magnitude  (i)  &
          &  =                   likelihoodSEDFitConstructor%magnitude  (i)  &
          &  -Filter_Vega_Offset(likelihoodSEDFitConstructor%filterIndex(i))
  end do
  ! Find stellar spectra postprocessing chain to use.
  likelihoodSEDFitConstructor%postprocessingChainIndex=Stellar_Population_Spectrum_Postprocess_Index(var_str('default'))
  ! Find the time at the observed redshift.
  thisElement => XML_Get_First_Element_By_Tag_Name(definition,"redshift")
  call extractDataContent(thisElement,redshift)
  cosmologyFunctions_ => cosmologyFunctions()
  likelihoodSEDFitConstructor%timeObserved                           &
       & =cosmologyFunctions_ %cosmicTime                 (          &
       &   cosmologyFunctions_%expansionFactorFromRedshift (         &
       &                                                    redshift &
       &                                                   )         &
       &                                                  )
  likelihoodSEDFitConstructor%redshift=redshift
  return
end function likelihoodSEDFitConstructor

double precision function likelihoodSEDFitEvaluate(self,simulationState,parameterMappings,simulationConvergence,temperature,logLikelihoodCurrent,logPriorCurrent,logPriorProposed)
  !% Return the log-likelihood for the SED fitting likelihood function.
  use, intrinsic :: ISO_C_Binding
  use FGSL
  use Constraints_State
  use Constraints_Convergence
  use Numerical_Integration
  use Galactic_Structure_Options
  use Abundances_Structure
  implicit none
  class           (likelihoodSEDFit          ), intent(inout)               :: self
  class           (state                     ), intent(in   )               :: simulationState
  type            (mappingList               ), intent(in   ), dimension(:) :: parameterMappings
  class           (convergence               ), intent(inout)               :: simulationConvergence
  double precision                            , intent(in   )               :: temperature                         , logLikelihoodCurrent, &
       &                                                                       logPriorCurrent                     , logPriorProposed
  double precision                            , parameter                   :: timeStart                     =0.0d0
  double precision                            , allocatable  , dimension(:) :: stateVector
  integer                                                                   :: iMagnitude                          , i
  double precision                                                          :: mass                                , timeScale           , &
       &                                                                       starFormationRateNormalization      , magnitude           , &
       &                                                                       luminosity                          , metallicity
  type            (abundances                )                              :: abundancesStars
  type            (fgsl_function             )                              :: integrandFunction
  type            (fgsl_integration_workspace)                              :: integrationWorkspace
  type            (c_ptr                     )                              :: parameterPointer
  logical                                                                   :: integrationReset

  ! Do not evaluate if the proposed prior is impossible.
  if (logPriorProposed <= logImpossible) then
     likelihoodSEDFitEvaluate=0.0d0
     return
  end if
  ! Get the simulation state.
  stateVector=simulationState%get()
  do i=1,size(stateVector)
     stateVector(i)=parameterMappings(i)%thisMapping%unmap(stateVector(i))
  end do
  mass       =stateVector(1)
  timeScale  =stateVector(2)
  metallicity=stateVector(3)
  ! Construct the abundances object.
  call abundancesStars%metallicitySet(metallicity,metallicityType=linearByMassSolar)
  ! Iterate over bands.
  likelihoodSEDFitEvaluate      =0.0d0
  starFormationRateNormalization=1.0d0
  do iMagnitude=0,size(self%magnitude)
     ! Integrate star formation history to get luminosity.
     integrationReset=.true.
     luminosity=Integrate(                               &
          &               timeStart                    , &
          &               self%timeObserved            , &
          &               luminosityIntegrand          , &
          &               parameterPointer             , &
          &               integrandFunction            , &
          &               integrationWorkspace         , &
          &               toleranceRelative     =1.0d-3, &
          &               reset=integrationReset         &
          &              )      
     call Integrate_Done(integrandFunction,integrationWorkspace)
     ! Determine if we're computing mass or luminosity.
     if (iMagnitude > 0) then
        ! Accumulate to likelihood.
        magnitude=-2.5d0*log10(luminosity)
        likelihoodSEDFitEvaluate=             &
             & +likelihoodSEDFitEvaluate      &
             & -0.5d0                         &
             & *(                             &
             &   (                            &
             &    +     magnitude             &
             &    -self%magnitude(iMagnitude) &
             &   )                            &
             &   / self%error    (iMagnitude) &
             &  )**2
     else
        ! Normalize the star formation rate.
        starFormationRateNormalization=mass/luminosity
     end if
  end do
  ! Destroy abundances object.
  call abundancesStars%destroy()
  return

contains
  
  function luminosityIntegrand(time,parameterPointer) bind(c)
    !% Star formation rate integrand.
    use, intrinsic :: ISO_C_Binding
    use Star_Formation_IMF
    use Stellar_Population_Luminosities
    implicit none
    real(c_double)        :: luminosityIntegrand
    real(c_double), value :: time
    type(c_ptr   ), value :: parameterPointer
    integer               :: imfIndex
    double precision      :: recycledFractionInstantaneous, starFormationRate
    
    ! Find stellar population age.
    self%age         =self%timeObserved-time
    ! Find star formation rate.
    starFormationRate=starFormationRateNormalization*exp(-(time-timeStart)/timeScale)
    ! Find the IMF.
    imfIndex         =IMF_Select(starFormationRate,abundancesStars,componentTypeDisk)
    ! Determine if we're computing luminosity or mass.
    if (iMagnitude > 0) then
       ! Find the mass-to-light ratio.
       self%massToLightRatio        =Stellar_Population_Luminosity      (                               &
            &                                                            self%luminosityIndex         , &
            &                                                            self%filterIndex             , &
            &                                                            self%postprocessingChainIndex, &
            &                                                            imfIndex                     , &
            &                                                            abundancesStars              , &
            &                                                            self%age                     , &
            &                                                            self%redshift                  &
            &                                                           )
       luminosityIntegrand=starFormationRate*self%massToLightRatio(iMagnitude)
    else
       ! Find the recycled fraction.
       recycledFractionInstantaneous=IMF_Recycled_Fraction_Instantaneous(starFormationRate,abundancesStars,componentTypeDisk)
       ! Find the contribution to final mass.
       luminosityIntegrand          =starFormationRate*(1.0d0-recycledFractionInstantaneous)
    end if
    return
  end function luminosityIntegrand
   
end function likelihoodSEDFitEvaluate

subroutine likelihoodSEDFitFunctionChanged(self)
  !% Respond to possible changes in the likelihood function.
  implicit none
  class(likelihoodSEDFit), intent(inout) :: self
  
  return
end subroutine likelihoodSEDFitFunctionChanged
