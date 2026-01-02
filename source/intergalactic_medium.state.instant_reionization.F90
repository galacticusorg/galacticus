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

  !!{
  An implementation of the intergalactic medium state class for a simplistic model of instantaneous and full reionization with some other class providing the pre-reionization state.
  !!}

  !![
  <intergalacticMediumState name="intergalacticMediumStateInstantReionization">
   <description>The intergalactic medium is assumed to be instantaneously and fully reionized at a fixed redshift, and heated to a fixed temperature. Prior to that, the reionization state is provided by some other class.</description>
  </intergalacticMediumState>
  !!]
  type, extends(intergalacticMediumStateClass) :: intergalacticMediumStateInstantReionization
     !!{
     An \gls{igm} state class for a model in which the \gls{igm} is assumed to be instantaneously and fully reionized at a
     fixed redshift, and heated to a fixed temperature. Prior to that, the reionization state is provided by some other class.
     !!}
     private
     class           (intergalacticMediumStateClass), pointer :: preReionizationState  => null()
     double precision                                         :: reionizationTime               , reionizationTemperature       , &
          &                                                      presentDayTemperature          , expansionFactorReionizationLog
   contains
     final     ::                                instantReionizationDestructor
     procedure :: electronFraction            => instantReionizationElectronFraction
     procedure :: temperature                 => instantReionizationTemperature
     procedure :: neutralHydrogenFraction     => instantReionizationNeutralHydrogenFraction
     procedure :: neutralHeliumFraction       => instantReionizationNeutralHeliumFraction
     procedure :: singlyIonizedHeliumFraction => instantReionizationSinglyIonizedHeliumFraction
     procedure :: descriptor                  => instantReionizationDescriptor
  end type intergalacticMediumStateInstantReionization

  interface intergalacticMediumStateInstantReionization
     !!{
     Constructors for the instantReionization intergalactic medium state class.
     !!}
     module procedure instantReionizationIGMConstructorParameters
     module procedure instantReionizationIGMConstructorInternal
  end interface intergalacticMediumStateInstantReionization

contains

  function instantReionizationIGMConstructorParameters(parameters) result (self)
    !!{
    Constructor for the instantReionization \gls{igm} state class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (intergalacticMediumStateInstantReionization)                :: self
    type            (inputParameters                            ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass                    ), pointer       :: cosmologyFunctions_
    class           (cosmologyParametersClass                   ), pointer       :: cosmologyParameters_
    class           (intergalacticMediumStateClass              ), pointer       :: preReionizationState
    double precision                                                             :: reionizationRedshift          , reionizationTemperature           , &
         &                                                                          electronScatteringOpticalDepth, presentDayTemperature
    logical                                                                      :: haveReionizationRedshift      , haveElectronScatteringOpticalDepth

    ! Check and read parameters.
    haveReionizationRedshift          =parameters%isPresent('reionizationRedshift'          )
    haveElectronScatteringOpticalDepth=parameters%isPresent('electronScatteringOpticalDepth')
    if (haveReionizationRedshift.and.haveElectronScatteringOpticalDepth) call Error_Report('only one of [reionizationRedshift] or [electronScatteringOpticalDepth] can be provided'//{introspection:location})
    if (haveElectronScatteringOpticalDepth) then
       !![
       <inputParameter>
         <name>electronScatteringOpticalDepth</name>
         <source>parameters</source>
         <variable>electronScatteringOpticalDepth</variable>
         <description>The optical depth to reionization in the instantReionization \gls{igm} state model.</description>
       </inputParameter>
       !!]
     else
       haveReionizationRedshift=.true.
       !![
       <inputParameter>
         <name>reionizationRedshift</name>
         <source>parameters</source>
         <variable>reionizationRedshift</variable>
         <defaultValue>9.97d0</defaultValue>
         <defaultSource>(\citealt{hinshaw_nine-year_2012}; CMB$+H_0+$BAO)</defaultSource>
         <description>The redshift of reionization in the instantReionization \gls{igm} state model.</description>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>reionizationTemperature</name>
      <source>parameters</source>
      <variable>reionizationTemperature</variable>
      <defaultValue>1.0d4</defaultValue>
      <description>The post-reionization temperature (in units of Kelvin) in the instantReionization \gls{igm} state model.</description>
    </inputParameter>
    <inputParameter>
      <name>presentDayTemperature</name>
      <source>parameters</source>
      <variable>presentDayTemperature</variable>
      <defaultValue>1.0d3</defaultValue>
      <description>The present day temperature (in units of Kelvin) in the instantReionization \gls{igm} state model.</description>
    </inputParameter>
    <objectBuilder class="intergalacticMediumState" name="preReionizationState" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_" source="parameters"/>
    !!]
    ! Construct the object.
    if (haveReionizationRedshift) then
       self=intergalacticMediumStateInstantReionization(cosmologyFunctions_,cosmologyParameters_,preReionizationState,reionizationTemperature,presentDayTemperature,reionizationRedshift          =reionizationRedshift          )
    else
       self=intergalacticMediumStateInstantReionization(cosmologyFunctions_,cosmologyParameters_,preReionizationState,reionizationTemperature,presentDayTemperature,electronScatteringOpticalDepth=electronScatteringOpticalDepth)
    end if
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="preReionizationState"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="cosmologyParameters_"/>
    !!]
    return
  end function instantReionizationIGMConstructorParameters

  function instantReionizationIGMConstructorInternal(cosmologyFunctions_,cosmologyParameters_,preReionizationState,reionizationTemperature,presentDayTemperature,reionizationRedshift,electronScatteringOpticalDepth) result(self)
    !!{
    Constructor for the instantReionization \gls{igm} state class.
    !!}
    use :: Error      , only : Error_Report
    use :: Root_Finder, only : rangeExpandMultiplicative, rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    type            (intergalacticMediumStateInstantReionization)                          :: self
    class           (cosmologyFunctionsClass                    ), intent(inout), target   :: cosmologyFunctions_
    class           (cosmologyParametersClass                   ), intent(inout), target   :: cosmologyParameters_
    class           (intergalacticMediumStateClass              ), intent(in   ), target   :: preReionizationState
    double precision                                             , intent(in   )           :: reionizationTemperature            , presentDayTemperature
    double precision                                             , intent(in   ), optional :: reionizationRedshift               , electronScatteringOpticalDepth
    ! Used as an initial guess when solving for the epoch of reionization.
    double precision                                             , parameter               :: reionizationRedshiftCanonical=6.0d0
    ! Used to set reionization epoch to a very early time while solving for the actual epoch.
    double precision                                             , parameter               :: reionizationRedshiftEarly    =1.0d2
    double precision                                                                       :: timePresent                        , timeReionizationGuess
    type            (rootFinder                                 )                          :: finder
    !![
    <constructorAssign variables="reionizationTemperature, presentDayTemperature, *preReionizationState, *cosmologyFunctions_, *cosmologyParameters_"/>
    !!]

    ! If reionization is defined by its redshift, simply use that to compute the time of reionization.
    if (present(reionizationRedshift)) then
       self%reionizationTime=cosmologyFunctions_ %cosmicTime                 (                      &
            &                 cosmologyFunctions_%expansionFactorFromRedshift (                     &
            &                                                                  reionizationRedshift &
            &                                                                 )                     &
            &                                                                )
    else if (present(electronScatteringOpticalDepth)) then
       ! Assert that we can not also have a reionization redshift specified.
       if (present(reionizationRedshift)) call Error_Report('only one of [reionizationRedshift] or [electronScatteringOpticalDepth] can be provided'//{introspection:location})
       ! Validate optical depth.
       if (electronScatteringOpticalDepth < 0.0d0) call Error_Report('electron scattering optical depth must be > 0'//{introspection:location})
       ! Solve for the redshift of reionization which gives the desired optical depth.
       timePresent          =cosmologyFunctions_ %cosmicTime                 (                               &
            &                 cosmologyFunctions_%expansionFactorFromRedshift (                              &
            &                                                                  0.0d0                         &
            &                                                                 )                              &
            &                                                                )
       self%reionizationTime=cosmologyFunctions_ %cosmicTime                 (                               &
            &                 cosmologyFunctions_%expansionFactorFromRedshift (                              &
            &                                                                  reionizationRedshiftEarly     &
            &                                                                 )                              &
            &                                                                )
       timeReionizationGuess=cosmologyFunctions_ %cosmicTime                 (                               &
            &                 cosmologyFunctions_%expansionFactorFromRedshift (                              &
            &                                                                  reionizationRedshiftCanonical &
            &                                                                 )                              &
            &                                                                )
       finder=rootFinder(                                                             &
            &            rootFunction                 =opticalDepth                 , &
            &            toleranceAbsolute            =1.0d-3                       , &
            &            toleranceRelative            =1.0d-3                       , &
            &            rangeExpandDownward          =0.5d0                        , &
            &            rangeExpandUpward            =2.0d0                        , &
            &            rangeExpandType              =rangeExpandMultiplicative    , &
            &            rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative, &
            &            rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive, &
            &            rangeUpwardLimit             =timePresent                    &
            &           )
       self%reionizationTime=finder%find(rootGuess=timeReionizationGuess)
       ! Unset electron scattering table initialization state to force it to be recomputed now that we have the correct
       ! reionization epoch.
       self%electronScatteringTableInitialized=.false.
    else
       call Error_Report('one of [reionizationRedshift] or [electronScatteringOpticalDepth] must be provided'//{introspection:location})
    end if
    ! Compute the expansion factor at reionization.
    self%expansionFactorReionizationLog=log(cosmologyFunctions_%expansionFactor(self%reionizationTime))
    return

  contains

    double precision function opticalDepth(time)
      !!{
      Compute the optical depth to electron scattering assuming instant reionization at the given time. Subtract from this the
      target optical depth so that we have a function whose root gives the epoch of instantaneous reionization.
      !!}
      implicit none
      double precision, intent(in   ) :: time

      ! The assumption here is that the optical depth to reionization is the excess optical depth over that found in a
      ! non-reionized universe. Furthermore, it is assumed that the pre-reionization IGM state object provides the ionization
      ! history of the universe in the case of no reionization. Therefore, the optical depth we must compute here is the excess
      ! optical depth in the reionized case compared to the unreionized case.
      opticalDepth=+self                     %electronScatteringOpticalDepth(time,assumeFullyIonized=.true.) &
           &       -self%preReionizationState%electronScatteringOpticalDepth(time                          ) &
           &       -                          electronScatteringOpticalDepth
      return
    end function opticalDepth

  end function instantReionizationIGMConstructorInternal

  subroutine instantReionizationDestructor(self)
    !!{
    Destructor for the instant reionization intergalactic medium state class.
    !!}
    implicit none
    type(intergalacticMediumStateInstantReionization), intent(inout) :: self

    !![
    <objectDestructor name="self%preReionizationState"/>
    <objectDestructor name="self%cosmologyParameters_"/>
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine instantReionizationDestructor

  double precision function instantReionizationElectronFraction(self,time)
    !!{
    Return the electron fraction of the \gls{igm} in the instantReionization model.
    !!}
    use :: Numerical_Constants_Astronomical, only : heliumByMassPrimordial, hydrogenByMassPrimordial
    use :: Numerical_Constants_Atomic      , only : atomicMassHelium      , atomicMassHydrogen
    implicit none
    class           (intergalacticMediumStateInstantReionization), intent(inout) :: self
    double precision                                             , intent(in   ) :: time

    if (time > self%reionizationTime) then
       instantReionizationElectronFraction=+(                                                    &
            &                                +      hydrogenByMassPrimordial/atomicMassHydrogen  &
            &                                +2.0d0*  heliumByMassPrimordial/atomicMassHelium    &
            &                               )                                                    &
            &                              /(       hydrogenByMassPrimordial/atomicMassHydrogen)
    else
       instantReionizationElectronFraction=self%preReionizationState%electronFraction(time)
    end if
    return
  end function instantReionizationElectronFraction

  double precision function instantReionizationNeutralHydrogenFraction(self,time)
    !!{
    Return the neutral hydrogen fraction of the \gls{igm} in the instantReionization model.
    !!}
    implicit none
    class           (intergalacticMediumStateInstantReionization), intent(inout) :: self
    double precision                                             , intent(in   ) :: time

    if (time > self%reionizationTime) then
       instantReionizationNeutralHydrogenFraction=0.0d0
    else
       instantReionizationNeutralHydrogenFraction=self%preReionizationState%neutralHydrogenFraction(time)
    end if
    return
  end function instantReionizationNeutralHydrogenFraction

  double precision function instantReionizationNeutralHeliumFraction(self,time)
    !!{
    Return the neutral helium fraction of the \gls{igm} in the instantReionization model.
    !!}
    implicit none
    class           (intergalacticMediumStateInstantReionization), intent(inout) :: self
    double precision                                             , intent(in   ) :: time

    if (time > self%reionizationTime) then
       instantReionizationNeutralHeliumFraction=0.0d0
    else
       instantReionizationNeutralHeliumFraction=self%preReionizationState%neutralHeliumFraction(time)
    end if
    return
  end function instantReionizationNeutralHeliumFraction

  double precision function instantReionizationSinglyIonizedHeliumFraction(self,time)
    !!{
    Return the singly-ionized helium fraction of the \gls{igm} in the instantReionization model.
    !!}
    implicit none
    class           (intergalacticMediumStateInstantReionization), intent(inout) :: self
    double precision                                             , intent(in   ) :: time

    if (time > self%reionizationTime) then
       instantReionizationSinglyIonizedHeliumFraction=0.0d0
    else
       instantReionizationSinglyIonizedHeliumFraction=self%preReionizationState%singlyIonizedHeliumFraction(time)
    end if
    return
  end function instantReionizationSinglyIonizedHeliumFraction

  double precision function instantReionizationTemperature(self,time)
    !!{
    Return the temperature of the \gls{igm} in the instantReionization model.
    !!}
    implicit none
    class           (intergalacticMediumStateInstantReionization), intent(inout) :: self
    double precision                                             , intent(in   ) :: time
    double precision                                                             :: expansioNFactor

    if (time > self%reionizationTime) then
       ! Interpolate temperature in expansion factor between reionization and present day.
       expansionFactor=self%cosmologyFunctions_%expansionFactor(time)
       instantReionizationTemperature=+exp(                                                              &
            &                              +log(                             self%presentDayTemperature) &
            &                              +log(self%reionizationTemperature/self%presentDayTemperature) &
            &                              *log(     expansionFactor                                   ) &
            &                              /    self%expansionFactorReionizationLog                      &
            &                             )
    else
       instantReionizationTemperature=self%preReionizationState%temperature(time)
    end if
    return
  end function instantReionizationTemperature

  subroutine instantReionizationDescriptor(self,descriptor,includeClass,includeFileModificationTimes)
    !!{
    Return an input parameter list descriptor which could be used to recreate this object.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    class    (intergalacticMediumStateInstantReionization), intent(inout)           :: self
    type     (inputParameters                            ), intent(inout)           :: descriptor
    logical                                               , intent(in   ), optional :: includeClass  , includeFileModificationTimes
    character(len=18                                     )                          :: parameterLabel
    type     (inputParameters                            )                          :: parameters

    if (.not.present(includeClass).or.includeClass) call descriptor%addParameter('intergalacticMediumState','instantReionization')
    parameters=descriptor%subparameters('intergalacticMediumState')
    write (parameterLabel,'(e17.10)') self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%reionizationTime       ))
    call parameters%addParameter('reionizationRedshift'   ,trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)')                                                                                               self%reionizationTemperature
    call parameters%addParameter('reionizationTemperature',trim(adjustl(parameterLabel)))
    write (parameterLabel,'(e17.10)')                                                                                               self%presentDayTemperature
    call parameters%addParameter('presentDayTemperature'  ,trim(adjustl(parameterLabel)))
    call self%cosmologyFunctions_ %descriptor(parameters,includeClass,includeFileModificationTimes)
    call self%cosmologyParameters_%descriptor(parameters,includeClass,includeFileModificationTimes)
    return
  end subroutine instantReionizationDescriptor
