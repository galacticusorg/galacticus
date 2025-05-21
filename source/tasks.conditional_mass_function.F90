!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

  use :: Conditional_Mass_Functions    , only : conditionalMassFunction   , conditionalMassFunctionClass
  use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
  use :: Geometry_Surveys              , only : surveyGeometry            , surveyGeometryClass
  use :: Halo_Mass_Functions           , only : haloMassFunction          , haloMassFunctionClass
  use :: Mass_Function_Incompletenesses, only : massFunctionIncompleteness, massFunctionIncompletenessClass

  !![
  <task name="taskConditionalMassFunction">
   <description>A task which computes the conditional mass function in bins of mass for a fixed halo mass.</description>
  </task>
  !!]
  type, extends(taskClass) :: taskConditionalMassFunction
     !!{
     Implementation of a task which computes and outputs the halo mass function and related quantities.
     !!}
     private
     class           (cosmologyFunctionsClass        ), pointer                   :: cosmologyFunctions_         => null()
     class           (conditionalMassFunctionClass   ), pointer                   :: conditionalMassFunction_    => null()
     class           (surveyGeometryClass            ), pointer                   :: surveyGeometry_             => null()
     class           (massFunctionIncompletenessClass), pointer                   :: massFunctionIncompleteness_ => null()
     class           (haloMassFunctionClass          ), pointer                   :: haloMassFunction_           => null()
     integer                                                                      :: countMass
     double precision                                                             :: massHalo                             , massMinimum                  , &
          &                                                                          massMaximum                          , timeMinimum                  , &
          &                                                                          timeMaximum                          , massHaloMinimum              , &
          &                                                                          massHaloMaximum                      , redshiftMinimum              , &
          &                                                                          redshiftMaximum
     logical                                                                      :: useSurveyLimits                      , integrateOverHaloMassFunction
     type            (varying_string                 )                            :: outputGroupName                      , massHalo_
     double precision                                 , allocatable, dimension(:) :: massBinCenters                       , massLogarithmDelta
   contains
     final     ::            conditionalMassFunctionDestructor
     procedure :: perform => conditionalMassFunctionPerform
  end type taskConditionalMassFunction

  interface taskConditionalMassFunction
     !!{
     Constructors for the \refClass{taskConditionalMassFunction} task.
     !!}
     module procedure conditionalMassFunctionConstructorParameters
     module procedure conditionalMassFunctionConstructorInternal
  end interface taskConditionalMassFunction

contains

  function conditionalMassFunctionConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{taskConditionalMassFunction} task class which takes a parameter set as input.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (taskConditionalMassFunction    )                              :: self
    type            (inputParameters                ), intent(inout)               :: parameters
    class           (cosmologyFunctionsClass        ), pointer                     :: cosmologyFunctions_
    class           (conditionalMassFunctionClass   ), pointer                     :: conditionalMassFunction_
    class           (surveyGeometryClass            ), pointer                     :: surveyGeometry_
    class           (massFunctionIncompletenessClass), pointer                     :: massFunctionIncompleteness_
    class           (haloMassFunctionClass          ), pointer                     :: haloMassFunction_
    double precision                                 , allocatable  , dimension(:) :: massBinCenters             , massLogarithmDelta
    integer                                                                        :: countMass
    double precision                                                               :: massHalo                   , massMinimum                  , &
         &                                                                            massMaximum                , redshiftMinimum              , &
         &                                                                            redshiftMaximum            , massHaloMinimum              , &
         &                                                                            massHaloMaximum            , timeMinimum                  , &
         &                                                                            timeMaximum
    logical                                                                        :: useSurveyLimits            , integrateOverHaloMassFunction
    type            (varying_string                 )                              :: outputGroupName            , massHalo_
    character       (len=32                         )                              :: text

    !![
    <inputParameter>
      <name>outputGroupName</name>
      <description>The name of the file to which the computed conditional mass function should be output.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMinimum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The minimum redshift for which to compute the conditional mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>redshiftMaximum</name>
      <defaultValue>0.0d0</defaultValue>
      <description>The maximum redshift for which to compute the conditional mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>useSurveyLimits</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether the limiting redshifts for integrating over the halo mass function should be limited by those of a galaxy survey.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    if (parameters%isPresent('massBinCenters').or.parameters%isPresent('massLogarithmDelta')) then
       if     (                                                 &
            &   .not.parameters%isPresent('massBinCenters'    ) &
            &  .or.                                             &
            &   .not.parameters%isPresent('massLogarithmDelta') &
            & ) call Error_Report('both [massBinCenters] and [massLogarithmDelta] must be specified if either is specified'   //{introspection:location})
       if     (                                                 &
            &        parameters%isPresent('massMinimum'       ) &
            &  .or.                                             &
            &        parameters%isPresent('massMaximum'       ) &
            &  .or.                                             &
            &        parameters%isPresent('countMass'         ) &
            & ) call Error_Report('ambigous mass specification'                                                               //{introspection:location})
       !![
       <inputParameter>
         <name>massBinCenters</name>
         <description>Logarithmic mass bins centers for conditional mass function calculations.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>massLogarithmDelta</name>
         <description>Logarithmic widths of mass bins for conditional mass function calculations.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    else
       if     (                                                 &
            &        parameters%isPresent('massBinCenters'    ) &
            &  .or.                                             &
            &        parameters%isPresent('massLogarithmDelta') &
            & ) call Error_Report('ambigous mass specification'                                                               //{introspection:location})
       if     (                                                 &
            &   .not.parameters%isPresent('massMinimum'       ) &
            &  .or.                                             &
            &   .not.parameters%isPresent('massMaximum'       ) &
            &  .or.                                             &
            &   .not.parameters%isPresent('countMass'         ) &
            & ) call Error_Report('all of [massMinimum], [massMaximum], and [countMass] must be specified if any is specified'//{introspection:location})
       !![
       <inputParameter>
         <name>massMinimum</name>
         <defaultValue>1.0d8</defaultValue>
         <description>The minimum mass for which to compute the conditional mass function.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>massMaximum</name>
         <defaultValue>1.0d12</defaultValue>
         <description>The maximum mass for which to compute the conditional mass function.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>countMass</name>
         <defaultValue>21</defaultValue>
         <description>The number of bins for which to compute the conditional mass function.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    end if
    !![
    <inputParameter>
      <name>massHalo</name>
      <defaultValue>var_str('all')</defaultValue>
      <description>The halo mass for which to compute the conditional mass function. A value of ``all'' will cause the conditional mass function to be integrated over the halo mass function, giving the mass function.</description>
      <source>parameters</source>
      <variable>massHalo_</variable>
    </inputParameter>
    !!]
    integrateOverHaloMassFunction=(massHalo_ == "all")
    if (integrateOverHaloMassFunction) then
       !![
       <inputParameter>
         <name>massHaloMinimum</name>
         <defaultValue>1.0d6</defaultValue>
         <description>The minimum halo mass to use when integrating over the halo mass function.</description>
         <source>parameters</source>
       </inputParameter>
       <inputParameter>
         <name>massHaloMaximum</name>
         <defaultValue>1.0d16</defaultValue>
         <description>The maximum halo mass to use when integrating over the halo mass function.</description>
         <source>parameters</source>
       </inputParameter>
       !!]
    else
       text=char(massHalo_)
       read (text,*) massHalo
    end if
    !![
    <objectBuilder class="cosmologyFunctions"         name="cosmologyFunctions_"         source="parameters"/>
    <objectBuilder class="conditionalMassFunction"    name="conditionalMassFunction_"    source="parameters"/>
    <objectBuilder class="surveyGeometry"             name="surveyGeometry_"             source="parameters"/>
    <objectBuilder class="massFunctionIncompleteness" name="massFunctionIncompleteness_" source="parameters"/>
    <objectBuilder class="haloMassFunction"           name="haloMassFunction_"           source="parameters"/>
    !!]
    ! Compute the time corresponding to the specified redshift.
    timeMinimum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMaximum))
    timeMaximum=cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshiftMinimum))
    ! Build the object.
    !![
    <conditionalCall>
    <call>
    self=conditionalMassFunctionConstructorInternal(outputGroupName,timeMinimum,timeMaximum,useSurveyLimits,cosmologyFunctions_,conditionalMassFunction_,surveyGeometry_,massFunctionIncompleteness_,haloMassFunction_{conditions})
    </call>
    <argument  name="massMinimum"        value="massMinimum"        parameterPresent="     parameters"                   />
    <argument  name="massMaximum"        value="massMaximum"        parameterPresent="     parameters"                   />
    <argument  name="countMass"          value="countMass"          parameterPresent="     parameters"                   />
    <argument  name="massBinCenters"     value="massBinCenters"     parameterPresent="     parameters"                   />
    <argument  name="massLogarithmDelta" value="massLogarithmDelta" parameterPresent="     parameters"                   />
    <argument  name="massHalo"           value="massHalo"           condition       =".not.integrateOverHaloMassFunction"/>
    <argument  name="massHaloMinimum"    value="massHaloMinimum"    condition       ="     integrateOverHaloMassFunction"/>
    <argument  name="massHaloMaximum"    value="massHaloMaximum"    condition       ="     integrateOverHaloMassFunction"/>
    </conditionalCall>
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"        />
    <objectDestructor name="conditionalMassFunction_"   />
    <objectDestructor name="surveyGeometry_"            />
    <objectDestructor name="massFunctionIncompleteness_"/>
    <objectDestructor name="haloMassFunction_"          />
    !!]
    return
  end function conditionalMassFunctionConstructorParameters

  function conditionalMassFunctionConstructorInternal(outputGroupName,timeMinimum,timeMaximum,useSurveyLimits,cosmologyFunctions_,conditionalMassFunction_,surveyGeometry_,massFunctionIncompleteness_,haloMassFunction_,countMass,massMinimum,massMaximum,massBinCenters,massLogarithmDelta,massHalo,massHaloMinimum,massHaloMaximum) result(self)
    !!{
    Constructor for the \refClass{taskConditionalMassFunction} task class which takes a parameter set as input.
    !!}
    use :: Error            , only : Error_Report
    use :: Numerical_Ranges , only : Make_Range   , rangeTypeLogarithmic
    implicit none
    class           (cosmologyFunctionsClass        ), intent(in   ), target                 :: cosmologyFunctions_
    type            (taskConditionalMassFunction    )                                        :: self
    class           (conditionalMassFunctionClass   ), intent(in   ), target                 :: conditionalMassFunction_
    class           (surveyGeometryClass            ), intent(in   ), target                 :: surveyGeometry_
    class           (massFunctionIncompletenessClass), intent(in   ), target                 :: massFunctionIncompleteness_
    class           (haloMassFunctionClass          ), intent(in   ), target                 :: haloMassFunction_
    integer                                          , intent(in   )              , optional :: countMass
    double precision                                 , intent(in   )              , optional :: massHalo                   , massMinimum       , &
         &                                                                                      massMaximum                , massHaloMinimum   , &
         &                                                                                      massHaloMaximum
    double precision                                 , intent(in   )                         :: timeMinimum                , timeMaximum
    logical                                          , intent(in   )                         :: useSurveyLimits
    double precision                                 , intent(in   ), dimension(:), optional :: massBinCenters             , massLogarithmDelta
    type            (varying_string                 ), intent(in   )                         :: outputGroupName
    character       (len=12                         )                                        :: label

    !![
    <constructorAssign variables="outputGroupName,massMinimum,massMaximum,countMass,timeMinimum,timeMaximum,useSurveyLimits,massHalo,massHaloMinimum,massHaloMaximum,massBinCenters,massLogarithmDelta,*cosmologyFunctions_,*conditionalMassFunction_,*surveyGeometry_,*massFunctionIncompleteness_, *haloMassFunction_"/>
    !!]

    if (present(massHalo)) then
       if (     present(massHaloMinimum).or.     present(massHaloMaximum)) call Error_Report('ambiguous halo mass selection'//{introspection:location})
    else
       if (.not.present(massHaloMinimum).or..not.present(massHaloMaximum)) call Error_Report('ambiguous halo mass selection'//{introspection:location})
    end if
    self%integrateOverHaloMassFunction=.not.present(massHalo)
    if (self%integrateOverHaloMassFunction) then
       self%massHalo_="all"
    else
       write (label,'(e12.6)') massHalo
       self%massHalo_=label
    end if
    if (present(massBinCenters)) then
       self%massLogarithmDelta=log(self%massLogarithmDelta)
    else
       allocate(self%massBinCenters    (countMass))
       allocate(self%massLogarithmDelta(countMass))
       self%massBinCenters    =Make_Range(massMinimum,massMaximum,countMass,rangeType=rangeTypeLogarithmic)
       self%massLogarithmDelta=log(self%massBinCenters(2)/self%massBinCenters(1))
    end if
    self%redshiftMinimum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMaximum))
    self%redshiftMaximum=self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(timeMinimum))
    return
  end function conditionalMassFunctionConstructorInternal

  subroutine conditionalMassFunctionDestructor(self)
    !!{
    Destructor for the \refClass{taskConditionalMassFunction} task class.
    !!}
    implicit none
    type(taskConditionalMassFunction), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"        />
    <objectDestructor name="self%surveyGeometry_"            />
    <objectDestructor name="self%conditionalMassFunction_"   />
    <objectDestructor name="self%massFunctionIncompleteness_"/>
    <objectDestructor name="self%haloMassFunction_"          />
    !!]
    return
  end subroutine conditionalMassFunctionDestructor

  subroutine conditionalMassFunctionPerform(self,status)
    !!{
    Compute and output the halo mass function.
    !!}
    use :: Display              , only : displayIndent, displayUnindent
    use :: Error                , only : Error_Report , errorStatusSuccess
    use :: Output_HDF5          , only : outputFile
    use :: IO_HDF5              , only : hdf5Object
    use :: ISO_Varying_String   , only : char         , var_str           , varying_string
    use :: Numerical_Integration, only : integrator
    use :: String_Handling      , only : operator(//)
    implicit none
    class           (taskConditionalMassFunction), intent(inout), target       :: self
    integer                                      , intent(  out), optional     :: status
    double precision                             , allocatable  , dimension(:) :: conditionalMassFunction    , conditionalMassFunctionIncomplete
    integer                                                                    :: iMass                      , iField                           , &
         &                                                                        fieldCount
    double precision                                                           :: volumeIntegrand            , massFunctionIntegrand            , &
         &                                                                        massBinMinimum             , massBinMaximum                   , &
         &                                                                        binTimeMinimum             , binTimeMaximum                   , &
         &                                                                        time                       , logHaloMassLower                 , &
         &                                                                        logHaloMassUpper           , distanceMaximum
    type            (integrator                )                               :: integratorMassHalo         , integratorTime                   , &
         &                                                                        integratorNormalizationTime
    type            (hdf5Object                )                               :: outputGroup
    type            (varying_string            )                               :: message
    character       (len=12                    )                               :: label

    call displayIndent('Begin task: conditional mass function' )
    if(present(status)) status=errorStatusSuccess
    allocate(conditionalMassFunction          (self%countMass))
    allocate(conditionalMassFunctionIncomplete(self%countMass))
    ! Find logarithmic limits for halo mass in integrations.
    logHaloMassLower=log10(self%massHaloMinimum)
    logHaloMassUpper=log10(self%massHaloMaximum)
    ! Build integrators.
    integratorMassHalo         =integrator(integrandMassHalo         ,toleranceRelative=1.0d-3)
    integratorTime             =integrator(integrandTime             ,toleranceRelative=1.0d-3)
    integratorNormalizationTime=integrator(integrandNormalizationTime,toleranceRelative=1.0d-3)
    ! Compute the conditional mass function and output to file.
    do iMass=1,self%countMass
       massBinMinimum=exp(log(self%massBinCenters(iMass))-0.5d0*self%massLogarithmDelta(iMass))
       massBinMaximum=exp(log(self%massBinCenters(iMass))+0.5d0*self%massLogarithmDelta(iMass))
       ! Branch on whether the conditional mass function is to be integrated over the halo mass function.
       if (self%integrateOverHaloMassFunction) then
          ! Branch on whether a range of redshifts was given.
          if (self%timeMaximum <= self%timeMinimum) then
             ! No range of redshifts given. Compute the mass function at the minimum redshift.
             time                          =self              %timeMaximum
             conditionalMassFunction(iMass)=integratorMassHalo%integrate  (logHaloMassLower,logHaloMassUpper)
          else
             ! Determine number of fields to integrate over.
             if (self%useSurveyLimits) then
                fieldCount=self%surveyGeometry_%fieldCount()
             else
                fieldCount=1
             end if
             massFunctionIntegrand=0.0d0
             volumeIntegrand      =0.0d0
             do iField=1,fieldCount
                if (self%useSurveyLimits) then
                   ! A survey geometry is imposed. Find the maximum distance at which a galaxy of the present
                   ! mass can be detected in this survey.
                   distanceMaximum=self%surveyGeometry_%distanceMaximum(sqrt(massBinMinimum*massBinMaximum),field=iField)
                   ! Set integration limits appropriately.
                   binTimeMaximum=        self                    %timeMaximum
                   binTimeMinimum=min(                                                                      &
                        &                 self                    %timeMaximum                            , &
                        &             max(                                                                  &
                        &                 self                    %timeMinimum                            , &
                        &                 self%cosmologyFunctions_%timeAtDistanceComoving(distanceMaximum)  &
                        &                )                                                                  &
                        &            )
                else
                   ! No survey geometry is imposed, so use the full range of specified redshifts.
                   binTimeMinimum=self%timeMinimum
                   binTimeMaximum=self%timeMaximum
                end if
                ! Range of redshifts was given, integrate the mass function over this time interval.
               massFunctionIntegrand =+massFunctionIntegrand                                                                 &
                     &                +self                       %surveyGeometry_%solidAngle(iField                       ) &
                     &                *integratorTime                             %integrate (binTimeMinimum,binTimeMaximum)
                volumeIntegrand      =+volumeIntegrand                                                                       &
                     &                +self                       %surveyGeometry_%solidAngle(iField                       ) &
                     &                *integratorNormalizationTime                %integrate (binTimeMinimum,binTimeMaximum)
             end do
             if (volumeIntegrand <= 0.0d0) then
                write (label,'(e12.6)') self%massBinCenters(iMass)
                message=var_str("zero volume for conditional mass function integral in mass bin ")//iMass//" (mass = "//trim(adjustl(label))//" M☉) - check your redshift interval"
                call Error_Report(message//{introspection:location})
             end if
             conditionalMassFunction(iMass)=massFunctionIntegrand/volumeIntegrand
          end if
       else
          conditionalMassFunction(iMass)=+self%conditionalMassFunction_%massFunction(                     &
               &                                                                     self%massHalo      , &
               &                                                                          massBinMinimum  &
               &                                                                    )                     &
               &                         -self%conditionalMassFunction_%massFunction(                     &
               &                                                                     self%massHalo      , &
               &                                                                          massBinMaximum  &
               &                                                                    )
       end if
    end do
    conditionalMassFunction=+conditionalMassFunction &
         &                  /self%massLogarithmDelta
    ! Find the mass function with incompleteness.
    do iMass=1,self%countMass
       conditionalMassFunctionIncomplete(iMass)=+conditionalMassFunction                                          (iMass)  &
            &                                   *self%massFunctionIncompleteness_%completeness(self%massBinCenters(iMass))
    end do
    ! Write the data to file.
    outputGroup=outputFile%openGroup(char(self%outputGroupName))
    call outputGroup%writeDataset(self%massBinCenters,"mass" ,comment="mass in units of M☉")
    if (self%integrateOverHaloMassFunction) then
       call outputGroup%writeDataset(conditionalMassFunction          ,"massFunction"          ,comment="Mass function in units of Mpc⁻³ per log(mass)."                 )
       call outputGroup%writeDataset(conditionalMassFunctionIncomplete,"massFunctionIncomplete",comment="Incomplete mass function in units of Mpc⁻³ per log(mass)."      )
    else
       call outputGroup%writeDataset(conditionalMassFunction          ,"massFunction"          ,comment="Conditional mass function in units of per log(mass)."           )
       call outputGroup%writeDataset(conditionalMassFunctionIncomplete,"massFunctionIncomplete",comment="Incomplete conditional mass function in units of per log(mass).")
    end if
    call outputGroup%close()
    call displayUnindent('Done task: conditional mass function' )
    return

  contains

    double precision function integrandTime(timePrime)
      !!{
      Integral over time.
      !!}
      implicit none
      double precision            , intent(in   ) :: timePrime
      type            (integrator)                :: integrator_

      time         =                                                           timePrime
      integrator_  = integrator                                               (integrandMassHalo,toleranceRelative=1.0d-3          )
      integrandTime=+integrator_                    %integrate                (logHaloMassLower ,                  logHaloMassUpper) &
           &        *self       %cosmologyFunctions_%comovingVolumeElementTime(timePrime                                           )
      return
    end function integrandTime

    double precision function integrandNormalizationTime(timePrime)
      !!{
      Normalization integral over time.
      !!}
      implicit none
      double precision, intent(in   ) :: timePrime

      integrandNormalizationTime=self%cosmologyFunctions_%comovingVolumeElementTime(timePrime)
      return
    end function integrandNormalizationTime

    double precision function integrandMassHalo(logMass)
      !!{
      Integral over halo mass function.
      !!}
      implicit none
      double precision, intent(in   ) :: logMass
      double precision                :: mass

      mass             =+10.0d0**logMass
      integrandMassHalo=+     self%haloMassFunction_       %differential(time,mass               )  &
           &            *                                                     mass                  &
           &            *log(10.0d0)                                                                &
           &            *max(                                                                       &
           &                 +0.0d0                                                               , &
           &                 +self%conditionalMassFunction_%massFunction(     mass,massBinMinimum)  &
           &                 -self%conditionalMassFunction_%massFunction(     mass,massBinMaximum)  &
           &                )
      return
    end function integrandMassHalo

  end subroutine conditionalMassFunctionPerform
