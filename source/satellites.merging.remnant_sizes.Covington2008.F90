!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.

  use :: Kind_Numbers                           , only : kind_int8
  use :: Dark_Matter_Halo_Scales                , only : darkMatterHaloScaleClass
  use :: Satellite_Merging_Progenitor_Properties, only : mergerProgenitorPropertiesClass

  !# <mergerRemnantSize name="mergerRemnantSizeCovington2008">
  !#  <description>A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.</description>
  !# </mergerRemnantSize>
  type, extends(mergerRemnantSizeClass) :: mergerRemnantSizeCovington2008
     !% A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.
     private
     class           (darkMatterHaloScaleClass       ), pointer :: darkMatterHaloScale_        => null()
     class           (mergerProgenitorPropertiesClass), pointer :: mergerProgenitorProperties_ => null()
     double precision                                           :: energyOrbital                        , efficiencyRadiative
     logical                                                    :: warningIssued
     integer         (kind=kind_int8                 )          :: lastUniqueID
     logical                                                    :: propertiesCalculated
     double precision                                           :: radius                               ,velocityCircular    , &
          &                                                        angularMomentumSpecific
   contains
     final     ::             covington2008Destructor
     procedure :: autoHook => covington2008AutoHook
     procedure :: get      => covington2008Get
  end type mergerRemnantSizeCovington2008

  interface mergerRemnantSizeCovington2008
     !% Constructors for the {\normalfont \ttfamily covington2008} merger remnant size class.
     module procedure covington2008ConstructorParameters
     module procedure covington2008ConstructorInternal
  end interface mergerRemnantSizeCovington2008

contains

  function covington2008ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily covington2008} merger remnant size class which takes a parameter list as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerRemnantSizeCovington2008 )                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    class           (mergerProgenitorPropertiesClass), pointer       :: mergerProgenitorProperties_
    double precision                                                 :: energyOrbital              , efficiencyRadiative

    !# <inputParameter>
    !#   <name>energyOrbital</name>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The orbital energy in units of the characteristic orbital energy.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>efficiencyRadiative</name>
    !#   <defaultSource>\citep{covington_predicting_2008}</defaultSource>
    !#   <defaultValue>2.75d0</defaultValue>
    !#   <description>The coefficient, $C_\mathrm{rad}$ energy used in the \cite{covington_predicting_2008} merger remnant size algorithm.</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale"        name="darkMatterHaloScale_"        source="parameters"/>
    !# <objectBuilder class="mergerProgenitorProperties" name="mergerProgenitorProperties_" source="parameters"/>
    self=mergerRemnantSizeCovington2008(energyOrbital, efficiencyRadiative,darkMatterHaloScale_,mergerProgenitorProperties_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"       />
    !# <objectDestructor name="mergerProgenitorProperties_"/>
    return
  end function covington2008ConstructorParameters

  function covington2008ConstructorInternal(energyOrbital, efficiencyRadiative,darkMatterHaloScale_,mergerProgenitorProperties_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily covington2008} merger remnant size class.
    implicit none
    type            (mergerRemnantSizeCovington2008 )                        :: self
    double precision                                 , intent(in   )         :: energyOrbital              , efficiencyRadiative
    class           (darkMatterHaloScaleClass       ), intent(in   ), target :: darkMatterHaloScale_
    class           (mergerProgenitorPropertiesClass), intent(in   ), target :: mergerProgenitorProperties_
    !# <constructorAssign variables="energyOrbital, efficiencyRadiative, *darkMatterHaloScale_, *mergerProgenitorProperties_"/>

    self%warningIssued=.false.
    return
  end function covington2008ConstructorInternal

  subroutine covington2008AutoHook(self)
    !% Attach to the calculation reset event.
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent, openMPThreadBindingAllLevels
    implicit none
    class(mergerRemnantSizeCovington2008), intent(inout) :: self

    call calculationResetEvent%attach(self,covington2008CalculationReset,openMPThreadBindingAllLevels                                                  )
    call satelliteMergerEvent %attach(self,covington2008GetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:remnantSizeCovington2008')
    return
  end subroutine covington2008AutoHook

  subroutine covington2008Destructor(self)
    !% Destructor for the {\normalfont \ttfamily covington2008} merger remnant size class.
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    type(mergerRemnantSizeCovington2008), intent(inout) :: self
    
    !# <objectDestructor name="self%darkMatterHaloScale_"       />
    !# <objectDestructor name="self%mergerProgenitorProperties_"/>
    call calculationResetEvent%detach(self,covington2008CalculationReset)
    call satelliteMergerEvent %detach(self,covington2008GetHook         )
    return
  end subroutine covington2008Destructor

  subroutine covington2008CalculationReset(self,node)
    !% Reset the dark matter profile calculation.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class(*       ), intent(inout) :: self
    type (treeNode), intent(inout) :: node

    select type (self)
    class is (mergerRemnantSizeCovington2008)
       self%propertiesCalculated=.false.
       self%lastUniqueID       =node%uniqueID()
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine covington2008CalculationReset

  subroutine covington2008GetHook(self,node)
    !% Hookable wrapper around the get function.
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (*       ), intent(inout)         :: self
    type            (treeNode), intent(inout), target :: node
    double precision                                  :: radius                 , velocityCircular, &
         &                                               angularMomentumSpecific

    select type (self)
    type is (mergerRemnantSizeCovington2008)
       call self%get(node,radius,velocityCircular,angularMomentumSpecific)
    class default
       call Galacticus_Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine covington2008GetHook
  
  subroutine covington2008Get(self,node,radius,velocityCircular,angularMomentumSpecific)
    !% Compute the size of the merger remnant for {\normalfont \ttfamily node} using the \cite{covington_predicting_2008} algorithm.
    use :: Galacticus_Display          , only : Galacticus_Display_Message     , Galacticus_Verbosity_Level, verbosityWarn
    use :: Galacticus_Error            , only : Galacticus_Error_Report        , Galacticus_Warn
    use :: Numerical_Comparison        , only : Values_Agree
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus
    use :: String_Handling             , only : operator(//)
    implicit none
    class           (mergerRemnantSizeCovington2008), intent(inout) :: self
    type            (treeNode                      ), intent(inout) :: node
    double precision                                , intent(  out) :: radius                                      , velocityCircular   , &
         &                                                             angularMomentumSpecific
    type            (treeNode                      ), pointer       :: nodeHost
    double precision                                , parameter     :: formFactorEnergyBinding             =0.5d+00
    double precision                                , parameter     :: toleranceMassAbsolute               =1.0d-06
    double precision                                , parameter     :: toleranceMassRelative               =1.0d-09
    double precision                                , parameter     :: fractionAngularMomentumSpecificSmall=1.0d-12
    double precision                                                :: factorAngularMomentum                       , energyFinal        , &
         &                                                             fractionGasInitial                          , massHost           , &
         &                                                             radiusHost                                  , massSpheroidHost   , &
         &                                                             massSpheroidHostPreMerger                   , energyProgenitors  , &
         &                                                             energyRadiated                              , radiusVirial       , &
         &                                                             massGasSpheroidRemnant                      , massSpheroidRemnant, &
         &                                                             massSatellite                               , radiusSatellite    , &
         &                                                             massSpheroidSatellite                       , velocityVirial
    character       (len= 5                        )                :: joinString
    character       (len=70                        )                :: dataString
    type            (varying_string                )                :: message
    logical                                                         :: errorCondition

    ! The calculation of remnant size is computed when first needed and then stored. This ensures that the results are determined
    ! by the properties of the merge target prior to any modification that will occur as node components are modified in response
    ! to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call covington2008CalculationReset(self,node)
    if (.not.self%propertiesCalculated) then
       self%propertiesCalculated=.true.
       nodeHost => node%mergesWith()
       call self%mergerProgenitorProperties_%get(                           &
            &                                    node                     , &
            &                                    nodeHost                 , &
            &                                    massSatellite            , &
            &                                    massHost                 , &
            &                                    massSpheroidSatellite    , &
            &                                    massSpheroidHost         , &
            &                                    massSpheroidHostPreMerger, &
            &                                    radiusSatellite          , &
            &                                    radiusHost               , &
            &                                    factorAngularMomentum    , &
            &                                    massSpheroidRemnant      , &
            &                                    massGasSpheroidRemnant     &
            &                                   )
       if (massSatellite <= 0.0d0 .and. Values_Agree(massSpheroidHost,massSpheroidHostPreMerger,relTol=toleranceMassRelative)) then
          self%radius                 =remnantNoChange
          self%velocityCircular       =remnantNoChange
          self%angularMomentumSpecific=remnantNoChange
       else
          ! Check that the properties of the galaxies are physically reasonable.
          errorCondition=.false.
          if     (                                                  &
               &   (                                                &
               &     radiusSatellite       <= +0.0d0                &
               &    .and.                                           &
               &     massSpheroidSatellite >  +0.0d0                &
               &   )                                                &
               &  .or.                                              &
               &     massSatellite         < -toleranceMassAbsolute &
               &  .or.                                              &
               &     massSpheroidSatellite < -toleranceMassAbsolute &
               & ) then
             write (dataString,'(3(e12.6,":",e12.6,":",e12.6))') radiusSatellite,massSatellite,massSpheroidSatellite
             message=var_str('Satellite galaxy [')//node%index()//'] has '
             joinString=""
             if (radiusSatellite       <= +0.0d0        ) then
                message=message//trim(joinString)//'non-positive radius'
                joinString=", "
             end if
             if (massSatellite         <  -toleranceMassAbsolute) then
                message=message//trim(joinString)//'negative mass'
                joinString=", "
             end if
             if (massSpheroidSatellite <  -toleranceMassAbsolute) then
                message=message//trim(joinString)//'negative spheroid mass'
                joinString=", "
             end if
             message=message//' (radius:mass:spheroidMass='//trim(dataString)//')'
             call Galacticus_Display_Message(message)
             errorCondition=.true.
          end if
          if     (                                             &
               &   (                                           &
               &     radiusHost       <= +0.0d0                &
               &    .and.                                      &
               &     massSpheroidHost >  +0.0d0                &
               &   )                                           &
               &  .or.                                         &
               &     massHost         < -toleranceMassAbsolute &
               &  .or.                                         &
               &     massSpheroidHost < -toleranceMassAbsolute &
               & ) then
             write (dataString,'(3(e12.6,":",e12.6,":",e12.6))') radiusHost,massHost,massSpheroidHost
             message=var_str('Host galaxy [')//nodeHost%index()//'] has '
             joinString=""
             if (radiusHost       <= +0.0d0        ) then
                message=message//trim(joinString)//'non-positive radius'
                joinString=", "
             end if
             if (massHost         <  -toleranceMassAbsolute) then
                message=message//trim(joinString)//'negative mass'
                joinString=", "
             end if
             if (massSpheroidHost <  -toleranceMassAbsolute) then
                message=message//trim(joinString)//'negative spheroid mass'
                joinString=", "
             end if
             message=message//' (radius:mass:spheroidMass='//trim(dataString)//')'
             call Galacticus_Display_Message(message)
             errorCondition=.true.
          end if
          if (errorCondition) then
             call node    %serializeASCII()
             call nodeHost%serializeASCII()
             call Galacticus_Error_Report('error condition detected'//{introspection:location})
          end if
          ! Apply the Covington et al. (2008) algorithm to compute the size of the new remnant.
          ! Check that remnant has finite mass.
          if (massSpheroidSatellite+massSpheroidHost > 0.0d0) then
             ! First calculate the energy of the progenitors.
             energyProgenitors=0.0d0
             if (+radiusHost                 > 0.0d0)                                                                           &
                  & energyProgenitors=energyProgenitors+                      massSpheroidHost**2/                  radiusHost
             if (            radiusSatellite > 0.0d0)                                                                           &
                  & energyProgenitors=energyProgenitors+massSpheroidSatellite                 **2/  radiusSatellite
             if (+radiusHost+radiusSatellite > 0.0d0)                                                                           &
                  & energyProgenitors=energyProgenitors+massSpheroidSatellite*massSpheroidHost   /(+radiusSatellite+radiusHost) &
                  &                                    *self%energyOrbital/formFactorEnergyBinding
             ! Compute the gas fraction in the remnant.
             fractionGasInitial=massGasSpheroidRemnant/massSpheroidRemnant
             ! Compute the energy lost through radiation.
             energyRadiated=self%efficiencyRadiative*fractionGasInitial*energyProgenitors
             ! Compute the final energy.
             energyFinal=energyProgenitors+energyRadiated
             if (energyFinal <= 0.0d0) then
                write (dataString,'(e12.6,":",e12.6)') energyProgenitors,energyRadiated
                message='remnant becomes unbound (energyProgenitors:energyRadiated='//trim(dataString)//')'
                call Galacticus_Error_Report(message//{introspection:location})
             end if
             ! Compute the remnant radius.
             self%radius=(massSpheroidSatellite+massSpheroidHost)**2/(energyProgenitors+energyRadiated)
             ! Also compute the specific angular momentum at the half-mass radius.
             self%velocityCircular=sqrt(gravitationalConstantGalacticus*(massSpheroidSatellite+massSpheroidHost)/self%radius)
             self%angularMomentumSpecific=self%radius*self%velocityCircular*factorAngularMomentum
             ! Check that the specific angular momentum is reasonable.
             if (.not.self%warningIssued.and.Galacticus_Verbosity_Level() >= verbosityWarn) then
                radiusVirial  =self%darkMatterHaloScale_%virialRadius  (nodeHost)
                velocityVirial=self%darkMatterHaloScale_%virialVelocity(nodeHost)
                if (angularMomentumSpecific < fractionAngularMomentumSpecificSmall*radiusVirial*velocityVirial) then
                   message='WARNING: the specific angular momentum for node '
                   message=message//nodeHost%index()//' has become very small'//char(10)
                   message=message//' --> this will likely lead to a crash soon'
                   message=message//'NOTE: this can happen with the covington2008 implementation of the mergerRemnantSizeMethod class'//char(10)
                   message=message//' --> an alternative choice (e.g. the cole2000 implementation) may avoid this problem'//{introspection:location}
                   call Galacticus_Warn(message)
                   self%warningIssued=.true.
                end if
             end if
          else
             ! Remnant has zero mass - don't do anything.
             self%radius                 =remnantNoChange
             self%velocityCircular       =remnantNoChange
             self%angularMomentumSpecific=remnantNoChange
          end if
       end if
    end if
    radius                 =self%radius
    velocityCircular       =self%velocityCircular
    angularMomentumSpecific=self%angularMomentumSpecific
    return
  end subroutine covington2008Get
