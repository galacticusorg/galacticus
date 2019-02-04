!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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
  
  use Satellite_Merging_Progenitor_Properties
  use Dark_Matter_Halo_Scales

  !# <mergerRemnantSize name="mergerRemnantSizeCovington2008">
  !#  <description>A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.</description>
  !# </mergerRemnantSize>
  type, extends(mergerRemnantSizeClass) :: mergerRemnantSizeCovington2008
     !% A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.
     private
     class           (darkMatterHaloScaleClass       ), pointer :: darkMatterHaloScale_
     class           (mergerProgenitorPropertiesClass), pointer :: mergerProgenitorProperties_
     double precision                                           :: energyOrbital              , efficiencyRadiative
     logical                                                    :: warningIssued
   contains
     final     ::        covington2008Destructor
     procedure :: get => covington2008Get
  end type mergerRemnantSizeCovington2008

  interface mergerRemnantSizeCovington2008
     !% Constructors for the {\normalfont \ttfamily covington2008} merger remnant size class.
     module procedure covington2008ConstructorParameters
     module procedure covington2008ConstructorInternal
  end interface mergerRemnantSizeCovington2008

contains

  function covington2008ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily covington2008} merger remnant size class which takes a parameter list as input.
    use Input_Parameters
    implicit none
    type            (mergerRemnantSizeCovington2008 )                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass       ), pointer       :: darkMatterHaloScale_
    class           (mergerProgenitorPropertiesClass), pointer       :: mergerProgenitorProperties_
    double precision                                                 :: energyOrbital              , efficiencyRadiative
    
    !# <inputParameter>
    !#   <name>energyOrbital</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The orbital energy in units of the characteristic orbital energy.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>efficiencyRadiative</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultSource>\citep{covington_predicting_2008}</defaultSource>
    !#   <defaultValue>2.75d0</defaultValue>
    !#   <description>The coefficient, $C_\mathrm{rad}$ energy used in the \cite{covington_predicting_2008} merger remnant size algorithm.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
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

  subroutine covington2008Destructor(self)
    !% Destructor for the {\normalfont \ttfamily covington2008} merger remnant size class.
    implicit none
    type(mergerRemnantSizeCovington2008), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"       />
    !# <objectDestructor name="self%mergerProgenitorProperties_"/>
    return
  end subroutine covington2008Destructor

  subroutine covington2008Get(self,node,radius,velocityCircular,angularMomentumSpecific)
    !% Compute the size of the merger remnant for {\normalfont \ttfamily node} using the \cite{covington_predicting_2008} algorithm.
    use Numerical_Constants_Physical
    use Numerical_Comparison
    use Galacticus_Error
    use String_Handling
    use Galacticus_Display
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
       radius                 =remnantNoChange
       velocityCircular       =remnantNoChange
       angularMomentumSpecific=remnantNoChange
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
          radius=(massSpheroidSatellite+massSpheroidHost)**2/(energyProgenitors+energyRadiated)
          ! Also compute the specific angular momentum at the half-mass radius.
          velocityCircular=sqrt(gravitationalConstantGalacticus*(massSpheroidSatellite+massSpheroidHost)/radius)
          angularMomentumSpecific=radius*velocityCircular*factorAngularMomentum
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
          radius                 =remnantNoChange
          velocityCircular       =remnantNoChange
          angularMomentumSpecific=remnantNoChange
       end if
    end if
    return
  end subroutine covington2008Get
