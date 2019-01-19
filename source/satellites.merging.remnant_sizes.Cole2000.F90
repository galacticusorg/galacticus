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

  !# <mergerRemnantSize name="mergerRemnantSizeCole2000">
  !#  <description>A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.</description>
  !# </mergerRemnantSize>
  type, extends(mergerRemnantSizeClass) :: mergerRemnantSizeCole2000
     !% A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.
     private
     class           (mergerProgenitorPropertiesClass), pointer :: mergerProgenitorProperties_
     double precision                                           :: energyOrbital
   contains
     final     ::        cole2000Destructor
     procedure :: get => cole2000Get
  end type mergerRemnantSizeCole2000

  interface mergerRemnantSizeCole2000
     !% Constructors for the {\normalfont \ttfamily cole2000} merger remnant size class.
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerRemnantSizeCole2000

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily cole2000} merger remnant size class which takes a parameter list as input.
    use Input_Parameters
    implicit none
    type            (mergerRemnantSizeCole2000      )                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (mergerProgenitorPropertiesClass), pointer       :: mergerProgenitorProperties_
    double precision                                                 :: energyOrbital

    !# <inputParameter>
    !#   <name>energyOrbital</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>1.0d0</defaultValue>
    !#   <description>The orbital energy used in the ``cole2000'' merger remnant sizes calculation in units of the characteristic orbital energy.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="mergerProgenitorProperties" name="mergerProgenitorProperties_" source="parameters"/>
    self=mergerRemnantSizeCole2000(energyOrbital,mergerProgenitorProperties_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(energyOrbital,mergerProgenitorProperties_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily cole2000} merger remnant size class.
    implicit none
    type            (mergerRemnantSizeCole2000      )                        :: self
    double precision                                 , intent(in   )         :: energyOrbital
    class           (mergerProgenitorPropertiesClass), intent(in   ), target :: mergerProgenitorProperties_
    !# <constructorAssign variables="energyOrbital, *mergerProgenitorProperties_"/>

    return
  end function cole2000ConstructorInternal

  subroutine cole2000Destructor(self)
    !% Destructor for the {\normalfont \ttfamily cole2000} merger remnant size class.
    implicit none
    type(mergerRemnantSizeCole2000), intent(inout) :: self

    !# <objectDestructor name="self%mergerProgenitorProperties_"/>
    return
  end subroutine cole2000Destructor

  subroutine cole2000Get(self,node,radius,velocityCircular,angularMomentumSpecific)
    !% Compute the size of the merger remnant for {\normalfont \ttfamily node} using the \cite{cole_hierarchical_2000} algorithm.
    use Numerical_Constants_Physical
    use Numerical_Comparison
    use Galacticus_Error
    use String_Handling
    use Galacticus_Display
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    implicit none
    class           (mergerRemnantSizeCole2000), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node
    double precision                           , intent(  out) :: radius                           , velocityCircular         , &
         &                                                        angularMomentumSpecific
    type            (treeNode                 ), pointer       :: nodeHost
    double precision                           , parameter     :: formFactorEnergyBinding   =0.5d+0
    double precision                           , parameter     :: toleranceMassAbsolute     =1.0d+0
    double precision                           , parameter     :: toleranceMassRelative     =1.0d-9
    double precision                                           :: factorAngularMomentum            , massDarkMatterHost       , &
         &                                                        massHost                         , radiusHost               , &
         &                                                        massSpheroidHost                 , massSpheroidHostPreMerger, &
         &                                                        massSpheroidHostTotal            , energyProgenitors        , &
         &                                                        massGasSpheroidRemnant           , massSpheroidRemnant      , &
         &                                                        massDarkMatterSatellite          , massSatellite            , &
         &                                                        radiusSatellite                  , massSpheroidSatellite    , &
         &                                                        massSpheroidTotalSatellite
    character       (len= 5                   )                :: joinString
    character       (len=40                   )                :: dataString
    type            (varying_string           )                :: message
    logical                                                    :: errorCondition

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
    if (massSpheroidSatellite <= 0.0d0 .and. Values_Agree(massSpheroidHost,massSpheroidHostPreMerger,relTol=toleranceMassRelative)) then
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
       ! Check if host has finite mass.
       if (massSpheroidSatellite+massSpheroidHost > 0.0d0) then
          ! Compute masses of dark matter within the host and satellite radii.
          massDarkMatterHost     =Galactic_Structure_Enclosed_Mass(nodeHost,radiusHost     ,massType=massTypeDark,haloLoaded=.true.)
          massDarkMatterSatellite=Galactic_Structure_Enclosed_Mass(node    ,radiusSatellite,massType=massTypeDark,haloLoaded=.true.)
          ! Combine baryonic and dark matter masses.
          massSpheroidHostTotal     =+massSpheroidHost     +2.0d0*massDarkMatterHost
          massSpheroidTotalSatellite=+massSpheroidSatellite+2.0d0*massDarkMatterSatellite
          ! Apply the Cole et al. (2000) algorithm to compute the size of the new remnant.
          energyProgenitors=0.0d0
          if (+radiusHost                 > 0.0d0)                                                                                      &
               & energyProgenitors=+energyProgenitors+                           massSpheroidHostTotal**2/                  radiusHost
          if (           +radiusSatellite > 0.0d0)                                                                                      &
               & energyProgenitors=+energyProgenitors+massSpheroidTotalSatellite                      **2/  radiusSatellite
          if (+radiusHost+radiusSatellite > 0.0d0)                                                                                      &
               & energyProgenitors=+energyProgenitors+massSpheroidTotalSatellite*massSpheroidHostTotal   /(+radiusSatellite+radiusHost) &
               &                                     *self%energyOrbital/formFactorEnergyBinding
          radius=(massSpheroidTotalSatellite+massSpheroidHostTotal)**2/energyProgenitors
          ! Also compute the specific angular momentum at the half-mass radius.
          velocityCircular       =sqrt(gravitationalConstantGalacticus*(massSpheroidSatellite+massSpheroidHost)/radius)
          angularMomentumSpecific=radius*velocityCircular*factorAngularMomentum
       else
          ! Remnant has zero mass - don't do anything.
          radius                 =remnantNoChange
          velocityCircular       =remnantNoChange
          angularMomentumSpecific=remnantNoChange
       end if
    end if
    return
  end subroutine cole2000Get
