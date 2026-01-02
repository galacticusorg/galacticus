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
  Implements a merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.
  !!}

  use :: Kind_Numbers                           , only : kind_int8
  use :: Satellite_Merging_Progenitor_Properties, only : mergerProgenitorPropertiesClass

  !![
  <mergerRemnantSize name="mergerRemnantSizeCole2000">
   <description>
    A merger remnant size class which uses uses the algorithm of \cite{cole_hierarchical_2000} to compute merger remnant
    spheroid sizes. Specifically
    \begin{equation}
    \frac{(M_1+M_2)^2}{ r_\mathrm{new}} =
    \frac{M_1^2}{r_1} + \frac{M_2^2}{r_2} + \frac{ f_\mathrm{orbit}}{c}
    \frac{M_1 M_2}{r_1+r_2},
    \end{equation}
    where $M_1$ and $M_2$ are the baryonic masses of the components of the merging galaxies that will end up in the spheroid
    \gls{component} of the remnant\footnote{Depending on the merging rules (see \protect\refPhysics{mergerMassMovements}) not
    all mass may be placed into the spheroid \gls{component} of the remnant.} and $r_1$ and $r_2$ are the half mass radii of
    those same components of the merging galaxies\footnote{In practice, \glc\ computes a weighted average of the disk and
    spheroid half-mass radii of each galaxy, with weights equal to the masses of each \gls{component} (disk and spheroid) which
    will become part of the spheroid \gls{component} of the remnant.}, $r_\mathrm{new}$ is the half mass radius of the
    spheroidal \gls{component} of the remnant galaxy and $c$ is a constant which depends on the distribution of the mass. For a
    Hernquist spheroid $c=0.40$ can be found by numerical integration while for a exponential disk $c=0.49$. For simplicity a
    value of $c=0.5$ is adopted for all components. The parameter $f_\mathrm{orbit}=${\normalfont \ttfamily energyOrbital}
    depends on the orbital parameters of the galaxy pair. For example, a value of $f_\mathrm{orbit} = 1$ corresponds to point
    mass galaxies in circular orbits about their center of mass.
    
    A subtlety arises because the above expression accounts for only the baryonic mass of material which becomes part of the
    spheroid \gls{component} of the remnant. In reality, there are additional terms in the energy equation due to the
    interaction of this material with any dark matter mass in each galaxy and any baryonic mass of each galaxy which does not
    become part of the spheroid \gls{component} of the remnant. To account for this additional matter, an effective boost
    factor, $f_\mathrm{boost}$, to the specific angular momentum of each \gls{component} of each merging galaxy is computed:
    \begin{equation}
     f_\mathrm{boost} = {j \over \sqrt{\mathrm{G} M r_{1/2}}},
    \end{equation}
    where $j$ is the specific angular momentum of the component, $M$ is its total baryonic mass and $r_\mathrm{1/2}$ is its
    half-mass radius. The mass-weighted mean boost factor is found by combining those of all components which will form part of
    the spheroid of the remnant. The final specific angular momentum of the remnant spheroid is then given by:
    \begin{equation}
     j_\mathrm{new} = \langle f_\mathrm{boost} \rangle r_\mathrm{new} V_\mathrm{new},
    \end{equation}
    where
    \begin{equation}
     V_\mathrm{new}^2 = {\mathrm{G} (M_1+M_2)\over r_\mathrm{new}}.
    \end{equation}
   </description>
  </mergerRemnantSize>
  !!]
  type, extends(mergerRemnantSizeClass) :: mergerRemnantSizeCole2000
     !!{
     A merger remnant size class which uses the \cite{cole_hierarchical_2000} algorithm.
     !!}
     private
     class           (mergerProgenitorPropertiesClass), pointer :: mergerProgenitorProperties_ => null()
     double precision                                           :: energyOrbital
     integer         (kind=kind_int8                 )          :: lastUniqueID
     logical                                                    :: propertiesCalculated                 , ignoreUnphysicalConditions
     double precision                                           :: radius                               , velocityCircular, &
          &                                                        angularMomentumSpecific
   contains
     final     ::             cole2000Destructor
     procedure :: autoHook => cole2000AutoHook
     procedure :: get      => cole2000Get
  end type mergerRemnantSizeCole2000

  interface mergerRemnantSizeCole2000
     !!{
     Constructors for the \refClass{mergerRemnantSizeCole2000} merger remnant size class.
     !!}
     module procedure cole2000ConstructorParameters
     module procedure cole2000ConstructorInternal
  end interface mergerRemnantSizeCole2000

contains

  function cole2000ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerRemnantSizeCole2000} merger remnant size class which takes a parameter list as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerRemnantSizeCole2000      )                :: self
    type            (inputParameters                ), intent(inout) :: parameters
    class           (mergerProgenitorPropertiesClass), pointer       :: mergerProgenitorProperties_
    double precision                                                 :: energyOrbital
    logical                                                          :: ignoreUnphysicalConditions
    
    !![
    <inputParameter>
      <name>energyOrbital</name>
      <defaultValue>1.0d0</defaultValue>
      <description>The orbital energy used in the ``cole2000'' merger remnant sizes calculation in units of the characteristic orbital energy.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>ignoreUnphysicalConditions</name>
      <defaultValue>.false.</defaultValue>
      <description>If true, ignore unphysical conditions (e.g. negative masses) and leave the size unchanged.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="mergerProgenitorProperties" name="mergerProgenitorProperties_" source="parameters"/>
    !!]
    self=mergerRemnantSizeCole2000(energyOrbital,ignoreUnphysicalConditions,mergerProgenitorProperties_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="mergerProgenitorProperties_"/>
    !!]
    return
  end function cole2000ConstructorParameters

  function cole2000ConstructorInternal(energyOrbital,ignoreUnphysicalConditions,mergerProgenitorProperties_) result(self)
    !!{
    Internal constructor for the \refClass{mergerRemnantSizeCole2000} merger remnant size class.
    !!}
    implicit none
    type            (mergerRemnantSizeCole2000      )                        :: self
    double precision                                 , intent(in   )         :: energyOrbital
    logical                                          , intent(in   )         :: ignoreUnphysicalConditions
    class           (mergerProgenitorPropertiesClass), intent(in   ), target :: mergerProgenitorProperties_
    !![
    <constructorAssign variables="energyOrbital, ignoreUnphysicalConditions, *mergerProgenitorProperties_"/>
    !!]

    self%propertiesCalculated   =.false.
    self%lastUniqueID           =-huge(0_kind_int8)
    self%radius                 =-huge(0.0d0      )
    self%velocityCircular       =-huge(0.0d0      )
    self%angularMomentumSpecific=-huge(0.0d0      )
  end function cole2000ConstructorInternal

  subroutine cole2000AutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels, satelliteMergerEvent
    implicit none
    class(mergerRemnantSizeCole2000), intent(inout) :: self

    call calculationResetEvent%attach(self,cole2000CalculationReset,openMPThreadBindingAllLevels,label='remnantStructure:remnantSizeCole2000')
    call satelliteMergerEvent %attach(self,cole2000GetHook         ,openMPThreadBindingAllLevels,label='remnantStructure:remnantSizeCole2000')
    return
  end subroutine cole2000AutoHook

  subroutine cole2000Destructor(self)
    !!{
    Destructor for the \refClass{mergerRemnantSizeCole2000} merger remnant size class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, satelliteMergerEvent
    implicit none
    type(mergerRemnantSizeCole2000), intent(inout) :: self

    !![
    <objectDestructor name="self%mergerProgenitorProperties_"/>
    !!]
    if (calculationResetEvent%isAttached(self,cole2000CalculationReset)) call calculationResetEvent%detach(self,cole2000CalculationReset)
    if (satelliteMergerEvent %isAttached(self,cole2000GetHook         )) call satelliteMergerEvent %detach(self,cole2000GetHook         )
    return
  end subroutine cole2000Destructor

  subroutine cole2000CalculationReset(self,node,uniqueID)
    !!{
    Reset the dark matter profile calculation.
    !!}
    use :: Error       , only : Error_Report
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (*        ), intent(inout) :: self
    type   (treeNode ), intent(inout) :: node
    integer(kind_int8), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    select type (self)
    class is (mergerRemnantSizeCole2000)
       self%propertiesCalculated=.false.
       self%lastUniqueID       =uniqueID
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine cole2000CalculationReset

  subroutine cole2000GetHook(self,node)
    !!{
    Hookable wrapper around the get function.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (*       ), intent(inout)         :: self
    type            (treeNode), intent(inout), target :: node
    double precision                                  :: radius                 , velocityCircular, &
         &                                               angularMomentumSpecific

    select type (self)
    type is (mergerRemnantSizeCole2000)
       call self%get(node,radius,velocityCircular,angularMomentumSpecific)
    class default
       call Error_Report('incorrect class'//{introspection:location})
    end select
    return
  end subroutine cole2000GetHook

  subroutine cole2000Get(self,node,radius,velocityCircular,angularMomentumSpecific)
    !!{
    Compute the size of the merger remnant for {\normalfont \ttfamily node} using the \cite{cole_hierarchical_2000} algorithm.
    !!}
    use :: Display                         , only : displayMessage                , verbosityLevelSilent
    use :: Galactic_Structure_Options      , only : massTypeDark
    use :: Error                           , only : Error_Report
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Comparison            , only : Values_Agree
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: String_Handling                 , only : operator(//)
    implicit none
    class           (mergerRemnantSizeCole2000), intent(inout) :: self
    type            (treeNode                 ), intent(inout) :: node
    double precision                           , intent(  out) :: radius                           , velocityCircular         , &
         &                                                        angularMomentumSpecific
    type            (treeNode                 ), pointer       :: nodeHost
    class           (massDistributionClass    ), pointer       :: massDistributionSatellite        , massDistributionHost
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

    ! The calculation of remnant size is computed when first needed and then stored. This ensures that the results are determined
    ! by the properties of the merge target prior to any modification that will occur as node components are modified in response
    ! to the merger.
    if (node%uniqueID() /= self%lastUniqueID) call cole2000CalculationReset(self,node,node%uniqueID())
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
       if (massSpheroidSatellite <= 0.0d0 .and. Values_Agree(massSpheroidHost,massSpheroidHostPreMerger,relTol=toleranceMassRelative)) then
          self%radius                 =remnantNoChange
          self%velocityCircular       =remnantNoChange
          self%angularMomentumSpecific=remnantNoChange
       else
          ! Check that the properties of the galaxies are physically reasonable.
          if (.not.self%ignoreUnphysicalConditions) then
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
                message=message//' (radius:mass:massSpheroid='//trim(dataString)//')'
                call displayMessage(message,verbosityLevelSilent)
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
                message=message//' (radius:mass:massSpheroid='//trim(dataString)//')'
                call displayMessage(message,verbosityLevelSilent)
                errorCondition=.true.
             end if
             if (errorCondition) then
                call node    %serializeASCII(verbosityLevelSilent)
                call nodeHost%serializeASCII(verbosityLevelSilent)
                call Error_Report('error condition detected'//{introspection:location})
             end if
          end if
          ! Check if host has finite mass.
          if (massSpheroidSatellite+massSpheroidHost > 0.0d0) then
             ! Compute masses of dark matter within the host and satellite radii.
             massDistributionHost      => nodeHost%massDistribution(massType=massTypeDark)
             massDistributionSatellite => node    %massDistribution(massType=massTypeDark)
             massDarkMatterHost        =  massDistributionHost     %massEnclosedBySphere(radiusHost     )
             massDarkMatterSatellite   =  massDistributionSatellite%massEnclosedBySphere(radiusSatellite)
             !![
	     <objectDestructor name="massDistributionHost"     />
	     <objectDestructor name="massDistributionSatellite"/>
	     !!]
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
             self%radius=(massSpheroidTotalSatellite+massSpheroidHostTotal)**2/energyProgenitors
             ! Also compute the specific angular momentum at the half-mass radius.
             self%velocityCircular       =sqrt(gravitationalConstant_internal*(massSpheroidSatellite+massSpheroidHost)/self%radius)
             self%angularMomentumSpecific=self%radius*self%velocityCircular*factorAngularMomentum
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
  end subroutine cole2000Get
