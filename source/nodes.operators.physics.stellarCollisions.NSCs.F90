!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
!!    Andrew Benson <abenson@carnegiescience.edu>
!!
!! This file is part of Galacticus.
!!
!!    Galacticus is free software: you can redistribute it and/or modify
!!    it under the terms of the GNU General Public License as published by
!!    the Free Software Foundation, e ither version 3 of the License, or
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
  Implements a node operator class that handle the collapse of nuclear star clusters into a black hole.
  !!}
  use :: Mass_Distributions, only : massDistributionClass  , kinematicsDistributionClass

  !![
  <nodeOperator name="nodeOperatorstellarCollisionsNSC">
   <description>
     A node operator class that handle the collapse of nuclear star clusters into a black hole. 
     Based on the model of \cite{vergara_global_2023} and \cite{escala_observational_2021}.
   </description>
  </nodeOperator>
  !!]

  type, extends(nodeOperatorClass) :: nodeOperatorstellarCollisionsNSC
     !!{
     A node operator class that handle the formation of a BH seed within nuclear star clusters due to runaway stellar collisions.
     !!}
     private
     double precision                        :: massSingleStar        , radiusSingleStar            , &
         &                                      massEfficiency        , radiusEfficiency            , &
         &                                      massThreshold 
     integer                                 :: stellarMassFormedNSCID, timeStellarMassFormedNSCID     
   contains
     procedure :: differentialEvolution => stellarCollisionsNSCDifferentialEvolution
  end type nodeOperatorstellarCollisionsNSC
  
  interface nodeOperatorstellarCollisionsNSC
     !!{
     Constructors for the {\normalfont \ttfamily stellarCollisionsNSC} node operator class.
     !!}
     module procedure stellarCollisionsNSCConstructorParameters
     module procedure stellarCollisionsNSCConstructorInternal
  end interface nodeOperatorstellarCollisionsNSC

  class (massDistributionClass ), pointer :: massDistribution_, massDistributionStellarNSC_ 
  !$omp threadprivate(massDistribution_,massDistributionStellarNSC_)

contains

  function stellarCollisionsNSCConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily stellarCollisionsNSC} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorstellarCollisionsNSC)                :: self
    type (inputParameters                 ), intent(inout) :: parameters
    double precision                                       :: massSingleStar    , radiusSingleStar, &
       &                                                      massEfficiency    , radiusEfficiency, &
       &                                                      massThreshold

    !![
    <inputParameter>
      <name>massSingleStar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the mass of a single star in solar units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusSingleStar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the radius of a single star in solar units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massEfficiency</name>
      <defaultValue>1.0d-1</defaultValue>
      <description>Specifies the efficiency of the mass converted into a BH seed.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusEfficiency</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the efficiency of the radius used to compute the critical mass.</description>
      <source>parameters</source>
    </inputParameter>
     <inputParameter>
      <name>massThreshold</name>
      <defaultValue>1.0d3</defaultValue>
      <description>Specifies the minimum stellar mass to apply the operator.</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=nodeOperatorstellarCollisionsNSC(massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold)

    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function stellarCollisionsNSCConstructorParameters
  
  function stellarCollisionsNSCConstructorInternal(massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily stellarCollisionsNSC} node operator class.
    !!}
    implicit none
    type            (nodeOperatorstellarCollisionsNSC)                :: self
    double precision                                  , intent(in   ) :: massSingleStar
    double precision                                  , intent(in   ) :: radiusSingleStar
    double precision                                  , intent(in   ) :: massEfficiency
    double precision                                  , intent(in   ) :: radiusEfficiency
    double precision                                  , intent(in   ) :: massThreshold

    !![
    <constructorAssign variables="massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold"/>
    !!]

    !![
    <addMetaProperty   component="NSC" name="agesStellarMassFormed"     id="self%stellarMassFormedNSCID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty   component="NSC" name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedNSCID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function stellarCollisionsNSCConstructorInternal

  subroutine stellarCollisionsNSCDifferentialEvolution(self,node,interrupt,functioninterrupt,propertyType)
      !!{
        Compute the nuclear star cluster collapse condition.
      !!}
    use :: Galacticus_Nodes                , only : interruptTask                  , nodeComponentNSC, nodeComponentNSCStandard, propertyInactive, &
        &                                           nodeComponentBasic             , treeNode        
    use :: Galactic_Structure_Options      , only : componentTypeNuclearStarCluster, massTypeStellar
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal , parsec          , megaParsec              , gigaYear
    use :: Numerical_Constants_Prefixes    , only : mega, kilo
    implicit none
    class(nodeOperatorstellarCollisionsNSC), intent(inout), target :: self
    type (treeNode                        ), intent(inout), target :: node
    logical                                , intent(inout)         :: interrupt
    procedure (interruptTask              ), intent(inout), pointer:: functioninterrupt
    integer                                , intent(in   )         :: propertyType
    class(nodeComponentNSC                )               , pointer:: NSC
    class(nodeComponentBasic              )               , pointer:: basic
    double precision                                               :: radiusNSC                  , velocityNSC       , &
        &                                                             massStellarNSC             , massCriticalNSC   , &
        &                                                             Theta                      , crossSectionNSC   , &
        &                                                             massFormedSeedNSC          , massTimeStellarNSC, &
        &                                                             ageNSC                     , time
    double precision                                               :: velocity        = 100.0d0                       !km s¯¹
    double precision                                               :: solarRadiusTopc = 2.2567d-8                     !pc

    ! Return immediately if inactive variables are requested.
    if (propertyInactive(propertyType)) return
    
    ! Get the nuclear star cluster component.
    NSC         => node%NSC  ()
    radiusNSC   =  self%radiusEfficiency*1.0d6*NSC%radius() !pc
    
    massDistributionStellarNSC_ => node%massDistribution(componentType=componentTypeNuclearStarCluster, massType=massTypeStellar)
    !Evaluates the radius in parsec.
    velocityNSC = massDistributionStellarNSC_%rotationCurve(radiusNSC*1.0e-6) 
    
    ! Detect nuclear star cluster component type.
    select type (NSC)
      type is (nodeComponentNSC)
          ! Nothing to do here.
          return

      class is (nodeComponentNSCStandard)
          ! Return for unphysical NSCs.
          if (NSC%massStellar()<=0.0d0.or.NSC%radius() <= 0.0d0) return

          massStellarNSC        = NSC%floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID)
          massTimeStellarNSC    = NSC%floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID)
          
          basic => node%basic()
          time  =  basic%time()

          !Trap cases after merger, when the age is reset
          if ( 0.0d0 < massStellarNSC) then 
              ageNSC= +time               &
                   &  -massTimeStellarNSC &
                   &  /massStellarNSC
          else 
              ageNSC= 0.0d0
          end if 

          call NSC%ageSet(ageNSC)      !Gyr

          if (ageNSC <= 0.0d0 .or. NSC%Collapse()) return 

          ! Safronov number defined by Binney & Tremaine (2008, https://ui.adsabs.harvard.edu/abs/2008gady.book.....B/abstract)
          Theta            = 9.54d0*(self%massSingleStar/self%radiusSingleStar)*(velocity/velocityNSC)**2.0d0
          
          ! Probabilistic mean free path defined as in Landau & Lifshitz (1980, https://ui.adsabs.harvard.edu/abs/1981PhT....34a..74L/abstract) and Shu (1991, https://ui.adsabs.harvard.edu/abs/1991pav..book.....S/abstract)
          CrossSectionNSC  = 16.0d0* sqrt(Pi)*(1+Theta)*(self%radiusSingleStar*solarRadiusTopc)**2.0d0

          ! Critical mass computation using equation (3) in the model of M.C. Vergara, A. Escala, D.R.G. Schleicher and B. Reinoso. (2023, https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.4224V/abstract)
          massCriticalNSC  = radiusNSC**(7.0d0/3.0d0)*((4.0d0*Pi*self%massSingleStar)/(3.0d0*CrossSectionNSC*ageNSC*sqrt((gravitationalConstant_internal*megaParsec*(kilo*gigaYear)**2.0d0)*parsec**-3.0d0)))**(2.0d0/3.0d0)
          massFormedSeedNSC= self%massEfficiency*NSC%massStellar()
          
          call NSC%massCriticalSet(massCriticalNSC)
          ! Generic type - interrupt and create a standard Black Hole if Nuclear Star Cluster mass is greater than the critical mass.
          if (0.0 <= massCriticalNSC .and. massCriticalNSC <= NSC%massStellar() .and. self%massThreshold <= NSC%massStellar()) then
            call NSC%massSeedSet   ( massFormedSeedNSC)
            interrupt=.true.
            functionInterrupt => BlackHoleStandardCreate
            call Collapse_Output   (node, radiusNSC, velocityNSC, NSC%massStellar(), NSC%massGas(), massCriticalNSC, ageNSC, massFormedSeedNSC)
            call NSC%massStellarset(                              NSC%massStellar()                                         -massFormedSeedNSC)
            call NSC%CollapseSet   (                                                                                                    .true.)
            return
          end if 
       !![
       <objectDestructor name="massDistributionStellarNSC_"/>
       !!]
    end select
    return
  end subroutine stellarCollisionsNSCDifferentialEvolution

  subroutine BlackHoleStandardCreate(node,timeEnd)
    !!{
    Creates a black hole component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode, nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    type (treeNode              ), intent(inout), target   :: node
    double precision             , intent(in   ), optional :: timeEnd
    class(nodeComponentBlackHole)               , pointer  :: blackHole
    class(nodeComponentNSC      )               , pointer  :: NSC
    !$GLC attributes unused :: timeEnd

    NSC       => node%NSC      (                 )
    ! Create the black hole.
    blackHole => node%blackHole(autoCreate=.true.)
    ! Set to the seed mass.
    call blackHole%          massSet(NSC%massSeed())
    call blackHole%          spinSet(         0.0d0)  
    call blackHole%radialPositionSet(         0.0d0)
    call blackHole%    NSCChannelSet(        .true.)
    return
  end subroutine BlackHoleStandardCreate
  
  subroutine Collapse_Output(node, radiusNSC, velocityNSC, massStellarNSC, massGasNSC, massCriticalNSC, ageNSC, massFormedSeedNSC)
    !!{
    Outputs properties of collapsing NSCs.
    !!}
    use :: Output_HDF5     , only : outputFile
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : hdf5Object
    implicit none
    type            (treeNode          ), intent(inout) :: node
    double precision                    , intent(in   ) :: radiusNSC        , velocityNSC, &
        &                                                  massStellarNSC   , massGasNSC , &
        &                                                  massCriticalNSC  , ageNSC     , &
        &                                                  massFormedSeedNSC
    class           (nodeComponentBasic), pointer       :: basic
    type            (hdf5Object        )                :: stellarCollisionsNSC

    ! Get the basic component.
    basic => node%basic()

    ! Open the group to which black hole formation should be written.
    !$ call hdf5Access%set()
    stellarCollisionsNSC=outputFile%openGroup("stellarCollisionsNSC","Data related to the black hole formation due to  stellar collisions in NSCs.")
    ! Append to the datasets.
    call    stellarCollisionsNSC%writeDataset([radiusNSC                 ],"radius"      ,"radius of the NSC."             ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([velocityNSC               ],"velocity"    ,"velocity mass of the NSC."      ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([massStellarNSC            ],"massStellar" ,"Stellar mass of the NSC."       ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([massGasNSC                ],"massGas"     ,"Gas mass of the NSC."           ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([massCriticalNSC           ],"massCritical","Critical mass of the NSC."      ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([ageNSC                    ],"Age"         ,"Age of the NSC."                ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([massFormedSeedNSC         ],"seedMass"    ,"Mass of the BH seed formed."    ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([basic%time()              ],"time"        ,"The time of the NSC collapse."  ,appendTo=.true.)
    call    stellarCollisionsNSC%writeDataset([node%hostTree%volumeWeight],"volumeWeight","The weight of the NSC."         ,appendTo=.true.)
    call    stellarCollisionsNSC%close       (                                                                                             )
    !$ call hdf5Access %unset       (                                                                                             )
    return
  end subroutine Collapse_Output
