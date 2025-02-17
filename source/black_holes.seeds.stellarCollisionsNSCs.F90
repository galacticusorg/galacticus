!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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
  Implements a black hole seed based on collapse of nuclear star clusters due to runaway stellar collisions.
  !!}

  use :: Mass_Distributions, only : massDistributionClass  , kinematicsDistributionClass

  !![
  <blackHoleSeeds name="blackHoleSeedsstellarCollisionsNSCs">
   <description>
    A model of black hole seeds in which seeds are formed due to the collapse of nuclear star clusters into a black hole. 
     Based on the model of \cite{Vergara_2023} and \cite{Escala_2021}.
   </description>
  </blackHoleSeeds>
  !!]

  type, extends(blackHoleSeedsClass) :: blackHoleSeedsstellarCollisionsNSCs
     !!{
     A node operator class that handle the collapse of nuclear star clusters into a black hole.
     !!}
     private
     double precision                        :: massSingleStar        , radiusSingleStar            , &
         &                                      massEfficiency        , radiusEfficiency            , &
         &                                      massThreshold         
     integer                                 :: stellarMassFormedNSCID, timeStellarMassFormedNSCID     

   contains
     procedure :: mass => stellarCollisionsNSCsSeedMass
     procedure :: spin => stellarCollisionsNSCsSeedSpin
  end type blackHoleSeedsstellarCollisionsNSCs
  
  interface blackHoleSeedsstellarCollisionsNSCs
     !!{
     Constructors for the {\normalfont \ttfamily stellarCollisionsNSCs} black hole seeds class.
     !!}
     module procedure stellarCollisionsNSCsConstructorParameters
     module procedure stellarCollisionsNSCsConstructorInternal
  end interface blackHoleSeedsstellarCollisionsNSCs

  class (massDistributionClass ), pointer :: massDistribution_, massDistributionStellarNSC_ 
  !$omp threadprivate(massDistribution_,massDistributionStellarNSC_)

contains

  function stellarCollisionsNSCsConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily stellarCollisionsNSCs} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (blackHoleSeedsstellarCollisionsNSCs)                :: self
    type (inputParameters                    ), intent(inout) :: parameters
    double precision                                          :: massSingleStar, radiusSingleStar, &
       &                                                         massEfficiency, radiusEfficiency, &
       &                                                         massThreshold 

    !![
    <inputParameter>
      <name>massSingleStar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the mass of a single star in Solar units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusSingleStar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the radius of a single star in Solar units.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>massEfficiency</name>
      <defaultValue>1.0d-1</defaultValue>
      <description>Specifies the efficiency of the mass converted into a black hole seed</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>radiusEfficiency</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the efficiency of the radius used to compute the critical mass</description>
      <source>parameters</source>
    </inputParameter>
     <inputParameter>
      <name>massThreshold</name>
      <defaultValue>1.0d3</defaultValue>
      <description>Specifies the minimum stellar mass to apply the operator</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=blackHoleSeedsstellarCollisionsNSCs(massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold)

    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function stellarCollisionsNSCsConstructorParameters
  
  function stellarCollisionsNSCsConstructorInternal(massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily stellarCollisionsNSCs} node operator class.
    !!}
    implicit none
    type            (blackHoleSeedsstellarCollisionsNSCs)                 :: self
    double precision                                     , intent(in   )  :: massSingleStar
    double precision                                     , intent(in   )  :: radiusSingleStar
    double precision                                     , intent(in   )  :: massEfficiency
    double precision                                     , intent(in   )  :: radiusEfficiency
    double precision                                     , intent(in   )  :: massThreshold

    !![
    <constructorAssign variables="massSingleStar, radiusSingleStar, massEfficiency, radiusEfficiency, massThreshold"/>
    !!]
    !![
    <addMetaProperty   component="NSC" name="agesStellarMassFormed"     id="self%stellarMassFormedNSCID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty   component="NSC" name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedNSCID" isEvolvable="yes" isCreator="no" />
    !!]
    return
  end function stellarCollisionsNSCsConstructorInternal

  double precision function stellarCollisionsNSCsSeedMass(self,node) 
      !!{
        Compute the nuclear star cluster collapse condition.
      !!}
    use :: Galacticus_Nodes                , only : nodeComponentNSC               , nodeComponentNSCStandard, treeNode  , nodeComponentBasic
    use :: Galactic_Structure_Options      , only : componentTypenuclearStarCluster, massTypeStellar
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus, parsec                  , megaParsec, gigaYear
    use :: Numerical_Constants_Prefixes    , only : mega, kilo
    implicit none
    class(blackHoleSeedsstellarCollisionsNSCs), intent(inout)          :: self
    type (treeNode                           ), intent(inout)          :: node
    class(nodeComponentNSC                   )               , pointer :: NSC
    class(nodeComponentBasic                 )               , pointer :: basic
    double precision                                                   :: radiusNSC        , velocityNSC       , &
        &                                                                 massStellarNSC   , massCriticalNSC   , &
        &                                                                 Theta            , crossSectionNSC   , &
        &                                                                 massFormedSeedNSC, massTimeStellarNSC, &
        &                                                                 ageNSC           , time
    double precision                                                   :: velocity        = 100.0d0              !km s¯¹
    double precision                                                   :: sun_rad_to_pc   = 2.2567d-8   
    
    ! Get the nuclear star cluster component.
    NSC         => node%NSC  ()
    radiusNSC   =  self%radiusEfficiency*NSC%radius()*1.0d6
    
    massDistributionStellarNSC_ => node%massDistribution(componentType=componentTypenuclearStarCluster, massType=massTypeStellar)
    velocityNSC = massDistributionStellarNSC_%rotationCurve(radiusNSC*1.0d-6)
    
    select type (NSC)
      type is (nodeComponentNSC)
          return

      class is (nodeComponentNSCStandard)
          if (NSC%massStellar()<= 0.0d0  .or.  NSC%radius() <= 0.0d0) return

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

          call NSC%AgeSet(ageNSC)      !Gyr

          if (ageNSC <= 0.0d0 .or. NSC%Collapse()) then
            stellarCollisionsNSCsSeedMass = 0.0d0
            return
          end if 

          ! Safronov number defined by Binney & Tremaine (2008, https://ui.adsabs.harvard.edu/abs/2008gady.book.....B/abstract)
          Theta            = 9.54d0*(self%massSingleStar/self%radiusSingleStar)*(velocity/velocityNSC)**2.0d0                   !Adimensional
          
          ! Probabilistic mean free path defined as in Landau & Lifshitz (1980, https://ui.adsabs.harvard.edu/abs/1981PhT....34a..74L/abstract) and Shu (1991, https://ui.adsabs.harvard.edu/abs/1991pav..book.....S/abstract)
          CrossSectionNSC  = 16.0d0* sqrt(Pi)*(1+Theta)*(self%radiusSingleStar*sun_rad_to_pc)**2.0d0      !pc²

          ! Critical mass computation using equation (3) in the model of M.C. Vergara, A. Escala, D.R.G. Schleicher and B. Reinoso. (2023, https://ui.adsabs.harvard.edu/abs/2023MNRAS.522.4224V/abstract)
          massCriticalNSC  = (radiusNSC)**(7.0d0/3.0d0)*((4.0d0*Pi*self%massSingleStar)/(3.0d0*CrossSectionNSC*ageNSC*sqrt((gravitationalConstantGalacticus*megaParsec*(kilo*gigaYear)**2.0d0)*parsec**-3.0d0)))**(2.0d0/3.0d0)
          massFormedSeedNSC= self%massEfficiency*NSC%massStellar()
          
          ! Generic type - interrupt and create a standard Black Hole if Nuclear Star Cluster mass is greater than the critical mass.
          if (0.0 <= massCriticalNSC .and. massCriticalNSC <= NSC%massStellar() .and. self%massThreshold <= NSC%massStellar()) then
            call NSC%massSeedSet   ( massFormedSeedNSC)
            call Collapse_Output   (node, radiusNSC, velocityNSC, NSC%massStellar()                  , NSC%massGas(), massCriticalNSC, ageNSC, massFormedSeedNSC)
            call NSC%massStellarset(                              NSC%massStellar()-massFormedSeedNSC                                                           )
            call NSC%CollapseSet   (                                                                                                                      .true.)
            stellarCollisionsNSCsSeedMass = massFormedSeedNSC
            return
          else
            stellarCollisionsNSCsSeedMass = 0.0d0
          end if 
       !![
       <objectDestructor name="massDistributionStellarNSC_"/>
       !!]
    end select
    return
  end function stellarCollisionsNSCsSeedMass

  double precision function stellarCollisionsNSCsSeedSpin(self,node) result(spin)
    !!{
    Compute the spin of the seed black hole.
    !!}
    implicit none
    class(blackHoleSeedsstellarCollisionsNSCs), intent(inout) :: self
    type (treeNode                           ), intent(inout) :: node
    !$GLC attributes unused :: node
    ! Assume zero spin.
    spin=0.0d0
    return
  end function stellarCollisionsNSCsSeedSpin
  
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
    double precision                    , intent(in   ) :: radiusNSC           , velocityNSC, &
        &                                                  massStellarNSC      , massGasNSC , &
        &                                                  massCriticalNSC     , ageNSC     , &
        &                                                  massFormedSeedNSC
    class           (nodeComponentBasic), pointer       :: basic
    type            (hdf5Object        )                :: stellarCollisionsNSCs

    ! Get the basic component.
    basic => node%basic()

    ! Open the group to which black hole formation should be written.
    !$ call hdf5Access%set()
    stellarCollisionsNSCs=outputFile%openGroup("stellarCollisionsNSCs","Data related to the nuclear star cluster collapse.")
    ! Append to the datasets.
    call    stellarCollisionsNSCs%writeDataset([radiusNSC                 ],"radius"      ,"radius of the nuclear star cluster."             ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([velocityNSC               ],"velocity"    ,"velocity mass of the nuclear star cluster."      ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([massStellarNSC            ],"massStellar" ,"Stellar mass of the nuclear star cluster."       ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([massGasNSC                ],"massGas"     ,"Gas mass of the nuclear star cluster."           ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([massCriticalNSC           ],"massCritical","Critical mass of the nuclear star cluster."      ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([ageNSC                    ],"Age"         ,"Age of the nuclear star cluster."                ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([massFormedSeedNSC         ],"Mseed"       ,"Mass of the black hole seed formed."             ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([basic%time()              ],"time"        ,"The time of the nuclear star cluster collapse."  ,appendTo=.true.)
    call    stellarCollisionsNSCs%writeDataset([node%hostTree%volumeWeight],"volumeWeight","The weight of the nuclear star cluster."         ,appendTo=.true.)
    call    stellarCollisionsNSCs%close       (                                                                                                              )
    !$ call hdf5Access %unset       (                                                                                             )
    return
  end subroutine Collapse_Output
