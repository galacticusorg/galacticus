!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Implements a node operator class that handle the collapse of nuclear star clusters into a black hole.
  !!}
  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  !![
  <nodeOperator name="nodeOperatorNSCCollapse">
   <description>A node operator class that handle the collapse of nuclear star clusters into a black hole.</description>
  </nodeOperator>
  !!]

  type, extends(nodeOperatorClass) :: nodeOperatorNSCCollapse
     !!{
     A node operator class that handle the collapse of nuclear star clusters into a black hole.
     !!}
     private
     class  (cosmologyFunctionsClass), pointer :: cosmologyFunctions_       => null()
     double precision                          :: Mstar                 , Rstar, Mefficiency, Refficiency, MThreshold
     integer                                   :: stellarMassFormedNSCID, timeStellarMassFormedNSCID     

   contains
     final     ::                          NSCCollapseDestructor 
     procedure :: differentialEvolution => NSCCollapseDifferentialEvolution
  end type nodeOperatorNSCCollapse
  
  interface nodeOperatorNSCCollapse
     !!{
     Constructors for the {\normalfont \ttfamily NSCCollapse} node operator class.
     !!}
     module procedure NSCCollapseConstructorParameters
     module procedure NSCCollapseConstructorInternal
  end interface nodeOperatorNSCCollapse


contains

  function NSCCollapseConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily NSCCollapse} node operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameters
    implicit none
    type (nodeOperatorNSCCollapse  )                :: self
    type (inputParameters          ), intent(inout) :: parameters
    class(cosmologyFunctionsClass  ), pointer       :: cosmologyFunctions_
    double precision                                :: Mstar, Rstar, Mefficiency, Refficiency, MThreshold

    !![
    <inputParameter>
      <name>Mstar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the mass of a single star in solar mass.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>Rstar</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the radius of a single star in solar radii.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>Mefficiency</name>
      <defaultValue>1.0d-1</defaultValue>
      <description>Specifies the efficiency of the mass converted into a BH</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>Refficiency</name>
      <defaultValue>1.0d0</defaultValue>
      <description>Specifies the efficiency of the radius used to compute the critcal mass</description>
      <source>parameters</source>
    </inputParameter>
     <inputParameter>
      <name>MThreshold</name>
      <defaultValue>1.0d3</defaultValue>
      <description>Specifies the minimum mass to apply the operator</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    self=nodeOperatorNSCCollapse(Mstar, Rstar, Mefficiency, Refficiency, MThreshold, cosmologyFunctions_)

    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function NSCCollapseConstructorParameters
  
  function NSCCollapseConstructorInternal(Mstar, Rstar, Mefficiency, Refficiency, MThreshold, cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily NSCCollapse} node operator class.
    !!}
    implicit none
    type            (nodeOperatorNSCCollapse)                        :: self
    class           (cosmologyFunctionsClass), intent(in   ), target :: cosmologyFunctions_
    double precision                         , intent(in   )         :: Mstar
    double precision                         , intent(in   )         :: Rstar
    double precision                         , intent(in   )         :: Mefficiency
    double precision                         , intent(in   )         :: Refficiency
    double precision                         , intent(in   )         :: MThreshold

    !![
    <addMetaProperty component="NSC"      name="agesStellarMassFormed"     id="self%stellarMassFormedNSCID"     isEvolvable="yes" isCreator="no" />
    <addMetaProperty component="NSC"      name="agesTimeStellarMassFormed" id="self%timeStellarMassFormedNSCID" isEvolvable="yes" isCreator="no" />
    <constructorAssign variables="Mstar, Rstar, Mefficiency, Refficiency, MThreshold, *cosmologyFunctions_"/>
    !!]
    return
  end function NSCCollapseConstructorInternal

  subroutine NSCCollapseDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily nodeOperatorNSCCollapse} node operator class
    !!} 
    implicit none
    type(nodeOperatorNSCCollapse), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"/>
    !!]
    return
  end subroutine NSCCollapseDestructor

  subroutine NSCCollapseDifferentialEvolution(self,node,interrupt,functioninterrupt,propertyType)
      !!{
        Compute the nuclear star cluster collapse condition.
      !!}
    use :: Galacticus_Nodes                , only : interruptTask    , nodeComponentNSC, nodeComponentNSCStandard, &
                                                &   propertyInactive , treeNode        , nodeComponentBasic
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstantGalacticus

    implicit none
    class(nodeOperatorNSCCollapse ), intent(inout), target :: self
    type (treeNode                ), intent(inout), target :: node
    logical                        , intent(inout)         :: interrupt
    procedure (interruptTask      ), intent(inout), pointer:: functioninterrupt
    integer                        , intent(in   )         :: propertyType
    class(nodeComponentNSC        )               , pointer:: NSC
    class(nodeComponentBasic      )               , pointer:: basic
    double precision                                       :: NSCradius, NSCvelocity, NSCStellar
    double precision                                       :: Mcrit    , Theta      , CrossSection, Mseed
    double precision                                       :: vel           = 100.0    !!km s^-1
    double precision                                       :: sun_rad_to_pc = 2.2567d-8  ! 
    double precision                                       :: massStellarNSC, massTimeStellarNSC, ageNSC, time
    ! Return immediately if inactive variables are requested.

    if (propertyInactive(propertyType)) return
    
    ! Get the nuclear star cluster component.
    NSC => node%NSC()

    NSCradius   = self%Refficiency*1.0d6*NSC%radius     () !pc
    NSCvelocity =                        NSC%velocity   ()
    NSCStellar  =                        NSC%massStellar() 

    ! Detect nuclear star cluster component type.
    select type (NSC)
      type is (nodeComponentNSC)
      return
      class is (nodeComponentNSCStandard)

          if (NSC%massStellar()<= 0.0d0  .or.  NSC%radius() <= 0.0d0 .or. NSC%Collapse() ) return
          massStellarNSC        =NSC    %floatRank0MetaPropertyGet(self%    stellarMassFormedNSCID)
          massTimeStellarNSC    =NSC    %floatRank0MetaPropertyGet(self%timeStellarMassFormedNSCID)
          basic => node %basic()
          time  =  basic%time ()
          ageNSC=+time                    &
               &      -massTimeStellarNSC &
               &      /massStellarNSC

          call NSC%AgeSet(ageNSC)                                       !pc
          if (ageNSC <= 0.0d0) return 
          Theta         = 9.54 *     (self%Mstar   / self%Rstar)*(vel/NSCvelocity)**2.0 !Adimensional
          CrossSection  = 16.0 * sqrt(Pi)*(1+Theta)*(self%Rstar *sun_rad_to_pc)   **2.0 !pc^2
          Mcrit         = NSCradius**(7.0/3.0)*((4.0*Pi*self%Mstar)/(3.0*CrossSection*ageNSC*sqrt(1.0d6*gravitationalConstantGalacticus*(3.2408d-14)**2/(3.1710d-17)**2)))**(2.0/3.0)
          Mseed         = self%Mefficiency*NSCStellar

          call NSC%CriticalMassSet(Mcrit)
          ! Generic type - interrupt and create a standard Black Hole if Nuclear Star Cluster mass is greater than the critical mass.
          if (0.0 <= Mcrit .and. Mcrit <= NSCStellar .and. self%MThreshold < NSCStellar) then
            call NSC%massSeedSet   ( Mseed)
            interrupt=.true.
            functionInterrupt => BlackHoleStandardCreate
            call NSC%massStellarset(NSCStellar-Mseed                                           )
            call NSC%CollapseSet   (.true.                                                     )
            call Collapse_Output   (node, NSCradius, NSCvelocity, NSCStellar, Mcrit, ageNSC, Mseed)
            return
          end if 
    end select
    return
  end subroutine NSCCollapseDifferentialEvolution

  subroutine BlackHoleStandardCreate(node)
    !!{
    Creates a black hole component for {\normalfont \ttfamily node}.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentBlackHole, treeNode, nodeComponentNSC, nodeComponentNSCStandard
    implicit none
    type (treeNode              ), intent(inout), target  :: node
    class(nodeComponentBlackHole)               , pointer :: blackHole
    class(nodeComponentNSC      )               , pointer :: NSC

    NSC       => node%NSC      (                 )
    ! Create the black hole.
    blackHole => node%blackHole(autoCreate=.true.)
    ! Set to the seed mass.
    call blackHole%          massSet(NSC      %massSeed())
    call blackHole%          spinSet(blackHole%spinSeed())  
    call blackHole%radialPositionSet(               0.0d0)
    call blackHole%    NSCChannelSet(              .true.)
    return
  end subroutine BlackHoleStandardCreate
  
  subroutine Collapse_Output(node, radius, velocity, massStellar, massCritical, ageNSC, Mseed)
    !!{
    Outputs properties of collapsing NSCs.
    !!}
    use :: Output_HDF5     , only : outputFile
    use :: Galacticus_Nodes, only : nodeComponentBasic, treeNode
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : hdf5Object
    implicit none
    type            (treeNode          ), intent(inout) :: node
    double precision                    , intent(in   ) :: radius     , velocity, massStellar, massCritical, ageNSC, MSeed
    class           (nodeComponentBasic), pointer       :: basic
    type            (hdf5Object        )                :: NSCCollapse

    ! Get the basic component.
    basic => node%basic()

    ! Open the group to which black hole mergers should be written.
    !$ call hdf5Access%set()
    NSCCollapse=outputFile%openGroup("NSCCollapse","Data related to the NSC collapse.")
    ! Append to the datasets.
    call    NSCCollapse%writeDataset([radius                    ],"radius"      ,"radius mass of the NSC."        ,appendTo=.true.)
    call    NSCCollapse%writeDataset([velocity                  ],"velocity"    ,"velocity mass of the NSC."      ,appendTo=.true.)
    call    NSCCollapse%writeDataset([massStellar               ],"massStellar" ,"Stellar mass of the NSC."       ,appendTo=.true.)
    call    NSCCollapse%writeDataset([massCritical              ],"massCritical","Critical mass of the NSC."      ,appendTo=.true.)
    call    NSCCollapse%writeDataset([ageNSC                    ],"Age"         ,"Age of the NSC."                ,appendTo=.true.)
    call    NSCCollapse%writeDataset([Mseed                     ],"Mseed"       ,"seed mass of the BH."           ,appendTo=.true.)
    call    NSCCollapse%writeDataset([basic%time()              ],"time"        ,"The time of the NSC collapse."  ,appendTo=.true.)
    call    NSCCollapse%writeDataset([node%hostTree%volumeWeight],"volumeWeight","The weight of the NSC."         ,appendTo=.true.)
    call    NSCCollapse%close       (                                                                                             )
    !$ call hdf5Access  %unset      (                                                                                             )
    return
  end subroutine Collapse_Output
