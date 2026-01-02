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
  An implementation of the merger tree importer class for \glc\ format merger tree files.
  !!}

  use :: Cosmological_Density_Field, only : cosmologicalMassVarianceClass
  use :: Cosmology_Functions       , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters      , only : cosmologyParametersClass
  use :: Halo_Mass_Functions       , only : haloMassFunctionClass
  use :: IO_HDF5                   , only : hdf5Object
  use :: Stateful_Types            , only : statefulDouble               , statefulInteger, statefulLogical

  ! Enumeration of particle epoch types.
  !![
  <enumeration>
   <name>galacticusParticleEpochType</name>
   <description>Particle epoch type enumerations.</description>
   <entry label="time"           />
   <entry label="expansionFactor"/>
   <entry label="redshift"       />
  </enumeration>
  !!]

  type, public, extends(nodeData) :: nodeDataGalacticus
     !!{
     Extension of the {\normalfont \ttfamily nodeData} class for \glc\ format merger trees. Stores particle indices and counts for nodes.
     !!}
     integer(c_size_t) :: particleIndexCount, particleIndexStart
  end type nodeDataGalacticus

  !![
  <mergerTreeImporter name="mergerTreeImporterGalacticus">
    <description>
    A merger tree importer class which imports trees from an HDF5 file. HDF5 file should follow the general purpose format
    described \href{https://github.com/galacticusorg/galacticus/wiki/Merger-Tree-File-Format}{here}. An example of how to
    construct such a file can be found in the {\normalfont \ttfamily tests/nBodyMergerTrees} folder. In that folder, the
    {\normalfont \ttfamily getMillenniumTrees.pl} script will retrieve a sample of merger trees from the
    \href{https://virgodb.dur.ac.uk:8443/Millennium/}{Millennium Simulation database} and use the {\normalfont \ttfamily
    Merger\_Tree\_File\_Maker.exe} code supplied with \glc\ to convert these into an HDF5 file suitable for reading into \glc. The
    {\normalfont \ttfamily getMillenniumTrees.pl} script requires you to have a username and password to access the Millennium
    Simulation database\footnote{If you do not have a username and password for the Millennium Simulation database you can request
    one from \href{mailto:contact@g-vo.org}{\normalfont \ttfamily contact@g-vo.org}.}. These can be entered manually or stored in
    a section of the
    \href{https://github.com/galacticusorg/galacticus/releases/download/bleeding-edge/Galacticus_Usage.pdf#sec.ConfigFile}{\normalfont
    \ttfamily galacticusConfig.xml} file as follows:   
    \begin{verbatim}
      &lt;millenniumDB>
        &lt;host>
          &lt;name>^myHost$&lt;/name>
          &lt;user>myUserName&lt;/user>
          &lt;passwordFrom>input&lt;/passwordFrom>
          &lt;treePath>/path/to/trees&lt;/treePath>
        &lt;/host>
        &lt;host>
          &lt;name>default&lt;/name>
          &lt;user>myUserName&lt;/user>
          &lt;password>myPassword&lt;/password>
        &lt;/host>
      &lt;/millenniumDB>
    \end{verbatim}
    Here, each {\normalfont \ttfamily host} section describes rules for a given computer (with ``default'' being used if no specific
    match to the regular expression give in {\normalfont \ttfamily name} is found). The {\normalfont \ttfamily user} element gives the
    user name to use, while the {\normalfont \ttfamily passwordFrom} element specifies how the password should be obtained. Currently
    the only allowed mechanism is ``input'', in which case the password is read from standard input. Alternatively, you can include a
    {\normalfont \ttfamily password} element which contains the password itself. Of course, this is insecure\ldots
    
    The optional {\normalfont \ttfamily treePath} element gives the location where merger trees from the Millennium Simulation can be
    stored. Some scripts will make use of this location so that Millennium Simulation merger trees can be shared between multiple
    scripts.
    </description>
  </mergerTreeImporter>
  !!]
  type, extends(mergerTreeImporterClass) :: mergerTreeImporterGalacticus
     !!{
     A merger tree importer class for \glc\ format merger tree files.
     !!}
     private
     class           (cosmologyFunctionsClass                   ), pointer                   :: cosmologyFunctions_       => null()
     class           (haloMassFunctionClass                     ), pointer                   :: haloMassFunction_         => null()
     class           (cosmologyParametersClass                  ), pointer                   :: cosmologyParameters_      => null()
     class           (cosmologicalMassVarianceClass             ), pointer                   :: cosmologicalMassVariance_ => null()
     type            (hdf5Object                                )                            :: file                               , forestHalos
     type            (statefulInteger                           )                            :: hasSubhalos                        , areSelfContained              , &
          &                                                                                     includesHubbleFlow                 , periodicPositions             , &
          &                                                                                     lengthStatus
     type            (statefulLogical                           )                            :: massesAreInclusive                 , angularMomentaAreInclusive
     type            (statefulDouble                            )                            :: length
     type            (importerUnits                             )                            :: massUnit                           , lengthUnit                    , &
          &                                                                                     timeUnit                           , velocityUnit
     logical                                                                                 :: fatalMismatches                    , forestIndicesRead             , &
          &                                                                                     angularMomentaIsScalar             , angularMomentaIsVector        , &
          &                                                                                     spinIsScalar                       , spinIsVector                  , &
          &                                                                                     reweightTrees                      , validateData
     integer                                                                                 :: forestsCount                       , formatVersion
     integer                                                     , allocatable, dimension(:) :: firstNodes                         , nodeCounts
     integer         (kind=kind_int8                            ), allocatable, dimension(:) :: forestIndices
     double precision                                            , allocatable, dimension(:) :: weights
     type            (hdf5Object                                )                            :: particles
     type            (enumerationGalacticusParticleEpochTypeType)                            :: particleEpochType
     type            (varying_string                            )                            :: particleEpochDataSetName
     character       (len=32                                    )                            :: forestHalosGroupName               , forestContainmentAttributeName, &
          &                                                                                     forestIndexGroupName               , forestIndexDatasetName        , &
          &                                                                                     forestWeightDatasetName
   contains
     final     ::                                  galacticusDestructor
     procedure :: open                          => galacticusOpen
     procedure :: close                         => galacticusClose
     procedure :: canReadSubsets                => galacticusCanReadSubsets
     procedure :: treesHaveSubhalos             => galacticusTreesHaveSubhalos
     procedure :: massesIncludeSubhalos         => galacticusMassesIncludeSubhalos
     procedure :: angularMomentaIncludeSubhalos => galacticusAngularMomentaIncludeSubhalos
     procedure :: treesAreSelfContained         => galacticusTreesAreSelfContained
     procedure :: velocitiesIncludeHubbleFlow   => galacticusVelocitiesIncludeHubbleFlow
     procedure :: positionsArePeriodic          => galacticusPositionsArePeriodic
     procedure :: cubeLength                    => galacticusCubeLength
     procedure :: treeWeight                    => galacticusTreeWeight
     procedure :: treeCount                     => galacticusTreeCount
     procedure :: treeIndex                     => galacticusTreeIndex
     procedure :: nodeCount                     => galacticusNodeCount
     procedure :: positionsAvailable            => galacticusPositionsAvailable
     procedure :: scaleRadiiAvailable           => galacticusScaleRadiiAvailable
     procedure :: particleCountAvailable        => galacticusParticleCountAvailable
     procedure :: velocityMaximumAvailable      => galacticusVelocityMaximumAvailable
     procedure :: velocityDispersionAvailable   => galacticusVelocityDispersionAvailable
     procedure :: angularMomentaAvailable       => galacticusAngularMomentaAvailable
     procedure :: angularMomenta3DAvailable     => galacticusAngularMomenta3DAvailable
     procedure :: spinAvailable                 => galacticusSpinAvailable
     procedure :: spin3DAvailable               => galacticusSpin3DAvailable
     procedure :: import                        => galacticusImport
     procedure :: subhaloTrace                  => galacticusSubhaloTrace
     procedure :: subhaloTraceCount             => galacticusSubhaloTraceCount
  end type mergerTreeImporterGalacticus

  interface mergerTreeImporterGalacticus
     !!{
     Constructors for the \glc\ format merger tree importer class.
     !!}
     module procedure galacticusConstructorParameters
     module procedure galacticusConstructorInternal
  end interface mergerTreeImporterGalacticus

contains

  function galacticusConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \glc\ format merger tree importer which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (mergerTreeImporterGalacticus )                :: self
    type   (inputParameters              ), intent(inout) :: parameters
    class  (cosmologyFunctionsClass      ), pointer       :: cosmologyFunctions_
    class  (haloMassFunctionClass        ), pointer       :: haloMassFunction_
    class  (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class  (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
    logical                                               :: fatalMismatches          , reweightTrees, &
         &                                                   validateData

    !![
    <inputParameter>
      <name>fatalMismatches</name>
      <defaultValue>.true.</defaultValue>
      <description>Specifies whether mismatches in cosmological parameter values between \glc\ and the merger tree file should be considered fatal.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>reweightTrees</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether merger tree weights should be recomputed from the halo mass function.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>validateData</name>
      <defaultValue>.false.</defaultValue>
      <description>If true perform some validation of imported data to identify possible problems.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"       name="cosmologyFunctions_"       source="parameters"/>
    <objectBuilder class="haloMassFunction"         name="haloMassFunction_"         source="parameters"/>
    <objectBuilder class="cosmologyParameters"      name="cosmologyParameters_"      source="parameters"/>
    <objectBuilder class="cosmologicalMassVariance" name="cosmologicalMassVariance_" source="parameters"/>
    !!]
    self=mergerTreeImporterGalacticus(fatalMismatches,reweightTrees,validateData,cosmologyFunctions_,haloMassFunction_,cosmologyParameters_,cosmologicalMassVariance_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"      />
    <objectDestructor name="haloMassFunction_"        />
    <objectDestructor name="cosmologyParameters_"     />
    <objectDestructor name="cosmologicalMassVariance_"/>
    !!]
    return
  end function galacticusConstructorParameters

  function galacticusConstructorInternal(fatalMismatches,reweightTrees,validateData,cosmologyFunctions_,haloMassFunction_,cosmologyParameters_,cosmologicalMassVariance_) result(self)
    !!{
    Internal constructor for the \glc\ format merger tree importer.
    !!}
    implicit none
    type   (mergerTreeImporterGalacticus )                        :: self
    logical                               , intent(in   )         :: fatalMismatches          , reweightTrees, &
         &                                                           validateData
    class  (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_
    class  (haloMassFunctionClass        ), intent(in   ), target :: haloMassFunction_
    class  (cosmologyParametersClass     ), intent(in   ), target :: cosmologyParameters_
    class  (cosmologicalMassVarianceClass), intent(in   ), target :: cosmologicalMassVariance_
    !![
    <constructorAssign variables="fatalMismatches, reweightTrees, validateData, *cosmologyFunctions_, *haloMassFunction_, *cosmologyParameters_, *cosmologicalMassVariance_"/>
    !!]

    self%hasSubhalos       %isSet=.false.
    self%massesAreInclusive%isSet=.false.
    self%areSelfContained  %isSet=.false.
    self%includesHubbleFlow%isSet=.false.
    self%periodicPositions %isSet=.false.
    self%length            %isSet=.false.
    self%forestIndicesRead       =.false.
    return
  end function galacticusConstructorInternal

  subroutine galacticusDestructor(self)
    !!{
    Destructor for the \glc\ format merger tree importer class.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    type(mergerTreeImporterGalacticus), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyFunctions_"      />
    <objectDestructor name="self%haloMassFunction_"        />
    <objectDestructor name="self%cosmologyParameters_"     />
    <objectDestructor name="self%cosmologicalMassVariance_"/>
    !!]
    !$ call hdf5Access%set()
    if (self%forestHalos%isOpen()) call self%forestHalos%close()
    if (self%file       %isOpen()) call self%file      %close()
    !$ call hdf5Access%unset()
    return
  end subroutine galacticusDestructor

  subroutine galacticusOpen(self,fileName)
    !!{
    Validate a \glc\ format merger tree file.
    !!}
    use :: Cosmology_Parameters, only : hubbleUnitsLittleH
    use :: Display             , only : displayMessage    , verbosityLevelWarn, displayMagenta, displayReset
    use :: Error               , only : Error_Report      , Warn
    use :: HDF5_Access         , only : hdf5Access
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    class           (mergerTreeImporterGalacticus ), intent(inout) :: self
    type            (varying_string               ), intent(in   ) :: fileName
    type            (hdf5Object                   )                :: cosmologicalParametersGroup, unitsGroup, angularMomentumDataset, spinDataset
    type            (varying_string               )                :: message
    character       (len=14                       )                :: valueString
    double precision                                               :: localLittleH0, localOmegaMatter, localOmegaDE, localOmegaBaryon, localSigma8, cosmologicalParameter

    ! Reset initialization status.
    self%hasSubhalos       %isSet=.false.
    self%massesAreInclusive%isSet=.false.
    self%areSelfContained  %isSet=.false.
    self%includesHubbleFlow%isSet=.false.
    self%periodicPositions %isSet=.false.
    self%length            %isSet=.false.
    self%forestIndicesRead       =.false.
    ! Get cosmological parameters. We do this in advance to avoid HDF5 thread conflicts.
    localLittleH0   =self%cosmologyParameters_     %HubbleConstant (hubbleUnitsLittleH)
    localOmegaMatter=self%cosmologyParameters_     %OmegaMatter    (                  )
    localOmegaDE    =self%cosmologyParameters_     %OmegaDarkEnergy(                  )
    localOmegaBaryon=self%cosmologyParameters_     %OmegaBaryon    (                  )
    localSigma8     =self%cosmologicalMassVariance_%sigma8         (                  )
    !$ call hdf5Access%set()
    ! Open the file.
    call self%file%openFile(char(fileName),readOnly=.true.)
    ! Get the file format version number.
    if (self%file%hasAttribute('formatVersion')) then
       call self%file%readAttribute('formatVersion',self%formatVersion,allowPseudoScalar=.true.)
    else
       self%formatVersion=1
    end if
    ! Validate format version number and set appropriate group names.
    select case (self%formatVersion)
    case (1)
       ! This version will be deprecated.
       !![
       <expiry version="1.0.0"/>
       !!]
       call Warn(displayMagenta()//'WARNING:'//displayReset()//' merger tree file format version is outdated - this format will soon be deprecated')
       self%forestHalosGroupName          ='haloTrees'
       self%forestContainmentAttributeName='treesAreSelfContained'
       self%forestIndexGroupName          ='treeIndex'
       self%forestIndexDatasetName        ='treeIndex'
       self%forestWeightDatasetName       ='treeWeight'
    case (2)
       ! This is the current version - no problems.
       self%forestHalosGroupName          ='forestHalos'
       self%forestContainmentAttributeName='forestsAreSelfContained'
       self%forestIndexGroupName          ='forestIndex'
       self%forestIndexDatasetName        ='forestIndex'
       self%forestWeightDatasetName       ='forestWeight'
    case default
       ! Unknown version - abort.
       call Error_Report('unknown file format version number'//{introspection:location})
    end select
    ! Open the merger trees group.
    self%forestHalos=self%file%openGroup(trim(self%forestHalosGroupName))
    ! Check that cosmological parameters are consistent with the internal ones.
    cosmologicalParametersGroup=self%file%openGroup("cosmology")
    if (cosmologicalParametersGroup%hasAttribute("OmegaMatter")) then
       call cosmologicalParametersGroup%readAttribute("OmegaMatter",cosmologicalParameter,allowPseudoScalar=.true.)
       if (Values_Differ(cosmologicalParameter,localOmegaMatter,absTol=0.001d0)) then
          message='Omega_Matter in merger tree file ['
          write (valueString,'(e14.8)') cosmologicalParameter
          message=message//trim(valueString)//'] differs from the internal value ['
          write (valueString,'(e14.8)') localOmegaMatter
          message=message//trim(valueString)//']'
          if (self%fatalMismatches) then
             call Error_Report(message//{introspection:location})
          else
             call displayMessage(message,verbosityLevelWarn)
          end if
       end if
    end if
    if (cosmologicalParametersGroup%hasAttribute("OmegaBaryon")) then
       call cosmologicalParametersGroup%readAttribute("OmegaBaryon",cosmologicalParameter,allowPseudoScalar=.true.)
       if (Values_Differ(cosmologicalParameter,localOmegaBaryon,absTol=0.001d0)) then
          message='Omega_b in merger tree file ['
          write (valueString,'(e14.8)') cosmologicalParameter
          message=message//trim(valueString)//'] differs from the internal value ['
          write (valueString,'(e14.8)') localOmegaBaryon
          message=message//trim(valueString)//']'
          if (self%fatalMismatches) then
             call Error_Report(message//{introspection:location})
          else
             call displayMessage(message,verbosityLevelWarn)
          end if
       end if
    end if
    if (cosmologicalParametersGroup%hasAttribute("OmegaLambda")) then
       call cosmologicalParametersGroup%readAttribute("OmegaLambda",cosmologicalParameter,allowPseudoScalar=.true.)
       if (Values_Differ(cosmologicalParameter,localOmegaDE,absTol=0.001d0)) then
          message='Omega_DE in merger tree file ['
          write (valueString,'(e14.8)') cosmologicalParameter
          message=message//trim(valueString)//'] differs from the internal value ['
          write (valueString,'(e14.8)') localOmegaDE
          message=message//trim(valueString)//']'
          if (self%fatalMismatches) then
             call Error_Report(message//{introspection:location})
          else
             call displayMessage(message,verbosityLevelWarn)
          end if
       end if
    end if
    if (cosmologicalParametersGroup%hasAttribute("HubbleParam")) then
       call cosmologicalParametersGroup%readAttribute("HubbleParam",cosmologicalParameter,allowPseudoScalar=.true.)
       if (Values_Differ(cosmologicalParameter,localLittleH0,absTol=0.00001d0)) then
          message='H_0 in merger tree file ['
          write (valueString,'(e14.8)') cosmologicalParameter
          message=message//trim(valueString)//'] differs from the internal value ['
          write (valueString,'(e14.8)') localLittleH0
          message=message//trim(valueString)//']'
          if (self%fatalMismatches) then
             call Error_Report(message//{introspection:location})
          else
             call displayMessage(message,verbosityLevelWarn)
          end if
       end if
    end if
    if (cosmologicalParametersGroup%hasAttribute("sigma_8")) then
       call cosmologicalParametersGroup%readAttribute("sigma_8",cosmologicalParameter,allowPseudoScalar=.true.)
       if (Values_Differ(cosmologicalParameter,localSigma8,absTol=0.00001d0)) then
          message='sigma_8 in merger tree file ['
          write (valueString,'(e14.8)') cosmologicalParameter
          message=message//trim(valueString)//'] differs from the internal value ['
          write (valueString,'(e14.8)') localSigma8
          message=message//trim(valueString)//'] - may not matter if sigma_8 is not used in other functions'
          call displayMessage(message)
       end if
    end if
    call cosmologicalParametersGroup%close()
    ! Read units.
    unitsGroup=self%file%openGroup("units")
    self%massUnit    %status=unitsGroup%hasAttribute("massUnitsInSI")
    if (self%massUnit    %status) then
       call unitsGroup%readAttribute("massUnitsInSI"              ,self%massUnit    %unitsInSI          ,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("massScaleFactorExponent"    ,self%massUnit    %scaleFactorExponent,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("massHubbleExponent"         ,self%massUnit    %hubbleExponent     ,allowPseudoScalar=.true.)
    end if
    self%lengthUnit  %status=unitsGroup%hasAttribute("lengthUnitsInSI")
    if (self%lengthUnit  %status) then
       call unitsGroup%readAttribute("lengthUnitsInSI"            ,self%lengthUnit  %unitsInSI          ,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("lengthScaleFactorExponent"  ,self%lengthUnit  %scaleFactorExponent,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("lengthHubbleExponent"       ,self%lengthUnit  %hubbleExponent     ,allowPseudoScalar=.true.)
    end if
    self%timeUnit    %status=unitsGroup%hasAttribute("timeUnitsInSI")
    if (self%timeUnit    %status) then
       call unitsGroup%readAttribute("timeUnitsInSI"              ,self%timeUnit    %unitsInSI          ,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("timeScaleFactorExponent"    ,self%timeUnit    %scaleFactorExponent,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("timeHubbleExponent"         ,self%timeUnit    %hubbleExponent     ,allowPseudoScalar=.true.)
    end if
    self%velocityUnit%status=unitsGroup%hasAttribute("velocityUnitsInSI")
    if (self%velocityUnit%status) then
       call unitsGroup%readAttribute("velocityUnitsInSI"          ,self%velocityUnit%unitsInSI          ,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("velocityScaleFactorExponent",self%velocityUnit%scaleFactorExponent,allowPseudoScalar=.true.)
       call unitsGroup%readAttribute("velocityHubbleExponent"     ,self%velocityUnit%hubbleExponent     ,allowPseudoScalar=.true.)
    end if
    call unitsGroup%close()
    ! Check for availability of particle data.
    if (self%file%hasGroup("particles")) then
       self%particles=self%file%openGroup("particles")
       if      (self%particles%hasDataset("time"           )) then
          self%particleEpochType       =galacticusParticleEpochTypeTime
          self%particleEpochDatasetName="time"
       else if (self%particles%hasDataset("expansionFactor")) then
          self%particleEpochType       =galacticusParticleEpochTypeExpansionFactor
          self%particleEpochDatasetName="expansionFactor"
       else if (self%particles%hasDataset("redshift"       )) then
          self%particleEpochType       =galacticusParticleEpochTypeRedshift
          self%particleEpochDatasetName="redshift"
       else
          call Error_Report("particles group must have one of time, redshift or expansionFactor datasets"//{introspection:location})
       end if
    end if
    ! Check for type of angular momenta data available.
    self%angularMomentaIsScalar=.false.
    self%angularMomentaIsVector=.false.
    if (self%forestHalos%hasDataset("angularMomentum")) then
       angularMomentumDataset=self%forestHalos%openDataset("angularMomentum")
       select case (angularMomentumDataset%rank())
       case (1)
          self%angularMomentaIsScalar=.true.
       case (2)
          if (angularMomentumDataset%size(1) /= 3) call Error_Report('2nd dimension of rank-2 angularMomentum dataset must be 3'//{introspection:location})
          self%angularMomentaIsVector=.true.
       case default
          call Error_Report('angularMomentum dataset must be rank 1 or 2'//{introspection:location})
       end select
       call angularMomentumDataset%close()
    end if
    ! Check for type of spin data available.
    self%spinIsScalar=.false.
    self%spinIsVector=.false.
    if (self%forestHalos%hasDataset("spin")) then
       spinDataset=self%forestHalos%openDataset("spin")
       select case (spinDataset%rank())
       case (1)
          self%spinIsScalar=.true.
       case (2)
          if (spinDataset%size(1) /= 3) call Error_Report('2nd dimension of rank-2 spin dataset must be 3'//{introspection:location})
          self%spinIsVector=.true.
       case default
          call Error_Report('spin dataset must be rank 1 or 2'//{introspection:location})
       end select
       call spinDataset%close()
    end if
    !$ call hdf5Access%unset()
    return
  end subroutine galacticusOpen

  subroutine galacticusClose(self)
    !!{
    Validate a \glc\ format merger tree file.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$ call hdf5Access%set()
    if (self%particles  %isOpen()) call self%particles  %close()
    if (self%forestHalos%isOpen()) call self%forestHalos%close()
    if (self%file       %isOpen()) call self%file       %close()
    !$ call hdf5Access%unset()
    return
  end subroutine galacticusClose

  logical function galacticusCanReadSubsets(self)
    !!{
    Return true since this format does permit reading of arbitrary subsets of halos from a forest.
    !!}
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    !$GLC attributes unused :: self

    galacticusCanReadSubsets=.true.
    return
  end function galacticusCanReadSubsets

  integer function galacticusTreesHaveSubhalos(self)
    !!{
    Return a Boolean integer specifying whether or not the trees have subhalos.
    !!}
    use :: HDF5_Access                , only : hdf5Access
    use :: Numerical_Constants_Boolean, only : booleanUnknown
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%hasSubhalos%isSet) then
       !$ call hdf5Access%set()
       if (self%forestHalos%hasAttribute("treesHaveSubhalos")) then
          call self%forestHalos%readAttribute("treesHaveSubhalos",self%hasSubhalos%value,allowPseudoScalar=.true.)
       else
          self%hasSubhalos%value=booleanUnknown
       end if
       !$ call hdf5Access%unset()
       self%hasSubhalos%isSet=.true.
    end if
    galacticusTreesHaveSubhalos=self%hasSubhalos%value
    return
  end function galacticusTreesHaveSubhalos

  logical function galacticusMassesIncludeSubhalos(self)
    !!{
    Return a Boolean specifying whether or not the halo masses include the contribution from subhalos.
    !!}
    use :: Error      , only : Error_Report
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                                              :: haloMassesIncludeSubhalosInteger

    if (.not.self%massesAreInclusive%isSet) then
       !$ call hdf5Access%set()
       if (self%forestHalos%hasAttribute("haloMassesIncludeSubhalos")) then
          call self%forestHalos%readAttribute("haloMassesIncludeSubhalos",haloMassesIncludeSubhalosInteger,allowPseudoScalar=.true.)
          self%massesAreInclusive%value=(haloMassesIncludeSubhalosInteger == 1)
       else
          call Error_Report('required attribute "haloMassesIncludeSubhalos" not present'//{introspection:location})
       end if
       !$ call hdf5Access%unset()
       self%massesAreInclusive%isSet=.true.
    end if
    galacticusMassesIncludeSubhalos=self%massesAreInclusive%value
    return
  end function galacticusMassesIncludeSubhalos

  logical function galacticusAngularMomentaIncludeSubhalos(self)
    !!{
    Return a Boolean specifying whether or not the halo momenta include the contribution from subhalos.
    !!}
    use :: Error      , only : Error_Report
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                                              :: haloAngularMomentaIncludeSubhalosInteger
    logical                                              :: attributeExists

    if (.not.self%angularMomentaAreInclusive%isSet) then
       !$ call hdf5Access%set()
       attributeExists=self%forestHalos%hasAttribute("haloAngularMomentaIncludeSubhalos")
       !$ call hdf5Access%unset()
       if (attributeExists) then
          !$ call hdf5Access%set()
          call self%forestHalos%readAttribute("haloAngularMomentaIncludeSubhalos",haloAngularMomentaIncludeSubhalosInteger,allowPseudoScalar=.true.)
          !$ call hdf5Access%unset()
          self%angularMomentaAreInclusive%value=(haloAngularMomentaIncludeSubhalosInteger == 1)
       else
          self%angularMomentaAreInclusive%value=self%massesIncludeSubhalos()
       end if
       self%angularMomentaAreInclusive%isSet=.true.
    end if
    galacticusAngularMomentaIncludeSubhalos=self%angularMomentaAreInclusive%value
    return
  end function galacticusAngularMomentaIncludeSubhalos

  integer function galacticusTreesAreSelfContained(self)
    !!{
    Return a Boolean integer specifying whether or not the trees are self-contained.
    !!}
    use :: HDF5_Access                , only : hdf5Access
    use :: Numerical_Constants_Boolean, only : booleanUnknown
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%areSelfContained%isSet) then
       !$ call hdf5Access%set()
       if (self%forestHalos%hasAttribute(trim(self%forestContainmentAttributeName))) then
          call self%forestHalos%readAttribute(trim(self%forestContainmentAttributeName),self%areSelfContained%value,allowPseudoScalar=.true.)
       else
          self%areSelfContained%value=booleanUnknown
       end if
       !$ call hdf5Access%unset()
       self%areSelfContained%isSet=.true.
    end if
    galacticusTreesAreSelfContained=self%areSelfContained%value
    return
  end function galacticusTreesAreSelfContained

  integer function galacticusVelocitiesIncludeHubbleFlow(self)
    !!{
    Return a Boolean integer specifying whether or not velocities include the Hubble flow.
    !!}
    use :: HDF5_Access                , only : hdf5Access
    use :: Numerical_Constants_Boolean, only : booleanUnknown
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%includesHubbleFlow%isSet) then
       !$ call hdf5Access%set()
       if (self%forestHalos%hasAttribute("velocitiesIncludeHubbleFlow")) then
          call self%forestHalos%readAttribute("velocitiesIncludeHubbleFlow",self%includesHubbleFlow%value,allowPseudoScalar=.true.)
       else
          self%includesHubbleFlow%value=booleanUnknown
       end if
       !$ call hdf5Access%unset()
       self%includesHubbleFlow%isSet=.true.
    end if
    galacticusVelocitiesIncludeHubbleFlow=self%includesHubbleFlow%value
    return
  end function galacticusVelocitiesIncludeHubbleFlow

  integer function galacticusPositionsArePeriodic(self)
    !!{
    Return a Boolean integer specifying whether or not positions are periodic.
    !!}
    use :: HDF5_Access                , only : hdf5Access
    use :: Numerical_Constants_Boolean, only : booleanUnknown
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%periodicPositions%isSet) then
       !$ call hdf5Access%set()
       if (self%forestHalos%hasAttribute("positionsArePeriodic")) then
          call self%forestHalos%readAttribute("positionsArePeriodic",self%periodicPositions%value,allowPseudoScalar=.true.)
       else
          self%periodicPositions%value=booleanUnknown
       end if
       !$ call hdf5Access%unset()
       self%periodicPositions%isSet=.true.
    end if
    galacticusPositionsArePeriodic=self%periodicPositions%value
    return
  end function galacticusPositionsArePeriodic

  double precision function galacticusCubeLength(self,time,status)
    !!{
    Return the length of the simulation cube.
    !!}
    use :: Error                           , only : Error_Report
    use :: HDF5_Access                     , only : hdf5Access
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Boolean     , only : booleanFalse, booleanTrue, booleanUnknown
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)           :: self
    double precision                              , intent(in   )           :: time
    integer                                       , intent(  out), optional :: status
    type            (hdf5Object                  )                          :: simulationGroup

    if (.not.self%length%isSet) then
       !$ call hdf5Access%set()
       if (self%file%hasGroup("simulation")) then
          simulationGroup=self%file%openGroup("simulation")
          if (simulationGroup%hasAttribute("boxSize")) then
             call simulationGroup%readAttribute("boxSize",self%length%value,allowPseudoScalar=.true.)
             if (self%length%value <= 0.0d0) call Error_Report('simulation box length must be positive'//{introspection:location})
             self%lengthStatus%value=booleanTrue
          else
             self%lengthStatus%value=booleanFalse
          end if
          call simulationGroup%close()
       else
          self%lengthStatus%value=booleanUnknown
       end if
       !$ call hdf5Access%unset()
       self%length      %isSet=.true.
       self%lengthStatus%isSet=.true.
    end if
    if (self%lengthStatus%value == booleanTrue) then
       galacticusCubeLength=importerUnitConvert(self%length%value,time,self%lengthUnit,megaParsec,self%cosmologyParameters_,self%cosmologyFunctions_)
    else
       galacticusCubeLength=0.0d0
    end if
    if (present(status)) then
       status=self%lengthStatus%value
    else
       if (self%lengthStatus%value == booleanFalse) call Error_Report('the boxSize attribute of the simulation group is required'//{introspection:location})
    end if
    return
  end function galacticusCubeLength

  function galacticusTreeCount(self)
    !!{
    Return a count of the number of trees available.
    !!}
    implicit none
    integer(c_size_t                    )                :: galacticusTreeCount
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    call galacticusForestIndicesRead(self)
    galacticusTreeCount=self%forestsCount
    return
  end function galacticusTreeCount

  integer(kind=kind_int8) function galacticusTreeIndex(self,i)
    !!{
    Return the index of the $i^\mathrm{th}$ tree.
    !!}
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call galacticusForestIndicesRead(self)
    galacticusTreeIndex=self%forestIndices(i)
    return
  end function galacticusTreeIndex

  function galacticusNodeCount(self,i)
    !!{
    Return a count of the number of nodes in the $i^\mathrm{th}$ tree.
    !!}
    implicit none
    integer(c_size_t                    )                :: galacticusNodeCount
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call galacticusForestIndicesRead(self)
    galacticusNodeCount=self%nodeCounts(i)
    return
  end function galacticusNodeCount

  subroutine galacticusForestIndicesRead(self)
    !!{
    Read the tree indices.
    !!}
    use :: Error                           , only : Error_Report
    use :: HDF5                            , only : HSIZE_T
    use :: HDF5_Access                     , only : hdf5Access
    use :: Numerical_Constants_Astronomical, only : gigaYear
    use :: Sorting                         , only : sortIndex
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)             :: self
    type            (hdf5Object                  )                            :: treeIndexGroup
    integer         (kind=kind_int8              ), allocatable, dimension(:) :: descendantIndex
    double precision                              , allocatable, dimension(:) :: nodeMass           , treeMass    , &
         &                                                                       nodeTime           , treeTime
    integer         (kind=HSIZE_T                )             , dimension(1) :: firstNodeIndex     , nodeCount
    integer         (kind=c_size_t               ), allocatable, dimension(:) :: sortOrder
    integer                                                                   :: i
    integer         (c_size_t                    )                            :: iNode
    double precision                                                          :: massMinimum        , massMaximum
    logical                                                                   :: hasForestWeights

    if (self%forestIndicesRead) return
    !$ call hdf5Access%set()
    if (.not.self%file%hasGroup(trim(self%forestIndexGroupName))) &
         & call Error_Report('merger tree file must contain the treeIndex group'//{introspection:location})
    treeIndexGroup=self%file%openGroup(trim(self%forestIndexGroupName))
    call treeIndexGroup%readDataset("firstNode"                      ,self%firstNodes   )
    call treeIndexGroup%readDataset("numberOfNodes"                  ,self%nodeCounts   )
    call treeIndexGroup%readDataset(trim(self%forestIndexDatasetName),self%forestIndices)
    hasForestWeights=treeIndexGroup%hasDataset(trim(self%forestWeightDatasetName))
    !$ call hdf5Access%unset()
    if (self%reweightTrees) then
       allocate(self%weights(size(self%firstNodes)))
       allocate(treeMass    (size(self%firstNodes)))
       allocate(treeTime    (size(self%firstNodes)))
       !$ call hdf5Access%set()
       do i=1,size(self%firstNodes)
          firstNodeIndex(1)=self%firstNodes(i)+1
          nodeCount     (1)=self%nodeCounts(i)
          ! Allocate the nodes array.
          allocate(descendantIndex(nodeCount(1)))
          allocate(nodeMass       (nodeCount(1)))
          allocate(nodeTime       (nodeCount(1)))
          call self%forestHalos%readDatasetStatic("descendantIndex",descendantIndex,firstNodeIndex,nodeCount)
          call self%forestHalos%readDatasetStatic("nodeMass"       ,nodeMass       ,firstNodeIndex,nodeCount)
          if      (self%forestHalos%hasDataset("time"           )) then
             ! Time is present, so read it.
             call self%forestHalos%readDatasetStatic("time"           ,nodeTime,firstNodeIndex,nodeCount)
             nodeTime=importerUnitConvert(nodeTime,nodeTime,self%timeUnit,gigaYear,self%cosmologyParameters_,self%cosmologyFunctions_)
          else if (self%forestHalos%hasDataset("expansionFactor")) then
             ! Expansion factor is present, read it instead.
             call self%forestHalos%readDatasetStatic("expansionFactor",nodeTime,firstNodeIndex,nodeCount)
             ! Convert expansion factors to times.
             do iNode=1,nodeCount(1)
                nodeTime(iNode)=self%cosmologyFunctions_%cosmicTime(nodeTime(iNode))
             end do
          else if (self%forestHalos%hasDataset("redshift"       )) then
             ! Redshift is present, read it instead.
             call self%forestHalos%readDatasetStatic("redshift"       ,nodeTime,firstNodeIndex,nodeCount)
             ! Convert redshifts to times.
             do iNode=1,nodeCount(1)
                nodeTime(iNode)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(nodeTime(iNode)))
             end do
          else
             call Error_Report("one of time, redshift or expansionFactor data sets must be present in forestHalos group"//{introspection:location})
          end if
          if (count(descendantIndex == -1) /= 1) call Error_Report('reweighting trees requires there to be only only root node'//{introspection:location})
          treeMass(i)=sum(nodeMass,mask=descendantIndex == -1)
          treeTime(i)=sum(nodeTime,mask=descendantIndex == -1)
          deallocate(descendantIndex)
          deallocate(nodeMass       )
          deallocate(nodeTime       )
       end do
       !$ call hdf5Access%unset()
       ! Sort the trees into mass order.
       sortOrder=sortIndex(treeMass)
       ! Abort if there is only one tree.
       if (size(self%firstNodes) <= 1) call Error_Report('reweighting trees requires there to be at least two forests'//{introspection:location})
       ! Compute the weight for each tree.
       do i=1,size(self%firstNodes)
          ! Get the minimum mass of the interval occupied by this tree.
          if (i == 1) then
             massMinimum=treeMass(sortOrder(i))*sqrt(treeMass(sortOrder(i))/treeMass(sortOrder(i+1)))
          else
             massMinimum=sqrt(treeMass(sortOrder(i))*treeMass(sortOrder(i-1)))
          end if
          ! Get the maximum mass of the interval occupied by this tree.
          if (i == size(self%firstNodes)) then
             massMaximum=treeMass(sortOrder(i))*sqrt(treeMass(sortOrder(i))/treeMass(sortOrder(i-1)))
          else
             massMaximum=sqrt(treeMass(sortOrder(i))*treeMass(sortOrder(i+1)))
          end if
          ! Get the integral of the halo mass function over this range.
          self%weights(sortOrder(i))=self%haloMassFunction_%integrated(treeTime(sortOrder(i)),massMinimum,massMaximum)
       end do
       !$ call hdf5Access%set()
       call treeIndexGroup%readDatasetStatic(trim(self%forestWeightDatasetName),self%weights)
       !$ call hdf5Access%unset()
       deallocate(treeMass)
       deallocate(treeTime)
    else if (hasForestWeights) then
       !$ call hdf5Access%set()
       call treeIndexGroup%readDataset(trim(self%forestWeightDatasetName),self%weights)
       !$ call hdf5Access%unset()
    end if
    !$ call hdf5Access%set()
    call treeIndexGroup%close()
    !$ call hdf5Access%unset()
    self%forestsCount=size(self%forestIndices)
    ! Reset first node indices to Fortran array standard.
    self%firstNodes=self%firstNodes+1
    self%forestIndicesRead=.true.
    return
  end subroutine galacticusForestIndicesRead

  double precision function galacticusTreeWeight(self,i)
    !!{
    Return the weight to assign to trees.
    !!}
    use :: Error                           , only : Error_Report
    use :: Numerical_Constants_Astronomical, only : megaParsec
    use :: Numerical_Constants_Boolean     , only : booleanTrue
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                                       , intent(in   ) :: i
    double precision                                              :: lengthSimulationBox, timePresent
    integer                                                       :: statusActual

    call galacticusForestIndicesRead(self)
    ! Determine the time at present.
    timePresent=self%cosmologyFunctions_%cosmicTime(1.0d0)
    ! Do we have an array of weights for trees?
    if (allocated(self%weights)) then
       ! We do, so simply return the appropriate weight.
       galacticusTreeWeight=importerUnitConvert(self%weights(i),timePresent,self%lengthUnit**(-3),1.0d0/megaParsec**3,self%cosmologyParameters_,self%cosmologyFunctions_)
    else
       ! We do not, so attempt to find the volume of the simulation cube.
       lengthSimulationBox=self%cubeLength(timePresent,statusActual)
       if (statusActual == booleanTrue) then
          ! Simulation cube length found, compute the inverse volume.
          galacticusTreeWeight=1.0d0/lengthSimulationBox**3
       else
          ! No method exists to determine the weight. Return unity.
          galacticusTreeWeight=1.0d0
       end if
    end if
    return
  end function galacticusTreeWeight

  logical function galacticusPositionsAvailable(self,positions,velocities)
    !!{
    Return true if positions and/or velocities are available.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    logical                              , intent(in   ) :: positions, velocities

    galacticusPositionsAvailable=.true.
    !$ call hdf5Access%set()
    if (positions .and..not.self%forestHalos%hasDataset("position")) galacticusPositionsAvailable=.false.
    if (velocities.and..not.self%forestHalos%hasDataset("velocity")) galacticusPositionsAvailable=.false.
    !$ call hdf5Access%unset()
    return
  end function galacticusPositionsAvailable

  logical function galacticusScaleRadiiAvailable(self)
    !!{
    Return true if scale radii are available.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$ call hdf5Access%set()
    galacticusScaleRadiiAvailable=                        &
         &  self%forestHalos%hasDataset("halfMassRadius") &
         & .or.                                           &
         &  self%forestHalos%hasDataset("scaleRadius"   )
    !$ call hdf5Access%unset()
    return
  end function galacticusScaleRadiiAvailable

  logical function galacticusParticleCountAvailable(self)
    !!{
    Return true if particle counts are available.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$ call hdf5Access%set()
    galacticusParticleCountAvailable=self%forestHalos%hasDataset("particleCount")
    !$ call hdf5Access%unset()
    return
  end function galacticusParticleCountAvailable

  logical function galacticusVelocityMaximumAvailable(self)
    !!{
    Return true if halo rotation curve velocity maxima are available.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$ call hdf5Access%set()
    galacticusVelocityMaximumAvailable=self%forestHalos%hasDataset("velocityMaximum")
    !$ call hdf5Access%unset()
    return
  end function galacticusVelocityMaximumAvailable

  logical function galacticusVelocityDispersionAvailable(self)
    !!{
    Return true if halo velocity dispersions are available.
    !!}
    use :: HDF5_Access, only : hdf5Access
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$ call hdf5Access%set()
    galacticusVelocityDispersionAvailable=self%forestHalos%hasDataset("velocityDispersion")
    !$ call hdf5Access%unset()
    return
  end function galacticusVelocityDispersionAvailable

  logical function galacticusAngularMomentaAvailable(self)
    !!{
    Return true if angular momenta are available.
    !!}
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    galacticusAngularMomentaAvailable=self%angularMomentaIsScalar.or.self%angularMomentaIsVector
    return
  end function galacticusAngularMomentaAvailable

  logical function galacticusAngularMomenta3DAvailable(self)
    !!{
    Return true if angular momenta vectors are available.
    !!}
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    galacticusAngularMomenta3DAvailable=self%angularMomentaIsVector
    return
  end function galacticusAngularMomenta3DAvailable

  logical function galacticusSpinAvailable(self)
    !!{
    Return true if spins are available.
    !!}
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    galacticusSpinAvailable=self%spinIsScalar.or.self%spinIsVector
    return
  end function galacticusSpinAvailable

  logical function galacticusSpin3DAvailable(self)
    !!{
    Return true if spins vectors are available.
    !!}
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    galacticusSpin3DAvailable=self%spinIsVector
    return
  end function galacticusSpin3DAvailable

  subroutine galacticusSubhaloTrace(self,node,time,position,velocity)
    !!{
    Returns a trace of subhalo position/velocity.
    !!}
    use :: Error                           , only : Error_Report
    use :: HDF5_Access                     , only : hdf5Access
    use :: Numerical_Constants_Astronomical, only : gigaYear    , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)                 :: self
    class           (nodeData                    ), intent(in   )                 :: node
    double precision                              , intent(  out), dimension(:  ) :: time
    double precision                              , intent(  out), dimension(:,:) :: position, velocity
    integer                                                                       :: i

    select type (node)
    type is (nodeDataGalacticus)
       ! Read epoch, position, and velocity data.
       !$ call hdf5Access%set()
       call self%particles%readDatasetStatic(char(self%particleEpochDatasetName),time    ,[            node%particleIndexStart+1],[            node%particleIndexCount])
       call self%particles%readDatasetStatic("position"                         ,position,[1_kind_int8,node%particleIndexStart+1],[3_kind_int8,node%particleIndexCount])
       call self%particles%readDatasetStatic("velocity"                         ,velocity,[1_kind_int8,node%particleIndexStart+1],[3_kind_int8,node%particleIndexCount])
       !$ call hdf5Access%unset()
       ! Convert epochs into times.
       select case (self%particleEpochType%ID)
       case (galacticusParticleEpochTypeTime           %ID)
          time=importerUnitConvert(time,time,self%timeUnit,gigaYear,self%cosmologyParameters_,self%cosmologyFunctions_)
       case (galacticusParticleEpochTypeExpansionFactor%ID)
          do i=1,size(time)
             time(i)=self%cosmologyFunctions_%cosmicTime(                                                     time(i) )
          end do
       case (galacticusParticleEpochTypeRedshift       %ID)
          do i=1,size(time)
             time(i)=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(time(i)))
          end do
       end select
       ! Convert units of position and velocity into Galacticus internal units.
       position=importerUnitConvert(position,time,self%lengthUnit  ,megaParsec,self%cosmologyParameters_,self%cosmologyFunctions_)
       velocity=importerUnitConvert(velocity,time,self%velocityUnit,kilo      ,self%cosmologyParameters_,self%cosmologyFunctions_)
    class default
       call Error_Report('node should be of type nodeDataGalacticus'//{introspection:location})
    end select
    return
  end subroutine galacticusSubhaloTrace

  function galacticusSubhaloTraceCount(self,node)
    !!{
    Returns the length of a subhalo trace.
    !!}
    use :: Error, only : Error_Report
    implicit none
    integer(c_size_t                    )                :: galacticusSubhaloTraceCount
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    class  (nodeData                    ), intent(in   ) :: node
    !$GLC attributes unused :: self

    select type (node)
    type is (nodeDataGalacticus)
       galacticusSubhaloTraceCount=node%particleIndexCount
    class default
       galacticusSubhaloTraceCount=0
       call Error_Report('node should be of type nodeDataGalacticus'//{introspection:location})
    end select
    return
  end function galacticusSubhaloTraceCount

  subroutine galacticusImport(self,i,nodes,nodeSubset,requireScaleRadii,requireAngularMomenta,requireAngularMomenta3D,requireSpin,requireSpin3D,requirePositions,structureOnly,requireNamedReals,requireNamedIntegers)
    !!{
    Import the $i^\mathrm{th}$ merger tree.
    !!}
    use :: Error                           , only : Error_Report       , Warn
    use :: HDF5                            , only : hsize_t
    use :: HDF5_Access                     , only : hdf5Access
    use :: Numerical_Constants_Astronomical, only : gigaYear           , massSolar      , megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)                              :: self
    integer                                       , intent(in   )                              :: i
    class           (nodeDataMinimal             ), intent(  out), allocatable, dimension(:  ) :: nodes
    integer         (c_size_t                    ), intent(in   ), optional   , dimension(:  ) :: nodeSubset
    logical                                       , intent(in   ), optional                    :: requireScaleRadii      , requireAngularMomenta, &
         &                                                                                        requireAngularMomenta3D, requirePositions     , &
         &                                                                                        structureOnly          , requireSpin          , &
         &                                                                                        requireSpin3D
    type            (varying_string              ), intent(in   ), optional   , dimension(:  ) :: requireNamedReals      , requireNamedIntegers
    integer         (hsize_t                     )                            , dimension(1  ) :: firstNodeIndex         , nodeCount
    integer         (c_size_t                    )                                             :: iNode
    integer         (c_size_t                    )               , allocatable, dimension(:  ) :: nodeSubsetOffset
    double precision                                             , allocatable, dimension(:,:) :: angularMomentum3D      , position             , &
         &                                                                                        velocity               , spin3D
    double precision                                             , allocatable, dimension(:  ) :: namedReal
    integer         (kind_int8                   )               , allocatable, dimension(:  ) :: namedInteger
    integer                                                                                    :: j
    logical                                                                                    :: timesAreInternal       , useNodeSubset

    ! Ensure tree indices have been read.
    call galacticusForestIndicesRead(self)
    ! Determine the first node index and the node count.
    firstNodeIndex(1)=self%firstNodes(i)
    nodeCount     (1)=self%nodeCounts(i)
    ! Handle node subsets.
    useNodeSubset=present(nodeSubset) .and. .not.nodeSubset(1) < 0
    if (useNodeSubset) then
       ! Check that all nodes are within range.
       if (any(nodeSubset > nodeCount(1))) call Error_Report('node subset lies outside of forest'//{introspection:location})
       ! Shift node subset to start of this forest.
       allocate(nodeSubsetOffset(size(nodeSubset)))
       nodeSubsetOffset=nodeSubset+firstNodeIndex(1)-1
       ! Reset size of node array to read.
       nodeCount(1)=size(nodeSubset)
    end if
    ! Allocate the nodes array.
    if (present(structureOnly).and.structureOnly) then
       allocate(nodeDataMinimal    :: nodes(nodeCount(1)))
    else
       allocate(nodeDataGalacticus :: nodes(nodeCount(1)))
    end if
    !$ call hdf5Access%set()
    if (useNodeSubset) then
       ! nodeIndex, hostIndex, parentNode
       call self%forestHalos%readDatasetStatic("nodeIndex"      ,nodes%nodeIndex                               ,readSelection=nodeSubsetOffset)
       ! hostIndex
       call self%forestHalos%readDatasetStatic("hostIndex"      ,nodes%hostIndex                               ,readSelection=nodeSubsetOffset)
       ! parentNode
       call self%forestHalos%readDatasetStatic("descendantIndex",nodes%descendantIndex                         ,readSelection=nodeSubsetOffset)
    else
       ! nodeIndex
       call self%forestHalos%readDatasetStatic("nodeIndex"      ,nodes%nodeIndex      ,firstNodeIndex,nodeCount                               )
       ! hostIndex
       call self%forestHalos%readDatasetStatic("hostIndex"      ,nodes%hostIndex      ,firstNodeIndex,nodeCount                               )
       ! parentNode
       call self%forestHalos%readDatasetStatic("descendantIndex",nodes%descendantIndex,firstNodeIndex,nodeCount                               )
    end if
    ! nodeTime
    timesAreInternal=.true.
    if      (self%forestHalos%hasDataset("time"           )) then
       ! Time is present, so read it.
       timesAreInternal=.false.
       if (useNodeSubset) then
          call self%forestHalos%readDatasetStatic("time"           ,nodes%nodeTime                         ,readSelection=nodeSubsetOffset)
       else
          call self%forestHalos%readDatasetStatic("time"           ,nodes%nodeTime,firstNodeIndex,nodeCount                               )
       end if
       if (self%validateData) then
          if (any(isNaN(nodes%nodeTime)        )) call Error_Report('"nodeTime" dataset contains NaN values'    //{introspection:location})
          if (any(      nodes%nodeTime <= 0.0d0)) call Error_Report('"nodeTime" dataset has non-positive values'//{introspection:location})
       end if
    else if (self%forestHalos%hasDataset("expansionFactor")) then
       ! Expansion factor is present, read it instead.
       if (useNodeSubset) then
          call self%forestHalos%readDatasetStatic("expansionFactor",nodes%nodeTime                         ,readSelection=nodeSubsetOffset)
       else
          call self%forestHalos%readDatasetStatic("expansionFactor",nodes%nodeTime,firstNodeIndex,nodeCount                               )
       end if
       ! Validate expansion factors.
       if (self%validateData) then
          if (any(isNaN(nodes%nodeTime)        )) call Error_Report('"expansionFactor" dataset contains NaN values'//{introspection:location})
          if (any(      nodes%nodeTime <= 0.0d0)) call Error_Report('"expansionFactor" dataset values must be >0'  //{introspection:location})
          if (any(      nodes%nodeTime >  1.0d0)) call Warn        ("WARNING: some expansion factors are in the future when importing merger tree")
       end if
       ! Convert expansion factors to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=self%cosmologyFunctions_%cosmicTime(nodes(iNode)%nodeTime)
       end do
    else if (self%forestHalos%hasDataset("redshift"       )) then
       ! Redshift is present, read it instead.
       if (useNodeSubset) then
          call self%forestHalos%readDatasetStatic("redshift"       ,nodes%nodeTime                         ,readSelection=nodeSubsetOffset)
       else
          call self%forestHalos%readDatasetStatic("redshift"       ,nodes%nodeTime,firstNodeIndex,nodeCount                               )
       end if
       ! Validate redshifts.
       if (self%validateData) then
          if (any(isNaN(nodes%nodeTime)         )) call Error_Report('"redshift" dataset contains NaN values'//{introspection:location})
          if (any(      nodes%nodeTime <= -1.0d0)) call Error_Report('"redshift" dataset values must be >-1' //{introspection:location})
          if (any(      nodes%nodeTime <   0.0d0)) call Warn        ("WARNING: some redshifts are in the future when importing merger tree")
       end if
       ! Convert redshifts to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(nodes(iNode)%nodeTime))
       end do
    else
       call Error_Report("one of time, redshift or expansionFactor data sets must be present in forestHalos group"//{introspection:location})
    end if
    ! nodeMass
    if (useNodeSubset) then
       call self%forestHalos%readDatasetStatic("nodeMass"       ,nodes%nodeMass                                ,readSelection=nodeSubsetOffset)
    else
       call self%forestHalos%readDatasetStatic("nodeMass"       ,nodes%nodeMass       ,firstNodeIndex,nodeCount                               )
    end if
    if (self%validateData) then
       if (any(isNaN(nodes%nodeMass))) call Error_Report('"nodeMass" dataset contains NaN values'//{introspection:location})
    end if
    !$ call hdf5Access%unset()
    ! If only structure is requested we are done.
    if (present(structureOnly).and.structureOnly) return
    select type (nodes)
    type is (nodeDataGalacticus)
       !$ call hdf5Access%set()
       ! Scale or half-mass radius.
       if (present(requireScaleRadii).and.requireScaleRadii) then
          if (self%forestHalos%hasDataset("scaleRadius")) then
             nodes%halfMassRadius=-1.0d0
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("scaleRadius"   ,nodes%scaleRadius                            ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("scaleRadius"   ,nodes%scaleRadius   ,firstNodeIndex,nodeCount                               )
             end if
             if (self%validateData) then
                if (any(isNaN(nodes%scaleRadius))) call Error_Report('"scaleRadius" dataset contains NaN values'//{introspection:location})
             end if
          else
             nodes%scaleRadius   =-1.0d0
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("halfMassRadius",nodes%halfMassRadius                         ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("halfMassRadius",nodes%halfMassRadius,firstNodeIndex,nodeCount                               )
             end if
             if (self%validateData) then
                if (any(isNaN(nodes%halfMassRadius))) call Error_Report('"halfMassRadius" dataset contains NaN values'//{introspection:location})
             end if
          end if
       end if
       ! Halo angular momenta.
       if (present(requireAngularMomenta).and.requireAngularMomenta) then
          if (self%angularMomentaIsVector) then
             if (useNodeSubset) then
                call self%forestHalos%readDataset("angularMomentum",angularMomentum3D                                                         ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDataset("angularMomentum",angularMomentum3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
             end if
             if (self%validateData) then
                if (any(isNaN(angularMomentum3D))) call Error_Report('"angularMomentum" dataset contains NaN values'//{introspection:location})
             end if
             ! Transfer to nodes.
             forall(iNode=1:nodeCount(1))
                nodes(iNode)%angularMomentum=Vector_Magnitude(angularMomentum3D(:,iNode))
             end forall
             deallocate(angularMomentum3D)
          else if (self%angularMomentaIsScalar) then
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("angularMomentum",nodes%angularMomentum                                   ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("angularMomentum",nodes%angularMomentum,[firstNodeIndex(1)],[nodeCount(1)]                               )
             end if
             if (self%validateData) then
                if (any(isNaN(nodes%angularMomentum))) call Error_Report('"angularMomentum" dataset contains NaN values'//{introspection:location})
             end if
          else
             call Error_Report("scalar angular momentum is not available"//{introspection:location})
          end if
       end if
       if (present(requireAngularMomenta3D).and.requireAngularMomenta3D) then
          if (.not.self%angularMomentaIsVector) call Error_Report("vector angular momentum is not available"//{introspection:location})
          if (useNodeSubset) then
             call self%forestHalos%readDataset("angularMomentum",angularMomentum3D                                                         ,readSelection=nodeSubsetOffset)
          else
             call self%forestHalos%readDataset("angularMomentum",angularMomentum3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
          end if
          if (self%validateData) then
             if (any(isNaN(angularMomentum3D))) call Error_Report('"angularMomentum" dataset contains NaN values'//{introspection:location})
          end if
       end if
       ! Halo spins.
       if (present(requireSpin).and.requireSpin) then
          if (self%spinIsVector) then
             if (useNodeSubset) then
                call self%forestHalos%readDataset("spin",spin3D                                                         ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDataset("spin",spin3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
             end if
             if (self%validateData) then
                if (any(isNaN(spin3D))) call Error_Report('"spin" dataset contains NaN values'//{introspection:location})
             end if
             ! Transfer to nodes.
             forall(iNode=1:nodeCount(1))
                nodes(iNode)%spin=Vector_Magnitude(spin3D(:,iNode))
             end forall
          deallocate(spin3D)
          else if (self%spinIsScalar) then
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("spin",nodes%spin                                   ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("spin",nodes%spin,[firstNodeIndex(1)],[nodeCount(1)]                               )
             end if
             if (self%validateData) then
                if (any(isNaN(nodes%spin))) call Error_Report('"spin" dataset contains NaN values'//{introspection:location})
             end if
          else
             call Error_Report("scalar spin is not available"//{introspection:location})
          end if
       end if
       if (present(requireSpin3D).and.requireSpin3D) then
          if (.not.self%spinIsVector) call Error_Report("vector spin is not available"//{introspection:location})
          if (useNodeSubset) then
             call self%forestHalos%readDataset("spin",spin3D                                                         ,readSelection=nodeSubsetOffset)
          else
             call self%forestHalos%readDataset("spin",spin3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
          end if
          if (self%validateData) then
             if (any(isNaN(spin3D))) call Error_Report('"spin" dataset contains NaN values'//{introspection:location})
          end if
          ! Transfer to nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%spin3D=spin3D(:,iNode)
          end forall
          deallocate(spin3D)
       end if
       ! Read arbitrary named real datasets.
       if (present(requireNamedReals   )) then
          do iNode=1,nodeCount(1)
             allocate(nodes(iNode)%reals   (size(requireNamedReals   )))
          end do
          do j=1,size(requireNamedReals   )
             if (.not.self%forestHalos%hasDataset(char(requireNamedReals   (j)))) &
                  & call Error_Report('named dataset "'//char(requireNamedReals   (j))//'" is not available'//{introspection:location})
             if (useNodeSubset) then
                call self%forestHalos%readDataset(char(requireNamedReals   (j)),namedReal                                      ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDataset(char(requireNamedReals   (j)),namedReal   ,[firstNodeIndex(1)],[nodeCount(1)]                               )
             end if
             ! Transfer to nodes.
             forall(iNode=1:nodeCount(1))
                nodes(iNode)%reals   (j)=namedReal   (iNode)
             end forall
             deallocate(namedReal   )
          end do
       end if
       ! Read arbitrary named integer datasets.
       if (present(requireNamedIntegers)) then
          do iNode=1,nodeCount(1)
             allocate(nodes(iNode)%integers(size(requireNamedIntegers)))
          end do
          do j=1,size(requireNamedIntegers)
             if (.not.self%forestHalos%hasDataset(char(requireNamedIntegers(j)))) &
                  & call Error_Report('named dataset "'//char(requireNamedIntegers(j))//'" is not available'//{introspection:location})
             if (useNodeSubset) then
                call self%forestHalos%readDataset(char(requireNamedIntegers(j)),namedInteger                                   ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDataset(char(requireNamedIntegers(j)),namedInteger,[firstNodeIndex(1)],[nodeCount(1)]                               )
             end if
             ! Transfer to nodes.
             forall(iNode=1:nodeCount(1))
                nodes(iNode)%integers(j)=namedInteger(iNode)
             end forall
             deallocate(namedInteger)
          end do
       end if
       ! Initialize particle data to null values.
       nodes%particleIndexStart=-1_c_size_t
       nodes%particleIndexCount=-1_c_size_t
       ! Positions (and velocities).
       if (present(requirePositions).and.requirePositions) then
          if (useNodeSubset) then
             ! position.
             call self%forestHalos%readDataset("position",position                                                         ,readSelection=nodeSubsetOffset)
             ! velocity.
             call self%forestHalos%readDataset("velocity",velocity                                                         ,readSelection=nodeSubsetOffset)
          else
             ! position.
             call self%forestHalos%readDataset("position",position,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
             ! velocity.
             call self%forestHalos%readDataset("velocity",velocity,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
          end if
          if (self%validateData) then
             if (any(isNaN(position))) call Error_Report('"position" dataset contains NaN values'//{introspection:location})
             if (any(isNaN(velocity))) call Error_Report('"velocity" dataset contains NaN values'//{introspection:location})
          end if
          ! If a set of most bound particle indices are present, read them.
          if (self%forestHalos%hasDataset("particleIndexStart").and.self%forestHalos%hasDataset("particleIndexCount")) then
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("particleIndexStart",nodes%particleIndexStart                         ,readSelection=nodeSubsetOffset)
                call self%forestHalos%readDatasetStatic("particleIndexCount",nodes%particleIndexCount                         ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("particleIndexStart",nodes%particleIndexStart,firstNodeIndex,nodeCount                               )
                call self%forestHalos%readDatasetStatic("particleIndexCount",nodes%particleIndexCount,firstNodeIndex,nodeCount                               )
             end if
          end if
       end if
       !$ call hdf5Access%unset()
       ! Unit conversion.
       if     (                                                                   &
            &   present(requirePositions)                                         &
            &  .and.                                                              &
            &           requirePositions                                          &
            &  .and.                                                              &
            &   .not.                                                             &
            &    (                                                                &
            &      self%lengthUnit%status                                         &
            &     .and.                                                           &
            &      self%velocityUnit%status                                       &
            &    )                                                                &
            & ) call Error_Report(                                                &
            &                     'length and velocity units must be given if '// &
            &                     'positions and velocities are to be read'    // &
            &                     {introspection:location}                        &
            &                    )
       if (self%timeUnit%status.and.self%timeUnit%scaleFactorExponent /= 0)                 &
            &   call Error_Report(                                                          &
            &                     'expect no scaling of time units with expansion factor'// &
            &                     {introspection:location}                                  &
            &                    )
       if (self%    massUnit%status.and.self%    massUnit%unitsInSI <= 0.0d0) call Error_Report('non-positive units for mass'    //{introspection:location})
       if (self%  lengthUnit%status.and.self%  lengthUnit%unitsInSI <= 0.0d0) call Error_Report('non-positive units for length'  //{introspection:location})
       if (self%velocityUnit%status.and.self%velocityUnit%unitsInSI <= 0.0d0) call Error_Report('non-positive units for velocity'//{introspection:location})
       if (self%    timeUnit%status.and.self%    timeUnit%unitsInSI <= 0.0d0) call Error_Report('non-positive units for time'    //{introspection:location})
       if (.not.timesAreInternal)                                                                                                                                         &
            & nodes%nodeTime       =importerUnitConvert(nodes%nodeTime       ,nodes%nodeTime,self%timeUnit                                  ,gigaYear                 ,self%cosmologyParameters_,self%cosmologyFunctions_)
       nodes       %nodeMass       =importerUnitConvert(nodes%nodeMass       ,nodes%nodeTime,                                  self%massUnit,                massSolar,self%cosmologyParameters_,self%cosmologyFunctions_)
       if (present(requireScaleRadii).and.requireScaleRadii) then
          nodes    %scaleRadius    =importerUnitConvert(nodes%scaleRadius    ,nodes%nodeTime,self%lengthUnit                                ,megaParsec               ,self%cosmologyParameters_,self%cosmologyFunctions_)
          nodes    %halfMassRadius =importerUnitConvert(nodes%halfMassRadius ,nodes%nodeTime,self%lengthUnit                                ,megaParsec               ,self%cosmologyParameters_,self%cosmologyFunctions_)
       end if
       if (present(requireAngularMomenta     ).and.requireAngularMomenta     )                                                                                                  &
            & nodes%angularMomentum   =importerUnitConvert(nodes%angularMomentum   ,nodes%nodeTime,self%lengthUnit*self%velocityUnit*self%massUnit,megaParsec*kilo*massSolar,self%cosmologyParameters_,self%cosmologyFunctions_)
       if (present(requireAngularMomenta3D).and.requireAngularMomenta3D) then
          angularmomentum3d=importerUnitConvert(angularmomentum3d,nodes%nodeTime,self%lengthUnit*self%velocityUnit*self%massUnit,megaParsec*kilo*massSolar,self%cosmologyParameters_,self%cosmologyFunctions_)
          ! Transfer to nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%angularMomentum3D=angularMomentum3D(:,iNode)
          end forall
       deallocate(angularMomentum3D)
       end if
       if (present(requirePositions).and.requirePositions) then
          position=importerUnitConvert(position,nodes%nodeTime,self%  lengthUnit,megaParsec,self%cosmologyParameters_,self%cosmologyFunctions_)
          velocity=importerUnitConvert(velocity,nodes%nodeTime,self%velocityUnit,kilo      ,self%cosmologyParameters_,self%cosmologyFunctions_)
          ! Transfer to the nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%position=position(:,iNode)
             nodes(iNode)%velocity=velocity(:,iNode)
          end forall
          deallocate(position)
          deallocate(velocity)
       end if
    end select
    return
  end subroutine galacticusImport
