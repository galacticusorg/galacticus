!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of the merger tree importer class for \glc\ format merger tree files.

  !# <mergerTreeImporter name="mergerTreeImporterGalacticus" description="Importer for \glc\ format merger tree files." />
  use IO_HDF5
  use Stateful_Types

  type, public, extends(nodeData) :: nodeDataGalacticus
     !% Extension of the {\normalfont \ttfamily nodeData} class for \glc\ format merger trees. Stores particle indices and counts for nodes.
     integer(c_size_t) :: particleIndexCount, particleIndexStart
  end type nodeDataGalacticus

  type, extends(mergerTreeImporterClass) :: mergerTreeImporterGalacticus
     !% A merger tree importer class for \glc\ format merger tree files.
     private
     type   (hdf5Object     )                            :: file                  , forestHalos
     type   (statefulInteger)                            :: hasSubhalos           , areSelfContained          , &
          &                                                 includesHubbleFlow    , periodicPositions         , &
          &                                                 lengthStatus
     type   (statefulLogical)                            :: massesAreInclusive    , angularMomentaAreInclusive
     type   (statefulDouble )                            :: length
     type   (importerUnits  )                            :: massUnit              , lengthUnit                , &
          &                                                 timeUnit              , velocityUnit
     logical                                             :: fatalMismatches       , forestIndicesRead         , &
          &                                                 angularMomentaIsScalar, angularMomentaIsVector    , &
          &                                                 spinIsScalar          , spinIsVector              , &
          &                                                 reweightTrees
     integer                                             :: forestsCount          , formatVersion
     integer                 , allocatable, dimension(:) :: firstNodes            , nodeCounts
     integer(kind=kind_int8 ), allocatable, dimension(:) :: forestIndices
     double precision        , allocatable, dimension(:) :: weights
     type   (hdf5Object     )                            :: particles
     integer                                             :: particleEpochType
     type   (varying_string )                            :: particleEpochDataSetName
     character(len=32)                                   :: forestHalosGroupName    , forestContainmentAttributeName, &
          &                                                 forestIndexGroupName    , forestIndexDatasetName        , &
          &                                                 forestWeightDatasetName
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
     !% Constructors for the \glc\ format merger tree importer class.
     module procedure galacticusDefaultConstructor
  end interface mergerTreeImporterGalacticus

  ! Record of implementation initialization state.
  logical            :: galacticusInitialized                     =.false.

  ! Default settings.
  logical            :: mergerTreeImportGalacticusMismatchIsFatal, mergerTreeImportGalacticusReweightTrees

  ! Particle epoch enumeration.
  integer, parameter :: galacticusParticleEpochTypeTime           =0
  integer, parameter :: galacticusParticleEpochTypeExpansionFactor=1
  integer, parameter :: galacticusParticleEpochTypeRedshift       =2

contains

  function galacticusDefaultConstructor()
    !% Default constructor for the \glc\ format merger tree importer.
    use Input_Parameters
    implicit none
    type (mergerTreeImporterGalacticus), target :: galacticusDefaultConstructor
    
    if (.not.galacticusInitialized) then
       !$omp critical (mergerTreeImporterGalacticusInitialize)
       if (.not.galacticusInitialized) then
          !@ <inputParameter>
          !@   <name>mergerTreeImportGalacticusMismatchIsFatal</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>true</defaultValue>
          !@   <description>
          !@     Specifies whether mismatches in cosmological parameter values between \glc\ and the merger tree file should be considered fatal.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportGalacticusMismatchIsFatal',mergerTreeImportGalacticusMismatchIsFatal,defaultValue=.true.)
          !@ <inputParameter>
          !@   <name>mergerTreeImportGalacticusReweightTrees</name>
          !@   <attachedTo>module</attachedTo>
          !@   <defaultValue>false</defaultValue>
          !@   <description>
          !@     Specifies whether merger tree weights should be recomputed from the halo mass function.
          !@   </description>
          !@   <type>boolean</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('mergerTreeImportGalacticusReweightTrees',mergerTreeImportGalacticusReweightTrees,defaultValue=.false.)
          galacticusInitialized=.true.
       end if
       !$omp end critical (mergerTreeImporterGalacticusInitialize)
    end if
    galacticusDefaultConstructor%hasSubhalos       %isSet=.false.
    galacticusDefaultConstructor%massesAreInclusive%isSet=.false.
    galacticusDefaultConstructor%areSelfContained  %isSet=.false.
    galacticusDefaultConstructor%includesHubbleFlow%isSet=.false.
    galacticusDefaultConstructor%periodicPositions %isSet=.false.
    galacticusDefaultConstructor%length            %isSet=.false.
    galacticusDefaultConstructor%forestIndicesRead       =.false.
    galacticusDefaultConstructor%fatalMismatches         =mergerTreeImportGalacticusMismatchIsFatal
    galacticusDefaultConstructor%reweightTrees           =mergerTreeImportGalacticusReweightTrees
    return
  end function galacticusDefaultConstructor

  subroutine galacticusDestructor(self)
    !% Destructor for the \glc\ format merger tree importer class.
    implicit none
    type(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    call self%file%close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine galacticusDestructor

  subroutine galacticusOpen(self,fileName)
    !% Validate a \glc\ format merger tree file.
    use Numerical_Comparison
    use Galacticus_Display
    use Galacticus_Error
    use Cosmological_Mass_Variance
    use Cosmology_Parameters
    implicit none
    class           (mergerTreeImporterGalacticus ), intent(inout) :: self
    type            (varying_string               ), intent(in   ) :: fileName
    class           (cosmologyParametersClass     ), pointer       :: cosmologyParameters_
    class           (cosmologicalMassVarianceClass), pointer       :: cosmologicalMassVariance_
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
    ! Get the default cosmology.
    cosmologyParameters_      => cosmologyParameters     ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    ! Get cosmological parameters. We do this in advance to avoid HDF5 thread conflicts.
    localLittleH0   =cosmologyParameters_     %HubbleConstant (hubbleUnitsLittleH)
    localOmegaMatter=cosmologyParameters_     %OmegaMatter    (                  )
    localOmegaDE    =cosmologyParameters_     %OmegaDarkEnergy(                  )
    localOmegaBaryon=cosmologyParameters_     %OmegaBaryon    (                  )
    localSigma8     =cosmologicalMassVariance_%sigma8         (                  )
    !$omp critical(HDF5_Access)
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
       !# <expiry version="1.0.0"/>
       call Galacticus_Warn('WARNING: merger tree file format version is outdated - this format will soon be deprecated')
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
       call Galacticus_Error_Report('galacticusOpen','unknown file format version number')
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
             call Galacticus_Error_Report('galacticusOpen',message)
          else
             call Galacticus_Display_Message(message,verbosityWarn)
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
             call Galacticus_Error_Report('galacticusOpen',message)
          else
             call Galacticus_Display_Message(message,verbosityWarn)
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
             call Galacticus_Error_Report('galacticusOpen',message)
          else
             call Galacticus_Display_Message(message,verbosityWarn)
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
             call Galacticus_Error_Report('galacticusOpen',message)
          else
             call Galacticus_Display_Message(message,verbosityWarn)
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
          message=message//trim(valueString)//'] - may not matter if sigma_8 is not use in other functions'
          call Galacticus_Display_Message(message)
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
          call Galacticus_Error_Report("galacticusOpen","particles group must have one of time, redshift or expansionFactor datasets")
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
          if (angularMomentumDataset%size(1) /= 3) call Galacticus_Error_Report('galacticusOpen','2nd dimension of rank-2 angularMomentum dataset must be 3')
          self%angularMomentaIsVector=.true.
       case default
          call Galacticus_Error_Report('galacticusOpen','angularMomentum dataset must be rank 1 or 2')
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
          if (spinDataset%size(1) /= 3) call Galacticus_Error_Report('galacticusOpen','2nd dimension of rank-2 spin dataset must be 3')
          self%spinIsVector=.true.
       case default
          call Galacticus_Error_Report('galacticusOpen','spin dataset must be rank 1 or 2')
       end select
       call spinDataset%close()
    end if
    !$omp end critical(HDF5_Access)
    return
  end subroutine galacticusOpen

  subroutine galacticusClose(self)
    !% Validate a \glc\ format merger tree file.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    
    !$omp critical(HDF5_Access)
    if (self%particles  %isOpen()) call self%particles  %close()
    if (self%forestHalos%isOpen()) call self%forestHalos%close()
    if (self%file       %isOpen()) call self%file       %close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine galacticusClose

  logical function galacticusCanReadSubsets(self)
    !% Return true since this format does permit reading of arbitrary subsets of halos from a forest.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    !GCC$ attributes unused :: self

    galacticusCanReadSubsets=.true.
    return
  end function galacticusCanReadSubsets

  integer function galacticusTreesHaveSubhalos(self)
    !% Return a Boolean integer specifying whether or not the trees have subhalos.
    use Numerical_Constants_Boolean
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%hasSubhalos%isSet) then
       !$omp critical(HDF5_Access)
       if (self%forestHalos%hasAttribute("treesHaveSubhalos")) then
          call self%forestHalos%readAttribute("treesHaveSubhalos",self%hasSubhalos%value,allowPseudoScalar=.true.)
       else
          self%hasSubhalos%value=booleanUnknown
       end if
       !$omp end critical(HDF5_Access)
       self%hasSubhalos%isSet=.true.
    end if
    galacticusTreesHaveSubhalos=self%hasSubhalos%value
    return
  end function galacticusTreesHaveSubhalos

  logical function galacticusMassesIncludeSubhalos(self)
    !% Return a Boolean specifying whether or not the halo masses include the contribution from subhalos.
    use Galacticus_Error
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                                              :: haloMassesIncludeSubhalosInteger

    if (.not.self%massesAreInclusive%isSet) then
       !$omp critical(HDF5_Access)
       if (self%forestHalos%hasAttribute("haloMassesIncludeSubhalos")) then
          call self%forestHalos%readAttribute("haloMassesIncludeSubhalos",haloMassesIncludeSubhalosInteger,allowPseudoScalar=.true.)
          self%massesAreInclusive%value=(haloMassesIncludeSubhalosInteger == 1)
       else
          call Galacticus_Error_Report('galacticusMassesIncludeSubhalos','required attribute "haloMassesIncludeSubhalos" not present')
       end if
       !$omp end critical(HDF5_Access)
       self%massesAreInclusive%isSet=.true.
    end if
    galacticusMassesIncludeSubhalos=self%massesAreInclusive%value
    return
  end function galacticusMassesIncludeSubhalos

  logical function galacticusAngularMomentaIncludeSubhalos(self)
    !% Return a Boolean specifying whether or not the halo momenta include the contribution from subhalos.
    use Galacticus_Error
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                                              :: haloAngularMomentaIncludeSubhalosInteger
    logical                                              :: attributeExists

    if (.not.self%angularMomentaAreInclusive%isSet) then
       !$omp critical(HDF5_Access)
       attributeExists=self%forestHalos%hasAttribute("haloAngularMomentaIncludeSubhalos")
       !$omp end critical(HDF5_Access)
       if (attributeExists) then
          !$omp critical(HDF5_Access)
          call self%forestHalos%readAttribute("haloAngularMomentaIncludeSubhalos",haloAngularMomentaIncludeSubhalosInteger,allowPseudoScalar=.true.)
          !$omp end critical(HDF5_Access)
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
    !% Return a Boolean integer specifying whether or not the trees are self-contained.
    use Numerical_Constants_Boolean
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%areSelfContained%isSet) then
       !$omp critical(HDF5_Access)
       if (self%forestHalos%hasAttribute(trim(self%forestContainmentAttributeName))) then
          call self%forestHalos%readAttribute(trim(self%forestContainmentAttributeName),self%areSelfContained%value,allowPseudoScalar=.true.)
       else
          self%areSelfContained%value=booleanUnknown
       end if
       !$omp end critical(HDF5_Access)
       self%areSelfContained%isSet=.true.
    end if
    galacticusTreesAreSelfContained=self%areSelfContained%value
    return
  end function galacticusTreesAreSelfContained

  integer function galacticusVelocitiesIncludeHubbleFlow(self)
    !% Return a Boolean integer specifying whether or not velocities include the Hubble flow.
    use Numerical_Constants_Boolean
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%includesHubbleFlow%isSet) then
       !$omp critical(HDF5_Access)
       if (self%forestHalos%hasAttribute("velocitiesIncludeHubbleFlow")) then
          call self%forestHalos%readAttribute("velocitiesIncludeHubbleFlow",self%includesHubbleFlow%value,allowPseudoScalar=.true.)
       else
          self%includesHubbleFlow%value=booleanUnknown
       end if
       !$omp end critical(HDF5_Access)
       self%includesHubbleFlow%isSet=.true.
    end if
    galacticusVelocitiesIncludeHubbleFlow=self%includesHubbleFlow%value
    return
  end function galacticusVelocitiesIncludeHubbleFlow

  integer function galacticusPositionsArePeriodic(self)
    !% Return a Boolean integer specifying whether or not positions are periodic.
    use Numerical_Constants_Boolean
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%periodicPositions%isSet) then
       !$omp critical(HDF5_Access)
       if (self%forestHalos%hasAttribute("positionsArePeriodic")) then
          call self%forestHalos%readAttribute("positionsArePeriodic",self%periodicPositions%value,allowPseudoScalar=.true.)
       else
          self%periodicPositions%value=booleanUnknown
       end if
       !$omp end critical(HDF5_Access)
       self%periodicPositions%isSet=.true.
    end if
    galacticusPositionsArePeriodic=self%periodicPositions%value
    return
  end function galacticusPositionsArePeriodic

  double precision function galacticusCubeLength(self,time,status)
    !% Return the length of the simulation cube.
    use Numerical_Constants_Boolean
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)           :: self
    double precision                              , intent(in   )           :: time
    integer                                       , intent(  out), optional :: status
    type            (hdf5Object                  )                          :: simulationGroup

    if (.not.self%length%isSet) then
       !$omp critical(HDF5_Access)
       if (self%file%hasGroup("simulation")) then
          simulationGroup=self%file%openGroup("simulation")
          if (simulationGroup%hasAttribute("boxSize")) then
             call simulationGroup%readAttribute("boxSize",self%length%value,allowPseudoScalar=.true.)
             if (self%length%value <= 0.0d0) call Galacticus_Error_Report('galacticusCubeLength','simulation box length must be positive')
             self%lengthStatus%value=booleanTrue
          else
             self%lengthStatus%value=booleanFalse
          end if
          call simulationGroup%close()
       else
          self%lengthStatus%value=booleanUnknown
       end if
       !$omp end critical(HDF5_Access)   
       self%length      %isSet=.true.
       self%lengthStatus%isSet=.true.
    end if
    if (self%lengthStatus%value == booleanTrue) then
       galacticusCubeLength=importerUnitConvert(self%length%value,time,self%lengthUnit,megaParsec)
    else
       galacticusCubeLength=0.0d0
    end if
    if (present(status)) then
       status=self%lengthStatus%value
    else
       if (self%lengthStatus%value == booleanFalse) call Galacticus_Error_Report('galacticusCubeLength','the boxSize attribute of the simulation group is required')
    end if
    return
  end function galacticusCubeLength

  function galacticusTreeCount(self)
    !% Return a count of the number of trees available.
    implicit none
    integer(c_size_t                    )                :: galacticusTreeCount
    class  (mergerTreeImporterGalacticus), intent(inout) :: self

    call galacticusForestIndicesRead(self)
    galacticusTreeCount=self%forestsCount
    return
  end function galacticusTreeCount

  integer(kind=kind_int8) function galacticusTreeIndex(self,i)
    !% Return the index of the $i^{\mathrm th}$ tree.
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call galacticusForestIndicesRead(self)
    galacticusTreeIndex=self%forestIndices(i)
    return
  end function galacticusTreeIndex

  function galacticusNodeCount(self,i)
    !% Return a count of the number of nodes in the $i^{\mathrm th}$ tree.
    implicit none
    integer(c_size_t                    )                :: galacticusNodeCount
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call galacticusForestIndicesRead(self)
    galacticusNodeCount=self%nodeCounts(i)
    return
  end function galacticusNodeCount

  subroutine galacticusForestIndicesRead(self)
    !% Read the tree indices.
    use Galacticus_Error
    use HDF5
    use Sort
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    use Halo_Mass_Functions
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)             :: self
    type            (hdf5Object                  )                            :: treeIndexGroup
    integer         (kind=kind_int8              ), allocatable, dimension(:) :: descendentIndex
    double precision                              , allocatable, dimension(:) :: nodeMass                 , treeMass    , &
         &                                                                       nodeTime                 , treeTime
    integer         (kind=HSIZE_T                )             , dimension(1) :: firstNodeIndex           , nodeCount
    integer         (kind=c_size_t               ), allocatable, dimension(:) :: sortOrder
    class           (cosmologyFunctionsClass     ), pointer                   :: cosmologyFunctionsDefault
    class           (haloMassFunctionClass       ), pointer                   :: haloMassFunction_
    integer                                                                   :: i
    integer         (c_size_t                    )                            :: iNode
    double precision                                                          :: massMinimum              , massMaximum
    logical                                                                   :: hasForestWeights
    
    if (self%forestIndicesRead) return
    !$omp critical(HDF5_Access)
    if (.not.self%file%hasGroup(trim(self%forestIndexGroupName))) &
         & call Galacticus_Error_Report('galacticusForestIndicesRead','merger tree file must contain the treeIndex group')
    treeIndexGroup=self%file%openGroup(trim(self%forestIndexGroupName))
    call treeIndexGroup%readDataset("firstNode"                      ,self%firstNodes   )
    call treeIndexGroup%readDataset("numberOfNodes"                  ,self%nodeCounts   )
    call treeIndexGroup%readDataset(trim(self%forestIndexDatasetName),self%forestIndices)
    hasForestWeights=treeIndexGroup%hasDataset(trim(self%forestWeightDatasetName))
    !$omp end critical(HDF5_Access)
    if (self%reweightTrees) then
       cosmologyFunctionsDefault => cosmologyFunctions()
       haloMassFunction_         => haloMassFunction  ()
       allocate(self%weights(size(self%firstNodes)))
       allocate(treeMass    (size(self%firstNodes)))
       allocate(treeTime    (size(self%firstNodes)))
       !$omp critical(HDF5_Access)
       do i=1,size(self%firstNodes)
          firstNodeIndex(1)=self%firstNodes(i)+1
          nodeCount     (1)=self%nodeCounts(i)
          ! Allocate the nodes array.
          allocate(descendentIndex(nodeCount(1)))
          allocate(nodeMass       (nodeCount(1)))
          allocate(nodeTime       (nodeCount(1)))
          call self%forestHalos%readDatasetStatic("descendentIndex",descendentIndex,firstNodeIndex,nodeCount)
          call self%forestHalos%readDatasetStatic("nodeMass"       ,nodeMass       ,firstNodeIndex,nodeCount)
          if      (self%forestHalos%hasDataset("time"           )) then
             ! Time is present, so read it.
             call self%forestHalos%readDatasetStatic("time"           ,nodeTime,firstNodeIndex,nodeCount)
             nodeTime=importerUnitConvert(nodeTime,nodeTime,self%timeUnit,gigaYear)
          else if (self%forestHalos%hasDataset("expansionFactor")) then
             ! Expansion factor is present, read it instead.
             call self%forestHalos%readDatasetStatic("expansionFactor",nodeTime,firstNodeIndex,nodeCount)
             ! Convert expansion factors to times.
             do iNode=1,nodeCount(1)
                nodeTime(iNode)=cosmologyFunctionsDefault%cosmicTime(nodeTime(iNode))
             end do
          else if (self%forestHalos%hasDataset("redshift"       )) then
             ! Redshift is present, read it instead.
             call self%forestHalos%readDatasetStatic("redshift"       ,nodeTime,firstNodeIndex,nodeCount)
             ! Convert redshifts to times.
             do iNode=1,nodeCount(1)
                nodeTime(iNode)=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(nodeTime(iNode)))
             end do
          else
             call Galacticus_Error_Report("galacticusImport","one of time, redshift or expansionFactor data sets must be present in forestHalos group")
          end if
          if (count(descendentIndex == -1) /= 1) call Galacticus_Error_Report('galacticusForestIndicesRead','reweighting trees requires there to be only only root node')
          treeMass(i)=sum(nodeMass,mask=descendentIndex == -1)
          treeTime(i)=sum(nodeTime,mask=descendentIndex == -1)
          deallocate(descendentIndex)
          deallocate(nodeMass       )
          deallocate(nodeTime       )
       end do
       !$omp end critical(HDF5_Access)
       ! Sort the trees into mass order.
       sortOrder=Sort_Index_Do(treeMass)
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
          self%weights(sortOrder(i))=haloMassFunction_%integrated(treeTime(sortOrder(i)),massMinimum,massMaximum)
       end do
       !$omp critical(HDF5_Access)
       call treeIndexGroup%readDatasetStatic(trim(self%forestWeightDatasetName),self%weights)
       !$omp end critical(HDF5_Access)
       deallocate(treeMass)
       deallocate(treeTime)
    else if (hasForestWeights) then
       !$omp critical(HDF5_Access)
       call treeIndexGroup%readDataset(trim(self%forestWeightDatasetName),self%weights)
       !$omp endcritical(HDF5_Access)
    end if
    !$omp critical(HDF5_Access)
    call treeIndexGroup%close()
    !$omp end critical(HDF5_Access)
    self%forestsCount=size(self%forestIndices)
    ! Reset first node indices to Fortran array standard.
    self%firstNodes=self%firstNodes+1
    self%forestIndicesRead=.true.
    return
  end subroutine galacticusForestIndicesRead

  double precision function galacticusTreeWeight(self,i)
    !% Return the weight to assign to trees.
    use Numerical_Constants_Boolean
    use Numerical_Constants_Astronomical
    use Galacticus_Error
    use Cosmology_Functions
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                                       , intent(in   ) :: i
    class           (cosmologyFunctionsClass     ), pointer       :: cosmologyFunctionsDefault
    double precision                                              :: lengthSimulationBox      , timePresent
    integer                                                       :: statusActual

    call galacticusForestIndicesRead(self)
    ! Determine the time at present.
    cosmologyFunctionsDefault => cosmologyFunctions()
    timePresent=cosmologyFunctionsDefault%cosmicTime(1.0d0)
    ! Do we have an array of weights for trees?
    if (allocated(self%weights)) then
       ! We do, so simply return the appropriate weight.
       galacticusTreeWeight=importerUnitConvert(self%weights(i),timePresent,self%lengthUnit**(-3),1.0d0/megaParsec**3)
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
    !% Return true if positions and/or velocities are available.
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    logical                              , intent(in   ) :: positions, velocities

    galacticusPositionsAvailable=.true.
    !$omp critical(HDF5_Access)
    if (positions .and..not.self%forestHalos%hasDataset("position")) galacticusPositionsAvailable=.false.
    if (velocities.and..not.self%forestHalos%hasDataset("velocity")) galacticusPositionsAvailable=.false.
    !$omp end critical(HDF5_Access)
    return
  end function galacticusPositionsAvailable

  logical function galacticusScaleRadiiAvailable(self)
    !% Return true if scale radii are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    galacticusScaleRadiiAvailable=                        &
         &  self%forestHalos%hasDataset("halfMassRadius") &
         & .or.                                           &
         &  self%forestHalos%hasDataset("position"      )
    !$omp end critical(HDF5_Access)
    return
  end function galacticusScaleRadiiAvailable

  logical function galacticusParticleCountAvailable(self)
    !% Return true if particle counts are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    
    !$omp critical(HDF5_Access)
    galacticusParticleCountAvailable=self%forestHalos%hasDataset("particleCount")
    !$omp end critical(HDF5_Access)
    return
  end function galacticusParticleCountAvailable

  logical function galacticusVelocityMaximumAvailable(self)
    !% Return true if halo rotation curve velocity maxima are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    galacticusVelocityMaximumAvailable=self%forestHalos%hasDataset("velocityMaximum")
    !$omp end critical(HDF5_Access)
    return
  end function galacticusVelocityMaximumAvailable

  logical function galacticusVelocityDispersionAvailable(self)
    !% Return true if halo velocity dispersions are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    galacticusVelocityDispersionAvailable=self%forestHalos%hasDataset("velocityDispersion")
    !$omp end critical(HDF5_Access)
    return
  end function galacticusVelocityDispersionAvailable

  logical function galacticusAngularMomentaAvailable(self)
    !% Return true if angular momenta are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    
    galacticusAngularMomentaAvailable=self%angularMomentaIsScalar.or.self%angularMomentaIsVector
    return
  end function galacticusAngularMomentaAvailable

  logical function galacticusAngularMomenta3DAvailable(self)
    !% Return true if angular momenta vectors are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    
    galacticusAngularMomenta3DAvailable=self%angularMomentaIsVector
    return
  end function galacticusAngularMomenta3DAvailable

  logical function galacticusSpinAvailable(self)
    !% Return true if spins are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    galacticusSpinAvailable=self%spinIsScalar.or.self%spinIsVector
    return
  end function galacticusSpinAvailable

  logical function galacticusSpin3DAvailable(self)
    !% Return true if spins vectors are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    galacticusSpin3DAvailable=self%spinIsVector
    return
  end function galacticusSpin3DAvailable

  subroutine galacticusSubhaloTrace(self,node,time,position,velocity)
    !% Returns a trace of subhalo position/velocity.
    use Galacticus_Error
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)                 :: self
    class           (nodeData                    ), intent(in   )                 :: node
    double precision                              , intent(  out), dimension(:  ) :: time
    double precision                              , intent(  out), dimension(:,:) :: position, velocity
    class           (cosmologyFunctionsClass     ), pointer                       :: cosmologyFunctionsDefault
    integer                                                                       :: i
    
    select type (node)
    type is (nodeDataGalacticus)
       ! Read epoch, position, and velocity data.
       !$omp critical(HDF5_Access)
       call self%particles%readDatasetStatic(char(self%particleEpochDatasetName),time    ,[            node%particleIndexStart+1],[            node%particleIndexCount])
       call self%particles%readDatasetStatic("position"                         ,position,[1_kind_int8,node%particleIndexStart+1],[3_kind_int8,node%particleIndexCount])
       call self%particles%readDatasetStatic("velocity"                         ,velocity,[1_kind_int8,node%particleIndexStart+1],[3_kind_int8,node%particleIndexCount])
       !$omp end critical(HDF5_Access)
       ! Convert epochs into times.
       cosmologyFunctionsDefault => cosmologyFunctions()
       select case (self%particleEpochType)
       case (galacticusParticleEpochTypeTime           )
          time=importerUnitConvert(time,time,self%timeUnit,gigaYear)
       case (galacticusParticleEpochTypeExpansionFactor)
          do i=1,size(time)
             time(i)=cosmologyFunctionsDefault%cosmicTime(                                                      time(i) )
          end do
       case (galacticusParticleEpochTypeRedshift       )
          do i=1,size(time)
             time(i)=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(time(i)))
          end do
       end select
       ! Convert units of position and velocity into Galacticus internal units.
       position=importerUnitConvert(position,time,self%lengthUnit  ,megaParsec)
       velocity=importerUnitConvert(velocity,time,self%velocityUnit,kilo      )
    class default
       call Galacticus_Error_Report('galacticusSubhaloTrace','node should be of type nodeDataGalacticus')
    end select
    return
  end subroutine galacticusSubhaloTrace

  function galacticusSubhaloTraceCount(self,node)
    !% Returns the length of a subhalo trace.
    use Galacticus_Error
    implicit none
    integer(c_size_t                    )                :: galacticusSubhaloTraceCount
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    class  (nodeData                    ), intent(in   ) :: node
    !GCC$ attributes unused :: self
    
    select type (node)
    type is (nodeDataGalacticus)
       galacticusSubhaloTraceCount=node%particleIndexCount
    class default
       galacticusSubhaloTraceCount=0
       call Galacticus_Error_Report('galacticusSubhaloTraceCount','node should be of type nodeDataGalacticus')
    end select
    return
  end function galacticusSubhaloTraceCount

  subroutine galacticusImport(self,i,nodes,nodeSubset,requireScaleRadii,requireAngularMomenta,requireAngularMomenta3D,requireSpin,requireSpin3D,requirePositions,requireParticleCounts,requireVelocityMaxima,requireVelocityDispersions,structureOnly)
    !% Import the $i^{\mathrm th}$ merger tree.
    use Memory_Management
    use Cosmology_Functions
    use HDF5
    use Galacticus_Error
    use Galacticus_Display
    use Vectors
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Prefixes
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)                              :: self
    integer                                       , intent(in   )                              :: i
    class           (nodeDataMinimal             ), intent(  out), allocatable, dimension(:  ) :: nodes
    integer         (c_size_t                    ), intent(in   ), optional   , dimension(:  ) :: nodeSubset
    logical                                       , intent(in   ), optional                    :: requireScaleRadii         , requireAngularMomenta, &
         &                                                                                        requireAngularMomenta3D   , requirePositions     , &
         &                                                                                        requireParticleCounts     , requireVelocityMaxima, &
         &                                                                                        requireVelocityDispersions, requireSpin          , &
         &                                                                                        requireSpin3D             , structureOnly
    class           (cosmologyFunctionsClass     ), pointer                                    :: cosmologyFunctionsDefault
    integer         (kind=HSIZE_T                )                            , dimension(1  ) :: firstNodeIndex            , nodeCount
    integer         (c_size_t                    )                                             :: iNode
    integer         (c_size_t                    )               , allocatable, dimension(:  ) :: nodeSubsetOffset
    double precision                                             , allocatable, dimension(:,:) :: angularMomentum3D         , position             , &
         &                                                                                        velocity                  , spin3D
    logical                                                                                    :: timesAreInternal          , useNodeSubset

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Ensure tree indices have been read.
    call galacticusForestIndicesRead(self)
    ! Determine the first node index and the node count.
    firstNodeIndex(1)=self%firstNodes(i)
    nodeCount     (1)=self%nodeCounts(i)
    ! Handle node subsets.
    useNodeSubset=present(nodeSubset) .and. .not.nodeSubset(1) < 0
    if (useNodeSubset) then
       ! Check that all nodes are within range.
       if (any(nodeSubset > nodeCount(1))) call Galacticus_Error_Report('galacticusImport','node subset lies outside of forest')
       ! Shift node subset to start of this forest.
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
    !# <workaround type="gfortran" PR="65889" url="https://gcc.gnu.org/bugzilla/show_bug.cgi?id=65889">
    select type (nodes)
    type is (nodeDataMinimal   )
       call Memory_Usage_Record(sizeof(nodes))
    type is (nodeDataGalacticus)
       call Memory_Usage_Record(sizeof(nodes))
    end select
    !# </workaround>
    !$omp critical(HDF5_Access)
    if (useNodeSubset) then
       ! nodeIndex
       call self%forestHalos%readDatasetStatic("nodeIndex"      ,nodes%nodeIndex                               ,readSelection=nodeSubsetOffset)
       ! hostIndex
       call self%forestHalos%readDatasetStatic("hostIndex"      ,nodes%hostIndex                               ,readSelection=nodeSubsetOffset)
       ! parentNode
       call self%forestHalos%readDatasetStatic("descendentIndex",nodes%descendentIndex                         ,readSelection=nodeSubsetOffset)
    else
       ! nodeIndex
       call self%forestHalos%readDatasetStatic("nodeIndex"      ,nodes%nodeIndex      ,firstNodeIndex,nodeCount                               )
       ! hostIndex
       call self%forestHalos%readDatasetStatic("hostIndex"      ,nodes%hostIndex      ,firstNodeIndex,nodeCount                               )
       ! parentNode
       call self%forestHalos%readDatasetStatic("descendentIndex",nodes%descendentIndex,firstNodeIndex,nodeCount                               )
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
    else if (self%forestHalos%hasDataset("expansionFactor")) then
       ! Expansion factor is present, read it instead.
       if (useNodeSubset) then
          call self%forestHalos%readDatasetStatic("expansionFactor",nodes%nodeTime                         ,readSelection=nodeSubsetOffset)
       else
          call self%forestHalos%readDatasetStatic("expansionFactor",nodes%nodeTime,firstNodeIndex,nodeCount                               )
       end if
       ! Validate expansion factors.
       if (any(nodes%nodeTime <= 0.0d0)) call Galacticus_Error_Report("galacticusImport","expansionFactor dataset values must be >0")
       if (any(nodes%nodeTime >  1.0d0)) call Galacticus_Warn        ("WARNING: some expansion factors are in the future when importing merger tree")
       ! Convert expansion factors to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=cosmologyFunctionsDefault%cosmicTime(nodes(iNode)%nodeTime)
       end do
    else if (self%forestHalos%hasDataset("redshift"       )) then
       ! Redshift is present, read it instead.
       if (useNodeSubset) then
          call self%forestHalos%readDatasetStatic("redshift"       ,nodes%nodeTime                         ,readSelection=nodeSubsetOffset)
       else
          call self%forestHalos%readDatasetStatic("redshift"       ,nodes%nodeTime,firstNodeIndex,nodeCount                               )
       end if
       ! Validate redshifts.
       if (any(nodes%nodeTime <= -1.0d0)) call Galacticus_Error_Report("galacticusImport","redshift dataset values must be >-1")
       if (any(nodes%nodeTime <   0.0d0)) call Galacticus_Warn        ("WARNING: some redshifts are in the future when importing merger tree")
       ! Convert redshifts to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(nodes(iNode)%nodeTime))
       end do
    else
       call Galacticus_Error_Report("galacticusImport","one of time, redshift or expansionFactor data sets must be present in forestHalos group")
    end if
    ! nodeMass
    if (useNodeSubset) then
       call self%forestHalos%readDatasetStatic("nodeMass"       ,nodes%nodeMass                                ,readSelection=nodeSubsetOffset)
    else
       call self%forestHalos%readDatasetStatic("nodeMass"       ,nodes%nodeMass       ,firstNodeIndex,nodeCount                               )
    end if
    !$omp end critical(HDF5_Access)
    ! If only structure is requested we are done.
    if (present(structureOnly).and.structureOnly) return
    select type (nodes)
    type is (nodeDataGalacticus)
       !$omp critical(HDF5_Access)
       ! Scale or half-mass radius.
       if (present(requireScaleRadii).and.requireScaleRadii) then
          if (self%forestHalos%hasDataset("scaleRadius")) then
             nodes%halfMassRadius=-1.0d0
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("scaleRadius"   ,nodes%scaleRadius                            ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("scaleRadius"   ,nodes%scaleRadius   ,firstNodeIndex,nodeCount                               )
             end if
          else
             nodes%scaleRadius   =-1.0d0
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("halfMassRadius",nodes%halfMassRadius                         ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("halfMassRadius",nodes%halfMassRadius,firstNodeIndex,nodeCount                               )
             end if
          end if
       end if
       ! Particle count.
       if (useNodeSubset) then
          if (present(requireParticleCounts     ).and.requireParticleCounts     )                                              &
               & call self%forestHalos%readDatasetStatic("particleCount"     ,nodes%particleCount                              ,readSelection=nodeSubsetOffset)
          ! Velocity maximum.
          if (present(requireVelocityMaxima     ).and.requireVelocityMaxima     )                                              &
               & call self%forestHalos%readDatasetStatic("velocityMaximum"   ,nodes%velocityMaximum                            ,readSelection=nodeSubsetOffset)
          ! Velocity dispersion.
          if (present(requireVelocityDispersions).and.requireVelocityDispersions)                                              &
               & call self%forestHalos%readDatasetStatic("velocityDispersion",nodes%velocityDispersion                         ,readSelection=nodeSubsetOffset)
       else
          if (present(requireParticleCounts     ).and.requireParticleCounts     )                                              &
               & call self%forestHalos%readDatasetStatic("particleCount"     ,nodes%particleCount     ,firstNodeIndex,nodeCount                               )
          ! Velocity maximum.
          if (present(requireVelocityMaxima     ).and.requireVelocityMaxima     )                                              &
               & call self%forestHalos%readDatasetStatic("velocityMaximum"   ,nodes%velocityMaximum   ,firstNodeIndex,nodeCount                               )
          ! Velocity dispersion.
          if (present(requireVelocityDispersions).and.requireVelocityDispersions)                                              &
               & call self%forestHalos%readDatasetStatic("velocityDispersion",nodes%velocityDispersion,firstNodeIndex,nodeCount                               )
       end if
       ! Halo angular momenta.
       if (present(requireAngularMomenta).and.requireAngularMomenta) then
          if (self%angularMomentaIsVector) then
             if (useNodeSubset) then
                call self%forestHalos%readDataset("angularMomentum",angularMomentum3D                                                         ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDataset("angularMomentum",angularMomentum3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
             end if
             ! Transfer to nodes.
             forall(iNode=1:nodeCount(1))
                nodes(iNode)%angularMomentum=Vector_Magnitude(angularMomentum3D(:,iNode))
             end forall
          call deallocateArray(angularMomentum3D)
          else if (self%angularMomentaIsScalar) then
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("angularMomentum",nodes%angularMomentum                                   ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("angularMomentum",nodes%angularMomentum,[firstNodeIndex(1)],[nodeCount(1)]                               )
             end if
          else
             call Galacticus_Error_Report("galacticusImport","scalar angular momentum is not available")
          end if
       end if
       if (present(requireAngularMomenta3D).and.requireAngularMomenta3D) then
          if (.not.self%angularMomentaIsVector) call Galacticus_Error_Report("galacticusImport","vector angular momentum is not available")
          if (useNodeSubset) then
             call self%forestHalos%readDataset("angularMomentum",angularMomentum3D                                                         ,readSelection=nodeSubsetOffset)
          else
             call self%forestHalos%readDataset("angularMomentum",angularMomentum3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
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
             ! Transfer to nodes.
             forall(iNode=1:nodeCount(1))
                nodes(iNode)%spin=Vector_Magnitude(spin3D(:,iNode))
             end forall
          call deallocateArray(spin3D)
          else if (self%spinIsScalar) then
             if (useNodeSubset) then
                call self%forestHalos%readDatasetStatic("spin",nodes%spin                                   ,readSelection=nodeSubsetOffset)
             else
                call self%forestHalos%readDatasetStatic("spin",nodes%spin,[firstNodeIndex(1)],[nodeCount(1)]                               )
             end if
          else
             call Galacticus_Error_Report("galacticusImport","scalar spin is not available")
          end if
       end if
       if (present(requireSpin3D).and.requireSpin3D) then
          if (.not.self%spinIsVector) call Galacticus_Error_Report("galacticusImport","vector spin is not available")
          if (useNodeSubset) then
             call self%forestHalos%readDataset("spin",spin3D                                                         ,readSelection=nodeSubsetOffset)
          else
             call self%forestHalos%readDataset("spin",spin3D,[1_c_size_t,firstNodeIndex(1)],[3_c_size_t,nodeCount(1)]                               )
          end if
          ! Transfer to nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%spin3D=spin3D(:,iNode)
          end forall
       call deallocateArray(spin3D)
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
       !$omp end critical(HDF5_Access)
       ! Unit conversion.
       if     (                                                                              &
            &   present(requirePositions)                                                    &
            &  .and.                                                                         &
            &           requirePositions                                                     &
            &  .and.                                                                         &
            &   .not.                                                                        &
            &    (                                                                           &
            &      self%lengthUnit%status                                                    &
            &     .and.                                                                      &
            &      self%velocityUnit%status                                                  &
            &    )                                                                           &
            & ) call Galacticus_Error_Report(                                                &
            &                                'galacticusImport'                           ,  &
            &                                'length and velocity units must be given if '// &
            &                                'positions and velocities are to be read'       &
            &                               )
       if (self%timeUnit%status.and.self%timeUnit%scaleFactorExponent /= 0)                           &
            &   call Galacticus_Error_Report(                                                         &
            &                                'galacticusImport'                                     , &
            &                                'expect no scaling of time units with expansion factor'  &
            &                               )
       if (self%    massUnit%status.and.self%    massUnit%unitsInSI <= 0.0d0) call Galacticus_Error_Report('galacticusImport','non-positive units for mass'    )
       if (self%  lengthUnit%status.and.self%  lengthUnit%unitsInSI <= 0.0d0) call Galacticus_Error_Report('galacticusImport','non-positive units for length'  )
       if (self%velocityUnit%status.and.self%velocityUnit%unitsInSI <= 0.0d0) call Galacticus_Error_Report('galacticusImport','non-positive units for velocity')
       if (self%    timeUnit%status.and.self%    timeUnit%unitsInSI <= 0.0d0) call Galacticus_Error_Report('galacticusImport','non-positive units for time'    )
       if (.not.timesAreInternal)                                                                                                                                         &
            & nodes%nodeTime       =importerUnitConvert(nodes%nodeTime       ,nodes%nodeTime,self%timeUnit                                  ,gigaYear                 )
       nodes       %nodeMass       =importerUnitConvert(nodes%nodeMass       ,nodes%nodeTime,                                  self%massUnit,                massSolar)
       if (present(requireScaleRadii).and.requireScaleRadii) then
          nodes    %scaleRadius    =importerUnitConvert(nodes%scaleRadius    ,nodes%nodeTime,self%lengthUnit                                ,megaParsec               )
          nodes    %halfMassRadius =importerUnitConvert(nodes%halfMassRadius ,nodes%nodeTime,self%lengthUnit                                ,megaParsec               )
       end if
       if (present(requireVelocityMaxima     ).and.requireVelocityMaxima     )                                                                                                  &
            &  nodes%velocityMaximum  =importerUnitConvert(nodes%velocityMaximum   ,nodes%nodeTime,                self%velocityUnit              ,           kilo          )
       if (present(requireVelocityDispersions).and.requireVelocityDispersions)                                                                                                  &
            & nodes%velocityDispersion=importerUnitConvert(nodes%velocityDispersion,nodes%nodeTime,                self%velocityUnit              ,           kilo          )
       if (present(requireAngularMomenta     ).and.requireAngularMomenta     )                                                                                                  &
            & nodes%angularMomentum   =importerUnitConvert(nodes%angularMomentum   ,nodes%nodeTime,self%lengthUnit*self%velocityUnit*self%massUnit,megaParsec*kilo*massSolar)
       if (present(requireAngularMomenta3D).and.requireAngularMomenta3D) then
          angularmomentum3d=importerUnitConvert(angularmomentum3d,nodes%nodeTime,self%lengthUnit*self%velocityUnit*self%massUnit,megaParsec*kilo*massSolar)
          ! Transfer to nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%angularMomentum3D=angularMomentum3D(:,iNode)
          end forall
       call deallocateArray(angularMomentum3D)
       end if
       if (present(requirePositions).and.requirePositions) then
          position=importerUnitConvert(position,nodes%nodeTime,self%  lengthUnit,megaParsec)
          velocity=importerUnitConvert(velocity,nodes%nodeTime,self%velocityUnit,kilo      )
          ! Transfer to the nodes.  
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%position=position(:,iNode)
             nodes(iNode)%velocity=velocity(:,iNode)
          end forall
       call deallocateArray(position)
       call deallocateArray(velocity)
       end if
    end select
    return
  end subroutine galacticusImport
