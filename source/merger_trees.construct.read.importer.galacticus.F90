!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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
     !% Extension of the {\tt nodeData} class for \glc\ format merger trees. Stores particle indices and counts for nodes.
     integer(kind=kind_int8) :: particleIndexCount, particleIndexStart
  end type nodeDataGalacticus

  type, extends(mergerTreeImporterClass) :: mergerTreeImporterGalacticus
     !% A merger tree importer class for \glc\ format merger tree files.
     private
     type   (hdf5Object     )                            :: file                  , haloTrees
     type   (statefulInteger)                            :: hasSubhalos           , areSelfContained          , &
          &                                                 includesHubbleFlow    , periodicPositions         , &
          &                                                 lengthStatus
     type   (statefulLogical)                            :: massesAreInclusive    , angularMomentaAreInclusive
     type   (statefulDouble )                            :: length
     type   (importerUnits  )                            :: massUnit              , lengthUnit                , &
          &                                                 timeUnit              , velocityUnit
     logical                                             :: fatalMismatches       , treeIndicesRead           , &
          &                                                 angularMomentaIsScalar, angularMomentaIsVector    , &
          &                                                 spinIsScalar          , spinIsVector              , &
          &                                                 reweightTrees
     integer                                             :: treesCount
     integer                 , allocatable, dimension(:) :: firstNodes            , nodeCounts
     integer(kind=kind_int8 ), allocatable, dimension(:) :: treeIndices
     double precision        , allocatable, dimension(:) :: weights
     type   (hdf5Object     )                            :: particles
     integer                                             :: particleEpochType
     type   (varying_string )                            :: particleEpochDataSetName
   contains
     !# <workaround type="gfortran" PR="58471 58470" url="http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58471 http://gcc.gnu.org/bugzilla/show_bug.cgi?id=58470">
     !# final     :: galacticusDestructor
     !# </workaround>
     procedure :: open                          => galacticusOpen
     procedure :: close                         => galacticusClose
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
    galacticusDefaultConstructor%treeIndicesRead         =.false.
    galacticusDefaultConstructor%fatalMismatches         =mergerTreeImportGalacticusMismatchIsFatal
    galacticusDefaultConstructor%reweightTrees           =mergerTreeImportGalacticusReweightTrees
    return
  end function galacticusDefaultConstructor

  subroutine galacticusDestructor(self)
    !% Destructor for the \glc\ format merger tree importer class.
    implicit none
    type(mergerTreeImporterGalacticus), intent(inout) :: self

    call self%file%close()
    return
  end subroutine galacticusDestructor

  subroutine galacticusOpen(self,fileName)
    !% Validate a \glc\ format merger tree file.
    use Numerical_Comparison
    use Galacticus_Display
    use Galacticus_Error
    use Power_Spectra
    use Cosmology_Parameters
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout) :: self
    type            (varying_string              ), intent(in   ) :: fileName
    class           (cosmologyParametersClass    ), pointer       :: thisCosmologyParameters
    type            (hdf5Object                  )                :: cosmologicalParametersGroup, unitsGroup, angularMomentumDataset, spinDataset
    type            (varying_string              )                :: message
    character       (len=14                      )                :: valueString
    double precision                                              :: localLittleH0, localOmegaMatter, localOmegaDE, localOmegaBaryon, localSigma8, cosmologicalParameter

    ! Get the default cosmology.
    thisCosmologyParameters => cosmologyParameters()
    ! Get cosmological parameters. We do this in advance to avoid HDF5 thread conflicts.
    localLittleH0   =thisCosmologyParameters%HubbleConstant (unitsLittleH)
    localOmegaMatter=thisCosmologyParameters%OmegaMatter    (            )
    localOmegaDE    =thisCosmologyParameters%OmegaDarkEnergy(            )
    localOmegaBaryon=thisCosmologyParameters%OmegaBaryon    (            )
    localSigma8     =sigma_8                                (            )
    !$omp critical(HDF5_Access)
    ! Open the file.
    call self%file%openFile(char(fileName),readOnly=.true.)
    ! Open the merger trees group.
    self%haloTrees=self%file%openGroup("haloTrees")
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
    else if (cosmologicalParametersGroup%hasAttribute("Omega0")) then
       ! <expiry>
       !  <label>Warning regarding use of Omega0 in merger tree files.</label>
       !  <date>25-Aug-2012</date>
       ! </expiry>
       call Galacticus_Display_Message('WARNING: Use of "Omega0" in merger tree files is deprecated - use OmegaMatter instead',verbosityWarn)
       call cosmologicalParametersGroup%readAttribute("Omega0",cosmologicalParameter,allowPseudoScalar=.true.)
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
    if (self%haloTrees%hasDataset("angularMomentum")) then     
       angularMomentumDataset=self%haloTrees%openDataset("angularMomentum")
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
    if (self%haloTrees%hasDataset("spin")) then     
       spinDataset=self%haloTrees%openDataset("spin")
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
    if (self%particles%isOpen()) call self%particles%close()
    if (self%haloTrees%isOpen()) call self%haloTrees%close()
    if (self%file     %isOpen()) call self%file     %close()
    !$omp end critical(HDF5_Access)
    return
  end subroutine galacticusClose

  integer function galacticusTreesHaveSubhalos(self)
    !% Return a Boolean integer specifying whether or not the trees have subhalos.
    use Numerical_Constants_Boolean
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    if (.not.self%hasSubhalos%isSet) then
       !$omp critical(HDF5_Access)
       if (self%haloTrees%hasAttribute("treesHaveSubhalos")) then
          call self%haloTrees%readAttribute("treesHaveSubhalos",self%hasSubhalos%value,allowPseudoScalar=.true.)
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
       if (self%haloTrees%hasAttribute("haloMassesIncludeSubhalos")) then
          call self%haloTrees%readAttribute("haloMassesIncludeSubhalos",haloMassesIncludeSubhalosInteger,allowPseudoScalar=.true.)
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
       attributeExists=self%haloTrees%hasAttribute("haloAngularMomentaIncludeSubhalos")
       !$omp end critical(HDF5_Access)
       if (attributeExists) then
          !$omp critical(HDF5_Access)
          call self%haloTrees%readAttribute("haloAngularMomentaIncludeSubhalos",haloAngularMomentaIncludeSubhalosInteger,allowPseudoScalar=.true.)
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
       if (self%haloTrees%hasAttribute("treesAreSelfContained")) then
          call self%haloTrees%readAttribute("treesAreSelfContained",self%areSelfContained%value,allowPseudoScalar=.true.)
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
       if (self%haloTrees%hasAttribute("velocitiesIncludeHubbleFlow")) then
          call self%haloTrees%readAttribute("velocitiesIncludeHubbleFlow",self%includesHubbleFlow%value,allowPseudoScalar=.true.)
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
       if (self%haloTrees%hasAttribute("positionsArePeriodic")) then
          call self%haloTrees%readAttribute("positionsArePeriodic",self%periodicPositions%value,allowPseudoScalar=.true.)
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
    double precision                                                        :: lengthSimulationBox

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
    if (self%lengthStatus%value == booleanTrue) galacticusCubeLength=importerUnitConvert(self%length%value,time,self%lengthUnit,megaParsec)
    if (present(status)) then
       status=self%lengthStatus%value
    else
    if (self%lengthStatus%value == booleanFalse) call Galacticus_Error_Report('galacticusCubeLength','the boxSize attribute of the simulation group is required')
    end if
    return
  end function galacticusCubeLength

  integer function galacticusTreeCount(self)
    !% Return a count of the number of trees available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    call galacticusTreeIndicesRead(self)
    galacticusTreeCount=self%treesCount
    return
  end function galacticusTreeCount

  integer(kind=kind_int8) function galacticusTreeIndex(self,i)
    !% Return the index of the $i^{\rm th}$ tree.
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call galacticusTreeIndicesRead(self)
    galacticusTreeIndex=self%treeIndices(i)
    return
  end function galacticusTreeIndex

  integer function galacticusNodeCount(self,i)
    !% Return a count of the number of nodes in the $i^{\rm th}$ tree.
    implicit none
    class  (mergerTreeImporterGalacticus), intent(inout) :: self
    integer                              , intent(in   ) :: i

    call galacticusTreeIndicesRead(self)
    galacticusNodeCount=self%nodeCounts(i)
    return
  end function galacticusNodeCount

  subroutine galacticusTreeIndicesRead(self)
    !% Read the tree indices.
    use Galacticus_Error
    use HDF5
    use Sort
    use Cosmology_Functions
    use Numerical_Constants_Astronomical
    use Halo_Mass_Function
    implicit none
    class           (mergerTreeImporterGalacticus), intent(inout)             :: self
    type            (hdf5Object                  )                            :: treeIndexGroup
    integer         (kind=kind_int8              ), allocatable, dimension(:) :: descendentIndex
    double precision                              , allocatable, dimension(:) :: nodeMass       , treeMass , nodeTime, treeTime
    integer         (kind=HSIZE_T                )             , dimension(1) :: firstNodeIndex , nodeCount
    integer         (kind=c_size_t               ), allocatable, dimension(:) :: sortOrder
    class           (cosmologyFunctionsClass     ), pointer                   :: cosmologyFunctionsDefault
    integer                                                                   :: i              , iNode
    double precision                                                          :: massMinimum    , massMaximum
    
    if (self%treeIndicesRead) return
    if (self%file%hasGroup("treeIndex")) then
       treeIndexGroup=self%file%openGroup("treeIndex")
       call treeIndexGroup%readDataset("firstNode"    ,self%firstNodes )
       call treeIndexGroup%readDataset("numberOfNodes",self%nodeCounts )
       call treeIndexGroup%readDataset("treeIndex"    ,self%treeIndices)
       if (self%reweightTrees) then
          cosmologyFunctionsDefault => cosmologyFunctions()
          allocate(self%weights(size(self%firstNodes)))
          allocate(treeMass    (size(self%firstNodes)))
          allocate(treeTime    (size(self%firstNodes)))
          do i=1,size(self%firstNodes)
             firstNodeIndex(1)=self%firstNodes(i)+1
             nodeCount     (1)=self%nodeCounts(i)
             ! Allocate the nodes array.
             allocate(descendentIndex(nodeCount(1)))
             allocate(nodeMass       (nodeCount(1)))
             allocate(nodeTime       (nodeCount(1)))
             !$omp critical(HDF5_Access)
             call self%haloTrees%readDatasetStatic("descendentIndex",descendentIndex,firstNodeIndex,nodeCount)
             call self%haloTrees%readDatasetStatic("nodeMass"       ,nodeMass       ,firstNodeIndex,nodeCount)
             if      (self%haloTrees%hasDataset("time"           )) then
                ! Time is present, so read it.
                call self%haloTrees%readDatasetStatic("time"           ,nodeTime,firstNodeIndex,nodeCount)
                nodeTime=importerUnitConvert(nodeTime,nodeTime,self%timeUnit,gigaYear)
             else if (self%haloTrees%hasDataset("expansionFactor")) then
                ! Expansion factor is present, read it instead.
                call self%haloTrees%readDatasetStatic("expansionFactor",nodeTime,firstNodeIndex,nodeCount)
                ! Convert expansion factors to times.
                do iNode=1,nodeCount(1)
                   nodeTime(iNode)=cosmologyFunctionsDefault%cosmicTime(nodeTime(iNode))
                end do
             else if (self%haloTrees%hasDataset("redshift"       )) then
                ! Redshift is present, read it instead.
                call self%haloTrees%readDatasetStatic("redshift"       ,nodeTime,firstNodeIndex,nodeCount)
                ! Convert redshifts to times.
                do iNode=1,nodeCount(1)
                   nodeTime(iNode)=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(nodeTime(iNode)))
                end do
             else
                call Galacticus_Error_Report("galacticusImport","one of time, redshift or expansionFactor data sets must be present in haloTrees group")
             end if
             !$omp end critical(HDF5_Access)
             if (count(descendentIndex == -1) /= 1) call Galacticus_Error_Report('galacticusTreeIndicesRead','reweighting trees requires there to be only only root node')
             treeMass(i)=sum(nodeMass,mask=descendentIndex == -1)
             treeTime(i)=sum(nodeTime,mask=descendentIndex == -1)
             deallocate(descendentIndex)
             deallocate(nodeMass       )
             deallocate(nodeTime       )
          end do
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
             self%weights(sortOrder(i))=Halo_Mass_Function_Integrated(treeTime(sortOrder(i)),massMinimum,massMaximum)
          end do
          call treeIndexGroup%readDatasetStatic("treeWeight",self%weights)
          deallocate(treeMass)
          deallocate(treeTime)
       else if (treeIndexGroup%hasDataset("treeWeight")) then
          call treeIndexGroup%readDataset("treeWeight",self%weights)
       end if
       call treeIndexGroup%close()
       self%treesCount=size(self%treeIndices)
       ! Reset first node indices to Fortran array standard.
       self%firstNodes=self%firstNodes+1
    else
       call Galacticus_Error_Report('galacticusTreeIndicesRead','merger tree file must contain the treeIndex group')
    end if
    self%treeIndicesRead=.true.
    return
  end subroutine galacticusTreeIndicesRead

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
    type            (hdf5Object                  )                :: simulationGroup
    double precision                                              :: lengthSimulationBox      , timePresent
    integer                                                       :: statusActual

    call galacticusTreeIndicesRead(self)
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
    if (positions .and..not.self%haloTrees%hasDataset("position")) galacticusPositionsAvailable=.false.
    if (velocities.and..not.self%haloTrees%hasDataset("velocity")) galacticusPositionsAvailable=.false.
    !$omp end critical(HDF5_Access)
    return
  end function galacticusPositionsAvailable

  logical function galacticusScaleRadiiAvailable(self)
    !% Return true if scale radii are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    galacticusScaleRadiiAvailable=                      &
         &  self%haloTrees%hasDataset("halfMassRadius") &
         & .or.                                         &
         &  self%haloTrees%hasDataset("position"      )
    !$omp end critical(HDF5_Access)
    return
  end function galacticusScaleRadiiAvailable

  logical function galacticusParticleCountAvailable(self)
    !% Return true if particle counts are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    
    !$omp critical(HDF5_Access)
    galacticusParticleCountAvailable=self%haloTrees%hasDataset("particleCount")
    !$omp end critical(HDF5_Access)
    return
  end function galacticusParticleCountAvailable

  logical function galacticusVelocityMaximumAvailable(self)
    !% Return true if halo rotation curve velocity maxima are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    galacticusVelocityMaximumAvailable=self%haloTrees%hasDataset("velocityMaximum")
    !$omp end critical(HDF5_Access)
    return
  end function galacticusVelocityMaximumAvailable

  logical function galacticusVelocityDispersionAvailable(self)
    !% Return true if halo velocity dispersions are available.
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self

    !$omp critical(HDF5_Access)
    galacticusVelocityDispersionAvailable=self%haloTrees%hasDataset("velocityDispersion")
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

  integer function galacticusSubhaloTraceCount(self,node)
    !% Returns the length of a subhalo trace.
    use Galacticus_Error
    implicit none
    class(mergerTreeImporterGalacticus), intent(inout) :: self
    class(nodeData                    ), intent(in   ) :: node
    
    select type (node)
    type is (nodeDataGalacticus)
       galacticusSubhaloTraceCount=node%particleIndexCount
    class default
       call Galacticus_Error_Report('galacticusSubhaloTraceCount','node should be of type nodeDataGalacticus')
    end select
    return
  end function galacticusSubhaloTraceCount

  subroutine galacticusImport(self,i,nodes,requireScaleRadii,requireAngularMomenta,requireAngularMomenta3D,requireSpin,requireSpin3D,requirePositions,requireParticleCounts,requireVelocityMaxima,requireVelocityDispersions)
    !% Import the $i^{\rm th}$ merger tree.
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
    class           (nodeData                    ), intent(  out), allocatable, dimension(:  ) :: nodes
    logical                                       , intent(in   ), optional                    :: requireScaleRadii        , requireAngularMomenta, &
         &                                                                                        requireAngularMomenta3D  , requirePositions     , &
         &                                                                                        requireParticleCounts    , requireVelocityMaxima, &
         &                                                                                        requireVelocityDispersions, requireSpin         , &
         &                                                                                        requireSpin3D
    class           (cosmologyFunctionsClass     ), pointer                                    :: cosmologyFunctionsDefault
    integer         (kind=HSIZE_T                )                            , dimension(1  ) :: firstNodeIndex           , nodeCount
    integer         (kind=kind_int8              )                                             :: iNode
    double precision                                             , allocatable, dimension(:,:) :: angularMomentum3D        , position             , &
         &                                                                                        velocity                 , spin3D
    logical                                                                                    :: timesAreInternal

    ! Get the default cosmology functions object.
    cosmologyFunctionsDefault => cosmologyFunctions()
    ! Ensure tree indices have been read.
    call galacticusTreeIndicesRead(self)
    ! Determine the first node index and the node count.
    firstNodeIndex(1)=self%firstNodes(i)
    nodeCount     (1)=self%nodeCounts(i)
    ! Allocate the nodes array.
    allocate(nodeDataGalacticus :: nodes(nodeCount(1)))
    call Memory_Usage_Record(sizeof(nodes))
    !$omp critical(HDF5_Access)
    ! nodeIndex
    call self%haloTrees%readDatasetStatic("nodeIndex"      ,nodes%nodeIndex      ,firstNodeIndex,nodeCount)
    ! hostIndex
    call self%haloTrees%readDatasetStatic("hostIndex"      ,nodes%hostIndex      ,firstNodeIndex,nodeCount)
    ! parentNode
    call self%haloTrees%readDatasetStatic("descendentIndex",nodes%descendentIndex,firstNodeIndex,nodeCount)
    ! nodeMass
    call self%haloTrees%readDatasetStatic("nodeMass"       ,nodes%nodeMass       ,firstNodeIndex,nodeCount)
    ! nodeTime
    timesAreInternal=.true. 
    if      (self%haloTrees%hasDataset("time"           )) then
       ! Time is present, so read it.
       timesAreInternal=.false.
       call self%haloTrees%readDatasetStatic("time"           ,nodes%nodeTime,firstNodeIndex,nodeCount)
    else if (self%haloTrees%hasDataset("expansionFactor")) then
       ! Expansion factor is present, read it instead.
       call self%haloTrees%readDatasetStatic("expansionFactor",nodes%nodeTime,firstNodeIndex,nodeCount)
       ! Validate expansion factors.
       if (any(nodes%nodeTime <= 0.0d0)) call Galacticus_Error_Report("galacticusImport","expansionFactor dataset values must be >0")
       if (any(nodes%nodeTime >  1.0d0)) call Galacticus_Display_Message("WARNING: some expansion factors are in the future when importing merger tree",verbosityWarn)
       ! Convert expansion factors to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=cosmologyFunctionsDefault%cosmicTime(nodes(iNode)%nodeTime)
       end do
    else if (self%haloTrees%hasDataset("redshift"       )) then
       ! Redshift is present, read it instead.
       call self%haloTrees%readDatasetStatic("redshift"       ,nodes%nodeTime,firstNodeIndex,nodeCount)
      ! Validate redshifts.
       if (any(nodes%nodeTime <= -1.0d0)) call Galacticus_Error_Report("galacticusImport","redshift dataset values must be >-1")
       if (any(nodes%nodeTime <   0.0d0)) call Galacticus_Display_Message("WARNING: some redshifts are in the future when importing merger tree",verbosityWarn)
       ! Convert redshifts to times.
       do iNode=1,nodeCount(1)
          nodes(iNode)%nodeTime=cosmologyFunctionsDefault%cosmicTime(cosmologyFunctionsDefault%expansionFactorFromRedshift(nodes(iNode)%nodeTime))
       end do
    else
       call Galacticus_Error_Report("galacticusImport","one of time, redshift or expansionFactor data sets must be present in haloTrees group")
    end if
    ! Scale or half-mass radius.
    if (present(requireScaleRadii).and.requireScaleRadii) then
       if (self%haloTrees%hasDataset("scaleRadius")) then
          nodes%halfMassRadius=-1.0d0
          call self%haloTrees%readDatasetStatic("scaleRadius"   ,nodes%scaleRadius   ,firstNodeIndex,nodeCount)
       else
          nodes%scaleRadius   =-1.0d0
          call self%haloTrees%readDatasetStatic("halfMassRadius",nodes%halfMassRadius,firstNodeIndex,nodeCount)
       end if
    end if
    ! Particle count.
    if (present(requireParticleCounts     ).and.requireParticleCounts     )                                              &
         & call self%haloTrees%readDatasetStatic("particleCount"     ,nodes%particleCount     ,firstNodeIndex,nodeCount)
    ! Velocity maximum.
    if (present(requireVelocityMaxima     ).and.requireVelocityMaxima     )                                              &
         & call self%haloTrees%readDatasetStatic("velocityMaximum"   ,nodes%velocityMaximum   ,firstNodeIndex,nodeCount)
    ! Velocity dispersion.
    if (present(requireVelocityDispersions).and.requireVelocityDispersions)                                              &
         & call self%haloTrees%readDatasetStatic("velocityDispersion",nodes%velocityDispersion,firstNodeIndex,nodeCount)
    ! Halo angular momenta.
    if (present(requireAngularMomenta).and.requireAngularMomenta) then
       if (self%angularMomentaIsVector) then
          call self%haloTrees%readDataset("angularMomentum",angularMomentum3D,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
          ! Transfer to nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%angularMomentum=Vector_Magnitude(angularMomentum3D(:,iNode))
          end forall
          call Dealloc_Array(angularMomentum3D)
       else if (self%angularMomentaIsScalar) then
          call self%haloTrees%readDatasetStatic("angularMomentum",nodes%angularMomentum,[firstNodeIndex(1)],[nodeCount(1)])
       else
          call Galacticus_Error_Report("galacticusImport","scalar angular momentum is not available")
       end if
    end if
    if (present(requireAngularMomenta3D).and.requireAngularMomenta3D) then
       if (.not.self%angularMomentaIsVector) call Galacticus_Error_Report("galacticusImport","vector angular momentum is not available")
       call self%haloTrees%readDataset("angularMomentum",angularMomentum3D,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
    end if
    ! Halo spins.
    if (present(requireSpin).and.requireSpin) then
       if (self%spinIsVector) then
          call self%haloTrees%readDataset("spin",spin3D,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
          ! Transfer to nodes.
          forall(iNode=1:nodeCount(1))
             nodes(iNode)%spin=Vector_Magnitude(spin3D(:,iNode))
          end forall
          call Dealloc_Array(spin3D)
       else if (self%spinIsScalar) then
          call self%haloTrees%readDatasetStatic("spin",nodes%spin,[firstNodeIndex(1)],[nodeCount(1)])
       else
          call Galacticus_Error_Report("galacticusImport","scalar spin is not available")
       end if
    end if
    if (present(requireSpin3D).and.requireSpin3D) then
       if (.not.self%spinIsVector) call Galacticus_Error_Report("galacticusImport","vector spin is not available")
       call self%haloTrees%readDataset("spin",spin3D,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
       ! Transfer to nodes.
       forall(iNode=1:nodeCount(1))
          nodes(iNode)%spin3D=spin3D(:,iNode)
       end forall
       call Dealloc_Array(spin3D)
    end if
    select type (nodes)
    type is (nodeDataGalacticus)
       ! Initialize particle data to null values.
       nodes%particleIndexStart=-1_kind_int8
       nodes%particleIndexCount=-1_kind_int8
       ! Positions (and velocities).
       if (present(requirePositions).and.requirePositions) then
          ! position.
          call self%haloTrees%readDataset("position",position,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
          ! velocity.
          call self%haloTrees%readDataset("velocity",velocity,[int(1,kind=kind_int8),firstNodeIndex(1)],[int(3,kind=kind_int8),nodeCount(1)])
          ! If a set of most bound particle indices are present, read them.
          if (self%haloTrees%hasDataset("particleIndexStart").and.self%haloTrees%hasDataset("particleIndexCount")) then
             call self%haloTrees%readDatasetStatic("particleIndexStart",nodes%particleIndexStart,firstNodeIndex,nodeCount)
             call self%haloTrees%readDatasetStatic("particleIndexCount",nodes%particleIndexCount,firstNodeIndex,nodeCount)
          end if
       end if
    class default
       call Galacticus_Error_Report('galacticusImport','nodes should be of type nodeDataGalacticus')
    end select
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
       call Dealloc_Array(angularMomentum3D)
    end if
    if (present(requirePositions).and.requirePositions) then
       position=importerUnitConvert(position,nodes%nodeTime,self%  lengthUnit,megaParsec)
       velocity=importerUnitConvert(velocity,nodes%nodeTime,self%velocityUnit,kilo      )
       ! Transfer to the nodes.  
       forall(iNode=1:self%nodeCounts(i))
          nodes(iNode)%position=position(:,iNode)
          nodes(iNode)%velocity=velocity(:,iNode)
       end forall
       call Dealloc_Array(position)
       call Dealloc_Array(velocity)
    end if
    return
  end subroutine galacticusImport
