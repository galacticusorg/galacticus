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

!% Contains a module which implements a merger tree operator which creates particle representations of \glc\ halos.

  use Cosmology_Parameters
  use Cosmology_Functions
  use ISO_Varying_String
  use Kind_Numbers

  !# <mergerTreeOperator name="mergerTreeOperatorParticulate" defaultThreadPrivate="yes">
  !#  <description>Provides a merger tree operator which create particle representations of \glc\ halos.</description>
  !# </mergerTreeOperator>
  type, extends(mergerTreeOperatorClass) :: mergerTreeOperatorParticulate
     !% A merger tree operator which create particle representations of \glc\ halos.
     private
     class           (cosmologyParametersClass), pointer :: cosmologyParameters_
     class           (cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     type            (varying_string          )          :: outputFileName
     double precision                                    :: massParticle        , radiusTruncateOverRadiusVirial, &
          &                                                 timeSnapshot
     logical                                             :: satelliteOffset     , nonCosmological               , &
          &                                                 positionOffset
     integer         (kind_int8               )          :: idMultiplier
   contains
     final     ::            particulateDestructor
     procedure :: operate => particulateOperate
  end type mergerTreeOperatorParticulate

  interface mergerTreeOperatorParticulate
     !% Constructors for the particulate merger tree operator class.
     module procedure particulateConstructorParameters
     module procedure particulateConstructorInternal
  end interface mergerTreeOperatorParticulate

contains

  function particulateConstructorParameters(parameters)
    !% Constructor for the particulate merger tree operator class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type(mergerTreeOperatorParticulate)                :: particulateConstructorParameters
    type(inputParameters              ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />
    
    call parameters%checkParameters(allowedParameterNames)
    !# <inputParameter>
    !#   <name>outputFileName</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%outputFileName</variable>
    !#   <description>Name of the file to which particle data should be written.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>idMultiplier</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0_kind_int8</defaultValue>
    !#   <variable>particulateConstructorParameters%idMultiplier</variable>
    !#   <description>If this parameter is greater than zero, particle IDs begin at {\normalfont \ttfamily nodeIndex\*[idMultiplier]} for each node. The multiplier should be chosen to be large enough that duplicate IDs can not occur.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%massParticle</variable>
    !#   <description>Mass of particles to be used to represent halos.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>timeSnapshot</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%timeSnapshot</variable>
    !#   <description>The time at which to snapshot the tree.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>radiusTruncateOverRadiusVirial</name>
    !#   <source>parameters</source>
    !#   <variable>particulateConstructorParameters%radiusTruncateOverRadiusVirial</variable>
    !#   <description>Radius (in units of the virial radius) at which to truncate halo profiles.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>satelliteOffset</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.true.</defaultValue>
    !#   <variable>particulateConstructorParameters%satelliteOffset</variable>
    !#   <description>If true, offset particle representations to the positions/velocities of satellites.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>positionOffset</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%positionOffset</variable>
    !#   <description>If true, offset particle representations to the positions/velocities of nodes.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>nonCosmological</name>
    !#   <source>parameters</source>
    !#   <defaultValue>.false.</defaultValue>
    !#   <variable>particulateConstructorParameters%nonCosmological</variable>
    !#   <description>If true, a non-cosmological snapshot file will be created.</description>
    !#   <type>boolean</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="cosmologyParameters" name="particulateConstructorParameters%cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="particulateConstructorParameters%cosmologyFunctions_"  source="parameters"/>
    return
  end function particulateConstructorParameters

  function particulateConstructorInternal(outputFileName,idMultiplier,massParticle,radiusTruncateOverRadiusVirial,timeSnapshot,satelliteOffset,positionOffset,nonCosmological,cosmologyParameters_,cosmologyFunctions_)
    !% Internal constructor for the particulate merger tree operator class.
    implicit none
    type            (mergerTreeOperatorParticulate)                        :: particulateConstructorInternal
    type            (varying_string               ), intent(in   )         :: outputFileName
    integer         (kind_int8                    ), intent(in   )         :: idMultiplier
    double precision                               , intent(in   )         :: massParticle                  , radiusTruncateOverRadiusVirial, &
         &                                                                    timeSnapshot
    logical                                        , intent(in   )         :: satelliteOffset               , nonCosmological               , &
         &                                                                    positionOffset
    class           (cosmologyParametersClass     ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass      ), intent(in   ), target :: cosmologyFunctions_

    particulateConstructorInternal%outputFileName                 =  outputFileName
    particulateConstructorInternal%idMultiplier                   =  idMultiplier
    particulateConstructorInternal%massParticle                   =  massParticle
    particulateConstructorInternal%radiusTruncateOverRadiusVirial =  radiusTruncateOverRadiusVirial
    particulateConstructorInternal%timeSnapshot                   =  timeSnapshot
    particulateConstructorInternal%satelliteOffset                =  satelliteOffset
    particulateConstructorInternal%positionOffset                 =  positionOffset
    particulateConstructorInternal%nonCosmological                =  nonCosmological
    particulateConstructorInternal%cosmologyParameters_           => cosmologyParameters_  
    particulateConstructorInternal%cosmologyFunctions_            => cosmologyFunctions_
   return
  end function particulateConstructorInternal

  subroutine particulateDestructor(self)
    !% Destructor for the particulate merger tree operator class.
    implicit none
    type(mergerTreeOperatorParticulate), intent(inout) :: self

    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_" />
    return
 end subroutine particulateDestructor

  subroutine particulateOperate(self,tree)
    !% Perform a particulation operation on a merger tree (i.e. create a particle representation of the tree).
    use IO_HDF5
    use FGSL
    use Numerical_Constants_Math
    use Pseudo_Random
    use Poisson_Random
    use Gaussian_Random
    use Dark_Matter_Halo_Scales
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Velocity_Dispersions
    use Galactic_Structure_Options
    use Galacticus_Calculations_Resets
    use Coordinates
    use Numerical_Comparison
    use Memory_Management
    implicit none
    class           (mergerTreeOperatorParticulate), intent(inout)               :: self
    type            (mergerTree                   ), intent(inout) , target      :: tree
    type            (treeNode                     ), pointer                     :: node
    class           (nodeComponentBasic           ), pointer                     :: basic
    class           (nodeComponentPosition        ), pointer                     :: position
    class           (nodeComponentSatellite       ), pointer                     :: satellite
    type            (mergerTree                   ), pointer                     :: treeCurrent
    class           (darkMatterHaloScaleClass     ), pointer                     :: darkMatterHaloScale_
    double precision                               , parameter                   :: tolerance           =1.0d-06
    double precision                               , parameter                   :: unitGadgetMass      =1.0d+10
    double precision                               , parameter                   :: unitGadgetLength    =1.0d-03
    double precision                               , parameter                   :: unitGadgetVelocity  =1.0d+00
    double precision                               , dimension(3  )              :: positionVector             , velocityVector
    integer                                        , dimension(6  )              :: particleCounts
    double precision                               , dimension(:,:), allocatable :: particlePosition           , particleVelocity
    integer         (kind_int8                    ), dimension(  :), allocatable :: particleIDs
    double precision                                                             :: particleCountMean          , velocityDispersion, &
         &                                                                          radiusVirial               , radiusTruncate    , &
         &                                                                          massFractionTruncate       , massTruncate
    integer                                                                      :: particleCountActual        , i
    type            (fgsl_rng                     )                              :: pseudoSequenceObject
    logical                                                                      :: pseudoSequenceReset =.true.
    type            (coordinateCartesian          )                              :: positionCartesian          , velocityCartesian
    type            (coordinateSpherical          )                              :: positionSpherical
    type            (hdf5Object                   )                              :: outputFile                 , header            , &
         &                                                                          particleGroup
    
    ! Get required objects.
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Open the HDF5 file for output.
    !$omp critical (HDF5_Access)
    call outputFile%openFile(char(self%outputFileName),overWrite=.false.,readOnly=.false.)
    ! Check if the Gadget header has already been created.
    if (.not.outputFile%hasGroup('Header')) then
       header=outputFile%openGroup('Header','Group containing Gadget metadata.')
       ! Particle properties.
       call header%writeAttribute(       1                                    ,'NumFilesPerSnapshot'   )
       call header%writeAttribute(spread(0                               ,1,6),'NumPart_ThisFile'      )
       call header%writeAttribute(spread(0                               ,1,6),'NumPart_Total_HighWord')
       call header%writeAttribute(spread(0                               ,1,6),'NumPart_Total'         )
       call header%writeAttribute(spread(self%massParticle/unitGadgetMass,1,6),'MassTable'             )
       ! Time.
       call header%writeAttribute(                                                                                              self%timeSnapshot  ,'Time'    )
       call header%writeAttribute(self%cosmologyFunctions_%redshiftFromExpansionFactor(self%cosmologyFunctions_%expansionFactor(self%timeSnapshot)),'Redshift')
       ! Cosmology.
       if (self%nonCosmological) then
          call header%writeAttribute(1.0d0                                                        ,'HubbleParam')
          call header%writeAttribute(0.0d0                                                        ,'Omega0'     )
          call header%writeAttribute(0.0d0                                                        ,'OmegaLambda')
       else
          call header%writeAttribute(self%cosmologyParameters_%HubbleConstant (hubbleUnitsLittleH),'HubbleParam')
          call header%writeAttribute(self%cosmologyParameters_%OmegaMatter    (                  ),'Omega0'     )
          call header%writeAttribute(self%cosmologyParameters_%OmegaDarkEnergy(                  ),'OmegaLambda')
       end if
       call header%writeAttribute(0.0d0,'BoxSize'         )
       ! Flags.
       call header%writeAttribute(0    ,'Flag_Cooling'    )
       call header%writeAttribute(0    ,'Flag_Sfr'        )
       call header%writeAttribute(0    ,'Flag_Feedback'   )
       call header%writeAttribute(0    ,'Flag_StellarAge' )
       call header%writeAttribute(0    ,'FlagMetals'      )
       call header%writeAttribute(0    ,'Flag_Entropy_ICs')
       call header%close()
    end if
    call outputFile%close()
    !$omp end critical (HDF5_Access)
    ! Iterate over trees.
    treeCurrent             => tree
    do while (associated(treeCurrent))
       ! Get root node of the tree.
       node  => treeCurrent%baseNode
       ! Walk the tree.
       do while (associated(node))
          call Galacticus_Calculations_Reset(node)
          ! Check if the node exists at the snapshot time.
          basic => node%basic()
          if (Values_Agree(basic%time(),self%timeSnapshot,relTol=tolerance)) then
             ! Determine the virial radius.
             radiusVirial=darkMatterHaloScale_%virialRadius(node)
             ! Determine the truncation radius.
             radiusTruncate=+self%radiusTruncateOverRadiusVirial &
                  &         *radiusVirial
             ! Determine the mass within the truncation radius in units of the virial radius.
             massTruncate=Galactic_Structure_Enclosed_Mass(                             &
                  &                                        node                       , &
                  &                                        radiusTruncate             , &
                  &                                        massType      =massTypeDark, &
                  &                                        haloLoaded    =.true.        &
                  &                                       ) 
             massFractionTruncate=+massTruncate      &
                  &               /basic%mass()
             ! Determine the mean number of particle required to represent this node.          
             particleCountMean  =+massTruncate      &
                  &              /self%massParticle
             ! Determine the actual number of particles to use to represent the node.
             particleCountActual=Poisson_Random_Get(pseudoSequenceObject,particleCountMean,pseudoSequenceReset)
             ! Allocate space for particle data.
             call Alloc_Array(particlePosition,[3,particleCountActual])
             call Alloc_Array(particleVelocity,[3,particleCountActual])
             call Alloc_Array(particleIDs     ,[  particleCountActual])
             ! Get required components.
             position  => node%position ()
             satellite => node%satellite()
             ! Iterate over particles.
             do i=1,particleCountActual
                ! Sample particle positions from the halo density distribution. Currently, we assume that halos are spherically
                ! symmetric.
                call positionSpherical%  phiSet(     2.0d0*Pi*Pseudo_Random_Get(pseudoSequenceObject,pseudoSequenceReset)       )
                call positionSpherical%thetaSet(acos(2.0d0   *Pseudo_Random_Get(pseudoSequenceObject,pseudoSequenceReset)-1.0d0))
                call positionSpherical%    rSet(                                                                                                                      &
                     &                          Galactic_Structure_Radius_Enclosing_Mass(                                                                             &
                     &                                                                   node                                                                       , &
                     &                                                                   fractionalMass=+massFractionTruncate                                         &
                     &                                                                                  *Pseudo_Random_Get(pseudoSequenceObject,pseudoSequenceReset), &
                     &                                                                   massType      =massTypeDark                                                , &
                     &                                                                   haloLoaded    =.true.                                                        &
                     &                                                                  )                                                                             &
                     &                         )
                ! Get the corresponding cartesian coordinates.
                positionCartesian=positionSpherical
               ! Find the velocity dispersion. Currently we assume isotropic velocity dispersion.
                velocityDispersion=Galactic_Structure_Velocity_Dispersion(                                             &
                     &                                                    node                                       , &
                     &                                                    positionSpherical%r()                      , &
                     &                                                    radiusOuter          =radiusTruncate       , &
                     &                                                    componentType        =componentTypeDarkHalo, &
                     &                                                    massType             =massTypeDark         , &
                     &                                                    haloLoaded           =.true.                 &
                     &                                                   )
                ! Select velocity components assuming isotropic velocity dispersion.
                call velocityCartesian%xSet(Gaussian_Random_Get(pseudoSequenceObject,velocityDispersion,pseudoSequenceReset))
                call velocityCartesian%ySet(Gaussian_Random_Get(pseudoSequenceObject,velocityDispersion,pseudoSequenceReset))
                call velocityCartesian%zSet(Gaussian_Random_Get(pseudoSequenceObject,velocityDispersion,pseudoSequenceReset))
                ! Offset position and velocity to the position and velocity of the node.
                if (self%positionOffset) then
                   positionVector=position%position()
                   velocityVector=position%velocity()
                   call positionCartesian%xSet(positionCartesian%x()+positionVector(1))
                   call positionCartesian%ySet(positionCartesian%y()+positionVector(2))
                   call positionCartesian%zSet(positionCartesian%z()+positionVector(3))
                   call velocityCartesian%xSet(velocityCartesian%x()+velocityVector(1))
                   call velocityCartesian%ySet(velocityCartesian%y()+velocityVector(2))
                   call velocityCartesian%zSet(velocityCartesian%z()+velocityVector(3))
                end if
                ! Offset position and velocity to position and velocity of satellite.
                if (self%satelliteOffset) then
                   positionVector=satellite%position()
                   velocityVector=satellite%velocity()
                   call positionCartesian%xSet(positionCartesian%x()+positionVector(1))
                   call positionCartesian%ySet(positionCartesian%y()+positionVector(2))
                   call positionCartesian%zSet(positionCartesian%z()+positionVector(3))
                   call velocityCartesian%xSet(velocityCartesian%x()+velocityVector(1))
                   call velocityCartesian%ySet(velocityCartesian%y()+velocityVector(2))
                   call velocityCartesian%zSet(velocityCartesian%z()+velocityVector(3))
                end if
                ! Store to particle data.
                particlePosition(:,i)=[positionCartesian%x(),positionCartesian%y(),positionCartesian%z()]
                particleVelocity(:,i)=[velocityCartesian%x(),velocityCartesian%y(),velocityCartesian%z()]
                particleIDs     (  i)=i-1
             end do
             ! Perform unit conversion.
             particlePosition=particlePosition/unitGadgetLength
             particleVelocity=particleVelocity/unitGadgetVelocity
             ! Accumulate the particle data to file.
             !$omp critical (HDF5_Access)
             call outputFile%openFile(char(self%outputFileName),overWrite=.false.,readOnly=.false.,objectsOverwritable=.true.)
             ! Get current count of particles in file.
             header=outputFile%openGroup('Header','Group containing Gadget metadata.')
             call header%readAttributeStatic('NumPart_Total',particleCounts)
             ! Offset particle IDs.
             if (self%idMultiplier > 0) then
                particleIDs=particleIDs+node%index()*self%idMultiplier
             else
                particleIDs=particleIDs+particleCounts(2)
             end if
             ! Write particle data.
             particleGroup=outputFile%openGroup('PartType1','Group containing particle data for halos')
             call particleGroup%writeDataset(particlePosition,'Coordinates','Particle coordinates',appendTo=.true.,appendDimension=2)
             call particleGroup%writeDataset(particleVelocity,'Velocities' ,'Particle velocities' ,appendTo=.true.,appendDimension=2)
             call particleGroup%writeDataset(particleIDs     ,'ParticleIDs','Particle IDs'        ,appendTo=.true.                  )
             call particleGroup%close()
             call Dealloc_Array(particlePosition)
             call Dealloc_Array(particleVelocity)
             call Dealloc_Array(particleIDs     )
             ! Update particle counts.
             particleCounts(2)=particleCounts(2)+particleCountActual
             call header%writeAttribute(particleCounts,'NumPart_ThisFile')
             call header%writeAttribute(particleCounts,'NumPart_Total'   )
             call header%close()
             call outputFile%close()
             !$omp end critical (HDF5_Access)
          end if
          ! Walk to the next node.
          node => node%walkTreeWithSatellites()
       end do
       ! Move to the next tree.
       treeCurrent => treeCurrent%nextTree
    end do
    return
  end subroutine particulateOperate
