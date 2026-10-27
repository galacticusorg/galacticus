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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a class of dark matter halo filaments placed at random, subject to a minimum separation.
  !!}

  !![
  <darkMatterHaloFilaments name="darkMatterHaloFilamentsRandomSequential" docformat="rst">
   <description>
   A dark matter halo filaments class in which each branch of the merger tree is fed by a number of filaments placed at random
   about the halo, subject to a minimum angular separation.

   The number of filaments, :math:`N`, has a mean of :math:`\langle N \rangle` (given by ``countFilamentsMean``) and the
   smallest variance possible for that mean: :math:`N` is :math:`\lfloor \langle N \rangle \rfloor + 1` with probability
   :math:`\langle N \rangle - \lfloor \langle N \rangle \rfloor`, and :math:`\lfloor \langle N \rangle \rfloor` otherwise. The
   first two filaments are a single filament passing through the halo, i.e. two opposite directions, oriented isotropically at
   random. Each further filament is placed at a direction drawn isotropically at random, accepted only if it is separated by
   at least ``angleSeparationMinimum`` from every filament already placed (random sequential placement). Since every further
   filament must clear both directions of the first filament, the separation must be less than :math:`90^\circ` if
   :math:`\langle N \rangle > 2`. (At exactly :math:`90^\circ` a further filament would have to be exactly perpendicular to
   the first, which a random direction never is.)

   The default, :math:`\langle N \rangle=3` (a filament passing through the halo plus one further filament) with a minimum
   separation of :math:`60^\circ`, is the configuration for which :galacticus-class:`virialOrbitFilamentary` is calibrated:
   with it, that class reproduces the distributions of the angles between the points at which subhalos enter their hosts and
   both the axes of slowest and fastest collapse of the velocity shear tensor measured by :cite:t:`libeskind_universal_2014`,
   if those axes are identified with the major and minor axes of the tensor
   :math:`\sum_i \hat{\mathbf{f}}_i \hat{\mathbf{f}}_i^\mathrm{T}` of the filament directions, :math:`\hat{\mathbf{f}}_i`.
   :math:`\langle N \rangle=3` agrees with the typical number of filaments connected to a halo, at all masses, in the peak
   theory calculations of :cite:t:`codis_connectivity_2018`, but exceeds the number crossing a sphere of
   :math:`1.5 R_\mathrm{200}` about groups in simulations, which is typically two :cite:p:`gouin_shape_2021`. It should
   therefore be regarded as a parameter of the model, rather than as a measured connectivity.

   The filaments of a branch are drawn the first time that any node on the branch requires them, and stored on every node of
   the branch, from its earliest progenitor up to the node of which it is a secondary progenitor. They are therefore the same
   whether a satellite's host is identified with a node before or after the promotion of its primary progenitor, and the same
   however the nodes of the tree are ordered when they are first required. This matters because virial orbits are assigned at
   several points: during the construction of a tree, by tree operators, and when nodes are initialized.
   </description>
  </darkMatterHaloFilaments>
  !!]
  type, extends(darkMatterHaloFilamentsClass) :: darkMatterHaloFilamentsRandomSequential
     !!{RST
     A dark matter halo filaments class in which filaments are placed at random, subject to a minimum separation.
     !!}
     private
     double precision :: countFilamentsMean     , angleSeparationMinimum, &
          &              cosineSeparationMaximum
     integer          :: filamentAxesID
   contains
     procedure :: axes => randomSequentialAxes
  end type darkMatterHaloFilamentsRandomSequential

  interface darkMatterHaloFilamentsRandomSequential
     !!{RST
     Constructors for the :galacticus-class:`darkMatterHaloFilamentsRandomSequential` dark matter halo filaments class.
     !!}
     module procedure randomSequentialConstructorParameters
     module procedure randomSequentialConstructorInternal
  end interface darkMatterHaloFilamentsRandomSequential

  ! The largest mean number of filaments permitted. Random sequential placement jams (i.e. no further filament can be placed)
  ! when the sphere becomes crowded, so the mean number of filaments is limited to well below the jamming limit for the default
  ! minimum separation.
  double precision, parameter :: countFilamentsMeanMaximum=6.0d0

  ! The number of attempts made to place each further filament before restarting the placement of all filaments of a branch,
  ! and the number of such restarts made before giving up.
  integer         , parameter :: countAttemptsMaximum     =1000 , countRestartsMaximum=100

contains

  function randomSequentialConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`darkMatterHaloFilamentsRandomSequential` dark matter halo filaments class which
    takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterHaloFilamentsRandomSequential)                :: self
    type            (inputParameters                        ), intent(inout) :: parameters
    double precision                                                         :: countFilamentsMean, angleSeparationMinimum

    !![
    <inputParameter docformat="rst">
      <name>countFilamentsMean</name>
      <defaultValue>3.0d0</defaultValue>
      <defaultSource>Calibration of :galacticus-class:`virialOrbitFilamentary` to :cite:t:`libeskind_universal_2014`.</defaultSource>
      <source>parameters</source>
      <description>
      The mean number of filaments, :math:`\langle N \rangle`, feeding each halo. Non-integer values are permitted, in which
      case the number of filaments of each halo is one of the two integers bracketing the mean.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>angleSeparationMinimum</name>
      <defaultValue>60.0d0</defaultValue>
      <defaultSource>Calibration of :galacticus-class:`virialOrbitFilamentary` to :cite:t:`libeskind_universal_2014`.</defaultSource>
      <source>parameters</source>
      <description>
      The minimum angle, in degrees, between any two filaments feeding a halo. Must be less than 90 degrees if
      [countFilamentsMean] exceeds 2.
      </description>
    </inputParameter>
    !!]
    self=darkMatterHaloFilamentsRandomSequential(countFilamentsMean,angleSeparationMinimum)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function randomSequentialConstructorParameters

  function randomSequentialConstructorInternal(countFilamentsMean,angleSeparationMinimum) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`darkMatterHaloFilamentsRandomSequential` dark matter halo filaments class.
    !!}
    use :: Error                   , only : Error_Report
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    type            (darkMatterHaloFilamentsRandomSequential)                :: self
    double precision                                         , intent(in   ) :: countFilamentsMean, angleSeparationMinimum
    !![
    <constructorAssign variables="countFilamentsMean, angleSeparationMinimum"/>
    !!]

    if (countFilamentsMean     <  1.0d0                                                     ) call Error_Report('[countFilamentsMean] must be at least 1'                   //{introspection:location})
    if (countFilamentsMean     >  countFilamentsMeanMaximum                                 ) call Error_Report('[countFilamentsMean] must be no larger than 6'             //{introspection:location})
    if (angleSeparationMinimum <  0.0d0                                                     ) call Error_Report('[angleSeparationMinimum] must be non-negative'             //{introspection:location})
    if (angleSeparationMinimum >  90.0d0                                                    ) call Error_Report('[angleSeparationMinimum] must be no larger than 90 degrees'//{introspection:location})
    if (angleSeparationMinimum >= 90.0d0                    .and. countFilamentsMean > 2.0d0) &
         & call Error_Report('[angleSeparationMinimum] must be less than 90 degrees if [countFilamentsMean] exceeds 2, so that further filaments can clear the first'//{introspection:location})
    self%cosineSeparationMaximum=cos(angleSeparationMinimum*Pi/180.0d0)
    !![
    <addMetaProperty component="basic" name="filamentAxes" id="self%filamentAxesID" rank="1" isEvolvable="no" isCreator="yes"/>
    !!]
    return
  end function randomSequentialConstructorInternal

  function randomSequentialAxes(self,node) result(axes)
    !!{RST
    Return the directions of the filaments feeding the halo of ``node``, drawing them for the branch of ``node`` if they do not
    yet exist.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentBasic
    implicit none
    double precision                                         , allocatable  , dimension(:,:) :: axes
    class           (darkMatterHaloFilamentsRandomSequential), intent(inout)                 :: self
    type            (treeNode                               ), intent(inout), target         :: node
    type            (treeNode                               ), pointer                       :: nodeBranch
    class           (nodeComponentBasic                     ), pointer                       :: basic
    double precision                                         , allocatable  , dimension(:  ) :: directionCosines
    double precision                                         , allocatable  , dimension(:,:) :: directions
    double precision                                                        , dimension(3  ) :: direction
    integer                                                                                  :: countFilaments  , countPlaced  , &
         &                                                                                      countAttempts   , countRestarts
    logical                                                                                  :: jammed

    basic            => node %basic                    (                   )
    directionCosines =  basic%floatRank1MetaPropertyGet(self%filamentAxesID)
    if (size(directionCosines) == 0) then
       ! No filaments exist for this node. Find the earliest progenitor on its branch, which holds the filaments of the branch if
       ! any node of the branch has been assigned them.
       nodeBranch => node
       do while (associated(nodeBranch%firstChild))
          nodeBranch => nodeBranch%firstChild
       end do
       basic            => nodeBranch%basic                    (                   )
       directionCosines =  basic     %floatRank1MetaPropertyGet(self%filamentAxesID)
       if (size(directionCosines) == 0) then
          ! The branch has no filaments, so draw them. First, choose the number of filaments, with the smallest variance
          ! consistent with the mean.
          countFilaments=int(self%countFilamentsMean)
          if (nodeBranch%hostTree%randomNumberGenerator_%uniformSample() < self%countFilamentsMean-dble(countFilaments)) &
               & countFilaments=countFilaments+1
          allocate(directions(3,countFilaments))
          countRestarts=0
          do while (.true.)
             ! Place a single filament passing through the halo, i.e. a pair of opposite directions, oriented isotropically.
             directions(:,1)=randomDirection()
             countPlaced    =1
             if (countFilaments > 1) then
                directions(:,2)=-directions(:,1)
                countPlaced    =2
             end if
             ! Place further filaments in isotropically-distributed directions, accepting each only if it is sufficiently
             ! separated from all filaments already placed.
             jammed=.false.
             do while (countPlaced < countFilaments)
                countAttempts=0
                do while (.true.)
                   direction=randomDirection()
                   if (all(matmul(direction,directions(:,1:countPlaced)) <= self%cosineSeparationMaximum)) exit
                   countAttempts=countAttempts+1
                   if (countAttempts >= countAttemptsMaximum) then
                      jammed=.true.
                      exit
                   end if
                end do
                if (jammed) exit
                countPlaced              =countPlaced+1
                directions(:,countPlaced)=direction
             end do
             if (.not.jammed) exit
             ! Placement jammed, so restart it.
             countRestarts=countRestarts+1
             if (countRestarts >= countRestartsMaximum) call Error_Report('unable to place filaments with the required separation - reduce [countFilamentsMean] or [angleSeparationMinimum]'//{introspection:location})
          end do
          directionCosines=reshape(directions,[3*countFilaments])
       end if
       ! Store the filaments on every node of the branch which lacks them, from its earliest progenitor upward.
       do while (associated(nodeBranch))
          basic => nodeBranch%basic()
          if (size(basic%floatRank1MetaPropertyGet(self%filamentAxesID)) == 0) &
               & call basic%floatRank1MetaPropertySet(self%filamentAxesID,directionCosines)
          if (.not.nodeBranch%isPrimaryProgenitor()) exit
          nodeBranch => nodeBranch%parent
       end do
    end if
    axes=reshape(directionCosines,[3,size(directionCosines)/3])
    return

  contains

    function randomDirection()
      !!{RST
      Return a direction drawn isotropically at random.
      !!}
      use :: Numerical_Constants_Math, only : Pi
      implicit none
      double precision, dimension(3) :: randomDirection
      double precision               :: cosineTheta    , phi

      cosineTheta    =+2.0d0   *nodeBranch%hostTree%randomNumberGenerator_%uniformSample()-1.0d0
      phi            =+2.0d0*Pi*nodeBranch%hostTree%randomNumberGenerator_%uniformSample()
      randomDirection=[                                     &
           &           sqrt(1.0d0-cosineTheta**2)*cos(phi), &
           &           sqrt(1.0d0-cosineTheta**2)*sin(phi), &
           &                      cosineTheta               &
           &          ]
      return
    end function randomDirection

  end function randomSequentialAxes
