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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a galactic inclination class deriving the inclination from the halo angular momentum vector.
  !!}

  use :: Geometry_Lightcones, only : geometryLightconeClass
  use :: ISO_Varying_String , only : varying_string

  !![
  <galacticInclination name="galacticInclinationSpinVector" docformat="rst">
   <description>
   A galactic inclination class in which the inclination is derived from the direction of the halo angular momentum
   vector, on the assumption that the galaxy's symmetry axis is aligned with it. The inclination is then

   .. math::

      i = \cos^{-1} \left| \hat{\mathbf{j}} \cdot \hat{\mathbf{d}} \right|,

   with :math:`\hat{\mathbf{j}}` the unit angular momentum vector and :math:`\hat{\mathbf{d}}` the direction to the
   observer. The absolute value folds the angle into :math:`[0,\pi/2]`: a disk seen from "above" and "below" is
   equally inclined.

   This requires no storage of its own, and it is the physically motivated choice: orientations come out correlated
   with the surrounding large-scale structure rather than drawn independently, which matters for lightcones and mock
   catalogs where neighboring galaxies should not have independent orientations. It requires the ``spin`` component
   to track the angular momentum *vector*---the ``vector`` implementation---and reports an error at construction if
   it does not.

   The direction to the observer is set by ``lineOfSight``:

   ``fixed``
     A direction fixed in the simulation box, given by the ``direction`` parameter. The default is the
     :math:`x`-axis. Appropriate for a snapshot, where there is no preferred observer.

   ``lightcone``
     The direction from the observer to the galaxy, taken from a :galacticus-class:`geometryLightconeClass` object as
     the position at which the galaxy crosses the lightcone. This is the correct choice when building a mock
     catalog, where the line of sight varies across the sky.

     Note that a galaxy may appear more than once in a lightcone through replication of the simulation box, and that
     each appearance strictly has its own line of sight. The lightcone crossing position is a single well-defined
     choice among these, and is used for every appearance; the alternative would require the replication instance to
     be threaded through every consumer of this class.
   </description>
  </galacticInclination>
  !!]

  !![
  <enumeration docformat="rst">
   <name>spinVectorLineOfSight</name>
   <description>
   Enumerates the ways in which the direction to the observer may be determined in the ``spinVector`` galactic
   inclination class.
   </description>
   <encodeFunction>yes</encodeFunction>
   <validator>yes</validator>
   <visibility>public</visibility>
   <entry label="fixed"    />
   <entry label="lightcone"/>
  </enumeration>
  !!]

  type, extends(galacticInclinationClass) :: galacticInclinationSpinVector
     !!{RST
     A galactic inclination class deriving the inclination from the halo angular momentum vector.
     !!}
     private
     class           (geometryLightconeClass                ), pointer                   :: geometryLightcone_ => null()
     type            (enumerationSpinVectorLineOfSightType  )                            :: lineOfSight
     double precision                                        , dimension(3)              :: direction
   contains
     final     ::                spinVectorDestructor
     procedure :: inclination => spinVectorInclination
  end type galacticInclinationSpinVector

  interface galacticInclinationSpinVector
     !!{RST
     Constructors for the :galacticus-class:`galacticInclinationSpinVector` galactic inclination class.
     !!}
     module procedure spinVectorConstructorParameters
     module procedure spinVectorConstructorInternal
  end interface galacticInclinationSpinVector

contains

  function spinVectorConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`galacticInclinationSpinVector` galactic inclination class which takes a
    parameter set as input.
    !!}
    use :: Input_Parameters   , only : inputParameter, inputParameters
    use :: ISO_Varying_String  , only : char          , var_str
    implicit none
    type            (galacticInclinationSpinVector)                :: self
    type            (inputParameters              ), intent(inout) :: parameters
    class           (geometryLightconeClass       ), pointer       :: geometryLightcone_
    type            (varying_string               )                :: lineOfSight
    double precision                               , dimension(3)  :: direction

    !![
    <inputParameter docformat="rst">
      <name>lineOfSight</name>
      <defaultValue>var_str('fixed')</defaultValue>
      <description>
      How the direction to the observer is determined: ``fixed`` for a direction fixed in the simulation box, or
      ``lightcone`` to take it from a lightcone geometry.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>direction</name>
      <defaultValue>[1.0d0,0.0d0,0.0d0]</defaultValue>
      <description>
      The direction to the observer, as a Cartesian vector in the coordinate system of the simulation box. It need
      not be normalized. Used only when ``lineOfSight`` is ``fixed``.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    nullify(geometryLightcone_)
    if (enumerationSpinVectorLineOfSightEncode(char(lineOfSight),includesPrefix=.false.) == spinVectorLineOfSightLightcone) then
       !![
       <objectBuilder class="geometryLightcone" name="geometryLightcone_" source="parameters"/>
       !!]
    end if
    self=galacticInclinationSpinVector(enumerationSpinVectorLineOfSightEncode(char(lineOfSight),includesPrefix=.false.),direction,geometryLightcone_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    if (associated(geometryLightcone_)) then
       !![
       <objectDestructor name="geometryLightcone_"/>
       !!]
    end if
    return
  end function spinVectorConstructorParameters

  function spinVectorConstructorInternal(lineOfSight,direction,geometryLightcone_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`galacticInclinationSpinVector` galactic inclination class.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : defaultSpinComponent
    implicit none
    type            (galacticInclinationSpinVector       )                                  :: self
    type            (enumerationSpinVectorLineOfSightType), intent(in   )                   :: lineOfSight
    double precision                                      , intent(in   ) , dimension(3)    :: direction
    class           (geometryLightconeClass              ), intent(in   ) , target, optional :: geometryLightcone_
    double precision                                                                        :: directionNorm
    !![
    <constructorAssign variables="lineOfSight, direction, *geometryLightcone_"/>
    !!]

    ! The angular momentum vector is needed, so the spin component must track it.
    if (.not.defaultSpinComponent%angularMomentumVectorIsGettable())                                                                     &
         & call Error_Report(                                                                                                            &
         &                   'the `spinVector` galactic inclination requires the angular momentum vector, so the `spin` component must'// &
         &                   ' be set to `vector`'                                                                                     // &
         &                   {introspection:location}                                                                                    &
         &                  )
    select case (self%lineOfSight%ID)
    case (spinVectorLineOfSightFixed    %ID)
       directionNorm=sqrt(sum(self%direction**2))
       if (directionNorm <= 0.0d0) call Error_Report('`direction` must be a non-zero vector'//{introspection:location})
       ! Store the direction normalized, so that the dot product below is a cosine directly.
       self%direction=+self%direction  &
            &         /     directionNorm
    case (spinVectorLineOfSightLightcone%ID)
       if (.not.associated(self%geometryLightcone_)) call Error_Report('a `geometryLightcone` is required when `lineOfSight` is `lightcone`'//{introspection:location})
    case default
       call Error_Report('unrecognized `lineOfSight`'//{introspection:location})
    end select
    return
  end function spinVectorConstructorInternal

  subroutine spinVectorDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`galacticInclinationSpinVector` galactic inclination class.
    !!}
    implicit none
    type(galacticInclinationSpinVector), intent(inout) :: self

    !![
    <objectDestructor name="self%geometryLightcone_"/>
    !!]
    return
  end subroutine spinVectorDestructor

  double precision function spinVectorInclination(self,node) result(inclination)
    !!{RST
    Return the inclination implied by the angular momentum vector of the halo.
    !!}
    use :: Error           , only : Error_Report
    use :: Galacticus_Nodes, only : nodeComponentSpin
    implicit none
    class           (galacticInclinationSpinVector), intent(inout)         :: self
    type            (treeNode                     ), intent(inout), target :: node
    class           (nodeComponentSpin            ), pointer               :: spin
    double precision                               , dimension(3)          :: angularMomentum, direction
    double precision                                                       :: normAngularMomentum, normDirection

    spin               => node%spin                 ()
    angularMomentum    =  spin%angularMomentumVector()
    normAngularMomentum=  sqrt(sum(angularMomentum**2))
    ! A halo with no angular momentum has no defined symmetry axis. Report it as face-on: this is the orientation
    ! which minimizes the effect of any flattened dust distribution, and so is the conservative choice.
    if (normAngularMomentum <= 0.0d0) then
       inclination=0.0d0
       return
    end if
    select case (self%lineOfSight%ID)
    case (spinVectorLineOfSightFixed    %ID)
       ! Already normalized at construction.
       direction=self%direction
    case (spinVectorLineOfSightLightcone%ID)
       direction    =self%geometryLightcone_%positionLightconeCrossing(node)
       normDirection=sqrt(sum(direction**2))
       ! A galaxy at the observer's position has no defined line of sight; treat it as face-on, as above.
       if (normDirection <= 0.0d0) then
          inclination=0.0d0
          return
       end if
       direction=+direction     &
            &    /normDirection
    case default
       direction=0.0d0
       call Error_Report('unrecognized `lineOfSight`'//{introspection:location})
    end select
    ! Fold into [0,pi/2]: only the angle between the axis and the line of sight matters, not its sign. The argument
    ! is clamped as protection against a dot product exceeding unity through round-off.
    inclination=+acos(                                     &
         &            min(                                 &
         &                +1.0d0                         , &
         &                +abs(                            &
         &                     +sum(                       &
         &                          +angularMomentum       &
         &                          *direction              &
         &                         )                       &
         &                     /normAngularMomentum         &
         &                    )                            &
         &               )                                 &
         &           )
    return
  end function spinVectorInclination
