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

  !!{RST
  Implements a property extractor class for the projected density at a set of radii.
  !!}
  use :: Dark_Matter_Halo_Scales             , only : darkMatterHaloScale, darkMatterHaloScaleClass
  use :: Galactic_Structure_Radii_Definitions, only : radiusDefinitions

  !![
  <nodePropertyExtractor name="nodePropertyExtractorProjectedDensity" docformat="rst">
   <description>
   A property extractor class for the projected density at a set of radii. The radii and types of projected density to output is specified by the ``radiusSpecifiers`` parameter. This parameter's value can contain multiple entries, each of which should be a valid :ref:`radius specifier &lt;manual-sec-radiusspecifiers&gt;`.

   A radius specifier can evaluate to zero---most often because the component on which it is based is absent or empty in the node in question. The line of sight then passes through the center of the mass distribution, where the projected density is finite only if the central logarithmic slope of the density profile exceeds :math:`-1` (an :term:`NFW` profile, for example, gives a logarithmically divergent integral). Where it is finite it is evaluated, by adding to the numerical integral an analytic estimate of the contribution from within some small radius, assuming that the density there follows the asymptotic central power law, and reducing that radius until the result no longer depends on it. Where it is not finite, a radius of zero is by default reported as a fatal error naming the offending specifier; setting ``zeroRadiusIsFatal`` to ``false`` instead writes the undefined-value sentinel in place of the projected density.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorArray) :: nodePropertyExtractorProjectedDensity
     !!{RST
     A property extractor class for the projected density at a set of radii.
     !!}
     private
     class  (darkMatterHaloScaleClass), pointer                   :: darkMatterHaloScale_ => null()
     integer                                                      :: radiiCount                    , elementCount_
     logical                                                      :: includeRadii                  , tolerateIntegrationFailures, &
          &                                                          zeroRadiusIsFatal
     type   (varying_string          ), allocatable, dimension(:) :: radiusSpecifiers
     type   (radiusDefinitions       )                            :: radii
   contains
     final     ::                       projectedDensityDestructor
     procedure :: columnDescriptions => projectedDensityColumnDescriptions
     procedure :: size               => projectedDensitySize
     procedure :: elementCount       => projectedDensityElementCount
     procedure :: extract            => projectedDensityExtract
     procedure :: names              => projectedDensityNames
     procedure :: descriptions       => projectedDensityDescriptions
     procedure :: unitsInSI          => projectedDensityUnitsInSI
     procedure :: units       => projectedDensityUnits
  end type nodePropertyExtractorProjectedDensity

  interface nodePropertyExtractorProjectedDensity
     !!{RST
     Constructors for the :galacticus-class:`nodePropertyExtractorProjectedDensity` property extractor class.
     !!}
     module procedure projectedDensityConstructorParameters
     module procedure projectedDensityConstructorInternal
  end interface nodePropertyExtractorProjectedDensity

  ! Module-scope variables used in integrands.
  double precision :: radius_
  !$omp threadprivate(radius_)

contains

  function projectedDensityConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`nodePropertyExtractorProjectedDensity` property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (nodePropertyExtractorProjectedDensity)                              :: self
    type   (inputParameters                      ), intent(inout)               :: parameters
    type   (varying_string                       ), allocatable  , dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass             ), pointer                     :: darkMatterHaloScale_
    logical                                                                     :: includeRadii        , tolerateIntegrationFailures, &
         &                                                                         zeroRadiusIsFatal

    allocate(radiusSpecifiers(parameters%count('radiusSpecifiers')))
    !![
    <inputParameter docformat="rst">
      <name>radiusSpecifiers</name>
      <description>
      A list of radius specifiers at which to output the projected density profile.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>includeRadii</name>
      <defaultValue>.false.</defaultValue>
      <description>
      Specifies whether or not the radii at which projected density data are output should also be included in the output file.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>tolerateIntegrationFailures</name>
      <defaultValue>.false.</defaultValue>
      <description>
      Specifies whether or not failures in integration of the projected density should be tolerated.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>zeroRadiusIsFatal</name>
      <defaultValue>.true.</defaultValue>
      <description>
      Specifies whether a radius specifier which evaluates to zero, in a node in which the projected density diverges at zero
      radius---that is, one in which the central logarithmic slope of the density profile is :math:`-1` or steeper---should be
      reported as a fatal error. If ``false``, the undefined-value sentinel is written in place of the projected density
      instead. A zero radius in a node in which the projected density is finite is always evaluated, and is unaffected by this
      parameter.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    !!]
    self=nodePropertyExtractorProjectedDensity(radiusSpecifiers,includeRadii,tolerateIntegrationFailures,zeroRadiusIsFatal,darkMatterHaloScale_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    !!]
    return
  end function projectedDensityConstructorParameters

  function projectedDensityConstructorInternal(radiusSpecifiers,includeRadii,tolerateIntegrationFailures,zeroRadiusIsFatal,darkMatterHaloScale_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`nodePropertyExtractorProjectedDensity` property extractor class.
    !!}
    implicit none
    type   (nodePropertyExtractorProjectedDensity)                              :: self
    type   (varying_string                       ), intent(in   ), dimension(:) :: radiusSpecifiers
    class  (darkMatterHaloScaleClass             ), intent(in   ), target       :: darkMatterHaloScale_
    logical                                       , intent(in   )               :: includeRadii        , tolerateIntegrationFailures, &
         &                                                                         zeroRadiusIsFatal
    !![
    <constructorAssign variables="radiusSpecifiers, includeRadii, tolerateIntegrationFailures, zeroRadiusIsFatal, *darkMatterHaloScale_"/>
    !!]

    if (includeRadii) then
       self%elementCount_=2
    else
       self%elementCount_=1
    end if
    self%radiiCount      =size(radiusSpecifiers)
    call self%radii%decode(radiusSpecifiers)
    return
  end function projectedDensityConstructorInternal

  subroutine projectedDensityDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`nodePropertyExtractorProjectedDensity` property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorProjectedDensity), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    !!]
    return
  end subroutine projectedDensityDestructor

  integer function projectedDensityElementCount(self,time)
    !!{RST
    Return the number of elements in the ``projectedDensity`` property extractors.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedDensityElementCount=self%elementCount_
    return
  end function projectedDensityElementCount

  function projectedDensitySize(self,time)
    !!{RST
    Return the number of array elements in the ``projectedDensity`` property extractors.
    !!}
    implicit none
    integer         (c_size_t                             )                :: projectedDensitySize
    class           (nodePropertyExtractorProjectedDensity), intent(inout) :: self
    double precision                                       , intent(in   ) :: time
    !$GLC attributes unused :: time

    projectedDensitySize=self%radiiCount
    return
  end function projectedDensitySize

  function projectedDensityExtract(self,node,time,instance) result(densityProjected)
    !!{RST
    Implement a ``projectedDensity`` property extractor.
    !!}
    use :: Galactic_Structure_Radii_Definitions, only : radiusResolver       , radiusUndefined
    use :: Galacticus_Nodes                    , only : treeNode
    use :: Numerical_Integration               , only : integrator           , GSL_Integ_Gauss15
    use :: Numerical_Comparison                , only : Values_Agree
    use :: Mass_Distributions                  , only : massDistributionClass
    use :: Coordinates                         , only : coordinateSpherical  , assignment(=)
    use :: Error                               , only : Error_Report
    implicit none
    double precision                                       , dimension(:,:), allocatable :: densityProjected
    class           (nodePropertyExtractorProjectedDensity), intent(inout) , target      :: self
    type            (treeNode                             ), intent(inout) , target      :: node
    double precision                                       , intent(in   )               :: time
    type            (multiCounter                         ), intent(inout) , optional    :: instance
    class           (massDistributionClass                ), pointer                     :: massDistribution_
    double precision                                       , parameter                   :: toleranceRelative        =1.0d-2, epsilonSingularity   =1.0d-3
    ! Parameters controlling the line of sight integral through the center of the distribution: the radius at which that
    ! integral begins, as a fraction of the outer radius, the factor by which that radius is reduced at each refinement, and the
    ! maximum number of such refinements permitted before convergence is declared to have failed.
    double precision                                       , parameter                   :: radiusInnerFractional    =1.0d-3, factorRadiusInner    =1.0d+1
    integer                                                , parameter                   :: refinementCountMaximum   =20
    type            (radiusResolver                       )                              :: resolver
    type            (integrator                           )                              :: integrator_
    integer                                                                              :: i                               , status                      , &
         &                                                                                  countRefinements                , countAgreements
    double precision                                                                     :: radiusVirial                    , radiusOuter                 , &
         &                                                                                  radiusInner                     , radiusInnerPrevious         , &
         &                                                                                  densityProjectedPrevious        , densityProjectedInner       , &
         &                                                                                  densityProjectedNumerical       , toleranceAbsolute
    logical                                                                              :: converged
    type            (coordinateSpherical                  )                              :: coordinates
    !$GLC attributes unused :: time, instance

    allocate(densityProjected(self%radiiCount,self%elementCount_))
    resolver                                                   =  radiusResolver(self%radii,node,self%darkMatterHaloScale_,radiusVirialRequired=.true.)
    radiusVirial                                               =  resolver%radiusVirial
    integrator_=integrator(projectedDensityIntegrand,toleranceRelative=1.0d-3,hasSingularities=.true.,integrationRule=GSL_Integ_Gauss15)
    do i=1,self%radiiCount
       call resolver%evaluate(i,radius_)
       if (radius_ < 0.0d0) then
          ! The radius is undefined in this node - report the sentinel. Note that this test must precede any arithmetic on the
          ! radius, since the sentinel would overflow (and so trap) under multiplication.
          densityProjected       (i,1)=radiusUndefined
          if (self%includeRadii)                       &
               & densityProjected(i,2)=radiusUndefined
          cycle
       end if
       massDistribution_ => node%massDistribution(self%radii%specifiers(i)%component,self%radii%specifiers(i)%mass)
       if     (                                                              &
            &   radius_                                            <=  0.0d0 &
            &  .and.                                                         &
            &   massDistribution_%densitySlopeLogarithmicCentral() <= -1.0d0 &
            & ) then
          ! The radius is zero, so the line of sight passes through the center of the mass distribution. The integral along that
          ! line of sight converges only if the central logarithmic slope of the density profile exceeds -1, which it does not
          ! here - an unknown central slope, being the most negative representable value, is rejected by the same test - so the
          ! projected density does not exist in this node.
          if (self%zeroRadiusIsFatal) &
               & call resolver%reportZeroRadius(i,'projectedDensity','the line of sight integral through the center of this mass distribution diverges')
          densityProjected       (i,1)=radiusUndefined
          if (self%includeRadii)                       &
               & densityProjected(i,2)=radius_
          !![
          <objectDestructor name="massDistribution_"/>
          !!]
          cycle
       end if
       radiusOuter              =  max(radius_*2.0d0,radiusVirial)
       ! Set an absolute tolerance scale for projected density convergence that is a small fraction of the mean halo density,
       ! integrated over a path length of 1 Mpc.
       toleranceAbsolute        =  +toleranceRelative                           &
            &                      *self%darkMatterHaloScale_%densityMean(node)
       if (radius_ > 0.0d0) then
          ! Cut out a small region round the coordinate singularity at the inner radius. This region will be integrated
          ! analytically assuming a constant density over this region. The region outside of this cut-out will be integrated
          ! numerically.
          radiusInner=min(radius_*(1.0d0+epsilonSingularity),radiusOuter)
       else
          ! The line of sight passes through the center of the distribution. There is no coordinate singularity in this case,
          ! but the integral is performed in the logarithm of radius and so can not begin at zero radius. Begin it instead at a
          ! small radius, and add an analytic estimate of the contribution from within that radius. As that radius is
          ! subsequently reduced until the result no longer depends on it, this initial choice need not resolve the scale radii
          ! of the distribution---it need only provide a starting point.
          radiusInner=radiusOuter*radiusInnerFractional
       end if
       !! Analytic integral within the innermost radius.
       densityProjectedInner    =  projectedDensityInner(radiusInner)
       !! Numerical integral outside of the innermost radius.
       densityProjectedNumerical  =0.0d0
       if (radiusInner < radiusOuter) then
          densityProjectedPrevious=0.0d0
          converged               =.false.
          do while (.not.converged)
             densityProjectedNumerical=integrator_%integrate(log(radiusInner),log(radiusOuter),status=status)
             if (status /= errorStatusSuccess .and. .not.self%tolerateIntegrationFailures) &
                  & call Error_Report('integration of projected density failed'//{introspection:location})
             converged                =Values_Agree(densityProjectedNumerical,densityProjectedPrevious,relTol=toleranceRelative,absTol=toleranceAbsolute)
             if (.not.converged) then
                radiusOuter             =2.0d0*radiusOuter
                densityProjectedPrevious=      densityProjectedNumerical
             end if
          end do
       end if
       densityProjected(i,1)=+densityProjectedInner     &
            &                +densityProjectedNumerical
       if (radius_ <= 0.0d0 .and. radiusInner > 0.0d0) then
          ! Reduce the radius at which the numerical integral begins, integrating each newly-exposed shell and re-estimating the
          ! contribution from within the new innermost radius, until the total is insensitive to that radius. This tests the
          ! accuracy of the analytic estimate directly, and so requires no knowledge of the scale radii of the distribution: the
          ! estimate becomes exact as the radius approaches zero, but converges as soon as the density has reached its
          ! asymptotic, central power law form. Two successive agreements are required, so that the slow variation of a poor
          ! estimate is not mistaken for convergence.
          countAgreements =0
          countRefinements=0
          do while (countAgreements < 2)
             radiusInnerPrevious      =                           radiusInner
             radiusInner              =                           radiusInner/factorRadiusInner
             densityProjectedNumerical=+densityProjectedNumerical                                                      &
                  &                    +integrator_%integrate(log(radiusInner),log(radiusInnerPrevious),status=status)
             if (status /= errorStatusSuccess .and. .not.self%tolerateIntegrationFailures) &
                  & call Error_Report('integration of projected density failed'//{introspection:location})
             densityProjectedPrevious =+densityProjected     (i          ,1)
             densityProjected(i,1)    =+projectedDensityInner(radiusInner   ) &
                  &                    +densityProjectedNumerical
             if (Values_Agree(densityProjected(i,1),densityProjectedPrevious,relTol=toleranceRelative,absTol=toleranceAbsolute)) then
                countAgreements=countAgreements+1
             else
                countAgreements=0
             end if
             countRefinements=countRefinements+1
             if (countRefinements >= refinementCountMaximum .and. countAgreements < 2) then
                if (.not.self%tolerateIntegrationFailures) &
                     & call Error_Report('projected density at zero radius failed to converge'//{introspection:location})
                exit
             end if
          end do
       end if
       if (self%includeRadii) densityProjected(i,2)=radius_
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
    end do
    return

  contains

    double precision function projectedDensityInner(radius)
      !!{RST
      Return the contribution to the projected density from within the given ``radius``---that is, the part of the line of sight
      integral which is not evaluated numerically.

      For a line of sight which passes to one side of the center the density is assumed to be constant, at its value at the
      radius of closest approach, over the small region cut out around the coordinate singularity there.

      For a line of sight which passes through the center there is no singularity to cut out, but the integral can not begin at
      zero radius. The contribution from within ``radius`` is instead estimated by assuming that the density follows its
      asymptotic, central power law, :math:`\rho \propto r^\gamma`, there:

      .. math::

         2 \int_0^r \rho(r^\prime) \mathrm{d}r^\prime = \frac{2 r \rho(r)}{1+\gamma},

      which is finite for :math:`\gamma > -1`---precisely the condition under which the full integral converges, and which is
      guaranteed to hold here as any distribution failing it has already been rejected.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      if (radius_ > 0.0d0) then
         coordinates          =[radius_,0.0d0,0.0d0]
         projectedDensityInner=+2.0d0                                  &
              &                *sqrt(                                  &
              &                      +radius **2                       &
              &                      -radius_**2                       &
              &                     )                                  &
              &                *massDistribution_%density(coordinates)
      else if (radius <= 0.0d0) then
         ! The distribution has no extent, so there is nothing to integrate.
         projectedDensityInner=+0.0d0
      else
         coordinates          =[radius ,0.0d0,0.0d0]
         projectedDensityInner=+2.0d0                                                           &
              &                *radius                                                          &
              &                *  massDistribution_%density                       (coordinates) &
              &                /(                                                               &
              &                  +1.0d0                                                         &
              &                  +massDistribution_%densitySlopeLogarithmicCentral(           ) &
              &                 )
      end if
      return
    end function projectedDensityInner

    double precision function projectedDensityIntegrand(radiusLogarithmic)
      !!{RST
      Integrand function used for computing projected densities.
      !!}
      implicit none
      double precision, intent(in   ) :: radiusLogarithmic
      double precision                :: radius

      radius=exp(radiusLogarithmic)
      if (radius <= radius_) then
         projectedDensityIntegrand=+0.0d0
      else
         coordinates=[radius,0.0d0,0.0d0]
         projectedDensityIntegrand=+2.0d0                                  &
              &                    *radius       **2                       &
              &                    /sqrt(                                  &
              &                          +radius **2                       &
              &                          -radius_**2                       &
              &                    )                                       &
              &                    *massDistribution_%density(coordinates)
      end if
      return
    end function projectedDensityIntegrand

  end function projectedDensityExtract

  subroutine projectedDensityNames(self,names,time)
    !!{RST
    Return the names of the ``projectedDensity`` properties.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: names
    !$GLC attributes unused :: time

    allocate(names(self%elementCount_))
    names       (1)="projectedDensity"
    if (self%includeRadii)                   &
         & names(2)="projectedDensityRadius"
    return
  end subroutine projectedDensityNames

  subroutine projectedDensityDescriptions(self,descriptions,time)
    !!{RST
    Return descriptions of the ``projectedDensity`` property.
    !!}
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout)                             :: self
    double precision                                       , intent(in   ), optional                   :: time
    type            (varying_string                       ), intent(inout), dimension(:) , allocatable :: descriptions
    !$GLC attributes unused :: time

    allocate(descriptions(self%elementCount_))
    descriptions       (1)="Projected density at a given radius [M☉/Mpc⁻²]."
    if (self%includeRadii)                                                      &
         & descriptions(2)="Radius at which projected density is output [Mpc]."
    return
  end subroutine projectedDensityDescriptions

  subroutine projectedDensityColumnDescriptions(self,descriptions,values,valuesDescription,valuesUnits,time)
    !!{RST
    Return column descriptions of the ``projectedDensity`` property.
    !!}
    use            :: Units_MetaData, only : unitType
    use, intrinsic :: ISO_C_Binding , only : c_int
    implicit none
    class           (nodePropertyExtractorProjectedDensity), intent(inout)                            :: self
    double precision                                       , intent(in   ), optional                  :: time
    type            (varying_string                       ), intent(inout), dimension(:), allocatable :: descriptions
    double precision                                       , intent(inout), dimension(:), allocatable :: values
    type            (varying_string                       ), intent(  out)                            :: valuesDescription
    type            (unitType                             ), intent(  out)                            :: valuesUnits
    !$GLC attributes unused :: time

    allocate(descriptions(self%radiiCount))
    allocate(values      (              0))
    valuesDescription=var_str('')
    valuesUnits      =unitType(1.0d0)
    descriptions     =self%radii%specifiers%name
    return
  end subroutine projectedDensityColumnDescriptions

  function projectedDensityUnitsInSI(self,time)
    !!{RST
    Return the units of the ``projectedDensity`` properties in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : massSolar, megaParsec
    implicit none
    double precision                                       , allocatable  , dimension(:) :: projectedDensityUnitsInSI
    class           (nodePropertyExtractorProjectedDensity), intent(inout)               :: self
    double precision                                       , intent(in   ), optional     :: time
    !$GLC attributes unused :: time

    allocate(projectedDensityUnitsInSI(self%elementCount_))
    projectedDensityUnitsInSI       (1)=massSolar/megaParsec**2
    if (self%includeRadii)                                       &
         & projectedDensityUnitsInSI(2)=          megaParsec
    return
  end function projectedDensityUnitsInSI

  function projectedDensityUnits(self,time) result(units)
    !!{RST
    Return the units of the projectedDensity properties.
    !!}
    use :: Units_MetaData, only : unitType
    implicit none
    type            (unitType                             ), dimension(:), allocatable :: units
    class           (nodePropertyExtractorProjectedDensity), intent(inout)             :: self
    double precision                                       , intent(in   ), optional   :: time
    double precision                                       , dimension(:), allocatable :: siValues
    integer                                                                            :: i

    siValues=self%unitsInSI(time)
    allocate(units(size(siValues)))
    do i=1,size(siValues)
       units(i)=unitType(siValues(i),description='M☉/Mpc²',quantity='solMass/Mpc**2')
    end do
    return
  end function projectedDensityUnits
