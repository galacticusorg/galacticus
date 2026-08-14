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

  !+    Contributions to this file made by: Andrew Benson. The pinning of the tabulation ranges to absolute lattices for issue
  !+    #1317 was drafted with assistance from Claude, and reviewed and verified by Andrew Benson.

  !!{RST
  Implementation of an abstract mass distribution class for tabulated spherically symmetric distributions.
  !!}

  use, intrinsic :: ISO_C_Binding   , only : c_size_t
  use            :: Numerical_Ranges, only : rangeLattice

  !![
  <massDistribution name="massDistributionSphericalTabulated" abstract="yes" docformat="rst">
   <description>
   An abstract mass distribution class for tabulated spherically symmetric distributions.
   </description>
  </massDistribution>
  !!]
  type, extends(massDistributionSpherical), abstract :: massDistributionSphericalTabulated
     !!{RST
     Implementation of an abstract mass distribution class for tabulated spherically symmetric distributions.
     !!}
     private
   contains
     !![
     <methods docformat="rst">
       <method method="parameters"           description="Return parameters of the current mass distribution."                  />
       <method method="factoryTabulation"    description="Return an instance of the class with the given tabulation parameters."/>
       <method method="suffix"               description="Return a suffix to append to table file names."                       />
       <method method="fileRead"             description="Read tabulation data from file."                                      />
       <method method="fileWrite"            description="Write tabulation data to file."                                       />
       <method method="tabulate"             description="(Re)tabulate the mass distribution."                                  />
       <method method="interpolate"          description="Interpolate in the mass distribution."                                />
       <method method="isTabulating"         description="Return true if the thread associated with the object is tabulating."  />
       <method method="velocityDispersion1D" description="Compute the 1D velocity dispersion at the given coordinates."         />
     </methods>
     !!]
     procedure(sphericalTabulatedParameters       ), deferred :: parameters
     procedure(sphericalTabulatedFactoryTabulation), deferred :: factoryTabulation
     procedure(sphericalTabulatedSuffix           ), deferred :: suffix
     procedure                                                :: fileRead                   => sphericalTabulatedFileRead
     procedure                                                :: fileWrite                  => sphericalTabulatedFileWrite
     procedure                                                :: tabulate                   => sphericalTabulatedTabulate
     procedure                                                :: interpolate                => sphericalTabulatedInterpolate
     procedure                                                :: isTabulating               => sphericalTabulatedIsTabulating
     procedure                                                :: massEnclosedBySphere       => sphericalTabulatedMassEnclosedBySphere
     procedure                                                :: radiusEnclosingDensity     => sphericalTabulatedRadiusEnclosingDensity
     procedure                                                :: potential                  => sphericalTabulatedPotential
     procedure                                                :: potentialDifference        => sphericalTabulatedPotentialDifference
     procedure                                                :: fourierTransform           => sphericalTabulatedFourierTransform
     procedure                                                :: energy                     => sphericalTabulatedEnergy
     procedure                                                :: densityRadialMoment        => sphericalTabulatedDensityRadialMoment
     procedure                                                :: radiusFreefall             => sphericalTabulatedRadiusFreefall
     procedure                                                :: radiusFreefallIncreaseRate => sphericalTabulatedRadiusFreefallIncreaseRate
     procedure                                                :: velocityDispersion1D       => sphericalTabulatedVelocityDispersion1D
  end type massDistributionSphericalTabulated

  !![
  <enumeration docformat="rst">
   <name>quantity</name>
   <description>
   Quantities for tabulation mass distributions.
   </description>
   <decodeFunction>yes</decodeFunction>
   <entry label="mass"                      />
   <entry label="radiusEnclosingDensity"    />
   <entry label="potential"                 />
   <entry label="energy"                    />
   <entry label="fourierTransform"          />
   <entry label="radiusFreefall"            />
   <entry label="radiusFreefallIncreaseRate"/>
   <entry label="densityRadialMoment0"      />
   <entry label="densityRadialMoment1"      />
   <entry label="densityRadialMoment2"      />
   <entry label="densityRadialMoment3"      />
   <entry label="velocityDispersion1D"      />
  </enumeration>
  !!]
  
  ! Every axis of these tabulations is pinned to a lattice of points per *octave*, with its bounds rounded outward to whole
  ! octaves. The octave is the natural interval here: the tabulated quantities are functions of scale-free radii and of
  ! dimensionless shape parameters, all of which the tabulation brackets by a factor of two on either side of a request - which
  ! is exactly one anchor interval, so the pinning costs nothing beyond the margin that was already being applied. It is also
  ! the interval on which these tabulations were already, in effect, being built: an accumulated cache of
  ! `SIDM_parametric_profile` tabulations has all four of its axes spanning whole numbers of octaves at exactly six points per
  ! octave, so pinning to this lattice reproduces those tabulations point for point rather than rebuilding them on a shifted
  ! grid. Anchoring to whole decades instead would inflate a three-dimensional table of that size fourfold.
  type :: massDistributionTabulation
     !!{RST
     Object used to store individual mass distribution tabulations.

     The extent of a tabulation is specified by the absolute lattices ``latticeRadius`` and ``latticeParameters``: these are
     the source of truth for it, and the minima, maxima and point counts are derived from them. Since a lattice is fixed by
     its gridding scheme and density of points alone, two tabulations of the same quantity have identical abscissae wherever
     they overlap, no matter what sequence of requests built each.
     !!}
     type            (enumerationQuantityType)                                   :: quantity
     logical                                                                     :: logTransform            =.false.   , isNegative              =.false.
     double precision                                                            :: radiusMinimum           =0.0d0     , radiusMaximum           =0.0d0
     double precision                          , allocatable, dimension(:      ) :: parametersMinimum                  , parametersMaximum
     ! Flags recording, per parameter, that the corresponding bound has reached the hard limit imposed on it, so that a request
     ! beyond that bound cannot be met by extending the tabulation and must instead be met by extrapolation.
     logical                                   , allocatable, dimension(:      ) :: parametersAtLimitMinimum           , parametersAtLimitMaximum
     integer         (c_size_t                )                                  :: radiusCountPer          =0_c_size_t, countRadii              =0_c_size_t
     integer         (c_size_t                ), allocatable, dimension(:      ) :: parametersCountPer                 , countParameters
     type            (rangeLattice            )                                  :: latticeRadius
     type            (rangeLattice            ), allocatable, dimension(:      ) :: latticeParameters
     double precision                          , allocatable, dimension(:,:,:,:) :: table
  end type massDistributionTabulation

  type :: massDistributionContainer
     !!{RST
     Object to store collections of mass distribution tabulations.
     !!}
     type            (varying_string            ), allocatable, dimension(:) :: nameParameters        , descriptionParameters
     double precision                            , allocatable, dimension(:) :: parametersMinimumLimit, parametersMaximumLimit
     ! Tabulations for individual quantities.
     type            (massDistributionTabulation)                            :: mass                      =massDistributionTabulation(quantity=quantityMass                      ,logTransform=.true. )
     type            (massDistributionTabulation)                            :: radiusEnclosingDensity    =massDistributionTabulation(quantity=quantityRadiusEnclosingDensity    ,logTransform=.false.)
     type            (massDistributionTabulation)                            :: potential                 =massDistributionTabulation(quantity=quantityPotential                 ,logTransform=.false.)
     type            (massDistributionTabulation)                            :: energy                    =massDistributionTabulation(quantity=quantityEnergy                    ,logTransform=.false.)
     type            (massDistributionTabulation)                            :: fourierTransform          =massDistributionTabulation(quantity=quantityFourierTransform          ,logTransform=.false.)
     type            (massDistributionTabulation)                            :: radiusFreefall            =massDistributionTabulation(quantity=quantityRadiusFreefall            ,logTransform=.false.)
     type            (massDistributionTabulation)                            :: radiusFreefallIncreaseRate=massDistributionTabulation(quantity=quantityRadiusFreefallIncreaseRate,logTransform=.false.)
     type            (massDistributionTabulation)                            :: densityRadialMoment0      =massDistributionTabulation(quantity=quantityDensityRadialMoment0      ,logTransform=.true. )
     type            (massDistributionTabulation)                            :: densityRadialMoment1      =massDistributionTabulation(quantity=quantityDensityRadialMoment1      ,logTransform=.true. )
     type            (massDistributionTabulation)                            :: densityRadialMoment2      =massDistributionTabulation(quantity=quantityDensityRadialMoment2      ,logTransform=.true. )
     type            (massDistributionTabulation)                            :: densityRadialMoment3      =massDistributionTabulation(quantity=quantityDensityRadialMoment3      ,logTransform=.true. )
     type            (massDistributionTabulation)                            :: velocityDispersion1D      =massDistributionTabulation(quantity=quantityVelocityDispersion1D      ,logTransform=.true. )
   contains
     !![
     <methods docformat="rst">
       <method method="initialize"      description="Initialize the container (specifically the number of parameters)."/>
       <method method="nameParameter"   description="Return the name of the index parameter for a given tabulation."   />
       <method method="countParameters" description="Return the number of parameters for a given tabulation."          />
     </methods>
     !!]
     procedure :: initialize      => massDistributionContainerInitialize
     procedure :: nameParameter   => massDistributionContainerNameParameter
     procedure :: countParameters => massDistributionContainerCountParameters
  end type massDistributionContainer
  
  abstract interface
     !!{RST
     Interfaces to deferred functions.
     !!}
     subroutine sphericalTabulatedParameters(self,densityNormalization,radiusNormalization,parameters,container)
       import massDistributionSphericalTabulated, massDistributionContainer
       class           (massDistributionSphericalTabulated), intent(inout)                              :: self
       double precision                                    , intent(  out)                              :: densityNormalization, radiusNormalization
       double precision                                    , intent(inout), allocatable, dimension(:  ) :: parameters
       type            (massDistributionContainer         ), intent(  out), pointer                     :: container
     end subroutine sphericalTabulatedParameters

     function sphericalTabulatedFactoryTabulation(self,parameters) result(instance)
       import massDistributionSphericalTabulated
       class           (massDistributionSphericalTabulated)               , pointer        :: instance
       class           (massDistributionSphericalTabulated), intent(inout)                 :: self
       double precision                                    , intent(in   ), dimension(:  ) :: parameters
     end function sphericalTabulatedFactoryTabulation
 
     function sphericalTabulatedSuffix(self) result(suffix)
       import massDistributionSphericalTabulated, varying_string
       type (varying_string                    )                :: suffix
       class(massDistributionSphericalTabulated), intent(inout) :: self
     end function sphericalTabulatedSuffix
  end interface

  ! Indicator used to record if tabulation is underway. When tabulating a property, *do not* use tabulated versions of other
  ! properties - we want a precise calculation when building the tables.
  logical :: tabulating=.false.
  !$omp threadprivate(tabulating)
  
contains

  double precision function sphericalTabulatedMassEnclosedBySphere(self,radius) result(mass)
    !!{RST
    Computes the mass enclosed within a sphere of given ``radius`` for spherically-symmetric mass distributions using a tabulation.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout) , target      :: self
    double precision                                    , intent(in   )               :: radius
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               radiusScaled

    if (radius <= 0.0d0) then
       mass=0.0d0
       return
    end if
    if (tabulating) then
       mass=self%massEnclosedBySphereNumerical(radius)
    else
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       radiusScaled=+radius              &
            &       /radiusNormalization
       call self%tabulate(radiusScaled,parameters,container,container%mass)
       ! Perform the interpolation.
       mass   =+self%interpolate         (radiusScaled,parameters,container%mass) &
            &  *     densityNormalization                                         &
            &  *     radiusNormalization **3
    end if
    return
  end function sphericalTabulatedMassEnclosedBySphere

  double precision function sphericalTabulatedRadiusEnclosingDensity(self,density,radiusGuess) result(radius)
    !!{RST
    Computes the radius enclosing the given ``radius`` for spherically-symmetric mass distributions using a tabulation.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout) , target      :: self
    double precision                                    , intent(in   )               :: density
    double precision                                    , intent(in   ) , optional    :: radiusGuess
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               densityScaled

    if (tabulating) then
       radius=self%radiusEnclosingDensityNumerical(density)
    else
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       densityScaled=+density              &
            &        /densityNormalization
       call self%tabulate(densityScaled,parameters,container,container%radiusEnclosingDensity)
       ! Perform the interpolation.
       radius=+self%interpolate        (densityScaled,parameters,container%radiusEnclosingDensity) &
            & *     radiusNormalization
    end if
    return
  end function sphericalTabulatedRadiusEnclosingDensity

  double precision function sphericalTabulatedPotential(self,coordinates,status) result(potential)
    !!{RST
    Compute the potential at the given ``coordinates`` in a spherical mass distribution using a tabulation.
    !!}
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout) , target      :: self
    class           (coordinate                        ), intent(in   )               :: coordinates
    type            (enumerationStructureErrorCodeType ), intent(  out) , optional    :: status
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               radiusScaled

    if (.not.tabulating) then
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       radiusScaled=+coordinates%rSpherical         () &
            &       /            radiusNormalization
       call self%tabulate(radiusScaled,parameters,container,container%potential)
       ! Perform the interpolation.
       potential=+self%interpolate         (radiusScaled,parameters,container%potential) &
            &    *     densityNormalization                                              &
            &    *     radiusNormalization **2
       ! Set return status.
       if (present(status)) status=structureErrorCodeSuccess
    else
       ! If tabulating then fall back to a numerical evaluation.
       potential=self%potentialNumerical(coordinates,status)
    end if
    return
  end function sphericalTabulatedPotential

  double precision function sphericalTabulatedPotentialDifference(self,coordinates1,coordinates2,status) result(potentialDifference)
    !!{RST
    Compute the potential difference between the given ``coordinates1`` and ``coordinates2`` in a spherical mass distribution using a tabulation.
    !!}
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout) , target   :: self
    class           (coordinate                        ), intent(in   )            :: coordinates1, coordinates2
    type            (enumerationStructureErrorCodeType ), intent(  out) , optional :: status
    double precision                                                               :: potential1  , potential2

    if (.not.tabulating) then
       ! Initialize to no potential difference.
       potentialDifference=+0.0d0
       ! Evaluate the potential at both coordinates.
       potential1=self%potential(coordinates1,status)
       if (present(status) .and. status /= structureErrorCodeSuccess) return
       potential2=self%potential(coordinates2,status)
       if (present(status) .and. status /= structureErrorCodeSuccess) return
       ! Re-evaluate potential1 now that we are sure that the potential is tabulated over a sufficient range of radii. This avoids
       ! any zero-point offset changes due to retabulation.
       potential1=self%potential(coordinates1,status)
       if (present(status) .and. status /= structureErrorCodeSuccess) return
       potentialDifference=+potential1 &
            &              -potential2
    else
       ! If tabulating then fall back to a numerical evaluation.
       potentialDifference=+self%potentialDifferenceNumerical(coordinates1,coordinates2,status)
    end if
    return
  end function sphericalTabulatedPotentialDifference

  double precision function sphericalTabulatedVelocityDispersion1D(self,coordinates) result(velocityDispersion)
    !!{RST
    Compute the 1D velocity dispersion at the given ``coordinates`` in a spherical mass distribution using a tabulation.
    !!}
    use :: Galactic_Structure_Options, only : structureErrorCodeSuccess
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout) , target      :: self
    class           (coordinate                        ), intent(in   )               :: coordinates
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               radiusScaled

    if (.not.tabulating) then
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       radiusScaled=+coordinates%rSpherical         () &
            &       /            radiusNormalization
       call self%tabulate(radiusScaled,parameters,container,container%velocityDispersion1D)
       ! Perform the interpolation.
       velocityDispersion=+self%interpolate         (radiusScaled,parameters,container%velocityDispersion1D) &
            &             *     densityNormalization**0.5d0                                                  &
            &             *     radiusNormalization
    else
       ! If tabulating then fall back to a numerical evaluation.
       velocityDispersion=+self%kinematicsDistribution_%velocityDispersion1D(coordinates,self,self)
    end if
    return
  end function sphericalTabulatedVelocityDispersion1D

  double precision function sphericalTabulatedEnergy(self,radiusOuter,massDistributionEmbedding) result(energy)
    !!{RST
    Compute the energy within a given ``radius`` in a spherical mass distribution using a tabulation.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout) , target      :: self
    double precision                                    , intent(in   )               :: radiusOuter
    class           (massDistributionClass             ), intent(inout) , target      :: massDistributionEmbedding
    double precision                                    , dimension(:  ), allocatable :: parameters
    class           (massDistributionClass             )                , pointer     :: self_
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization     , radiusNormalization, &
         &                                                                               radiusOuterScaled

    self_ => self
    if (.not.tabulating.and.associated(self_,massDistributionEmbedding)) then
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       radiusOuterScaled=+radiusOuter         &
            &            /radiusNormalization
       call self%tabulate(radiusOuterScaled,parameters,container,container%energy)
       ! Perform the interpolation.
       energy =+self%interpolate         (radiusOuterScaled,parameters,container%energy) &
            &  *     densityNormalization**2                                             &
            &  *     radiusNormalization **5
    else
       ! If embedded in some other mass distribution, we must fall back to a fully numerical calculation.
       energy=self%energyNumerical(radiusOuter,massDistributionEmbedding)
    end if
    return
  end function sphericalTabulatedEnergy

  double precision function sphericalTabulatedFourierTransform(self,radiusOuter,wavenumber) result(fourierTransform)
    !!{RST
    Compute the Fourier transform of the density profile at the given ``wavenumber`` in a spherical mass distribution using a tabulation.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout)               :: self
    double precision                                    , intent(in   )               :: radiusOuter         , wavenumber
    double precision                                    , dimension(:  ), allocatable :: parameters          , parametersExtended
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               radiusOuterScaled   , wavenumberScaled

    if (.not.tabulating) then
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       wavenumberScaled =+wavenumber          &
            &            *radiusNormalization
       radiusOuterScaled=+radiusOuter         &
            &            /radiusNormalization
       ! Make the scaled outer radius an extra parameter.
       allocate(parametersExtended(size(parameters)+1))
       parametersExtended(1:size(parameters)  )=parameters
       parametersExtended(  size(parameters)+1)=radiusOuterScaled
       call self%tabulate(wavenumberScaled,parametersExtended,container,container%fourierTransform)
       ! Perform the interpolation.
       fourierTransform=+self%interpolate(wavenumberScaled,parametersExtended,container%fourierTransform)
    else
       ! If already tabulating fall back to a numerical calculation.
       fourierTransform=self%fourierTransformNumerical(radiusOuter,wavenumber)
    end if
    return
  end function sphericalTabulatedFourierTransform

  double precision function sphericalTabulatedRadiusFreefall(self,time) result(radius)
    !!{RST
    Compute the freefall radius at a given ``time`` in a spherical mass distribution using a tabulation.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout)               :: self
    double precision                                    , intent(in   )               :: time
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               timeScaled          , timeNormalization

    if (tabulating) then
       radius=self%radiusFreefallNumerical(time)
    else
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       timeNormalization=+1.0d0                      &
            &            /sqrt(densityNormalization)
       timeScaled      =+time              &
            &           /timeNormalization
       call self%tabulate(timeScaled,parameters,container,container%radiusFreefall)
       ! Perform the interpolation.
       radius =+self%interpolate        (timeScaled,parameters,container%radiusFreefall) &
            &  *     radiusNormalization
    end if
    return
  end function sphericalTabulatedRadiusFreefall

  double precision function sphericalTabulatedRadiusFreefallIncreaseRate(self,time) result(radiusIncreaseRate)
    !!{RST
    Compute the rate of increase of freefall radius at a given ``time`` in a spherical mass distribution using a tabulation.
    !!}
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout)               :: self
    double precision                                    , intent(in   )               :: time
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               timeScaled          , timeNormalization

    if (tabulating) then
       radiusIncreaseRate=self%radiusFreefallIncreaseRateNumerical(time)
    else
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       timeNormalization=+1.0d0                      &
            &            /sqrt(densityNormalization)
       timeScaled      =+time              &
            &           /timeNormalization
       call self%tabulate(timeScaled,parameters,container,container%radiusFreefallIncreaseRate)
       ! Perform the interpolation.
       radiusIncreaseRate=+self%interpolate        (timeScaled,parameters,container%radiusFreefallIncreaseRate) &
            &             *     radiusNormalization                                                             &
            &             /     timeNormalization
    end if
    return
  end function sphericalTabulatedRadiusFreefallIncreaseRate

  double precision function sphericalTabulatedDensityRadialMoment(self,moment,radiusMinimum,radiusMaximum,isInfinite) result(densityRadialMoment)
    !!{RST
    Computes radial density moments for spherically-symmetric mass distributions using a tabulation.
    !!}
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout)               :: self
    double precision                                    , intent(in   )               :: moment
    double precision                                    , intent(in   ) , optional    :: radiusMinimum       , radiusMaximum
    logical                                             , intent(  out) , optional    :: isInfinite
    double precision                                    , dimension(:  ), allocatable :: parameters
    type            (massDistributionContainer         )                , pointer     :: container
    double precision                                                                  :: densityNormalization, radiusNormalization, &
         &                                                                               radiusMinimumScaled , radiusMaximumScaled
    logical                                                                           :: useTabulation
    integer                                                                           :: moment_

    useTabulation=.not.tabulating.and.present(radiusMaximum) ! Do not tabulate for infinite outer radius or if already tabulating some other property.
    ! Detect cases for which we can use a tabulation - we tabulate only for the most commonly-used moments.
    if (useTabulation) then
       if      (Values_Agree(moment,0.0d0,absTol=1.0d-6)) then
          moment_      =0
       else if (Values_Agree(moment,1.0d0,absTol=1.0d-6)) then
          moment_      =1
       else if (Values_Agree(moment,2.0d0,absTol=1.0d-6)) then
          moment_      =2
       else if (Values_Agree(moment,3.0d0,absTol=1.0d-6)) then
          moment_      =3
       else
          moment_      =-huge(0)
          useTabulation=.false.
       end if
    end if
    if (useTabulation) then
       ! Get the tabulation properties.
       call self%parameters(densityNormalization,radiusNormalization,parameters,container)
       if (present(radiusMinimum) .and. radiusMinimum > 0.0d0) then
          radiusMinimumScaled=+radiusMinimum       &
               &              /radiusNormalization
       end if
       radiusMaximumScaled=+radiusMaximum       &
            &              /radiusNormalization
       select case (moment_)
       case (0)
          call    self%tabulate(radiusMaximumScaled,parameters,container,container%densityRadialMoment0)
          densityRadialMoment   =                    +self%interpolate(radiusMaximumScaled,parameters,container%densityRadialMoment0)
          if (present(radiusMinimum) .and. radiusMinimum > 0.0d0) then
             call self%tabulate(radiusMinimumScaled,parameters,container,container%densityRadialMoment0)
             densityRadialMoment=+densityRadialMoment-self%interpolate(radiusMinimumScaled,parameters,container%densityRadialMoment0)
          end if
       case (1)
          call    self%tabulate(radiusMaximumScaled,parameters,container,container%densityRadialMoment1)
          densityRadialMoment   =                    +self%interpolate(radiusMaximumScaled,parameters,container%densityRadialMoment1)
          if (present(radiusMinimum) .and. radiusMinimum > 0.0d0) then
             call self%tabulate(radiusMinimumScaled,parameters,container,container%densityRadialMoment1)
             densityRadialMoment=+densityRadialMoment-self%interpolate(radiusMinimumScaled,parameters,container%densityRadialMoment1)
          end if
       case (2)
          call    self%tabulate(radiusMaximumScaled,parameters,container,container%densityRadialMoment2)
          densityRadialMoment   =                    +self%interpolate(radiusMaximumScaled,parameters,container%densityRadialMoment2)
          if (present(radiusMinimum) .and. radiusMinimum > 0.0d0) then
             call self%tabulate(radiusMinimumScaled,parameters,container,container%densityRadialMoment2)
             densityRadialMoment=+densityRadialMoment-self%interpolate(radiusMinimumScaled,parameters,container%densityRadialMoment2)
          end if
       case (3)
          call    self%tabulate(radiusMaximumScaled,parameters,container,container%densityRadialMoment3)
          densityRadialMoment   =                    +self%interpolate(radiusMaximumScaled,parameters,container%densityRadialMoment3)
          if (present(radiusMinimum) .and. radiusMinimum > 0.0d0) then
             call self%tabulate(radiusMinimumScaled,parameters,container,container%densityRadialMoment3)
             densityRadialMoment=+densityRadialMoment-self%interpolate(radiusMinimumScaled,parameters,container%densityRadialMoment3)
          end if
       end select
       ! Perform the interpolation.
       densityRadialMoment=+densityRadialMoment                  &
            &              *densityNormalization                 &
            &              *radiusNormalization **(moment+1.0d0)
    else
       ! Fall back to a numerical calculation.
       densityRadialMoment=self%densityRadialMomentNumerical(moment,radiusMinimum,radiusMaximum,isInfinite)
    end if
    return
  end function sphericalTabulatedDensityRadialMoment

  subroutine sphericalTabulatedTabulate(self,radiusScaled,parameters,container,tabulation)
    !!{RST
    Tabulate the mass distribution.
    !!}
    use :: Coordinates                 , only : coordinateSpherical , assignment(=)
    use :: Input_Paths                 , only : inputPath           , pathTypeDataDynamic
    use :: File_Utilities              , only : File_Lock           , File_Path          , File_Unlock   , lockDescriptor       , &
         &                                      File_Exists         , Directory_Make
    use :: Multi_Counters              , only : multiCounter
    use :: Display                     , only : displayIndent       , displayUnindent    , displayMessage, verbosityLevelWorking, &
         &                                      displayCounter      , displayCounterClear
    use :: Numerical_Constants_Prefixes, only : siFormat
    use :: Numerical_Ranges            , only : Range_Lattice_Offset
    implicit none
    class           (massDistributionSphericalTabulated )                    , intent(inout) :: self
    double precision                                                         , intent(in   ) :: radiusScaled
    double precision                                     , dimension(:      ), intent(in   ) :: parameters
    type            (massDistributionContainer          )                    , intent(inout) :: container
    type            (massDistributionTabulation         )                    , intent(inout) :: tabulation
    class           (massDistributionSphericalTabulated ), save              , pointer       :: instance
    type            (kinematicsDistributionCollisionless), save              , pointer       :: instanceKinematicsDistribution
    !$omp threadprivate(instance,instanceKinematicsDistribution)
    double precision                                     , dimension(:      ), allocatable   :: parameters_                   , parametersReduced_    , &
         &                                                                                      valuesRadius
    double precision                                     , dimension(:,:    ), allocatable   :: valuesParameters
    double precision                                     , dimension(:,:,:,:), allocatable   :: tablePrevious
    integer         (c_size_t                           ), dimension(:      ), allocatable   :: iParameters
    type            (rangeLattice                       ), dimension(:      ), allocatable   :: latticeParametersNew
    type            (rangeLattice                       )                                    :: latticeRadiusNew
    integer         (c_size_t                           ), dimension(4      )                :: countsNew                     , countsPrevious        , &
         &                                                                                      offsets                       , iTable
    integer         (c_size_t                           )                                    :: lengthMaximum                 , iRadius               , &
         &                                                                                      iterationCount                , iterationCountTotal   , &
         &                                                                                      i                             , countParameters_      , &
         &                                                                                      iParameter
    double precision                                                                         :: radius_                       , quantity_             , &
         &                                                                                      time_                         , wavenumber_           , &
         &                                                                                      radiusOuter_                  , density_
    logical                                                                                  :: workRemains                   , carryOver             , &
         &                                                                                      carryOverRadiusComplete       , isCarriedParameters
    type            (multiCounter                       )                                    :: counter
    character       (len= 8                             )                                    :: labelLower                    , labelUpper
    character       (len=64                             )                                    :: labelSize

    ! Test if within current tabulation range.
    if (retabulate()) then
       block
         type(varying_string                     ) :: message     , fileName            , &
              &                                       quantityName
         type(coordinateSpherical                ) :: coordinates , coordinatesZeroPoint
         type(lockDescriptor                     ) :: fileLock
         
         ! Generate the file name.
         quantityName=enumerationQuantityDecode(tabulation%quantity)
         fileName    =inputPath(pathTypeDataDynamic)//'massDistributions/'//self%objectType()//'_'//quantityName//'_'//self%suffix()//'.hdf5'
         ! Restore tabulation from file if necessary.
         call Directory_Make(File_Path(fileName))
         if (File_Exists(fileName)) then
            call File_Lock(fileName,fileLock,lockIsShared=.true.)
            call self%fileRead(fileName,quantityName,container,tabulation)
            call File_Unlock(fileLock)
         end if
         if (retabulate()) then
            tabulating=.true.
            call File_Lock(char(fileName),fileLock,lockIsShared=.false.)
            if (File_Exists(fileName)) then
               call self%fileRead(fileName,quantityName,container,tabulation)
            end if
            ! Test if within current tabulation range.
            if (retabulate()) then
               call displayIndent("tabulating "//enumerationQuantityDecode(tabulation%quantity)//" profile for '"//char(self%objectType())//"'",verbosityLevelWorking)
               ! Construct the radius and parameter ranges, pinning each to an absolute lattice so that the points evaluated -
               ! and therefore every value interpolated between them - depend only on which lattice points are spanned, and not
               ! on the sequence of values which happened to be requested. Each request is passed as the target and the range
               ! already tabulated is unioned in through `latticeCurrent`; folding the latter into the target instead - as the
               ! `min`/`max` against the current bounds formerly did - would apply the factor-of-two margin to an already
               ! margined bound and so ratchet the range outward on every retabulation.
               countParameters_=size(tabulation%countParameters,kind=c_size_t)
               if (countParameters_ > 3_c_size_t) call Error_Report('rank not supported'//{introspection:location})
               latticeRadiusNew=sphericalTabulatedLatticeAxis(radiusScaled,tabulation%radiusCountPer,-huge(0.0d0),+huge(0.0d0),tabulation%latticeRadius)
               allocate(latticeParametersNew(countParameters_))
               do i=1,countParameters_
                  latticeParametersNew(i)=sphericalTabulatedLatticeAxis(                                      &
                       &                                                parameters                       (i), &
                       &                                                tabulation%parametersCountPer    (i), &
                       &                                                container %parametersMinimumLimit(i), &
                       &                                                container %parametersMaximumLimit(i), &
                       &                                                tabulation%latticeParameters     (i)  &
                       &                                               )
               end do
               ! Record where the tabulation already in hand sits within the extended one, so that the values it holds can be
               ! carried over. Every offset is found in exact integer arithmetic from the lattice indices, so no abscissa is
               ! compared.
               carryOver=allocated(tabulation%table) .and. tabulation%latticeRadius%isDefined()
               do i=1,countParameters_
                  carryOver=carryOver .and. tabulation%latticeParameters(i)%isDefined()
               end do
               offsets                        =0_c_size_t
               countsPrevious                 =1_c_size_t
               countsNew                      =1_c_size_t
               countsNew(1                   )=latticeRadiusNew                        %count
               countsNew(2:1+countParameters_)=latticeParametersNew(1:countParameters_)%count
               if (carryOver) then
                  offsets       (1)     =int(Range_Lattice_Offset(tabulation%latticeRadius       ,latticeRadiusNew       ),kind=c_size_t)
                  countsPrevious(1)     =tabulation%latticeRadius%count
                  do i=1,countParameters_
                     offsets       (i+1)=int(Range_Lattice_Offset(tabulation%latticeParameters(i),latticeParametersNew(i)),kind=c_size_t)
                     countsPrevious(i+1)=tabulation%latticeParameters(i)%count
                  end do
                  call Move_Alloc(tabulation%table,tablePrevious)
               end if
               ! Adopt the new lattices, and recover from them every quantity which describes the extent of the tabulation, so
               ! that a tabulation reached by extension cannot come to be described differently from one built in a single pass.
               tabulation%latticeRadius    =latticeRadiusNew
               tabulation%latticeParameters=latticeParametersNew
               tabulation%radiusMinimum    =latticeRadiusNew%minimum()
               tabulation%radiusMaximum    =latticeRadiusNew%maximum()
               tabulation%countRadii       =countsNew(1)
               do i=1,countParameters_
                  tabulation%parametersMinimum(i)=latticeParametersNew(i)%minimum()
                  tabulation%parametersMaximum(i)=latticeParametersNew(i)%maximum()
                  tabulation%countParameters  (i)=countsNew(i+1)
               end do
               call sphericalTabulatedLimitsFlag(container,tabulation)
               ! Take the abscissae from the lattices. They must come from there, and never from an exponentiation open-coded
               ! here: the lattice evaluates them through a single, deliberately un-inlined path, so that a given lattice point
               ! is bit-identical between one tabulation and another regardless of how many points each spans.
               valuesRadius=latticeRadiusNew%values()
               allocate(valuesParameters(maxval(countsNew(2:4)),max(countParameters_,1_c_size_t)))
               valuesParameters=0.0d0
               do i=1,countParameters_
                  valuesParameters(1:countsNew(i+1),i)=latticeParametersNew(i)%values()
               end do
               ! Allocate the table, and carry over the block of values already computed. The trailing dimensions of the table
               ! are of extent 1 where they are unused, in both the previous table and the new one, so the assignment below is
               ! correct at every rank without a rank-by-rank case.
               allocate(tabulation%table(countsNew(1),countsNew(2),countsNew(3),countsNew(4)))
               tabulation%table=0.0d0
               if (carryOver) then
                  tabulation%table(                                           &
                       &           offsets(1)+1:offsets(1)+countsPrevious(1), &
                       &           offsets(2)+1:offsets(2)+countsPrevious(2), &
                       &           offsets(3)+1:offsets(3)+countsPrevious(3), &
                       &           offsets(4)+1:offsets(4)+countsPrevious(4)  &
                       &          )=tablePrevious
                  deallocate(tablePrevious)
               end if
               carryOverRadiusComplete= carryOver                               &
                    &                  .and. offsets       (1) == 0_c_size_t    &
                    &                  .and. countsPrevious(1) == countsNew (1)
               ! Report on tabulation.
               lengthMaximum=max(12,maxval(len(container%nameParameters)))
               write (labelLower,'(e8.2)') tabulation%radiusMinimum
               write (labelUpper,'(e8.2)') tabulation%radiusMaximum
               message=labelLower//" ≤ radiusScaled"//repeat(" ",lengthMaximum-12)//" ≤ "//labelUpper
               call displayMessage(message,verbosityLevelWorking)
               do i=1,size(tabulation%countParameters)
                  write (labelLower,'(e8.2)') tabulation%parametersMinimum(i)
                  write (labelUpper,'(e8.2)') tabulation%parametersMaximum(i)
                  message=labelLower//" ≤ "//container%nameParameter(i,tabulation)//repeat(" ",lengthMaximum-len(container%nameParameter(i,tabulation)))//" ≤ "//labelUpper
                  call displayMessage(message,verbosityLevelWorking)
               end do
               labelSize=siFormat(dble(sizeof(tabulation%table)),'f8.2,1x')
               message="tabulation size = "//trim(adjustl(labelSize))//"B"
               call displayMessage(message,verbosityLevelWorking)
               ! Iterate over parameters. Where the whole radius axis is carried over, a parameter state which also lies within
               ! the carried block does no work at all, and so is not counted here.
               iterationCount     =0_c_size_t
               if (carryOverRadiusComplete) then
                  iterationCountTotal=product(countsNew(2:4))-product(countsPrevious(2:4))
               else
                  iterationCountTotal=product(countsNew(2:4))
               end if
               iterationCountTotal=max(iterationCountTotal,1_c_size_t)
               ! Tabulate in parallel.
               !$omp parallel private(iRadius,radius_,time_,radiusOuter_,wavenumber_,density_,coordinates,coordinatesZeroPoint,quantity_,iParameters,iParameter,iTable,parameters_,parametersReduced_,workRemains,counter,isCarriedParameters)
               ! This is a new thread, so mark it as tabulating.
               tabulating=.true.
               ! Initialize the counter and iterate over parameter states. Note that the counter is private to each thread but
               ! is incremented in lockstep by all of them, so that the decision to skip a parameter state below is uniform
               ! across the team - as it must be, since the loop over radii which follows it is a worksharing construct.
               counter=multiCounter(tabulation%countParameters)
               do while (.true.)
                  !$omp barrier
                  workRemains=counter%increment()
                  if (.not.workRemains) exit
                  iParameters   =counter%states()
                  ! Determine whether this parameter state lies within the block of values carried over from the tabulation
                  ! already in hand, and skip it entirely if it does and the radius axis was carried over complete.
                  isCarriedParameters=carryOver
                  do iParameter=1,countParameters_
                     isCarriedParameters=       isCarriedParameters                                                           &
                          &               .and. iParameters(iParameter) >  offsets(iParameter+1)                              &
                          &               .and. iParameters(iParameter) <= offsets(iParameter+1)+countsPrevious(iParameter+1)
                  end do
                  if (isCarriedParameters .and. carryOverRadiusComplete) cycle
                  !$omp masked
                  call displayCounter(int(100.0d0*dble(iterationCount)/dble(iterationCountTotal)),iterationCount==0,verbosityLevelWorking)
                  iterationCount=iterationCount+1_c_size_t
                  !$omp end masked
                  if (.not.allocated(parameters_)) allocate(parameters_(countParameters_))
                  do iParameter=1,countParameters_
                     parameters_(iParameter)=valuesParameters(iParameters(iParameter),iParameter)
                  end do
                  ! Record the position of this parameter state within the table. The trailing indices are unity where the
                  ! corresponding dimension is unused.
                  iTable        =1_c_size_t
                  iTable(2:1+countParameters_)=iParameters(1:countParameters_)
                  ! Call the factory function in the child class to get an instance built with the current parameters.
                  select case (tabulation%quantity%ID)
                  case (quantityFourierTransform%ID)
                     if (.not.allocated(parametersReduced_)) allocate(parametersReduced_(size(parameters_-1)))
                     parametersReduced_ =  parameters_(1:size(parameters_)-1)
                     radiusOuter_       =  parameters_(  size(parameters_)  )
                     instance           => self%factoryTabulation(parametersReduced_)
                  case default
                     instance           => self%factoryTabulation(parameters_       )
                  end select
                  allocate(instanceKinematicsDistribution)
                  !![
		  <referenceConstruct object="instanceKinematicsDistribution">
		    <constructor>
		      kinematicsDistributionCollisionless(                                                                                                                    &amp;
                       &amp;                              toleranceRelativeVelocityDispersion       =self%kinematicsDistribution_%toleranceRelativeVelocityDispersion       , &amp;
                       &amp;                              toleranceRelativeVelocityDispersionMaximum=self%kinematicsDistribution_%toleranceRelativeVelocityDispersionMaximum  &amp;
                       &amp;                             )
		    </constructor>
		  </referenceConstruct>
                  !!]
                  call instance%setKinematicsDistribution(instanceKinematicsDistribution)
                  !![
		  <objectDestructor name="instanceKinematicsDistribution"/>
                  !!]
                  ! Iterate over scaled radii.
                  !$omp do schedule(dynamic)
                  do iRadius=1,tabulation%countRadii
                     ! Skip the values carried over from the tabulation already in hand - evaluating them again would merely
                     ! reproduce them, at the cost of a numerical solution apiece.
                     if     (                                              &
                          &        isCarriedParameters                     &
                          &  .and. iRadius >  offsets(1)                   &
                          &  .and. iRadius <= offsets(1)+countsPrevious(1) &
                          & ) cycle
                     radius_             =valuesRadius(iRadius)
                     time_               = radius_
                     wavenumber_         = radius_
                     density_            = radius_
                     coordinates         =[radius_,0.0d0,0.0d0]
                     coordinatesZeroPoint=[1.0d0  ,0.0d0,0.0d0]
                     ! Compute the quantity numerically.
                     select case (tabulation%quantity%ID)
                     case (quantityMass                      %ID)
                        quantity_=+instance%massEnclosedBySphereNumerical      (             radius_                      )
                     case (quantityRadiusEnclosingDensity    %ID)
                        quantity_=+instance%radiusEnclosingDensityNumerical    (             density_                     )
                     case (quantityEnergy                    %ID)
                        quantity_=+instance%energyNumerical                    (             radius_             ,instance)
                     case (quantityPotential                 %ID)
                        ! Potential is always referenced to a zero-point at the normalization radius.
                        quantity_=+instance%potentialNumerical                 (             coordinates                  ) &
                             &    -instance%potentialNumerical                 (             coordinatesZeroPoint         )
                     case (quantityVelocityDispersion1D      %ID)
                        quantity_=+instance%velocityDispersion1D               (             coordinates                  )
                     case (quantityFourierTransform          %ID)
                        quantity_=+instance%fourierTransformNumerical          (radiusOuter_,wavenumber_                  )
                     case (quantityRadiusFreefall            %ID)
                        quantity_=+instance%radiusFreefallNumerical            (             time_                        )
                     case (quantityRadiusFreefallIncreaseRate%ID)
                        quantity_=+instance%radiusFreefallIncreaseRateNumerical(             time_                        )
                     case (quantityDensityRadialMoment0      %ID)
                        quantity_=+instance% densityRadialMomentNumerical      (0.0d0,0.0d0 ,radius_                      )
                     case (quantityDensityRadialMoment1      %ID)
                        quantity_=+instance% densityRadialMomentNumerical      (1.0d0,0.0d0 ,radius_                      )
                     case (quantityDensityRadialMoment2      %ID)
                        quantity_=+instance% densityRadialMomentNumerical      (2.0d0,0.0d0 ,radius_                      )
                     case (quantityDensityRadialMoment3      %ID)
                        quantity_=+instance% densityRadialMomentNumerical      (3.0d0,0.0d0 ,radius_                      )
                     case default
                        quantity_=+0.0d0
                        call Error_Report('unknown quantity'//{introspection:location})
                     end select
                     ! For quantities that are negative, invert the sign prior to any logarithmic transform.
                     if (tabulation%isNegative  ) quantity_=-    quantity_
                     ! If logarithmic interpolation is requested, log-transform now.
                     if (tabulation%logTransform) quantity_=+log(quantity_)
                     ! Store the quantity.
                     iTable(1)=iRadius
                     tabulation%table(iTable(1),iTable(2),iTable(3),iTable(4))=quantity_
                  end do
                  !$omp end do
                  deallocate(instance)
                  nullify   (instance)
               end do
               !$omp masked
               call displayCounterClear(verbosityLevelWorking)
               !$omp end masked
               tabulating=.false.
               !$omp end parallel
               ! Store tabulation to file.
               call self%fileWrite(fileName,quantityName,container,tabulation)
               call displayUnindent('done',verbosityLevelWorking)
            end if
            tabulating=.false.
            call File_Unlock(fileLock)
         end if
       end block
    end if
    return

  contains

    logical function retabulate()
      !!{RST
      Test if the mass profile must be retabulated.

      A request which lies beyond a parameter bound calls for a retabulation only if that bound can actually be moved. A bound
      which has reached its hard limit cannot: the limit is imposed by snapping the bound *inward* to a lattice point, so the
      bound stops marginally short of the limit and a request in between would otherwise ask for a retabulation which
      reproduced exactly the same range - on every evaluation, each taking the file lock and rewriting the stored tabulation.
      Such requests are met by the flat extrapolation which the interpolation performs beyond the ends of the table.
      !!}
      implicit none

      retabulate=.not.allocated(tabulation%table)
      if (.not.retabulate)                                                                                                &
           &  retabulate=     radiusScaled < tabulation%radiusMinimum                                                     &
           &             .or.                                                                                             &
           &                  radiusScaled > tabulation%radiusMaximum                                                     &
           &             .or.                                                                                             &
           &              any(parameters   < tabulation%parametersMinimum .and. .not.tabulation%parametersAtLimitMinimum) &
           &             .or.                                                                                             &
           &              any(parameters   > tabulation%parametersMaximum .and. .not.tabulation%parametersAtLimitMaximum)
      return
    end function retabulate

  end subroutine sphericalTabulatedTabulate

  double precision function sphericalTabulatedInterpolate(self,radiusScaled,parameters,tabulation) result(interpolated)
    !!{RST
    Interpolate in the tabulated mass distribution.
    !!}
    use :: Multi_Counters, only : multiCounter
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout)                 :: self
    double precision                                    , intent(in   )                 :: radiusScaled
    double precision                                    , intent(in   ), dimension(:  ) :: parameters
    type            (massDistributionTabulation        ), intent(in   )                 :: tabulation
    double precision                                                   , dimension(2  ) :: hRadius
    integer         (c_size_t                          ), allocatable  , dimension(:  ) :: iParameters , jParameters
    double precision                                    , allocatable  , dimension(:,:) :: hParameters
    integer         (c_size_t                          )                                :: iRadius     , jRadius
    type            (multiCounter                      )                                :: counter

    ! Compute interpolating factors. Along each axis the position of the value is found as its coordinate on the absolute
    ! lattice - a quantity which depends only on the value and on the density of lattice points, never on which part of the
    ! lattice this particular tabulation happens to span - and is split there into the index of the lattice point below it and
    ! the fractional position within the interval. Only then is the index of the first tabulated point subtracted, in exact
    ! integer arithmetic.
    !
    ! Doing it in that order is what makes the interpolation invariant under extension of the tabulation. Forming the position
    ! relative to the first tabulated point *before* taking its fractional part would not be: the subtraction is exact, but the
    ! fractional part is then extracted from a number whose magnitude depends on where the range begins, and so is rounded to a
    ! coarser grid the further from the lattice origin that is. Extending the range downward would then perturb every
    ! interpolated value in the last few bits, which is precisely the sequence dependence this tabulation exists to remove.
    !
    ! Beyond the ends of the table - which can happen for a parameter whose bound has reached its hard limit - the index is
    ! confined to the interior range [1,count-1], so that the pair of nodes spanning the interval always exists, and the weight
    ! is set to the nearer end, giving fixed extrapolation from the edge value. That is continuous across the table boundary,
    ! which matters for iterative solvers that evaluate the distribution as they converge.
    allocate(hParameters(2,size(parameters)))
    allocate(iParameters(  size(parameters)))
    hRadius    (2  )=log(radiusScaled)/log(2.0d0)*dble(tabulation%radiusCountPer    )
    hParameters(2,:)=log(parameters  )/log(2.0d0)*dble(tabulation%parametersCountPer)
    iRadius         =floor(hRadius    (2  ),kind=c_size_t)
    iParameters     =floor(hParameters(2,:),kind=c_size_t)
    hRadius    (2  )=hRadius    (2  )-dble(iRadius    )
    hParameters(2,:)=hParameters(2,:)-dble(iParameters)
    iRadius         =iRadius    -int(tabulation%latticeRadius    %indexMinimum,kind=c_size_t)+1_c_size_t
    iParameters     =iParameters-int(tabulation%latticeParameters%indexMinimum,kind=c_size_t)+1_c_size_t
    if (iRadius <  1_c_size_t                          ) hRadius(2  )=0.0d0
    if (iRadius > tabulation%countRadii     -1_c_size_t) hRadius(2  )=1.0d0
    where (iParameters <  1_c_size_t                          ) hParameters(2,:)=0.0d0
    where (iParameters > tabulation%countParameters-1_c_size_t) hParameters(2,:)=1.0d0
    iRadius         =min(max(iRadius    ,1_c_size_t),tabulation%countRadii     -1_c_size_t)
    iParameters     =min(max(iParameters,1_c_size_t),tabulation%countParameters-1_c_size_t)
    hRadius    (1  )=1.0d0-hRadius    (2  )
    hParameters(1,:)=1.0d0-hParameters(2,:)
    ! Perform the interpolation.
    interpolated=0.0d0
    counter     =multiCounter(spread(2_c_size_t,1,size(parameters)))
    select case(size(parameters))
    case (0)
          do jRadius=1,2
             interpolated=+interpolated                                                 &
                  &       +tabulation %table(                                           &
                  &                          iRadius       +jRadius       -1_c_size_t,  &
                  &                                                        1_c_size_t,  &
                  &                                                        1_c_size_t,  &
                  &                                                        1_c_size_t   &
                  &                         )                                           &
                  &       *hRadius          (                                           &
                  &                                         jRadius                     &
                  &                         )
          end do
    case (1)
       do while (counter%increment())
          jParameters=counter%states()
          do jRadius=1,2
             interpolated=+interpolated                                                 &
                  &       +tabulation %table(                                           &
                  &                          iRadius       +jRadius       -1_c_size_t,  &
                  &                          iParameters(1)+jParameters(1)-1_c_size_t,  &
                  &                                                        1_c_size_t,  &
                  &                                                        1_c_size_t   &
                  &                         )                                           &
                  &       *hRadius          (                                           &
                  &                                         jRadius                     &
                  &                         )                                           &
                  &       *hParameters      (                                           &
                  &                                         jParameters(1)           ,1 &
                  &                         )
          end do
       end do
    case (2)
       do while (counter%increment())
          jParameters=counter%states()
          do jRadius=1,2
             interpolated=+interpolated                                                 &
                  &       +tabulation %table(                                           &
                  &                          iRadius       +jRadius       -1_c_size_t,  &
                  &                          iParameters(1)+jParameters(1)-1_c_size_t,  &
                  &                          iParameters(2)+jParameters(2)-1_c_size_t,  &
                  &                                                        1_c_size_t   &
                  &                         )                                           &
                  &       *hRadius          (                                           &
                  &                                         jRadius                     &
                  &                         )                                           &
                  &       *hParameters      (                                           &
                  &                                         jParameters(1)           ,1 &
                  &                         )                                           &
                  &       *hParameters      (                                           &
                  &                                         jParameters(2)           ,2 &
                  &                         )
          end do
       end do
    case (3)
       do while (counter%increment())
          jParameters=counter%states()
          do jRadius=1,2
             interpolated=+interpolated                                                 &
                  &       +tabulation %table(                                           &
                  &                          iRadius       +jRadius       -1_c_size_t,  &
                  &                          iParameters(1)+jParameters(1)-1_c_size_t,  &
                  &                          iParameters(2)+jParameters(2)-1_c_size_t,  &
                  &                          iParameters(3)+jParameters(3)-1_c_size_t   &
                  &                         )                                           &
                  &       *hRadius          (                                           &
                  &                                         jRadius                     &
                  &                         )                                           &
                  &       *hParameters      (                                           &
                  &                                         jParameters(1)           ,1 &
                  &                         )                                           &
                  &       *hParameters      (                                           &
                  &                                         jParameters(2)           ,2 &
                  &                         )                                           &
                  &       *hParameters      (                                           &
                  &                                         jParameters(3)           ,3 &
                  &                         )
          end do
       end do
    case default
       call Error_Report('rank not supported'//{introspection:location})
    end select
    ! If logarithmic interpolation was used, inverse transform now.
    if (tabulation%logTransform) interpolated=+exp(interpolated)
    ! For quantities that are negative, invert the sign.
    if (tabulation%isNegative  ) interpolated=-    interpolated
    return
  end function sphericalTabulatedInterpolate
  
  logical function sphericalTabulatedIsTabulating(self) result(isTabulating)
    !!{RST
    Return true if this thread is currently tabulating.
    !!}
    implicit none
    class(massDistributionSphericalTabulated), intent(inout) :: self
    !$GLC attributes unused :: self
    
    isTabulating=tabulating
    return
  end function sphericalTabulatedIsTabulating
  
  subroutine sphericalTabulatedFileRead(self,fileName,quantityName,container,tabulation)
    !!{RST
    Read tabulated data from file.

    The stored tabulation is adopted only if the file records, for every axis, a lattice which is self-consistent and which
    uses the density of points that this object would use, and if the array stored alongside those lattices has the extent
    they imply. A file written before the lattices were recorded, or with a different grid density, therefore leaves the
    tabulation already in hand untouched rather than being misread.
    !!}
    use :: HDF5_Access     , only : hdf5Access
    use :: IO_HDF5         , only : hdf5File
    use :: String_Handling , only : String_Upper_Case_First
    use :: Display         , only : displayMessage          , verbosityLevelWorking
    use :: Numerical_Ranges, only : gridSchemePerOctave
    use :: Table_Caches    , only : Table_Cache_Lattice_Read
    implicit none
    class           (massDistributionSphericalTabulated), intent(inout)                     :: self
    type            (varying_string                    ), intent(in   )                     :: fileName         , quantityName
    type            (massDistributionContainer         ), intent(inout)                     :: container
    type            (massDistributionTabulation        ), intent(inout)                     :: tabulation
    type            (hdf5File                          )                                    :: file
    type            (rangeLattice                      )                                    :: latticeRadius
    type            (rangeLattice                      ), allocatable  , dimension(:      ) :: latticeParameters
    double precision                                    , allocatable  , dimension(:,:,:,:) :: table
    integer         (c_size_t                          )                                    :: i                , countParameters_
    logical                                                                                 :: isUsable

    call displayMessage("reading tabulated "//enumerationQuantityDecode(tabulation%quantity)//" profile from '"//char(fileName)//"'",verbosityLevelWorking)
    countParameters_=container%countParameters(tabulation)
    if (.not.allocated(tabulation%parametersMinimum       )) allocate(tabulation%parametersMinimum       (countParameters_))
    if (.not.allocated(tabulation%parametersMaximum       )) allocate(tabulation%parametersMaximum       (countParameters_))
    if (.not.allocated(tabulation%parametersAtLimitMinimum)) allocate(tabulation%parametersAtLimitMinimum(countParameters_))
    if (.not.allocated(tabulation%parametersAtLimitMaximum)) allocate(tabulation%parametersAtLimitMaximum(countParameters_))
    if (.not.allocated(tabulation%countParameters         )) allocate(tabulation%countParameters         (countParameters_))
    if (.not.allocated(tabulation%latticeParameters       )) allocate(tabulation%latticeParameters       (countParameters_))
    allocate(latticeParameters(countParameters_))
    !$ call hdf5Access%set()
    file=hdf5File(fileName,readOnly=.true.)
    ! Recover the lattices on which the stored tabulation was built.
    call    Table_Cache_Lattice_Read(file,char(quantityName)//'Radius'                                                            ,gridSchemePerOctave,int(tabulation%radiusCountPer       ),latticeRadius       )
    isUsable=latticeRadius%isDefined()
    do i=1,countParameters_
       call Table_Cache_Lattice_Read(file,char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation))),gridSchemePerOctave,int(tabulation%parametersCountPer(i)),latticeParameters(i))
       isUsable=isUsable .and. latticeParameters(i)%isDefined()
    end do
    if (isUsable) call file%readDataset(char(quantityName),table)
    !$ call hdf5Access%unset()
    ! Check that the array stored alongside the lattices has the extent which they imply - including that the dimensions which
    ! are unused, this tabulation having fewer than three parameters, are of extent unity.
    if (isUsable) then
       isUsable=size(table,dim=1,kind=c_size_t) == int(latticeRadius%count,kind=c_size_t)
       do i=1,countParameters_
          isUsable=isUsable .and. size(table,dim=int(i)+1,kind=c_size_t) == int(latticeParameters(i)%count,kind=c_size_t)
       end do
       do i=countParameters_+1,3_c_size_t
          isUsable=isUsable .and. size(table,dim=int(i)+1,kind=c_size_t) == 1_c_size_t
       end do
    end if
    ! Decline a stored tabulation which does not contain the one already in hand: since this is called only when the latter has
    ! been found insufficient, adopting a narrower tabulation in its place would discard values which must then be computed
    ! again. The range in hand is extended in its place, and the file rewritten from it.
    if (isUsable .and. allocated(tabulation%table)) then
       if (tabulation%latticeRadius%isDefined()) then
          isUsable=latticeRadius%covers(tabulation%latticeRadius)
          do i=1,countParameters_
             isUsable=isUsable .and. latticeParameters(i)%covers(tabulation%latticeParameters(i))
          end do
       end if
    end if
    if (.not.isUsable) return
    ! Adopt the stored tabulation, recovering every quantity which describes its extent from the lattices rather than reading
    ! them from the file, so that a restored tabulation cannot come to be described differently from a freshly built one.
    if (allocated(tabulation%table)) deallocate(tabulation%table)
    call Move_Alloc(table,tabulation%table)
    tabulation%latticeRadius    =latticeRadius
    tabulation%latticeParameters=latticeParameters
    tabulation%radiusMinimum    =latticeRadius%minimum()
    tabulation%radiusMaximum    =latticeRadius%maximum()
    tabulation%countRadii       =latticeRadius%count
    do i=1,countParameters_
       tabulation%parametersMinimum(i)=latticeParameters(i)%minimum()
       tabulation%parametersMaximum(i)=latticeParameters(i)%maximum()
       tabulation%countParameters  (i)=latticeParameters(i)%count
    end do
    call sphericalTabulatedLimitsFlag(container,tabulation)
    return
  end subroutine sphericalTabulatedFileRead

  subroutine sphericalTabulatedFileWrite(self,fileName,quantityName,container,tabulation)
    !!{RST
    Write tabulated data to file.
    !!}
    use :: HDF5_Access    , only : hdf5Access
    use :: IO_HDF5        , only : hdf5File
    use :: String_Handling, only : String_Upper_Case_First
    use :: Display        , only : displayMessage           , verbosityLevelWorking
    use :: Table_Caches   , only : Table_Cache_Lattice_Write
    implicit none
    class  (massDistributionSphericalTabulated), intent(inout) :: self
    type   (varying_string                    ), intent(in   ) :: fileName  , quantityName
    type   (massDistributionContainer         ), intent(inout) :: container
    type   (massDistributionTabulation        ), intent(in   ) :: tabulation
    type   (hdf5File                          )                :: file
    integer(c_size_t                          )                :: i

    call displayMessage("writing tabulated "//char(quantityName)//" profile to '"//char(fileName)//"'",verbosityLevelWorking)
    !$ call hdf5Access%set()
    file=hdf5File(fileName,overWrite=.true.)
    ! Record the lattices on which the axes are built. The limits and inverse steps formerly stored alongside them are not:
    ! each is a function of the lattices, and is recovered from them when the file is read.
    call    Table_Cache_Lattice_Write(file,char(quantityName)//'Radius'                                                            ,tabulation%latticeRadius       )
    do i=1,container%countParameters(tabulation)
       call Table_Cache_Lattice_Write(file,char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation))),tabulation%latticeParameters(i))
    end do
    call file%writeDataset(tabulation%table,char(quantityName),'Tabulated '//char(quantityName)//' profile.')
    !$ call hdf5Access%unset()
    return
  end subroutine sphericalTabulatedFileWrite

  function sphericalTabulatedLatticeAxis(valueTarget,countPer,limitMinimum,limitMaximum,latticeCurrent) result(lattice)
    !!{RST
    Return the pinned ``rangeLattice`` covering ``valueTarget`` for a single axis of a tabulation, taking the union with the
    lattice ``latticeCurrent`` on which that axis is currently built.

    A hard lower limit is imposed only where it is positive: a logarithmic lattice cannot express a limit of zero or below, and
    such a limit is in any case inert here since the request is itself always positive and the safety margin is a factor. Note
    that where the limit *is* imposed the request is first clamped up to it, rather than the limit being added to the set of
    target values: adding it would make it the smallest target and so would drag the lower bound of every tabulation down to
    the limit instead of leaving it where the request puts it.
    !!}
    use :: Numerical_Ranges, only : Range_Pinned, gridSchemePerOctave
    implicit none
    type            (rangeLattice)                :: lattice
    double precision              , intent(in   ) :: valueTarget   , limitMinimum, &
         &                                           limitMaximum
    integer         (c_size_t    ), intent(in   ) :: countPer
    type            (rangeLattice), intent(in   ) :: latticeCurrent
    double precision                              :: valueTarget_
    integer                                       :: pointsPer

    pointsPer   =int(countPer)
    valueTarget_=valueTarget
    ! The bounds are pinned to whole octaves - the default anchor interval - which is exactly the factor of two by which the
    ! request is bracketed, so no tabulation is made wider than the margin already made it.
    if (limitMinimum > 0.0d0) then
       valueTarget_=max(valueTarget_,limitMinimum)
       lattice=Range_Pinned(                                     &
            &                              [valueTarget_ ]     , &
            &                               pointsPer          , &
            &                               gridSchemePerOctave, &
            &               marginFactor  = 2.0d0              , &
            &               limitMinimum  = limitMinimum       , &
            &               limitMaximum  = limitMaximum       , &
            &               latticeCurrent= latticeCurrent       &
            &              )
    else
       lattice=Range_Pinned(                                     &
            &                              [valueTarget_ ]     , &
            &                               pointsPer          , &
            &                               gridSchemePerOctave, &
            &               marginFactor  = 2.0d0              , &
            &               limitMaximum  = limitMaximum       , &
            &               latticeCurrent= latticeCurrent       &
            &              )
    end if
    return
  end function sphericalTabulatedLatticeAxis

  subroutine sphericalTabulatedLimitsFlag(container,tabulation)
    !!{RST
    Record, for each parameter of ``tabulation``, whether the corresponding bound has reached the hard limit imposed on it -
    that is, whether the next lattice point beyond that bound would lie outside the limit. Where it has, no retabulation can
    move that bound, and a request beyond it must be met by extrapolation instead.
    !!}
    implicit none
    type            (massDistributionContainer ), intent(inout) :: container
    type            (massDistributionTabulation), intent(inout) :: tabulation
    integer         (c_size_t                  )                :: i
    double precision                                            :: stepFactor

    do i=1,container%countParameters(tabulation)
       stepFactor                            = 2.0d0**(1.0d0/dble(tabulation%parametersCountPer(i)))
       tabulation%parametersAtLimitMinimum(i)=tabulation%parametersMinimum(i)/stepFactor < container%parametersMinimumLimit(i)
       tabulation%parametersAtLimitMaximum(i)=tabulation%parametersMaximum(i)*stepFactor > container%parametersMaximumLimit(i)
    end do
    return
  end subroutine sphericalTabulatedLimitsFlag

  subroutine massDistributionContainerInitialize(self,countParameters)
    !!{RST
    Initialize the container to the specified number of parameters.
    !!}
    implicit none
    class  (massDistributionContainer), intent(inout) :: self
    integer                           , intent(in   ) :: countParameters

    ! Allocate the internal arrays.
    allocate(self%descriptionParameters (countParameters  ))
    allocate(self%nameParameters        (countParameters  ))
    allocate(self%parametersMinimumLimit(countParameters+1))
    allocate(self%parametersMaximumLimit(countParameters+1))
    ! Initialize each tabulation. The Fourier transform tabulation carries the scaled outer radius as an extra parameter.
    call massDistributionTabulationInitialize(self%mass                      ,countParameters  )
    call massDistributionTabulationInitialize(self%radiusEnclosingDensity    ,countParameters  )
    call massDistributionTabulationInitialize(self%energy                    ,countParameters  )
    call massDistributionTabulationInitialize(self%potential                 ,countParameters  )
    call massDistributionTabulationInitialize(self%velocityDispersion1D      ,countParameters  )
    call massDistributionTabulationInitialize(self%radiusFreefall            ,countParameters  )
    call massDistributionTabulationInitialize(self%radiusFreefallIncreaseRate,countParameters  )
    call massDistributionTabulationInitialize(self%densityRadialMoment0      ,countParameters  )
    call massDistributionTabulationInitialize(self%densityRadialMoment1      ,countParameters  )
    call massDistributionTabulationInitialize(self%densityRadialMoment2      ,countParameters  )
    call massDistributionTabulationInitialize(self%densityRadialMoment3      ,countParameters  )
    call massDistributionTabulationInitialize(self%fourierTransform          ,countParameters+1)
    ! Initialize limits on the parameters. By default no limit is imposed.
    self                           %parametersMinimumLimit=-huge(0.0d0)
    self                           %parametersMaximumLimit=+huge(0.0d0)
    return
  end subroutine massDistributionContainerInitialize

  subroutine massDistributionTabulationInitialize(tabulation,countParameters)
    !!{RST
    Initialize a single tabulation to the specified number of parameters. The minima and maxima of the tabulation ranges are
    chosen to ensure that the first tabulation will force them to be reset.
    !!}
    implicit none
    type   (massDistributionTabulation), intent(inout) :: tabulation
    integer                            , intent(in   ) :: countParameters

    allocate(tabulation%parametersMinimum       (countParameters))
    allocate(tabulation%parametersMaximum       (countParameters))
    allocate(tabulation%parametersAtLimitMinimum(countParameters))
    allocate(tabulation%parametersAtLimitMaximum(countParameters))
    allocate(tabulation%parametersCountPer      (countParameters))
    allocate(tabulation%countParameters         (countParameters))
    allocate(tabulation%latticeParameters       (countParameters))
    tabulation%radiusMinimum           =+huge(0.0d0)
    tabulation%radiusMaximum           =-huge(0.0d0)
    tabulation%parametersMinimum       =+huge(0.0d0)
    tabulation%parametersMaximum       =-huge(0.0d0)
    tabulation%parametersAtLimitMinimum=.false.
    tabulation%parametersAtLimitMaximum=.false.
    return
  end subroutine massDistributionTabulationInitialize

  function massDistributionContainerNameParameter(self,indexParameter,tabulation) result(nameParameter)
    !!{RST
    Return the name of the indexed parameter.
    !!}
    implicit none
    type   (varying_string            )                :: nameParameter
    class  (massDistributionContainer ), intent(inout) :: self
    integer(c_size_t                  ), intent(in   ) :: indexParameter
    type   (massDistributionTabulation), intent(in   ) :: tabulation
    integer(c_size_t                  )                :: countParameters

    countParameters=self%countParameters(tabulation)
    if     (                                  &
         &   indexParameter < 1               &
         &  .or.                              &
         &   indexParameter > countParameters &
         & ) call Error_Report('`indexParameter` is out of range'//{introspection:location})
    select case (tabulation%quantity%ID)
    case (quantityFourierTransform%ID)
       if (indexParameter == countParameters) then
          nameParameter='radiusOuter'
       else
          nameParameter=self%nameParameters(indexParameter)
       end if
    case default
       nameParameter   =self%nameParameters(indexParameter)
    end select
    return
  end function massDistributionContainerNameParameter
  
  function massDistributionContainerCountParameters(self,tabulation) result(countParameters)
    !!{RST
    Return the name of the indexed parameter.
    !!}
    implicit none
    integer(c_size_t                  )                :: countParameters
    class  (massDistributionContainer ), intent(inout) :: self
    type   (massDistributionTabulation), intent(in   ) :: tabulation

    select case (tabulation%quantity%ID)
    case (quantityFourierTransform%ID)
       countParameters=size(self%nameParameters)+1_c_size_t
    case default
       countParameters=size(self%nameParameters)
    end select
    return
  end function massDistributionContainerCountParameters
