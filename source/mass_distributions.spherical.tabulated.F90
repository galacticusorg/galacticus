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
  Implementation of an abstract mass distribution class for tabulated spherically symmetric distributions.
  !!}

  use, intrinsic :: ISO_C_Binding, only : c_size_t

  !![
  <massDistribution name="massDistributionSphericalTabulated" abstract="yes">
   <description>An abstract mass distribution class for tabulated spherically symmetric distributions.</description>
  </massDistribution>
  !!]
  type, extends(massDistributionSpherical), abstract :: massDistributionSphericalTabulated
     !!{
     Implementation of an abstract mass distribution class for tabulated spherically symmetric distributions.
     !!}
     private
   contains
     !![
     <methods>
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
  <enumeration>
   <name>quantity</name>
   <description>Quantities for tabulation mass distributions.</description>
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
  
  type :: massDistributionTabulation
     !!{
     Object used to store individual mass distribution tabulations.
     !!}
     type            (enumerationQuantityType)                                   :: quantity
     logical                                                                     :: logTransform         , isNegative
     double precision                                                            :: radiusMinimum        , radiusMaximum    , &
          &                                                                         radiusInverseStep
     double precision                          , allocatable, dimension(:      ) :: parametersMinimum    , parametersMaximum, &
          &                                                                         parametersInverseStep
     integer         (c_size_t                )                                  :: radiusCountPer       , countRadii
     integer         (c_size_t                ), allocatable, dimension(:      ) :: parametersCountPer   , countParameters
     double precision                          , allocatable, dimension(:,:,:,:) :: table
  end type massDistributionTabulation

  type :: massDistributionContainer
     !!{
     Object to store collections of mass distribution tabulations.
     !!}
     type            (varying_string            ), allocatable, dimension(:) :: nameParameters        , descriptionParameters
     double precision                            , allocatable, dimension(:) :: parametersMinimumLimit, parametersMaximumLimit
     ! Tabulations for individual quantities.
     type            (massDistributionTabulation)                            :: mass                      =massDistributionTabulation(quantityMass                      ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: radiusEnclosingDensity    =massDistributionTabulation(quantityRadiusEnclosingDensity    ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: potential                 =massDistributionTabulation(quantityPotential                 ,.false.,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: energy                    =massDistributionTabulation(quantityEnergy                    ,.false.,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: fourierTransform          =massDistributionTabulation(quantityFourierTransform          ,.false.,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: radiusFreefall            =massDistributionTabulation(quantityRadiusFreefall            ,.false.,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: radiusFreefallIncreaseRate=massDistributionTabulation(quantityRadiusFreefallIncreaseRate,.false.,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: densityRadialMoment0      =massDistributionTabulation(quantityDensityRadialMoment0      ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: densityRadialMoment1      =massDistributionTabulation(quantityDensityRadialMoment1      ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: densityRadialMoment2      =massDistributionTabulation(quantityDensityRadialMoment2      ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: densityRadialMoment3      =massDistributionTabulation(quantityDensityRadialMoment3      ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
     type            (massDistributionTabulation)                            :: velocityDispersion1D      =massDistributionTabulation(quantityVelocityDispersion1D      ,.true. ,.false.,0.0d0,0.0d0,0.0d0,null(),null(),null(),0_c_size_t,0_c_size_t,null(),null(),null())
   contains
     !![
     <methods>
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
     !!{
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
    !!{
    Computes the mass enclosed within a sphere of given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions using a tabulation.
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
    !!{
    Computes the radius enclosing the given {\normalfont \ttfamily radius} for spherically-symmetric mass
    distributions using a tabulation.
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
    !!{
    Compute the potential at the given {\normalfont \ttfamily coordinates} in a spherical mass distribution using a tabulation.
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
    !!{
    Compute the potential difference between the given {\normalfont \ttfamily coordinates1} and {\normalfont \ttfamily coordinates2} in a spherical mass distribution using a tabulation.
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
    !!{
    Compute the 1D velocity dispersion at the given {\normalfont \ttfamily coordinates} in a spherical mass distribution using a tabulation.
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
    !!{
    Compute the energy within a given {\normalfont \ttfamily radius} in a spherical mass distribution using a tabulation.
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
    !!{    
    Compute the Fourier transform of the density profile at the given {\normalfont \ttfamily wavenumber} in a spherical mass
    distribution using a tabulation.
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
    !!{
    Compute the freefall radius at a given {\normalfont \ttfamily time} in a spherical mass distribution using a tabulation.
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
    !!{
    Compute the rate of increase of freefall radius at a given {\normalfont \ttfamily time} in a spherical mass distribution using a tabulation.
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
    !!{
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
    !!{
    Tabulate the mass distribution.
    !!}
    use :: Coordinates                 , only : coordinateSpherical, assignment(=)
    use :: Input_Paths                 , only : inputPath          , pathTypeDataDynamic
    use :: File_Utilities              , only : File_Lock          , File_Path          , File_Unlock   , lockDescriptor       , &
         &                                      File_Exists        , Directory_Make
    use :: Multi_Counters              , only : multiCounter
    use :: Display                     , only : displayIndent      , displayUnindent    , displayMessage, verbosityLevelWorking, &
         &                                      displayCounter     , displayCounterClear
    use :: Numerical_Constants_Prefixes, only : siFormat
    implicit none
    class           (massDistributionSphericalTabulated )                , intent(inout) :: self
    double precision                                                     , intent(in   ) :: radiusScaled
    double precision                                     , dimension(:  ), intent(in   ) :: parameters
    type            (massDistributionContainer          )                , intent(inout) :: container
    type            (massDistributionTabulation         )                , intent(inout) :: tabulation
    class           (massDistributionSphericalTabulated ), save          , pointer       :: instance
    type            (kinematicsDistributionCollisionless), save          , pointer       :: instanceKinematicsDistribution
    !$omp threadprivate(instance,instanceKinematicsDistribution)
    double precision                                     , dimension(:  ), allocatable   :: parameters_                   , parametersReduced_
    integer         (c_size_t                           ), dimension(:  ), allocatable   :: iParameters
    integer         (c_size_t                           )                                :: lengthMaximum                 , iRadius             , &
         &                                                                                  iterationCount                , iterationCountTotal , &
         &                                                                                  i
    double precision                                                                     :: radius_                       , quantity_           , &
         &                                                                                  time_                         , wavenumber_         , &
         &                                                                                  radiusOuter_                  , density_
    logical                                                                              :: workRemains
    type            (multiCounter                       )                                :: counter
    character       (len= 8                             )                                :: labelLower                    , labelUpper
    character       (len=64                             )                                :: labelSize

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
         call Directory_Make(char(File_Path(char(fileName))))
         if (File_Exists(fileName)) then
            call File_Lock(char(fileName),fileLock,lockIsShared=.true.)
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
               ! Construct radius and parameter ranges.
               tabulation%radiusMinimum        =min(tabulation%radiusMinimum    ,0.5d0*radiusScaled)
               tabulation%radiusMaximum        =max(tabulation%radiusMaximum    ,2.0d0*radiusScaled)
               tabulation%parametersMinimum    =max(min(tabulation%parametersMinimum,0.5d0*parameters  ),container%parametersMinimumLimit)
               tabulation%parametersMaximum    =min(max(tabulation%parametersMaximum,2.0d0*parameters  ),container%parametersMaximumLimit)
               tabulation%countRadii           =int(log10(    tabulation%radiusMaximum/tabulation%    radiusMinimum)*dble(    tabulation%radiusCountPer),kind=c_size_t)+1_c_size_t
               tabulation%countParameters      =int(log10(tabulation%parametersMaximum/tabulation%parametersMinimum)*dble(tabulation%parametersCountPer),kind=c_size_t)+1_c_size_t
               tabulation%radiusInverseStep    =dble(tabulation%countRadii     -1_c_size_t)/log(tabulation%    radiusMaximum/tabulation%    radiusMinimum)
               tabulation%parametersInverseStep=dble(tabulation%countParameters-1_c_size_t)/log(tabulation%parametersMaximum/tabulation%parametersMinimum)
               if (allocated(tabulation%table)) deallocate(tabulation%table)
               select case(size(tabulation%countParameters))
               case (1)
                  allocate(tabulation%table(tabulation%countRadii,tabulation%countParameters(1),1                            ,1                            ))
               case (2)
                  allocate(tabulation%table(tabulation%countRadii,tabulation%countParameters(1),tabulation%countParameters(2),1                            ))
               case (3)
                  allocate(tabulation%table(tabulation%countRadii,tabulation%countParameters(1),tabulation%countParameters(2),tabulation%countParameters(3)))
               case default
                  call Error_Report('rank not supported'//{introspection:location})
               end select
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
               ! Iterate over parameters.
               iterationCount     =0_c_size_t
               iterationCountTotal=product(tabulation%countParameters)
               ! Tabulate in parallel.
               !$omp parallel private(iRadius,radius_,time_,radiusOuter_,wavenumber_,density_,coordinates,coordinatesZeroPoint,quantity_,iParameters,parameters_,parametersReduced_,workRemains,counter)
               ! This is a new thread, so mark it as tabulating.
               tabulating         =.true.
               ! Initialize the counter and iterate over parameter states.
               counter            =multiCounter(tabulation%countParameters)
               do while (.true.)
                  !$omp barrier
                  workRemains=counter%increment()
                  if (.not.workRemains) exit
                  !$omp master
                  call displayCounter(int(100.0d0*dble(iterationCount)/dble(iterationCountTotal)),iterationCount==0,verbosityLevelWorking)
                  iterationCount=iterationCount+1_c_size_t
                  !$omp end master
                  iParameters   =counter%states()
                  parameters_   =exp(log(tabulation%parametersMinimum)+dble(iParameters-1_c_size_t)/tabulation%parametersInverseStep)
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
                     radius_             =exp(log(tabulation%radiusMinimum)+dble(iRadius-1_c_size_t)/tabulation%radiusInverseStep)
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
                     select case(size(tabulation%countParameters))
                     case (1)
                        tabulation%table(iRadius,iParameters(1),1             ,1             )=quantity_
                     case (2)
                        tabulation%table(iRadius,iParameters(1),iParameters(2),1             )=quantity_
                     case (3)
                        tabulation%table(iRadius,iParameters(1),iParameters(2),iParameters(3))=quantity_
                     case default
                        call Error_Report('rank not supported'//{introspection:location})
                     end select
                  end do
                  !$omp end do
                  deallocate(instance)
                  nullify   (instance)
               end do
               !$omp master
               call displayCounterClear(verbosityLevelWorking)
               !$omp end master
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
      !!{
      Test if the mass profile must be retabulated.
      !!}
      implicit none
      
      retabulate=.not.allocated(tabulation%table)
      if (.not.retabulate)                                                                                                                       &
           &  retabulate=     radiusScaled < tabulation%radiusMinimum                                                                            &
           &             .or.                                                                                                                    &
           &                  radiusScaled > tabulation%radiusMaximum                                                                            &
           &             .or.                                                                                                                    &
           &              any(parameters   < tabulation%parametersMinimum .and. tabulation%parametersMinimum > container%parametersMinimumLimit) &
           &             .or.                                                                                                                    &
           &              any(parameters   > tabulation%parametersMaximum .and. tabulation%parametersMaximum < container%parametersMaximumLimit)
      return
    end function retabulate

  end subroutine sphericalTabulatedTabulate

  double precision function sphericalTabulatedInterpolate(self,radiusScaled,parameters,tabulation) result(interpolated)
    !!{
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

    ! Compute interpolating factors.
    allocate(hParameters(3,size(parameters)))
    hRadius    (2  )=min(log(radiusScaled/tabulation%radiusMinimum    )*tabulation%radiusInverseStep    +1.0d0,dble(tabulation%countRadii     -1_c_size_t))
    hParameters(2,:)=min(log(parameters  /tabulation%parametersMinimum)*tabulation%parametersInverseStep+1.0d0,dble(tabulation%countParameters-1_c_size_t))
    iRadius         =int(hRadius    (2  ),kind=c_size_t)
    iParameters     =int(hParameters(2,:),kind=c_size_t)
    hRadius    (2  )=hRadius    (2  )-dble(iRadius    )
    hParameters(2,:)=hParameters(2,:)-dble(iParameters)
    hRadius    (1  )=1.0d0-hRadius    (2  )
    hParameters(1,:)=1.0d0-hParameters(2,:)    
    ! Perform the interpolation.
    interpolated=0.0d0
    counter     =multiCounter(spread(2_c_size_t,1,size(parameters)))
    select case(size(parameters))
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
    !!{
    Return true if this thread is currently tabulating.
    !!}
    implicit none
    class(massDistributionSphericalTabulated), intent(inout) :: self
    !$GLC attributes unused :: self
    
    isTabulating=tabulating
    return
  end function sphericalTabulatedIsTabulating
  
  subroutine sphericalTabulatedFileRead(self,fileName,quantityName,container,tabulation)
    !!{
    Read tabulated data from file.
    !!}
    use :: HDF5_Access    , only : hdf5Access
    use :: IO_HDF5        , only : hdf5Object
    use :: String_Handling, only : String_Upper_Case_First
    use :: Display        , only : displayMessage         , verbosityLevelWorking
    implicit none
    class  (massDistributionSphericalTabulated), intent(inout) :: self
    type   (varying_string                    ), intent(in   ) :: fileName  , quantityName
    type   (massDistributionContainer         ), intent(inout) :: container
    type   (massDistributionTabulation        ), intent(inout) :: tabulation
    type   (hdf5Object                        )                :: file
    integer(c_size_t                          )                :: i

    if (allocated(tabulation%table)) deallocate(tabulation%table)
    call displayMessage("reading tabulated "//enumerationQuantityDecode(tabulation%quantity)//" profile from '"//char(fileName)//"'",verbosityLevelWorking)
    if (.not.allocated(tabulation%parametersMinimum    )) allocate(tabulation%parametersMinimum    (container%countParameters(tabulation)))
    if (.not.allocated(tabulation%parametersMaximum    )) allocate(tabulation%parametersMaximum    (container%countParameters(tabulation)))
    if (.not.allocated(tabulation%parametersInverseStep)) allocate(tabulation%parametersInverseStep(container%countParameters(tabulation)))
    if (.not.allocated(tabulation%countParameters      )) allocate(tabulation%countParameters      (container%countParameters(tabulation)))
    !$ call hdf5Access%set()
    file=hdf5Object(char(fileName),readOnly=.true.)
    call    file%readAttribute(char(quantityName)//'RadiusMinimum'                                                                    ,tabulation%radiusMinimum           )
    call    file%readAttribute(char(quantityName)//'RadiusMaximum'                                                                    ,tabulation%radiusMaximum           )
    call    file%readAttribute(char(quantityName)//'RadiusInverseStep'                                                                ,tabulation%radiusInverseStep       )
    do i=1,container%countParameters(tabulation)
       call file%readAttribute(char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation)))//'Minimum'    ,tabulation%parametersMinimum    (i))
       call file%readAttribute(char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation)))//'Maximum'    ,tabulation%parametersMaximum    (i))
       call file%readAttribute(char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation)))//'InverseStep',tabulation%parametersInverseStep(i))
    end do
    call    file%readDataset  (char(quantityName)                                                                                     ,tabulation%table                   )
    !$ call hdf5Access%unset()
    tabulation   %countRadii        =size(tabulation%table,dim=  1)
    do i=1,container%countParameters(tabulation)
       tabulation%countParameters(i)=size(tabulation%table,dim=i+1)
    end do
    return
  end subroutine sphericalTabulatedFileRead

  subroutine sphericalTabulatedFileWrite(self,fileName,quantityName,container,tabulation)
    !!{
    Read tabulated data from file.
    !!}
    use :: HDF5_Access    , only : hdf5Access
    use :: IO_HDF5        , only : hdf5Object
    use :: String_Handling, only : String_Upper_Case_First
    use :: Display        , only : displayMessage         , verbosityLevelWorking
    implicit none
    class  (massDistributionSphericalTabulated), intent(inout) :: self
    type   (varying_string                    ), intent(in   ) :: fileName  , quantityName
    type   (massDistributionContainer         ), intent(inout) :: container
    type   (massDistributionTabulation        ), intent(in   ) :: tabulation
    type   (hdf5Object                        )                :: file
    integer(c_size_t                          )                :: i

    call displayMessage("writing tabulated "//char(quantityName)//" profile to '"//char(fileName)//"'",verbosityLevelWorking)
    !$ call hdf5Access%set()
   file=hdf5Object(char(fileName),overWrite=.true.)
    call    file%writeAttribute(tabulation%radiusMinimum           ,char(quantityName)//'RadiusMinimum'                                                                                                                  )
    call    file%writeAttribute(tabulation%radiusMaximum           ,char(quantityName)//'RadiusMaximum'                                                                                                                  )
    call    file%writeAttribute(tabulation%radiusInverseStep       ,char(quantityName)//'RadiusInverseStep'                                                                                                              )
    do i=1,container%countParameters(tabulation)
       call file%writeAttribute(tabulation%parametersMinimum    (i),char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation)))//'Minimum'                                                  )
       call file%writeAttribute(tabulation%parametersMaximum    (i),char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation)))//'Maximum'                                                  )
       call file%writeAttribute(tabulation%parametersInverseStep(i),char(quantityName)//String_Upper_Case_First(char(container%nameParameter(i,tabulation)))//'InverseStep'                                              )
    end do
    call    file%writeDataset  (tabulation%table                   ,char(quantityName)                                                                                     ,'Tabulated '//char(quantityName)//' profile.')
    !$ call hdf5Access%unset()
    return
  end subroutine sphericalTabulatedFileWrite

  subroutine massDistributionContainerInitialize(self,countParameters)
    !!{
    Initialize the container to the specified number of parameters.
    !!}
    implicit none
    class  (massDistributionContainer), intent(inout) :: self
    integer                           , intent(in   ) :: countParameters

    ! Allocate the internal arrays.
    allocate(self                           %descriptionParameters (countParameters  ))
    allocate(self                           %nameParameters        (countParameters  ))
    allocate(self                           %parametersMinimumLimit(countParameters+1))
    allocate(self                           %parametersMaximumLimit(countParameters+1))
    allocate(self%mass                      %parametersMinimum     (countParameters  ))
    allocate(self%mass                      %parametersMaximum     (countParameters  ))
    allocate(self%mass                      %parametersCountPer    (countParameters  ))
    allocate(self%mass                      %countParameters       (countParameters  ))
    allocate(self%radiusEnclosingDensity    %parametersMinimum     (countParameters  ))
    allocate(self%radiusEnclosingDensity    %parametersMaximum     (countParameters  ))
    allocate(self%radiusEnclosingDensity    %parametersCountPer    (countParameters  ))
    allocate(self%radiusEnclosingDensity    %countParameters       (countParameters  ))
    allocate(self%energy                    %parametersMinimum     (countParameters  ))
    allocate(self%energy                    %parametersMaximum     (countParameters  ))
    allocate(self%energy                    %parametersCountPer    (countParameters  ))
    allocate(self%energy                    %countParameters       (countParameters  ))
    allocate(self%potential                 %parametersMinimum     (countParameters  ))
    allocate(self%potential                 %parametersMaximum     (countParameters  ))
    allocate(self%potential                 %parametersCountPer    (countParameters  ))
    allocate(self%potential                 %countParameters       (countParameters  ))
    allocate(self%velocityDispersion1D      %parametersMinimum     (countParameters  ))
    allocate(self%velocityDispersion1D      %parametersMaximum     (countParameters  ))
    allocate(self%velocityDispersion1D      %parametersCountPer    (countParameters  ))
    allocate(self%velocityDispersion1D      %countParameters       (countParameters  ))
    allocate(self%radiusFreefall            %parametersMinimum     (countParameters  ))
    allocate(self%radiusFreefall            %parametersMaximum     (countParameters  ))
    allocate(self%radiusFreefall            %parametersCountPer    (countParameters  ))
    allocate(self%radiusFreefall            %countParameters       (countParameters  ))
    allocate(self%radiusFreefallIncreaseRate%parametersMinimum     (countParameters  ))
    allocate(self%radiusFreefallIncreaseRate%parametersMaximum     (countParameters  ))
    allocate(self%radiusFreefallIncreaseRate%parametersCountPer    (countParameters  ))
    allocate(self%radiusFreefallIncreaseRate%countParameters       (countParameters  ))
    allocate(self%densityRadialMoment0      %parametersMinimum     (countParameters  ))
    allocate(self%densityRadialMoment0      %parametersMaximum     (countParameters  ))
    allocate(self%densityRadialMoment0      %parametersCountPer    (countParameters  ))
    allocate(self%densityRadialMoment0      %countParameters       (countParameters  ))
    allocate(self%densityRadialMoment1      %parametersMinimum     (countParameters  ))
    allocate(self%densityRadialMoment1      %parametersMaximum     (countParameters  ))
    allocate(self%densityRadialMoment1      %parametersCountPer    (countParameters  ))
    allocate(self%densityRadialMoment1      %countParameters       (countParameters  ))
    allocate(self%densityRadialMoment2      %parametersMinimum     (countParameters  ))
    allocate(self%densityRadialMoment2      %parametersMaximum     (countParameters  ))
    allocate(self%densityRadialMoment2      %parametersCountPer    (countParameters  ))
    allocate(self%densityRadialMoment2      %countParameters       (countParameters  ))
    allocate(self%densityRadialMoment3      %parametersMinimum     (countParameters  ))
    allocate(self%densityRadialMoment3      %parametersMaximum     (countParameters  ))
    allocate(self%densityRadialMoment3      %parametersCountPer    (countParameters  ))
    allocate(self%densityRadialMoment3      %countParameters       (countParameters  ))
    allocate(self%fourierTransform          %parametersMinimum     (countParameters+1))
    allocate(self%fourierTransform          %parametersMaximum     (countParameters+1))
    allocate(self%fourierTransform          %parametersCountPer    (countParameters+1))
    allocate(self%fourierTransform          %countParameters       (countParameters+1))
    ! Initialize minima/maxima of tabulation ranges. These are chosen to ensure that the first tabulation will force these to
    ! be reset.
    self                           %parametersMinimumLimit=-huge(0.0d0)
    self                           %parametersMaximumLimit=+huge(0.0d0)
    self%mass                      %radiusMinimum         =+huge(0.0d0)
    self%mass                      %radiusMaximum         =-huge(0.0d0)
    self%radiusEnclosingDensity    %radiusMinimum         =+huge(0.0d0)
    self%radiusEnclosingDensity    %radiusMaximum         =-huge(0.0d0)
    self%energy                    %radiusMinimum         =+huge(0.0d0)
    self%energy                    %radiusMaximum         =-huge(0.0d0)
    self%potential                 %radiusMinimum         =+huge(0.0d0)
    self%potential                 %radiusMaximum         =-huge(0.0d0)
    self%velocityDispersion1D      %radiusMinimum         =+huge(0.0d0)
    self%velocityDispersion1D      %radiusMaximum         =-huge(0.0d0)
    self%radiusFreefall            %radiusMinimum         =+huge(0.0d0)
    self%radiusFreefall            %radiusMaximum         =-huge(0.0d0)
    self%radiusFreefallIncreaseRate%radiusMinimum         =+huge(0.0d0)
    self%radiusFreefallIncreaseRate%radiusMaximum         =-huge(0.0d0)
    self%densityRadialMoment0      %radiusMinimum         =+huge(0.0d0)
    self%densityRadialMoment0      %radiusMaximum         =-huge(0.0d0)
    self%densityRadialMoment1      %radiusMinimum         =+huge(0.0d0)
    self%densityRadialMoment1      %radiusMaximum         =-huge(0.0d0)
    self%densityRadialMoment2      %radiusMinimum         =+huge(0.0d0)
    self%densityRadialMoment2      %radiusMaximum         =-huge(0.0d0)
    self%densityRadialMoment3      %radiusMinimum         =+huge(0.0d0)
    self%densityRadialMoment3      %radiusMaximum         =-huge(0.0d0)
    self%fourierTransform          %radiusMinimum         =+huge(0.0d0)
    self%fourierTransform          %radiusMaximum         =-huge(0.0d0)
    self%mass                      %parametersMinimum     =+huge(0.0d0)
    self%mass                      %parametersMaximum     =-huge(0.0d0)
    self%radiusEnclosingDensity    %parametersMinimum     =+huge(0.0d0)
    self%radiusEnclosingDensity    %parametersMaximum     =-huge(0.0d0)
    self%energy                    %parametersMinimum     =+huge(0.0d0)
    self%energy                    %parametersMaximum     =-huge(0.0d0)
    self%potential                 %parametersMinimum     =+huge(0.0d0)
    self%potential                 %parametersMaximum     =-huge(0.0d0)
    self%velocityDispersion1D      %parametersMinimum     =+huge(0.0d0)
    self%velocityDispersion1D      %parametersMaximum     =-huge(0.0d0)
    self%radiusFreefall            %parametersMinimum     =+huge(0.0d0)
    self%radiusFreefall            %parametersMaximum     =-huge(0.0d0)
    self%radiusFreefallIncreaseRate%parametersMinimum     =+huge(0.0d0)
    self%radiusFreefallIncreaseRate%parametersMaximum     =-huge(0.0d0)
    self%densityRadialMoment0      %parametersMinimum     =+huge(0.0d0)
    self%densityRadialMoment0      %parametersMaximum     =-huge(0.0d0)
    self%densityRadialMoment1      %parametersMinimum     =+huge(0.0d0)
    self%densityRadialMoment1      %parametersMaximum     =-huge(0.0d0)
    self%densityRadialMoment2      %parametersMinimum     =+huge(0.0d0)
    self%densityRadialMoment2      %parametersMaximum     =-huge(0.0d0)
    self%densityRadialMoment3      %parametersMinimum     =+huge(0.0d0)
    self%densityRadialMoment3      %parametersMaximum     =-huge(0.0d0)
    self%fourierTransform          %parametersMinimum     =+huge(0.0d0)
    self%fourierTransform          %parametersMaximum     =-huge(0.0d0)
    return
  end subroutine massDistributionContainerInitialize

  function massDistributionContainerNameParameter(self,indexParameter,tabulation) result(nameParameter)
    !!{
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
    !!{
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
