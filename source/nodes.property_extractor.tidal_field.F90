!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!+    Contributions to this file made by: Andrew Benson, Charles Gannon.

!!{
Implements a tidal field property extractor class.
!!}

  !![
  <nodePropertyExtractor name="nodePropertyExtractorTidalField">
   <description> 
    A property extractor class which extracts the radial component of the tidal tensor assuming spherical symmetry in units of $\mathrm{Gyr^{-2}}$.
    This is calculated using the equation:
    \begin{equation}
     \mathcal{F} = {\mathrm{G} M_\mathrm{host}(&lt;r) \over r^3} - 4 \pi \mathrm{G}
     \rho_\mathrm{host}(r) + \omega^2,
    \end{equation}
    where $r$ is the current orbital radius. $M_\mathrm{host}(&lt;r)$ is the mass of the host halo enclosed within a sphere
    of radius $r$, $\rho_\mathrm{host}(r)$ is the host density at radius $r$, and $\omega$ is the orbital angular velocity.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorTidalField
     !!{
      A property extractor class which extracts the radial component of the tidal tensor assuming spherical symmetry in units of $\mathrm{Gyr^{-2}}$.
     !!}
   contains
     final     ::                tidalFieldDestructor
     procedure :: extract     => tidalFieldExtract
     procedure :: name        => tidalFieldName
     procedure :: description => tidalFieldDescription
     procedure :: unitsInSI   => tidalFieldUnitsInSI
  end type nodePropertyExtractorTidalField

  interface nodePropertyExtractorTidalField
     !!{
     Constructors for the \refClass{nodePropertyExtractorTidalField} output analysis class.
     !!}
     module procedure tidalFieldConstructorParameters
     module procedure tidalFieldConstructorInternal
  end interface nodePropertyExtractorTidalField

contains

  function tidalFieldConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nodePropertyExtractorTidalField} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorTidalField)                :: self
    type (inputParameters                ), intent(inout) :: parameters
    !$GLC attributes unused :: parameters

    self=nodePropertyExtractorTidalField()
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function tidalFieldConstructorParameters

  function tidalFieldConstructorInternal() result(self)
    !!{
    Internal constructor for the \refClass{nodePropertyExtractorTidalField} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorTidalField  ) :: self
    !$GLC attributes unused :: self

    return
  end function tidalFieldConstructorInternal

  subroutine tidalFieldDestructor(self)
    !!{
    Destructor for the \refClass{nodePropertyExtractorTidalField} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorTidalField), intent(inout) :: self
    !$GLC attributes unused :: self

    return
  end subroutine tidalFieldDestructor

  double precision function tidalFieldExtract(self,node,instance)
    !!{
    Return the radial part of the tidal tensor for satellite halos assuming spherical symmetry of the host.
    !!}
    use :: coordinates                     , only : coordinatespherical                     , assignment(=)
    use :: galacticus_nodes                , only : nodecomponentsatellite        , treenode
    use :: mass_distributions              , only : massdistributionclass
    use :: numerical_constants_math        , only : pi
    use :: numerical_constants_astronomical, only : gravitationalconstant_internal, gigaYear,               &
         &                                          megaParsec
    use :: Numerical_Constants_Prefixes    , only : kilo 
    use :: vectors                         , only : vector_magnitude
    implicit none
    class           (nodePropertyExtractorTidalField), intent(inout), target   :: self
    type            (treeNode                       ), intent(inout), target   :: node
    type            (multiCounter                   ), intent(inout), optional :: instance 
    type            (treeNode                       ), pointer                 :: nodeHost
    class           (nodeComponentSatellite         ), pointer                 :: satellite
    class           (massDistributionClass          ), pointer                 :: massDistribution_
    double precision                                                           :: densityHost       , enclosedMassHost, &
         &                                                                        radiusOrbital     , velocityOrbital
    type            (coordinateSpherical            )                          :: coordinatesOrbital
    !$glc attributes unused :: self, instance

    ! For isolated halos, always return zero tidal field.
    if (node%isSatellite()) then
       ! Find the host node.
       nodeHost  => node%parent
       ! Get the satellite component.
       satellite => node%satellite()
       ! Get the postition components

       ! Compute orbital radius and velocity
       radiusOrbital      =  +Vector_Magnitude(satellite%position())
       velocityOrbital    =  +Vector_Magnitude(satellite%velocity())

       coordinatesOrbital =  [radiusOrbital,0.0d0,0.0d0]
       massDistribution_  => nodeHost         %massDistribution    (                  )
       densityHost        =  massDistribution_%density             (coordinatesOrbital)
       enclosedMassHost   =  massDistribution_%massEnclosedBySphere(     radiusOrbital)
       !![
       <objectDestructor name="massDistribution_"/>
       !!]
       ! Compute the tidal field.
       tidalFieldExtract  = +gravitationalConstant_internal*enclosedMassHost/radiusOrbital **3 &
            &                   -4.0d0*Pi*gravitationalConstant_internal*densityHost           &
            &                   +(velocityOrbital/radiusOrbital)**2

       ! Convert to Gyr⁻² units.
       tidalFieldExtract  = +tidalFieldExtract                                                 &
            &                   *gigaYear**2                                                   &
            &                   *kilo**2                                                       &
            &                   /megaParsec**2        
    else
       tidalFieldExtract  = +0.0d0
    end if
    return
  end function tidalFieldExtract

  function tidalFieldName(self)
    !!{
    Return the name of the tidal radius property.
    !!}
    implicit none
    type (varying_string                  )                         :: tidalFieldName
    class(nodePropertyExtractorTidalField), intent(inout)           :: self
    !$GLC attributes unused :: self

    tidalFieldName=var_str('satelliteTidalField')
    return
  end function tidalFieldName

  function tidalFieldDescription(self)
    !!{
    Return a description of the tidal radius property.
    !!}
    implicit none
    type (varying_string                 )                :: tidalFieldDescription
    class(nodePropertyExtractorTidalField), intent(inout) :: self
    !$GLC attributes unused :: self

    tidalFieldDescription=var_str('Tidal field in the halo Gyr⁻².')
    return
  end function tidalFieldDescription

  double precision function tidalFieldUnitsInSI(self)
    !!{
    Return the units of the tidal radius property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorTidalField), intent(inout) :: self
    !$GLC attributes unused :: self

    tidalFieldUnitsInSI= 1 / gigaYear**2
    return
  end function tidalFieldUnitsInSI


