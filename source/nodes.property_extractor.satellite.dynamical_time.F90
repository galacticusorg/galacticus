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

!+    Contributions to this file made by: Charles Gannon
  
!!{
Provides a class that implements a satellite dynamical time extractor.
!!}

  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusClass

  !![
  <nodePropertyExtractor name="nodePropertyExtractorSatelliteDynamicalTime">
   <description>
     A satellite dynamical time extractor class. Extracts the satellite dynamical time in units of Gyr, following the definition
     from \cite{binney_galactic_2008}:     
     \begin{equation}
     \tau_\mathrm{dyn} = \sqrt{\frac{3 \pi}{16 \mathrm{G} \rho_\mathrm{tidal}}} = \sqrt{\frac{\pi^2 r_\mathrm{tidal}^3}{4 \mathrm{G} M_\mathrm{tidal}}},
     \end{equation}     
     with $\rho_\mathrm{tidal}$ being the density within the tidal radius, $r_\mathrm{tidal}$, of the satellite which encloses a
     mass $M_\mathrm{tidal}$.
   </description>
  </nodePropertyExtractor>
  !!]
  type, extends(nodePropertyExtractorScalar) :: nodePropertyExtractorSatelliteDynamicalTime
     !!{
     A dynamical time extractor class.
     !!}
     private
     class(satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
   contains
     final     ::                dynamicalTimeDestructor
     procedure :: extract     => dynamicalTimeExtract
     procedure :: name        => dynamicalTimeName
     procedure :: description => dynamicalTimeDescription
     procedure :: unitsInSI   => dynamicalTimeUnitsInSI
  end type nodePropertyExtractorSatelliteDynamicalTime

  interface nodePropertyExtractorSatelliteDynamicalTime
     !!{
     Constructors for the {\normalfont \ttfamily satelliteDynamicalTime} output analysis class.
     !!}
     module procedure dynamicalTimeConstructorParameters
     module procedure dynamicalTimeConstructorInternal
  end interface nodePropertyExtractorSatelliteDynamicalTime

contains

  function dynamicalTimeConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily satelliteDynamicalTime} property extractor class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type (nodePropertyExtractorSatelliteDynamicalTime)                :: self
    type (inputParameters                            ), intent(inout) :: parameters
    class(satelliteTidalStrippingRadiusClass         ), pointer       :: satelliteTidalStrippingRadius_

    !![
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
    !!]
    self=nodePropertyExtractorSatelliteDynamicalTime(satelliteTidalStrippingRadius_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    !!]
    return
  end function dynamicalTimeConstructorParameters

  function dynamicalTimeConstructorInternal(satelliteTidalStrippingRadius_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily satelliteDynamicalTime} property extractor class.
    !!}
    implicit none
    type (nodePropertyExtractorSatelliteDynamicalTime)                        :: self
    class(satelliteTidalStrippingRadiusClass         ), intent(in   ), target :: satelliteTidalStrippingRadius_
    !![
    <constructorAssign variables="*satelliteTidalStrippingRadius_"/>
    !!]

    return
  end function dynamicalTimeConstructorInternal

  subroutine dynamicalTimeDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily satelliteDynamicalTime} property extractor class.
    !!}
    implicit none
    type(nodePropertyExtractorSatelliteDynamicalTime), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    !!]
    return
  end subroutine dynamicalTimeDestructor

  double precision function dynamicalTimeExtract(self,node,instance)
    !!{
    Implement a dynamical time property extractor.
    !!}
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, MpcPerKmPerSToGyr
    implicit none
    class           (nodePropertyExtractorSatelliteDynamicalTime), intent(inout), target   :: self
    type            (treeNode                                   ), intent(inout), target   :: node
    type            (multiCounter                               ), intent(inout), optional :: instance
    class           (massDistributionClass                      )               , pointer  :: massDistribution_
    double precision                                                                       :: radiusTidal      , massTidal
    !$GLC attributes unused :: instance

    dynamicalTimeExtract=-1.0d0
    radiusTidal         =self%satelliteTidalStrippingRadius_%radius(node)
    if (radiusTidal <= 0.0d0) return
    massDistribution_ => node             %massDistribution    (           )
    massTidal         =  massDistribution_%massEnclosedBySphere(radiusTidal)
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    if (massTidal <= 0.0d0) return
    dynamicalTimeExtract  =+sqrt(                                   &  
         &                       +Pi                            **2 &
         &                       /4.0d0                             &
         &                       *radiusTidal**3                    &
         &                       /gravitationalConstant_internal    &
         &                       /massTidal                         &
         &                       *MpcPerKmPerSToGyr             **2 &
         &                      )
    return
  end function dynamicalTimeExtract

  function dynamicalTimeName(self)
    !!{
    Return the name of the satellite dynamical time property.
    !!}
    implicit none
    type (varying_string                             )                :: dynamicalTimeName
    class(nodePropertyExtractorSatelliteDynamicalTime), intent(inout) :: self
    !$GLC attributes unused :: self

    dynamicalTimeName=var_str('satelliteDynamicalTime')
    return
  end function dynamicalTimeName

  function dynamicalTimeDescription(self)
    !!{
    Return a description of the satellite dynamical time property.
    !!}
    implicit none
    type (varying_string                             )                :: dynamicalTimeDescription
    class(nodePropertyExtractorSatelliteDynamicalTime), intent(inout) :: self
    !$GLC attributes unused :: self

    dynamicalTimeDescription=var_str('Dynamical time of satellite [Gyr].')
    return
  end function dynamicalTimeDescription

  double precision function dynamicalTimeUnitsInSI(self)
    !!{
    Return the units of the satellite dynamical time property in the SI system.
    !!}
    use :: Numerical_Constants_Astronomical, only : gigaYear
    implicit none
    class(nodePropertyExtractorSatelliteDynamicalTime), intent(inout) :: self
    !$GLC attributes unused :: self

    dynamicalTimeUnitsInSI=gigaYear
    return
  end function dynamicalTimeUnitsInSI




