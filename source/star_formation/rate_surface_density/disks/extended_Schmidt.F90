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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson, Claude.

  !!{RST
  Implementation of the extended Schmidt star formation rate surface density law of :cite:t:`shi_extended_2011` for galactic disks.
  !!}
  
  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksExtendedSchmidt" docformat="rst">
   <description>
   A star formation rate surface density class implementing the extended Schmidt law :cite:p:`shi_extended_2011`:

   .. math::

      \dot{\Sigma}_\star = A \left({\Sigma_\mathrm{gas}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)^{N_1}
      \left({\Sigma_{\star}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)^{N_2}

   where :math:`A=`\ ``[normalization]``, :math:`N_1=`\ ``[exponentGas]`` and :math:`N_2=`\ ``[exponentStars]`` are parameters. Note that :math:`\Sigma_\mathrm{gas}` here is the *total* gas surface density, not that of hydrogen alone, since the gas masses of :cite:t:`shi_extended_2011` include a factor of 1.36 to account for elements heavier than hydrogen.
   </description>
  </starFormationRateSurfaceDensityDisks>
  !!]
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksExtendedSchmidt
     !!{RST
     Implementation of the extended Schmidt star formation rate surface density law of :cite:t:`shi_extended_2011` for galactic disks.
     !!}
     private
     double precision :: normalization, exponentGas, &
          &              exponentStars
   contains
     procedure :: rate => extendedSchmidtRate
  end type starFormationRateSurfaceDensityDisksExtendedSchmidt

  interface starFormationRateSurfaceDensityDisksExtendedSchmidt
     !!{RST
     Constructors for the :galacticus-class:`starFormationRateSurfaceDensityDisksExtendedSchmidt` star formation surface density rate in disks class.
     !!}
     module procedure extendedSchmidtConstructorParameters
     module procedure extendedSchmidtConstructorInternal
  end interface starFormationRateSurfaceDensityDisksExtendedSchmidt

contains

  function extendedSchmidtConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`starFormationRateSurfaceDensityDisksExtendedSchmidt` star formation surface density rate in disks class which takes a parameter set as input.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (starFormationRateSurfaceDensityDisksExtendedSchmidt)                :: self
    type            (inputParameters                                    ), intent(inout) :: parameters
    double precision                                                                     :: normalization, exponentGas, &
         &                                                                                  exponentStars

    !![
    <inputParameter docformat="rst">
      <name>normalization</name>
      <defaultSource>
      :cite:p:`shi_extended_2011`
      </defaultSource>
      <defaultValue>0.5248d-1</defaultValue>
      <description>
      The normalization of the extended Schmidt star formation law [:math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}`\ pc\ :math:`^{-2}`].
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentGas</name>
      <defaultSource>
      :cite:p:`shi_extended_2011`
      </defaultSource>
      <defaultValue>1.0000d+0</defaultValue>
      <description>
      The exponent of gas surface density in the extended Schmidt star formation law.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>exponentStars</name>
      <defaultSource>
      :cite:p:`shi_extended_2011`
      </defaultSource>
      <defaultValue>0.4800d+0</defaultValue>
      <description>
      The exponent of stellar surface density in the extended Schmidt star formation law.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateSurfaceDensityDisksExtendedSchmidt(normalization,exponentGas,exponentStars)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function extendedSchmidtConstructorParameters

  function extendedSchmidtConstructorInternal(normalization,exponentGas,exponentStars) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`starFormationRateSurfaceDensityDisksExtendedSchmidt` star formation surface density rate in disks class.
    !!}
    use :: Numerical_Constants_Prefixes, only : mega
    implicit none
    type            (starFormationRateSurfaceDensityDisksExtendedSchmidt)                :: self
    double precision                                                     , intent(in   ) :: normalization, exponentGas, &
         &                                                                                  exponentStars
    !![
    <constructorAssign variables="normalization, exponentGas, exponentStars"/>
    !!]

    ! Renormalize the relation to internal units.
    self%normalization=+self%normalization                    &
         &             *(mega**2)                             & ! Convert from M☉/pc²/Gyr to M☉/Mpc²/Gyr
         &             *((1.0d0/mega**2)**self%exponentStars) & ! Unit conversion for stars.
         &             *((1.0d0/mega**2)**self%exponentGas  )   ! Unit conversion for gas.
    return
  end function extendedSchmidtConstructorInternal

  double precision function extendedSchmidtRate(self,node,radius)
    !!{RST
    Returns the star formation rate surface density (in :math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}` Mpc\ :math:`^{-2}`) for star formation in the galactic disk of ``node``. The disk is assumed to obey the extended Schmidt law of :cite:t:`shi_extended_2011`:

    .. math::

       \dot{\Sigma}_\star = A \left({\Sigma_\mathrm{gas}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)
       ^{N_1} \left({\Sigma_{\star}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)^{N_2},

    where :math:`A=`\ ``[normalization]``, :math:`N_1=`\ ``[exponentGas]``, and :math:`N_2=`\ ``[exponentStars]``.
    !!}
    use :: Coordinates               , only : coordinateCylindrical, assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeGaseous, massTypeStellar
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (starFormationRateSurfaceDensityDisksExtendedSchmidt), intent(inout) :: self
    type            (treeNode                                           ), intent(inout) :: node
    double precision                                                     , intent(in   ) :: radius
    class           (massDistributionClass                              ), pointer       :: massDistributionGaseous, massDistributionStellar
    double precision                                                                     :: surfaceDensityGas      , surfaceDensityStellar
    type            (coordinateCylindrical                              )                :: coordinates

    ! Return zero rate for non-positive radius.
    if (radius <= 0.0d0) then
       extendedSchmidtRate=0.0d0
       return
    end if
    ! Get stellar and gas surface densities.
    coordinates             =  [radius,0.0d0,0.0d0]
    massDistributionGaseous => node                   %massDistribution(componentType=componentTypeDisk,massType=massTypeGaseous)
    massDistributionStellar => node                   %massDistribution(componentType=componentTypeDisk,massType=massTypeStellar)
    surfaceDensityGas       =  massDistributionGaseous%surfaceDensity  (              coordinates                               )
    surfaceDensityStellar   =  massDistributionStellar%surfaceDensity  (              coordinates                               )
    !![
    <objectDestructor name="massDistributionGaseous"/>
    <objectDestructor name="massDistributionStellar"/>
    !!]
    ! Compute the star formation rate surface density.
    ! Note that the gas surface densities of :cite:t:`shi_extended_2011` include a factor of 1.36 to account for elements heavier
    ! than hydrogen, so the law is applied here to the total gas surface density, not to that of hydrogen alone.
    extendedSchmidtRate=+self%normalization                        & ! Normalization of the star formation rate.
         &              *surfaceDensityGas    **self%exponentGas   &
         &              *surfaceDensityStellar**self%exponentStars
    return
  end function extendedSchmidtRate
