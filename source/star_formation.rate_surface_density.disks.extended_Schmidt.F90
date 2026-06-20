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

!+    Contributions to this file made by: Arya Farahi, Andrew Benson.

  !!{RST
  Implementation of the extended Schmidt star formation rate surface density law of :cite:t:`shi_extended_2011` for galactic disks.
  !!}
  
  !![
  <starFormationRateSurfaceDensityDisks name="starFormationRateSurfaceDensityDisksExtendedSchmidt" docformat="rst">
   <description>
   A star formation rate surface density class implementing the extended Schmidt law :cite:p:`shi_extended_2011`:

   .. math::

      \dot{\Sigma}_\star = A \left(x_\mathrm{H} {\Sigma_\mathrm{gas}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)^{N_1}
      \left({\Sigma_{\star}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)^{N_2}

   where :math:`A=`\ ``[normalization]``, :math:`N_1=`\ ``[exponentGas]`` and :math:`N_2=`\ ``[exponentStars]`` are parameters.
   </description>
  </starFormationRateSurfaceDensityDisks>
  !!]
  type, extends(starFormationRateSurfaceDensityDisksClass) :: starFormationRateSurfaceDensityDisksExtendedSchmidt
     !!{RST
     Implementation of the extended Schmidt star formation rate surface density law of :cite:t:`shi_extended_2011` for galactic disks.
     !!}
     private
     integer         (kind_int8)          :: lastUniqueID
     logical                              :: factorsComputed
     double precision                     :: normalization  , exponentGas         , &
          &                                  exponentStars  , hydrogenMassFraction
   contains
     !![
     <methods docformat="rst">
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                     extendedSchmidtDestructor
     procedure :: autoHook         => extendedSchmidtAutoHook
     procedure :: calculationReset => extendedSchmidtCalculationReset
     procedure :: rate             => extendedSchmidtRate
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
      <defaultValue>0.5248d-10</defaultValue>
      <description>
      The normalization of the extended Schmidt star formation law [:math:`\mathrm{M}_\odot` yr\ :math:`^{-1}`\ pc\ :math:`^{-2}`].
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
    use :: Numerical_Constants_Prefixes, only : giga, mega
    implicit none
    type            (starFormationRateSurfaceDensityDisksExtendedSchmidt)                :: self
    double precision                                                     , intent(in   ) :: normalization, exponentGas, &
         &                                                                                  exponentStars
    !![
    <constructorAssign variables="normalization, exponentGas, exponentStars"/>
    !!]

    self%lastUniqueID   =-1_kind_int8
    self%factorsComputed=.false.
    ! Renormalize the relation to internal units.
    self%normalization=+self%normalization                    &
         &             *(mega**2)*giga                        & ! Convert from M☉/pc²/yr to M☉/Mpc²/Gyr
         &             *((1.0d0/mega**2)**self%exponentStars) & ! Unit conversion for stars.
         &             *((1.0d0/mega**2)**self%exponentGas  )   ! Hydrogen fraction and unit conversion for gas.
    return
  end function extendedSchmidtConstructorInternal

  subroutine extendedSchmidtAutoHook(self)
    !!{RST
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(starFormationRateSurfaceDensityDisksExtendedSchmidt), intent(inout) :: self

    call calculationResetEvent%attach(self,extendedSchmidtCalculationReset,openMPThreadBindingAllLevels,label='starFormationRateSurfaceDensityDisksExtendedS\chmidt')
    return
  end subroutine extendedSchmidtAutoHook

  subroutine extendedSchmidtDestructor(self)
    !!{RST
    Destructor for the extendedSchmidt cooling radius class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(starFormationRateSurfaceDensityDisksExtendedSchmidt), intent(inout) :: self

    if (calculationResetEvent%isAttached(self,extendedSchmidtCalculationReset)) call calculationResetEvent%detach(self,extendedSchmidtCalculationReset)
    return
  end subroutine extendedSchmidtDestructor

  subroutine extendedSchmidtCalculationReset(self,node,uniqueID)
    !!{RST
    Reset the Kennicutt-Schmidt relation calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (starFormationRateSurfaceDensityDisksExtendedSchmidt), intent(inout) :: self
    type   (treeNode                                           ), intent(inout) :: node
    integer(kind_int8                                          ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node
    
    self%factorsComputed=.false.
    self%lastUniqueID   =uniqueID
    return
  end subroutine extendedSchmidtCalculationReset

  double precision function extendedSchmidtRate(self,node,radius)
    !!{RST
    Returns the star formation rate surface density (in :math:`\mathrm{M}_\odot` Gyr\ :math:`^{-1}` Mpc\ :math:`^{-2}`) for star formation in the galactic disk of ``node``. The disk is assumed to obey the extended Schmidt law of :cite:t:`shi_extended_2011`:

    .. math::

       \dot{\Sigma}_\star = A \left(x_\mathrm{H} {\Sigma_\mathrm{gas}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)
       ^{N_1} \left({\Sigma_{\star}\over \mathrm{M}_\odot \hbox{pc}^{-2}}\right)^{N_2},

    where :math:`A=`\ ``[normalization]``, :math:`N_1=`\ ``[exponentGas]``, and :math:`N_2=`\ ``[exponentStars]``.
    !!}
    use :: Abundances_Structure      , only : abundances
    use :: Coordinates               , only : coordinateCylindrical, assignment(=)
    use :: Galactic_Structure_Options, only : componentTypeDisk    , massTypeGaseous, massTypeStellar
    use :: Galacticus_Nodes          , only : nodeComponentDisk
    use :: Mass_Distributions        , only : massDistributionClass
    implicit none
    class           (starFormationRateSurfaceDensityDisksExtendedSchmidt), intent(inout) :: self
    type            (treeNode                                           ), intent(inout) :: node
    double precision                                                     , intent(in   ) :: radius
    class           (nodeComponentDisk                                  ), pointer       :: disk
    class           (massDistributionClass                              ), pointer       :: massDistributionGaseous, massDistributionStellar
    type            (abundances                                         ), save          :: abundancesFuel
    !$omp threadprivate(abundancesFuel)
    double precision                                                                     :: massGas                , surfaceDensityGas, &
         &                                                                                  surfaceDensityStellar
    type            (coordinateCylindrical                              )                :: coordinates

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if factors have been precomputed.
    if (.not.self%factorsComputed) then
       ! Get the disk properties.
       disk    => node%disk   ()
       massGas =  disk%massGas()
       ! Find the hydrogen fraction in the disk gas of the fuel supply.
       abundancesFuel=disk%abundancesGas()
       call abundancesFuel%massToMassFraction(massGas)
       self%hydrogenMassFraction=abundancesFuel%hydrogenMassFraction()
       ! Record that factors have now been computed.
       self%factorsComputed=.true.
    end if
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
    extendedSchmidtRate=+  self%normalization                                                    & ! Normalization of the star formation rate.
         &              *((self%hydrogenMassFraction*surfaceDensityGas    )**self%exponentGas  ) &
         &              *(                           surfaceDensityStellar **self%exponentStars)
    return
  end function extendedSchmidtRate
