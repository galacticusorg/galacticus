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
  Implements a dust properties class with a universal dust-to-metals ratio and a fixed opacity.
  !!}

  !![
  <dustProperties name="dustPropertiesSimple" docformat="rst">
   <description>
   Dust properties with a universal dust-to-metals ratio, :math:`f_\mathrm{dust:metals}` (``dustToMetalsRatio``), and a
   fixed :math:`V`-band extinction opacity per unit mass of dust, :math:`\kappa_\mathrm{V}` (``opacityExtinctionV``).

   The default ratio, :math:`f_\mathrm{dust:metals}=0.44`, is approximately correct for the Milky Way
   (:cite:t:`popping_dust_2017`). The default opacity is the one for which dust with that ratio reproduces the relation
   between column density and reddening of :cite:t:`savage_observed_1979`,

   .. math::

      \kappa_\mathrm{V} f_\mathrm{dust:metals} = \frac{A_\mathrm{V}/E(B-V)}{N_\mathrm{H}/E(B-V)} \, \frac{X_\odot}{m_\mathrm{u} Z_\mathrm{ISM}} \, \frac{1}{2.5 \log_{10} \mathrm{e}} \approx 1.048\times 10^4\,\hbox{cm}^2\,\hbox{g}^{-1},

   with :math:`A_\mathrm{V}/E(B-V)=3.1`, :math:`N_\mathrm{H}/E(B-V)=5.8\times10^{21}\,\hbox{atoms cm}^{-2}\,\hbox{mag}^{-1}`,
   and :math:`Z_\mathrm{ISM}=0.02` the metallicity of the local interstellar medium, giving
   :math:`\kappa_\mathrm{V} \approx 2.38\times 10^4\,\hbox{cm}^2\,\hbox{g}^{-1}`. With both at their defaults, dust
   attenuation models taking their normalization from this class reproduce that Milky Way calibration.

   The two parameters are physically distinct: the ratio sets how much dust a galaxy has, and the opacity how strongly
   that dust extinguishes light. Changing the ratio alone changes the mass of dust, and every optical depth derived
   from it, in proportion; the opacity is not adjusted to compensate.
   </description>
  </dustProperties>
  !!]
  type, extends(dustPropertiesClass) :: dustPropertiesSimple
     !!{RST
     A dust properties class with a universal dust-to-metals ratio and a fixed opacity.
     !!}
     private
     double precision :: dustToMetalsRatio_, opacityExtinctionV_
   contains
     procedure :: dustToMetalsRatio  => simpleDustToMetalsRatio
     procedure :: opacityExtinctionV => simpleOpacityExtinctionV
  end type dustPropertiesSimple

  interface dustPropertiesSimple
     !!{RST
     Constructors for the :galacticus-class:`dustPropertiesSimple` dust properties class.
     !!}
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface dustPropertiesSimple

contains

  function simpleConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`dustPropertiesSimple` dust properties class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (dustPropertiesSimple)                :: self
    type            (inputParameters     ), intent(inout) :: parameters
    double precision                                      :: dustToMetalsRatio_, opacityExtinctionV_

    !![
    <inputParameter docformat="rst">
      <name>dustToMetalsRatio</name>
      <variable>dustToMetalsRatio_</variable>
      <defaultValue>dustToMetalsRatioMilkyWay</defaultValue>
      <defaultSource>Approximately correct for the Milky Way (e.g. :cite:t:`popping_dust_2017`).</defaultSource>
      <description>
      The fraction of the mass of metals which is in dust.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>opacityExtinctionV</name>
      <variable>opacityExtinctionV_</variable>
      <defaultValue>depthOpticalVPerSurfaceDensityMetalsMilkyWay/dustToMetalsRatioMilkyWay</defaultValue>
      <defaultSource>Chosen so that, with the default ``dustToMetalsRatio``, the relation between column density and reddening of :cite:t:`savage_observed_1979` is reproduced.</defaultSource>
      <description>
      The :math:`V`-band extinction opacity per unit mass of dust, in cm² g⁻¹.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=dustPropertiesSimple(dustToMetalsRatio_,opacityExtinctionV_)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(dustToMetalsRatio_,opacityExtinctionV_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`dustPropertiesSimple` dust properties class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (dustPropertiesSimple)                :: self
    double precision                      , intent(in   ) :: dustToMetalsRatio_, opacityExtinctionV_
    !![
    <constructorAssign variables="dustToMetalsRatio_, opacityExtinctionV_"/>
    !!]

    if (self%dustToMetalsRatio_  < 0.0d0 .or. self%dustToMetalsRatio_ > 1.0d0) &
         & call Error_Report('`dustToMetalsRatio` must lie between zero and one'//{introspection:location})
    if (self%opacityExtinctionV_ < 0.0d0                                     ) &
         & call Error_Report('`opacityExtinctionV` must be non-negative'         //{introspection:location})
    return
  end function simpleConstructorInternal

  double precision function simpleDustToMetalsRatio(self,node,componentType) result(ratio)
    !!{RST
    Return the dust-to-metals ratio, which is the same for every component of every galaxy.
    !!}
    implicit none
    class(dustPropertiesSimple        ), intent(inout)         :: self
    type (treeNode                    ), intent(inout), target :: node
    type (enumerationComponentTypeType), intent(in   )         :: componentType
    !$GLC attributes unused :: node, componentType

    ratio=self%dustToMetalsRatio_
    return
  end function simpleDustToMetalsRatio

  double precision function simpleOpacityExtinctionV(self) result(opacity)
    !!{RST
    Return the :math:`V`-band extinction opacity per unit mass of dust, in cm² g⁻¹.
    !!}
    implicit none
    class(dustPropertiesSimple), intent(inout) :: self

    opacity=self%opacityExtinctionV_
    return
  end function simpleOpacityExtinctionV
