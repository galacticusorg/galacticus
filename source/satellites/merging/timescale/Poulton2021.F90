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

  !+    Contributions to this file made by: Claude.

  !!{RST
  Implements calculations of satellite merging times using the :cite:t:`poulton_extracting_2021` method.
  !!}

  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass

  !![
  <satelliteMergingTimescales name="satelliteMergingTimescalesPoulton2021" docformat="rst">
   <description>
   A satellite merging timescales class which computes merging timescales using the fitting formula of :cite:t:`poulton_extracting_2021`. The merging timescale is:

   .. math::

      T_\mathrm{merge} = A \, T_\mathrm{dyn}(r) \left(\frac{r}{R_\mathrm{vir}}\right)^b \left(\frac{R_\mathrm{peri}}{R_\mathrm{vir}}\right)^c,

   where :math:`T_\mathrm{dyn}(r) = r/v_\mathrm{circ}(r)` is the local dynamical time, :math:`R_\mathrm{peri}` is the pericentric distance computed using the locally Keplerian approximation with the enclosed host mass (for unbound as well as bound orbits, so that unbound satellites also receive a finite merging time), and

   .. math::

      b = \begin{cases} b_\mathrm{interior} &amp; r &lt; R_\mathrm{vir} \\ b_\mathrm{exterior} &amp; r \geq R_\mathrm{vir} \end{cases}.

   Parameters are set via :math:`A=`\ ``A``, :math:`c=`\ ``[c]``, :math:`b_\mathrm{interior}=`\ ``[bInterior]``, and :math:`b_\mathrm{exterior}=`\ ``[bExterior]``. The default parameter values :math:`A=5.5`, :math:`c=0.2`, :math:`b_\mathrm{interior}=-0.5`, and :math:`b_\mathrm{exterior}=-1.0` are those favored by :cite:t:`poulton_extracting_2021`.
   </description>
  </satelliteMergingTimescales>
  !!]
  type, extends(satelliteMergingTimescalesClass) :: satelliteMergingTimescalesPoulton2021
     !!{RST
     A class implementing the :cite:t:`poulton_extracting_2021` method for satellite merging timescales.
     !!}
     private
     class           (darkMatterHaloScaleClass ), pointer :: darkMatterHaloScale_  => null()
     class           (darkMatterProfileDMOClass), pointer :: darkMatterProfileDMO_ => null()
     double precision                                     :: A                              , c        , &
          &                                                  bInterior                      , bExterior
   contains
     final     ::                     poulton2021Destructor
     procedure :: timeUntilMerging => poulton2021TimeUntilMerging
  end type satelliteMergingTimescalesPoulton2021

  interface satelliteMergingTimescalesPoulton2021
     !!{RST
     Constructors for the :galacticus-class:`satelliteMergingTimescalesPoulton2021` satellite merging timescale class.
     !!}
     module procedure poulton2021ConstructorParameters
     module procedure poulton2021ConstructorInternal
  end interface satelliteMergingTimescalesPoulton2021

contains

  function poulton2021ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :cite:t:`poulton_extracting_2021` merging timescale class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteMergingTimescalesPoulton2021)                :: self
    type            (inputParameters                      ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass             ), pointer       :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass            ), pointer       :: darkMatterProfileDMO_
    double precision                                                       :: A                    , c        , &
         &                                                                    bInterior            , bExterior

    !![
    <componentPropertyAssert class="basic" properties="mass" require="gettable"/>
    <inputParameter docformat="rst">
      <name>A</name>
      <defaultValue>5.5d0</defaultValue>
      <description>
      The normalization constant, :math:`A`, in the :cite:t:`poulton_extracting_2021` satellite merging timescale fitting formula.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>c</name>
      <defaultValue>0.2d0</defaultValue>
      <description>
      The exponent, :math:`c`, on the pericentric distance ratio in the :cite:t:`poulton_extracting_2021` satellite merging timescale fitting formula.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>bInterior</name>
      <defaultValue>-0.5d0</defaultValue>
      <description>
      The radial exponent, :math:`b`, in the :cite:t:`poulton_extracting_2021` satellite merging timescale fitting formula for satellites inside the virial radius (:math:`r &lt; R_\mathrm{vir}`).
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>bExterior</name>
      <defaultValue>-1.0d0</defaultValue>
      <description>
      The radial exponent, :math:`b`, in the :cite:t:`poulton_extracting_2021` satellite merging timescale fitting formula for satellites outside the virial radius (:math:`r \geq R_\mathrm{vir}`).
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteMergingTimescalesPoulton2021(A,c,bInterior,bExterior,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function poulton2021ConstructorParameters

  function poulton2021ConstructorInternal(A,c,bInterior,bExterior,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{RST
    Constructor for the :cite:t:`poulton_extracting_2021` merging timescale class.
    !!}
    implicit none
    type            (satelliteMergingTimescalesPoulton2021)                        :: self
    double precision                                       , intent(in   )         :: A                    , c        , &
         &                                                                            bInterior            , bExterior
    class           (darkMatterHaloScaleClass             ), intent(in   ), target :: darkMatterHaloScale_
    class           (darkMatterProfileDMOClass            ), intent(in   ), target :: darkMatterProfileDMO_
    !![
    <constructorAssign variables="A, c, bInterior, bExterior, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function poulton2021ConstructorInternal

  subroutine poulton2021Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteMergingTimescalesPoulton2021` satellite merging timescale class.
    !!}
    implicit none
    type(satelliteMergingTimescalesPoulton2021), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine poulton2021Destructor

  double precision function poulton2021TimeUntilMerging(self,node,orbit) result(timeUntilMerging)
    !!{RST
    Return the timescale for merging satellites using the :cite:t:`poulton_extracting_2021` method. All orbital quantities
    (energy, angular momentum, eccentricity, pericentric distance, and dynamical time) are computed, as in
    :cite:t:`poulton_extracting_2021` (their eqns. 2--7), from the host mass enclosed within the current orbital radius,
    :math:`M_\mathrm{encl}(r)`, which is taken to be the virial mass for :math:`r \ge R_\mathrm{vir}`, and the reduced mass
    :math:`\mu=M_\mathrm{sat}M_\mathrm{encl}/(M_\mathrm{sat}+M_\mathrm{encl})`. Unbound orbits (:math:`e \ge 1`) receive the
    same fitted, finite, merging time.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , treeNode
    use :: Kepler_Orbits                   , only : keplerOrbit
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal, MpcPerKmPerSToGyr
    implicit none
    class           (satelliteMergingTimescalesPoulton2021), intent(inout) :: self
    type            (treeNode                             ), intent(inout) :: node
    type            (keplerOrbit                          ), intent(inout) :: orbit
    type            (treeNode                             ), pointer       :: nodeHost
    class           (nodeComponentBasic                   ), pointer       :: basic
    class           (massDistributionClass                ), pointer       :: massDistribution_
    double precision                                                       :: radiusOrbital    , radiusVirial    , &
         &                                                                    massEnclosed     , timeDynamical   , &
         &                                                                    massReduced      , radiusPericenter, &
         &                                                                    massSatellite    , bExponent       , &
         &                                                                    energyOrbital    , angularMomentum , &
         &                                                                    eccentricity

    ! Find the host node.
    if (node%isSatellite()) then
       nodeHost => node%parent
    else
       nodeHost => node%parent%firstChild
    end if
    ! Get the orbital radius and host virial radius.
    radiusOrbital=orbit                     %radius      (        )
    radiusVirial =self %darkMatterHaloScale_%radiusVirial(nodeHost)
    ! Get the host mass enclosed within the orbital radius (Poulton et al. 2021, eqn. 2), which is the virial mass outside the
    ! virial radius.
    massDistribution_ => self             %darkMatterProfileDMO_%get                 (nodeHost                       )
    massEnclosed      =  massDistribution_                      %massEnclosedBySphere(min(radiusOrbital,radiusVirial))
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Return instantaneous merging if the enclosed mass vanishes (e.g. r → 0).
    if (massEnclosed <= 0.0d0) then
       timeUntilMerging=0.0d0
       return
    end if
    ! Compute the orbital energy, angular momentum, and eccentricity (Poulton et al. 2021, eqns. 3-5), using the reduced mass of
    ! the satellite and the enclosed host mass.
    basic            =>  node%basic()
    massSatellite    =   basic%mass()
    massReduced      =  +massSatellite   &
         &              *massEnclosed    &
         &              /(               &
         &                +massSatellite &
         &                +massEnclosed  &
         &               )
    energyOrbital    =  +0.5d0                           &
         &              *massReduced                     &
         &              *(                               &
         &                +orbit%velocityRadial    ()**2 &
         &                +orbit%velocityTangential()**2 &
         &               )                               &
         &              -gravitationalConstant_internal  &
         &              *massEnclosed*massSatellite      &
         &              /radiusOrbital
    angularMomentum  =  +massReduced                &
         &              *radiusOrbital              &
         &              *orbit%velocityTangential()
    eccentricity     =  +sqrt(                                                                                                                        &
         &                    max(                                                                                                                    &
         &                        +0.0d0                                                                                                            , &
         &                        +1.0d0                                                                                                              &
         &                        +2.0d0*energyOrbital*angularMomentum**2/(gravitationalConstant_internal*massEnclosed*massSatellite)**2/massReduced  &
         &                       )                                                                                                                    &
         &                   )
    ! Pericentric radius via the locally-Keplerian approximation (Poulton et al. 2021, eqn. 6). This holds for unbound
    ! (e ≥ 1) as well as bound orbits, so the fitting formula is applied to both.
    radiusPericenter =  +angularMomentum**2             &
         &              /(1.0d0+eccentricity)           &
         &              /gravitationalConstant_internal &
         &              /massEnclosed                   &
         &              /massSatellite                  &
         &              /massReduced
    ! Local dynamical time (Poulton et al. 2021, eqn. 7), converted from Mpc/(km/s) to Gyr.
    timeDynamical    =  +sqrt(                                &
         &                    +radiusOrbital**3               &
         &                    /gravitationalConstant_internal &
         &                    /massEnclosed                   &
         &                   )                                &
         &              *MpcPerKmPerSToGyr
    ! Select radial exponent b depending on whether satellite is inside or outside R_vir.
    if (radiusOrbital < radiusVirial) then
       bExponent=self%bInterior
    else
       bExponent=self%bExterior
    end if
    ! Evaluate the Poulton et al. (2021) fitting formula.
    timeUntilMerging=+self%A                                          &
         &           *     timeDynamical                              &
         &           *(radiusOrbital   /radiusVirial)**     bExponent &
         &           *(radiusPericenter/radiusVirial)**self%c
    return
  end function poulton2021TimeUntilMerging
