!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !!{
  Implementation of a satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.
  !!}

  use :: Satellite_Tidal_Stripping_Radii, only : satelliteTidalStrippingRadiusClass
  use :: Galactic_Structure             , only : galacticStructureClass

  !![
  <satelliteTidalStripping name="satelliteTidalStrippingZentner2005">
   <description>
    A satellite tidal stripping class which uses the formalism of \cite{zentner_physics_2005} to compute the mass loss rate
    $\dot{M}_\mathrm{sat}$:
    \begin{equation}
    \dot{M}_\mathrm{sat}=-\alpha \frac{M_\mathrm{sat}(>r_\mathrm{tidal})}{T_\mathrm{orb}},
    \end{equation}
    where $\alpha=${\normalfont \ttfamily [efficiency]},
    \begin{equation}
    T_\mathrm{orb} = {1 \over \hbox{max}(\omega/2\pi,v_\mathrm{r}/r)},
    \end{equation}
    where $\omega$ is the angular velocity of the satellite, $v_\mathrm{r}$ is the radial velocity, $r$ is the orbital radius,
    and $r_\mathrm{tidal}$ is the tidal radius of the satellite, given by the \cite{king_structure_1962} formula:
    \begin{equation}
    r_\mathrm{tidal}=\left(\frac{GM_\mathrm{sat}}{\omega^2-d^2\Phi/dr^2}\right)^{1/3},
    \end{equation}
    where $\omega$ is the orbital angular velocity of the satellite, and $\Phi(r)$ is the gravitational potential due to the
    host.
   </description>
  </satelliteTidalStripping>
  !!]
  type, extends(satelliteTidalStrippingClass) :: satelliteTidalStrippingZentner2005
     !!{
     Implementation of a satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.
     !!}
     private
     class           (satelliteTidalStrippingRadiusClass), pointer :: satelliteTidalStrippingRadius_ => null()
     class           (galacticStructureClass            ), pointer :: galacticStructure_             => null()
     double precision                                              :: efficiency
   contains
     final     ::                 zentner2005Destructor
     procedure :: massLossRate => zentner2005MassLossRate
  end type satelliteTidalStrippingZentner2005

  interface satelliteTidalStrippingZentner2005
     !!{
     Constructors for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
     !!}
     module procedure zentner2005ConstructorParameters
     module procedure zentner2005ConstructorInternal
  end interface satelliteTidalStrippingZentner2005

contains

  function zentner2005ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (satelliteTidalStrippingZentner2005)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (satelliteTidalStrippingRadiusClass), pointer       :: satelliteTidalStrippingRadius_
    class           (galacticStructureClass            ), pointer       :: galacticStructure_
    double precision                                                    :: efficiency

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultValue>2.5d0</defaultValue>
      <description>The dimensionless rate coefficient appearing in the \cite{zentner_physics_2005} expression for the tidal mass loss rate from subhalos.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteTidalStrippingRadius" name="satelliteTidalStrippingRadius_" source="parameters"/>
    <objectBuilder class="galacticStructure"             name="galacticStructure_"             source="parameters"/>
    !!]
    self=satelliteTidalStrippingZentner2005(efficiency,satelliteTidalStrippingRadius_,galacticStructure_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStrippingRadius_"/>
    <objectDestructor name="galacticStructure_"            />
    !!]
    return
  end function zentner2005ConstructorParameters

  function zentner2005ConstructorInternal(efficiency,satelliteTidalStrippingRadius_,galacticStructure_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
    !!}
    implicit none
    type            (satelliteTidalStrippingZentner2005)                        :: self
    class           (satelliteTidalStrippingRadiusClass), intent(in   ), target :: satelliteTidalStrippingRadius_
    class           (galacticStructureClass            ), intent(in   ), target :: galacticStructure_
    double precision                                    , intent(in)            :: efficiency
    !![
    <constructorAssign variables="efficiency, *satelliteTidalStrippingRadius_, *galacticStructure_"/>
    !!]

    return
  end function zentner2005ConstructorInternal

  subroutine zentner2005Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
    !!}
    implicit none
    type(satelliteTidalStrippingZentner2005), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStrippingRadius_"/>
    <objectDestructor name="self%galacticStructure_"            />
    !!]
    return
  end subroutine zentner2005Destructor

  double precision function zentner2005MassLossRate(self,node)
    !!{
    Return a mass loss rate for satellites due to tidal stripping using the formulation of \cite{zentner_physics_2005}.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentSatellite, treeNode
    use :: Numerical_Constants_Astronomical, only : gigaYear              , megaParsec    , gravitationalConstantGalacticus
    use :: Numerical_Constants_Math        , only : Pi
    use :: Numerical_Constants_Prefixes    , only : kilo
    use :: Vectors                         , only : Vector_Magnitude      , Vector_Product
    implicit none
    class           (satelliteTidalStrippingZentner2005), intent(inout)  :: self
    type            (treeNode                          ), intent(inout)  :: node
    class           (nodeComponentSatellite            ), pointer        :: satellite
    double precision                                    , dimension(3  ) :: position          , velocity
    double precision                                                     :: massSatellite     , frequencyAngular, &
         &                                                                  periodOrbital     , radius          , &
         &                                                                  massOuterSatellite, frequencyRadial

    ! Get required quantities from the satellite node.
    satellite          =>  node     %satellite (        )
    massSatellite      =   satellite%boundMass (        )
    position           =   satellite%position  (        )
    velocity           =   satellite%velocity  (        )
    radius             =   Vector_Magnitude    (position)
    ! Compute the orbital period.
    frequencyAngular   =  +Vector_Magnitude(Vector_Product(position,velocity)) &
         &                /radius**2                                           &
         &                *kilo                                                &
         &                *gigaYear                                            &
         &                /megaParsec
    frequencyRadial    =  +abs             (   Dot_Product(position,velocity)) &
         &                /radius**2                                           &
         &                *kilo                                                &
         &                *gigaYear                                            &
         &                /megaParsec
    ! Find the orbital period. We use the larger of the angular and radial frequencies to avoid numerical problems for purely
    ! radial or purely circular orbits.
    periodOrbital      =  +2.0d0                 &
         &                *Pi                    &
         &                /max(                  &
         &                     frequencyAngular, &
         &                     frequencyRadial   &
         &                    )
    ! Compute the mass of the satellite outside of the tidal radius.
    massOuterSatellite =  max(                                                                                              &
         &                    +massSatellite                                                                                &
         &                    -self%galacticStructure_%massEnclosed(node,self%satelliteTidalStrippingRadius_%radius(node)), &
         &                    +0.0d0                                                                                        &
         &                   )
    ! Compute the rate of mass loss.
    zentner2005MassLossRate=-self%efficiency    &
         &                  *massOuterSatellite &
         &                  /periodOrbital
    return
  end function zentner2005MassLossRate
