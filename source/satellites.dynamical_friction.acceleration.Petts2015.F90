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
  Implementation of a satellite dynamical friction class which uses the model of :cite:t:`petts_semi-analytic_2015` for the Coulomb logarithm.
  !!}

  use :: Cosmology_Parameters, only : cosmologyParametersClass

  !![
  <satelliteDynamicalFriction name="satelliteDynamicalFrictionPetts2015" docformat="rst">
   <description>
   A satellite dynamical friction class which computes the Coulomb logarithm following the model of :cite:t:`petts_semi-analytic_2015`. The minimum impact parameter is taken to be :math:`b_\mathrm{min} = \max(r_{1/2,\mathrm{sat}}, G M_\mathrm{sat} / v_\mathrm{orbital}^2)`, where :math:`r_{1/2,\mathrm{sat}}` is the half-mass radius of the dark matter component of the satellite, and the maximum is :math:`b_\mathrm{max} = r_\mathrm{orbital}/|\gamma|` for :math:`|\gamma| &gt; 1` (where :math:`\gamma` is the logarithmic density slope of the host at the satellite position), or :math:`r_\mathrm{orbital}` otherwise. The ``[logarithmCoulombApproximate]`` parameter controls whether the Coulomb logarithm is evaluated as :math:`\ln\Lambda` (for :math:`\Lambda \geq 1`, else zero) or the approximate form :math:`\frac{1}{2}\ln(1+\Lambda^2)`.
   </description>
  </satelliteDynamicalFriction>
  !!]
  type, extends(satelliteDynamicalFrictionChandrasekhar1943) :: satelliteDynamicalFrictionPetts2015
     !!{RST
     Implementation of a satellite dynamical friction class which uses the model of :cite:t:`petts_semi-analytic_2015` for the Coulomb logarithm.
     !!}
     private
     class  (cosmologyParametersClass), pointer :: cosmologyParameters_        => null()
     logical                                    :: logarithmCoulombApproximate
   contains
     final     ::                     petts2015Destructor
     procedure :: coulombLogarithm => petts2015CoulombLogarithm
  end type satelliteDynamicalFrictionPetts2015

  interface satelliteDynamicalFrictionPetts2015
     !!{RST
     Constructors for the petts2015 satellite dynamical friction class.
     !!}
     module procedure petts2015ConstructorParameters
     module procedure petts2015ConstructorInternal
  end interface satelliteDynamicalFrictionPetts2015

contains

  function petts2015ConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`satelliteDynamicalFrictionPetts2015` satellite dynamical friction class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (satelliteDynamicalFrictionPetts2015)                :: self
    type   (inputParameters                    ), intent(inout) :: parameters
    class  (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass          ), pointer       :: darkMatterProfileDMO_
    class  (cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    logical                                                     :: logarithmCoulombApproximate
    
    !![
    <inputParameter docformat="rst">
      <name>logarithmCoulombApproximate</name>
      <defaultValue>.false.</defaultValue>
      <description>
      If true, the Coulomb logarithm term is evaluated as :math:`\frac{1}{2}\log(1+\Lambda^2)`. Otherwise, it is evaluated as :math:`\log\Lambda` for :math:`\Lambda \ge 1` and set to zero otherwise.
      </description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="cosmologyParameters"  name="cosmologyParameters_"  source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=satelliteDynamicalFrictionPetts2015(logarithmCoulombApproximate,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    <objectDestructor name="cosmologyParameters_" />
    !!]
    return
  end function petts2015ConstructorParameters

  function petts2015ConstructorInternal(logarithmCoulombApproximate,cosmologyParameters_,darkMatterHaloScale_,darkMatterProfileDMO_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`satelliteDynamicalFrictionPetts2015` satellite dynamical friction class.
    !!}
    implicit none
    type   (satelliteDynamicalFrictionPetts2015)                        :: self
    class  (cosmologyParametersClass           ), intent(in   ), target :: cosmologyParameters_
    class  (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass          ), intent(in   ), target :: darkMatterProfileDMO_
    logical                                     , intent(in   )         :: logarithmCoulombApproximate
    !![
    <constructorAssign variables="logarithmCoulombApproximate, *cosmologyParameters_, *darkMatterHaloScale_, *darkMatterProfileDMO_"/>
    !!]

    return
  end function petts2015ConstructorInternal

  subroutine petts2015Destructor(self)
    !!{RST
    Destructor for the :galacticus-class:`satelliteDynamicalFrictionPetts2015` satellite dynamical friction class.
    !!}
    implicit none
    type(satelliteDynamicalFrictionPetts2015), intent(inout) :: self

    !![
    <objectDestructor name="self%cosmologyParameters_" />
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    return
  end subroutine petts2015Destructor
  
  double precision function petts2015CoulombLogarithm(self,node) result(coulombLogarithm)
    !!{RST
    Evaluate the Coulomb logarithm for the :cite:t:`petts_semi-analytic_2015` dynamical friction model.
    !!}
    use :: Coordinates                     , only : coordinateSpherical           , assignment(=)
    use :: Galacticus_Nodes                , only : nodeComponentBasic            , nodeComponentSatellite, treeNode
    use :: Galactic_Structure_Options      , only : componentTypeAll              , massTypeDark
    use :: Mass_Distributions              , only : massDistributionClass
    use :: Numerical_Constants_Astronomical, only : gravitationalConstant_internal
    use :: Vectors                         , only : Vector_Magnitude
    implicit none
    class           (satelliteDynamicalFrictionPetts2015), intent(inout) :: self
    type            (treeNode                           ), intent(inout) :: node
    class           (nodeComponentSatellite             ), pointer       :: satellite
    class           (nodeComponentBasic                 ), pointer       :: basic
    type            (treeNode                           ), pointer       :: nodeHost
    class           (massDistributionClass              ), pointer       :: massDistribution_
    double precision                                     , dimension(3)  :: position               , velocity    
    double precision                                                     :: speedOrbital           , radiusOrbital          , &
         &                                                                  impactParameterMinimum , impactParameterMaximum , &
         &                                                                  massSatellite          , densitySlopeLogarithmic, &
         &                                                                  radiusHalfMassSatellite, fractionDarkMatter     , &
         &                                                                  massHalfSatellite
    type            (coordinateSpherical                )                :: coordinates

    nodeHost                =>  node     %mergesWith(        )
    satellite               =>  node     %satellite (        )
    basic                   =>  node     %basic     (        )
    massSatellite           =   satellite%boundMass (        )
    position                =   satellite%position  (        )
    velocity                =   satellite%velocity  (        )
    radiusOrbital           =   Vector_Magnitude    (position )
    speedOrbital            =   Vector_Magnitude    (velocity )
    fractionDarkMatter      =  +(                                          &
         &                        +self%cosmologyParameters_%OmegaMatter() &
         &                        -self%cosmologyParameters_%OmegaBaryon() &
         &                       )                                         &
         &                      /  self%cosmologyParameters_%OmegaMatter()
    massHalfSatellite       =  +0.50d0             &
         &                     *fractionDarkMatter &
         &                     *min(               &
         &                          massSatellite, &
         &                          basic%mass()   &
         &                         )
    massDistribution_       =>  node             %massDistribution   (componentType=componentTypeAll ,massType=massTypeDark)
    radiusHalfMassSatellite =   massDistribution_%radiusEnclosingMass(mass         =massHalfSatellite                      )
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    coordinates             =  [radiusOrbital,0.0d0,0.0d0]
    massDistribution_       =>     self             %darkMatterProfileDMO_%get                  (nodeHost                      )
    densitySlopeLogarithmic =  abs(massDistribution_                      %densityGradientRadial(coordinates,logarithmic=.true.))
    !![
    <objectDestructor name="massDistribution_"/>
    !!]
    ! Evaluate the minimum and maximum impact parameters.
    if (densitySlopeLogarithmic > 1.0d0) then
       impactParameterMaximum=radiusOrbital/densitySlopeLogarithmic
    else
       impactParameterMaximum=radiusOrbital
    end if
    impactParameterMinimum=max(radiusHalfMassSatellite,gravitationalConstant_internal*massSatellite/speedOrbital**2)
    ! Evaluate the Coulomb logarithm, using either the approximate or full expression.
    if (impactParameterMinimum <= 0.0d0) then
       coulombLogarithm=0.0d0
    else
       if (self%logarithmCoulombApproximate) then
          coulombLogarithm=      log(max(1.0d0, impactParameterMaximum/impactParameterMinimum)   )
       else
          coulombLogarithm=0.5d0*log(    1.0d0+(impactParameterMaximum/impactParameterMinimum)**2)
       end if
    end if
    return
  end function petts2015CoulombLogarithm
