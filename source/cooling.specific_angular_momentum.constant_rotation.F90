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

  !!{
  Implementation of a specific angular momentum of cooling gas class assuming a constant rotation velocity as a function of
  radius.
  !!}

  use :: Dark_Matter_Profiles_DMO, only : darkMatterProfileDMOClass
  use :: Dark_Matter_Halo_Scales , only : darkMatterHaloScaleClass
  use :: Kind_Numbers            , only : kind_int8

  ! Enumeration for angular momentum source.
  !![
  <enumeration>
   <name>angularMomentumSource</name>
   <description>Enumeration specifying the origin of angular momentum of cooling gas in the constant rotation class.</description>
   <encodeFunction>yes</encodeFunction>
   <entry label="darkMatter"/>
   <entry label="hotGas"    />
  </enumeration>
  !!]

  !![
  <coolingSpecificAngularMomentum name="coolingSpecificAngularMomentumConstantRotation">
   <description>
    A cooling specific angular momentum class which assumes a constant rotation velocity as a function of radius. The specific
    angular momentum of cooling gas is given either by
    \begin{equation}
     j_\mathrm{cool} = \langle j \rangle r_\mathrm{cool} A,
    \end{equation}
    where $r_\mathrm{cool}$ is the cooling radius, $A$ is the rotation normalization and $\langle j \rangle$ is the mean
    specific angular momentum of the cooling gas, if {\normalfont \ttfamily [useInteriorMean]}$=${\normalfont \ttfamily false},
    or by
    \begin{equation}
     j_\mathrm{cool} = \langle j \rangle {I_3(r_\mathrm{cool})/I_2(r_\mathrm{cool})} A,
    \end{equation}
    where $I_n(r)$ is the $n^\mathrm{th}$ radial moment of the hot gas density profile from $0$ to $r$ (this therefore gives
    the mean specific angular momentum interior to radius $r$), if {\normalfont \ttfamily [useInteriorMean]}$=${\normalfont
    \ttfamily true}.
  
    If {\normalfont \ttfamily [sourceAngularMomentumSpecificMean]}$=${\normalfont \ttfamily darkMatter} then $\langle j \rangle$
    is the mean specific angular momentum of the dark matter halo, while if {\normalfont \ttfamily
    [sourceAngularMomentumSpecificMean]}$=${\normalfont \ttfamily hotGas} then $\langle j \rangle$ is equal to the mean specific
    angular momentum of gas currently in the hot gas reservoir. If {\normalfont \ttfamily
    [sourceNormalizationRotation]}$=${\normalfont \ttfamily darkMatter} then the rotation normalization $A$ is computed using the
    dark matter density profile, while if {\normalfont \ttfamily [sourceNormalizationRotation]}$=${\normalfont \ttfamily hotGas}
    it is computed using the density profile of the hot gas reservoir.
   </description>
  </coolingSpecificAngularMomentum>
  !!]
  type, extends(coolingSpecificAngularMomentumClass) :: coolingSpecificAngularMomentumConstantRotation
     !!{
     Implementation of the specific angular momentum of cooling gas class which assumes a constant rotation velocity as a function of radius.
     !!}
     private
     class           (darkMatterHaloScaleClass            ), pointer :: darkMatterHaloScale_              => null()
     class           (darkMatterProfileDMOClass           ), pointer :: darkMatterProfileDMO_             => null()
     integer         (kind=kind_int8                      )          :: lastUniqueID
     logical                                                         :: angularMomentumSpecificComputed
     double precision                                                :: angularMomentumSpecificPrevious
     type            (enumerationAngularMomentumSourceType)          :: sourceAngularMomentumSpecificMean          , sourceNormalizationRotation
     logical                                                         :: useInteriorMean
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                            constantRotationDestructor
     procedure :: autoHook                => constantRotationAutoHook
     procedure :: calculationReset        => constantRotationCalculationReset
     procedure :: angularMomentumSpecific => constantRotationAngularMomentumSpecific
  end type coolingSpecificAngularMomentumConstantRotation

  interface coolingSpecificAngularMomentumConstantRotation
     !!{
     Constructors for the constantRotation specific angular momentum of cooling gas class.
     !!}
     module procedure constantRotationConstructorParameters
     module procedure constantRotationConstructorInternal
  end interface coolingSpecificAngularMomentumConstantRotation

contains

  function constantRotationConstructorParameters(parameters) result(self)
    !!{
    Constructor for the constantRotation freefall radius class which builds the object from a parameter set.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type   (coolingSpecificAngularMomentumConstantRotation)                :: self
    type   (inputParameters                               ), intent(inout) :: parameters
    class  (darkMatterProfileDMOClass                     ), pointer       :: darkMatterProfileDMO_
    class  (darkMatterHaloScaleClass                      ), pointer       :: darkMatterHaloScale_
    logical                                                                :: useInteriorMean
    type   (varying_string                                )                :: sourceAngularMomentumSpecificMean, sourceNormalizationRotation

    !![
    <inputParameter>
      <name>sourceAngularMomentumSpecificMean</name>
      <defaultValue>var_str('hotGas')</defaultValue>
      <description>
       The component (``{\normalfont \ttfamily hotGas}'' or ``{\normalfont \ttfamily darkMatter}'') from which the mean specific angular momentum should be computed for
       calculations of cooling gas specific angular momentum.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>sourceNormalizationRotation</name>
      <defaultValue>var_str('hotGas')</defaultValue>
      <description>
       The component (``{\normalfont \ttfamily hotGas}'' or ``{\normalfont \ttfamily darkMatter}'') from which the constant rotation speed should be computed for
       calculations of cooling gas specific angular momentum.
      </description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>useInteriorMean</name>
      <defaultValue>.false.</defaultValue>
      <description>Specifies whether to use the specific angular momentum at the cooling radius, or the mean specific angular momentum interior to that radius.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale"  name="darkMatterHaloScale_"  source="parameters"/>
    <objectBuilder class="darkMatterProfileDMO" name="darkMatterProfileDMO_" source="parameters"/>
    !!]
    self=coolingSpecificAngularMomentumConstantRotation(                                                                                                        &
         &                                              darkMatterHaloScale_                                                                                  , &
         &                                              darkMatterProfileDMO_                                                                                 , &
         &                                              enumerationAngularMomentumSourceEncode(char(sourceAngularMomentumSpecificMean),includesPrefix=.false.), &
         &                                              enumerationAngularMomentumSourceEncode(char(sourceNormalizationRotation      ),includesPrefix=.false.), &
         &                                              useInteriorMean                                                                                         &
         &                                             )
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_" />
    <objectDestructor name="darkMatterProfileDMO_"/>
    !!]
    return
  end function constantRotationConstructorParameters

  function constantRotationConstructorInternal(darkMatterHaloScale_,darkMatterProfileDMO_,sourceAngularMomentumSpecificMean,sourceNormalizationRotation,useInteriorMean) result(self)
    !!{
    Internal constructor for the darkMatterHalo freefall radius class.
    !!}
    implicit none
    type   (coolingSpecificAngularMomentumConstantRotation)                        :: self
    class  (darkMatterHaloScaleClass                      ), intent(in   ), target :: darkMatterHaloScale_
    class  (darkMatterProfileDMOClass                     ), intent(in   ), target :: darkMatterProfileDMO_
    type   (enumerationAngularMomentumSourceType          ), intent(in   )         :: sourceAngularMomentumSpecificMean, sourceNormalizationRotation
    logical                                                , intent(in   )         :: useInteriorMean
    !![
    <constructorAssign variables="*darkMatterHaloScale_, *darkMatterProfileDMO_, sourceAngularMomentumSpecificMean, sourceNormalizationRotation, useInteriorMean"/>
    !!]

    self%lastUniqueID                   =-1_kind_int8
    self%angularMomentumSpecificComputed=.false.
    return
  end function constantRotationConstructorInternal

  subroutine constantRotationAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(coolingSpecificAngularMomentumConstantRotation), intent(inout) :: self
    call calculationResetEvent%attach(self,constantRotationCalculationReset,openMPThreadBindingAllLevels,label='coolingSpecificAngularMomentumConstantRotation')
    return
  end subroutine constantRotationAutoHook

  subroutine constantRotationDestructor(self)
    !!{
    Destructor for the constant rotation specific angular momentum of cooling gas class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(coolingSpecificAngularMomentumConstantRotation), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_" />
    <objectDestructor name="self%darkMatterProfileDMO_"/>
    !!]
    if (calculationResetEvent%isAttached(self,constantRotationCalculationReset)) call calculationResetEvent%detach(self,constantRotationCalculationReset)
    return
  end subroutine constantRotationDestructor

  double precision function constantRotationAngularMomentumSpecific(self,node,radius)
    !!{
    Return the specific angular momentum of cooling gas in the constantRotation model.
    !!}
    use :: Error                     , only : Error_Report
    use :: Galacticus_Nodes          , only : nodeComponentBasic   , nodeComponentHotHalo, nodeComponentSpin, treeNode
    use :: Mass_Distributions        , only : massDistributionClass
    use :: Galactic_Structure_Options, only : componentTypeHotHalo , massTypeGaseous
    implicit none
    class           (coolingSpecificAngularMomentumConstantRotation), intent(inout) :: self
    type            (treeNode                                      ), intent(inout) :: node
    double precision                                                , intent(in   ) :: radius
    class           (nodeComponentBasic                            ), pointer       :: basic
    class           (nodeComponentSpin                             ), pointer       :: spin
    class           (nodeComponentHotHalo                          ), pointer       :: hotHalo
    class           (massDistributionClass                         ), pointer       :: massDistribution_
    double precision                                                                :: angularMomentumSpecificMean, normalizationRotation

    ! Check if node differs from previous one for which we performed calculations.
    if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node,node%uniqueID())
    ! Check if specific angular momentum of cooling gas is already computed.
    if (.not.self%angularMomentumSpecificComputed) then
       ! Flag that cooling radius is now computed.
       self%angularMomentumSpecificComputed=.true.
       ! Compute the mean specific angular momentum.
       select case (self%sourceAngularMomentumSpecificMean%ID)
       case (angularMomentumSourceDarkMatter%ID)
          ! Compute mean specific angular momentum of the dark matter halo.
          basic                       =>  node               %basic()
          spin                        =>  node               %spin ()
          angularMomentumSpecificMean =  +spin%angularMomentum     () &
               &                         /basic              %mass ()
       case (angularMomentumSourceHotGas    %ID)
          ! Compute mean specific angular momentum from the hot halo component.
          hotHalo => node%hotHalo()
          angularMomentumSpecificMean=+hotHalo%angularMomentum() &
               &                      /hotHalo%mass           ()
       case default
          angularMomentumSpecificMean=0.0d0
          call Error_Report('unknown profile type'//{introspection:location})
       end select
       ! Compute the rotation normalization.
       select case (self%sourceNormalizationRotation      %ID)
       case (angularMomentumSourceDarkMatter%ID)
          massDistribution_     =>  self             %darkMatterProfileDMO_%get                (node                                                                                         )
          normalizationRotation =  +massDistribution_                      %densityRadialMoment(moment       =2.0d0               ,radiusMaximum=self%darkMatterHaloScale_%radiusVirial(node)) &
               &                   /massDistribution_                      %densityRadialMoment(moment       =3.0d0               ,radiusMaximum=self%darkMatterHaloScale_%radiusVirial(node))
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]          
       case (angularMomentumSourceHotGas    %ID)
          hotHalo               =>  node                                   %hotHalo            (                                                                                             )
          massDistribution_     =>  node                                   %massDistribution   (componentType=componentTypeHotHalo,massType     =massTypeGaseous                             )
          normalizationRotation =  +massDistribution_                      %densityRadialMoment(moment       =2.0d0               ,radiusMaximum=hotHalo                  %outerRadius (    )) &
               &                   /massDistribution_                      %densityRadialMoment(moment       =3.0d0               ,radiusMaximum=hotHalo                  %outerRadius (    ))
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]          
       case default
          normalizationRotation=0.0d0
          call Error_Report('unknown profile type'//{introspection:location})
       end select
       ! Compute the specific angular momentum of the cooling gas.
       self%angularMomentumSpecificPrevious=+normalizationRotation       &
            &                               *angularMomentumSpecificMean
    end if
    ! Check that the radius is positive.
    if (radius > 0.0d0) then
       ! Return the computed value.
       if (self%useInteriorMean) then
          ! Find the specific angular momentum interior to the specified radius.
          massDistribution_                       =>  node             %massDistribution               (componentType=componentTypeHotHalo,massType     =massTypeGaseous)
          constantRotationAngularMomentumSpecific =  +self             %angularMomentumSpecificPrevious                                                                   &
               &                                     *massDistribution_%densityRadialMoment            (moment       =3.0d0               ,radiusMaximum=radius         ) &
               &                                     /massDistribution_%densityRadialMoment            (moment       =2.0d0               ,radiusMaximum=radius         )
          !![
	  <objectDestructor name="massDistribution_"/>
          !!]          
       else
          ! Find the specific angular momentum at the specified radius.
          constantRotationAngularMomentumSpecific =  +self             %angularMomentumSpecificPrevious                                                                   &
               &                                     *                                                                                                   radius
       end if
    else
       ! Radius is non-positive - return zero.
       constantRotationAngularMomentumSpecific    =  +0.0d0
    end if
    return
  end function constantRotationAngularMomentumSpecific

  subroutine constantRotationCalculationReset(self,node,uniqueID)
    !!{
    Reset the specific angular momentum of cooling gas calculation.
    !!}
    use :: Kind_Numbers, only : kind_int8
    implicit none
    class  (coolingSpecificAngularMomentumConstantRotation), intent(inout) :: self
    type   (treeNode                                      ), intent(inout) :: node
    integer(kind_int8                                     ), intent(in   ) :: uniqueID
    !$GLC attributes unused :: node

    self%angularMomentumSpecificComputed=.false.
    self%lastUniqueID                   =uniqueID
    return
  end subroutine constantRotationCalculationReset
