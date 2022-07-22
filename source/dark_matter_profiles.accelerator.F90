!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  An accelerator class for non-dark-matter-only dark matter halo profiles.
  !!}
  
  use :: Binary_Search_Trees, only : binaryTree
  use :: Kind_Numbers       , only : kind_int8
  
  !![
  <darkMatterProfile name="darkMatterProfileAccelerator">
   <description>An accelerator class for non-dark-matter-only dark matter halo profiles.</description>
  </darkMatterProfile>
  !!]
  type, extends(darkMatterProfileClass) :: darkMatterProfileAccelerator
     !!{
     An accelerator class for non-dark-matter-only dark matter halo profiles.
     !!}
     private
     integer         (kind_int8             ), dimension(2) :: uniqueIDPrevious
     class           (darkMatterProfileClass), pointer      :: darkMatterProfile_             => null()
     type            (binaryTree            ), dimension(2) :: treeMassEnclosed
     integer                                                :: treePrevious
     double precision                                       :: toleranceRelative                       , factorRadiusMaximum, &
          &                                                    factorRadiusLogarithmicMaximum
   contains
     !![
     <methods>
       <method description="Reset memoized calculations." method="calculationReset" />
     </methods>
     !!]
     final     ::                                      acceleratorDestructor
     procedure :: autoHook                          => acceleratorAutoHook
     procedure :: calculationReset                  => acceleratorCalculationReset
     procedure :: density                           => acceleratorDensity
     procedure :: densityLogSlope                   => acceleratorDensityLogSlope
     procedure :: radiusEnclosingDensity            => acceleratorRadiusEnclosingDensity
     procedure :: radiusEnclosingMass               => acceleratorRadiusEnclosingMass
     procedure :: radialMoment                      => acceleratorRadialMoment
     procedure :: enclosedMass                      => acceleratorEnclosedMass
     procedure :: potential                         => acceleratorPotential
     procedure :: circularVelocity                  => acceleratorCircularVelocity
     procedure :: circularVelocityMaximum           => acceleratorCircularVelocityMaximum
     procedure :: radialVelocityDispersion          => acceleratorRadialVelocityDispersion
     procedure :: radiusFromSpecificAngularMomentum => acceleratorRadiusFromSpecificAngularMomentum
     procedure :: rotationNormalization             => acceleratorRotationNormalization
     procedure :: energy                            => acceleratorEnergy
     procedure :: kSpace                            => acceleratorKSpace
     procedure :: freefallRadius                    => acceleratorFreefallRadius
     procedure :: freefallRadiusIncreaseRate        => acceleratorFreefallRadiusIncreaseRate
  end type darkMatterProfileAccelerator

  interface darkMatterProfileAccelerator
     !!{
     Constructors for the {\normalfont \ttfamily accelerator} non-dark-matter-only dark matter halo profile class.
     !!}
     module procedure acceleratorConstructorParameters
     module procedure acceleratorConstructorInternal
  end interface darkMatterProfileAccelerator

contains

  function acceleratorConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily accelerator} non-dark-matter-only dark matter halo profile class which takes
    a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (darkMatterProfileAccelerator)                :: self
    type            (inputParameters             ), intent(inout) :: parameters
    class           (darkMatterProfileClass      ), pointer       :: darkMatterProfile_
    class           (darkMatterHaloScaleClass    ), pointer       :: darkMatterHaloScale_
    double precision                                              :: toleranceRelative   , factorRadiusMaximum

    !![
    <inputParameter>
      <name>toleranceRelative</name>
      <defaultValue>1.0d-2</defaultValue>
      <source>parameters</source>
      <description>The toleranceRelative with which to accept accelerated estimates.</description>
    </inputParameter>
    <inputParameter>
      <name>factorRadiusMaximum</name>
      <defaultValue>3.0d0</defaultValue>
      <source>parameters</source>
      <description>The maximum factor by which to interpolate in radius.</description>
    </inputParameter>
    <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    <objectBuilder class="darkMatterProfile"   name="darkMatterProfile_"   source="parameters"/>
    !!]
    self=darkMatterProfileAccelerator(toleranceRelative,factorRadiusMaximum,darkMatterHaloScale_,darkMatterProfile_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="darkMatterHaloScale_"/>
    <objectDestructor name="darkMatterProfile_"  />
    !!]
    return
  end function acceleratorConstructorParameters

  function acceleratorConstructorInternal(toleranceRelative,factorRadiusMaximum,darkMatterHaloScale_,darkMatterProfile_) result(self)
    !!{
    Generic constructor for the {\normalfont \ttfamily accelerator} dark matter profile class.
    !!}
    implicit none
    type            (darkMatterProfileAccelerator)                        :: self
    class           (darkMatterProfileClass      ), intent(in   ), target :: darkMatterProfile_
    class           (darkMatterHaloScaleClass    ), intent(in   ), target :: darkMatterHaloScale_
    double precision                              , intent(in)            :: toleranceRelative   , factorRadiusMaximum
    !![
    <constructorAssign variables="toleranceRelative, factorRadiusMaximum, *darkMatterHaloScale_, *darkMatterProfile_"/>
    !!]

    self%uniqueIDPrevious              =-1_kind_int8
    self%treePrevious                  =+1
    self%factorRadiusLogarithmicMaximum=+log(sqrt(factorRadiusMaximum))
    return
  end function acceleratorConstructorInternal

  subroutine acceleratorAutoHook(self)
    !!{
    Attach to the calculation reset event.
    !!}
    use :: Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(darkMatterProfileAccelerator), intent(inout) :: self

    call calculationResetEvent%attach(self,acceleratorCalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine acceleratorAutoHook

  subroutine acceleratorDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily accelerator} dark matter halo profile class.
    !!}
    use :: Events_Hooks, only : calculationResetEvent
    implicit none
    type(darkMatterProfileAccelerator), intent(inout) :: self

    !![
    <objectDestructor name="self%darkMatterHaloScale_"/>
    <objectDestructor name="self%darkMatterProfile_"  />
    !!]
    if (calculationResetEvent%isAttached(self,acceleratorCalculationReset)) call calculationResetEvent%detach(self,acceleratorCalculationReset)
    return
  end subroutine acceleratorDestructor

  subroutine acceleratorCalculationReset(self,node)
    !!{
    Reset the dark matter profile calculation.
    !!}
    implicit none
    class  (darkMatterProfileAccelerator), intent(inout) :: self
    type   (treeNode                    ), intent(inout) :: node
    integer                                              :: i

    ! Trees are maintained for two nodes - this is often advantageous as queries are often made for satellite and host nodes
    ! together. If the current node is one for which we currently have a tree, invalidate that tree. Otherwise, if the current
    ! node is a satellite in the current host node, place it into the second tree, otherwise, into the first.
    if      (node%uniqueID() == self%uniqueIDPrevious(1)) then
       i      =+1
    else if (node%uniqueID() == self%uniqueIDPrevious(2)) then
       i      =+2
    else
       if (node%isSatellite() .and. node%parent%uniqueID() == self%uniqueIDPrevious(1)) then
          i=2
       else
          i=1
       end if
    end if
    self%uniqueIDPrevious(i)=node%uniqueID()
    self%treePrevious       =i
    if (associated(self%treeMassEnclosed(i)%root)) deallocate(self%treeMassEnclosed(i)%root)
    return
  end subroutine acceleratorCalculationReset

  double precision function acceleratorDensity(self,node,radius)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    acceleratorDensity=self%darkMatterProfile_%density(node,radius)
    return
  end function acceleratorDensity

  double precision function acceleratorDensityLogSlope(self,node,radius)
    !!{
    Returns the logarithmic slope of the density in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    acceleratorDensityLogSlope=self%darkMatterProfile_%densityLogSlope(node,radius)
    return
  end function acceleratorDensityLogSlope

  double precision function acceleratorRadiusEnclosingDensity(self,node,density)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily density} (given in units of $M_\odot/$Mpc$^{-3}$).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: density

    acceleratorRadiusEnclosingDensity=self%darkMatterProfile_%radiusEnclosingDensity(node,density)
    return
  end function acceleratorRadiusEnclosingDensity

  double precision function acceleratorRadiusEnclosingMass(self,node,mass)
    !!{
    Returns the radius (in Mpc) in the dark matter profile of {\normalfont \ttfamily node} which encloses the given
    {\normalfont \ttfamily mass} (given in units of $M_\odot$).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout), target :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: mass

    acceleratorRadiusEnclosingMass=self%darkMatterProfile_%radiusEnclosingMass(node,mass)
    return
  end function acceleratorRadiusEnclosingMass

  double precision function acceleratorRadialMoment(self,node,moment,radiusMinimum,radiusMaximum)
    !!{
    Returns the density (in $M_\odot$ Mpc$^{-3}$) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout)           :: self
    type            (treeNode                    ), intent(inout)           :: node
    double precision                              , intent(in   )           :: moment
    double precision                              , intent(in   ), optional :: radiusMinimum, radiusMaximum

    acceleratorRadialMoment=self%darkMatterProfile_%radialMoment(node,moment,radiusMinimum,radiusMaximum)
    return
  end function acceleratorRadialMoment

  double precision function acceleratorEnclosedMass(self,node,radius)
    !!{
    Returns the enclosed mass (in $M_\odot$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont \ttfamily radius} (given in
    units of Mpc).
    !!}
    use :: Binary_Search_Trees, only : binaryTreeNode
    use :: Numerical_Comparison, only : Values_Agree
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius
    type            (binaryTreeNode              ), pointer       :: left1            , left2        , &
         &                                                           right1           , right2
    double precision                                              :: massEnclosed1    , massEnclosed2, &
         &                                                           radiusLogarithmic
    logical                                                       :: found
    integer                                                       :: i
    
    if      (node%uniqueID() == self%uniqueIDPrevious(1)) then
       i=1
    else if (node%uniqueID() == self%uniqueIDPrevious(2)) then
       i=2
    else
       call self%calculationReset(node)
       i=self%treePrevious
    end if
    found=.false.
    radiusLogarithmic=log(radius)
    call self%treeMassEnclosed(i)%bracket(radiusLogarithmic,left1,right1)
    if (associated(left1).and.associated(right1)) then
       if (associated(left1,right1)) then
          acceleratorEnclosedMass=exp(left1%value)
          found                  =.true.
       else
          if     (                                                                     &
               &   +radiusLogarithmic- left1%key < self%factorRadiusLogarithmicMaximum &
               &  .and.                                                                &
               &   -radiusLogarithmic+right1%key < self%factorRadiusLogarithmicMaximum &
               & ) then
             left2  =>  left1%predecessor()
             right2 => right1%  successor()
             if (associated(left2).and.associated(right2)) then
                massEnclosed1=(radiusLogarithmic-left1%key)*(right1%value-left1%value)/(right1%key-left1%key)+left1%value
                massEnclosed2=(radiusLogarithmic-left2%key)*(right2%value-left2%value)/(right2%key-left2%key)+left2%value
                if (Values_Agree(massEnclosed1,massEnclosed2,relTol=self%toleranceRelative)) then
                   acceleratorEnclosedMass=exp(massEnclosed1)
                   found                  =.true.
                end if
             end if
          end if
       end if
    end if
    if (.not.found) then
       acceleratorEnclosedMass=self%darkMatterProfile_%enclosedMass(node,radius)
       call self%treeMassEnclosed(i)%insert(radiusLogarithmic,log(acceleratorEnclosedMass))
    end if
    return
  end function acceleratorEnclosedMass

  double precision function acceleratorPotential(self,node,radius,status)
    !!{
    Returns the potential (in (km/s)$^2$) in the dark matter profile of {\normalfont \ttfamily node} at the given {\normalfont
    \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator     ), intent(inout), target   :: self
    type            (treeNode                         ), intent(inout), target   :: node
    double precision                                   , intent(in   )           :: radius
    type            (enumerationStructureErrorCodeType), intent(  out), optional :: status

    acceleratorPotential=self%darkMatterProfile_%potential(node,radius,status)
    return
  end function acceleratorPotential

  double precision function acceleratorCircularVelocity(self,node,radius)
    !!{
    Returns the circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                       ), intent(inout) :: node
    double precision                                 , intent(in   ) :: radius

    acceleratorCircularVelocity=self%darkMatterProfile_%circularVelocity(node,radius)
    return
  end function acceleratorCircularVelocity

  double precision function acceleratorCircularVelocityMaximum(self,node)
    !!{
    Returns the maximum circular velocity (in km/s) in the dark matter profile of {\normalfont \ttfamily node}.
    !!}
    implicit none
    class(darkMatterProfileAccelerator), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    acceleratorCircularVelocityMaximum=self%darkMatterProfile_%circularVelocityMaximum(node)
    return
  end function acceleratorCircularVelocityMaximum

  double precision function acceleratorRadialVelocityDispersion(self,node,radius)
    !!{
    Returns the radial velocity dispersion (in km/s) in the dark matter profile of {\normalfont \ttfamily node} at the given
    {\normalfont \ttfamily radius} (given in units of Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: radius

    acceleratorRadialVelocityDispersion=self%darkMatterProfile_%radialVelocityDispersion(node,radius)
    return
  end function acceleratorRadialVelocityDispersion

  double precision function acceleratorRadiusFromSpecificAngularMomentum(self,node,specificAngularMomentum)
    !!{
    Returns the radius (in Mpc) in {\normalfont \ttfamily node} at which a circular orbit has the given {\normalfont \ttfamily specificAngularMomentum} (given
    in units of km s$^{-1}$ Mpc).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout), target :: self
    type            (treeNode                    ), intent(inout)         :: node
    double precision                              , intent(in   )         :: specificAngularMomentum

    acceleratorRadiusFromSpecificAngularMomentum=self%darkMatterProfile_%radiusFromSpecificAngularMomentum(node,specificAngularMomentum)
    return
  end function acceleratorRadiusFromSpecificAngularMomentum

  double precision function acceleratorRotationNormalization(self,node)
    !!{
    Return the normalization of the rotation velocity vs. specific angular momentum relation.
    !!}
    implicit none
    class(darkMatterProfileAccelerator), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    acceleratorRotationNormalization=self%darkMatterProfile_%rotationNormalization(node)
    return
  end function acceleratorRotationNormalization

  double precision function acceleratorEnergy(self,node)
    !!{
    Return the energy of the dark matter halo density profile.
    !!}
    implicit none
    class(darkMatterProfileAccelerator), intent(inout) :: self
    type (treeNode                    ), intent(inout) :: node

    acceleratorEnergy=self%darkMatterProfile_%energy(node)
    return
  end function acceleratorEnergy

  double precision function acceleratorKSpace(self,node,waveNumber)
    !!{
    Returns the Fourier transform of the dark matter halo density profile at the specified {\normalfont \ttfamily waveNumber}
    (given in Mpc$^{-1}$).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout)         :: self
    type            (treeNode                    ), intent(inout), target :: node
    double precision                              , intent(in   )         :: waveNumber

    acceleratorKSpace=+self%darkMatterProfile_%kSpace(node,waveNumber)
    return
  end function acceleratorKSpace

  double precision function acceleratorFreefallRadius(self,node,time)
    !!{
    Returns the freefall radius in the dark matter halo density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: time

    acceleratorFreefallRadius=self%darkMatterProfile_%freefallRadius(node,time)
    return
  end function acceleratorFreefallRadius

  double precision function acceleratorFreefallRadiusIncreaseRate(self,node,time)
    !!{
    Returns the freefall radius in the dark matter halo density profile at the specified {\normalfont \ttfamily time} (given in
    Gyr).
    !!}
    implicit none
    class           (darkMatterProfileAccelerator), intent(inout) :: self
    type            (treeNode                    ), intent(inout) :: node
    double precision                              , intent(in   ) :: time

    acceleratorFreefallRadiusIncreaseRate=self%darkMatterProfile_%freefallRadiusIncreaseRate(node,time)
    return
  end function acceleratorFreefallRadiusIncreaseRate
