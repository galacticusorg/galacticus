!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo scales based on virial density contrast.

  !# <darkMatterHaloScale name="darkMatterHaloScaleVirialDensityContrastDefinition">
  !#  <description>Dark matter halo scales derived from virial density contrasts.</description>
  !# </darkMatterHaloScale>
  use Kind_Numbers
  use Tables
  use Virial_Density_Contrast

  type, extends(darkMatterHaloScaleClass) :: darkMatterHaloScaleVirialDensityContrastDefinition
     !% A dark matter halo scale contrast class using virial density contrasts.
     private
     ! Virial density contrast object.
     class           (virialDensityContrastClass), pointer   :: virialDensityContrastDefinition => null()
     ! Record of unique ID of node which we last computed results for.
     integer         (kind=kind_int8            )            :: lastUniqueID
     ! Record of whether or not halo scales have already been computed for this node.
     logical                                                 :: dynamicalTimescaleComputed               , virialRadiusComputed            , &
          &                                                     virialTemperatureComputed                , virialVelocityComputed
     ! Stored values of halo scales.
     double precision                                        :: dynamicalTimescaleStored                 , virialRadiusStored              , &
          &                                                     virialTemperatureStored                  , virialVelocityStored            , &
          &                                                     timePrevious                             , densityGrowthRatePrevious       , &
          &                                                     massPrevious
     ! Table for fast lookup of the mean density of halos.
     double precision                                        :: meanDensityTimeMaximum                   , meanDensityTimeMinimum   =-1.0d0
     type            (table1DLogarithmicLinear  )            :: meanDensityTable
     logical                                                 :: resetMeanDensityTable
   contains
     final     ::                             virialDensityContrastDefinitionDestructor
     procedure :: dynamicalTimescale       => virialDensityContrastDefinitionDynamicalTimescale
     procedure :: virialVelocity           => virialDensityContrastDefinitionVirialVelocity
     procedure :: virialVelocityGrowthRate => virialDensityContrastDefinitionVirialVelocityGrowthRate
     procedure :: virialTemperature        => virialDensityContrastDefinitionVirialTemperature
     procedure :: virialRadius             => virialDensityContrastDefinitionVirialRadius
     procedure :: virialRadiusGrowthRate   => virialDensityContrastDefinitionVirialRadiusGrowthRate
     procedure :: meanDensity              => virialDensityContrastDefinitionMeanDensity
     procedure :: meanDensityGrowthRate    => virialDensityContrastDefinitionMeanDensityGrowthRate
     procedure :: calculationReset         => virialDensityContrastDefinitionCalculationReset
     procedure :: stateStore               => virialDensityContrastDefinitionStateStore
     procedure :: stateRestore             => virialDensityContrastDefinitionStateRestore
  end type darkMatterHaloScaleVirialDensityContrastDefinition

  interface darkMatterHaloScaleVirialDensityContrastDefinition
     !% Constructors for the {\tt virialDensityContrastDefinition} dark matter halo scales class.
     module procedure virialDensityContrastDefinitionDefaultConstructor
     module procedure virialDensityContrastDefinitionConstructor
  end interface darkMatterHaloScaleVirialDensityContrastDefinition

  integer, parameter :: virialDensityContrastDefinitionMeanDensityTablePointsPerDecade=100

contains

  function virialDensityContrastDefinitionDefaultConstructor()
    !% Default constructor for the {\tt virialDensityContrastDefinition} dark matter halo scales class.
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition), target  :: virialDensityContrastDefinitionDefaultConstructor
    class(virialDensityContrastClass              ), pointer :: virialDensityContrast_

    virialDensityContrast_                            => virialDensityContrast                     (                      )
    virialDensityContrastDefinitionDefaultConstructor =  virialDensityContrastDefinitionConstructor(virialDensityContrast_)
    return
  end function virialDensityContrastDefinitionDefaultConstructor

  function virialDensityContrastDefinitionConstructor(virialDensityContrastDefinition)
    !% Default constructor for the {\tt virialDensityContrastDefinition} dark matter halo scales class.
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition)               , target :: virialDensityContrastDefinitionConstructor
    class(virialDensityContrastClass                        ), intent(in   ), target :: virialDensityContrastDefinition

    virialDensityContrastDefinitionConstructor%virialDensityContrastDefinition => virialDensityContrastDefinition
    virialDensityContrastDefinitionConstructor%lastUniqueID                    =  -1_kind_int8
    virialDensityContrastDefinitionConstructor%dynamicalTimescaleComputed      =  .false.
    virialDensityContrastDefinitionConstructor%virialRadiusComputed            =  .false.
    virialDensityContrastDefinitionConstructor%virialTemperatureComputed       =  .false.
    virialDensityContrastDefinitionConstructor%virialVelocityComputed          =  .false.
    virialDensityContrastDefinitionConstructor%meanDensityTimeMaximum          =  -1.0d0
    virialDensityContrastDefinitionConstructor%meanDensityTimeMinimum          =  -1.0d0
    virialDensityContrastDefinitionConstructor%resetMeanDensityTable           =  .false.
    virialDensityContrastDefinitionConstructor%timePrevious                    =  -1.0d0
    virialDensityContrastDefinitionConstructor%massPrevious                    =  -1.0d0
    return
  end function virialDensityContrastDefinitionConstructor

  subroutine virialDensityContrastDefinitionDestructor(self)
    !% Destructir for the {\tt virialDensityContrastDefinition} dark matter halo scales class.
    implicit none
    type (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self

    if (associated(self%virialDensityContrastDefinition).and.self%virialDensityContrastDefinition%isFinalizable()) deallocate(self%virialDensityContrastDefinition)
    if (self%meanDensityTimeMinimum >= 0.0d0) call self%meanDensityTable%destroy()
    return
  end subroutine virialDensityContrastDefinitionDestructor

  subroutine virialDensityContrastDefinitionCalculationReset(self,thisNode)
    !% Reset the halo scales calculation.
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode

    self%virialRadiusComputed      =.false.
    self%virialTemperatureComputed =.false.
    self%virialVelocityComputed    =.false.
    self%dynamicalTimescaleComputed=.false.
    self%lastUniqueID              =thisNode%uniqueID()
    return
  end subroutine virialDensityContrastDefinitionCalculationReset

  double precision function virialDensityContrastDefinitionDynamicalTimescale(self,thisNode)
    !% Returns the dynamical timescale for {\tt thisNode}.
    use Numerical_Constants_Astronomical
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(thisNode)
    ! Check if halo dynamical timescale is already computed. Compute and store if not.
    if (.not.self%dynamicalTimescaleComputed) then
       self%dynamicalTimescaleComputed= .true.
       self%dynamicalTimescaleStored  = self%virialRadius  (thisNode) &
            &                          /self%virialVelocity(thisNode) &
            &                          *(megaParsec/kilo/gigaYear)
     end if
    ! Return the stored timescale.
    virialDensityContrastDefinitionDynamicalTimescale=self%dynamicalTimescaleStored
    return
  end function virialDensityContrastDefinitionDynamicalTimescale

  double precision function virialDensityContrastDefinitionVirialVelocity(self,thisNode)
    !% Returns the virial velocity scale for {\tt thisNode}.
    use Numerical_Constants_Physical
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic                                )               , pointer :: thisBasic

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(thisNode)
    ! Check if virial velocity is already computed. Compute and store if not.
    if (.not.self%virialVelocityComputed) then
       ! Get the basic component.
       thisBasic => thisNode%basic()
       ! Compute the virial velocity.
       self%virialVelocityStored=sqrt(gravitationalConstantGalacticus*thisBasic%mass() &
            &/self%virialRadius(thisNode))
       ! Record that virial velocity has now been computed.
       self%virialVelocityComputed=.true.
    end if
    ! Return the stored virial velocity.
    virialDensityContrastDefinitionVirialVelocity=self%virialVelocityStored
    return
  end function virialDensityContrastDefinitionVirialVelocity

  double precision function virialDensityContrastDefinitionVirialVelocityGrowthRate(self,thisNode)
    !% Returns the growth rate of the virial velocity scale for {\tt thisNode}.
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic                                )               , pointer :: thisBasic

    ! Get the basic component.
    thisBasic                                              => thisNode%basic()
    virialDensityContrastDefinitionVirialVelocityGrowthRate= &
         & +0.5d0                                            &
         & *  self%virialVelocity        (thisNode)          &
         & *(                                                &
         &    thisBasic%accretionRate    (        )          &
         &   /thisBasic%mass             (        )          &
         &   -self%virialRadiusGrowthRate(thisNode)          &
         &   /self%virialRadius          (thisNode)          &
         &  )
    return
  end function virialDensityContrastDefinitionVirialVelocityGrowthRate

  double precision function virialDensityContrastDefinitionVirialTemperature(self,thisNode)
    !% Returns the virial temperature (in Kelvin) for {\tt thisNode}.
    use Numerical_Constants_Physical
    use Numerical_Constants_Astronomical
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(thisNode)
    ! Check if virial temperature is already computed. Compute and store if not.
    if (.not.self%virialTemperatureComputed) then
       self%virialTemperatureComputed=.true.
       self%virialTemperatureStored=0.5d0*atomicMassUnit*meanAtomicMassPrimordial*((kilo&
            &*self%virialVelocity(thisNode))**2)/boltzmannsConstant
    end if
    ! Return the stored temperature.
    virialDensityContrastDefinitionVirialTemperature=self%virialTemperatureStored
    return
  end function virialDensityContrastDefinitionVirialTemperature

  double precision function virialDensityContrastDefinitionVirialRadius(self,thisNode)
    !% Returns the virial radius scale for {\tt thisNode}.
    use Numerical_Constants_Math
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic                                )               , pointer :: thisBasic

    ! Check if node differs from previous one for which we performed calculations.
    if (thisNode%uniqueID() /= self%lastUniqueID) call self%calculationReset(thisNode)
    ! Check if virial radius is already computed. Compute and store if not.
    if (.not.self%virialRadiusComputed) then
       ! Get the basic component.
       thisBasic => thisNode%basic()
       ! Compute the virial radius.
       self%virialRadiusStored=(3.0d0*thisBasic%mass()/4.0d0/Pi/self%meanDensity(thisNode))**(1.0d0/3.0d0)
       ! Record that the virial radius has been computed.
       self%virialRadiusComputed=.true.
    end if
    ! Return the stored value.
    virialDensityContrastDefinitionVirialRadius=self%virialRadiusStored
    return
  end function virialDensityContrastDefinitionVirialRadius

  double precision function virialDensityContrastDefinitionVirialRadiusGrowthRate(self,thisNode)
    !% Returns the growth rate of the virial radius scale for {\tt thisNode}.
    implicit none
    class(darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type (treeNode                                          ), intent(inout), pointer :: thisNode
    class(nodeComponentBasic)                                               , pointer :: thisBasic

    ! Get the basic component.
    thisBasic => thisNode%basic()
    virialDensityContrastDefinitionVirialRadiusGrowthRate=(1.0d0/3.0d0)*self%virialRadius(thisNode)&
         &*(thisBasic%accretionRate()/thisBasic%mass()-self%meanDensityGrowthRate(thisNode)&
         &/virialDensityContrastDefinitionMeanDensity(self,thisNode))
    return
  end function virialDensityContrastDefinitionVirialRadiusGrowthRate

  double precision function virialDensityContrastDefinitionMeanDensity(self,thisNode)
    !% Returns the mean density for {\tt thisNode}.
    use Cosmology_Parameters
    use Cosmology_Functions
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type            (treeNode                                          ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic                                )               , pointer :: thisBasic
    class           (cosmologyParametersClass                          )               , pointer :: thisCosmologyParameters
    class           (cosmologyFunctionsClass                           )               , pointer :: cosmologyFunctionsDefault
    integer                                                                                      :: i                        , meanDensityTablePoints
    double precision                                                                             :: time

    ! Get the basic component.
    thisBasic => thisNode%basic()
    ! Get the time at which this halo was last an isolated halo.
    time=thisBasic%timeLastIsolated()
    if (time <= 0.0d0) time=thisBasic%time()
    ! For mass-dependent virial density contrasts we must always recompute the result.
    if (self%virialDensityContrastDefinition%isMassDependent()) then
       ! Get default objects.
       thisCosmologyParameters   => cosmologyParameters  ()
       cosmologyFunctionsDefault => cosmologyFunctions   ()
       virialDensityContrastDefinitionMeanDensity=&
            &+self%virialDensityContrastDefinition%densityContrast(thisBasic%mass(),time)&
            &*thisCosmologyParameters             %OmegaMatter    (                                           )     &
                  &  *thisCosmologyParameters             %densityCritical(                                           )     &
                  &  /cosmologyFunctionsDefault           %expansionFactor(                 time)**3


    else
       ! For non-mass-dependent virial density contrasts we can tabulate as a function of time.
       ! Retabulate the mean density vs. time if necessary.
       if (self%resetMeanDensityTable .or. time < self%meanDensityTimeMinimum .or. time > self%meanDensityTimeMaximum) then
          self%resetMeanDensityTable=.false.
          if (self%meanDensityTimeMinimum <= 0.0d0) then
             self%meanDensityTimeMinimum=                                time/10.0d0
             self%meanDensityTimeMaximum=                                time* 2.0d0
          else
             self%meanDensityTimeMinimum=min(self%meanDensityTimeMinimum,time/10.0d0)
             self%meanDensityTimeMaximum=max(self%meanDensityTimeMaximum,time* 2.0d0)
          end if
          meanDensityTablePoints=int(log10(self%meanDensityTimeMaximum/self%meanDensityTimeMinimum)*dble(virialDensityContrastDefinitionMeanDensityTablePointsPerDecade))+1
          call self%meanDensityTable%destroy()
          call self%meanDensityTable%create(self%meanDensityTimeMinimum,self%meanDensityTimeMaximum,meanDensityTablePoints)
          ! Get default objects.
          thisCosmologyParameters   => cosmologyParameters  ()
          cosmologyFunctionsDefault => cosmologyFunctions   ()
          do i=1,meanDensityTablePoints
             call self%meanDensityTable%populate                                                                            &
                  & (                                                                                                       &
                  &  +self%virialDensityContrastDefinition%densityContrast(thisBasic%mass(),self%meanDensityTable%x(i))     &
                  &  *thisCosmologyParameters             %OmegaMatter    (                                           )     &
                  &  *thisCosmologyParameters             %densityCritical(                                           )     &
                  &  /cosmologyFunctionsDefault           %expansionFactor(                 self%meanDensityTable%x(i))**3, &
                  &  i                                                                                                      &
                  & )
          end do
       end if
       ! Return the stored value.
       virialDensityContrastDefinitionMeanDensity=self%meanDensityTable%interpolate(time)
    end if
    return
  end function virialDensityContrastDefinitionMeanDensity

  double precision function virialDensityContrastDefinitionMeanDensityGrowthRate(self,thisNode)
    !% Returns the growth rate of the mean density for {\tt thisNode}.
    use Cosmology_Functions
    implicit none
    class           (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout)          :: self
    type            (treeNode                                          ), intent(inout), pointer :: thisNode
    class           (nodeComponentBasic                                )               , pointer :: thisBasic
    class           (cosmologyFunctionsClass                           )               , pointer :: cosmologyFunctionsDefault
    double precision                                                                             :: aExpansion               , time

    if (thisNode%isSatellite()) then
       ! Satellite halo is not growing, return zero rate.
       virialDensityContrastDefinitionMeanDensityGrowthRate=0.0d0
    else
       ! Get the basic component.
       thisBasic => thisNode%basic()
       ! Get the time at which this halo was last an isolated halo.
       time=thisBasic%timeLastIsolated()
       ! Check if the time is different from that one previously used.
       if (time /= self%timePrevious .or. thisBasic%mass() /= self%massPrevious) then
          ! It is not, so recompute the density growth rate.
          self%timePrevious=time
          self%massPrevious=thisBasic%mass()
          ! Get default objects.
          cosmologyFunctionsDefault => cosmologyFunctions   ()
          ! Get the expansion factor at this time.
          aExpansion=cosmologyFunctionsDefault%expansionFactor(time)
          ! Compute growth rate of its mean density based on mean cosmological density and overdensity of a collapsing halo.
          self%densityGrowthRatePrevious=                                                                         &
               & self%meanDensity(thisNode)                                                                       &
               & *(                                                                                               &
               &   +self%virialDensityContrastDefinition%densityContrastRateOfChange(thisBasic%mass(),time      ) &
               &   /self%virialDensityContrastDefinition%densityContrast            (thisBasic%mass(),time      ) &
               &   -3.0d0                                                                                         &
               &   *cosmologyFunctionsDefault           %expansionRate              (aExpansion                 ) &
               &  )
       end if
       ! Return the stored value.
       virialDensityContrastDefinitionMeanDensityGrowthRate=self%densityGrowthRatePrevious
    end if
    return
  end function virialDensityContrastDefinitionMeanDensityGrowthRate

  subroutine virialDensityContrastDefinitionStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    class  (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    integer                                          , intent(in   ) :: stateFile
    type   (fgsl_file                               ), intent(in   ) :: fgslStateFile

    write (stateFile) self%meanDensityTimeMinimum,self%meanDensityTimeMaximum
    return
  end subroutine virialDensityContrastDefinitionStateStore

  subroutine virialDensityContrastDefinitionStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    class  (darkMatterHaloScaleVirialDensityContrastDefinition), intent(inout) :: self
    integer                                                    , intent(in   ) :: stateFile
    type   (fgsl_file                                         ), intent(in   ) :: fgslStateFile

    read (stateFile) self%meanDensityTimeMinimum,self%meanDensityTimeMaximum
    ! Ensure that interpolation objects will get reset.
    self%resetMeanDensityTable=.true.
    return
  end subroutine virialDensityContrastDefinitionStateRestore
