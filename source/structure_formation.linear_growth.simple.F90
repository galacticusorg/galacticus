!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017
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

  !% An implementation of linear growth of cosmological structure in simple cosomologies. Ignores pressure terms for the growth of
  !% baryons and has no wavenumber dependence. Also assumes no growth of radiation perturbations.

  !# <linearGrowth name="linearGrowthSimple" defaultThreadPrivate="yes">
  !#  <description>Linear growth of cosmological structure in simple cosomologies. Ignores pressure terms for the growth of baryons and has no wavenumber dependence. Also assumes no growth of radiation perturbations.</description>
  !# </linearGrowth>
  use Tables
  use Cosmology_Parameters, only : cosmologyParameters, cosmologyParametersClass, hubbleUnitsTime
  use Cosmology_Functions , only : cosmologyFunctions , cosmologyFunctionsClass

  type, extends(linearGrowthClass) :: linearGrowthSimple
     !% A linear growth of cosmological structure contrast class in simple cosomologies.
     private
     logical                                                 :: tableInitialized
     double precision                                        :: tableTimeMinimum            , tableTimeMaximum, &
          &                                                     normalizationMatterDominated
     class           (table1D                 ), allocatable :: growthFactor
     class           (cosmologyParametersClass), pointer     :: cosmologyParameters_
     class           (cosmologyFunctionsClass ), pointer     :: cosmologyFunctions_
   contains
     !@ <objectMethods>
     !@   <object>linearGrowthSimple</object>
     !@   <objectMethod>
     !@     <method>retabulate</method>
     !@     <type>void</type>
     !@     <arguments>\doublezero\ time\argin</arguments>
     !@     <description>Tabulate linear growth factor.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                         simpleDestructor
     procedure :: stateStore                           => simpleStateStore
     procedure :: stateRestore                         => simpleStateRestore
     procedure :: descriptor                           => simpleDescriptor
     procedure :: value                                => simpleValue
     procedure :: logarithmicDerivativeExpansionFactor => simpleLogarithmicDerivativeExpansionFactor
     procedure :: retabulate                           => simpleRetabulate
  end type linearGrowthSimple

  interface linearGrowthSimple
     !% Constructors for the {\normalfont \ttfamily simple} linear growth class.
     module procedure simpleConstructorParameters
     module procedure simpleConstructorInternal
  end interface linearGrowthSimple

contains

  function simpleConstructorParameters(parameters)
    !% Constructor for the {\normalfont \ttfamily simple} linear growth class which takes a parameter set as input.
    use Input_Parameters2
    implicit none
    type (linearGrowthSimple      )                :: simpleConstructorParameters
    type (inputParameters         ), intent(inout) :: parameters
    class(cosmologyParametersClass), pointer       :: cosmologyParameters_    
    class(cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_    

    !# <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    !# <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    simpleConstructorParameters=simpleConstructorInternal(cosmologyParameters_,cosmologyFunctions_)
    return
  end function simpleConstructorParameters

  function simpleConstructorInternal(cosmologyParameters_,cosmologyFunctions_)
    !% Internal constructor for the {\normalfont \ttfamily simple} linear growth class.
    implicit none
    type (linearGrowthSimple      )                        :: simpleConstructorInternal
    class(cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_    
    class(cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_    

    simpleConstructorInternal%cosmologyParameters_ => cosmologyParameters_
    simpleConstructorInternal%cosmologyFunctions_  => cosmologyFunctions_
    simpleConstructorInternal%tableInitialized     =  .false.
    simpleConstructorInternal%tableTimeMinimum     =   1.0d0
    simpleConstructorInternal%tableTimeMaximum     =  20.0d0
    return
  end function simpleConstructorInternal

  subroutine simpleDestructor(self)
    !% Destructor for the {\normalfont \ttfamily simple} linear growth class.
    implicit none
    type (linearGrowthSimple), intent(inout) :: self
    
    !# <objectDestructor name="self%cosmologyParameters_"/>
    !# <objectDestructor name="self%cosmologyFunctions_"/>
    if (self%tableInitialized) then
       call self%growthFactor%destroy()
       deallocate(self%growthFactor)
    end if
    return
  end subroutine simpleDestructor

  subroutine simpleRetabulate(self,time)
    !% Returns the linear growth factor $D(a)$ for expansion factor {\normalfont \ttfamily aExpansion}, normalized such that
    !% $D(1)=1$ for a simple matter plus cosmological constant cosmology.
    use Tables
    use Table_Labels
    use ODE_Solver
    implicit none
    class           (linearGrowthSimple      ), intent(inout) :: self
    double precision                          , intent(in   ) :: time
    double precision                          , parameter     :: dominateFactor               =   1.0d+04
    double precision                          , parameter     :: odeToleranceAbsolute         =   1.0d-10, odeToleranceRelative     =1.0d-10
    integer                                   , parameter     :: growthTablePointsPerDecade   =1000
    double precision                          , dimension(2)  :: growthFactorODEVariables
    logical                                                   :: remakeTable
    integer                                                   :: i
    double precision                                          :: expansionFactorMatterDominant           , growthFactorDerivative           , &
         &                                                       timeNow                                 , linearGrowthFactorPresent        , &
         &                                                       timeMatterDominant                      , timePresent
    integer                                                   :: growthTableNumberPoints    
    type            (fgsl_odeiv_step         )                :: odeStepper
    type            (fgsl_odeiv_control      )                :: odeController
    type            (fgsl_odeiv_evolve       )                :: odeEvolver
    type            (fgsl_odeiv_system       )                :: odeSystem
    logical                                                   :: odeReset                     =.true.
    
    ! Check if we need to recompute our table.
    if (self%tableInitialized) then
       remakeTable=(                              &
            &        time < self%tableTimeMinimum &
            &       .or.                          &
            &        time > self%tableTimeMaximum &
            &      )
    else
       remakeTable=.true.
    end if
    if (remakeTable) then
       ! Find the present-day epoch.
       timePresent                 =     self%cosmologyFunctions_%cosmicTime           (                        1.0d0)
       ! Find epoch of matter-dark energy equality.
       expansionFactorMatterDominant=min(                                                                               &
            &                            self%cosmologyFunctions_%expansionFactor      (        self%tableTimeMinimum), &
            &                            self%cosmologyFunctions_%dominationEpochMatter(               dominateFactor)  &
            &                           )
       timeMatterDominant           =    self%cosmologyFunctions_%cosmicTime           (expansionFactorMatterDominant)       
       ! Find minimum and maximum times to tabulate.
       self%tableTimeMinimum=min(self%tableTimeMinimum,min(timePresent,min(time/2.0,timeMatterDominant)      ))
       self%tableTimeMaximum=max(self%tableTimeMaximum,max(timePresent,max(time    ,timeMatterDominant)*2.0d0))
       ! Determine number of points to tabulate.
       growthTableNumberPoints=int(log10(self%tableTimeMaximum/self%tableTimeMinimum)*dble(growthTablePointsPerDecade))       
       ! Destroy current table.
       if (allocated(self%growthFactor)) then
          call self%growthFactor%destroy()
          deallocate(self%growthFactor)
       end if
       ! Create table.
       allocate(table1DLogarithmicLinear :: self%growthFactor)
       select type (growthFactor => self%growthFactor)
       type is (table1DLogarithmicLinear)
          call growthFactor%create(self%tableTimeMinimum,self%tableTimeMaximum,growthTableNumberPoints)
          ! Solve ODE to get corresponding expansion factors. Initialize with solution for matter dominated phase.
          call growthFactor%populate(1.0d0,1)
          growthFactorDerivative=self%cosmologyFunctions_ %expansionRate  (                   &
               &                  self%cosmologyFunctions_%expansionFactor (                  &
               &                                                            growthFactor%x(1) &
               &                                                           )                  &
               &                                                          )
          do i=2,growthTableNumberPoints
             timeNow                    =growthFactor          %x(i-1)
             growthFactorODEVariables(1)=growthFactor          %y(i-1)
             growthFactorODEVariables(2)=growthFactorDerivative
             call ODE_Solve(                                   &
                  &         odeStepper                       , &
                  &         odeController                    , &
                  &         odeEvolver                       , &
                  &         odeSystem                        , &
                  &         timeNow                          , &
                  &         growthFactor%x(i)                , &
                  &         2                                , &
                  &         growthFactorODEVariables         , &
                  &         growthFactorODEs                 , &
                  &         odeToleranceAbsolute             , &
                  &         odeToleranceRelative             , &
                  &         reset                   =odeReset  &
                  &        )
             call growthFactor%populate(growthFactorODEVariables(1),i)
             growthFactorDerivative=growthFactorODEVariables(2)
          end do
          call ODE_Solver_Free(odeStepper,odeController,odeEvolver,odeSystem)
          ! Normalize to growth factor of unity at present day.
          linearGrowthFactorPresent=growthFactor%interpolate(timePresent)
          call growthFactor%populate(reshape(growthFactor%ys(),[growthTableNumberPoints])/linearGrowthFactorPresent)
          ! Compute relative normalization factor such that growth factor behaves as expansion factor at early times.
          self%normalizationMatterDominated=+(                                                           &
               &                              +9.0d0                                                     &
               &                              *self%cosmologyParameters_%OmegaMatter   (               ) &
               &                              /4.0d0                                                     &
               &                             )**(1.0d0/3.0d0)                                            &
               &                            *(                                                           &
               &                              +self%cosmologyParameters_%HubbleConstant(hubbleUnitsTime) &
               &                              *growthFactor             %x             (              1) &
               &                             )**(2.0d0/3.0d0)                                            &
               &                            /  growthFactor             %y             (              1)
          self%tableInitialized=.true.
       end select
    end if
    return
    
  contains
    
    integer function growthFactorODEs(time,values,derivatives)
      !% System of differential equations to solve for the growth factor.
      double precision              , intent(in   ) :: time
      double precision, dimension(:), intent(in   ) :: values
      double precision, dimension(:), intent(  out) :: derivatives
      double precision                              :: expansionFactor
      
      expansionFactor   =+self%cosmologyFunctions_%expansionFactor   (                           time)
      derivatives    (1)=+values                                     (                              2)
      derivatives    (2)=+1.5d0                                                                           &
           &             *self%cosmologyFunctions_%expansionRate     (                expansionFactor)**2 &
           &             *self%cosmologyFunctions_%omegaMatterEpochal(expansionFactor=expansionFactor)    &
           &             *values                                     (                              1)    &
           &             -2.0d0                                                                           &
           &             *self%cosmologyFunctions_%expansionRate     (                expansionFactor)    &
           &             *values                                     (                              2)
      growthFactorODEs  = FGSL_Success
    end function growthFactorODEs

  end subroutine simpleRetabulate

  double precision function simpleValue(self,time,expansionFactor,collapsing,normalize,component,wavenumber)
    !% Return the linear growth factor at the given epoch.
    implicit none
    class           (linearGrowthSimple     ), intent(inout)           :: self
    double precision                         , intent(in   ), optional :: time               , expansionFactor
    logical                                  , intent(in   ), optional :: collapsing
    integer                                  , intent(in   ), optional :: normalize          , component
    double precision                         , intent(in   ), optional :: wavenumber
    double precision                                                   :: time_
    !# <optionalArgument name="normalize" defaultsTo="normalizePresentDay" />
    !GCC$ attributes unused :: component, wavenumber

    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    simpleValue=self%growthFactor%interpolate(time_)
    ! Normalize.
    select case (normalize_)
    case (normalizeMatterDominated)
       simpleValue=simpleValue*self%normalizationMatterDominated
    end select
    return
  end function simpleValue

  double precision function simpleLogarithmicDerivativeExpansionFactor(self,time,expansionFactor,collapsing,component,wavenumber)
    !% Return the logarithmic gradient of linear growth factor with respect to expansion factor at the given epoch.
    use Galacticus_Error
    implicit none
    class           (linearGrowthSimple     ), intent(inout)           :: self
    double precision                         , intent(in   ), optional :: time               , expansionFactor
    logical                                  , intent(in   ), optional :: collapsing
    integer                                  , intent(in   ), optional :: component 
    double precision                         , intent(in   ), optional :: wavenumber 
    double precision                                                   :: time_              , expansionFactor_
    !GCC$ attributes unused :: component, wavenumber
    
    ! Determine cosmological time.
    call self%cosmologyFunctions_%epochValidate(time,expansionFactor,collapsing,timeOut=time_,expansionFactorOut=expansionFactor_)
    ! Remake the table if necessary.
    call self%retabulate(time_)
    ! Interpolate to get the expansion factor.
    simpleLogarithmicDerivativeExpansionFactor=+self%growthFactor       %interpolateGradient(time_           ) &
         &                                     /self%growthFactor       %interpolate        (time_           ) &
         &                                     /self%cosmologyFunctions_%expansionRate      (expansionFactor_)
    return
  end function simpleLogarithmicDerivativeExpansionFactor

  subroutine simpleStateStore(self,stateFile,fgslStateFile)
    !% Store the tabulation.
    use Galacticus_Display
    use FGSL
    implicit none
    class  (linearGrowthSimple), intent(inout) :: self
    integer                    , intent(in   ) :: stateFile
    type   (fgsl_file         ), intent(in   ) :: fgslStateFile
    !GCC$ attributes unused :: fgslStateFile
    
    call Galacticus_Display_Message('Storing state for: linearGrowth -> simple',verbosity=verbosityInfo)
    write (stateFile) self%tableTimeMinimum,self%tableTimeMaximum
    return
  end subroutine simpleStateStore
  
  subroutine simpleStateRestore(self,stateFile,fgslStateFile)
    !% Reset the tabulation if state is to be retrieved.
    use Galacticus_Display
    use FGSL
    implicit none
    class  (linearGrowthSimple), intent(inout) :: self
    integer                    , intent(in   ) :: stateFile
    type   (fgsl_file         ), intent(in   ) :: fgslStateFile
    !GCC$ attributes unused :: fgslStateFile

    call Galacticus_Display_Message('Retrieving state for: linearGrowth -> simple',verbosity=verbosityInfo)
    read (stateFile) self%tableTimeMinimum,self%tableTimeMaximum
    self%tableInitialized=.false.
    call self%retabulate(sqrt(self%tableTimeMinimum*self%tableTimeMaximum))
    return
  end subroutine simpleStateRestore

  subroutine simpleDescriptor(self,descriptor)
    !% Add parameters to an input parameter list descriptor which could be used to recreate this object.
    use Input_Parameters2
    use FoX_DOM
    implicit none
    class(linearGrowthSimple), intent(inout) :: self
    type (inputParameters   ), intent(inout) :: descriptor
    type (inputParameters   )                :: subParameters
    
    call descriptor%addParameter("linearGrowthMethod","simple")
    subParameters=descriptor%subparameters("linearGrowthMethod")
    call self%cosmologyParameters_%descriptor(subParameters)
    call self%cosmologyFunctions_ %descriptor(subParameters)
    return
  end subroutine simpleDescriptor
