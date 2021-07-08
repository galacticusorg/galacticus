!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021
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
  An accretion flow class which models the accretion flow using the model of \cite{lam_modeling_2013}.
  !!}

  use :: Cosmology_Functions            , only : cosmologyFunctionsClass
  use :: Cosmological_Density_Field     , only : criticalOverdensityClass
  use :: Dark_Matter_Halo_Scales        , only : darkMatterHaloScaleClass
  use :: Correlation_Functions_Two_Point, only : correlationFunctionTwoPointClass
  use :: Dark_Matter_Halo_Biases        , only : darkMatterHaloBiasClass
  use :: Linear_Growth                  , only : linearGrowthClass

  
  !![
  <accretionFlows name="accretionFlowsLam2013">
   <description>
    An accretion flow class using the model of \cite{lam_modeling_2013}.
   </description>
  </accretionFlows>
  !!]
  type, extends(accretionFlowsClass) :: accretionFlowsLam2013
     !!{
     An accretion flow class which models the accretion flow using the model of \cite{lam_modeling_2013}.
     !!}
     private
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_          => null()
     class           (criticalOverdensityClass        ), pointer :: criticalOverdensity_         => null()
     class           (darkMatterHaloBiasClass         ), pointer :: darkMatterHaloBias_          => null()
     class           (darkMatterHaloScaleClass        ), pointer :: darkMatterHaloScale_         => null()
     class           (correlationFunctionTwoPointClass), pointer :: correlationFunctionTwoPoint_ => null()
     class           (linearGrowthClass               ), pointer :: linearGrowth_                => null()
     double precision                                            :: scaleFactorVelocity
   contains
     final     ::             lam2013Destructor
     procedure :: density  => lam2013Density
     procedure :: velocity => lam2013Velocity
  end type accretionFlowsLam2013

  interface accretionFlowsLam2013
     !!{
     Constructors for the {\normalfont \ttfamily lam2013} accretion flows class.
     !!}
     module procedure lam2013ConstructorParameters
     module procedure lam2013ConstructorInternal
   end interface accretionFlowsLam2013

contains
  
  function lam2013ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily lam2013} accretion flow class that takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (accretionFlowsLam2013)                           :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (criticalOverdensityClass        ), pointer       :: criticalOverdensity_
    class           (darkMatterHaloBiasClass         ), pointer       :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass        ), pointer       :: darkMatterHaloScale_
    class           (correlationFunctionTwoPointClass), pointer       :: correlationFunctionTwoPoint_
    class           (linearGrowthClass               ), pointer       :: linearGrowth_
    double precision                                                  :: scaleFactorVelocity

    !![
    <inputParameter>
      <name>scaleFactorVelocity</name>
      <source>parameters</source>
      <defaultValue>1.0d0</defaultValue>
      <description>A scale factor to be applied to inflow velocities.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions"          name="cosmologyFunctions_"          source="parameters"/>
    <objectBuilder class="criticalOverdensity"         name="criticalOverdensity_"         source="parameters"/>
    <objectBuilder class="linearGrowth"                name="linearGrowth_"                source="parameters"/>
    <objectBuilder class="darkMatterHaloScale"         name="darkMatterHaloScale_"         source="parameters"/>
    <objectBuilder class="darkMatterHaloBias"          name="darkMatterHaloBias_"          source="parameters"/>
    <objectBuilder class="correlationFunctionTwoPoint" name="correlationFunctionTwoPoint_" source="parameters"/>
    !!]
    self=accretionFlowsLam2013(scaleFactorVelocity,cosmologyFunctions_,criticalOverdensity_,darkMatterHaloBias_,darkMatterHaloScale_,correlationFunctionTwoPoint_,linearGrowth_)
    !![
    <inputParametersValidate source="parameters"         />
    <objectDestructor name="darkMatterHaloScale_"        />
    <objectDestructor name="darkMatterHaloBias_"         />
    <objectDestructor name="correlationFunctionTwoPoint_"/>
    <objectDestructor name="cosmologyFunctions_"         />
    <objectDestructor name="linearGrowth_"               />
    <objectDestructor name="criticalOverdensity_"        />
    !!]
    return
  end function lam2013ConstructorParameters

  function lam2013ConstructorInternal(scaleFactorVelocity,cosmologyFunctions_,criticalOverdensity_,darkMatterHaloBias_,darkMatterHaloScale_,correlationFunctionTwoPoint_,linearGrowth_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily lam2013} accretion flows class.
    !!}
    implicit none
    type            (accretionFlowsLam2013           )                        :: self
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (criticalOverdensityClass        ), intent(in   ), target :: criticalOverdensity_
    class           (darkMatterHaloBiasClass         ), intent(in   ), target :: darkMatterHaloBias_
    class           (darkMatterHaloScaleClass        ), intent(in   ), target :: darkMatterHaloScale_
    class           (correlationFunctionTwoPointClass), intent(in   ), target :: correlationFunctionTwoPoint_
    class           (linearGrowthClass               ), intent(in   ), target :: linearGrowth_
    double precision                                  , intent(in   )         :: scaleFactorVelocity
   !![
   <constructorAssign variables="scaleFactorVelocity, *cosmologyFunctions_, *criticalOverdensity_, *darkMatterHaloBias_, *darkMatterHaloScale_, *correlationFunctionTwoPoint_, *linearGrowth_"/>
   !!]

    return
  end function lam2013ConstructorInternal

  subroutine lam2013Destructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily lam2013} accretion flows class.
    !!}
    implicit none
    type(accretionFlowsLam2013), intent(inout) :: self
    
    !![
    <objectDestructor name="self%cosmologyFunctions_"         />
    <objectDestructor name="self%criticalOverdensity_"        />
    <objectDestructor name="self%darkMatterHaloBias_"         />
    <objectDestructor name="self%darkMatterHaloScale_"        />
    <objectDestructor name="self%correlationFunctionTwoPoint_"/>
    <objectDestructor name="self%linearGrowth_"               />
    !!]
   return
  end subroutine lam2013Destructor
  
  double precision function lam2013Density(self,node,radius)
    !!{
    Compute the density of the accretion flow at the given radius.
    !!}
    use :: Galacticus_Error, only : Galacticus_Error_Report
    implicit none
    class           (accretionFlowsLam2013), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: radius
    !$GLC attributes unused :: self, node, radius
    
    lam2013Density=0.0d0
    call Galacticus_Error_Report('density is unsupported'//{introspection:location})
    return
  end function lam2013Density

  double precision function lam2013Velocity(self,node,radius)
    !!{
    Compute the mean radial velocity of the accretion flow at the given radius.
    !!}
    use :: Galacticus_Nodes        , only : nodeComponentBasic
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class           (accretionFlowsLam2013), intent(inout) :: self
    type            (treeNode             ), intent(inout) :: node
    double precision                       , intent(in   ) :: radius
    class           (nodeComponentBasic   ), pointer       :: basic
    double precision                                       :: time                    , overdensityCritical, &
         &                                                    radiusVirial            , massShell          , &
         &                                                    densityContrastNonLinear

    basic                    =>  node                      %basic               (    )
    time                     =   basic                     %time                (    )
    overdensityCritical      =   self %criticalOverdensity_%value               (time)
    radiusVirial             =   self %darkMatterHaloScale_%virialRadius        (node)
    ! Evaluate the mass in the shell outside the halo virial radius using equation (B4) of Lam et al. (2013).
    if (radius > radiusVirial) then
       massShell             =  +self %cosmologyFunctions_ %matterDensityEpochal(time)                                                                                                     &
            &                   *4.0d0                                                                                                                                                     &
            &                   *Pi                                                                                                                                                        &
            &                   /3.0d0                                                                                                                                                     &
            &                   *(                                                                                                                                                         &
            &                     +radius      **3*(1.0d0+self%darkMatterHaloBias_%bias(node,radius      )*self%correlationFunctionTwoPoint_%correlationVolumeAveraged(radius      ,time)) &
            &                     -radiusVirial**3*(1.0d0+self%darkMatterHaloBias_%bias(node,radiusVirial)*self%correlationFunctionTwoPoint_%correlationVolumeAveraged(radiusVirial,time)) &
            &                    )
    else
       massShell=0.0d0
    end if
    ! Compute the nonlinear density contrast using equation (B1) of Lam et al. (2013).
    densityContrastNonlinear =  +(                                                                                                                                                         &
         &                        +basic%mass     ()                                                                                                                                       &
         &                        +      massShell                                                                                                                                         &
         &                       )                                                                                                                                                         &
         &                      /self%cosmologyFunctions_%matterDensityEpochal                (time)                                                                                       &
         &                      /(                                                                                                                                                         &
         &                        +4.0d0                                                                                                                                                   &
         &                        *Pi                                                                                                                                                      &
         &                        /3.0d0                                                                                                                                                   &
         &                        *radius**3                                                                                                                                               &
         &                       )
    ! Evaluate the inflow velocity in the spherical collapse model using equation (B2) of Lam et al. (2013).
    lam2013Velocity          =  -self%scaleFactorVelocity                                                                                                                                  &
         &                      *self%cosmologyFunctions_%hubbleParameterEpochal              (time)                                                                                       &
         &                      *radius                                                                                                                                                    &
         &                      *self%cosmologyFunctions_%expansionFactor                     (time)                                                                                       &
         &                      *self%linearGrowth_      %logarithmicDerivativeExpansionFactor(time)                                                                                       &
         &                      /3.0d0                                                                                                                                                     &
         &                      *                                   overdensityCritical                                                                                                    &
         &                      *(                                                                                                                                                         &
         &                        +densityContrastNonLinear**(1.0d0/overdensityCritical)                                                                                                   &
         &                        -1.0d0                                                                                                                                                   &
         &                      )
    return
  end function lam2013Velocity
