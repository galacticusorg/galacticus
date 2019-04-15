!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

!% Contains a module which handles outputting of velocity dispersion data to the \glc\ output file.

module Galacticus_Output_Trees_Velocity_Dispersion
  !% Handles outputting of velocity dispersion data to the \glc\ output file.
  use ISO_Varying_String
  use Galactic_Structure_Radii_Definitions
  use Galacticus_Nodes                    , only : treeNode
  implicit none
  private
  public :: Galacticus_Output_Tree_Velocity_Dispersion, Galacticus_Output_Tree_Velocity_Dispersion_Property_Count,&
       & Galacticus_Output_Tree_Velocity_Dispersion_Names

  ! Flag indicating whether or not velocity dispersion information is to be output.
  logical :: outputVelocityDispersionData

  ! Flag indicating whether or not this module has been initialized.
  logical :: outputVelocityDispersionDataInitialized=.false.

  ! Radii at which velocity dispersion is to be output.
  integer                                                      :: radiiCount
  type            (radiusSpecifier), allocatable, dimension(:) :: radii
  type            (varying_string ), allocatable, dimension(:) :: outputVelocityDispersionRadii
  logical                                                      :: darkMatterScaleRadiusIsNeeded, diskIsNeeded        , &
       &                                                          spheroidIsNeeded             , virialRadiusIsNeeded

  ! Module scope variables used in integrations.
  type            (treeNode       ), pointer                   :: activeNode
  integer                                                      :: componentType                , massType            , &
       &                                                          weightBy                     , weightIndex
  double precision                                             :: radiusImpact                 , radiusOuter
  !$omp threadprivate(activeNode,radiusOuter,massType,componentType,weightBy,weightIndex,radiusImpact)
  
contains

  subroutine Galacticus_Output_Tree_Velocity_Dispersion_Initialize
    !% Initializes the module by determining whether or not velocity dispersion should be output.
    use Input_Parameters
    use Galacticus_Error
    implicit none

    if (.not.outputVelocityDispersionDataInitialized) then
       !$omp critical(Galacticus_Output_Tree_Velocity_Dispersion_Initialize)
       if (.not.outputVelocityDispersionDataInitialized) then
          !# <inputParameter>
          !#   <name>outputVelocityDispersionData</name>
          !#   <cardinality>1</cardinality>
          !#   <defaultValue>.false.</defaultValue>
          !#   <description>Specifies whether or not velocity dispersion data should be incldued in the output file.</description>
          !#   <group>output</group>
          !#   <source>globalParameters</source>
          !#   <type>boolean</type>
          !# </inputParameter>
          ! Read and parse parameter controlling radii of output.
          if (outputVelocityDispersionData) then
             if (.not.globalParameters%isPresent('outputVelocityDispersionRadii'))                                                  &
                  & call Galacticus_Error_Report(                                                                                   &
                  &                              'outputVelocityDispersionRadii must be specified for velocity dispersion output'// &
                  &                              {introspection:location}                                                           &
                  &                             )
             radiiCount=globalParameters%count('outputVelocityDispersionRadii')
             allocate(outputVelocityDispersionRadii(radiiCount))
             !# <inputParameter>
             !#   <name>outputVelocityDispersionRadii</name>
             !#   <cardinality>1..*</cardinality>
             !#   <description>Specifies the radii at which the velocity dispersion will be output.</description>
             !#   <group>output</group>
             !#   <source>globalParameters</source>
             !#   <type>string</type>
             !# </inputParameter>
             ! Parse the radii definitions.
             call Galactic_Structure_Radii_Definition_Decode(outputVelocityDispersionRadii,radii,diskIsNeeded,spheroidIsNeeded,virialRadiusIsNeeded,darkMatterScaleRadiusIsNeeded)
             deallocate(outputVelocityDispersionRadii)
          end if
          ! Flag that module is now initialized.
          outputVelocityDispersionDataInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Tree_Velocity_Dispersion_Initialize)
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion_Initialize

  !# <mergerTreeOutputNames>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Dispersion_Names</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Dispersion</sortName>
  !# </mergerTreeOutputNames>
  subroutine Galacticus_Output_Tree_Velocity_Dispersion_Names(thisNode,integerProperty,integerPropertyNames&
       &,integerPropertyComments,integerPropertyUnitsSI,doubleProperty ,doublePropertyNames,doublePropertyComments&
       &,doublePropertyUnitsSI,time)
    !% Set the names of velocity dispersion properties to be written to the \glc\ output file.
    use Numerical_Constants_Astronomical
    implicit none
    type            (treeNode)              , intent(inout) :: thisNode
    double precision                        , intent(in   ) :: time
    integer                                 , intent(inout) :: doubleProperty         , integerProperty
    character       (len=*   ), dimension(:), intent(inout) :: doublePropertyComments , doublePropertyNames   , &
         &                                                     integerPropertyComments, integerPropertyNames
    double precision          , dimension(:), intent(inout) :: doublePropertyUnitsSI  , integerPropertyUnitsSI
    integer                                                 :: i
    !GCC$ attributes unused :: thisNode, time, integerProperty, integerPropertyNames, integerPropertyComments, integerPropertyUnitsSI
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Dispersion_Initialize

    ! Return property names if we are outputting velocity dispersion data.
    if (outputVelocityDispersionData) then
       do i=1,radiiCount
          doubleProperty=doubleProperty+1
          doublePropertyNames   (doubleProperty)='velocityDispersion:'//char(radii(i)%name)
          doublePropertyComments(doubleProperty)='Velocity dispersion at a given radius'
          doublePropertyUnitsSI (doubleProperty)=kilo
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion_Names

  !# <mergerTreeOutputPropertyCount>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Dispersion_Property_Count</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Dispersion</sortName>
  !# </mergerTreeOutputPropertyCount>
  subroutine Galacticus_Output_Tree_Velocity_Dispersion_Property_Count(thisNode,integerPropertyCount,doublePropertyCount,time)
    !% Account for the number of velocity dispersion properties to be written to the \glc\ output file.
    implicit none
    type            (treeNode), intent(inout) :: thisNode
    double precision          , intent(in   ) :: time
    integer                   , intent(inout) :: doublePropertyCount, integerPropertyCount
    !GCC$ attributes unused :: thisNode, integerPropertyCount, time
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Dispersion_Initialize

    ! Increment property count if we are outputting velocity dispersion data.
    if (outputVelocityDispersionData) doublePropertyCount=doublePropertyCount+radiiCount
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion_Property_Count

  !# <mergerTreeOutputTask>
  !#  <unitName>Galacticus_Output_Tree_Velocity_Dispersion</unitName>
  !#  <sortName>Galacticus_Output_Tree_Velocity_Dispersion</sortName>
  !# </mergerTreeOutputTask>
  subroutine Galacticus_Output_Tree_Velocity_Dispersion(thisNode,integerProperty,integerBufferCount,integerBuffer,doubleProperty&
       &,doubleBufferCount,doubleBuffer,time,instance)
    !% Store velocity dispersion properties in the \glc\ output file buffers.
    use Kind_Numbers
    use Galactic_Structure_Velocity_Dispersions
    use Dark_Matter_Halo_Scales
    use FGSL                                   , only : fgsl_function    , fgsl_integration_workspace
    use Numerical_Integration
    use Galactic_Structure_Options
    use Galactic_Structure_Enclosed_Masses
    use Multi_Counters
    use Galacticus_Nodes                       , only : nodeComponentDisk, nodeComponentSpheroid     , nodeComponentDarkMatterProfile
    implicit none
    double precision                                , intent(in   )         :: time
    type            (treeNode                      ), intent(inout), target :: thisNode
    integer                                         , intent(inout)         :: doubleBufferCount                , doubleProperty          , &
         &                                                                     integerBufferCount               , integerProperty
    integer         (kind=kind_int8                ), intent(inout)         :: integerBuffer        (:,:)
    double precision                                , intent(inout)         :: doubleBuffer         (:,:)
    type            (multiCounter                  ), intent(inout)         :: instance
    class           (nodeComponentDisk             ), pointer               :: thisDisk
    class           (nodeComponentSpheroid         ), pointer               :: thisSpheroid
    class           (nodeComponentDarkMatterProfile), pointer               :: thisDarkMatterProfile
    class           (darkMatterHaloScaleClass      ), pointer               :: darkMatterHaloScale_
    double precision                                , parameter             :: outerRadiusMultiplier     =10.0d0
    type            (fgsl_function                 )                        :: integrandFunction
    type            (fgsl_integration_workspace    )                        :: integrationWorkspace
    integer                                                                 :: i
    logical                                                                 :: scaleIsZero
    double precision                                                        :: densityIntegrand                 , radius                  , &
         &                                                                     radiusFromFraction               , radiusVirial            , &
         &                                                                     radiusZero                       , velocityDensityIntegrand, &
         &                                                                     numerator                        , denominator             , &
         &                                                                     massDisk                         , massSpheroid
    !GCC$ attributes unused :: time, integerProperty, integerBufferCount, integerBuffer, instance
    
    ! Initialize the module.
    call Galacticus_Output_Tree_Velocity_Dispersion_Initialize
    ! Store property data if we are outputting velocity dispersion data.
    if (outputVelocityDispersionData) then
       ! Compute required quantities.
       darkMatterHaloScale_ => darkMatterHaloScale()
       radiusVirial         =  0.0d0
       if (         virialRadiusIsNeeded) radiusVirial          =  darkMatterHaloScale_%virialRadius(thisNode                    )
       if (                 diskIsNeeded) thisDisk              =>                                   thisNode%disk             ()
       if (             spheroidIsNeeded) thisSpheroid          =>                                   thisNode%spheroid         ()
       if (darkMatterScaleRadiusIsNeeded) thisDarkMatterProfile =>                                   thisNode%darkMatterProfile()
       do i=1,radiiCount
          ! Find the radius.
          scaleIsZero=.false.
          radius     =radii(i)%value
          select case (radii(i)%type)
          case   (radiusTypeRadius                )
             radiusOuter=    radius                                        *outerRadiusMultiplier
          case   (radiusTypeVirialRadius          )
             radius     =    radius*radiusVirial
             radiusOuter=max(radius,radiusVirial                          )*outerRadiusMultiplier
          case   (radiusTypeDarkMatterScaleRadius )
             radius     =    radius*thisDarkMatterProfile%         scale()
             radiusOuter=max(radius,thisDarkMatterProfile%         scale())*outerRadiusMultiplier
          case   (radiusTypeDiskRadius            )
             radius     =    radius*thisDisk             %        radius()
             radiusOuter=max(radius,thisDisk             %        radius())*outerRadiusMultiplier
             scaleIsZero=(thisDisk                       %        radius() <= 0.0d0)
          case   (radiusTypeSpheroidRadius        )
             radius     =    radius*thisSpheroid         %        radius()
             radiusOuter=max(radius,thisSpheroid         %        radius())*outerRadiusMultiplier
             scaleIsZero=(thisSpheroid                   %        radius() <= 0.0d0)
          case   (radiusTypeDiskHalfMassRadius    )
             radius     =    radius*thisDisk             %halfMassRadius()
             radiusOuter=max(radius,thisDisk             %halfMassRadius())*outerRadiusMultiplier
             scaleIsZero=(thisDisk                       %halfMassRadius() <= 0.0d0)
          case   (radiusTypeSpheroidHalfMassRadius)
             radius     =    radius*thisSpheroid         %halfMassRadius()
             radiusOuter=max(radius,thisSpheroid         %halfMassRadius())*outerRadiusMultiplier
             scaleIsZero=(thisSpheroid                   %halfMassRadius() <= 0.0d0)
          case   (radiusTypeGalacticMassFraction ,  &
               &  radiusTypeGalacticLightFraction )
             radiusFromFraction=                              &
                  &  Galactic_Structure_Radius_Enclosing_Mass &
                  &  (                                        &
                  &   thisNode                             ,  &
                  &   fractionalMass=radii(i)%fraction     ,  &
                  &   massType      =massTypeGalactic      ,  &
                  &   componentType =componentTypeAll      ,  &
                  &   weightBy      =radii(i)%weightBy     ,  &
                  &   weightIndex   =radii(i)%weightByIndex   &
                  &  )
             radius     =    radius*radiusFromFraction
             radiusOuter=max(radius,radiusFromFraction)*outerRadiusMultiplier
          end select
          ! Store the rotation curve.
          doubleProperty=doubleProperty+1
          if (scaleIsZero) then
             ! Do not compute dispersions if the component scale is zero.
             doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
          else
             select case (radii(i)%direction)
             case (directionRadial                    )
                ! Radial velocity dispersion.
                doubleBuffer(doubleBufferCount,doubleProperty)=                                                                   &
                     &                                   Galactic_Structure_Velocity_Dispersion(                                  &
                     &                                                                          thisNode                        , &
                     &                                                                          radius                          , &
                     &                                                                          radiusOuter                     , &
                     &                                                                          massType     =radii(i)%mass     , &
                     &                                                                          componentType=radii(i)%component  &
                     &                                                                         )
             case (directionLineOfSight               )
                ! Line-of-sight velocity dispersion.
                activeNode   => thisNode
                massType     =  radii(i)%mass
                componentType=  radii(i)%component
                weightBy     =  radii(i)%integralWeightBy
                weightIndex  =  radii(i)%integralWeightByIndex
                radiusImpact =  radius
                doubleBuffer(doubleBufferCount,doubleProperty)=Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion(radius)
             case (directionLineOfSightInteriorAverage)
                ! Average over the line-of-sight velocity dispersion within the radius.
                activeNode   => thisNode
                massType     =  radii(i)%mass
                componentType=  radii(i)%component
                weightBy     =  radii(i)%integralWeightBy
                weightIndex  =  radii(i)%integralWeightByIndex
                radiusZero   =  0.0d0
                radiusImpact =  radius
                velocityDensityIntegrand=Integrate(radiusZero,radiusOuter&
                     &,Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd ,integrandFunction&
                     &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                densityIntegrand        =Integrate(radiusZero,radiusOuter&
                     &,Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd ,integrandFunction&
                     &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
                call Integrate_Done(integrandFunction,integrationWorkspace)
                if (velocityDensityIntegrand <= 0.0d0) then
                   doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
                else
                   doubleBuffer(doubleBufferCount,doubleProperty)=sqrt(velocityDensityIntegrand/densityIntegrand)
                end if
             case (directionLambdaR                   )
                ! The "lambdaR" parameter of Cappellari et al. (2007; MNRAS; 379; 418)
                activeNode   => thisNode
                massType     =  radii(i)%mass
                componentType=  radii(i)%component
                weightBy     =  radii(i)%integralWeightBy
                weightIndex  =  radii(i)%integralWeightByIndex
                ! Check the total masses of the disk and spheroid components. If either is zero we can use the solutions for the
                ! appropriate limiting case.
                massSpheroid=Galactic_Structure_Enclosed_Mass(                                 &
                     &                                        thisNode                       , &
                     &                                        radiusLarge                    , &
                     &                                        massType     =massType         , &
                     &                                        componentType=componentType    , &
                     &                                        weightBy     =weightBy         , &
                     &                                        weightIndex  =weightIndex        &
                     &                                       )
                massDisk    =Galactic_Structure_Enclosed_Mass(                                 &
                     &                                        thisNode                       , &
                     &                                        radiusLarge                    , &
                     &                                        massType     =massTypeStellar  , &
                     &                                        componentType=componentTypeDisk, &
                     &                                        weightBy     =weightBy         , &
                     &                                        weightIndex  =weightIndex        &
                     &                                       )                             
                if (massDisk <= 0.0d0) then
                   doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
                else if (massSpheroid <= 0.0d0) then
                   doubleBuffer(doubleBufferCount,doubleProperty)=1.0d0
                else
                   ! Full calculation is required.
                   radiusZero=0.0d0
                   numerator=Integrate(                                                                 &
                        &                radiusZero                                                   , &
                        &                radius                                                       , &
                        &                Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd2        , &
                        &                integrandFunction                                            , &
                        &                integrationWorkspace                                         , &
                        &                toleranceAbsolute                                     =0.0d0 , &
                        &                toleranceRelative                                     =1.0d-2  &
                        &               )
                   call Integrate_Done(integrandFunction,integrationWorkspace)
                   denominator=Integrate(                                                               &
                        &                radiusZero                                                   , &
                        &                radius                                                       , &
                        &                Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd1        , &
                        &                integrandFunction                                            , &
                        &                integrationWorkspace                                         , &
                        &                toleranceAbsolute                                     =0.0d0 , &
                        &                toleranceRelative                                     =1.0d-2  &
                        &               )
                   call Integrate_Done(integrandFunction,integrationWorkspace)                      
                   if (denominator <= 0.0d0) then
                      doubleBuffer(doubleBufferCount,doubleProperty)=0.0d0
                   else
                      doubleBuffer(doubleBufferCount,doubleProperty)=numerator/denominator
                   end if
                end if
             end select
          end if
       end do
    end if
    return
  end subroutine Galacticus_Output_Tree_Velocity_Dispersion

  double precision function Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd1(radius)
    !% Integrand function used for integrating the $\lambda_\mathrm{R}$ statistic of \cite{cappellari_sauron_2007}. In this case we
    !% want to evaluate
    !% \begin{equation}
    !% \int_0^r 2 \pi r^\prime \Sigma(r^\prime) \sqrt{\sigma^2(r^\prime)+V^2(r^\prime)} \mathrm{d}r^\prime,
    !% \end{equation}
    !% where $\Sigma(r)$ is the projected surface density (in mass or light) of the galaxy at radius $r$, $\sigma^2(r)$ is the
    !% measured velocity dispersion and $V(r)$ the measured rotation speed. Assuming that the selected component is purely
    !% dispersion dominated with velocity dispersion $\sigma_\mathrm{s}(r)$, and that rotation is present in only the disk component
    !% with rotation curve $V_\mathrm{d}(r)$ then we can model the velocity distribution, $P(V)$, at $r$ as the sum of a Gaussian of
    !% width $\sigma_\mathrm{s}(r)$ and normalized area $\Sigma_\mathrm{s}(r)$, and a delta function at $V_\mathrm{d}(r)$ with normalized
    !% area $\Sigma_\mathrm{d}(r)$. The measured rotation speed is then:
    !% \begin{equation}
    !% V(r) = \left. \int_{-\infty}^{+\infty} P(V) V \mathrm{d}V \right/ \int_{-\infty}^{+\infty} P(V) \mathrm{d}V = {\Sigma_\mathrm{d}(r) V_\mathrm{d}(r) \over [\Sigma_\mathrm{d}(r)+\Sigma_\mathrm{s}(r)]},
    !% \end{equation}
    !% and the measured velocity dispersion is:
    !% \begin{equation}
    !% \sigma^2(r) = \left. \int_{-\infty}^{+\infty} P(V) [V-V(r)]^2 \mathrm{d}V \right/ \int_{-\infty}^{+\infty} P(V) \mathrm{d}V = {  \Sigma_\mathrm{s}(r) [\sigma_\mathrm{s}^2(r)] + \Sigma_\mathrm{d}(r) [V_\mathrm{d}(r)-V(r)]^2  \over [\Sigma_\mathrm{d}(r)+\Sigma_\mathrm{s}(r)]}.
    !% \end{equation}
    use Galactic_Structure_Options
    use Galactic_Structure_Densities
    use Galactic_Structure_Velocity_Dispersions
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Surface_Densities
    use Numerical_Constants_Math
    use Numerical_Integration
    use FGSL                                   , only : fgsl_function, fgsl_integration_workspace
    implicit none
    double precision                            , intent(in   ) :: radius
    double precision                            , parameter     :: fractionSmall=1.0d-3
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: sigmaLineOfSightSquaredSpheroidDensity, densitySpheroid        , &
         &                                                         densityDisk                           , velocityDisk           , &
         &                                                         velocityMean                          , sigmaLineOfSightSquared

    if (radius <= 0.0d0) then
       Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd1=0.0d0
    else
       radiusImpact=radius
       densitySpheroid                                                                         &
            & =Integrate(                                                                      &
            &            radius                                                              , &
            &            radiusOuter                                                         , &
            &            Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand       , &
            &            integrandFunction                                                   , &
            &            integrationWorkspace                                                , &
            &            toleranceAbsolute                                            =0.0d0 , &
            &            toleranceRelative                                            =1.0d-2  &
            &           )
       call Integrate_Done(integrandFunction,integrationWorkspace)       
       densityDisk=Galactic_Structure_Surface_Density            (                                        &
            &                                                     activeNode                            , &
            &                                                     [radius,0.0d0,0.0d0]                  , &
            &                                                     massType            =massTypeStellar  , &
            &                                                     componentType       =componentTypeDisk, &
            &                                                     weightBy            =weightBy         , &
            &                                                     weightIndex         =weightIndex        &
            &                                                    )                     
       velocityDisk=Galactic_Structure_Rotation_Curve            (                                        &
            &                                                     activeNode                            , &
            &                                                     radius                                , &
            &                                                     massType            =massTypeAll      , &
            &                                                     componentType       =componentTypeAll   &
            &                                                    )
       ! Test if the spheroid density is significant....
       if (densitySpheroid < fractionSmall*densityDisk) then
          ! ...it is not, so we can avoid computing the spheroid velocity dispersion.
          Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd1 &
               & =2.0d0                                         &
               & *Pi                                            &
               & *radius                                        &
               & *densityDisk                                   &
               & *velocityDisk
       else
          ! ...it is, so we must do the full calculation.
          sigmaLineOfSightSquaredSpheroidDensity                                                  &
               & =Integrate(                                                                      &
               &            radius                                                              , &
               &            radiusOuter                                                         , &
               &            Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd        , &
               &            integrandFunction                                                   , &
               &            integrationWorkspace                                                , &
               &            toleranceAbsolute                                            =0.0d0 , &
               &            toleranceRelative                                            =1.0d-2  &
               &           )
          call Integrate_Done(integrandFunction,integrationWorkspace)
          velocityMean                                          &
               & =densityDisk                                   &
               & *velocityDisk                                  &
               & /(densityDisk+densitySpheroid)
          sigmaLineOfSightSquared                               &
               & =(                                             &
               &    sigmaLineOfSightSquaredSpheroidDensity      &
               &   +densitySpheroid                             &
               &   *velocityMean**2                             &
               &   +densityDisk                                 &
               &   *(velocityDisk-velocityMean)**2              &
               &  )                                             &
               & /(densityDisk+densitySpheroid)
          Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd1 &
               & =2.0d0                                         &
               & *Pi                                            &
               & *radius                                        &
               & *(densityDisk+densitySpheroid)                 &
               & *sqrt(                                         &
               &        sigmaLineOfSightSquared                 &
               &       +velocityMean**2                         &
               &      )
       end if
    end if
    return
  end function Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd1
  
  double precision function Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd2(radius)
    !% Integrand function used for integrating the $\lambda_\mathrm{R}$ statistic of \cite{cappellari_sauron_2007}. In this case we
    !% want to evaluate
    !% \begin{equation}
    !% \int_0^r 2 \pi r^\prime \Sigma(r^\prime) V(r^\prime) \mathrm{d}r^\prime,
    !% \end{equation}
    !% where $\Sigma(r)$ is the projected surface density (in mass or light) of the galaxy at radius $r$, and $V(r)$ the measured
    !% rotation speed. Assuming that the selected component is purely dispersion dominated with velocity dispersion $\sigma_\mathrm{
    !% s}(r)$, and that rotation is present in only the disk component with rotation curve $V_\mathrm{d}(r)$ then we can model the
    !% velocity distribution, $P(V)$, at $r$ as the sum of a Gaussian of width $\sigma_\mathrm{s}(r)$ and normalized area
    !% $\Sigma_\mathrm{s}(r)$, and a delta function at $V_\mathrm{d}(r)$ with normalized area $\Sigma_\mathrm{d}(r)$. The measured rotation
    !% speed is then:
    !% \begin{equation}
    !% V(r) = \left. \int_{-\infty}^{+\infty} P(V) V \mathrm{d}V \right/ \int_{-\infty}^{+\infty} P(V) \mathrm{d}V = {\Sigma_\mathrm{d}(r) V_\mathrm{d}(r) \over [\Sigma_\mathrm{d}(r)+\Sigma_\mathrm{s}(r)]}.
    !% \end{equation}
    use Galactic_Structure_Options
    use Galactic_Structure_Rotation_Curves
    use Galactic_Structure_Surface_Densities
    use Numerical_Integration
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: densityDisk, velocityDisk

    if (radius <= 0.0d0) then
       Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd2=0.0d0
    else
       densityDisk=Galactic_Structure_Surface_Density(                                        &
            &                                         activeNode                            , &
            &                                         [radius,0.0d0,0.0d0]                  , &
            &                                         massType            =massTypeStellar  , &
            &                                         componentType       =componentTypeDisk, &
            &                                         weightBy            =weightBy         , &
            &                                         weightIndex         =weightIndex        &
            &                                        )                     
       velocityDisk=Galactic_Structure_Rotation_Curve(                                        &
            &                                         activeNode                            , &
            &                                         radius                                , &
            &                                         massType            =massTypeAll      , &
            &                                         componentType       =componentTypeAll   &
            &                                        )           
       Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd2=2.0d0*Pi*radius*densityDisk*velocityDisk
    end if
    return
  end function Galacticus_Output_Trees_Vlcty_Dsprsn_LambdaR_Intgrnd2

  double precision function Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd(radius)
    !% Integrand function used for integrating line-of-sight velocity dispersion over surface density.
    use Galactic_Structure_Densities
    use Galactic_Structure_Velocity_Dispersions
    implicit none
    double precision, intent(in   ) :: radius

    if (radius <= 0.0d0) then
       Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd=0.0d0
    else
       Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd=                                    &
            &                        Spherical_Shell_Solid_Angle_In_Cylcinder(radius)                    &
            &                       *                                         radius **2                 &
            &                       *Galactic_Structure_Density            (                             &
            &                                                               activeNode                 , &
            &                                                               [radius,0.0d0,0.0d0]       , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               weightBy     =weightBy     , &
            &                                                               weightIndex  =weightIndex    &
            &                                                              )                             &
            &                       *Galactic_Structure_Velocity_Dispersion(                             &
            &                                                               activeNode                 , &
            &                                                               radius                     , &
            &                                                               radiusOuter                , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType  &
            &                                                              )**2
    end if
   return
  end function Galacticus_Output_Trees_Vlcty_Dsprsn_Vlcty_Dnsty_Srfc_Intgrnd

  double precision function Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd(radius)
    !% Integrand function used for integrating line-of-sight surface density dispersion over area.
    use Galactic_Structure_Densities
    implicit none
    double precision, intent(in   ) :: radius

    if (radius <= 0.0d0) then
       Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd=0.0d0
    else
       Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd=                               &
            &                         Spherical_Shell_Solid_Angle_In_Cylcinder(radius)        &
            &                        *                                         radius **2     &
            &                        *Galactic_Structure_Density(                             &
            &                                                    activeNode                 , &
            &                                                    [radius,0.0d0,0.0d0]       , &
            &                                                    massType     =massType     , &
            &                                                    componentType=componentType, &
            &                                                    weightBy     =weightBy     , &
            &                                                    weightIndex  =weightIndex    &
            &                                                   )
    end if
    return
  end function Galacticus_Output_Trees_Vlcty_Dsprsn_Dnsty_Srfc_Intgrnd

  double precision function Spherical_Shell_Solid_Angle_In_Cylcinder(radius)
    !% Computes the solid angle of a spherical shelll of given {\normalfont \ttfamily radius} that lies within a cylinder of radius {\normalfont \ttfamily
    !% radiusImpact}.
    use Numerical_Constants_Math
    implicit none
    double precision, intent(in   ) :: radius

    ! Test size of sphere relative to cylinder.
    if (radius <= radiusImpact) then
       ! Sphere is entirely within the cylinder.
       Spherical_Shell_Solid_Angle_In_Cylcinder=4.0d0*Pi
    else
       ! Some of sphere lies outside of cylinder.
       Spherical_Shell_Solid_Angle_In_Cylcinder=4.0d0*Pi*(1.0d0-cos(asin(radiusImpact/radius)))
    end if
    return
  end function Spherical_Shell_Solid_Angle_In_Cylcinder

  double precision function Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion(radius)
    !% Compute the line-of-sight velocity dispersion at the given {\normalfont \ttfamily radius}.
    use FGSL                 , only : fgsl_function, fgsl_integration_workspace
    use Numerical_Integration
    implicit none
    double precision                            , intent(in   ) :: radius
    type            (fgsl_function             )                :: integrandFunction
    type            (fgsl_integration_workspace)                :: integrationWorkspace
    double precision                                            :: densityIntegral     , velocityDensityIntegral

    velocityDensityIntegral=Integrate(radius,radiusOuter&
         &,Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd ,integrandFunction&
         &,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    densityIntegral        =Integrate(radius,radiusOuter,Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand&
         & ,integrandFunction,integrationWorkspace,toleranceAbsolute=0.0d0,toleranceRelative=1.0d-3)
    call Integrate_Done(integrandFunction,integrationWorkspace)
    if (velocityDensityIntegral <= 0.0d0) then
       Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion=0.0d0
    else
       Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion=sqrt(velocityDensityIntegral/densityIntegral)
    end if
    return
  end function Galacticus_Output_Trees_Line_of_Sight_Velocity_Dispersion

  double precision function Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand(radius)
    !% Integrand function used for computing line-of-sight velocity dispersions.
    use Galactic_Structure_Densities
    implicit none
    double precision, intent(in   ) :: radius

    if (radius == 0.0d0) then
       Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand=0.0d0
    else
       Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand=                                    &
            &                        Galactic_Structure_Density           (                              &
            &                                                               activeNode                 , &
            &                                                               [radius,0.0d0,0.0d0]       , &
            &                                                               massType     =massType     , &
            &                                                               componentType=componentType, &
            &                                                               weightBy     =weightBy     , &
            &                                                               weightIndex  =weightIndex    &
            &                                                              )                             &
            &                       *     radius                                                         &
            &                       /sqrt(radius**2-radiusImpact**2)
    end if
    return
  end function Galacticus_Output_Trees_Velocity_Dispersion_Density_Integrand

  double precision function Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd(radius)
    !% Integrand function used for computing line-of-sight velocity dispersions. Specifically, we wish to evaluate the integral:
    !% \begin{equation}
    !% \int_{r_\mathrm{i}}^{r_\mathrm{o}} \sigma^2(r) \rho(r) {r \over \sqrt{r^2-r_\mathrm{i}^2}} \mathrm{d}r,
    !% \label{eq:velocityDispersionDensityIntegral}
    !% \end{equation}
    !% where $r_\mathrm{i}$ is the impact parameter, $r_\mathrm{o}$ is an outer radius at which we assume $\rho(r_\mathrm{
    !% o})\sigma^2(r_\mathrm{o}) = 0$ (i.e. it is the radius at which we begin integrating the Jeans equation), $\rho(r)$ is density,
    !% and $\sigma(r)$ is the velocity dispersion at radius $r$. Assuming spherical symmetry and isotropic velocity dispersion,
    !% the Jeans equation tells us
    !% \begin{equation}
    !% \rho(r) \sigma^2(r) = \int^{r_\mathrm{o}}_r {\mathrm{G} M(<r^\prime) \over r^{\prime 2}} \rho(r^\prime) \mathrm{d}r^\prime,
    !% \label{eq:sphericalIsotropicJeans}
    !% \end{equation}
    !% where $\mathrm{G}$ is the gravitational constant, and $M(<r)$ is the total mass contained within radius
    !% $r$. Equation~(\ref{eq:velocityDispersionDensityIntegral}) can then be simplified using integration by parts to give:
    !% \begin{equation}
    !% \left[ \sigma^2(r)\rho(r)\sqrt{r^2-r_\mathrm{i}^2}\right]_{r_\mathrm{i}}^{r_\mathrm{o}} + \int_{r_\mathrm{i}}^{r_\mathrm{o}} {\mathrm{d}\over \mathrm{d}r} \left[ \sigma^2(r) \rho(r) \right] \sqrt{r^2-r_\mathrm{i}^2} \mathrm{d}r.
    !% \end{equation}
    !% The first term is zero at both limits (due to the constraint $\rho(r_\mathrm{o})\sigma^2(r_\mathrm{o}) = 0$ at $r_\mathrm{o}$ and
    !% due to $sqrt{r^2-r_\mathrm{i}^2}=0$ at $r_\mathrm{i}$), and the second term can be simplified using
    !% eqn.~(\ref{eq:sphericalIsotropicJeans}) to give
    !% \begin{equation}
    !% \int_{r_\mathrm{i}}^{r_\mathrm{o}} {\mathrm{G} M(<r) \over r^2} \rho(r) \sqrt{r^2-r_\mathrm{i}^2} \mathrm{d}r.
    !% \end{equation}
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    implicit none
    double precision, intent(in   ) :: radius

    if (radius <= radiusImpact) then
       Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd=0.0d0
    else
       Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd=                                        &
            &                        gravitationalConstantGalacticus                                        &
            &                       *Galactic_Structure_Density           (                                 &
            &                                                               activeNode                    , &
            &                                                               [radius,0.0d0,0.0d0]          , &
            &                                                               massType     =massType        , &
            &                                                               componentType=componentType   , &
            &                                                               weightBy     =weightBy        , &
            &                                                               weightIndex  =weightIndex       &
            &                                                              )                                &
            &                       *Galactic_Structure_Enclosed_Mass(                                      &
            &                                                               activeNode                    , &
            &                                                               radius                        , &
            &                                                               massType     =massTypeAll     , &
            &                                                               componentType=componentTypeAll  &
            &                                                              )                                &
            &                       /     radius**2                                                         &
            &                       *sqrt(radius**2-radiusImpact**2)
    end if
    return
  end function Galacticus_Output_Trees_Vlcty_Dispersion_Vlcty_Dnsty_Intgrnd

end module Galacticus_Output_Trees_Velocity_Dispersion
