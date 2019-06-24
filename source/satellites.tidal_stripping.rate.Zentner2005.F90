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

!+    Contributions to this file made by:  Anthony Pullen, Andrew Benson.

  !% Implementation of a satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.
  
  use Dark_Matter_Halo_Scales
  use Kind_Numbers  
  
  !# <satelliteTidalStripping name="satelliteTidalStrippingZentner2005">
  !#  <description>A satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.</description>
  !# </satelliteTidalStripping>
  type, extends(satelliteTidalStrippingClass) :: satelliteTidalStrippingZentner2005
     !% Implementation of a satellite tidal stripping class which follows the model of \cite{zentner_physics_2005}.
     private
     class           (darkMatterHaloScaleClass), pointer :: darkMatterHaloScale_ => null()
     double precision                                    :: efficiency          , expandMultiplier, &
          &                                                 radiusTidalPrevious
     integer         (kind_int8               )          :: lastUniqueID
   contains
     !@ <objectMethods>
     !@   <object>satelliteTidalStrippingZentner2005</object>
     !@   <objectMethod>
     !@     <method>calculationReset</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(table)\textgreater} node\arginout</arguments>
     !@     <description>Reset memoized calculations.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                     zentner2005Destructor
     procedure :: autoHook         => zentner2005AutoHook
     procedure :: calculationReset => zentner2005CalculationReset
     procedure :: massLossRate     => zentner2005MassLossRate
  end type satelliteTidalStrippingZentner2005

  interface satelliteTidalStrippingZentner2005
     !% Constructors for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
     module procedure zentner2005ConstructorParameters
     module procedure zentner2005ConstructorInternal
  end interface satelliteTidalStrippingZentner2005

  ! Module-scope objects used for root finding.
  type            (treeNode), pointer :: zentner2005Node
  double precision                    :: zentner2005TidalPull
  !$omp threadprivate(zentner2005Node,zentner2005TidalPull)
  
contains
  
  function zentner2005ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily zentner2005} satellite tidaql stripping class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (satelliteTidalStrippingZentner2005)                :: self
    type            (inputParameters                   ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass          ), pointer       :: darkMatterHaloScale_
    double precision                                                    :: efficiency

    !# <inputParameter>
    !#   <name>efficiency</name>
    !#   <cardinality>1</cardinality>
    !#   <defaultValue>2.5d0</defaultValue>
    !#   <description>The dimensionless rate coefficient apeparing in the \cite{zentner_physics_2005} expression for the tidal mass loss rate from subhalos.</description>
    !#   <group></group>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale" name="darkMatterHaloScale_" source="parameters"/>
    self=satelliteTidalStrippingZentner2005(efficiency,darkMatterHaloScale_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"/>
    return
  end function zentner2005ConstructorParameters

  function zentner2005ConstructorInternal(efficiency,darkMatterHaloScale_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
    implicit none
    type            (satelliteTidalStrippingZentner2005)                        :: self
    class           (darkMatterHaloScaleClass          ), intent(in   ), target :: darkMatterHaloScale_
    double precision                                    , intent(in)            :: efficiency
    !# <constructorAssign variables="efficiency, *darkMatterHaloScale_"/>

    self%expandMultiplier=2.0d0
    return
  end function zentner2005ConstructorInternal
  
  subroutine zentner2005AutoHook(self)
    !% Attach to the calculation reset event.
    use Events_Hooks, only : calculationResetEvent, openMPThreadBindingAllLevels
    implicit none
    class(satelliteTidalStrippingZentner2005), intent(inout) :: self

    call calculationResetEvent%attach(self,zentner2005CalculationReset,openMPThreadBindingAllLevels)
    return
  end subroutine zentner2005AutoHook
  
  subroutine zentner2005Destructor(self)
    !% Destructor for the {\normalfont \ttfamily zentner2005} satellite tidal stripping class.
    use Events_Hooks, only : calculationResetEvent
    implicit none
    type(satelliteTidalStrippingZentner2005), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"/>
    call calculationResetEvent%detach(self,zentner2005CalculationReset)
    return
  end subroutine zentner2005Destructor

  double precision function zentner2005MassLossRate(self,node)
    !% Return a mass loss rate for satellites due to tidal stripping using the formulation of \cite{zentner_physics_2005}.
    use Galacticus_Nodes                  , only : nodeComponentSatellite
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Numerical_Constants_Math
    use Error_Functions
    use Galactic_Structure_Densities
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Vectors
    use Dark_Matter_Halo_Scales
    use Root_Finder
    implicit none
    class           (satelliteTidalStrippingZentner2005), intent(inout)         :: self
    type            (treeNode                          ), intent(inout), target :: node
    type            (treeNode                          ), pointer               :: nodeHost
    class           (nodeComponentSatellite            ), pointer               :: satellite
    double precision                                    , dimension(3)          :: position                      , velocity
    double precision                                    , parameter             :: toleranceAbsolute      =0.0d0 , toleranceRelative=1.0d-3
    double precision                                    , parameter             :: radiusZero             =0.0d0
    double precision                                    , parameter             :: radiusTidalTinyFraction=1.0d-6
    type            (rootFinder                        ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                                                            :: massSatellite                 , densityHost             , &
         &                                                                         massEnclosedHost              , frequencyAngular        , &
         &                                                                         periodOrbital                 , radius                  , &
         &                                                                         tidalTensor                   , radiusTidal             , &
         &                                                                         massOuterSatellite            , frequencyRadial

    ! Get required quantities from satellite and host nodes.
    nodeHost          => node     %mergesWith()
    satellite         => node     %satellite ()
    massSatellite     =  satellite%boundMass ()
    position          =  satellite%position  ()
    velocity          =  satellite%velocity  ()
    radius            =  Vector_Magnitude                (         position                          )
    densityHost       =  Galactic_Structure_Density      (nodeHost,position,coordinateSystemCartesian)
    massEnclosedHost  =  Galactic_Structure_Enclosed_Mass(nodeHost,radius                            )
    ! Compute the orbital period.
    frequencyAngular  = +Vector_Magnitude(Vector_Product(position,velocity)) &
         &              /radius**2                                           &
         &              *kilo                                                &
         &              *gigaYear                                            &
         &              /megaParsec
    frequencyRadial   = +abs             (   Dot_Product(position,velocity)) &
         &              /radius**2                                           &
         &              *kilo                                                &
         &              *gigaYear                                            &
         &              /megaParsec
    ! Find the orbital period. We use the larger of the angular and radial frequencies to avoid numerical problems for purely
    ! radial or purely circular orbits.
    periodOrbital     = 2.0d0*Pi/max(frequencyAngular,frequencyRadial)
    ! Find the tidal tensor.
    tidalTensor       = -gravitationalConstantGalacticus &
         &              *(                               &
         &                 2.0d0                         &
         &                *massEnclosedHost              &
         &                /radius**3                     &
         &                -4.0d0                         &
         &                *Pi                            &
         &                *densityHost                   &
         &               )                               &
         &              *(kilo*gigaYear/megaParsec)**2
    ! If the tidal force is stretching (not compressing), compute the tidal radius.
    if     (                                                                 &
         &   frequencyAngular**2                               > tidalTensor &
         &  .and.                                                            &
         &   massSatellite                                     >  0.0d0      &
         &  .and.                                                            &
         &   Galactic_Structure_Enclosed_Mass(node,radiusZero) >= 0.0d0      &
         & ) then
       ! Check if node differs from previous one for which we performed calculations.
       if (node%uniqueID() /= self%lastUniqueID) call self%calculationReset(node)
       ! Initial estimate of the tidal radius.
       zentner2005TidalPull=  frequencyAngular**2-tidalTensor
       if (self%radiusTidalPrevious <= 0.0d0) then
          self%radiusTidalPrevious=+(                                 &
               &                     +gravitationalConstantGalacticus &
               &                     *massSatellite                   &
               &                     /zentner2005TidalPull            &
               &                     *(kilo*gigaYear/megaParsec)**2   &
               &                    )**(1.0d0/3.0d0)
          self%expandMultiplier   =+2.0d0
       end if
       ! Check if tidal radius will lie outside of current boundary.
       if (Galactic_Structure_Enclosed_Mass(node,self%radiusTidalPrevious) >= massSatellite) then
          ! Tidal radius lies outside current boundary, so no additional stripping will occur.
          zentner2005MassLossRate=0.0d0
       else
          ! Find the tidal radius in the dark matter profile.
          if (.not.finder%isInitialized()) then
             call finder%rootFunction(zentner2005TidalRadiusSolver)
             call finder%tolerance   (toleranceAbsolute,toleranceRelative)
          end if
        call finder%rangeExpand(                                                             &
               &                rangeExpandUpward            =1.0d0*self%expandMultiplier  , &
               &                rangeExpandDownward          =1.0d0/self%expandMultiplier  , &
               &                rangeExpandDownwardSignExpect=rangeExpandSignExpectNegative, &
               &                rangeExpandUpwardSignExpect  =rangeExpandSignExpectPositive, &
               &                rangeExpandType              =rangeExpandMultiplicative      &
               &               )
          zentner2005Node => node
          ! Check for extremes.
          if (zentner2005TidalRadiusSolver(radiusTidalTinyFraction*self%radiusTidalPrevious) >  0.0d0) then
             ! Complete stripping.
             radiusTidal       =0.0d0
             massOuterSatellite=massSatellite
          else
             ! Find the tidal radius, using the previous result as an initial guess.
             radiusTidal       =finder%find(rootGuess=self%radiusTidalPrevious)
             massOuterSatellite=max(                                                                  &
                  &                 massSatellite-Galactic_Structure_Enclosed_Mass(node,radiusTidal), &
                  &                 0.0d0                                                             &
                  &                )
             self%radiusTidalPrevious  =radiusTidal             
             self%expandMultiplier=1.2d0
          end if
          zentner2005MassLossRate=-self%efficiency    &
               &                  *massOuterSatellite &
               &                  /periodOrbital
       end if
    else
       zentner2005MassLossRate=0.0d0
    end if      
    return
  end function zentner2005MassLossRate
  
  double precision function zentner2005TidalRadiusSolver(radius)
    !% Root function used to find the tidal radius within a subhalo.
    use Numerical_Constants_Prefixes
    use Numerical_Constants_Astronomical
    use Numerical_Constants_Physical
    use Galactic_Structure_Enclosed_Masses
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: enclosedMass

    ! Get the satellite component.
    enclosedMass                =+Galactic_Structure_Enclosed_Mass(zentner2005Node,radius)
    zentner2005TidalRadiusSolver=+zentner2005TidalPull               &
         &                       -gravitationalConstantGalacticus    &
         &                       *enclosedMass                       &
         &                       /radius                         **3 &
         &                       *(                                  &
         &                         +kilo                             &
         &                         *gigaYear                         &
         &                         /megaParsec                       &
         &                        )                              **2
    return
  end function zentner2005TidalRadiusSolver

  subroutine zentner2005CalculationReset(self,node)
    !% Reset the stored tidal radii.
    implicit none
    class(satelliteTidalStrippingZentner2005), intent(inout) :: self
    type (treeNode                          ), intent(inout) :: node

    self%radiusTidalPrevious=-1.0d0
    self%lastUniqueID       =node%uniqueID()
    return
  end subroutine zentner2005CalculationReset
