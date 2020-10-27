!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

  !% Implements a class for ram pressure stripping of hot halos based on the methods of \cite{font_colours_2008}.

  use :: Dark_Matter_Halo_Scales     , only : darkMatterHaloScaleClass
  use :: Hot_Halo_Mass_Distributions , only : hotHaloMassDistributionClass
  use :: Hot_Halo_Ram_Pressure_Forces, only : hotHaloRamPressureForce     , hotHaloRamPressureForceClass
  use :: Kind_Numbers                , only : kind_int8

  !# <hotHaloRamPressureStripping name="hotHaloRamPressureStrippingFont2008">
  !#  <description>
  !#   A hot halo ram pressure stripping class based on the methods of \cite{font_colours_2008}. Specifically, the radius,
  !#   $r_\mathrm{rp}$, is computed as the solution of
  !#   \begin{equation}
  !#   \alpha_\mathrm{rp} {\mathrm{G} M_\mathrm{satellite}(r_\mathrm{rp}) \rho_\mathrm{hot, satellite}(r_\mathrm{rp}) \over
  !#   r_\mathrm{rp} } = \mathcal{F}_\mathrm{ram, hot, host},
  !#   \end{equation}
  !#   where $M_\mathrm{satellite}(r)$ is the total mass of the satellite within radius $r$, $\mathcal{F}_\mathrm{ram, hot, host}$
  !#   is the ram pressure force due to the hot halo (computed using the selected hot halo ram pressure force method; see
  !#   \S\ref{sec:HotHaloRamPressureForce}). The parameter $\alpha_\mathrm{rp}=${\normalfont \ttfamily
  !#   [formFactor]} is a geometric factor of order unity.
  !#  </description>
  !# </hotHaloRamPressureStripping>
  type, extends(hotHaloRamPressureStrippingClass) :: hotHaloRamPressureStrippingFont2008
     !% Implementation of a hot halo ram pressure stripping class based on the methods of \cite{font_colours_2008}.
     private
     class           (darkMatterHaloScaleClass    ), pointer :: darkMatterHaloScale_     => null()
     class           (hotHaloRamPressureForceClass), pointer :: hotHaloRamPressureForce_ => null()
     class           (hotHaloMassDistributionClass), pointer :: hotHaloMassDistribution_ => null()
     double precision                                        :: formFactor
     integer         (kind_int8                   )          :: uniqueIDLast             =  -1
     double precision                                        :: radiusLast               =  -1.0d0
   contains
     final     ::                   font2008Destructor
     procedure :: radiusStripped => font2008RadiusStripped
  end type hotHaloRamPressureStrippingFont2008

  interface hotHaloRamPressureStrippingFont2008
     !% Constructors for the {\normalfont \ttfamily font2008} hot halo ram pressure timescale class.
     module procedure font2008ConstructorParameters
     module procedure font2008ConstructorInternal
  end interface hotHaloRamPressureStrippingFont2008

  ! Global variables used in root finding.
  class           (hotHaloRamPressureStrippingFont2008), pointer :: font2008Self
  type            (treeNode                           ), pointer :: font2008Node
  double precision                                               :: font2008ForceRamPressure
  !$omp threadprivate(font2008Self,font2008Node,font2008ForceRamPressure)

contains

  function font2008ConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily font2008} hot halo ram pressure stripping class which builds the object from a parameter set.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (hotHaloRamPressureStrippingFont2008)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (darkMatterHaloScaleClass           ), pointer       :: darkMatterHaloScale_
    class           (hotHaloRamPressureForceClass       ), pointer       :: hotHaloRamPressureForce_
    class           (hotHaloMassDistributionClass       ), pointer       :: hotHaloMassDistribution_
    double precision                                                     :: formFactor

    !# <inputParameter>
    !#   <name>formFactor</name>
    !#   <defaultValue>2.0d0</defaultValue>
    !#   <description>The form factor appearing in the gravitational binding force (per unit area) in the ram pressure stripping model
    !#      of \citeauthor{font_colours_2008}~(\citeyear{font_colours_2008}; their eqn.~4).</description>
    !#   <source>parameters</source>
    !# </inputParameter>
    !# <objectBuilder class="darkMatterHaloScale"     name="darkMatterHaloScale_"     source="parameters"/>
    !# <objectBuilder class="hotHaloRamPressureForce" name="hotHaloRamPressureForce_" source="parameters"/>
    !# <objectBuilder class="hotHaloMassDistribution" name="hotHaloMassDistribution_" source="parameters"/>
    self=hotHaloRamPressureStrippingFont2008(formFactor,darkMatterHaloScale_,hotHaloRamPressureForce_,hotHaloMassDistribution_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="darkMatterHaloScale_"    />
    !# <objectDestructor name="hotHaloRamPressureForce_"/>
    !# <objectDestructor name="hotHaloMassDistribution_"/>
    return
  end function font2008ConstructorParameters

  function font2008ConstructorInternal(formFactor,darkMatterHaloScale_,hotHaloRamPressureForce_,hotHaloMassDistribution_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily font2008} hot halo ram pressure stripping class.
    implicit none
    type            (hotHaloRamPressureStrippingFont2008)                        :: self
    class           (darkMatterHaloScaleClass           ), intent(in   ), target :: darkMatterHaloScale_
    class           (hotHaloRamPressureForceClass       ), intent(in   ), target :: hotHaloRamPressureForce_
    class           (hotHaloMassDistributionClass       ), intent(in   ), target :: hotHaloMassDistribution_
    double precision                                     , intent(in   )         :: formFactor
    !# <constructorAssign variables="formFactor, *darkMatterHaloScale_, *hotHaloRamPressureForce_, *hotHaloMassDistribution_"/>

    return
  end function font2008ConstructorInternal

  subroutine font2008Destructor(self)
    !% Destructor for the {\normalfont \ttfamily font2008} hot halo ram pressure stripping class.
    implicit none
    type(hotHaloRamPressureStrippingFont2008), intent(inout) :: self

    !# <objectDestructor name="self%darkMatterHaloScale_"    />
    !# <objectDestructor name="self%hotHaloRamPressureForce_"/>
    !# <objectDestructor name="self%hotHaloMassDistribution_"/>
    return
  end subroutine font2008Destructor

  double precision function font2008RadiusStripped(self,node)
    !% Return the ram pressure stripping radius due to the hot halo using the model of \cite{font_colours_2008}.
    use :: Galacticus_Display, only : Galacticus_Display_Message, verbositySilent
    use :: Galacticus_Error  , only : Galacticus_Error_Report   , errorStatusSuccess
    use :: Root_Finder       , only : rangeExpandMultiplicative , rangeExpandSignExpectNegative, rangeExpandSignExpectPositive, rootFinder
    implicit none
    class           (hotHaloRamPressureStrippingFont2008), intent(inout), target :: self
    type            (treeNode                           ), intent(inout), target :: node
    double precision                                     , parameter             :: toleranceAbsolute             =0.0d+0, toleranceRelative=1.0d-3
    double precision                                     , parameter             :: radiusSmallestOverRadiusVirial=1.0d-6
    type            (rootFinder                         ), save                  :: finder
    !$omp threadprivate(finder)
    double precision                                                             :: radiusVirial                         , radiusVirialRoot        , &
         &                                                                          radiusSmallRoot
    integer                                                                      :: status
    type            (varying_string                     )                        :: message
    character       (len=16                             )                        :: label

    ! Get the virial radius of the satellite.
    radiusVirial=self%darkMatterHaloScale_%virialRadius(node)
    ! Test whether node is a satellite.
    if (node%isSatellite()) then
       ! Set a pointer to the satellite node.
       font2008Self             => self
       font2008Node             =>                                     node
       ! Get the ram pressure force due to the hot halo.
       font2008ForceRamPressure =  self%hotHaloRamPressureForce_%force(node)
       ! Find the radial range within which the ram pressure radius must lie.
       radiusVirialRoot=font2008RadiusSolver(radiusVirial)
       if      (radiusVirialRoot >= 0.0d0) then
          ! The ram pressure force is not sufficiently strong to strip even at the satellite virial radius - simply return the
          ! virial radius as the stripping radius in this case.
          font2008RadiusStripped=radiusVirial
       else
          radiusSmallRoot=font2008RadiusSolver(radiusSmallestOverRadiusVirial*radiusVirial)
          if (radiusSmallRoot <= 0.0d0) then
             ! The ram pressure force can strip to (essentially) arbitrarily small radii.
             font2008RadiusStripped=0.0d0
          else
             ! Solver for the ram pressure stripping radius.
             if (.not.finder%isInitialized()) then
                call finder%rootFunction(font2008RadiusSolver)
                call finder%tolerance   (toleranceAbsolute,toleranceRelative)
             end if
             ! If we have a previously found radius, and if the node is the same as the previous node for which this function was
             ! called, then use that previous radius as a guess for the new solution. Note that we do not reset the previous
             ! radius between ODE steps (i.e. we do not make use of "calculationReset" events to reset the radius) as we
             ! specifically want to retain knowledge from the previous step.
             if (self%radiusLast > 0.0d0 .and. node%uniqueID() == self%uniqueIDLast) then
                call finder%rangeExpand(                                                                           &
                     &                  rangeExpandDownward          =0.9d0                                      , &
                     &                  rangeExpandUpward            =1.1d0                                      , &
                     &                  rangeExpandType              =rangeExpandMultiplicative                  , &
                     &                  rangeDownwardLimit           =radiusSmallestOverRadiusVirial*radiusVirial, &
                     &                  rangeUpwardLimit             =                               radiusVirial, &
                     &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive              , &
                     &                  rangeExpandUpwardSignExpect  =rangeExpandSignExpectNegative                &
                     &                 )
             else
                call finder%rangeExpand(                                                                           &
                     &                  rangeExpandDownward          =0.5d0                                      , &
                     &                  rangeExpandType              =rangeExpandMultiplicative                  , &
                     &                  rangeDownwardLimit           =radiusSmallestOverRadiusVirial*radiusVirial, &
                     &                  rangeExpandDownwardSignExpect=rangeExpandSignExpectPositive                &
                     &                 )
                self%radiusLast  =radiusVirial
                self%uniqueIDLast=node%uniqueID()
             end if
             font2008RadiusStripped=finder%find(rootGuess=min(self%radiusLast,radiusVirial),status=status)
             if (status /= errorStatusSuccess) then
                message='virial radius / root function at virial radius = '
                write (label,'(e12.6)') radiusVirial
                message=message//trim(adjustl(label))
                write (label,'(e12.6)') radiusVirialRoot
                message=message//" / "//trim(adjustl(label))
                write (label,'(e12.6)') font2008RadiusSolver(radiusVirial)
                message=message//" : "//trim(adjustl(label))
                call Galacticus_Display_Message(message,verbositySilent)
                message='small radius / root function at small radius = '
                write (label,'(e12.6)') radiusSmallestOverRadiusVirial*radiusVirial
                message=message//trim(adjustl(label))
                write (label,'(e12.6)') radiusSmallRoot
                message=message//" / "//trim(adjustl(label))
                write (label,'(e12.6)') font2008RadiusSolver(radiusSmallestOverRadiusVirial*radiusVirial)
                message=message//" : "//trim(adjustl(label))
                call Galacticus_Display_Message(message,verbositySilent)
                call Galacticus_Error_Report('root finding failed'//{introspection:location})
             end if
             self%radiusLast=font2008RadiusStripped
          end if
       end if
    else
       ! If node is not a satellite, return a ram pressure stripping radius equal to the virial radius.
       font2008RadiusStripped=radiusVirial
    end if
    return
  end function font2008RadiusStripped

  double precision function font2008RadiusSolver(radius)
    !% Root function used in finding the ram pressure stripping radius.
    use :: Galactic_Structure_Enclosed_Masses, only : Galactic_Structure_Enclosed_Mass
    use :: Galactic_Structure_Options        , only : componentTypeAll                , massTypeAll
    use :: Numerical_Constants_Astronomical      , only : gravitationalConstantGalacticus
    implicit none
    double precision, intent(in   ) :: radius
    double precision                :: massEnclosed  , forceBindingGravitational, &
         &                             densityHotHalo

    ! Get the hot halo mass distribution.
    massEnclosed             =+Galactic_Structure_Enclosed_Mass(font2008Node,radius,massType=massTypeAll,componentType=componentTypeAll)
    densityHotHalo           =+font2008Self%hotHaloMassDistribution_%density(font2008Node,radius)
    forceBindingGravitational=+font2008Self%formFactor         &
         &                    *gravitationalConstantGalacticus &
         &                    *massEnclosed                    &
         &                    *densityHotHalo                  &
         &                    /radius
    if (forceBindingGravitational >= 0.0d0) then
       font2008RadiusSolver=forceBindingGravitational-font2008ForceRamPressure
    else
       font2008RadiusSolver=                         -font2008ForceRamPressure
    end if
    return
  end function font2008RadiusSolver
