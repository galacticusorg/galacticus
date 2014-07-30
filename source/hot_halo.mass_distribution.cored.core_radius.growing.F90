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

!% An implementation of the hot halo mass distribution core radius class in which the core grows as the hot halo content is depleted.
    
  use Tables

  integer         , parameter :: growingCoreRadiusTablePointsPerDecade=100
  double precision            :: growingCoreRadiusOverScaleRadius             , growingCoreRadiusOverVirialRadiusMaximum
  logical                     :: growingInitialized                   =.false.

  !# <hotHaloMassDistributionCoreRadius name="hotHaloMassDistributionCoreRadiusGrowing">
  !#  <description>Provides an implementation of the hot halo mass distribution core radius class in which the core grows as the hot halo content is depleted.</description>
  !# </hotHaloMassDistributionCoreRadius>
  type, extends(hotHaloMassDistributionCoreRadiusClass) :: hotHaloMassDistributionCoreRadiusGrowing
     !% An implementation of the hot halo mass distribution core radius class in which the core grows as the hot halo content is depleted.
     private
     double precision                                        :: coreRadiusOverScaleRadius      , coreRadiusOverVirialRadiusMaximum
     double precision                                        :: coreRadiusMaximum              , coreRadiusMinimum
     double precision                                        :: hotGasFractionSaved            , coreRadiusOverVirialRadiusInitialSaved, &
          &                                                     coreRadiusOverVirialRadiusSaved
     integer                                                 :: coreRadiusTableCount
     logical                                                 :: coreRadiusTableInitialized
     type            (table1DLogarithmicLinear)              :: coreRadiusTable
     class           (table1D                 ), allocatable :: coreRadiusTableInverse
   contains
     final     ::                 growingDestructor
     procedure :: radius       => growingRadius
     procedure :: stateStore   => growingStateStore
     procedure :: stateRestore => growingStateRestore
  end type hotHaloMassDistributionCoreRadiusGrowing

  interface hotHaloMassDistributionCoreRadiusGrowing
     !% Constructors for the {\tt growing} hot halo mass distribution core radius class.
     module procedure growingDefaultConstructor
     module procedure growingConstructor
  end interface hotHaloMassDistributionCoreRadiusGrowing

contains

  function growingDefaultConstructor()
    !% Default constructor for the {\tt growing} hot halo mass distribution core radius class.
    use Input_Parameters
    implicit none
    type(hotHaloMassDistributionCoreRadiusGrowing) :: growingDefaultConstructor

    if (.not.growingInitialized) then
       !$omp critical (hotHaloMassDistributionCoreRadiusGrowingInitialize)
       if (.not.growingInitialized) then
          !@ <inputParameter>
          !@   <name>hotHaloCoreRadiusOverScaleRadius</name>
          !@   <defaultValue>0.1</defaultValue>
          !@   <description>
          !@     The core radius in the hot halo density profile in units of the dark matter profile scale radius.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloCoreRadiusOverScaleRadius',growingCoreRadiusOverScaleRadius,defaultValue=0.1d0)
          !@ <inputParameter>
          !@   <name>hotHaloCoreRadiusOverVirialRadiusMaximum</name>
          !@   <defaultValue>10</defaultValue>
          !@   <description>
          !@     The maximum core radius in the ``cored isothermal'' hot halo density profile in units of the virial radius.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('hotHaloCoreRadiusOverVirialRadiusMaximum',growingCoreRadiusOverVirialRadiusMaximum,defaultValue=10.0d0)
          ! Record that this implementation is now initialized.
          growingInitialized=.true.
       end if
       !$omp end critical (hotHaloMassDistributionCoreRadiusGrowingInitialize)
    end if
    growingDefaultConstructor=growingConstructor(growingCoreRadiusOverScaleRadius,growingCoreRadiusOverVirialRadiusMaximum)
    return
  end function growingDefaultConstructor

  function growingConstructor(coreRadiusOverScaleRadius,coreRadiusOverVirialRadiusMaximum)
    !% Default constructor for the {\tt growing} hot halo mass distribution core radius class.
    use Input_Parameters
    use Galacticus_Error
    implicit none
    type            (hotHaloMassDistributionCoreRadiusGrowing)                :: growingConstructor
    double precision                                          , intent(in   ) :: coreRadiusOverScaleRadius,coreRadiusOverVirialRadiusMaximum

    growingConstructor%coreRadiusOverScaleRadius        =coreRadiusOverScaleRadius
    growingConstructor%coreRadiusOverVirialRadiusMaximum=coreRadiusOverVirialRadiusMaximum
    growingConstructor%coreRadiusTableInitialized       =.false.
    ! Ensure that the dark matter profile supports the scale property.
    if (.not.defaultDarkMatterProfileComponent%scaleIsGettable())                                                          &
         & call Galacticus_Error_Report                                                                                    &
         &      (                                                                                                          &
         &       'growingConstructor'                                                                                    , &
         &       'method requires a dark matter profile component that provides a gettable "scale" property.'//            &
         &       Galacticus_Component_List(                                                                                &
         &                                 'darkMatterProfile'                                                           , &
         &                                  defaultDarkMatterProfileComponent%scaleAttributeMatch(requireGettable=.true.)  &
         &                                )                                                                                &
         &      )    
   return
  end function growingConstructor

  subroutine growingDestructor(self)
    !% Destructor for the {\tt growing} hot halo mass distribution class.
    implicit none
    type(hotHaloMassDistributionCoreRadiusGrowing), intent(inout) :: self

    call self%coreRadiusTable%destroy()
    return
  end subroutine growingDestructor

  double precision function growingRadius(self,node)
    !% Return the core radius of the hot halo mass distribution.
    use Cosmology_Parameters
    use Dark_Matter_Halo_Scales
    implicit none
    class           (hotHaloMassDistributionCoreRadiusGrowing), intent(inout)          :: self
    type            (treeNode                                ), intent(inout), pointer :: node
    class           (nodeComponentBasic                      )               , pointer :: basicComponent
    class           (nodeComponentHotHalo                    )               , pointer :: hotHaloComponent
    class           (nodeComponentDarkMatterProfile          )               , pointer :: darkMatterProfileComponent
    class           (cosmologyParametersClass                )               , pointer :: cosmologyParameters_
    class           (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                                                   :: hotGasFraction                   , coreRadiusOverVirialRadius, &
         &                                                                                coreRadiusOverVirialRadiusInitial, targetValue
    logical                                                                            :: makeTable

    ! Get components.
    basicComponent             => node%basic            ()
    hotHaloComponent           => node%hotHalo          ()
    darkMatterProfileComponent => node%darkMatterProfile()
    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Find the fraction of gas in the hot halo relative to that expected from the universal baryon fraction.
    hotGasFraction=(hotHaloComponent%mass()/basicComponent%mass())*(cosmologyParameters_%OmegaMatter()/cosmologyParameters_%OmegaBaryon())
    ! Return an arbitrary value for empty halos.
    if (hotGasFraction <= 0.0d0) then
       growingRadius=self%coreRadiusOverScaleRadius
       return
    end if
    ! Comptue the desired core radius (in units of the virial radius) for a fully populated halo.
    coreRadiusOverVirialRadiusInitial=          &
         &  self%coreRadiusOverScaleRadius      &
         & *darkMatterProfileComponent%scale()  &
         & /darkMatterHaloScale_%virialRadius(node)
    ! Check if the initial core radius and hot gas fraction equal the previously stored values.
    if     (                                                                                        &
         &  .not.                                                                                   &
         &       (                                                                                  &
         &         coreRadiusOverVirialRadiusInitial == self%coreRadiusOverVirialRadiusInitialSaved &
         &        .and.                                                                             &
         &         hotGasFraction                    == self%hotGasFractionSaved                    &
         &       )                                                                                  &
         & ) then
       ! Create a tabulation of core radius vs. virial density factor if necessary.
       makeTable=.not.self%coreRadiusTableInitialized
       if (.not.makeTable) makeTable=(coreRadiusOverVirialRadiusInitial < self%coreRadiusTable%x(-1))
       if (makeTable) then
          self%coreRadiusMinimum   =min(self%coreRadiusOverScaleRadius,coreRadiusOverVirialRadiusInitial)
          self%coreRadiusMaximum   =self%coreRadiusOverVirialRadiusMaximum
          self%coreRadiusTableCount=int(log10(self%coreRadiusMaximum/self%coreRadiusMinimum)*dble(growingCoreRadiusTablePointsPerDecade))+1
          call self%coreRadiusTable%destroy (                                                                       )
          call self%coreRadiusTable%create  (self%coreRadiusMaximum,self%coreRadiusMinimum,self%coreRadiusTableCount)
          call self%coreRadiusTable%populate(growingCoreVirialDensityFunction(self%coreRadiusTable%xs())            )
          call self%coreRadiusTable%reverse (self%coreRadiusTableInverse                                            )
          self%coreRadiusTableInitialized=.true.
       end if
       ! Compute the target value of the function giving the density at the virial radius per unit gas mass.
       targetValue=growingCoreVirialDensityFunction(coreRadiusOverVirialRadiusInitial)*hotGasFraction
       ! Interpolate to get the required core radius.
       if      (hotGasFraction >= 1.0d0                    ) then
          coreRadiusOverVirialRadius=coreRadiusOverVirialRadiusInitial
       else if (targetValue    <= self%coreRadiusTable%y(1)) then
          coreRadiusOverVirialRadius=self%coreRadiusTable%x(1)
       else
          coreRadiusOverVirialRadius=self%coreRadiusTableInverse%interpolate(targetValue)
       end if
       self%coreRadiusOverVirialRadiusInitialSaved=coreRadiusOverVirialRadiusInitial
       self%hotGasFractionSaved                   =hotGasFraction
       self%coreRadiusOverVirialRadiusSaved       =coreRadiusOverVirialRadius
    end if
    ! Compute the resulting core radius.
    growingRadius=self%coreRadiusOverVirialRadiusSaved*darkMatterHaloScale_%virialRadius(node)
    return
  end function growingRadius

  elemental double precision function growingCoreVirialDensityFunction(radiusOverVirialRadius)
    !% Returns the function $(1+r_{\rm c}^2)[1-r_{\rm c} \tan^{-1}(1/r_{\rm c}]$ which is proportional to the density at the
    !% virial radius of a cored isothermal profile with core radius $r_{\rm c}$ (in units of the virial radius) per unit mass.
    implicit none
    double precision, intent(in   ) :: radiusOverVirialRadius

    growingCoreVirialDensityFunction=(1.0d0+radiusOverVirialRadius**2)*(1.0d0-radiusOverVirialRadius*atan(1.0d0/radiusOverVirialRadius))
    return
  end function growingCoreVirialDensityFunction

  subroutine growingStateStore(self,stateFile,fgslStateFile)
    !% Write the tablulation state to file.
    use FGSL
    implicit none
    class  (hotHaloMassDistributionCoreRadiusGrowing), intent(inout) :: self
    integer                                          , intent(in   ) :: stateFile
    type   (fgsl_file                               ), intent(in   ) :: fgslStateFile

    write (stateFile) self%coreRadiusMinimum,self%coreRadiusMaximum
    return
  end subroutine growingStateStore

  subroutine growingStateRestore(self,stateFile,fgslStateFile)
    !% Retrieve the tabulation state from the file.
    use FGSL
    implicit none
    class  (hotHaloMassDistributionCoreRadiusGrowing), intent(inout) :: self
    integer                                          , intent(in   ) :: stateFile
    type   (fgsl_file                               ), intent(in   ) :: fgslStateFile

    ! Read the minimum and maximum tabulated times.
    read (stateFile) self%coreRadiusMinimum,self%coreRadiusMaximum
    ! Force retabulation on next evaluation.
    self%coreRadiusTableInitialized=.false.
    return
  end subroutine growingStateRestore
