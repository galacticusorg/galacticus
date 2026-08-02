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
  An implementation of dark matter halo virial density contrasts based on the percolation analysis of :cite:t:`more_overdensity_2011`.
  !!}

  use :: Cosmology_Functions, only : cosmologyFunctionsClass
  use :: Kind_Numbers       , only : kind_int8
  use :: Tables             , only : table2DLogLogLin
  use :: Resource_Manager   , only : resourceManager

  !![
  <virialDensityContrast name="virialDensityContrastPercolation" recursive="yes" docformat="rst">
   <description>
   A dark matter halo virial density contrast class based on the percolation analysis of :cite:t:`more_overdensity_2011`. The virial density contrast is computed to be consistent with a given friends-of-friends algorithm linking length using the percolation-theory-motivated calibration of :cite:t:`more_overdensity_2011`. Specifically, the friends-of-friends algorithm is assumed to link together particles forming an isodensity surface of density :math:`\rho = \bar{\rho} n_\mathrm{c}/b^3`, where :math:`\bar{\rho}` is the mean density of the universe, :math:`n_\mathrm{c}=0.652960` is a critical density for percolation as given by :cite:t:`more_overdensity_2011`, and :math:`b` is the linking length. Given this bounding density, the virial density contrast is found by requiring that the halo contain the required mass within such a bounding density, given the halo density profile.
   </description>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastPercolation
     !!{RST
     A dark matter halo virial density contrast class based on the percolation analysis of :cite:t:`more_overdensity_2011`.
     !!}
     private
     double precision                                   :: linkingLength
     type            (varying_string         )          :: fileName
     class           (cosmologyFunctionsClass), pointer :: cosmologyFunctions_             => null()
     class           (*                      ), pointer :: percolationObjects_             => null()
     type            (resourceManager        )          :: percolationObjectsManager
     ! Tabulation of density contrast vs. time and mass.
     double precision                                   :: densityContrastTableTimeMinimum           , densityContrastTableTimeMaximum
     double precision                                   :: densityContrastTableMassMinimum           , densityContrastTableMassMaximum
     integer                                            :: densityContrastTableMassCount             , densityContrastTableTimeCount
     logical                                            :: densityContrastTableInitialized =  .false.
     type            (table2DLogLogLin       )          :: densityContrastTable
     integer                                            :: densityContrastTableRemakeCount
   contains
     !![
     <methods docformat="rst">
       <method description="Tabulate the virial density contrast as a function of mass and time." method="tabulate"    />
     </methods>
     !!]
     final     ::                                percolationDestructor
     procedure :: densityContrast             => percolationDensityContrast
     procedure :: densityContrastRateOfChange => percolationDensityContrastRateOfChange
     procedure :: isMassDependent             => percolationIsMassDepdendent
     procedure :: tabulate                    => percolationTabulate
     procedure :: deepCopy                    => percolationDeepCopy
     procedure :: deepCopyReset               => percolationDeepCopyReset
     procedure :: deepCopyFinalize            => percolationDeepCopyFinalize
  end type virialDensityContrastPercolation

  interface virialDensityContrastPercolation
     !!{RST
     Constructors for the :galacticus-class:`virialDensityContrastPercolation` dark matter halo virial density contrast class.
     !!}
     module procedure percolationConstructorParameters
     module procedure percolationConstructorInternal
  end interface virialDensityContrastPercolation

  ! Granularity parameters for tabulations.
  integer                    , parameter :: densityContrastTableTimePointsPerDecade=5
  integer                    , parameter :: densityContrastTableMassPointsPerDecade=5

  ! Intervals, measured in lattice steps, to which the bounds of each tabulated axis are pinned. The time axis is anchored to
  ! half decades since a whole decade of cosmic time spans most of the history of the universe; the mass axis, which already
  ! spans eleven decades by default, is anchored to whole decades.
  ! Ranges tabulated irrespective of the mass and time requested.
  double precision           , parameter :: densityContrastTableTimeMinimumDefault = 1.0d-03, densityContrastTableTimeMaximumDefault=2.0d+01
  double precision           , parameter :: densityContrastTableMassMinimumDefault = 4.0d+05, densityContrastTableMassMaximumDefault=1.0d+16

  integer                    , parameter :: densityContrastTableTimeAnchorEvery    =densityContrastTableTimePointsPerDecade/2
  integer                    , parameter :: densityContrastTableMassAnchorEvery    =densityContrastTableMassPointsPerDecade

  ! Module-scope record of state used when solving for the percolation density contrast.
  logical                                :: solving                                =.false.
  double precision                       :: densityContrastCurrent                 =-1.0d0
  !$omp threadprivate(densityContrastCurrent,solving)

contains

  recursive function percolationConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`virialDensityContrastPercolation` dark matter halo virial density contrast class that takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctions                                      , cosmologyFunctionsClass
    use :: Functions_Global   , only : Virial_Density_Contrast_Percolation_Objects_Constructor_
    use :: Input_Parameters   , only : inputParameter                                          , inputParameters
    implicit none
    type            (virialDensityContrastPercolation)                :: self
    type            (inputParameters                 ), intent(inout) :: parameters
    class           (cosmologyFunctionsClass         ), pointer       :: cosmologyFunctions_
    class           (*                               ), pointer       :: percolationObjects_
    double precision                                                  :: linkingLength

    !![
    <inputParameter docformat="rst">
     <name>linkingLength</name>
     <source>parameters</source>
     <defaultValue>0.2d0</defaultValue>
     <description>
     The friends-of-friends linking length to use in computing virial density contrasts with the percolation analysis of :cite:t:`more_overdensity_2011`.
     </description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    ! Build a pointer to a container object which stores all objects needed by the percolation density solver.
    percolationObjects_ => Virial_Density_Contrast_Percolation_Objects_Constructor_(parameters)
    self=virialDensityContrastPercolation(linkingLength,cosmologyFunctions_,percolationObjects_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    return
  end function percolationConstructorParameters

  recursive function percolationConstructorInternal(linkingLength,cosmologyFunctions_,percolationObjects_) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`virialDensityContrastPercolation` dark matter halo virial density contrast class.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : operator(//) , char
    use :: Input_Parameters  , only : inputParameter, inputParameters
    use :: Table_Caches      , only : Table_Cache_File_Name
    implicit none
    type            (virialDensityContrastPercolation)                        :: self
    double precision                                  , intent(in   )         :: linkingLength
    class           (cosmologyFunctionsClass         ), intent(in   ), target :: cosmologyFunctions_
    class           (*                               ), intent(in   ), target :: percolationObjects_
    !![
    <constructorAssign variables="linkingLength, *cosmologyFunctions_, *percolationObjects_"/>
    !!]

    ! Add management to the shared percolationObjects resource.
    self%percolationObjectsManager=resourceManager(percolationObjects_)
    ! File name for tabulation.
    self%fileName=Table_Cache_File_Name(                                                                                          &
         &                              subDirectory    ='darkMatterHalos'                                                      , &
         &                              objectType      =char(self%objectType      (                                          )), &
         &                              hashedDescriptor=char(self%hashedDescriptor(includeSourceDigest          =.true.        , &
         &                                                                          includeFileModificationTimes=.true.        ))  &
         &                             )
    ! Initialize tabulations.
    self%densityContrastTableTimeMinimum=densityContrastTableTimeMinimumDefault
    self%densityContrastTableTimeMaximum=densityContrastTableTimeMaximumDefault
    self%densityContrastTableMassMinimum=densityContrastTableMassMinimumDefault
    self%densityContrastTableMassMaximum=densityContrastTableMassMaximumDefault
    self%densityContrastTableInitialized=.false.
    self%densityContrastTableRemakeCount= 0
    return
  end function percolationConstructorInternal

  subroutine percolationDestructor(self)
    !!{RST
    Destructor for the :galacticus-class:`virialDensityContrastPercolation` dark matter halo virial density contrast class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type(virialDensityContrastPercolation), intent(inout) :: self

    if (self%densityContrastTableInitialized) call self%densityContrastTable%destroy()
    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    return
  end subroutine percolationDestructor

  subroutine percolationTabulate(self,mass,time)
    !!{RST
    Tabulate virial density contrast as a function of mass and time for the ``percolation`` density contrast class.
    !!}
    use :: Display         , only : displayCounter                             , displayCounterClear, displayIndent      , displayUnindent, &
          &                         verbosityLevelWorking
    use :: Functions_Global, only : Virial_Density_Contrast_Percolation_Solver_
    use :: Error           , only : Error_Report
    use :: Numerical_Ranges, only : Range_Pinned                               , rangeLattice       , gridSchemePerDecade
    use :: Table_Caches    , only : Table_Cache_Restore                        , Table_Cache_Store
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)               :: self
    double precision                                  , intent(in   )               :: mass       , time
    type            (rangeLattice                    )                              :: latticeMass, latticeTime
    logical                                           , allocatable, dimension(:,:) :: isComputed
    integer                                                                         :: iMass      , iTime      , &
         &                                                                             iCount     , status
    double precision                                                                :: tableMass  , tableTime

    ! Find the ranges of mass and time to tabulate, pinning each axis independently to an absolute lattice so that the
    ! tabulation - and hence every value interpolated from it - is independent of the mass and time at which it was first
    ! requested, and so that the table can be extended without changing or recomputing any previously computed value.
    latticeMass=Range_Pinned(                                                                              &
         &                                  mass                                                         , &
         &                                  densityContrastTableMassPointsPerDecade                      , &
         &                                  gridSchemePerDecade                                          , &
         &                   rangeCurrent  =[densityContrastTableMassMinimumDefault,densityContrastTableMassMaximumDefault], &
         &                   latticeCurrent=self%densityContrastTable%latticeX                            , &
         &                   anchorEvery   =densityContrastTableMassAnchorEvery                             &
         &                  )
    latticeTime=Range_Pinned(                                                                              &
         &                                  time                                                         , &
         &                                  densityContrastTableTimePointsPerDecade                      , &
         &                                  gridSchemePerDecade                                          , &
         &                   rangeCurrent  =[densityContrastTableTimeMinimumDefault,densityContrastTableTimeMaximumDefault], &
         &                   latticeCurrent=self%densityContrastTable%latticeY                            , &
         &                   anchorEvery   =densityContrastTableTimeAnchorEvery                             &
         &                  )
    if     (                                                          &
         &   self%densityContrastTable%latticeX%covers(latticeMass)   &
         &  .and.                                                     &
         &   self%densityContrastTable%latticeY%covers(latticeTime)   &
         & ) return
    ! Merge in any tabulation already cached on disk, then re-evaluate the ranges required in the light of it.
    call Table_Cache_Restore(self%densityContrastTable,self%fileName,status)
    latticeMass=Range_Pinned(                                                                              &
         &                                  mass                                                         , &
         &                                  densityContrastTableMassPointsPerDecade                      , &
         &                                  gridSchemePerDecade                                          , &
         &                   rangeCurrent  =[densityContrastTableMassMinimumDefault,densityContrastTableMassMaximumDefault], &
         &                   latticeCurrent=self%densityContrastTable%latticeX                            , &
         &                   anchorEvery   =densityContrastTableMassAnchorEvery                             &
         &                  )
    latticeTime=Range_Pinned(                                                                              &
         &                                  time                                                         , &
         &                                  densityContrastTableTimePointsPerDecade                      , &
         &                                  gridSchemePerDecade                                          , &
         &                   rangeCurrent  =[densityContrastTableTimeMinimumDefault,densityContrastTableTimeMaximumDefault], &
         &                   latticeCurrent=self%densityContrastTable%latticeY                            , &
         &                   anchorEvery   =densityContrastTableTimeAnchorEvery                             &
         &                  )
    call self%densityContrastTable%extend(latticeMass,latticeTime,isComputed)
    if (any(.not.isComputed)) then
       ! Check that we have a pointer to the required objects.
       if (.not.associated(self%percolationObjects_)) call Error_Report('no percolationObjects available'//{introspection:location})
       ! Increment the number of table remakes.
       self%densityContrastTableRemakeCount=self%densityContrastTableRemakeCount+1
       ! Record that we are in the solving phase of calculation, so we will avoid recursive calls to this function.
       solving=.true.
       ! Tabulate the density contrast at those points whose values were not preserved by the extension.
       call displayIndent('Tabulating virial density contrasts for percolation class',verbosity=verbosityLevelWorking)
       iCount=0
       do iMass=1,latticeMass%count
          tableMass=self%densityContrastTable%x(iMass)
          do iTime=1,latticeTime%count
             if (isComputed(iMass,iTime)) cycle
             tableTime=self%densityContrastTable%y(iTime)
             iCount=iCount+1
             call displayCounter(int(100.0d0*dble(iCount)/dble(count(.not.isComputed))),isNew=(iCount==1),verbosity=verbosityLevelWorking)
             call self%densityContrastTable%populate(Virial_Density_Contrast_Percolation_Solver_(tableMass,tableTime,self%linkingLength,densityContrastCurrent,self%percolationObjects_,self),iMass,iTime)
          end do
       end do
       call displayCounterClear(verbosity=verbosityLevelWorking)
       call displayUnindent('done',verbosity=verbosityLevelWorking)
       ! Merge with, and write back, the cache file.
       call Table_Cache_Store(self%densityContrastTable,self%fileName)
       ! Solving phase is finished.
       solving=.false.
    end if
    ! Record the tabulated ranges.
    self%densityContrastTableInitialized=.true.
    self%densityContrastTableMassMinimum=self%densityContrastTable%latticeX%minimum()
    self%densityContrastTableMassMaximum=self%densityContrastTable%latticeX%maximum()
    self%densityContrastTableTimeMinimum=self%densityContrastTable%latticeY%minimum()
    self%densityContrastTableTimeMaximum=self%densityContrastTable%latticeY%maximum()
    self%densityContrastTableMassCount  =self%densityContrastTable%latticeX%count
    self%densityContrastTableTimeCount  =self%densityContrastTable%latticeY%count
    return
  end subroutine percolationTabulate

  double precision function percolationDensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{RST
    Return the virial density contrast at the given epoch, based on the percolation algorithm of :cite:t:`more_overdensity_2011`.
    !!}
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)            :: self
    double precision                                  , intent(in   )            :: mass
    double precision                                  , intent(in   ) , optional :: time      , expansionFactor
    logical                                           , intent(in   ) , optional :: collapsing
    double precision                                                             :: timeActual

    ! Get the time to use.
    if (.not.solving) timeActual=self%cosmologyFunctions_%epochTime(time,expansionFactor,collapsing)
    ! Determine how to compute density contrast.
    if (solving) then
       ! Currently solving for solutions - return the current guess.
       percolationDensityContrast=densityContrastCurrent
    else
       ! Use our tabulated solution.
       call self%tabulate(mass,timeActual)
       percolationDensityContrast=self%densityContrastTable%interpolate(mass,timeActual)
    end if
    return
  end function percolationDensityContrast

  double precision function percolationDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !!{RST
    Return the rate of change of the virial density contrast at the given epoch, based on the percolation algorithm of :cite:t:`more_overdensity_2011`.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)            :: self
    double precision                                  , intent(in   )            :: mass
    double precision                                  , intent(in   ) , optional :: time          , expansionFactor
    logical                                           , intent(in   ) , optional :: collapsing
    double precision                                                             :: timeActual

    ! Get the time to use.
    timeActual=self%cosmologyFunctions_%epochTime(time,expansionFactor,collapsing)
    ! Compute the solution.
    call self%tabulate(mass,timeActual)
    percolationDensityContrastRateOfChange=self%densityContrastTable%interpolateGradient(mass,timeActual,dim=2)
    return
  end function percolationDensityContrastRateOfChange

  logical function percolationIsMassDepdendent(self)
    !!{RST
    Specify that the ``percolation`` virial density contrast class is mass-dependent.
    !!}
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self
    !$GLC attributes unused :: self

    percolationIsMassDepdendent=.true.
    return
  end function percolationIsMassDepdendent

  subroutine percolationDeepCopyReset(self)
    !!{RST
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : percolationObjectsDeepCopyReset_
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self

    self%copiedSelf => null()
    if (associated(self%cosmologyfunctions_)) call self%cosmologyfunctions_%deepCopyReset()
    if (associated(self%percolationObjects_)) call percolationObjectsDeepCopyReset_(self%percolationObjects_)
    return
  end subroutine percolationDeepCopyReset

  subroutine percolationDeepCopyFinalize(self)
    !!{RST
    Finalize a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : percolationObjectsDeepCopyFinalize_
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self

    if (associated(self%cosmologyfunctions_)) call self%cosmologyfunctions_%deepCopyFinalize()
    if (associated(self%percolationObjects_)) call percolationObjectsDeepCopyFinalize_(self%percolationObjects_)
    return
  end subroutine percolationDeepCopyFinalize

  subroutine percolationDeepCopy(self,destination)
    !!{RST
    Perform a deep copy of the object.
    !!}
    use :: Functions_Global  , only : percolationObjectsDeepCopy_
    use :: Error             , only : Error_Report
    implicit none
    class(virialDensityContrastPercolation), intent(inout), target :: self
    class(virialDensityContrastClass      ), intent(inout)         :: destination

    call self%virialDensityContrastClass%deepCopy(destination)
    select type (destination)
    type is (virialDensityContrastPercolation)
       destination%linkingLength                  =self%linkingLength
       destination%fileName                       =self%fileName
       destination%densityContrastTableTimeMinimum=self%densityContrastTableTimeMinimum
       destination%densityContrastTableTimeMaximum=self%densityContrastTableTimeMaximum
       destination%densityContrastTableMassMinimum=self%densityContrastTableMassMinimum
       destination%densityContrastTableMassMaximum=self%densityContrastTableMassMaximum
       destination%densityContrastTableMassCount  =self%densityContrastTableMassCount
       destination%densityContrastTableTimeCount  =self%densityContrastTableTimeCount
       destination%densityContrastTableInitialized=self%densityContrastTableInitialized
       destination%densityContrastTable           =self%densityContrastTable
       destination%densityContrastTableRemakeCount=self%densityContrastTableRemakeCount
       if (associated(self%cosmologyFunctions_)) then
          if (associated(self%cosmologyFunctions_%copiedSelf)) then
             select type(s => self%cosmologyFunctions_%copiedSelf)
                class is (cosmologyFunctionsClass)
                destination%cosmologyFunctions_ => s
             class default
                call Error_Report('copiedSelf has incorrect type'//{introspection:location})
             end select
             call self%cosmologyFunctions_%copiedSelf%referenceCountIncrement()
          else
             allocate(destination%cosmologyFunctions_,mold=self%cosmologyFunctions_)
             call self%cosmologyFunctions_%deepCopy(destination%cosmologyFunctions_)
             self%cosmologyFunctions_%copiedSelf => destination%cosmologyFunctions_
             call destination%cosmologyFunctions_%autoHook()
          end if
       end if
       nullify(destination%percolationobjects_)
       if (associated(self%percolationobjects_)) then
          allocate(destination%percolationobjects_,mold=self%percolationobjects_)
          call percolationObjectsDeepCopy_(self%percolationobjects_,destination%percolationobjects_)
       end if
       call destination%densitycontrasttable%deepCopyActions()
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine percolationDeepCopy


