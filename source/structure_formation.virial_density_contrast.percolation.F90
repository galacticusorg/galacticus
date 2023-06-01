!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023
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
  An implementation of dark matter halo virial density contrasts based on the percolation analysis of \cite{more_overdensity_2011}.
  !!}

  use    :: Cosmology_Functions, only : cosmologyFunctionsClass
  use    :: Kind_Numbers       , only : kind_int8
  !$ use :: OMP_Lib            , only : OMP_Destroy_Lock       , OMP_Init_Lock, OMP_Set_Lock, OMP_Unset_Lock, &
  !$      &                             omp_lock_kind
  use    :: Tables             , only : table2DLogLogLin

  !![
  <virialDensityContrast name="virialDensityContrastPercolation" recursive="yes">
   <description>
    A dark matter halo virial density contrast class based on the percolation analysis of \cite{more_overdensity_2011}. The
    virial density contrast is computed to be consistent with a given friends-of-friends algorithm linking length using the
    percolation-theory-motivated calibration of \cite{more_overdensity_2011}. Specifically, the friends-of-friends algorithm is
    assumed to link together particles forming an isodensity surface of density $\rho = \bar{\rho} n_\mathrm{c}/b^3$, where
    $\bar{\rho}$ is the mean density of the universe, $n_\mathrm{c}=0.652960$ is a critical density for percolation as given by
    \cite{more_overdensity_2011}, and $b$ is the linking length. Given this bounding density, the virial density contrast is
    found by requiring that the halo contain the required mass within such a bounding density, given the halo density profile.
   </description>
   <deepCopy>
    <ignore   variables="recursiveSelf"/>
   </deepCopy>
   <stateStorable>
     <exclude variables="recursiveSelf"/>
   </stateStorable>
  </virialDensityContrast>
  !!]
  type, extends(virialDensityContrastClass) :: virialDensityContrastPercolation
     !!{
     A dark matter halo virial density contrast class based on the percolation analysis of \cite{more_overdensity_2011}.
     !!}
     private
     double precision                                            :: linkingLength
     type            (varying_string                  )          :: fileName
     class           (cosmologyFunctionsClass         ), pointer :: cosmologyFunctions_             => null()
     logical                                                     :: isRecursive                               , parentDeferred
     class           (virialDensityContrastPercolation), pointer :: recursiveSelf                   => null()
     class           (*                               ), pointer :: percolationObjects_             => null()
     integer         (kind_int8                       )          :: selfID
     ! Tabulation of density contrast vs. time and mass.
     double precision                                            :: densityContrastTableTimeMinimum           , densityContrastTableTimeMaximum
     double precision                                            :: densityContrastTableMassMinimum           , densityContrastTableMassMaximum
     integer                                                     :: densityContrastTableMassCount             , densityContrastTableTimeCount
     logical                                                     :: densityContrastTableInitialized =  .false.
     logical                                                     :: densityContrastLockInitialized  =  .false.
     type            (table2DLogLogLin                )          :: densityContrastTable
     !$ integer      (omp_lock_kind                   )          :: densityContrastTableLock
     integer                                                     :: densityContrastTableRemakeCount
   contains
     !![
     <methods>
       <method description="Tabulate the virial density contrast as a function of mass and time." method="tabulate"    />
       <method description="Restore a tabulated solution from file."                              method="restoreTable"/>
       <method description="Store a tabulated solution to file."                                  method="storeTable"  />
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
     procedure :: restoreTable                => percolationRestoreTable
     procedure :: storeTable                  => percolationStoreTable
  end type virialDensityContrastPercolation

  interface virialDensityContrastPercolation
     !!{
     Constructors for the {\normalfont \ttfamily percolation} dark matter halo virial density contrast class.
     !!}
     module procedure percolationConstructorParameters
     module procedure percolationConstructorInternal
  end interface virialDensityContrastPercolation

  ! Granularity parameters for tabulations.
  integer                    , parameter :: densityContrastTableTimePointsPerDecade=5
  integer                    , parameter :: densityContrastTableMassPointsPerDecade=5

  ! Module-scope record of state used when solving for the percolation density contrast.
  logical                                :: solving                                =.false.
  double precision                       :: densityContrastCurrent                 =-1.0d0
  !$omp threadprivate(densityContrastCurrent,solving)

  ! Unique ID used to identify parent-child recursive relationships.
  integer         (kind_int8)            :: ID                                     =0_kind_int8, IDRecursive=-1_kind_int8
  !$omp threadprivate(IDRecursive)

contains

  recursive function percolationConstructorParameters(parameters,recursiveConstruct,recursiveSelf) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily percolation} dark matter halo virial density contrast class that takes a parameter set as input.
    !!}
    use :: Cosmology_Functions, only : cosmologyFunctions                                      , cosmologyFunctionsClass
    use :: Functions_Global   , only : Virial_Density_Contrast_Percolation_Objects_Constructor_
    use :: Input_Parameters   , only : inputParameter                                          , inputParameters
    implicit none
    type            (virialDensityContrastPercolation)                          :: self
    type            (inputParameters                 ), intent(inout)           :: parameters
    logical                                           , intent(in   ), optional :: recursiveConstruct
    class           (virialDensityContrastClass      ), intent(in   ), optional :: recursiveSelf
    class           (cosmologyFunctionsClass         ), pointer                 :: cosmologyFunctions_
    class           (*                               ), pointer                 :: percolationObjects_
    double precision                                                            :: linkingLength
    !![
    <optionalArgument name="recursiveConstruct" defaultsTo=".false." />
    !!]

    !![
    <inputParameter>
     <name>linkingLength</name>
     <source>parameters</source>
     <defaultValue>0.2d0</defaultValue>
     <description>The friends-of-friends linking length to use in computing virial density contrasts with the percolation analysis of \cite{more_overdensity_2011}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyFunctions" name="cosmologyFunctions_" source="parameters"/>
    !!]
    ! Build a pointer to a container object which stores all objects needed by the percolation density solver. If this construct
    ! has been called recursively, then instead set this to null - recursive construction can happen as some objects used in
    ! solving for percolation density (e.g. dark matter profiles) themselves require a virial density contrast object. By setting
    ! this pointer to null we avoid circular dependencies between functionClass objects which lead to memory leaks.
    if (recursiveConstruct_) then
       percolationObjects_ => null()
    else
       percolationObjects_ => Virial_Density_Contrast_Percolation_Objects_Constructor_(parameters)
    end if
    self=virialDensityContrastPercolation(linkingLength,cosmologyFunctions_,percolationObjects_,recursiveSelf)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyFunctions_"/>
    !!]
    nullify(percolationObjects_)
    return
  end function percolationConstructorParameters

  recursive  function percolationConstructorInternal(linkingLength,cosmologyFunctions_,percolationObjects_,recursiveSelf) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily percolation} dark matter halo virial density contrast class.
    !!}
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath     , pathTypeDataDynamic
    use :: ISO_Varying_String, only : operator(//)
    use :: Input_Parameters  , only : inputParameter, inputParameters
    implicit none
    type            (virialDensityContrastPercolation)                                  :: self
    double precision                                  , intent(in   )                   :: linkingLength
    class           (cosmologyFunctionsClass         ), intent(in   ), target           :: cosmologyFunctions_
    class           (*                               ), intent(in   ), target           :: percolationObjects_
    class           (virialDensityContrastClass      ), intent(in   ), target, optional :: recursiveSelf
    !![
    <constructorAssign variables="linkingLength, *cosmologyFunctions_, *percolationObjects_"/>
    !!]

    ! File name for tabulation.
    self%fileName=inputPath(pathTypeDataDynamic)                   // &
         &        'darkMatterHalos/'                               // &
         &        self%objectType      (                          )// &
         &        '_'                                              // &
         &        self%hashedDescriptor(includeSourceDigest=.true.)// &
         &        '.hdf5'
    ! Initialize tabulations.
    self%densityContrastTableTimeMinimum= 1.0d-03
    self%densityContrastTableTimeMaximum=20.0d+00
    self%densityContrastTableMassMinimum= 4.0d+05
    self%densityContrastTableMassMaximum= 1.0d+16
    self%densityContrastTableInitialized=.false.
    self%densityContrastTableRemakeCount= 0
    self%densityContrastLockInitialized =.true.
    !$ call OMP_Init_Lock(self%densityContrastTableLock)
    ! Check if we have been given a percolation objects container. If not, we expect that we are being built recursively, in which
    ! case we require a pointer to the originating object under construction.
    if (.not.associated(self%percolationObjects_)) then
       if (.not.present(recursiveSelf)) call Error_Report('for recursive construction a pointer to the originating object is required'//{introspection:location})
       select type (recursiveSelf)
       class is (virialDensityContrastPercolation)
          self%recursiveSelf => recursiveSelf
          self%isRecursive   =  .true.
       class default
          call Error_Report('originating object is of incorrect class'//{introspection:location})
       end select
       ! Generate a new ID.
       !$omp critical (percolationIDIncrement)
       ID                =ID+1_kind_int8
       self                  %selfID=ID
       if (ID == -1_kind_int8) call Error_Report('ran out of IDs for percolation class'//{introspection:location})
       !$omp end critical (percolationIDIncrement)
       IDRecursive       =self%selfID
    else
       self%recursiveSelf => null()
       self%isRecursive   =  .false.
       ! Get a unique ID.
       if (IDRecursive == -1_kind_int8) then
          ! Generate a new ID.
          !$omp critical (percolationIDIncrement)
          ID                =ID+1_kind_int8
          self                  %selfID=ID
          if (ID == -1_kind_int8) call Error_Report('ran out of IDs for percolation class'//{introspection:location})
          !$omp end critical (percolationIDIncrement)
       else
          ! An ID was already set during recursive construction.
          self                  %selfID=IDRecursive
          IDRecursive       =-1_kind_int8
       end if
    end if
    self%parentDeferred=.false.
    return
  end function percolationConstructorInternal

  subroutine percolationDestructor(self)
    !!{
    Destructor for the {\normalfont \ttfamily percolation} virial density contrast class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type(virialDensityContrastPercolation), intent(inout) :: self

    if    (self%densityContrastTableInitialized) call                  self%densityContrastTable    %destroy()
    !$ if (self%densityContrastLockInitialized ) call OMP_Destroy_Lock(self%densityContrastTableLock          )
    !![
    <objectDestructor name="self%cosmologyFunctions_" />
    !!]
    if (associated(self%percolationObjects_)) deallocate(self%percolationObjects_)
    return
  end subroutine percolationDestructor

  subroutine percolationTabulate(self,mass,time,mustRetabulate)
    !!{
    Tabulate virial density contrast as a function of mass and time for the {\normalfont \ttfamily percolation} density contrast class.
    !!}
    use :: Display         , only : displayCounter                             , displayCounterClear, displayIndent, displayUnindent, &
          &                         verbosityLevelWorking
    use :: Functions_Global, only : Virial_Density_Contrast_Percolation_Solver_
    use :: Error           , only : Error_Report
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass          , time
    logical                                           , intent(  out), optional :: mustRetabulate
    integer                                                                     :: iMass         , iTime    , &
         &                                                                         iCount
    logical                                                                     :: makeTable
    double precision                                                            :: tableMass     , tableTime

    ! Always check if we need to make the table.
    makeTable=.true.
    do while (makeTable)
       ! Assume table does not need remaking.
       makeTable=.false.
       ! Check for uninitialized table.
       if (.not.self%densityContrastTableInitialized) call self%restoreTable()
       if (.not.self%densityContrastTableInitialized) then
          makeTable=.true.
       ! Check for mass out of range.
       else if (                                            &
            &   mass < self%densityContrastTableMassMinimum &
            &    .or.                                       &
            &   mass > self%densityContrastTableMassMaximum &
            &  ) then
          makeTable=.true.
          ! Compute the range of tabulation and number of points to use.
          self%densityContrastTableMassMinimum=min(self%densityContrastTableMassMinimum,0.5d0*mass)
          self%densityContrastTableMassMaximum=max(self%densityContrastTableMassMaximum,2.0d0*mass)
       ! Check for time out of range.
       else if (                                            &
            &   time < self%densityContrastTableTimeMinimum &
            &    .or.                                       &
            &   time > self%densityContrastTableTimeMaximum &
            &  ) then
          makeTable=.true.
          self%densityContrastTableTimeMinimum=min(self%densityContrastTableTimeMinimum,0.5d0*time)
          self%densityContrastTableTimeMaximum=max(self%densityContrastTableTimeMaximum,2.0d0*time)
       end if
       ! If we are just being asked to determine remake status, return that now.
       if (present(mustRetabulate)) then
          mustRetabulate=makeTable
          return
       end if
       ! Remake the table is necessary.
       if (makeTable) then
          ! Check that we have a pointer to the required objects.
          if (.not.associated(self%percolationObjects_)) call Error_Report('no percolationObjects available'//{introspection:location})
          ! Increment the number of table remakes.
          self%densityContrastTableRemakeCount=self%densityContrastTableRemakeCount+1
          ! Record that we are in the solving phase of calculation, so we will avoid recursive calls to this function.
          solving=.true.
          ! Allocate arrays to the appropriate sizes.
          self%densityContrastTableMassCount=int(log10(self%densityContrastTableMassMaximum/self%densityContrastTableMassMinimum)*dble(densityContrastTableMassPointsPerDecade))+1
          self%densityContrastTableTimeCount=int(log10(self%densityContrastTableTimeMaximum/self%densityContrastTableTimeMinimum)*dble(densityContrastTableTimePointsPerDecade))+1
          ! Create the table.
          call self%densityContrastTable%create(                                      &
               &                                self%densityContrastTableMassMinimum, &
               &                                self%densityContrastTableMassMaximum, &
               &                                self%densityContrastTableMassCount  , &
               &                                self%densityContrastTableTimeMinimum, &
               &                                self%densityContrastTableTimeMaximum, &
               &                                self%densityContrastTableTimeCount    &
               &                               )
          ! Tabulate the density contrast.
          call displayIndent('Tabulating virial density contrasts for percolation class',verbosity=verbosityLevelWorking)
          iCount=0
          do iMass=1,self%densityContrastTableMassCount
             tableMass=self%densityContrastTable%x(iMass)
             do iTime=1,self%densityContrastTableTimeCount
                tableTime=self%densityContrastTable%y(iTime)
                iCount=iCount+1
                call displayCounter(int(100.0d0*dble(iCount)/dble(self%densityContrastTableMassCount*self%densityContrastTableTimeCount)),isNew=(iCount==1),verbosity=verbosityLevelWorking)
                call self%densityContrastTable%populate(Virial_Density_Contrast_Percolation_Solver_(tableMass,tableTime,self%linkingLength,densityContrastCurrent,self%percolationObjects_,self),iMass,iTime)
             end do
          end do
          call displayCounterClear(verbosity=verbosityLevelWorking)
          call displayUnindent('done',verbosity=verbosityLevelWorking)
          ! Flag that the table is now initialized.
          self%densityContrastTableInitialized=.true.
          ! Store the table.
          call self%storeTable()
          ! Solving phase is finished.
          solving=.false.
       end if
    end do
    return
  end subroutine percolationTabulate

  double precision function percolationDensityContrast(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)            :: self
    double precision                                  , intent(in   )            :: mass
    double precision                                  , intent(in   ) , optional :: time            , expansionFactor
    logical                                           , intent(in   ) , optional :: collapsing
    double precision                                                             :: timeActual
    logical                                                                      :: useSolverCurrent, useTable       , &
         &                                                                          mustRetabulate

    ! Get the time to use.
    if (.not.solving) timeActual=self%cosmologyFunctions_%epochTime(time,expansionFactor,collapsing)
    ! Initialize choices for how we will compute the density contrast.
    useSolverCurrent=.false.
    useTable        =.false.
    ! Determine how to compute density contrast.
    if (self%isRecursive) then
       if (solving) then
          ! Currently solving for solutions - return the current guess.
          if (self%selfID /= self%recursiveSelf%selfID) call Error_Report('recursively-constructed percolation class object ID does not match that of actively solving object'//{introspection:location})
          useSolverCurrent=.true.
       else
          ! Not solving for solutions - if our table is sufficient (ensure it is up to date first), use it, otherwise request that
          ! the parent self performs the calculation (in which case it will update its table, and we can later copy it).
          call percolationCopyTable(self)
          !$ call OMP_Set_Lock(self%densityContrastTableLock)
          call self%tabulate(mass,timeActual,mustRetabulate)
          !$ call OMP_Unset_Lock(self%densityContrastTableLock)
          if (mustRetabulate) then
             percolationDensityContrast=self%recursiveSelf%densityContrast(mass,time,expansionFactor,collapsing)
             return
          else
             useTable=.true.
          end if
       end if
    else if (solving) then
       ! Currently solving for solutions - return the current guess.
       useSolverCurrent=.true.
    else
       ! Use our tabulated solution.
       useTable        =.true.
    end if
    ! Provide a result based on the chosen method.
    if (useSolverCurrent) then
       percolationDensityContrast=densityContrastCurrent
    else if (useTable) then
       !$ call OMP_Set_Lock(self%densityContrastTableLock)
       call self%tabulate(mass,timeActual)
       percolationDensityContrast=self%densityContrastTable%interpolate(mass,timeActual)
       !$ call OMP_Unset_Lock(self%densityContrastTableLock)
    else
       call Error_Report('no method selected to compute density contrast'//{introspection:location})
    end if
    return
  end function percolationDensityContrast

  double precision function percolationDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !!{
    Return the rate of change of the virial density contrast at the given epoch, based on the percolation algorithm of
    \cite{more_overdensity_2011}.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)            :: self
    double precision                                  , intent(in   )            :: mass
    double precision                                  , intent(in   ) , optional :: time          , expansionFactor
    logical                                           , intent(in   ) , optional :: collapsing
    double precision                                                             :: timeActual
    logical                                                                      :: mustRetabulate

    ! Get the time to use.
    timeActual=self%cosmologyFunctions_%epochTime(time,expansionFactor,collapsing)
    ! Determine how to compute density contrast.
    if (self%isRecursive) then
       call percolationCopyTable (self)
       ! If our table is sufficient, use it, otherwise request that the parent self performs the calculation (in which case it
       ! will update its table, and we can later copy it).
       !$ call OMP_Set_Lock(self%densityContrastTableLock)
       call self%tabulate(mass,timeActual,mustRetabulate)
       !$ call OMP_Unset_Lock(self%densityContrastTableLock)
       if (mustRetabulate) then
          percolationDensityContrastRateOfChange=self%recursiveSelf%densityContrastRateOfChange(mass,time,expansionFactor,collapsing)
          return
       end if
    end if
    ! Compute the solution.
    !$ call OMP_Set_Lock(self%densityContrastTableLock)
    call self%tabulate(mass,timeActual)
    percolationDensityContrastRateOfChange=self%densityContrastTable%interpolateGradient(mass,timeActual,dim=2)
    !$ call OMP_Unset_Lock(self%densityContrastTableLock)
    return
  end function percolationDensityContrastRateOfChange

  logical function percolationIsMassDepdendent(self)
    !!{
    Specify that the {\normalfont \ttfamily percolation} virial density contrast class is mass-dependent.
    !!}
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self
    !$GLC attributes unused :: self

    percolationIsMassDepdendent=.true.
    return
  end function percolationIsMassDepdendent

  subroutine percolationDeepCopyReset(self)
    !!{
    Perform a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : percolationObjectsDeepCopyReset_
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self

    self                           %   copiedSelf => null()
    if (.not.self%isRecursive) self%recursiveSelf => null()
    if (associated(self%cosmologyfunctions_)) call self%cosmologyfunctions_%deepCopyReset()
    if (associated(self%percolationObjects_)) call percolationObjectsDeepCopyReset_(self%percolationObjects_)
    return
  end subroutine percolationDeepCopyReset

  subroutine percolationDeepCopyFinalize(self)
    !!{
    Finalize a deep copy reset of the object.
    !!}
    use :: Functions_Global, only : percolationObjectsDeepCopyFinalize_
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self

    if (self%isRecursive) call percolationFindParent(self)
    if (associated(self%cosmologyfunctions_)) call self%cosmologyfunctions_%deepCopyFinalize()
    if (associated(self%percolationObjects_)) call percolationObjectsDeepCopyFinalize_(self%percolationObjects_)
    return
  end subroutine percolationDeepCopyFinalize

  subroutine percolationDeepCopy(self,destination)
    !!{
    Perform a deep copy of the object.
    !!}
    use :: Functions_Global  , only : percolationObjectsDeepCopy_
    use :: Error             , only : Error_Report
#ifdef OBJECTDEBUG
    use :: MPI_Utilities     , only : mpiSelf
#endif
#ifdef OBJECTDEBUG
    use :: Function_Classes  , only : debugReporting
#endif
#ifdef OBJECTDEBUG
    use :: Display           , only : displayMessage             , verbosityLevelSilent
#endif
#ifdef OBJECTDEBUG
    use :: ISO_Varying_String, only : operator(//)               , var_str
#endif
#ifdef OBJECTDEBUG
    use :: String_Handling   , only : operator(//)
#endif
    implicit none
    class(virialDensityContrastPercolation), intent(inout), target :: self
    class(virialDensityContrastClass      ), intent(inout)         :: destination

    call self%virialDensityContrastClass%deepCopy(destination)
    select type (destination)
    type is (virialDensityContrastPercolation)
       destination%linkingLength                  =self%linkingLength
       destination%fileName                       =self%fileName
       destination%isRecursive                    =self%isRecursive
       destination%densityContrastTableTimeMinimum=self%densityContrastTableTimeMinimum
       destination%densityContrastTableTimeMaximum=self%densityContrastTableTimeMaximum
       destination%densityContrastTableMassMinimum=self%densityContrastTableMassMinimum
       destination%densityContrastTableMassMaximum=self%densityContrastTableMassMaximum
       destination%densityContrastTableMassCount  =self%densityContrastTableMassCount
       destination%densityContrastTableTimeCount  =self%densityContrastTableTimeCount
       destination%densityContrastTableInitialized=self%densityContrastTableInitialized
       destination%densityContrastTable           =self%densityContrastTable
       destination%densityContrastTableRemakeCount=self%densityContrastTableRemakeCount
       destination%selfID                         =self%selfID
       destination%parentDeferred                 =.false.
       if (self%isRecursive) then
          if (associated(self%recursiveSelf%recursiveSelf)) then
             ! If the parent self's recursiveSelf pointer is set, it indicates that it was deep-copied, and the pointer points to
             ! that copy. In that case we set the parent self of our destination to that copy.
             destination%recursiveSelf  => self%recursiveSelf%recursiveSelf
          else
             ! The parent self does not appear to have been deep-copied yet. Retain the same parent self pointer in our copy, but
             ! indicate that we need to look for the new parent later.
             destination%recursiveSelf  => self%recursiveSelf
             destination%parentDeferred =  .true.
          end if
       else
          ! This is a parent of a recursively-constructed object. Record the location of our copy so that it can be used to set
          ! the parent in deep copies of the child object.
          call percolationDeepCopyAssign(self,destination)
          destination%recursiveSelf                     => null()
       end if
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
#ifdef OBJECTDEBUG
          if (debugReporting.and.mpiSelf%isMaster()) call displayMessage(var_str('functionClass[own] (class : ownerName : ownerLoc : objectLoc : sourceLoc): cosmologyfunctions_ : [destination] : ')//loc(destination)//' : '//loc(destination%cosmologyFunctions_)//' : '//{introspection:location:compact},verbosityLevelSilent)
#endif
       end if
       nullify(destination%percolationobjects_)
       if (associated(self%percolationobjects_)) then
          allocate(destination%percolationobjects_,mold=self%percolationobjects_)
          call percolationObjectsDeepCopy_(self%percolationobjects_,destination%percolationobjects_)
       end if
       call destination%densitycontrasttable%deepCopyActions()
       destination%densityContrastLockInitialized=.true.
       !$ call OMP_Init_Lock(destination%densityContrastTableLock)
    class default
       call Error_Report('destination and source types do not match'//{introspection:location})
    end select
    return
  end subroutine percolationDeepCopy

  subroutine percolationDeepCopyAssign(self,destination)
    !!{
    Perform pointer assignment during a deep copy of the object.
    !!}
    implicit none
    class(virialDensityContrastPercolation), intent(inout)         :: self
    class(virialDensityContrastClass      ), intent(inout), target :: destination

    select type (destination)
    type is (virialDensityContrastPercolation)
       self%recursiveSelf => destination
    end select
    return
  end subroutine percolationDeepCopyAssign

  subroutine percolationFindParent(self)
    !!{
    Find the deep-copied parent of a recursive child.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self

    if (self%parentDeferred) then
       if (associated(self%recursiveSelf%recursiveSelf)) then
          self%recursiveSelf => self%recursiveSelf%recursiveSelf
       else
        call Error_Report("recursive child's parent was not copied"//{introspection:location})
       end if
       self%parentDeferred=.false.
    end if
    return
  end subroutine percolationFindParent

  subroutine percolationCopyTable(self)
    !!{
    Copy the table from a recursive child's parent.
    !!}
    use :: Error, only : Error_Report
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self

    if (associated(self%recursiveSelf)) then
       if (self%densityContrastTableRemakeCount /= self%recursiveSelf%densityContrastTableRemakeCount) then
          call self%densityContrastTable%destroy()
          self%densityContrastTableTimeMinimum=self%recursiveSelf%densityContrastTableTimeMinimum
          self%densityContrastTableTimeMaximum=self%recursiveSelf%densityContrastTableTimeMaximum
          self%densityContrastTableMassMinimum=self%recursiveSelf%densityContrastTableMassMinimum
          self%densityContrastTableMassMaximum=self%recursiveSelf%densityContrastTableMassMaximum
          self%densityContrastTableMassCount  =self%recursiveSelf%densityContrastTableMassCount
          self%densityContrastTableTimeCount  =self%recursiveSelf%densityContrastTableTimeCount
          self%densityContrastTableInitialized=self%recursiveSelf%densityContrastTableInitialized
          self%densityContrastTable           =self%recursiveSelf%densityContrastTable
          self%densityContrastTableRemakeCount=self%recursiveSelf%densityContrastTableRemakeCount
       end if
    else
       call Error_Report("recursive child has no parent"//{introspection:location})
    end if
    return
  end subroutine percolationCopyTable

  subroutine percolationStoreTable(self)
    !!{
    Store the table to file.
    !!}
    use :: File_Utilities    , only : Directory_Make, File_Lock     , File_Path, File_Unlock, &
          &                           lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char          , varying_string
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self
    type (hdf5Object                      )                :: file
    type (lockDescriptor                  )                :: fileLock

    call Directory_Make(char(File_Path(char(self%fileName))))
    ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
    call File_Lock     (char(self%fileName),fileLock,lockIsShared=.false.)
    !$ call hdf5Access%set()
    file=hdf5Object(char(self%fileName),overWrite=.true.,readOnly=.false.)
    call file%writeAttribute(        self%densityContrastTableTimeMinimum                                                                 ,'timeMinimum'    )
    call file%writeAttribute(        self%densityContrastTableTimeMaximum                                                                 ,'timeMaximum'    )
    call file%writeAttribute(        self%densityContrastTableMassMinimum                                                                 ,'massMinimum'    )
    call file%writeAttribute(        self%densityContrastTableMassMaximum                                                                 ,'massMaximum'    )
    call file%writeAttribute(        self%densityContrastTableTimeCount                                                                   ,'timeCount'      )
    call file%writeAttribute(        self%densityContrastTableMassCount                                                                   ,'massCount'      )
    call file%writeDataset  (        self%densityContrastTable%xs()                                                                       ,'mass'           )
    call file%writeDataset  (        self%densityContrastTable%ys()                                                                       ,'time'           )
    call file%writeDataset  (reshape(self%densityContrastTable%zs(),[self%densityContrastTable%size(1),self%densityContrastTable%size(2)]),'densityContrast')
    !$ call hdf5Access%unset()
    call File_Unlock(fileLock)
    return
  end subroutine percolationStoreTable

  subroutine percolationRestoreTable(self)
    !!{
    Attempt to restore a table from file.
    !!}
    use :: File_Utilities    , only : File_Exists, File_Lock     , File_Unlock, lockDescriptor
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char       , varying_string
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)               :: self
    double precision                                  , allocatable, dimension(:  ) :: timeTable           , massTable
    double precision                                  , allocatable, dimension(:,:) :: densityContrastTable
    type            (hdf5Object                      )                              :: file
    type            (lockDescriptor                  )                              :: fileLock

    if (File_Exists(self%fileName)) then
       ! Always obtain the file lock before the hdf5Access lock to avoid deadlocks between OpenMP threads.
       call File_Lock(char(self%fileName),fileLock,lockIsShared=.true.)
       !$ call hdf5Access%set()
       file=hdf5Object(char(self%fileName))
       call file%readAttribute('timeMinimum'    ,self%densityContrastTableTimeMinimum)
       call file%readAttribute('timeMaximum'    ,self%densityContrastTableTimeMaximum)
       call file%readAttribute('massMinimum'    ,self%densityContrastTableMassMinimum)
       call file%readAttribute('massMaximum'    ,self%densityContrastTableMassMaximum)
       call file%readAttribute('timeCount'      ,self%densityContrastTableTimeCount  )
       call file%readAttribute('massCount'      ,self%densityContrastTableMassCount  )
       call file%readDataset  ('mass'           ,massTable                           )
       call file%readDataset  ('time'           ,timeTable                           )
       call file%readDataset  ('densityContrast',densityContrastTable                )
       !$ call hdf5Access%unset()
       call File_Unlock(fileLock)
       call self%densityContrastTable%create  (                                       &
            &                                  massTable           (1              ), &
            &                                  massTable           (size(massTable)), &
            &                                                       size(massTable) , &
            &                                  timeTable           (1              ), &
            &                                  timeTable           (size(timeTable)), &
            &                                                       size(timeTable)   &
            &                                 )
       call self%densityContrastTable%populate(                                       &
            &                                  densityContrastTable                   &
            &                                 )
       self%densityContrastTableInitialized=.true.
       self%densityContrastTableRemakeCount=0
    end if
    return
  end subroutine percolationRestoreTable
