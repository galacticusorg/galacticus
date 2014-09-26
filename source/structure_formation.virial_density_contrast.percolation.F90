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

  !% An implementation of dark matter halo virial density contrasts based on the percolation analysis of \cite{more_overdensity_2011}.

  use FGSL

  !# <virialDensityContrast name="virialDensityContrastPercolation">
  !#  <description>Dark matter halo virial density contrasts based on the percolation analysis of \cite{more_overdensity_2011}.</description>
  !# </virialDensityContrast>

  type, extends(virialDensityContrastClass) :: virialDensityContrastPercolation
     !% A dark matter halo virial density contrast class based on the percolation analysis of \cite{more_overdensity_2011}.
     private
     double precision                                                   :: linkingLength
     logical                                                            :: solving
     double precision                   , pointer                       :: densityContrastCurrent
     ! Tabulation of density contrast vs. time and mass.
     double precision                                                   :: densityContrastTableTimeMinimum                    
     double precision                                                   :: densityContrastTableTimeMaximum                    
     double precision                                                   :: densityContrastTableMassMinimum                     
     double precision                                                   :: densityContrastTableMassMaximum                     
     logical                                                            :: densityContrastTableInitialized                      
     integer                                                            :: densityContrastTableMassCount                       , densityContrastTableTimeCount
     double precision                   , allocatable, dimension(:  )   :: densityContrastTableMass                            , densityContrastTableTime
     double precision                   , allocatable, dimension(:,:)   :: densityContrastTable
     type            (fgsl_interp      )                                :: densityContrastTableTimeInterpolationObject
     type            (fgsl_interp_accel)                                :: densityContrastTableMassInterpolationAccelerator    , densityContrastTableTimeInterpolationAccelerator
     logical                                                            :: densityContrastTableMassInterpolationReset          , densityContrastTableTimeInterpolationReset            
   contains
     final     ::                                percolationDestructor
     procedure :: densityContrast             => percolationDensityContrast
     procedure :: densityContrastRateOfChange => percolationDensityContrastRateOfChange
     procedure :: isMassDependent             => percolationIsMassDepdendent
     procedure :: tabulate                    => percolationTabulate
  end type virialDensityContrastPercolation

  interface virialDensityContrastPercolation
     !% Constructors for the {\tt percolation} dark matter halo virial density contrast class.
     module procedure percolationDefaultConstructor
     module procedure percolationConstructor
  end interface virialDensityContrastPercolation

  ! Initialization state.
  logical          :: percolationInitialized=.false.

  ! Default value of the linking length and density ratio parameters.
  double precision :: virialDensityContrastPercolationLinkingLength, virialDensityContrastPercolationDensityRatio

  ! Granularity parameters for tabulations.
  integer, parameter :: percolationDensityContrastTableTimePointsPerDecade=10
  integer, parameter :: percolationDensityContrastTableMassPointsPerDecade=10

contains

  function percolationDefaultConstructor()
    !% Default constructor for the {\tt percolation} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type (virialDensityContrastPercolation), target  :: percolationDefaultConstructor
    
    if (.not.percolationInitialized) then
       !$omp critical(virialDensityContrastPercolationInitialize)
       if (.not.percolationInitialized) then
          ! Get the linking length to use.
          !@ <inputParameter>
          !@   <name>virialDensityContrastPercolationLinkingLength</name>
          !@   <defaultValue>0.2</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The friends-of-friends linking length to use in computing virial density contrasts with the percolation analysis of \cite{more_overdensity_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("virialDensityContrastPercolationLinkingLength",virialDensityContrastPercolationLinkingLength,defaultValue=0.2d0)
          ! Record initialization.
          percolationInitialized=.true.
       end if
       !$omp end critical(virialDensityContrastPercolationInitialize)
    end if
    percolationDefaultConstructor=percolationConstructor(virialDensityContrastPercolationLinkingLength)
    return
  end function percolationDefaultConstructor

  function percolationConstructor(linkingLength)
    !% Generic constructor for the {\tt percolation} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type            (virialDensityContrastPercolation), target        :: percolationConstructor
    double precision                                  , intent(in   ) :: linkingLength

    allocate(percolationConstructor%densityContrastCurrent)
    percolationConstructor%densityContrastCurrent=-1.0d0
    percolationConstructor%solving               =.false.
    percolationConstructor%linkingLength         =linkingLength
    ! Initialize tabulations.
    percolationConstructor%densityContrastTableTimeMinimum           = 1.0d-3
    percolationConstructor%densityContrastTableTimeMaximum           =20.0d+0
    percolationConstructor%densityContrastTableMassMinimum           = 1.0d8
    percolationConstructor%densityContrastTableMassMaximum           = 1.0d16
    percolationConstructor%densityContrastTableInitialized           =.false.
    percolationConstructor%densityContrastTableMassInterpolationReset=.true.
    percolationConstructor%densityContrastTableTimeInterpolationReset=.true.
    return
  end function percolationConstructor

  subroutine percolationDestructor(self)
    !% Destructor for the {\tt percolation} virial density contrast class.
    use Numerical_Interpolation
    implicit none
    type(virialDensityContrastPercolation), intent(inout) :: self

    call Interpolate_Done(                                                                                &
         &                interpolationObject     =self%densityContrastTableTimeInterpolationObject     , &
         &                interpolationAccelerator=self%densityContrastTableTimeInterpolationAccelerator, &
         &                reset                   =self%densityContrastTableTimeInterpolationReset        &
         &               )
    call Interpolate_Done(                                                                                &
         &                interpolationAccelerator=self%densityContrastTableMassInterpolationAccelerator, &
         &                reset                   =self%densityContrastTableMassInterpolationReset        &
         &                )
    return
  end subroutine percolationDestructor

  subroutine percolationTabulate(self,mass,time)
    !% Tabulate virial density contrast as a function of mass and time for the {\tt percolation} density contrast class.
    use Functions_Global, only : Virial_Density_Contrast_Percolation_Solver_
    use Numerical_Interpolation
    use Numerical_Ranges
    use Memory_Management
    use Galacticus_Display
    implicit none
    class           (virialDensityContrastPercolation), intent(inout) :: self
    double precision                                  , intent(in   ) :: mass     , time
    integer                                                           :: iMass    , iTime    , &
         &                                                               iCount
    logical                                                           :: makeTable
    double precision                                                  :: tableMass, tableTime

    ! Always check if we need to make the table.
    makeTable=.true.
    do while (makeTable)
       ! Assume table does not need remaking.
       makeTable=.false.
       ! Check for uninitialized table.
       if (.not.self%densityContrastTableInitialized) then
          makeTable=.true.
          ! Check for mass out of range.
       else if (                                                                         &
            &   mass < self%densityContrastTableMass(                                 1) &
            &    .or.                                                                    &
            &   mass > self%densityContrastTableMass(self%densityContrastTableMassCount) &
            &  ) then
          makeTable=.true.
          ! Compute the range of tabulation and number of points to use.
          self%densityContrastTableMassMinimum=min(self%densityContrastTableMassMinimum,0.5d0*mass)
          self%densityContrastTableMassMaximum=max(self%densityContrastTableMassMaximum,2.0d0*mass)
       ! Check for time out of range.
       else if (                                                                         &
            &   time < self%densityContrastTableTime(                                 1) &
            &    .or.                                                                    &
            &   time > self%densityContrastTableTime(self%densityContrastTableTimeCount) &
            &  ) then
          makeTable=.true.
          self%densityContrastTableTimeMinimum=min(self%densityContrastTableTimeMinimum,0.5d0*time)
          self%densityContrastTableTimeMaximum=max(self%densityContrastTableTimeMaximum,2.0d0*time)
       end if
       ! Remake the table if necessary.
       if (makeTable) then
          ! Record that we are in the solving phase of calculation, so we will avoid recursive calls to this function.
          self%solving=.true.       
          ! Allocate arrays to the appropriate sizes.
          self%densityContrastTableMassCount=int(log10(self%densityContrastTableMassMaximum/self%densityContrastTableMassMinimum) &
               &*dble(percolationDensityContrastTableMassPointsPerDecade))+1
          self%densityContrastTableTimeCount=int(log10(self%densityContrastTableTimeMaximum/self%densityContrastTableTimeMinimum) &
               &*dble(percolationDensityContrastTableTimePointsPerDecade))+1
          if (allocated(self%densityContrastTableMass)) call Dealloc_Array(self%densityContrastTableMass)
          if (allocated(self%densityContrastTableTime)) call Dealloc_Array(self%densityContrastTableTime)
          if (allocated(self%densityContrastTable    )) call Dealloc_Array(self%densityContrastTable    )
          call Alloc_Array(self%densityContrastTableMass,[                                   self%densityContrastTableMassCount])
          call Alloc_Array(self%densityContrastTableTime,[self%densityContrastTableTimeCount                                   ])
          call Alloc_Array(self%densityContrastTable    ,[self%densityContrastTableTimeCount,self%densityContrastTableMassCount])
          ! Create ranges of mass and time.
          self%densityContrastTableMass=Make_Range(self%densityContrastTableMassMinimum,self%densityContrastTableMassMaximum &
               &,self%densityContrastTableMassCount,rangeType=rangeTypeLogarithmic)
          self%densityContrastTableTime=Make_Range(self%densityContrastTableTimeMinimum,self%densityContrastTableTimeMaximum &
               &,self%densityContrastTableTimeCount,rangeType=rangeTypeLogarithmic)
          ! Tabulate the density contrast.
          call Galacticus_Display_Indent('Tabulating virial density contrasts for percolation class',verbosity=verbosityWorking)
          iCount=0
          do iMass=1,self%densityContrastTableMassCount
             tableMass=self%densityContrastTableMass(iMass)
             do iTime=1,self%densityContrastTableTimeCount
                tableTime=self%densityContrastTableTime(iTime)
                iCount=iCount+1
                call Galacticus_Display_Counter(int(100.0d0*dble(iCount)/dble(self%densityContrastTableMassCount*self%densityContrastTableTimeCount)),isNew=(iCount==1),verbosity=verbosityWorking)
               self%densityContrastTable(iTime,iMass)=Virial_Density_Contrast_Percolation_Solver_(tableMass,tableTime,self%linkingLength,self%densityContrastCurrent)
            end do
          end do
          call Galacticus_Display_Counter_Clear(verbosity=verbosityWorking)
          call Galacticus_Display_Unindent('done',verbosity=verbosityWorking)
          ! Reset interpolators.
          call Interpolate_Done(self%densityContrastTableTimeInterpolationObject,self%densityContrastTableTimeInterpolationAccelerator &
               &,self%densityContrastTableTimeInterpolationReset)
          call Interpolate_Done(interpolationAccelerator=self%densityContrastTableMassInterpolationAccelerator &
               &,reset=self%densityContrastTableMassInterpolationReset)
          self%densityContrastTableTimeInterpolationReset=.true.
          self%densityContrastTableMassInterpolationReset=.true.
          ! Flag that the table is now initialized.
          self%densityContrastTableInitialized=.true.
          ! Solving phase is finished.
          self%solving=.false.
       end if
    end do
    return
  end subroutine percolationTabulate
  
  double precision function percolationDensityContrast(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    use Cosmology_Functions
    use Numerical_Interpolation
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)            :: self
    double precision                                  , intent(in   )            :: mass
    double precision                                  , intent(in   ) , optional :: time               , expansionFactor
    logical                                           , intent(in   ) , optional :: collapsing
    class           (cosmologyFunctionsClass         ), pointer                  :: cosmologyFunctions_
    integer                                           , dimension(0:1)           :: jMass
    double precision                                  , dimension(0:1)           :: hMass
    double precision                                                             :: timeActual
    integer                                                                      :: iMass

    if (self%solving) then
       percolationDensityContrast=self%densityContrastCurrent
    else
       ! Get the time to use.
       cosmologyFunctions_ => cosmologyFunctions           (                               )
       timeActual          =  cosmologyFunctions_%epochTime(time,expansionFactor,collapsing)
       ! Ensure tabulation is built.
       call self%tabulate(mass,timeActual)  
       ! Get interpolating factors in mass.
       jMass(0)=Interpolate_Locate(self%densityContrastTableMassCount,self%densityContrastTableMass&
            &,self%densityContrastTableMassInterpolationAccelerator,mass,reset=self%densityContrastTableMassInterpolationReset)
       jMass(1)=jMass(0)+1
       hMass=Interpolate_Linear_Generate_Factors(self%densityContrastTableMassCount,self%densityContrastTableMass,jMass(0),mass)       
       ! Interpolate in time to get density contrast.
       percolationDensityContrast=0.0d0
       do iMass=0,1
          percolationDensityContrast=                                                &
               &  percolationDensityContrast                                         &
               & +Interpolate(self%densityContrastTableTimeCount                   , &
               &              self%densityContrastTableTime                        , &
               &              self%densityContrastTable(:,jMass(iMass))            , &
               &              self%densityContrastTableTimeInterpolationObject     , &
               &              self%densityContrastTableTimeInterpolationAccelerator, &
               &              timeActual                                           , &
               &              reset=self%densityContrastTableTimeInterpolationReset  &
               &             )                                                       &
               & *hMass(iMass)
       end do
    end if
    return
  end function percolationDensityContrast

  double precision function percolationDensityContrastRateOfChange(self,mass,time,expansionFactor,collapsing)
    !% Return the virial density contrast at the given epoch, based on the percolation algorithm of \cite{more_overdensity_2011}.
    use Cosmology_Functions
    use Numerical_Interpolation
    implicit none
    class           (virialDensityContrastPercolation), intent(inout)           :: self
    double precision                                  , intent(in   )           :: mass
    double precision                                  , intent(in   ), optional :: time               , expansionFactor
    logical                                           , intent(in   ), optional :: collapsing
    class           (cosmologyFunctionsClass         ), pointer                 :: cosmologyFunctions_
    integer                                           , dimension(0:1)          :: jMass
    double precision                                  , dimension(0:1)          :: hMass
    double precision                                                            :: timeActual
    integer                                                                     :: iMass

    ! Get the time to use.
    cosmologyFunctions_ => cosmologyFunctions           (                               )
    timeActual          =  cosmologyFunctions_%epochTime(time,expansionFactor,collapsing)
    ! Ensure tabulation is built.
    call self%tabulate(mass,timeActual)  
    ! Get interpolating factors in mass.
    jMass(0)=Interpolate_Locate(self%densityContrastTableMassCount,self%densityContrastTableMass&
         &,self%densityContrastTableMassInterpolationAccelerator,mass,reset=self%densityContrastTableMassInterpolationReset)
    jMass(1)=jMass(0)+1
    hMass=Interpolate_Linear_Generate_Factors(self%densityContrastTableMassCount,self%densityContrastTableMass,jMass(0),mass)       
    ! Interpolate in time to get density contrast.
    percolationDensityContrastRateOfChange=0.0d0
    do iMass=0,1
       percolationDensityContrastRateOfChange=                                               &
            &  percolationDensityContrastRateOfChange                                        &
            & +Interpolate_Derivative(self%densityContrastTableTimeCount                   , &
            &                         self%densityContrastTableTime                        , &
            &                         self%densityContrastTable(:,jMass(iMass))            , &
            &                         self%densityContrastTableTimeInterpolationObject     , &
            &                         self%densityContrastTableTimeInterpolationAccelerator, &
            &                         timeActual                                           , &
            &                         reset=self%densityContrastTableTimeInterpolationReset  &
            &                        )                                                       &
            & *hMass(iMass)
    end do
    return
  end function percolationDensityContrastRateOfChange

  logical function percolationIsMassDepdendent(self)
    !% Specify that the {\tt percolation} virial density contrast class is mass-dependent.
    implicit none
    class(virialDensityContrastPercolation), intent(inout) :: self
    
    percolationIsMassDepdendent=.true.
    return
  end function percolationIsMassDepdendent
  
