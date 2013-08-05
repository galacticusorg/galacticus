!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which performs analysis to compute a variety of mass functions.

module Galacticus_Output_Analyses_Mass_Functions
  !% Performs analysis to compute a variety of mass functions. Currently supported mass functions include:
  !% \begin{itemize}
  !% \item The $z\approx 0.07$ stellar mass function in the SDSS measured by \cite{li_distribution_2009}. Assumes $0.07$
  !% dex random errors on stellar masses. This is approximate, but motivated by the discussion of
  !% \cite{behroozi_comprehensive_2010}.
  !% \end{itemize}
  use Galacticus_Nodes
  use Galactic_Structure_Options
  implicit none
  private
  public :: Galacticus_Output_Analysis_Mass_Functions, Galacticus_Output_Analysis_Mass_Functions_Output

  ! Record of module initialization.
  logical                                                   :: moduleInitialized=.false.

  ! Record of whether this analysis is active.
  logical                                                   :: analysisActive

  ! Number of supported mass functions.
  integer          , parameter                              :: massFunctionsSupportedCount=9

  ! Labels for supported mass functions.
  character(len=32), dimension(massFunctionsSupportedCount) :: massFunctionLabels= &
       & [                                                                         &
       &  'sdssStellarMassFunctionZ0.07   ',                                       &
       &  'alfalfaHiMassFunctionZ0.00     ',                                       &
       &  'primusStellarMassFunctionZ0.100',                                       &
       &  'primusStellarMassFunctionZ0.250',                                       &
       &  'primusStellarMassFunctionZ0.350',                                       &
       &  'primusStellarMassFunctionZ0.450',                                       &
       &  'primusStellarMassFunctionZ0.575',                                       &
       &  'primusStellarMassFunctionZ0.725',                                       &
       &  'primusStellarMassFunctionZ0.900'                                        &
       & ]

  ! Interface for mass mapping functions.
  abstract interface
     double precision function Map_Mass(mass,thisNode)
       import treeNode
       double precision          , intent(in   )          :: mass
       type            (treeNode), intent(inout), pointer :: thisNode
     end function Map_Mass
  end interface

  ! Type for descriptors of mass functions.
  type :: massFunctionDescriptor
     double precision                    :: redshift
     double precision                    :: systematicLogM0
     double precision                    :: randomError
     double precision                    :: massLogarithmicMinimum
     integer                             :: systematicCoefficientCount
     integer                             :: massType
     character       (len= 32 )          :: label
     character       (len=128 )          :: comment
     procedure       (Map_Mass), pointer :: mapMass
  end type massFunctionDescriptor

  ! Mass function descriptors.
  type(massFunctionDescriptor), dimension(massFunctionsSupportedCount), target :: massFunctionDescriptors=    &
       & [ &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.07d00                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   2                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'sdssStellarMassFunctionZ0.07'            ,        &
       &                                                   'SDSS stellar mass function at z=0.07'    ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.000d0                                   ,        &
       &                                                  10.000d0                                   ,        &
       &                                                   0.080d0                                   ,        &
       &                                                   4.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeGaseous                           ,        &
       &                                                   'alfalfaHiMassFunctionZ0.00'              ,        &
       &                                                   'ALFALFA HI mass function at z=0.00'      ,        &
       &                                                   null()                                             &
       &                                                 )                                           ,        &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.100d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.100'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.100' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.250d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.250'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.250' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.350d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.350'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.350' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.450d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.450'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.450' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.575d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.575'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.575' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.725d0                                   ,        &
       &                                                  11.300d0                                   ,        &
       &                                                   0.070d0                                   ,        &
       &                                                   6.500d0                                   ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.725'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.725' ,        &
       &                                                   null()                                             &
       &                                                 )                                                  , &
       &                           massFunctionDescriptor(                                                    &
       &                                                   0.90d0                                    ,        &
       &                                                  11.30d0                                    ,        &
       &                                                   0.07d0                                    ,        &
       &                                                   6.50d0                                    ,        &
       &                                                   0                                         ,        &
       &                                                   massTypeStellar                           ,        &
       &                                                   'primusStellarMassFunctionZ0.900'         ,        &
       &                                                   'PRMIUS stellar mass function at z=0.900' ,        &
       &                                                   null()                                             &
       &                                                 )                                                    &
       & ]

  ! Type to store mass functions.
  type :: massFunction
     ! Copy of the mass function descriptor for this mass function.
     type            (massFunctionDescriptor), pointer                     :: descriptor
     ! Parameters for the systematic error model.
     double precision                        , allocatable, dimension(:  ) :: systematicCoefficients
     ! The index of the output corresponding to the required redshift.
     integer                                                               :: outputNumber
     ! The number of mass bins.
     integer                                                               :: massesCount
     ! Arrays for the masses and mass function.
     double precision                        , allocatable, dimension(:  ) :: masses,massesLogarithmic,massFunction,massesLogarithmicMinimum&
          &,massesLogarithmicMaximum
     ! Array for the covariance matrix.
     double precision                        , allocatable, dimension(:,:) :: massFunctionCovariance
  end type massFunction

  ! Mass functions.
  type(massFunction), allocatable, dimension(:) :: massFunctions

  ! Type for storing temporary mass functions during cumulation.
  type :: massFunctionWork
     double precision, allocatable, dimension(:  ) :: massFunction
     double precision, allocatable, dimension(:,:) :: covariance
  end type massFunctionWork

  ! Work array.
  type(massFunctionWork), allocatable, dimension(:) :: thisGalaxy
  !$omp threadprivate(thisGalaxy)

contains

  !# <mergerTreeAnalysisTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Functions</unitName>
  !# </mergerTreeAnalysisTask>
  subroutine Galacticus_Output_Analysis_Mass_Functions(thisTree,thisNode,iOutput,mergerTreeAnalyses)
    !% Construct a mass functions to compare to various observational determinations.
    use Galacticus_Nodes
    use Galacticus_Input_Paths
    use IO_HDF5
    use ISO_Varying_String
    use Memory_Management
    use Cosmological_Parameters
    use Galactic_Structure_Enclosed_Masses
    use Input_Parameters
    use Galacticus_Output_Times
    use Galacticus_Error
    use Cosmology_Functions
    use Numerical_Comparison
    use String_Handling
    use FoX_dom
    implicit none
    type            (mergerTree    ), intent(in   )                       :: thisTree
    type            (treeNode      ), intent(inout), pointer              :: thisNode
    integer                         , intent(in   )                       :: iOutput
    type            (varying_string), intent(in   ), dimension(:  )       :: mergerTreeAnalyses
    type            (node          ), pointer                             :: doc,massFunctionElement,columnElement,massElement,hubbleElement,hubbleExponentElement,datum
    type            (nodeList      ), pointer                             :: dataList
    type            (hdf5Object    )                                      :: dataFile,massDataset,parameters
    integer                                                               :: i,j,k,currentAnalysis,activeAnalysisCount,iDatum,ioErr
    double precision                                                      :: massHubbleExponent,dataHubbleParameter,hubbleRatio &
         &,mass,massLogarithmic
    type            (varying_string)                                      :: parameterName

    ! Initialize the module if necessary.
    if (.not.moduleInitialized) then
       !$omp critical(Galacticus_Output_Analysis_Mass_Functions_Initialize)
       if (.not.moduleInitialized) then
          ! Establish mass mapping functions for mass function descriptors.
          massFunctionDescriptors(2)%mapMass => Map_Mass_ALFALFA_HI_Mass_Function_Z0_00
          ! Determine how many supported mass functions are requested.
          activeAnalysisCount=0
          do i=1,massFunctionsSupportedCount
             if (any(trim(mergerTreeAnalyses) == trim(massFunctionLabels(i)))) activeAnalysisCount=activeAnalysisCount+1
          end do
          ! Allocate mass function arrays and populate with required data.
          if (activeAnalysisCount <= 0) then
             analysisActive=.false.
          else
             analysisActive=.true.
             currentAnalysis=0
             allocate(massFunctions(activeAnalysisCount))
             do i=1,size(mergerTreeAnalyses)
                do j=1,massFunctionsSupportedCount
                   if (mergerTreeAnalyses(i) == trim(massFunctionLabels(j))) then
                      currentAnalysis=currentAnalysis+1
                      ! Copy the descriptor for this mass function.
                      massFunctions(currentAnalysis)%descriptor => massFunctionDescriptors(j)
                      ! Read parameters of the systematic error model.
                      if (massFunctionDescriptors(j)%systematicCoefficientCount > 0) then
                         allocate(massFunctions(currentAnalysis)%systematicCoefficients(massFunctionDescriptors(j)%systematicCoefficientCount))
                         do k=1,massFunctionDescriptors(j)%systematicCoefficientCount
                            parameterName=trim(massFunctionLabels(j))//'MassSystematic'
                            parameterName=parameterName//(k-1)
                            call Get_Input_Parameter(char(parameterName),massFunctions(currentAnalysis)%systematicCoefficients(k),defaultValue=0.0d0)
                         end do
                      end if
                      ! Find which output number corresponds to the required redshift.
                      massFunctions(currentAnalysis)%outputNumber=-1
                      do k=1,Galacticus_Output_Time_Count()
                         if     (                                                                                        &
                              &  Values_Agree(                                                                           &
                              &               massFunctionDescriptors(j)%redshift                                      , &
                              &               Redshift_From_Expansion_Factor(                                            &
                              &                                              Expansion_Factor(                           &
                              &                                                               Galacticus_Output_Time(k)  &
                              &                                                              )                           &
                              &                                             )                                          , &
                              &               absTol=0.001d0                                                             &
                              &              )                                                                           &
                              & ) then
                            massFunctions(currentAnalysis)%outputNumber=k
                            exit
                         end if
                      end do
                      if (massFunctions(currentAnalysis)%outputNumber < 0) call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Functions','unable to find required redshift in outputs for mass function '//trim(massFunctionLabels(j)))
                      ! Read the appropriate observational data definition.
                      select case (trim(massFunctionLabels(j)))
                      case ('sdssStellarMassFunctionZ0.07')
                         ! SDSS z=0.07 stellar mass function.
                         !$omp critical(HDF5_Access)
                         call dataFile%openFile(char(Galacticus_Input_Path())//'/data/observations/massFunctionsStellar/Stellar_Mass_Function_Li_White_2009.hdf5',readOnly=.true.)
                         call dataFile   %readDataset  ('mass'          ,massFunctions(currentAnalysis)%masses             )
                         massDataset=dataFile%openDataset('mass'      )
                         call massDataset%readAttribute('hubbleExponent',massHubbleExponent                                )
                         call massDataset%close()
                         parameters =dataFile%openGroup  ('Parameters')
                         call parameters %readAttribute('H_0'           ,dataHubbleParameter                               )
                         call parameters %close()
                         call dataFile   %close()
                         !$omp end critical(HDF5_Access)
                         ! Adjust masses to current Hubble parameter.
                         hubbleRatio=H_0()/dataHubbleParameter
                         massFunctions(currentAnalysis)%masses=massFunctions(currentAnalysis)%masses*hubbleRatio**massHubbleExponent
                         ! Construct mass function array.
                         massFunctions(currentAnalysis)%massesCount=size(massFunctions(currentAnalysis)%masses)
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmic       ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMinimum,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMaximum,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunction            ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunctionCovariance  ,[massFunctions(currentAnalysis)%massesCount,massFunctions(currentAnalysis)%massesCount])
                         massFunctions(currentAnalysis)%massesLogarithmic     =log10(massFunctions(currentAnalysis)%masses)
                         massFunctions(currentAnalysis)%massFunction          =0.0d0
                         massFunctions(currentAnalysis)%massFunctionCovariance=0.0d0
                         do k=1,massFunctions(currentAnalysis)%massesCount
                            if (k ==                                          1) then
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)-0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)-massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k-1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                            if (k == massFunctions(currentAnalysis)%massesCount) then
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)+0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k  )-massFunctions(currentAnalysis)%massesLogarithmic(k-1))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                         end do
                      case ('alfalfaHiMassFunctionZ0.00')
                         !$omp critical (FoX_DOM_Access)
                         doc => parseFile(char(Galacticus_Input_Path())//"data/observations/massFunctionsHI/HI_Mass_Function_ALFALFA_2010.xml",iostat=ioErr)
                         if (ioErr /= 0) call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Functions','Unable to find data file')
                         massFunctionElement   => item(getElementsByTagname(doc                ,"massFunctions" ),0)
                         columnElement         => item(getElementsByTagname(massFunctionElement,"columns"       ),0)
                         massElement           => item(getElementsByTagname(      columnElement,"mass"          ),0)
                         hubbleElement         => item(getElementsByTagname(      columnElement,"hubble"        ),0)
                         hubbleExponentElement => item(getElementsByTagname(      columnElement,"hubbleExponent"),0)
                         dataList              =>      getElementsByTagname(        massElement,"data"          ) 
                         call extractDataContent(hubbleExponentElement,massHubbleExponent )
                         call extractDataContent(hubbleElement        ,dataHubbleParameter)
                         ! Construct mass function array.
                         massFunctions(currentAnalysis)%massesCount=getLength(dataList)
                         call Alloc_Array(massFunctions(currentAnalysis)%masses                  ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmic       ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMinimum,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massesLogarithmicMaximum,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunction            ,[massFunctions(currentAnalysis)%massesCount                                           ])
                         call Alloc_Array(massFunctions(currentAnalysis)%massFunctionCovariance  ,[massFunctions(currentAnalysis)%massesCount,massFunctions(currentAnalysis)%massesCount])
                         do iDatum=0,getLength(dataList)-1
                            datum => item(dataList,iDatum)
                            call extractDataContent(datum,massFunctions(currentAnalysis)%masses(iDatum+1))
                         end do
                         ! Destroy the document.
                         call destroy(doc)
                         !$omp end critical (FoX_DOM_Access)
                         ! Adjust masses to current Hubble parameter.
                         hubbleRatio                          =H_0()/dataHubbleParameter
                         massFunctions(currentAnalysis)%masses=massFunctions(currentAnalysis)%masses*hubbleRatio**massHubbleExponent
                         massFunctions(currentAnalysis)%massesLogarithmic     =log10(massFunctions(currentAnalysis)%masses)
                         massFunctions(currentAnalysis)%massFunction          =0.0d0
                         massFunctions(currentAnalysis)%massFunctionCovariance=0.0d0
                         do k=1,massFunctions(currentAnalysis)%massesCount
                            if (k ==                                          1) then
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)-0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)-massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMinimum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k-1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                            if (k == massFunctions(currentAnalysis)%massesCount) then
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=massFunctions(currentAnalysis)%massesLogarithmic(k)+0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k  )-massFunctions(currentAnalysis)%massesLogarithmic(k-1))
                            else
                               massFunctions(currentAnalysis)%massesLogarithmicMaximum(k)=                                                   +0.5d0*(massFunctions(currentAnalysis)%massesLogarithmic(k+1)+massFunctions(currentAnalysis)%massesLogarithmic(k  ))
                            end if
                         end do
                      case ('primusStellarMassFunctionZ0.100')
                         call Load_PRIMUS_Mass_Function(0,massFunctions(currentAnalysis))
                      case ('primusStellarMassFunctionZ0.250')
                         call Load_PRIMUS_Mass_Function(1,massFunctions(currentAnalysis))
                      case ('primusStellarMassFunctionZ0.350')
                         call Load_PRIMUS_Mass_Function(2,massFunctions(currentAnalysis))
                      case ('primusStellarMassFunctionZ0.450')
                         call Load_PRIMUS_Mass_Function(3,massFunctions(currentAnalysis))
                      case ('primusStellarMassFunctionZ0.575')
                         call Load_PRIMUS_Mass_Function(4,massFunctions(currentAnalysis))
                      case ('primusStellarMassFunctionZ0.725')
                         call Load_PRIMUS_Mass_Function(5,massFunctions(currentAnalysis))
                      case ('primusStellarMassFunctionZ0.900')
                         call Load_PRIMUS_Mass_Function(6,massFunctions(currentAnalysis))
                      case default
                         call Galacticus_Error_Report('Galacticus_Output_Analysis_Mass_Functions','unknown mass function')
                      end select
                      exit
                   end if
                end do
             end do
          end if
          ! Record that module is initialized.
          moduleInitialized=.true.
       end if
       !$omp end critical(Galacticus_Output_Analysis_Mass_Functions_Initialize)
    end if
    ! Return if this analysis is not active.
    if (.not.analysisActive) return
    ! Allocate work arrays.
    if (.not.allocated(thisGalaxy)) allocate(thisGalaxy(size(massFunctions)))
    ! Iterate over active analyses.
    do i=1,size(massFunctions)
       ! Return if this analysis is not active, or if this is not the correct output.
       if (iOutput /= massFunctions(i)%outputNumber) cycle
       ! Allocate workspace.
       if (.not.allocated(thisGalaxy(i)%massFunction)) call Alloc_Array(thisGalaxy(i)%massFunction,[massFunctions(i)%massesCount                             ])
       if (.not.allocated(thisGalaxy(i)%covariance  )) call Alloc_Array(thisGalaxy(i)%covariance  ,[massFunctions(i)%massesCount,massFunctions(i)%massesCount])
       ! Get the galactic stellar mass.
       mass=                                                                                                                &
            &  Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeDisk    ,massType=massFunctions(i)%descriptor%massType) &
            & +Galactic_Structure_Enclosed_Mass(thisNode,radiusLarge,componentType=componentTypeSpheroid,massType=massFunctions(i)%descriptor%massType)
       if (mass            <=                  0.0d0) return
       if (associated(massFunctions(i)%descriptor%mapMass)) mass=massFunctions(i)%descriptor%mapMass(mass,thisNode)
       massLogarithmic=log10(mass)
       do j=1,massFunctions(i)%descriptor%systematicCoefficientCount
          massLogarithmic=massLogarithmic+massFunctions(i)%systematicCoefficients(j)*(log10(mass)-massFunctions(i)%descriptor%systematicLogM0)**(j-1)
       end do
       if (massLogarithmic <  massFunctions(i)%descriptor%massLogarithmicMinimum) return
       ! Compute contributions to each bin.
       thisGalaxy(i)%massFunction=(                                                                                                                      &
            &                      +erf((massFunctions(i)%massesLogarithmicMaximum-massLogarithmic)/massFunctions(i)%descriptor%randomError/sqrt(2.0d0)) &
            &                      -erf((massFunctions(i)%massesLogarithmicMinimum-massLogarithmic)/massFunctions(i)%descriptor%randomError/sqrt(2.0d0)) &
            &                     )                                                                                                                      &
            &                     /2.0d0                                                                                                                 &
            &                     *thisTree%volumeWeight
       forall(j=1:massFunctions(i)%massesCount)
          forall(k=j:massFunctions(i)%massesCount)
             thisGalaxy(i)%covariance(j,k)=thisGalaxy(i)%massFunction(j)*thisGalaxy(i)%massFunction(k)
             thisGalaxy(i)%covariance(k,j)=thisGalaxy(i)%covariance(j,k)
          end forall
       end forall
       ! Accumulate mass function and covariance.
       !$omp critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
       massFunctions(i)%massFunction          =massFunctions(i)%massFunction          +thisGalaxy(i)%massFunction
       massFunctions(i)%massFunctionCovariance=massFunctions(i)%massFunctionCovariance+thisGalaxy(i)%covariance
       !$omp end critical (Galacticus_Output_Analysis_Mass_Functions_Accumulate)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Functions

  !# <hdfPreCloseTask>
  !#  <unitName>Galacticus_Output_Analysis_Mass_Functions_Output</unitName>
  !# </hdfPreCloseTask>
  subroutine Galacticus_Output_Analysis_Mass_Functions_Output
    !% Outputs SDSS $z\approx 0.07$ stellar mass function to file.
    use Galacticus_HDF5
    use Numerical_Constants_Astronomical
    implicit none
    integer             :: i,j,k
    type   (hdf5Object) :: analysisGroup,massFunctionGroup,thisDataset

    ! Return immediately if this analysis is not active.
    if (.not.analysisActive) return
    ! Iterate over mass functions.
    do k=1,size(massFunctions)
       ! Convert model mass function to differential per log(M).
       forall(i=1:massFunctions(k)%massesCount)
          massFunctions(k)%massFunction             (i  )=  massFunctions(k)%massFunction          (i  )                              &
               &                         /(massFunctions(k)%massesLogarithmicMaximum(i)-massFunctions(k)%massesLogarithmicMinimum(i)) &
               &                         / log(10.0d0)
          forall(j=1:massFunctions(k)%massesCount)
             massFunctions(k)%massFunctionCovariance(i,j)=  massFunctions(k)%massFunctionCovariance(i,j)                              &
                  &                      /(massFunctions(k)%massesLogarithmicMaximum(i)-massFunctions(k)%massesLogarithmicMinimum(i)) &
                  &                      /(massFunctions(k)%massesLogarithmicMaximum(j)-massFunctions(k)%massesLogarithmicMinimum(j)) &
                  &                      / log(10.0d0)**2
          end forall
       end forall
       ! Output the mass function.
       !$omp critical(HDF5_Access)
       analysisGroup    =galacticusOutputFile%openGroup('analysis','Model analysis')
       massFunctionGroup=analysisGroup       %openGroup(trim(massFunctions(k)%descriptor%label),trim(massFunctions(k)%descriptor%comment))
       call massFunctionGroup%writeDataset  (massFunctions(k)%masses                ,'mass'                  ,'Mass'             ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(massSolar             ,'unitsInSI'                                                                            )
       call thisDataset      %close()
       call massFunctionGroup%writeDataset  (massFunctions(k)%massFunction          ,'massFunction'          ,'Mass function'           ,datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(1.0d0/megaParsec**3   ,'unitsInSI'                                                                            )
       call thisDataset      %close()
       call massFunctionGroup%writeDataset  (massFunctions(k)%massFunctionCovariance,'massFunctionCovariance','Mass function covariance',datasetReturned=thisDataset)
       call thisDataset      %writeAttribute(1.0d0/megaParsec**6   ,'unitsInSI'                                                                            )
       call thisDataset      %close()
       call massFunctionGroup%close()
       call analysisGroup    %close()
       !$omp end critical(HDF5_Access)
    end do
    return
  end subroutine Galacticus_Output_Analysis_Mass_Functions_Output

  double precision function Map_Mass_ALFALFA_HI_Mass_Function_Z0_00(mass,thisNode)
    !% Maps gas masses into HI masses for the ALFALFA survey analysis. Assumes a constant gas to HI mass conversion factor of
    !% $0.54$ \citep{power_redshift_2010}.
    double precision          , intent(in   )          :: mass
    type            (treeNode), intent(inout), pointer :: thisNode
    ! Constant gas to HI mass conversion factor.
    double precision          , parameter              :: hiFractionByMass=0.54d0

    Map_Mass_ALFALFA_HI_Mass_Function_Z0_00=hiFractionByMass*mass
    return
  end function Map_Mass_ALFALFA_HI_Mass_Function_Z0_00

  subroutine Load_PRIMUS_Mass_Function(massFunctionIndex,thisMassFunction)
    !% Load the specified mass function from the PRIMUS stellar mass function dataset.
    use ISO_Varying_String
    use Galacticus_Error
    use Cosmological_Parameters
    use Galacticus_Input_Paths
    use Memory_Management
    use FoX_dom
    implicit none
    integer                         , intent(in   ) :: massFunctionIndex
    type            (massFunction  ), intent(inout) :: thisMassFunction
    type            (node          ), pointer       :: doc,massFunctionElement,columnElement,massElement,hubbleElement,hubbleExponentElement,datum
    type            (nodeList      ), pointer       :: dataList
    integer                                         :: k,iDatum,ioErr
    double precision                                :: hubbleRatio,dataHubbleParameter,massHubbleExponent

    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(Galacticus_Input_Path())//"data/observations/massFunctionsStellar/Stellar_Mass_Function_PRIMUS_2013.xml",iostat=ioErr)
    if (ioErr /= 0) call Galacticus_Error_Report('Load_PRIMUS_Mass_Function','Unable to find data file')
    massFunctionElement   => item(getElementsByTagname(doc                ,"stellarMassFunction" ),massFunctionIndex)
    columnElement         => item(getElementsByTagname(massFunctionElement,"columns"             ),                0)
    massElement           => item(getElementsByTagname(      columnElement,"stellarMass"         ),                0)
    hubbleElement         => item(getElementsByTagname(      columnElement,"hubble"              ),                0)
    hubbleExponentElement => item(getElementsByTagname(      columnElement,"hubbleExponent"      ),                0)
    dataList              =>      getElementsByTagname(        massElement,"datum"               ) 
    call extractDataContent(hubbleExponentElement,massHubbleExponent )
    call extractDataContent(hubbleElement        ,dataHubbleParameter)
    ! Construct mass function array.
    thisMassFunction%massesCount=getLength(dataList)
    call Alloc_Array(thisMassFunction%masses                  ,[thisMassFunction%massesCount                             ])
    call Alloc_Array(thisMassFunction%massesLogarithmic       ,[thisMassFunction%massesCount                             ])
    call Alloc_Array(thisMassFunction%massesLogarithmicMinimum,[thisMassFunction%massesCount                             ])
    call Alloc_Array(thisMassFunction%massesLogarithmicMaximum,[thisMassFunction%massesCount                             ])
    call Alloc_Array(thisMassFunction%massFunction            ,[thisMassFunction%massesCount                             ])
    call Alloc_Array(thisMassFunction%massFunctionCovariance  ,[thisMassFunction%massesCount,thisMassFunction%massesCount])
    do iDatum=0,getLength(dataList)-1
       datum => item(dataList,iDatum)
       call extractDataContent(datum,thisMassFunction%masses(iDatum+1))
    end do
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    ! Adjust masses to current Hubble parameter.
    hubbleRatio                        =H_0()/dataHubbleParameter
    thisMassFunction%masses                =(10.0d0**thisMassFunction%masses)*(hubbleRatio**massHubbleExponent)
    thisMassFunction%massesLogarithmic     =log10(thisMassFunction%masses)
    thisMassFunction%massFunction          =0.0d0
    thisMassFunction%massFunctionCovariance=0.0d0
    do k=1,thisMassFunction%massesCount
       if (k ==                            1) then
          thisMassFunction%massesLogarithmicMinimum(k)=thisMassFunction%massesLogarithmic(k)-0.5d0*(thisMassFunction%massesLogarithmic(k+1)-thisMassFunction%massesLogarithmic(k  ))
       else
          thisMassFunction%massesLogarithmicMinimum(k)=                                 +0.5d0*(thisMassFunction%massesLogarithmic(k-1)+thisMassFunction%massesLogarithmic(k  ))
       end if
       if (k == thisMassFunction%massesCount) then
          thisMassFunction%massesLogarithmicMaximum(k)=thisMassFunction%massesLogarithmic(k)+0.5d0*(thisMassFunction%massesLogarithmic(k  )-thisMassFunction%massesLogarithmic(k-1))
       else
          thisMassFunction%massesLogarithmicMaximum(k)=                                 +0.5d0*(thisMassFunction%massesLogarithmic(k+1)+thisMassFunction%massesLogarithmic(k  ))
       end if
    end do
    return
  end subroutine Load_PRIMUS_Mass_Function

end module Galacticus_Output_Analyses_Mass_Functions
