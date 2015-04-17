!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo profile concentrations using the \cite{klypin_multidark_2014} algorithm.
  
  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationKlypin2015">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{klypin_multidark_2014}.</description>
  !# </darkMatterProfileConcentration>

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationKlypin2015
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{klypin_multidark_2014}.
     private
     integer                 :: virialDensityContrast, fittingFunction
     type   (table1DGeneric) :: fitParameters
   contains
     procedure :: concentration               => klypin2015Concentration
     procedure :: densityContrastDefinition   => klypin2015DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => klypin2015DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationKlypin2015
  
  interface darkMatterProfileConcentrationKlypin2015
     !% Constructors for the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.
     module procedure klypin2015DefaultConstructor
     module procedure klypin2015Constructor
  end interface darkMatterProfileConcentrationKlypin2015

  ! Labels for virial density contrast definition.
  integer, parameter :: klypin2015VirialDensityContrastFixed       = 0
  integer, parameter :: klypin2015VirialDensityContrastVirial      = 1

  ! Labels for fitting function type.
  integer, parameter :: klypin2015FittingFunctionEqn24             = 0
  integer, parameter :: klypin2015FittingFunctionEqn25             = 1

  ! Labels for sample selection.
  integer, parameter :: klypin2015SamplePlanck200CritRelaxedMass   = 0
  integer, parameter :: klypin2015SamplePlanck200CritAllMass       = 1
  integer, parameter :: klypin2015SamplePlanck200CritRelaxedVmax   = 2
  integer, parameter :: klypin2015SamplePlanck200CritAllVmax       = 3
  integer, parameter :: klypin2015SamplePlanckVirialRelaxedMass    = 4
  integer, parameter :: klypin2015SamplePlanckVirialAllMass        = 5
  integer, parameter :: klypin2015SamplePlanckVirialRelaxedVmax    = 6
  integer, parameter :: klypin2015SamplePlanckVirialAllVmax        = 7
  integer, parameter :: klypin2015SampleWMAP7200CritRelaxedMass    = 8
  integer, parameter :: klypin2015SampleWMAP7200CritAllMass        = 9
  integer, parameter :: klypin2015SampleWMAP7200CritRelaxedVmax    =10
  integer, parameter :: klypin2015SampleWMAP7VirialRelaxedMass     =11
  integer, parameter :: klypin2015SampleWMAP7VirialAllMass         =12
  integer, parameter :: klypin2015SampleWMAP7VirialRelaxedVmax     =13
  integer, parameter :: klypin2015SamplePlanck200CritAllMassUni    =14
  integer, parameter :: klypin2015SamplePlanck200CritRelaxedMassUni=15
  integer, parameter :: klypin2015SamplePlanckVirialAllMassUni     =16
  integer, parameter :: klypin2015SamplePlanckVirialRelaxedMassUni =17

  ! Default sample.
  integer            :: klypin2015ConcentrationSample

  ! Initialization status.
  logical            :: klypin2015Initialized                      =.false.
  
contains
  
  function klypin2015DefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.
    use Galacticus_Error
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationKlypin2015), target  :: klypin2015DefaultConstructor
    type(varying_string                          )          :: klypin2015ConcentrationSampleText
    
    if (.not.klypin2015Initialized) then
       !$omp critical(klypin2015DefaultInitialize)
       if (.not.klypin2015Initialized) then
          ! Get parameters of the model.
          !@ <inputParameter>
          !@   <name>klypin2015ConcentrationSample</name>
          !@   <defaultValue>planck200Crit</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The sample to use for the halo concentration algorithm of \cite{klypin_multidark_2014}.
          !@   </description>
          !@   <type>string</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("klypin2015ConcentrationSample",klypin2015ConcentrationSampleText,defaultValue='planck200Crit')
          select case (char(klypin2015ConcentrationSampleText))
          case ('planck200CritRelaxedMass'         )
             klypin2015ConcentrationSample=klypin2015SamplePlanck200CritRelaxedMass
          case ('planck200CritAllMass'             )
             klypin2015ConcentrationSample=klypin2015SamplePlanck200CritAllMass
          case ('planck200CritRelaxedVmax'         )
             klypin2015ConcentrationSample=klypin2015SamplePlanck200CritRelaxedVmax
          case ('planck200CritAllVmax'             )
             klypin2015ConcentrationSample=klypin2015SamplePlanck200CritAllVmax
          case ('planckVirialRelaxedMass'          )
             klypin2015ConcentrationSample=klypin2015SamplePlanckVirialRelaxedMass
          case ('planckVirialAllMass'              )
             klypin2015ConcentrationSample=klypin2015SamplePlanckVirialAllMass
          case ('planckVirialRelaxedVmax'          )
             klypin2015ConcentrationSample=klypin2015SamplePlanckVirialRelaxedVmax
          case ('planckVirialAllVmax'              )
             klypin2015ConcentrationSample=klypin2015SamplePlanckVirialAllVmax
          case ('wmap7200CritRelaxedMass'          )
             klypin2015ConcentrationSample=klypin2015SampleWMAP7200CritRelaxedMass
          case ('wmap7200CritAllMass'              )
             klypin2015ConcentrationSample=klypin2015SampleWMAP7200CritAllMass
          case ('wmap7200CritRelaxedVmax'          )
             klypin2015ConcentrationSample=klypin2015SampleWMAP7200CritRelaxedVmax
          case ('wmap7VirialRelaxedMass'           )
             klypin2015ConcentrationSample=klypin2015SampleWMAP7VirialRelaxedMass
          case ('wmap7VirialAllMass'               )
             klypin2015ConcentrationSample=klypin2015SampleWMAP7VirialAllMass
          case ('wmap7VirialRelaxedVmax'           )
             klypin2015ConcentrationSample=klypin2015SampleWMAP7VirialRelaxedVmax
          case ('planck200CritAllMassUniversal'    )
             klypin2015ConcentrationSample=klypin2015SamplePlanck200CritAllMassUni
          case ('planck200CritRelaxedMassUniversal')
             klypin2015ConcentrationSample=klypin2015SamplePlanck200CritRelaxedMassUni
          case ('planckVirialAllVmaxUniversal'     )
             klypin2015ConcentrationSample=klypin2015SamplePlanckVirialAllMassUni
          case ('planckVirialRelaxedVmaxUniversal' )
             klypin2015ConcentrationSample=klypin2015SamplePlanckVirialRelaxedMassUni
          case default
             call Galacticus_Error_Report('klypin2015DefaultConstructor','unrecognized sample')
          end select
          ! Record that method is now initialized.
          klypin2015Initialized=.true.
       end if
       !$omp end critical(klypin2015DefaultInitialize)
    end if
    ! Construct the object.
    klypin2015DefaultConstructor=klypin2015Constructor(klypin2015ConcentrationSample)
    return
  end function klypin2015DefaultConstructor
  
  function klypin2015Constructor(sample)
    !% Constructor for the {\normalfont \ttfamily klypin2015} dark matter halo profile concentration class.
    implicit none
    type   (darkMatterProfileConcentrationKlypin2015)                :: klypin2015Constructor
    integer                                          , intent(in   ) :: sample
    
    select case (sample)
    case (klypin2015SamplePlanck200CritRelaxedMass)
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.750d0,6.700d0,6.250d0,5.020d0,4.190d0,3.300d0,3.000d0,2.720d0,2.400d0,2.100d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.000d-1,0.950d-1,0.920d-1,0.880d-1,0.850d-1,0.830d-1,0.800d-1,0.800d-1,0.800d-1,0.800d-1], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[4.500d5,2.000d4,8.000d3,7.800d2,1.600d2,2.700d1,1.400d1,6.800d0,1.600d0,0.300d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritAllMass    )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.400d0,6.250d0,5.650d0,4.300d0,3.530d0,2.700d0,2.420d0,2.200d0,1.920d0,1.650d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.120d0,0.117d0,0.115d0,0.110d0,0.095d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.500d5,1.000d5,2.000d4,9.000d2,3.000d2,4.200d1,1.700d1,8.500d0,2.000d0,0.300d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritRelaxedVmax)
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[8.000d0,6.820d0,6.400d0,5.200d0,4.350d0,3.500d0,3.120d0,2.850d0,2.550d0,2.160d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.000d-1,0.950d-1,0.920d-1,0.880d-1,0.850d-1,0.800d-1,0.800d-1,0.800d-1,0.800d-1,0.800d-1], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d5,9.000d3,4.500d3,6.000d2,1.500d2,2.700d1,1.100d1,5.500d0,1.500d0,0.220d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritAllVmax    )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.750d0,6.500d0,5.950d0,4.550d0,3.680d0,2.750d0,2.500d0,2.250d0,2.050d0,1.760d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.115d0,0.115d0,0.115d0,0.110d0,0.105d0,1.000d-1,0.950d-1,0.900d-1,0.800d-1,0.800d-1], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.500d5,1.800d4,6.000d3,6.000d2,1.500d2,2.000d1,1.000d1,5.000d0,1.500d0,0.250d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialRelaxedMass )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.020d1,7.850d0,7.160d0,5.450d0,4.550d0,3.550d0,3.240d0,2.920d0,2.600d0,2.300d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.100d0,0.095d0,0.092d0,0.088d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.000d5,1.200d4,5.500d3,7.000d2,1.800d2,3.000d1,1.500d1,7.000d0,1.900d0,0.360d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialAllMass     )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.750d0,7.250d0,6.500d0,4.750d0,3.800d0,3.000d0,2.650d0,2.420d0,2.100d0,1.860d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.110d0,0.107d0,0.105d0,0.100d0,0.095d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.000d5,2.200d4,1.000d4,1.000d3,2.100d2,4.300d1,1.800d1,9.000d0,1.900d0,0.420d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialRelaxedVmax )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.070d1,8.100d0,7.330d0,5.650d0,4.650d0,3.700d0,2.350d0,2.980d0,2.700d0,2.350d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.110d0,0.100d0,0.100d0,0.088d0,0.085d0,0.080d0,0.080d0,0.080d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.400d4,5.000d3,2.200d3,5.200d2,1.200d2,2.500d1,1.200d1,5.000d0,1.400d0,0.260d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialAllVmax     )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.350d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0,5.400d0], &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.030d1,7.600d0,6.830d0,4.960d0,3.960d0,3.000d0,2.730d0,2.450d0,2.240d0,2.030d0], &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.115d0,0.115d0,0.115d0,0.110d0,0.105d0,0.100d0,0.095d0,0.090d0,0.080d0,0.080d0], &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[4.800d4,5.000d3,2.700d3,3.900d2,1.100d2,1.800d1,1.000d1,5.000d0,1.400d0,0.360d0], &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7200CritRelaxedMass  )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[6.900d0,5.700d0,4.550d0,3.750d0,2.900d0,2.600d0,2.400d0,2.200d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.090d0,0.088d0,0.086d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[5.500d5,6.000d3,5.000d2,1.000d2,2.000d1,1.000d1,6.000d0,3.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7200CritAllMass      )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[6.600d0,5.250d0,3.850d0,3.000d0,2.100d0,1.800d0,1.600d0,1.400d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.110d0,0.105d0,0.103d0,0.097d0,0.095d0,0.095d0,0.095d0,0.095d0]                , &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d6,6.000d4,8.000d2,1.100d2,1.300d1,6.000d0,3.000d0,1.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7200CritRelaxedVmax  )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[7.200d0,5.900d0,4.700d0,3.850d0,3.000d0,2.700d0,2.500d0,2.300d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.090d0,0.088d0,0.086d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d5,4.000d3,4.000d2,8.000d1,1.300d1,7.000d0,3.500d0,2.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7VirialRelaxedMass   )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.500d0,6.750d0,5.000d0,4.050d0,3.100d0,2.800d0,2.450d0,2.200d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.090d0,0.088d0,0.086d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[3.000d5,5.000d3,4.500d2,9.000d1,1.500d1,8.000d0,3.500d0,1.500d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7VirialAllMass       )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.000d0,6.000d0,4.300d0,3.300d0,2.300d0,2.100d0,1.850d0,1.700d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.100d0,0.100d0,0.100d0,0.100d0,0.095d0,0.095d0,0.095d0,0.095d0]                , &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[2.000d6,7.000d3,5.500d2,9.000d1,1.100d1,6.000d0,2.500d0,1.000d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SampleWMAP7VirialRelaxedVmax    )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn24
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.500d0,1.000d0,1.440d0,2.150d0,2.500d0,2.900d0,4.100d0]                , &
            &            tableCount       =3                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[9.750d0,7.020d0,5.230d0,4.250d0,3.200d0,2.900d0,2.500d0,2.350d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.085d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0,0.085d0]                , &
            &            table            =2                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[1.300d5,4.000d3,4.000d2,8.000d1,1.100d1,6.000d0,2.500d0,1.200d0]                , &
            &            table            =3                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritAllMassUni    )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn25
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,2.890d0,5.410d0], &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.400d0,0.650d0,0.820d0,1.080d0,1.230d0,1.600d0,1.680d0,1.700d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.278d0,0.375d0,0.411d0,0.436d0,0.426d0,0.375d0,0.360d0,0.351d0]                , &
            &            table            =2                                                                                  &
            &           )
    case (klypin2015SamplePlanck200CritRelaxedMassUni)
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastFixed
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn25
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,2.890d0,5.410d0]                , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.950d0,1.060d0,1.150d0,1.280d0,1.390d0,1.660d0,1.700d0,1.720d0]                , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.522d0,0.550d0,0.562d0,0.562d0,0.541d0,0.480d0,0.464d0,0.450d0]                , &
            &            table            =2                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialAllMassUni     )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn25
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,5.500d0]                        , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.750d0,0.900d0,0.970d0,1.120d0,1.280d0,1.520d0,1.620d0]                        , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.567d0,0.541d0,0.529d0,0.496d0,0.474d0,0.421d0,0.393d0]                        , &
            &            table            =2                                                                                  &
            &           )
    case (klypin2015SamplePlanckVirialRelaxedMassUni )
       klypin2015Constructor%virialDensityContrast=klypin2015VirialDensityContrastVirial
       klypin2015Constructor%fittingFunction      =klypin2015FittingFunctionEqn25
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %create  (                                                                                                     &
            &            x                =[0.000d0,0.380d0,0.500d0,1.000d0,1.440d0,2.500d0,5.500d0]                        , &
            &            tableCount       =2                                                                                , &
            &            extrapolationType=extrapolationTypeFix                                                               &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.990d0,1.100d0,1.160d0,1.290d0,1.410d0,1.650d0,1.720d0]                        , &
            &            table            =1                                                                                  &
            &           )
       call klypin2015Constructor                                                                                             &
            & %fitParameters                                                                                                  &
            &  %populate(                                                                                                     &
            &            y                =[0.716d0,0.673d0,0.660d0,0.615d0,0.585d0,0.518d0,0.476d0]                        , &
            &            table            =2                                                                                  &
            &           )
    end select
    return
  end function klypin2015Constructor
  
  double precision function klypin2015Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the \cite{klypin_multidark_2014} algorithm.
    use Power_Spectra
    use Cosmology_Parameters
    use Cosmology_Functions
    implicit none
    class           (darkMatterProfileConcentrationKlypin2015), intent(inout)          :: self
    type            (treeNode                                ), intent(inout), pointer :: node
    class           (nodeComponentBasic                      )               , pointer :: basic
    class           (cosmologyParametersClass                )               , pointer :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 )               , pointer :: cosmologyFunctions_
    double precision                                                                   :: massLittleH         , concentration0, &
         &                                                                                mass0               , gamma         , &
         &                                                                                redshift            , a0            , &
         &                                                                                b0                  , sigma

    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    cosmologyFunctions_  => cosmologyFunctions ()
    ! Get the basic component, and find the halo mass in "little-h" units.
    basic       =>  node                %basic         (                  )
    massLittleH =  +basic               %mass          (                  ) &
         &         *cosmologyParameters_%HubbleConstant(units=hubbleUnitsLittleH)
    ! Find redshift.
    redshift=cosmologyFunctions_ %redshiftFromExpansionFactor(    &
         &    cosmologyFunctions_%expansionFactor             (   &
         &     basic             %time                         () &
         &                                                    )   &
         &                                                   )
    ! Determine which fitting function to use.
    select case (self%fittingFunction)
    case (klypin2015FittingFunctionEqn24)
       ! Evaluate fitting function parameters.
       concentration0=self%fitParameters%interpolate(redshift,table=1)
       mass0         =self%fitParameters%interpolate(redshift,table=2)
       gamma         =self%fitParameters%interpolate(redshift,table=3)
       ! Evaluate the concentration.
       klypin2015Concentration=+concentration0                &
            &                  *(                             &
            &                    +1.0d0                       &
            &                    +(massLittleH/mass0 )**0.4d0 &
            &                  )                              &
            &                  /  (massLittleH/1.0d12)**gamma
    case (klypin2015FittingFunctionEqn25)
       ! Find sigma.
       sigma         =Cosmological_Mass_Root_Variance(basic%mass())
       ! Evaluate fitting function parameters.
       a0            =self%fitParameters%interpolate(redshift,table=1)
       b0            =self%fitParameters%interpolate(redshift,table=2)
       ! Evaluate the concentration.
       klypin2015Concentration=+b0          &
            &                  *(           &
            &                    +1.00d0    &
            &                    +7.37d0    &
            &                    *(         &
            &                      +sigma   &
            &                      /a0      &
            &                     )**0.75d0 &
            &                   )           &
            &                  *(           &
            &                    +1.00d0    &
            &                    +0.14d0    &
            &                    /(         &
            &                      +sigma   &
            &                      /a0      &
            &                     )**2      &
            &                   )
    end select
    return
  end function klypin2015Concentration
  
  function klypin2015DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the \cite{klypin_multidark_2014} algorithm.
    implicit none
    class(virialDensityContrastClass              ), pointer       :: klypin2015DensityContrastDefinition
    class(darkMatterProfileConcentrationKlypin2015), intent(inout) :: self

    select case (self%virialDensityContrast)
    case (klypin2015VirialDensityContrastFixed)
       allocate(virialDensityContrastFixed :: klypin2015DensityContrastDefinition)
       select type (klypin2015DensityContrastDefinition)
       type is (virialDensityContrastFixed)
          klypin2015DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
       end select
    case (klypin2015VirialDensityContrastVirial)
       allocate(virialDensityContrastSphericalCollapseMatterLambda :: klypin2015DensityContrastDefinition)
       select type (klypin2015DensityContrastDefinition)
       type is (virialDensityContrastSphericalCollapseMatterLambda)
          klypin2015DensityContrastDefinition=virialDensityContrastSphericalCollapseMatterLambda()
       end select
    end select
    return
  end function klypin2015DensityContrastDefinition

  function klypin2015DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{klypin_multidark_2014} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: klypin2015DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationKlypin2015           ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: klypin2015DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition       )
    select type (klypin2015DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition       =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          klypin2015DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function klypin2015DarkMatterProfileDefinition
