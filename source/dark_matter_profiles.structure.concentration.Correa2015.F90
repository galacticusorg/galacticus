!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

  !% An implementation of dark matter halo profile concentrations using the \cite{correa_accretion_2015} algorithm.

  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationCorrea2015">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{correa_accretion_2015}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationCorrea2015
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{correa_accretion_2015}.
     private
     double precision :: A
   contains
     final     ::                                correa2015Destructor
     procedure :: concentration               => correa2015Concentration
     procedure :: densityContrastDefinition   => correa2015DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => correa2015DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationCorrea2015

  interface darkMatterProfileConcentrationCorrea2015
     !% Constructors for the {\normalfont \ttfamily correa2015} dark matter halo profile concentration class.
     module procedure correa2015ConstructorParameters
     module procedure correa2015ConstructorInternal
  end interface darkMatterProfileConcentrationCorrea2015

contains

  function correa2015ConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily correa2015} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationCorrea2015)                :: correa2015ConstructorParameters
    type(inputParameters                         ), intent(in   ) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>A</name>
    !#   <source>parameters</source>
    !#   <variable>correa2015ConstructorParameters%A</variable>
    !#   <defaultValue>887.0d0</defaultValue>
    !#   <defaultSource>\cite{correa_accretion_2015}</defaultSource>
    !#   <description>The parameter $A$ appearin in eqn.~(17) of \cite{correa_accretion_2015}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    return
  end function correa2015ConstructorParameters
  
  function correa2015ConstructorInternal(A)
    !% Constructor for the {\normalfont \ttfamily correa2015} dark matter halo profile concentration class.
    implicit none
    type            (darkMatterProfileConcentrationCorrea2015)                :: correa2015ConstructorInternal
    double precision                                          , intent(in   ) :: A

    correa2015ConstructorInternal%A=A
    return
  end function correa2015ConstructorInternal
  
  subroutine correa2015Destructor(self)
    !% Destructor for the {\normalfont \ttfamily correa2015} dark matter halo profile concentration
    !% class.
    implicit none
    type(darkMatterProfileConcentrationCorrea2015), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine correa2015Destructor

  double precision function correa2015Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the
    !% \cite{correa_accretion_2015} algorithm.
    use Cosmology_Functions
    use Cosmology_Parameters
    use Dark_Matter_Halos_Correa2015
    use Root_Finder
    implicit none
    class           (darkMatterProfileConcentrationCorrea2015), intent(inout)          :: self
    type            (treeNode                                ), intent(inout), pointer :: node
    class           (nodeComponentBasic                      )               , pointer :: basic
    class           (cosmologyParametersClass                )               , pointer :: cosmologyParameters_
    class           (cosmologyFunctionsClass                 )               , pointer :: cosmologyFunctions_
    double precision                                          , parameter              :: toleranceRelative   =1.0d-6
    double precision                                          , parameter              :: toleranceAbsolute   =0.0d+0
    type            (rootFinder                              ), save                   :: finder
    !$omp threadprivate(finder)
    double precision                                                                   :: mass                       , time           , &
         &                                                                                aTilde                     , bTilde         , &
         &                                                                                redshift                   , expansionFactor

    ! Get properties of the node.
    basic => node %basic()
    mass  =  basic%mass ()
    time  =  basic%time ()
    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    cosmologyFunctions_  => cosmologyFunctions ()
    ! Determine the base redshift. 
    expansionFactor=cosmologyFunctions_%expansionFactor            (time           )
    redshift       =cosmologyFunctions_%redshiftFromExpansionFactor(expansionFactor)
    ! Find the a~ and b~ parameters.
    call Dark_Matter_Halo_Correa2015_Fit_Parameters(mass,expansionFactor,aTilde,bTilde)
    ! Solve for the redshift corresponding to this mass.
    if (.not.finder%isInitialized()) then
       call finder%rootFunction(concentrationSolver                  )
       call finder%tolerance   (toleranceAbsolute  ,toleranceRelative)
       call finder%rangeExpand (                                                   &
            &                   rangeExpandDownward=0.5d0                        , &
            &                   rangeExpandUpward  =2.0d0                        , &
            &                   rangeExpandType    =rangeExpandMultiplicative      &
            &                  )
    end if
    correa2015Concentration=finder%find(rootGuess=6.0d0)
    return

  contains

    double precision function concentrationSolver(concentration)
      !% Solver used in finding halo concentration using the \cite{correa_accretion_2015} algorithm.
      implicit none
      double precision, intent(in   ) :: concentration
      double precision                :: Y1           , Yc               , &
           &                             f1           , f2               , &
           &                             densityInner , redshiftFormation

      ! Evaluate left-hand side of eqn. 18 of Correa et al. (2015).      
      Y1                 =log(2.0d0              )-0.5d0
      Yc                 =log(1.0d0+concentration)-concentration/(1.0d0+concentration)
      f1                 =+Y1 &
           &              /Yc
      ! Evaluate mean inner density.
      densityInner       =+200.0d0          &
           &              *concentration**3 &
           &              *Y1               &
           &              /Yc
      ! Evaluate eqn. 17 of Correa et al. (2015), rearranged to give formation redshift.
      redshiftFormation  =+(                                                         &
           &                +(                                                       &
           &                  +(                                                     &
           &                    +1.0d0                                               &
           &                    +redshift                                            &
           &                   )                                     **3             &
           &                  +cosmologyParameters_%omegaDarkEnergy()                &
           &                  /cosmologyParameters_%omegaMatter    ()                &
           &                 )                                                       &
           &                *densityInner                                            &
           &                /self%A                                                  &
           &                -  cosmologyParameters_%omegaDarkEnergy()                &
           &                /  cosmologyParameters_%omegaMatter    ()                &
           &               )                                         **(1.0d0/3.0d0) &
           &              -1.0d0
      ! Evaluate right-hand side of eqn. 19 of Correa et al. (2015).
      f2                 =+(                        &
           &                +1.0d0                  &
           &                +redshiftFormation      &
           &                -redshift               &
           &               )**aTilde                &
           &              *exp(                     &
           &                   +(                   &
           &                     +redshiftFormation &
           &                     -redshift          &
           &                    )                   &
           &                   *bTilde              &
           &                  )
      ! Left- and right-hand sides will be equal for correct concentration.
      concentrationSolver=+f1 &
           &              -f2
      return
    end function concentrationSolver

  end function correa2015Concentration

  function correa2015DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the
    !% \cite{correa_accretion_2015} algorithm.
    implicit none
    class(virialDensityContrastClass              ), pointer       :: correa2015DensityContrastDefinition
    class(darkMatterProfileConcentrationCorrea2015), intent(inout) :: self
    
    allocate(virialDensityContrastFixed :: correa2015DensityContrastDefinition)
    select type (correa2015DensityContrastDefinition)
    type is (virialDensityContrastFixed)
       correa2015DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
    end select    
    return
  end function correa2015DensityContrastDefinition

  function correa2015DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{correa_accretion_2015} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: correa2015DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationCorrea2015          ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: correa2015DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition        )
    select type (correa2015DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition        =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          correa2015DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function correa2015DarkMatterProfileDefinition
