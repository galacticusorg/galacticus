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

  !% An implementation of dark matter halo profile concentrations using the \cite{prada_halo_2011} algorithm.
  
  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationPrada2011">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{prada_halo_2011}.</description>
  !# </darkMatterProfileConcentration>
  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationPrada2011
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{prada_halo_2011}.
     private
     double precision :: A, B, C, D, C0, C1, X0, X1, inverseSigma0, inverseSigma1, alpha, beta
   contains
     final     ::                                prada2011Destructor
     procedure :: concentration               => prada2011Concentration
     procedure :: densityContrastDefinition   => prada2011DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => prada2011DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationPrada2011
  
  interface darkMatterProfileConcentrationPrada2011
     !% Constructors for the {\normalfont \ttfamily prada2011} dark matter halo profile concentration class.
     module procedure prada2011ConstructorParameters
     module procedure prada2011ConstructorInternal
  end interface darkMatterProfileConcentrationPrada2011
  
contains
  
  function prada2011ConstructorParameters(parameters)
    !% Default constructor for the {\normalfont \ttfamily prada2011} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationPrada2011)                :: prada2011ConstructorParameters
    type(inputParameters                        ), intent(inout) :: parameters
    !# <inputParameterList label="allowedParameterNames" />

    ! Check and read parameters.
    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>A</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%A</variable>
    !#   <defaultValue>2.881d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $A$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>B</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%B</variable>
    !#   <defaultValue>1.257d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $b$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>C</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%C</variable>
    !#   <defaultValue>1.022d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $c$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>D</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%D</variable>
    !#   <defaultValue>0.060d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $d$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>C0</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%C0</variable>
    !#   <defaultValue>3.681d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $c_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>C1</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%C1</variable>
    !#   <defaultValue>5.033d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $c_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>X0</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%X0</variable>
    !#   <defaultValue>0.424d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $x_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>X1</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%X1</variable>
    !#   <defaultValue>0.526d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $x_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>inverseSigma0</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%inverseSigma0</variable>
    !#   <defaultValue>1.047d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $\sigma^{-1}_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>inverseSigma1</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%inverseSigma1</variable>
    !#   <defaultValue>1.646d0</defaultValue>
    !#   <attachedTo>module</attachedTo>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $\sigma^{-1}_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%alpha</variable>
    !#   <defaultValue>6.948d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <source>parameters</source>
    !#   <variable>prada2011ConstructorParameters%beta</variable>
    !#   <defaultValue>7.386d0</defaultValue>
    !#   <defaultSource>\cite{prada_halo_2011}</defaultSource>
    !#   <description>The parameter $\beta$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.</description>
    !#   <type>real</type>
    !#   <cardinality>1</cardinality>
    !# </inputParameter>
    return
  end function prada2011ConstructorParameters
  
  function prada2011ConstructorInternal(A,B,C,D,C0,C1,X0,X1,inverseSigma0,inverseSigma1,alpha,beta)
    !% Constructor for the {\normalfont \ttfamily prada2011} dark matter halo profile concentration class.
    implicit none
    type            (darkMatterProfileConcentrationPrada2011)                :: prada2011ConstructorInternal
    double precision                                         , intent(in   ) :: A                           , B            , &
         &                                                                      C                           , D            , &
         &                                                                      C0                          , C1           , &
         &                                                                      X0                          , inverseSigma0, &
         &                                                                      inverseSigma1               , alpha        , &
         &                                                                      beta                        , X1
    
    prada2011ConstructorInternal%A            =A
    prada2011ConstructorInternal%B            =B
    prada2011ConstructorInternal%C            =C
    prada2011ConstructorInternal%D            =D
    prada2011ConstructorInternal%C0           =C0
    prada2011ConstructorInternal%C1           =C1
    prada2011ConstructorInternal%X0           =X0
    prada2011ConstructorInternal%X1           =X1
    prada2011ConstructorInternal%inverseSigma0=inverseSigma0
    prada2011ConstructorInternal%inverseSigma1=inverseSigma1
    prada2011ConstructorInternal%alpha        =alpha
    prada2011ConstructorInternal%beta         =beta
    return
  end function prada2011ConstructorInternal
  
  subroutine prada2011Destructor(self)
    !% Destructor for the {\normalfont \ttfamily prada2011} dark matter halo profile concentration class.
    implicit none
    type(darkMatterProfileConcentrationPrada2011), intent(inout) :: self

    ! Nothing to do.
    return
  end subroutine prada2011Destructor

  double precision function prada2011Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\normalfont \ttfamily node} using the \cite{prada_halo_2011} algorithm.
    use Cosmological_Mass_Variance
    use Cosmology_Functions
    use Cosmology_Parameters
    use Linear_Growth
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout)          :: self
    type            (treeNode                               ), intent(inout), pointer :: node
    class           (nodeComponentBasic                     )               , pointer :: basic
    class           (cosmologyParametersClass               )               , pointer :: cosmologyParameters_
    class           (cosmologyFunctionsClass                )               , pointer :: cosmologyFunctions_
    class           (linearGrowthClass                      )               , pointer :: linearGrowth_
    class           (cosmologicalMassVarianceClass          )               , pointer :: cosmologicalMassVariance_
    double precision                                                                  :: massNode            , sigmaPrime, &
         &                                                                               timeNode            , x
    
    ! Get required objects.
    cosmologyFunctions_       => cosmologyFunctions      ()
    cosmologyParameters_      => cosmologyParameters     ()
    linearGrowth_             => linearGrowth            ()
    cosmologicalMassVariance_ => cosmologicalMassVariance()
    ! Compute concentration.
    basic                 => node %basic()
    massNode              =  basic%mass ()
    timeNode              =  basic%time ()
    x                     =  (cosmologyParameters_%OmegaDarkEnergy()/cosmologyParameters_%OmegaMatter())**(1.0d0/3.0d0)*cosmologyFunctions_%expansionFactor(timeNode)
    sigmaPrime            =  prada2011B1(self,x)*cosmologicalMassVariance_%rootVariance(massNode)*linearGrowth_%value(timeNode)
    prada2011Concentration=  prada2011B0(self,x)*prada2011C(self,sigmaPrime)
    return
  end function prada2011Concentration
  
  double precision function prada2011B0(self,x)
    !% The function $B_0(x)$ as defined in eqn.~(18) of \cite{prada_halo_2011}.
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x
    
    prada2011B0=prada2011cMin(self,x)/prada2011cMin(self,1.393d0)
  end function prada2011B0
  
  double precision function prada2011B1(self,x)
    !% The function $B_1(x)$ as defined in eqn.~(18) of \cite{prada_halo_2011}.
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x
    
    prada2011B1=prada2011inverseSigmaMin(self,x)/prada2011inverseSigmaMin(self,1.393d0)
  end function prada2011B1
  
  double precision function prada2011cMin(self,x)
    !% The function $c_{\mathrm min}(x)$ as defined in eqn.~(19) of \cite{prada_halo_2011}.
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x
    
    prada2011cMin=self%C0+(self%C1-self%C0)*(atan(self%alpha*(x-self%X0))/Pi+0.5d0)
  end function prada2011cMin
  
  double precision function prada2011inverseSigmaMin(self,x)
    !% The function $\sigma^{-1}_{\mathrm min}(x)$ as defined in eqn.~(20) of \cite{prada_halo_2011}.
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x
    
    prada2011inverseSigmaMin=self%inverseSigma0+(self%inverseSigma1-self%inverseSigma0)*(atan(self%beta*(x-self%X1))/Pi+0.5d0)
  end function prada2011inverseSigmaMin
  
  double precision function prada2011C(self,sigmaPrime)
    !% The function $\mathcal{C}(\sigma^\prime)$ as defined in eqn.~(17) of \cite{prada_halo_2011}.
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: sigmaPrime
    
    prada2011C=self%A*((sigmaPrime/self%B)**self%C+1.0d0)*exp(self%D/sigmaPrime**2)
  end function prada2011C

  function prada2011DensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of concentration in the \cite{prada_halo_2011} algorithm.
    implicit none
    class(virialDensityContrastClass             ), pointer       :: prada2011DensityContrastDefinition
    class(darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    
    allocate(virialDensityContrastFixed :: prada2011DensityContrastDefinition)
    select type (prada2011DensityContrastDefinition)
    type is (virialDensityContrastFixed)
      prada2011DensityContrastDefinition=virialDensityContrastFixed(200.0d0,virialDensityContrastFixedDensityTypeCritical)
    end select
    return
  end function prada2011DensityContrastDefinition

  function prada2011DarkMatterProfileDefinition(self)
    !% Return a dark matter density profile object defining that used in the definition of concentration in the
    !% \cite{prada_halo_2011} algorithm.
    use Dark_Matter_Halo_Scales
    implicit none
    class(darkMatterProfileClass                            ), pointer       :: prada2011DarkMatterProfileDefinition
    class(darkMatterProfileConcentrationPrada2011           ), intent(inout) :: self
    class(darkMatterHaloScaleVirialDensityContrastDefinition), pointer       :: darkMatterHaloScaleDefinition

    allocate(darkMatterProfileNFW                               :: prada2011DarkMatterProfileDefinition)
    allocate(darkMatterHaloScaleVirialDensityContrastDefinition :: darkMatterHaloScaleDefinition       )
    select type (prada2011DarkMatterProfileDefinition)
    type is (darkMatterProfileNFW)
       select type (darkMatterHaloScaleDefinition)
       type is (darkMatterHaloScaleVirialDensityContrastDefinition)
          darkMatterHaloScaleDefinition       =darkMatterHaloScaleVirialDensityContrastDefinition(self%densityContrastDefinition())
          prada2011DarkMatterProfileDefinition=darkMatterProfileNFW                              (darkMatterHaloScaleDefinition   )
       end select
    end select
    return
  end function prada2011DarkMatterProfileDefinition
