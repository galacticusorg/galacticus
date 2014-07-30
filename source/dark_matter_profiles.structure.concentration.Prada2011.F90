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
  
  !% An implementation of dark matter halo profile concentrations using the \cite{prada_halo_2011} algorithm.
  
  !# <darkMatterProfileConcentration name="darkMatterProfileConcentrationPrada2011">
  !#  <description>Dark matter halo concentrations are computed using the algorithm of \cite{prada_halo_2011}.</description>
  !# </darkMatterProfileConcentration>
  
  ! Parameters of the fit.
  logical          :: prada2011Initialized               =.false.
  double precision :: prada2011ConcentrationA                    , prada2011ConcentrationB            , &
       &              prada2011ConcentrationC                    , prada2011ConcentrationD            , &
       &              prada2011ConcentrationC0                   , prada2011ConcentrationC1           , &
       &              prada2011ConcentrationX0                   , prada2011ConcentrationInverseSigma0, &
       &              prada2011ConcentrationInverseSigma1        , prada2011ConcentrationAlpha        , &
       &              prada2011ConcentrationBeta                 , prada2011ConcentrationX1

  type, extends(darkMatterProfileConcentrationClass) :: darkMatterProfileConcentrationPrada2011
     !% A dark matter halo profile concentration class implementing the algorithm of \cite{prada_halo_2011}.
     private
     double precision :: A, B, C, D, C0, C1, X0, X1, InverseSigma0, InverseSigma1, Alpha, Beta
   contains
     procedure :: concentration               => prada2011Concentration
     procedure :: densityContrastDefinition   => prada2011DensityContrastDefinition
     procedure :: darkMatterProfileDefinition => prada2011DarkMatterProfileDefinition
  end type darkMatterProfileConcentrationPrada2011
  
  interface darkMatterProfileConcentrationPrada2011
     !% Constructors for the {\tt prada2011} dark matter halo profile concentration class.
     module procedure prada2011DefaultConstructor
     module procedure prada2011Constructor
  end interface darkMatterProfileConcentrationPrada2011
  
contains
  
  function prada2011DefaultConstructor()
    !% Default constructor for the {\tt prada2011} dark matter halo profile concentration class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileConcentrationPrada2011), target  :: prada2011DefaultConstructor
    
    if (.not.prada2011Initialized) then
       !$omp critical(prada2011DefaultInitialize)
       if (.not.prada2011Initialized) then
          ! Get parameters of the model.
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationA</name>
          !@   <defaultValue>2.881 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $A$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationA",prada2011ConcentrationA,defaultValue=2.881d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationB</name>
          !@   <defaultValue>1.257 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $b$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationB",prada2011ConcentrationB,defaultValue=1.257d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationC</name>
          !@   <defaultValue>1.022 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $c$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationC",prada2011ConcentrationC,defaultValue=1.022d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationD</name>
          !@   <defaultValue>0.060 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $d$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationD",prada2011ConcentrationD,defaultValue=0.060d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationC0</name>
          !@   <defaultValue>3.681 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $c_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationC0",prada2011ConcentrationC0,defaultValue=3.681d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationC1</name>
          !@   <defaultValue>5.033 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $c_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationC1",prada2011ConcentrationC1,defaultValue=5.033d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationX0</name>
          !@   <defaultValue>0.424 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $x_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationX0",prada2011ConcentrationX0,defaultValue=0.424d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationX1</name>
          !@   <defaultValue>0.526 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $x_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationX1",prada2011ConcentrationX1,defaultValue=0.526d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationInverseSigma0</name>
          !@   <defaultValue>1.047 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\sigma^{-1}_0$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationInverseSigma0",prada2011ConcentrationInverseSigma0,defaultValue=1.047d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationInverseSigma1</name>
          !@   <defaultValue>1.646 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\sigma^{-1}_1$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationInverseSigma1",prada2011ConcentrationInverseSigma1,defaultValue=1.646d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationAlpha</name>
          !@   <defaultValue>6.948 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\alpha$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationAlpha",prada2011ConcentrationAlpha,defaultValue=6.948d0)
          !@ <inputParameter>
          !@   <name>prada2011ConcentrationBeta</name>
          !@   <defaultValue>7.386 \cite{prada_halo_2011}</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The parameter $\beta$ appearing in the halo concentration algorithm of \cite{prada_halo_2011}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter("prada2011ConcentrationBeta",prada2011ConcentrationBeta,defaultValue=7.386d0)
          ! Record that method is now initialized.
          prada2011Initialized=.true.
       end if
       !$omp end critical(prada2011DefaultInitialize)
    end if
    ! Construct the object.
    prada2011DefaultConstructor=&
         & prada2011Constructor(&
         &                      prada2011ConcentrationA            , &
         &                      prada2011ConcentrationB            , &
         &                      prada2011ConcentrationC            , &
         &                      prada2011ConcentrationD            , &
         &                      prada2011ConcentrationC0           , &
         &                      prada2011ConcentrationC1           , &
         &                      prada2011ConcentrationX0           , &
         &                      prada2011ConcentrationX1           , &
         &                      prada2011ConcentrationInverseSigma0, &
         &                      prada2011ConcentrationInverseSigma1, &
         &                      prada2011ConcentrationAlpha        , &
         &                      prada2011ConcentrationBeta           &
         &                     )
    return
  end function prada2011DefaultConstructor
  
  function prada2011Constructor(A,B,C,D,C0,C1,X0,X1,InverseSigma0,InverseSigma1,Alpha,Beta)
    !% Constructor for the {\tt prada2011} dark matter halo profile concentration class.
    implicit none
    type            (darkMatterProfileConcentrationPrada2011)                :: prada2011Constructor
    double precision                                         , intent(in   ) :: A                   , B            , &
         &                                                                      C                   , D            , &
         &                                                                      C0                  , C1           , &
         &                                                                      X0                  , InverseSigma0, &
         &                                                                      InverseSigma1       , Alpha        , &
         &                                                                      Beta                , X1
    
    prada2011Constructor%A            =A
    prada2011Constructor%B            =B
    prada2011Constructor%C            =C
    prada2011Constructor%D            =D
    prada2011Constructor%C0           =C0
    prada2011Constructor%C1           =C1
    prada2011Constructor%X0           =X0
    prada2011Constructor%X1           =X1
    prada2011Constructor%InverseSigma0=InverseSigma0
    prada2011Constructor%InverseSigma1=InverseSigma1
    prada2011Constructor%Alpha        =Alpha
    prada2011Constructor%Beta         =Beta
    return
  end function prada2011Constructor
  
  double precision function prada2011Concentration(self,node)
    !% Return the concentration of the dark matter halo profile of {\tt node} using the \cite{prada_halo_2011} algorithm.
    use Power_Spectra
    use Cosmology_Functions
    use Cosmology_Parameters
    use Linear_Growth
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout)          :: self
    type            (treeNode                               ), intent(inout), pointer :: node
    class           (nodeComponentBasic                     )               , pointer :: basic
    class           (cosmologyParametersClass               )               , pointer :: cosmologyParameters_
    class           (cosmologyFunctionsClass                )               , pointer :: cosmologyFunctions_
    double precision                                                                  :: massNode            , sigmaPrime, &
         &                                                                               timeNode            , x
    
    ! Get the default cosmology functions object.
    cosmologyFunctions_   => cosmologyFunctions ()
    ! Get the default cosmology.
    cosmologyParameters_  => cosmologyParameters()
    basic                 => node %basic()
    massNode              =  basic%mass ()
    timeNode              =  basic%time ()
    x                     =  (cosmologyParameters_%OmegaDarkEnergy()/cosmologyParameters_%OmegaMatter())**(1.0d0/3.0d0)*cosmologyFunctions_%expansionFactor(timeNode)
    sigmaPrime            =  prada2011B1(self,x)*Cosmological_Mass_Root_Variance(massNode)*Linear_Growth_Factor(timeNode)
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
    
    prada2011B1=prada2011InverseSigmaMin(self,x)/prada2011InverseSigmaMin(self,1.393d0)
  end function prada2011B1
  
  double precision function prada2011cMin(self,x)
    !% The function $c_{\rm min}(x)$ as defined in eqn.~(19) of \cite{prada_halo_2011}.
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x
    
    prada2011cMin=self%C0+(self%C1-self%C0)*(atan(self%Alpha*(x-self%X0))/Pi+0.5d0)
  end function prada2011cMin
  
  double precision function prada2011InverseSigmaMin(self,x)
    !% The function $\sigma^{-1}_{\rm min}(x)$ as defined in eqn.~(20) of \cite{prada_halo_2011}.
    use Numerical_Constants_Math
    implicit none
    class           (darkMatterProfileConcentrationPrada2011), intent(inout) :: self
    double precision                                         , intent(in   ) :: x
    
    prada2011InverseSigmaMin=self%InverseSigma0+(self%InverseSigma1-self%InverseSigma0)*(atan(self%Beta*(x-self%X1))/Pi+0.5d0)
  end function prada2011InverseSigmaMin
  
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
