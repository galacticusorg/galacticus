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

!% Contains a module which implements a molecular ratio class that assumes the model of \cite{obreschkow_simulation_2009}.

  !# <outputAnalysisMolecularRatio name="outputAnalysisMolecularRatioObreschkow2009">
  !#  <description>A high-pass filter analysis property operator class.</description>
  !# </outputAnalysisMolecularRatio>
  type, extends(outputAnalysisMolecularRatioClass) :: outputAnalysisMolecularRatioObreschkow2009
     !% A multiplication property operator class.
     private
     double precision :: K     , fSigma , &
          &              A1    , A2     , &
          &              alpha1, alpha2 , &
          &              beta  , scatter
   contains
     procedure :: ratio        => obreschkow2009Ratio
     procedure :: ratioScatter => obreschkow2009RatioScatter
  end type outputAnalysisMolecularRatioObreschkow2009

  interface outputAnalysisMolecularRatioObreschkow2009
     !% Constructors for the ``obreschkow2009'' output analysis class.
     module procedure obreschkow2009ConstructorParameters
     module procedure obreschkow2009ConstructorInternal
  end interface outputAnalysisMolecularRatioObreschkow2009

contains

  function obreschkow2009ConstructorParameters(parameters) result(self)
    !% Constructor for the ``obreschkow2009'' output analysis property operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(outputAnalysisMolecularRatioObreschkow2009)                :: self
    type(inputParameters                           ), intent(inout) :: parameters
    double precision                                                :: K         , fSigma , &
         &                                                             A1        , A2     , &
         &                                                             alpha1    , alpha2 , &
         &                                                             beta      , scatter

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>K</name>
    !#   <defaultValue>11.3d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The parameter, $K$ (in units of m$^4$ kg$^{-2}$), appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>fSigma</name>
    !#   <defaultValue>0.4d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The parameter, $\langle f_\sigma \rangle$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>A1</name>
    !#   <defaultValue>3.44d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The parameter, $A_1$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>A2</name>
    !#   <defaultValue>4.82d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>
    !#    The parameter, $A_2$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
    !#   </description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha1</name>
    !#   <defaultValue>-0.506d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The parameter, $\alpha_1$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha2</name>
    !#   <defaultValue>-1.054d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The parameter, $\alpha_2$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <defaultValue>0.8d0</defaultValue>
    !#   <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The parameter, $\beta$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>scatter</name>
    !#   <defaultValue>0.4d0</defaultValue>
    !#   <defaultSource>(Obsreschkow, private communication)</defaultSource>
    !#   <source>parameters</source>
    !#   <description>The scatter (in dex) in the molecular ratio $\log_{10}R_\mathrm{mol}$ of \cite{obreschkow_simulation_2009} compared to observational data.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    self=outputAnalysisMolecularRatioObreschkow2009(K,fSigma,A1,A2,alpha1,alpha2,beta,scatter)
    !# <inputParametersValidate source="parameters"/>
    return
  end function obreschkow2009ConstructorParameters

  function obreschkow2009ConstructorInternal(K,fSigma,A1,A2,alpha1,alpha2,beta,scatter) result (self)
    !% Internal constructor for the ``obreschkow2009'' output analysis distribution operator class.
    implicit none
    type            (outputAnalysisMolecularRatioObreschkow2009)                :: self
    double precision                                            , intent(in   ) :: K     , fSigma , &
         &                                                                         A1    , A2     , &
         &                                                                         alpha1, alpha2 , &
         &                                                                         beta  , scatter
    !# <constructorAssign variables="K,fSigma,A1,A2,alpha1,alpha2,beta,scatter"/>

    return
  end function obreschkow2009ConstructorInternal

  double precision function obreschkow2009Ratio(self,massISM,node)
    !% Compute the molecular fraction in the {\normalfont \ttfamily obreschkow2009} class.
    use :: Galacticus_Nodes                , only : nodeComponentDisk, treeNode
    use :: Numerical_Constants_Astronomical, only : massSolar        , megaParsec
    implicit none
    class           (outputAnalysisMolecularRatioObreschkow2009), intent(inout) :: self
    double precision                                            , intent(in   ) :: massISM
    type            (treeNode                                  ), intent(inout) :: node
    class           (nodeComponentDisk                         ), pointer       :: disk
    double precision                                                            :: massStellar          , diskRadius, &
         &                                                                         molecularRatioCentral

    ! Get disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    diskRadius  =  disk%radius     ()
    ! Get the molecular to atomic mass ratio (H2/HI).
    if (diskRadius <= 0.0d0 .or. massStellar < 0.0d0) then
       obreschkow2009Ratio=1.0d0
    else
       molecularRatioCentral   &
            & =(               &
            &   massSolar  **2 &
            &   /megaParsec**4 &
            &   *self%K        &
            &   *massISM       &
            &   *(             &
            &     +massISM     &
            &     +self%fSigma &
            &     *massStellar &
            &    )             &
            &   /diskRadius**4 &
            &  )**self%beta
       if (molecularRatioCentral > 0.0d0) then
          obreschkow2009Ratio=+1.0d0                                &
               &              /(                                    &
               &                +                       self%A1     &
               &                *molecularRatioCentral**self%alpha1 &
               &                +                       self%A2     &
               &                *molecularRatioCentral**self%alpha2 &
               &               )
       else
          obreschkow2009Ratio=0.0d0
       end if
    end if
    return
  end function obreschkow2009Ratio

  double precision function obreschkow2009RatioScatter(self,massISM,node)
    !% Compute the scatter in molecular fraction in the {\normalfont \ttfamily obreschkow2009} class.
    implicit none
    class           (outputAnalysisMolecularRatioObreschkow2009), intent(inout) :: self
    double precision                                            , intent(in   ) :: massISM
    type            (treeNode                                  ), intent(inout) :: node
    !$GLC attributes unused :: massISM, node

    obreschkow2009RatioScatter=self%scatter
    return
  end function obreschkow2009RatioScatter
