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

!!{
Implements a molecular ratio class that assumes the model of \cite{obreschkow_simulation_2009}.
!!}

  !![
  <outputAnalysisMolecularRatio name="outputAnalysisMolecularRatioObreschkow2009">
   <description>
    A molecular ratio class which computes the molecular ratio. The class assumes that only the total \gls{ism} mass of each
    galaxy is available, along with the disk radius (assuming an exponential disk). To infer the HI mass the model of
    \cite{obreschkow_simulation_2009} is used. Specifically, the molecular ratio, $R_\mathrm{mol}\equiv
    M_\mathrm{H_2}/M_\mathrm{HI}$, is given by:
    \begin{equation}
     R_\mathrm{mol} = \left( A_1 R_\mathrm{mol}^\mathrm{c\,\alpha_1} + A_2 R_\mathrm{mol}^\mathrm{c\,\alpha_2} \right)^{-1},
     \label{eq:HIMassSystematic}
    \end{equation}
    where the ratio at the disk center is given by    
    \begin{equation}
     R_\mathrm{mol}^\mathrm{c} = [ K r_\mathrm{disk}^{-4} M_\mathrm{gas} (M_\mathrm{gas} + \langle f_\sigma \rangle M_\star)]^\beta.
    \end{equation}
      
    Here, $R_\mathrm{mol}$ is the mass ratio of H$_2$ to HI, $M_\star$ is the stellar mass of the disk, $r_\mathrm{disk}$ is the
    disk exponential scale length, $\langle f_\sigma \rangle$ is the average ratio of the vertical velocity dispersions of gas to
    stars, and $K=\mathrm{G}/(8\pi P_\star)$. The HI mass is then determined from:
    \begin{equation}
     M_\mathrm{HI} = X_\mathrm{H} M_\mathrm{gas} / ( 1 + R_\mathrm{mol} ),
    \end{equation}
    where $X_\mathrm{H}=0.778$ is the primordial hydrogen fraction by mass. In the above $K=${\normalfont \ttfamily [K]},
    $\langle f_\sigma \rangle=${\normalfont \ttfamily [fSigma]}, $A_1=${\normalfont \ttfamily [A1]}, $A_2=${\normalfont \ttfamily
    [A2]}, $\alpha_1=${\normalfont \ttfamily [alpha1]}, $\alpha_2=${\normalfont \ttfamily [alpha2]}, and $\beta=${\normalfont
    \ttfamily [beta]}. Default values for these parameters are taken from \cite{obreschkow_simulation_2009}. According to
    Obreschkow (private communication), there remains significant scatter of $\sigma_{R_\mathrm{mol}}=0.4$~dex between the
    predicted $R_\mathrm{mol}$ from this model and that observed. This is accounted for in when constructing the mass function
    (see below).
   </description>
  </outputAnalysisMolecularRatio>
  !!]
  type, extends(outputAnalysisMolecularRatioClass) :: outputAnalysisMolecularRatioObreschkow2009
     !!{
     A multiplication property operator class.
     !!}
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
     !!{
     Constructors for the \refClass{outputAnalysisMolecularRatioObreschkow2009} output analysis class.
     !!}
     module procedure obreschkow2009ConstructorParameters
     module procedure obreschkow2009ConstructorInternal
  end interface outputAnalysisMolecularRatioObreschkow2009

contains

  function obreschkow2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisMolecularRatioObreschkow2009} output analysis property operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(outputAnalysisMolecularRatioObreschkow2009)                :: self
    type(inputParameters                           ), intent(inout) :: parameters
    double precision                                                :: K         , fSigma , &
         &                                                             A1        , A2     , &
         &                                                             alpha1    , alpha2 , &
         &                                                             beta      , scatter

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>K</name>
      <defaultValue>11.3d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>The parameter, $K$ (in units of m$^4$ kg$^{-2}$), appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    </inputParameter>
    <inputParameter>
      <name>fSigma</name>
      <defaultValue>0.4d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>The parameter, $\langle f_\sigma \rangle$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    </inputParameter>
    <inputParameter>
      <name>A1</name>
      <defaultValue>3.44d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>The parameter, $A_1$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    </inputParameter>
    <inputParameter>
      <name>A2</name>
      <defaultValue>4.82d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>
       The parameter, $A_2$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.
      </description>
    </inputParameter>
    <inputParameter>
      <name>alpha1</name>
      <defaultValue>-0.506d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>The parameter, $\alpha_1$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    </inputParameter>
    <inputParameter>
      <name>alpha2</name>
      <defaultValue>-1.054d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>The parameter, $\alpha_2$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <defaultValue>0.8d0</defaultValue>
      <defaultSource>\citep{obreschkow_simulation_2009}</defaultSource>
      <source>parameters</source>
      <description>The parameter, $\beta$, appearing in the model for the H$_2$/HI ratio in galaxies from \cite{obreschkow_simulation_2009}.</description>
    </inputParameter>
    <inputParameter>
      <name>scatter</name>
      <defaultValue>0.4d0</defaultValue>
      <defaultSource>(Obsreschkow, private communication)</defaultSource>
      <source>parameters</source>
      <description>The scatter (in dex) in the molecular ratio $\log_{10}R_\mathrm{mol}$ of \cite{obreschkow_simulation_2009} compared to observational data.</description>
    </inputParameter>
    !!]
    self=outputAnalysisMolecularRatioObreschkow2009(K,fSigma,A1,A2,alpha1,alpha2,beta,scatter)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function obreschkow2009ConstructorParameters

  function obreschkow2009ConstructorInternal(K,fSigma,A1,A2,alpha1,alpha2,beta,scatter) result (self)
    !!{
    Internal constructor for the \refClass{outputAnalysisMolecularRatioObreschkow2009} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisMolecularRatioObreschkow2009)                :: self
    double precision                                            , intent(in   ) :: K     , fSigma , &
         &                                                                         A1    , A2     , &
         &                                                                         alpha1, alpha2 , &
         &                                                                         beta  , scatter
    !![
    <constructorAssign variables="K,fSigma,A1,A2,alpha1,alpha2,beta,scatter"/>
    !!]

    return
  end function obreschkow2009ConstructorInternal

  double precision function obreschkow2009Ratio(self,massISM,node)
    !!{
    Compute the molecular fraction in the {\normalfont \ttfamily obreschkow2009} class.
    !!}
    use :: Galacticus_Nodes                , only : nodeComponentDisk, treeNode
    use :: Numerical_Constants_Astronomical, only : massSolar        , megaParsec
    implicit none
    class           (outputAnalysisMolecularRatioObreschkow2009), intent(inout) :: self
    double precision                                            , intent(in   ) :: massISM
    type            (treeNode                                  ), intent(inout) :: node
    class           (nodeComponentDisk                         ), pointer       :: disk
    double precision                                                            :: massStellar          , radiusDisk    , &
         &                                                                         molecularRatioCentral

    ! Get disk properties.
    disk        => node%disk       ()
    massStellar =  disk%massStellar()
    radiusDisk  =  disk%radius     ()
    ! Get the molecular to atomic mass ratio (H₂/HI).
    if (radiusDisk <= 0.0d0 .or. massStellar < 0.0d0) then
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
            &   /radiusDisk**4 &
            &  )**self%beta
       if (molecularRatioCentral > 0.0d0) then
          ! Find the H₂/HI ratio using eqn. (15) of Obreschkow et al. (2009).
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
    !!{
    Compute the scatter in molecular fraction in the {\normalfont \ttfamily obreschkow2009} class.
    !!}
    implicit none
    class           (outputAnalysisMolecularRatioObreschkow2009), intent(inout) :: self
    double precision                                            , intent(in   ) :: massISM
    type            (treeNode                                  ), intent(inout) :: node
    !$GLC attributes unused :: massISM, node

    obreschkow2009RatioScatter=self%scatter
    return
  end function obreschkow2009RatioScatter
