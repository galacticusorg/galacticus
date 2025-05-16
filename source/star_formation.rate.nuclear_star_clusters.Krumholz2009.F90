!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024
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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implementation of the \cite{antonini_coevolution_2015} star formation rate law for galactic \glspl{nsc}.
  !!}
  use :: Abundances_Structure, only : abundances

  !![
  <starFormationRateNuclearStarClusters name="starFormationRateNuclearStarClustersKrumholz2009">
   <description>
    A star formation rate implementing the model of \citep{antonini_coevolution_2015} for galactic \glspl{nsc}. This model
    uses the \citep{krumholz_star_2009} star formation rule:
    \begin{equation}
     \dot{M}_\star^\mathrm{NSC} = f_c\frac{M_\mathrm{gas}^\mathrm{gas}}{t_{SF}},
    \end{equation}
    where $f_c$ is the fraction of cold gas available for star formation given by
    \begin{equation}
     f_c = 1 - \left( 1 + \left[ { 3 s \over 4 (1+\delta)} \right]^{-5} \right)^{-1/5},
    \end{equation}
     if $f_c > 0.02$ and $f_c = 0.02 $ otherwise, with 
    \begin{equation}
     \delta = 0.0712 \left[ 0.1 s^{-1} + 0.675 \right]^{-2.8},
    \end{equation}
    and
    \begin{equation}
     s = {\ln(1+0.6\chi+0.01\chi) \over 0.04 \Sigma_1 Z^\prime},
    \end{equation}
    with
    \begin{equation}
     \chi = 0.77 \left[ 1 + 3.1 Z^{\prime 0.365} \right],
    \end{equation}
    and $\Sigma_1= \Sigma_\mathrm{gas}^\mathrm{NSC}/M_\odot \hbox{pc}^{-2}$ where $\Sigma_\mathrm{gas}^\mathrm{NSC}=\frac{M_\mathrm{gas}^{NSC}}{4\pi r^\mathrm{NSC}}$
    is the surface density of the NSC gas reservoir. The timescale is given by 
    \begin{equation}
    t_\mathrm{SF}^{-1} = (2.36~\mathrm{Gyr})^{-1}\times \left\{ \begin{array}{cc} \left(\frac{\Sigma_\mathrm{res}}{\Sigma_\mathrm{th}} \right) ^{-0.33}, &amp;
    \Sigma_\mathrm{res} \le \Sigma_\mathrm{th} \\  \left(\frac{\Sigma_\mathrm{res}}{\Sigma_\mathrm{th}} \right) ^{0.33}, &amp; \Sigma_\mathrm{res} &gt; \Sigma_\mathrm{th} \end{array}  \right. ,
    \end{equation}
    with $\Sigma_\mathrm{th}=85\mathrm{M}_\odot\,\hbox{pc}^{-2}$
   </description>
  </starFormationRateNuclearStarClusters>
  !!]
  type, extends(starFormationRateNuclearStarClustersClass) :: starFormationRateNuclearStarClustersKrumholz2009
     !!{
     Implementation of the \cite{krumholz_star_2009} star formation rate law for galactic \glspl{nsc}.
     !!}
     private
     double precision :: frequencyStarFormation                       
     contains
     procedure :: rate => krumholz2009Rate
  end type starFormationRateNuclearStarClustersKrumholz2009

  interface starFormationRateNuclearStarClustersKrumholz2009
     !!{
     Constructors for the \refClass{starFormationRateNuclearStarClustersKrumholz2009} star formation rate law for galactic \glspl{nsc}.
     !!}
     module procedure krumholz2009ConstructorParameters
     module procedure krumholz2009ConstructorInternal
  end interface starFormationRateNuclearStarClustersKrumholz2009
    
contains

  function krumholz2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{starFormationRateNuclearStarClustersKrumholz2009} star formation rate law for galactic \glspl{nsc} which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationRateNuclearStarClustersKrumholz2009)                :: self
    type            (inputParameters                                 ), intent(inout) :: parameters
    double precision                                                                  :: frequencyStarFormation         

    !![
    <inputParameter>
      <name>frequencyStarFormation</name>
      <defaultSource>\citep{krumholz_star_2009}</defaultSource>
      <defaultValue>2.36d0</defaultValue>
      <description>The star formation frequency (in units of Gyr).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateNuclearStarClustersKrumholz2009(frequencyStarFormation)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function krumholz2009ConstructorParameters

  function krumholz2009ConstructorInternal(frequencyStarFormation) result(self)
    !!{
    Internal constructor for the \refClass{starFormationRateNuclearStarClustersKrumholz2009} star formation rate law for galactic \glspl{nsc}.
    !!}
    implicit none
    type            (starFormationRateNuclearStarClustersKrumholz2009)                 :: self
    double precision                                                  , intent(in   )  :: frequencyStarFormation
    !![
    <constructorAssign variables="frequencyStarFormation"/>
    !!]
    return
  end function krumholz2009ConstructorInternal

  double precision function krumholz2009Rate(self,node)
    !!{
    Returns the star formation rate (in $M_\odot$ Gyr$^{-1}$) for star formation in the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The \gls{nsc} is assumed to obey the \cite{krumholz_star_2009} star formation rule.
    !!}
    use :: Galacticus_Nodes                          , only : nodeComponentNSC
    use :: Abundances_Structure                      , only : metallicityTypeLinearByMassSolar
    use :: Star_Formation_Rate_Krumholz2009_Utilities, only : krumholz2009MolecularFractionSlow
    implicit none
    class           (starFormationRateNuclearStarClustersKrumholz2009), intent(inout), target :: self
    type            (treeNode                                        ), intent(inout)         :: node
    double precision                                                  , parameter             :: surfaceDensityThreshold    =85.0d0 ! M☉/pc²
    double precision                                                  , parameter             :: surfaceDensityNormalization= 1.0d0 ! M☉/pc²
    class           (nodeComponentNSC                                ), pointer               :: nuclearStarCluster              
    type            (abundances                                      ), save                  :: abundancesFuel
    !$omp threadprivate(abundancesFuel)
    double precision                                                                          :: molecularGasFraction                      , radiusNuclearStarCluster                 , &
         &                                                                                       massGasNuclearStarCluster                 , timescaleStarFormation                   , &
         &                                                                                       surfaceDensityGasNuclearStarCluster       , surfaceDensityGasNuclearStarClusterScaled, &
         &                                                                                       metallicityRelativeToSolar                , chi                                      , &
         &                                                                                       s

    nuclearStarCluster        => node              %NSC    ()
    massGasNuclearStarCluster =  nuclearStarCluster%massGas()
    radiusNuclearStarCluster  =  nuclearStarCluster%radius ()
    if     (                                    &
         &   massGasNuclearStarCluster <= 0.0d0 &
         &  .or.                                &
         &   radiusNuclearStarCluster  <= 0.0d0 &
         & ) then
       ! Unphysical nclear star cluster, so return zero rate.
       krumholz2009Rate=0.0d0
       return
    else 
       surfaceDensityGasNuclearStarCluster      =+surfaceDensityGas(radiusNuclearStarCluster,massGasNuclearStarCluster)
       surfaceDensityGasNuclearStarClusterScaled=+surfaceDensityGasNuclearStarCluster &
            &                                    /surfaceDensityNormalization
       ! Find the hydrogen fraction in the nuclear star cluster gas fuel supply.
       abundancesFuel=nuclearStarCluster%abundancesGas()
       call abundancesFuel%massToMassFraction(massGasNuclearStarCluster)
       ! Get the metallicity in Solar units, and related quantities.
       metallicityRelativeToSolar=abundancesFuel%metallicity(metallicityTypeLinearByMassSolar)
       if (metallicityRelativeToSolar /= 0.0d0) then
           chi              =+0.77d0*(1.0d0+3.1d0*metallicityRelativeToSolar**0.365d0)
           s                =+log(1.0d0+0.6d0*chi)/(0.04d0*metallicityRelativeToSolar*surfaceDensityGasNuclearStarClusterScaled)
       else
            krumholz2009Rate=+0.0d0
       end if 
       ! Computation of the timescale (in units of Gyr⁻¹) given by Krumholz et al. (2009;
       ! https://ui.adsabs.harvard.edu/abs/2009ApJ...699..850K/abstract) and Sesana et al. (2014,
       ! https://ui.adsabs.harvard.edu/abs/2014ApJ...794..104S/abstract).
       if (surfaceDensityGasNuclearStarCluster > surfaceDensityThreshold ) then 
           timescaleStarFormation=(1.0d0/self%frequencyStarFormation)*(surfaceDensityGasNuclearStarCluster/surfaceDensityThreshold)**(-0.33d0)  
       else
           timescaleStarFormation=(1.0d0/self%frequencyStarFormation)*(surfaceDensityGasNuclearStarCluster/surfaceDensityThreshold)**(+0.34d0)
       end if 
       ! Compute the molecular fraction.
       molecularGasFraction=max(0.02d0,krumholz2009MolecularFractionSlow(s))
       !Compute the star formation rate density.
       krumholz2009Rate    =+massGasNuclearStarCluster &
            &               *molecularGasFraction      &
            &               *timescaleStarFormation
    end if
    return
  end function krumholz2009Rate
  
  double precision function surfaceDensityGas(radiusNuclearStarCluster,massGasNuclearStarCluster)
    !!{
    Compute surface density of the nuclear star cluster.
    !!}
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Prefixes, only : mega
    implicit none
    double precision, intent(in   ) :: radiusNuclearStarCluster, massGasNuclearStarCluster
    
    ! Return the surface gas density in units of M☉/pc² of the nuclear star cluster.
    surfaceDensityGas=+massGasNuclearStarCluster          &
         &            /4.0d0                              &
         &            /Pi                                 &
         &            /(mega*radiusNuclearStarCluster)**2
    return
  end function surfaceDensityGas
