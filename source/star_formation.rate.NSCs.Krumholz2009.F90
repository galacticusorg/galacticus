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

  !!{
  Implementation of the \cite{Antonini_2015} star formation rate law for galactic NSCs.
  !!}
  use :: Abundances_Structure, only : abundances

  !![
  <starFormationRateNSCs name="starFormationRateNSCsKrumholz2009">
   <description>
    A star formation rate implementing the model of \citep{Antonini_2015} for galactic NSCs. This model
    uses the \citep{krumholz_star_2009} star formation rule, with minor modifications.
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
    and $\Sigma_1= \Sigma_\mathrm{gas}^\mathrm{NSC}/M_\odot \hbox{pc}^{-2}$ where $\Sigma_\mathrm{gas}^\mathrm{NSC}=\frac{M_\mathrm{gas}^{NSC}}{2\pi r^\mathrm{NSC}$
    is the surface density of the NSC gas reservoir. The timescale is given by 
    \begin{equation}
    t_\mathrm{SF}^{-1} = (2.6~\mathrm{Gyr})^{-1}\times \left\{ \begin{array}{cc} \left(\frac{\Sigma_\mathrm{res}}{\Sigma_\mathrm{th}} \right) ^{-0.33}, &amp;
    \Sigma_\mathrm{res} \le \Sigma_\mathrm{th} \\  \left(\frac{\Sigma_\mathrm{res}}{\Sigma_\mathrm{th}} \right) ^{0.34}, &amp; \Sigma_\mathrm{res} &gt; \Sigma_\mathrm{th} \end{array}  \right. ,
    \end{equation}
    with $\Sigma_\mathrm{th}=85M_\odot\,\box{pc}^{-2}$
   </description>
  </starFormationRateNSCs>
  !!]
  type, extends(starFormationRateNSCsClass) :: starFormationRateNSCsKrumholz2009
     !!{
     Implementation of the \cite{krumholz_star_2009} star formation rate surface density law for galactic NSCs.
     !!}
     private
     double precision   ::  metallicityRelativeToSolar, frequencyStarFormation, &
          &                 s                         , chi                            
     contains
     procedure :: rate                  => krumholz2009Rate
  end type starFormationRateNSCsKrumholz2009

  interface starFormationRateNSCsKrumholz2009
     !!{
     Constructors for the {\normalfont \ttfamily krumholz2009} star formation surface density rate in NSCs class.
     !!}
     module procedure krumholz2009ConstructorParameters
     module procedure krumholz2009ConstructorInternal
  end interface starFormationRateNSCsKrumholz2009
    
contains

  function krumholz2009ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily krumholz2009} star formation surface density rate in NSCs class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (starFormationRateNSCsKrumholz2009)                :: self
    type            (inputParameters                  ), intent(inout) :: parameters
    double precision                                                   :: frequencyStarFormation         

    !![
    <inputParameter>
      <name>frequencyStarFormation</name>
      <defaultSource>\citep{krumholz_star_2009}</defaultSource>
      <defaultValue>2.36d0</defaultValue>
      <description>The star formation frequency (in units of Gyr).</description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=starFormationRateNSCsKrumholz2009(frequencyStarFormation)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function krumholz2009ConstructorParameters

  function krumholz2009ConstructorInternal(frequencyStarFormation) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily krumholz2009} star formation surface density rate from NSCs class.
    !!}
    implicit none
    type            (starFormationRateNSCsKrumholz2009)                 :: self
    double precision                                  , intent(in   )  :: frequencyStarFormation
    !![
    <constructorAssign variables="frequencyStarFormation"/>
    !!]
    return
  end function krumholz2009ConstructorInternal

  double precision function krumholz2009Rate(self,node)
    !!{
    Returns the star formation rate (in $M_\odot$ Gyr$^{-1}$) for star formation
    in the galactic NSC of {\normalfont \ttfamily node}. The NSC is assumed to obey the
    \cite{krumholz_star_2009} star formation rule.
    !!}
    use :: Galacticus_Nodes    , only : nodeComponentNSC                , treeNode
    use :: Abundances_Structure, only : metallicityTypeLinearByMassSolar
    implicit none
    class           (starFormationRateNSCsKrumholz2009), intent(inout), target :: self
    type            (treeNode                         ), intent(inout)         :: node
    class           (nodeComponentNSC                 ), pointer               :: NSC
    type            (abundances                       ), save                  :: abundancesFuel
    double precision                                                           :: molecularGasFraction, radiusNSC   , &
         &                                                                        massGasNSC          , timescale_SF                       
    double precision                                   , parameter             :: Sigma_th = 85.0d0                   !M⊙ pc⁻²
    double precision                                                           :: surfaceDensityGasNSC, Sigma_1
                        
    !$omp threadprivate(abundancesFuel)


    NSC       => node%NSC    ()
    massGasNSC=  NSC %massGas()
    radiusNSC =  NSC %radius ()*1.0d6 !pc

    if     (                                             &
         &   massGasNSC                         <= 0.0d0 &
         &  .or.                                         &
         &   radiusNSC                          <= 0.0d0 &
         & ) then
       ! It is not, so return zero rate.
       krumholz2009Rate=0.0d0
       return
    else 
       surfaceDensityGasNSC = SurfaceDensityGas(radiusNSC,massGasNSC)
       Sigma_1              = surfaceDensityGasNSC/1.0d0

       ! Find the hydrogen fraction in the NSC gas of the fuel supply.
       abundancesFuel=NSC%abundancesGas()
       call abundancesFuel%massToMassFraction(massGasNSC)

       ! Get the metallicity in Solar units, and related quantities.
       self%metallicityRelativeToSolar=abundancesFuel%metallicity(metallicityTypeLinearByMassSolar)
       if (self%metallicityRelativeToSolar /= 0.0d0) then
           self%chi                         = 0.77d0*(1.0d0+3.1d0*self%metallicityRelativeToSolar**0.365d0)
           self%s                           = log(1.0d0+0.6d0*self%chi)/(0.04d0*self%metallicityRelativeToSolar*Sigma_1)
       else
            krumholz2009Rate=0.0d0
       end if 
       
       !Computation of the timescale given by Krumholz et al. (2009; https://ui.adsabs.harvard.edu/abs/2009ApJ...699..850K/abstract)
       !and Sesana et al. (2014, https://ui.adsabs.harvard.edu/abs/2014ApJ...794..104S/abstract) This variable is in units of Gyr⁻¹
       if (surfaceDensityGasNSC > Sigma_th ) then 
           timescale_SF = (1.0d0/self%frequencyStarFormation)*(surfaceDensityGasNSC/Sigma_th)**(-0.33d0)  
       else
           timescale_SF = (1.0d0/self%frequencyStarFormation)*(surfaceDensityGasNSC/Sigma_th)**0.34d0 
       end if 
      ! Compute the molecular fraction.
       molecularGasFraction = MolecularFraction(self%s)
       !Compute the star formation rate density.
       krumholz2009Rate     = +massGasNSC*molecularGasFraction*timescale_SF
    end if
    return
  end function krumholz2009Rate
  
  double precision function SurfaceDensityGas(radiusNSC,massGasNSC)
    !!{
    Compute surface density of the NSC.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    double precision               , intent(in   ) :: radiusNSC, massGasNSC
    ! Return the surface gas density in units of M☉/pc² of the NSC
    SurfaceDensityGas = massGasNSC / (4.0d0*Pi*radiusNSC**2.0d0)
    return
  end function SurfaceDensityGas

  double precision function MolecularFraction(s)
    !!{
    Slow (but more accurate at low molecular fraction) fitting function from \cite{krumholz_star_2009} for the molecular
    hydrogen fraction.
    !!}
    implicit none
    double precision, intent(in   ) :: s
    double precision, parameter     :: sTiny        =1.000000d-06
    double precision, parameter     :: sHuge        =1.000000d+10
    double precision, parameter     :: deltaInfinity=0.214008d+00 ! The value of δ for s → ∞.
    double precision, parameter     :: sMaximum     =10.0d0 
    double precision                :: delta

    if      (s <  sTiny   ) then
       ! Series expansion for very small s.
       MolecularFraction=1.0d0-0.75d0*s
    else if (s >= sHuge   ) then
       ! Truncate to zero for extremely large s.
       MolecularFraction=0.0d0
    else if (s >= sMaximum) then
       ! Simplified form for very large s.
       MolecularFraction=1.0d0/(0.75d0/(1.0d0+deltaInfinity))**5.0d0/5.0d0/s**5.0d0
    else
       ! Full expression.
       delta            =0.0712d0/((0.1d0/s+0.675d0)**2.8d0)
       MolecularFraction=1.0d0-1.0d0/((1.0d0+(((1.0d0+delta)/0.75d0/s)**5.0d0))**0.2d0)
    end if

    if (MolecularFraction > 0.02d0) then
        return
    else 
       MolecularFraction = 0.02d0
    end if 
    return
  end function MolecularFraction
