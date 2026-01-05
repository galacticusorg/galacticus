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
  Implements a random error output analysis distribution operator class providing errors in HI mass for
  the ALFALFA survey.
  !!}

  use :: Output_Analysis_Molecular_Ratios, only : outputAnalysisMolecularRatioClass

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomErrorALFLF">
   <description>
    A random error output analysis distribution operator class providing errors in HI mass for the ALFALFA survey. To account for
    both observational errors and scatter in $R_\mathrm{mol}$, the HI mass of each galaxy is modeled as a Gaussian in
    $\log_{10}M_\mathrm{HI}$ when constructing the mass function. Observational random errors on HI mass, including those arising
    from flux density uncertainties and errors in the assumed distance to each source, are taken from Fig.~19 of
    \cite{haynes_arecibo_2011}. The magnitude of the error as a function of HI mass is fit using a functional form:
    \begin{equation}
     \sigma_\mathrm{obs} = a + \exp\left(-{\log_{10}(M_\mathrm{HI}/M_\odot)-b\over c}\right),
    \end{equation}
     where $\sigma_\mathrm{obs}$ is the error on $\log_{10}(M_\mathrm{HI}/M_\odot)$. We find a reasonable fit using
     values\footnote{This should not be regarded as a formal good fit. Error estimates are approximate---we have simply found a
     functional form that roughly describes them, along with conservative errors on the parameters of this function which are
     included in the priors.} of $a=${\normalfont \ttfamily a}$=0.100 \pm 0.010$, $b=${\normalfont \ttfamily b}$=5.885 \pm
     0.100$, and $c=${\normalfont \ttfamily c}$=0.505 \pm 0.020$ as shown in Fig.~\ref{fig:ALFALFAErrorModel}. The total random
     error on the logarithm of each galaxy mass is given by $\sigma^2 = \sigma_{R_\mathrm{mol}}^2+\sigma_\mathrm{obs}^2$, and is
     used as the width of the Gaussian kernel when applying each galaxy to the mass function histogram (as described above).
  
    \begin{figure}
     \begin{center}
     \includegraphics[width=85mm,trim=0mm 0mm 0mm 4mm,clip]{Plots/DataAnalysis/alfalfaHIMassErrorModel.pdf}
     \caption{The observational random error in galaxy HI mass as a function of HI mass for the ALFALFA survey. Points show the errors reported by \protect\cite{haynes_arecibo_2011}, while the line shows a simple functional form fit to these errors.}
     \end{center}
     \label{fig:ALFALFAErrorModel}
    \end{figure}
   </description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRandomErrorALFLF
     !!{
     A random error output distribution operator class providing errors in HI mass for the ALFALFA survey.
     !!}
     private
     double precision                                             :: a                            , b, &
          &                                                          c
     class           (outputAnalysisMolecularRatioClass), pointer :: outputAnalysisMolecularRatio_ => null()
   contains
     final     ::                 randomErrorHIALFALFADestructor
     procedure :: rootVariance => randomErrorHIALFALFARootVariance
  end type outputAnalysisDistributionOperatorRandomErrorALFLF

  interface outputAnalysisDistributionOperatorRandomErrorALFLF
     !!{
     Constructors for the \refClass{outputAnalysisDistributionOperatorRandomErrorALFLF} output analysis distribution operator class.
     !!}
     module procedure randomErrorHIALFALFAConstructorParameters
     module procedure randomErrorHIALFALFAConstructorInternal
  end interface outputAnalysisDistributionOperatorRandomErrorALFLF

contains

  function randomErrorHIALFALFAConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{outputAnalysisDistributionOperatorRandomErrorALFLF} output analysis distribution operator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorALFLF)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (outputAnalysisMolecularRatioClass                 ), pointer       :: outputAnalysisMolecularRatio_
    double precision                                                                    :: a                            , b, &
         &                                                                                 c

    ! Check and read parameters.
    !![
    <inputParameter>
      <name>a</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <variable>a</variable>
      <description>Parameter $a$ in the ALFALFA HI mass error model.</description>
    </inputParameter>
    <inputParameter>
      <name>b</name>
      <source>parameters</source>
      <defaultValue>5.885d0</defaultValue>
      <variable>b</variable>
      <description>Parameter $b$ in the ALFALFA HI mass error model.</description>
    </inputParameter>
    <inputParameter>
      <name>c</name>
      <source>parameters</source>
      <defaultValue>0.505d0</defaultValue>
      <variable>c</variable>
      <description>Parameter $c$ in the ALFALFA HI mass error model.</description>
    </inputParameter>
    <objectBuilder class="outputAnalysisMolecularRatio" name="outputAnalysisMolecularRatio_" source="parameters"/>
    !!]
    ! Construct the object.
    self=outputAnalysisDistributionOperatorRandomErrorALFLF(a,b,c,outputAnalysisMolecularRatio_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="outputAnalysisMolecularRatio_"/>
    !!]
    return
  end function randomErrorHIALFALFAConstructorParameters

  function randomErrorHIALFALFAConstructorInternal(a,b,c,outputAnalysisMolecularRatio_) result(self)
    !!{
    Internal constructor for the \refClass{outputAnalysisDistributionOperatorRandomErrorALFLF} output analysis distribution operator class.
    !!}
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorALFLF)                        :: self
    double precision                                                    , intent(in   )         :: a                            , b, &
         &                                                                                         c
    class           (outputAnalysisMolecularRatioClass                 ), intent(in   ), target :: outputAnalysisMolecularRatio_
    !![
    <constructorAssign variables="a, b, c, *outputAnalysisMolecularRatio_"/>
    !!]

    return
  end function randomErrorHIALFALFAConstructorInternal

  subroutine randomErrorHIALFALFADestructor(self)
    !!{
    Destructor for the \refClass{outputAnalysisDistributionOperatorRandomErrorALFLF} output analysis distribution operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorRandomErrorALFLF), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysisMolecularRatio_"/>
    !!]
    return
  end subroutine randomErrorHIALFALFADestructor

  double precision function randomErrorHIALFALFARootVariance(self,propertyValue,node)
    !!{
    Computes errors on $\log_{10}($HI masses$)$ for the ALFALFA survey analysis. Uses a simple fitting function. See
    {\normalfont \ttfamily constraints/dataAnalysis/hiMassFunction\_ALFALFA\_z0.00/alfalfaHIMassErrorModel.pl} for details.
    !!}
    implicit none
    class           (outputAnalysisDistributionOperatorRandomErrorALFLF), intent(inout) :: self
    double precision                                                    , intent(in   ) :: propertyValue
    type            (treeNode                                          ), intent(inout) :: node
    double precision                                                                    :: exponentialArgument, exponentialTerm, &
         &                                                                                 molecularFraction  , mass

    ! Compute the random error on the mass.
    exponentialArgument=+(                          &
         &                +max(propertyValue,6.0d0) &
         &                -self%b                   &
         &               )                          &
         &              /  self%c
    if (exponentialArgument > 0.0d0) then
       exponentialTerm=exp(-exponentialArgument)
    else
       exponentialTerm=0.0d0
    end if
    randomErrorHIALFALFARootVariance=+self%a          &
         &                           +exponentialTerm
    ! Add in quadrature a term accounting for the scatter in the molecular ratio model.
    mass                            =10.0d0**propertyValue
    molecularFraction               =        self%outputAnalysisMolecularRatio_%ratio       (mass,node)
    randomErrorHIALFALFARootVariance=sqrt(                                                                 &
         &                                +randomErrorHIALFALFARootVariance                            **2 &
         &                                +(                                                               &
         &                                  +self%outputAnalysisMolecularRatio_%ratioScatter(mass,node)    &
         &                                  *molecularFraction                                             &
         &                                  /(                                                             &
         &                                    +1.0                                                         &
         &                                    +molecularFraction                                           &
         &                                   )                                                             &
         &                                 )                                                           **2 &
         &                               )
    return
  end function randomErrorHIALFALFARootVariance
