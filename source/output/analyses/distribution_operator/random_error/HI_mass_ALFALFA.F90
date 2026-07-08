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

  !!{RST
  Implements a random error output analysis distribution operator class providing errors in HI mass for the ALFALFA survey.
  !!}

  use :: Output_Analysis_Molecular_Ratios, only : outputAnalysisMolecularRatioClass

  !![
  <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomErrorALFLF" docformat="rst">
   <description>
   A random error output analysis distribution operator class providing errors in HI mass for the ALFALFA survey. To account for both observational errors and scatter in :math:`R_\mathrm{mol}`, the HI mass of each galaxy is modeled as a Gaussian in :math:`\log_{10}M_\mathrm{HI}` when constructing the mass function. Observational random errors on HI mass, including those arising from flux density uncertainties and errors in the assumed distance to each source, are taken from Fig. 19 of :cite:t:`haynes_arecibo_2011`. The magnitude of the error as a function of HI mass is fit using a functional form:

   .. math::

      \sigma_\mathrm{obs} = a + \exp\left(-{\log_{10}(M_\mathrm{HI}/\mathrm{M}_\odot)-b\over c}\right),

   where :math:`\sigma_\mathrm{obs}` is the error on :math:`\log_{10}(M_\mathrm{HI}/\mathrm{M}_\odot)`. We find a reasonable fit using values\footnote{This should not be regarded as a formal good fit. Error estimates are approximate---we have simply found a functional form that roughly describes them, along with conservative errors on the parameters of this function which are included in the priors.} of :math:`a=`\ ``a``\ :math:`=0.100 \pm 0.010`, :math:`b=`\ ``b``\ :math:`=5.885 \pm 0.100`, and :math:`c=`\ ``c``\ :math:`=0.505 \pm 0.020` as shown in Fig. :numref:`{number} &lt;fig-ALFALFAErrorModel&gt;`. The total random error on the logarithm of each galaxy mass is given by :math:`\sigma^2 = \sigma_{R_\mathrm{mol}}^2+\sigma_\mathrm{obs}^2`, and is used as the width of the Gaussian kernel when applying each galaxy to the mass function histogram (as described above).

   .. figure:: Plots/DataAnalysis/alfalfaHIMassErrorModel.pdf
      :name: fig-ALFALFAErrorModel

      The observational random error in galaxy HI mass as a function of HI mass for the ALFALFA survey. Points show the errors reported by :cite:t:`haynes_arecibo_2011`, while the line shows a simple functional form fit to these errors.
   </description>
  </outputAnalysisDistributionOperator>
  !!]
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRandomErrorALFLF
     !!{RST
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
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisDistributionOperatorRandomErrorALFLF` output analysis distribution operator class.
     !!}
     module procedure randomErrorHIALFALFAConstructorParameters
     module procedure randomErrorHIALFALFAConstructorInternal
  end interface outputAnalysisDistributionOperatorRandomErrorALFLF

contains

  function randomErrorHIALFALFAConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisDistributionOperatorRandomErrorALFLF` output analysis distribution operator class which takes a parameter set as input.
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
    <inputParameter docformat="rst">
      <name>a</name>
      <source>parameters</source>
      <defaultValue>0.1d0</defaultValue>
      <variable>a</variable>
      <description>
      Parameter :math:`a` in the ALFALFA HI mass error model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>b</name>
      <source>parameters</source>
      <defaultValue>5.885d0</defaultValue>
      <variable>b</variable>
      <description>
      Parameter :math:`b` in the ALFALFA HI mass error model.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>c</name>
      <source>parameters</source>
      <defaultValue>0.505d0</defaultValue>
      <variable>c</variable>
      <description>
      Parameter :math:`c` in the ALFALFA HI mass error model.
      </description>
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
    !!{RST
    Internal constructor for the :galacticus-class:`outputAnalysisDistributionOperatorRandomErrorALFLF` output analysis distribution operator class.
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
    !!{RST
    Destructor for the :galacticus-class:`outputAnalysisDistributionOperatorRandomErrorALFLF` output analysis distribution operator class.
    !!}
    implicit none
    type(outputAnalysisDistributionOperatorRandomErrorALFLF), intent(inout) :: self

    !![
    <objectDestructor name="self%outputAnalysisMolecularRatio_"/>
    !!]
    return
  end subroutine randomErrorHIALFALFADestructor

  double precision function randomErrorHIALFALFARootVariance(self,propertyValue,node)
    !!{RST
    Computes errors on :math:`\log_{10}(`\ HI masses\ :math:`)` for the ALFALFA survey analysis. Uses a simple fitting function. See ``constraints/dataAnalysis/hiMassFunction_ALFALFA_z0.00/alfalfaHIMassErrorModel.py`` for details.
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
