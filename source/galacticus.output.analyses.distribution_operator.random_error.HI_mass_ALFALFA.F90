!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019
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

  !% Contains a module which implements a random error output analysis distribution operator class providing errors in HI mass for
  !% the ALFALFA survey.

  use :: Output_Analysis_Molecular_Ratios, only : outputAnalysisMolecularRatioClass

  !# <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRandomErrorALFLF">
  !#  <description>A random error output analysis distribution operator class providing errors in HI mass for the ALFALFA survey. Specifically, $\sigma_\mathrm{obs} = a + \exp\left(-{\log_{10}(M_\mathrm{HI}/M_\odot)-b\over c}\right)$.</description>
  !# </outputAnalysisDistributionOperator>
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRandomErrorALFLF
     !% A random error output distribution operator class providing errors in HI mass for the ALFALFA survey.
     private
     double precision                                             :: a                            , b, &
          &                                                          c
     class           (outputAnalysisMolecularRatioClass), pointer :: outputAnalysisMolecularRatio_ => null()
   contains
     final     ::                 randomErrorHIALFALFADestructor
     procedure :: rootVariance => randomErrorHIALFALFARootVariance
  end type outputAnalysisDistributionOperatorRandomErrorALFLF

  interface outputAnalysisDistributionOperatorRandomErrorALFLF
     !% Constructors for the ``randomErrorHIALFALFA'' output analysis distribution operator class.
     module procedure randomErrorHIALFALFAConstructorParameters
     module procedure randomErrorHIALFALFAConstructorInternal
  end interface outputAnalysisDistributionOperatorRandomErrorALFLF

contains

  function randomErrorHIALFALFAConstructorParameters(parameters) result(self)
    !% Constructor for the ``randomErrorHIALFALFA'' output analysis distribution operator class which takes a parameter set as input.
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorALFLF)                :: self
    type            (inputParameters                                   ), intent(inout) :: parameters
    class           (outputAnalysisMolecularRatioClass                 ), pointer       :: outputAnalysisMolecularRatio_
    double precision                                                                    :: a                            , b, &
         &                                                                                 c

    ! Check and read parameters.
    !# <inputParameter>
    !#   <name>a</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.1d0</defaultValue>
    !#   <variable>a</variable>
    !#   <description>Parameter $a$ in the ALFALFA HI mass error model.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>b</name>
    !#   <source>parameters</source>
    !#   <defaultValue>5.885d0</defaultValue>
    !#   <variable>b</variable>
    !#   <description>Parameter $b$ in the ALFALFA HI mass error model.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>c</name>
    !#   <source>parameters</source>
    !#   <defaultValue>0.505d0</defaultValue>
    !#   <variable>c</variable>
    !#   <description>Parameter $c$ in the ALFALFA HI mass error model.</description>
    !#   <type>float</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisMolecularRatio" name="outputAnalysisMolecularRatio_" source="parameters"/>
    ! Construct the object.
    self=outputAnalysisDistributionOperatorRandomErrorALFLF(a,b,c,outputAnalysisMolecularRatio_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="outputAnalysisMolecularRatio_"/>
    return
  end function randomErrorHIALFALFAConstructorParameters

  function randomErrorHIALFALFAConstructorInternal(a,b,c,outputAnalysisMolecularRatio_) result(self)
    !% Internal constructor for the ``randomErrorHIALFALFA'' output analysis distribution operator class.
    implicit none
    type            (outputAnalysisDistributionOperatorRandomErrorALFLF)                        :: self
    double precision                                                    , intent(in   )         :: a                            , b, &
         &                                                                                         c
    class           (outputAnalysisMolecularRatioClass                 ), intent(in   ), target :: outputAnalysisMolecularRatio_
    !# <constructorAssign variables="a, b, c, *outputAnalysisMolecularRatio_"/>

    return
  end function randomErrorHIALFALFAConstructorInternal

  subroutine randomErrorHIALFALFADestructor(self)
    !% Destructor for the ``randomErrorHIALFALFA'' output analysis distribution operator class.
    implicit none
    type(outputAnalysisDistributionOperatorRandomErrorALFLF), intent(inout) :: self

    !# <objectDestructor name="self%outputAnalysisMolecularRatio_"/>
    return
  end subroutine randomErrorHIALFALFADestructor

  double precision function randomErrorHIALFALFARootVariance(self,propertyValue,node)
    !% Computes errors on $\log_{10}($HI masses$)$ for the ALFALFA survey analysis. Uses a simple fitting function. See
    !% {\normalfont \ttfamily constraints/dataAnalysis/hiMassFunction\_ALFALFA\_z0.00/alfalfaHIMassErrorModel.pl} for details.
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
