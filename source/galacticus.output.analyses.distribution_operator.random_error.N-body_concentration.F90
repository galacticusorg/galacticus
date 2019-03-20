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

  !% Contains a module which implements a random error output analysis distribution operator class providing errors in $\log_{10}$
  !% of N-body halo concentration.

  use Output_Analysis_Property_Extractions
  
  !# <outputAnalysisDistributionOperator name="outputAnalysisDistributionOperatorRndmErrNbdyCnc">
  !#  <description>A random error output analysis distribution operator class providing errors in $\log_{10}$ of N-body halo concentration.</description>
  !# </outputAnalysisDistributionOperator>
  type, extends(outputAnalysisDistributionOperatorRandomError) :: outputAnalysisDistributionOperatorRndmErrNbdyCnc
     !% A random error output distribution operator class providing errors in $\log_{10}$ of N-body halo mass.
     private
     class           (outputAnalysisPropertyExtractorClass), pointer                   :: outputAnalysisPropertyExtractor_ => null()
     double precision                                      , allocatable, dimension(:) :: a
     double precision                                                                  :: b                               , massParticle
   contains
     final     ::                 randomErrorNbdyCncDestructor
     procedure :: rootVariance => randomErrorNbdyCncRootVariance
  end type outputAnalysisDistributionOperatorRndmErrNbdyCnc

  interface outputAnalysisDistributionOperatorRndmErrNbdyCnc
     !% Constructors for the ``randomErrorNbdyCnc'' output analysis distribution operator class.
     module procedure randomErrorNbdyCncConstructorParameters
     module procedure randomErrorNbdyCncConstructorInternal
  end interface outputAnalysisDistributionOperatorRndmErrNbdyCnc

contains

  function randomErrorNbdyCncConstructorParameters(parameters) result(self)
    !% Constructor for the ``randomErrorNbdyCnc'' output analysis distribution operator class which takes a parameter set as input.
    use Input_Parameters
    implicit none
    type            (outputAnalysisDistributionOperatorRndmErrNbdyCnc)                              :: self
    type            (inputParameters                                 ), intent(inout)               :: parameters
    class           (outputAnalysisPropertyExtractorClass            ), pointer                     :: outputAnalysisPropertyExtractor_
    double precision                                                  , allocatable  , dimension(:) :: a
    double precision                                                                                :: b                               , massParticle

    allocate(a(parameters%count('a')))
    !# <inputParameter>
    !#   <name>a</name>
    !#   <source>parameters</source>
    !#   <description>Coefficients of the polynomial in concentration in the concentration error model.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>b</name>
    !#   <source>parameters</source>
    !#   <description>The exponent of particle number in the concentration error model.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massParticle</name>
    !#   <source>parameters</source>
    !#   <description>The mass of the particle in the N-body simulation.</description>
    !#   <type>real</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    !# <objectBuilder class="outputAnalysisPropertyExtractor" name="outputAnalysisPropertyExtractor_" source="parameters"/>
    self=outputAnalysisDistributionOperatorRndmErrNbdyCnc(a,b,massParticle,outputAnalysisPropertyExtractor_)
    !# <inputParametersValidate source="parameters"/>
    !# <objectDestructor name="outputAnalysisPropertyExtractor_"/>
    return
  end function randomErrorNbdyCncConstructorParameters

  function randomErrorNbdyCncConstructorInternal(a,b,massParticle,outputAnalysisPropertyExtractor_) result(self)
    !% Internal constructor for the ``randomErrorNbdyCnc'' output analysis distribution operator class.
    implicit none
    type            (outputAnalysisDistributionOperatorRndmErrNbdyCnc)                              :: self
    double precision                                                  , intent(in   ), dimension(:) :: a
    double precision                                                  , intent(in   )               :: b                               , massParticle
    class           (outputAnalysisPropertyExtractorClass            ), intent(in   ), target       :: outputAnalysisPropertyExtractor_
    !# <constructorAssign variables="a, b, massParticle, *outputAnalysisPropertyExtractor_"/>

    return
  end function randomErrorNbdyCncConstructorInternal

  subroutine randomErrorNbdyCncDestructor(self)
    !% Destructor for  the ``nbodyConcentration'' output analysis distribution operator class.
    type(outputAnalysisDistributionOperatorRndmErrNbdyCnc), intent(inout) :: self
    
    !# <objectDestructor name="self%outputAnalysisPropertyExtractor_" />
    return
  end subroutine randomErrorNbdyCncDestructor

  double precision function randomErrorNbdyCncRootVariance(self,propertyValue,node)
    !% Computes errors on $\log_{10}($halo concentration$)$ for N-body halos.
    implicit none
    class           (outputAnalysisDistributionOperatorRndmErrNbdyCnc), intent(inout) :: self
    double precision                                                  , intent(in   ) :: propertyValue
    type            (treeNode                                        ), intent(inout) :: node
    double precision                                                                  :: nbodyMassPropertyValue
    integer                                                                           :: nbodyMassPropertyType , i

    nbodyMassPropertyType         =+self%outputAnalysisPropertyExtractor_%type   (    )
    nbodyMassPropertyValue        =+self%outputAnalysisPropertyExtractor_%extract(node)
    randomErrorNbdyCncRootVariance=+self%b                        &
         &                         *log10(                        &
         &                                +nbodyMassPropertyValue &
         &                                /self%massParticle      &
         &                               )
    do i=1,size(self%a)
       randomErrorNbdyCncRootVariance=randomErrorNbdyCncRootVariance+self%a(i)*propertyValue**(i-1)
    end do
    randomErrorNbdyCncRootVariance=10.0d0**randomErrorNbdyCncRootVariance    
    return
  end function randomErrorNbdyCncRootVariance
