!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018
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

  !% Implementation of a merger tree halo mass function sampling class optimized to minimize variance in the model stellar mass function.

  use Halo_Mass_Functions

  !# <mergerTreeHaloMassFunctionSampling name="mergerTreeHaloMassFunctionSamplingStllrMssFnctn" defaultThreadPrivate="yes">
  !#  <description>A merger tree halo mass function sampling class optimized to minimize variance in the model stellar mass function.</description>
  !# </mergerTreeHaloMassFunctionSampling>
  type, extends(mergerTreeHaloMassFunctionSamplingClass) :: mergerTreeHaloMassFunctionSamplingStllrMssFnctn
     !% Implementation of merger tree halo mass function sampling class optimized to minimize variance in the model stellar mass function.
     private
     class           (haloMassFunctionClass), pointer :: haloMassFunction_
     double precision                                 :: alpha             , beta               , &
          &                                              constant          , binWidthLogarithmic, &
          &                                              massMinimum       , massMaximum        , &
          &                                              massCharacteristic, normalization
   contains
     final     ::           stellarMassFunctionDestructor
     procedure :: sample => stellarMassFunctionSample
  end type mergerTreeHaloMassFunctionSamplingStllrMssFnctn

  interface mergerTreeHaloMassFunctionSamplingStllrMssFnctn
     !% Constructors for the {\normalfont \ttfamily stellarMassFunction} merger tree halo mass function sampling class.
     module procedure stellarMassFunctionConstructorParameters
     module procedure stellarMassFunctionConstructorInternal
  end interface mergerTreeHaloMassFunctionSamplingStllrMssFnctn

contains

  function stellarMassFunctionConstructorParameters(parameters) result(self)
    !% Constructor for the {\normalfont \ttfamily stellarMassFunction} merger tree halo mass function sampling class which builds the object from a parameter set.
    use Input_Parameters
    implicit none
    type            (mergerTreeHaloMassFunctionSamplingStllrMssFnctn)                :: self
    type            (inputParameters                                      ), intent(inout) :: parameters
    class           (haloMassFunctionClass                                ), pointer       :: haloMassFunction_
    double precision                                                                       :: alpha             , beta               , &
          &                                                                                   constant          , binWidthLogarithmic, &
          &                                                                                   massMinimum       , massMaximum        , &
          &                                                                                   massCharacteristic, normalization
    
    !# <inputParameter>
    !#   <name>normalization</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The value $\phi_0$ in a Schechter function, $\sigma(M) = \phi_0 (M/M_\star)^\alpha \exp(-[M/M_\star]^\beta)$, describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>alpha</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The value $\alpha$ in a Schechter function describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>beta</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The value $\beta$ in a Schechter function describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massCharacteristic</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The value $M_\star$ in a Schechter function describing the errors on the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>constant</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The constant error contribution to the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>binWidthLogarithmic</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The logarithmic width of bins in the stellar mass function to be assumed when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMinimum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum stellar mass to consider when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <inputParameter>
    !#   <name>massMaximum</name>
    !#   <cardinality>1</cardinality>
    !#   <description>The minimum stellar mass to consider when computing the optimal sampling density function for tree masses.</description>
    !#   <source>parameters</source>
    !#   <type>real</type>
    !# </inputParameter>
    !# <objectBuilder class="haloMassFunction" name="haloMassFunction_" source="parameters"/>
    self=mergerTreeHaloMassFunctionSamplingStllrMssFnctn(alpha,beta,constant,binWidthLogarithmic,massMinimum,massMaximum,massCharacteristic,normalization,haloMassFunction_)
    !# <inputParametersValidate source="parameters"/>
    return
  end function stellarMassFunctionConstructorParameters

  function stellarMassFunctionConstructorInternal(alpha,beta,constant,binWidthLogarithmic,massMinimum,massMaximum,massCharacteristic,normalization,haloMassFunction_) result(self)
    !% Internal constructor for the {\normalfont \ttfamily stellarMassFunction} merger tree halo mass function sampling class.
   implicit none
    type            (mergerTreeHaloMassFunctionSamplingStllrMssFnctn)                        :: self
    class           (haloMassFunctionClass                          ), intent(in   ), target :: haloMassFunction_
    double precision                                                 , intent(in   )         :: alpha             , beta               , &
         &                                                                                      constant          , binWidthLogarithmic, &
         &                                                                                      massMinimum       , massMaximum        , &
         &                                                                                      massCharacteristic, normalization
    !# <constructorAssign variables="alpha, beta, constant, binWidthLogarithmic, massMinimum, massMaximum, massCharacteristic, normalization, *haloMassFunction_"/>
    
    return
  end function stellarMassFunctionConstructorInternal

  subroutine stellarMassFunctionDestructor(self)
    !% Destructor for the {\normalfont \ttfamily stellarMassFunction} merger tree halo mass sampling class.
    implicit none
    type(mergerTreeHaloMassFunctionSamplingStllrMssFnctn), intent(inout) :: self

    !# <objectDestructor name="self%haloMassFunction_"/>
    return
  end subroutine stellarMassFunctionDestructor

  double precision function stellarMassFunctionSample(self,mass,time,massMinimum,massMaximum)
    !% Computes the halo mass function sampling rate optimized to minimize errors in the stellar mass function.
    use FGSL
    use Galacticus_Meta_Compute_Times
    use Numerical_Integration
    implicit none
    class           (mergerTreeHaloMassFunctionSamplingStllrMssFnctn), intent(inout) :: self
    double precision                                                 , intent(in   ) :: mass                               , massMaximum                 , &
         &                                                                              massMinimum                        , time
    double precision                                                 , parameter     :: toleranceAbsolute           =1.0d-3, toleranceRelative    =1.0d-2
    double precision                                                                 :: haloMassFunctionDifferential       , logStellarMassMaximum       , &
         &                                                                              logStellarMassMinimum              , treeComputeTime             , &
         &                                                                              xi                                 , xiIntegral                  , &
         &                                                                              massHalo
    type            (fgsl_function                                  )                :: integrandFunction
    type            (fgsl_integration_workspace                     )                :: integrationWorkspace
    !GCC$ attributes unused :: massMinimum, massMaximum

    ! Get the halo mass function, defined per logarithmic interval in halo mass.
    haloMassFunctionDifferential=+                                         mass  &
         &                       *self%haloMassFunction_%differential(time,mass)
    ! Compute the integral that appears in the "xi" function.
    massHalo             =mass
    logStellarMassMinimum=log10(self%massMinimum)
    logStellarMassMaximum=log10(self%massMaximum)
    xiIntegral           =Integrate(                                         &
         &                                            logStellarMassMinimum, &
         &                                            logStellarMassMaximum, &
         &                                            xiIntegrand          , &
         &                                            integrandFunction    , &
         &                                            integrationWorkspace , &
         &                          toleranceAbsolute=toleranceAbsolute    , &
         &                          toleranceRelative=toleranceRelative      &
         &                         )
    call Integrate_Done(integrandFunction,integrationWorkspace)
    ! Compute the "xi" function.
    xi              =+haloMassFunctionDifferential**2 &
         &           *xiIntegral
    ! Get the time taken to compute a tree of this mass.
    treeComputeTime =Galacticus_Time_Per_Tree(mass)
    ! Compute the optimal weighting for trees of this mass.
    stellarMassFunctionSample=sqrt(xi/treeComputeTime)
    return

  contains

    double precision function xiIntegrand(logStellarMass)
      !% The integrand appearing in the $\xi$ function.
      use Conditional_Mass_Functions
      implicit none
      double precision                              , intent(in   ) :: logStellarMass
      class           (conditionalMassFunctionClass), pointer       :: conditionalMassFunction_
      double precision                                              :: conditionalMassFunctionVariance , stellarMass       , &
           &                                                           stellarMassFunctionObservedError, stellarMassMaximum, &
           &                                                           stellarMassMinimum

      ! Compute the stellar mass and range corresponding to data bins.
      stellarMass       =10.0d0** logStellarMass
      stellarMassMinimum=10.0d0**(logStellarMass-0.5d0*self%binWidthLogarithmic)
      stellarMassMaximum=10.0d0**(logStellarMass+0.5d0*self%binWidthLogarithmic)
      ! Compute the variance in the model conditional stellar mass function.
      conditionalMassFunction_        => conditionalMassFunction()
      conditionalMassFunctionVariance =  conditionalMassFunction_%massFunctionVariance(massHalo,stellarMassMinimum,stellarMassMaximum)
      ! Compute the error in the observed stellar mass. We use a simple Schechter function (plus minimum error) fit.
      stellarMassFunctionObservedError=+self%normalization                                      &
           &                           *exp(-(stellarMass/self%massCharacteristic)**self%beta ) &
           &                           *     (stellarMass/self%massCharacteristic)**self%alpha  &
           &                           +self%constant

      ! Compute the integrand for the xi function integral.
      xiIntegrand=conditionalMassFunctionVariance/stellarMassFunctionObservedError**2
      return
    end function xiIntegrand

  end function stellarMassFunctionSample
