!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022
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
  Contains a module which implements a transfer function class for fuzzy dark matter using the fitting function of
  \cite{murgia_non-cold_2017}.
  !!}

  !![
  <transferFunction name="transferFunctionFuzzyDM">
   <description>A transfer function class for fuzzy dark matter using the fitting function of \cite{murgia_non-cold_2017}.</description>
  </transferFunction>
  !!]
  type, extends(transferFunctionMurgia2017) :: transferFunctionFuzzyDM
     !!{
     A transfer function class for fuzzy dark matter using the fitting function of \cite{murgia_non-cold_2017}.
     !!}
     private
     double precision :: m22
  end type transferFunctionFuzzyDM
   
  interface transferFunctionFuzzyDM
     !!{
     Constructors for the {\normalfont \ttfamily fuzzyDM} transfer function class.
     !!}
     module procedure fuzzyDMConstructorParameters
     module procedure fuzzyDMConstructorInternal
  end interface transferFunctionFuzzyDM

contains

  function fuzzyDMConstructorParameters(parameters) result(self)
    !!{
    Constructor for the {\normalfont \ttfamily fuzzyDM} transfer function class which takes a parameter set as input.
    !!}
    use :: Cosmology_Functions           , only : cosmologyFunctions        , cosmologyFunctionsClass
    use :: Cosmology_Functions_Parameters, only : requestTypeExpansionFactor
    use :: Input_Parameters              , only : inputParameter            , inputParameters
    implicit none
    type            (transferFunctionFuzzyDM )                :: self
    type            (inputParameters         ), intent(inout) :: parameters
    class           (transferFunctionClass   ), pointer       :: transferFunctionCDM
    class           (cosmologyParametersClass), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), pointer       :: cosmologyFunctions_
    double precision                                          :: m22                 , beta    , &
         &                                                       gamma               , redshift

    ! Validate parameters.
    if (.not.parameters%isPresent('transferFunction')) call Galacticus_Error_Report("an explicit 'transferFunction' must be given"//{introspection:location})
    ! Read parameters.
    !![
    <inputParameter>
      <name>m22</name>
      <source>parameters</source>
      <description>The mass of the fuzzy dark matter particle in units of $10^{-22}$~eV, $m_{22}$.</description>
    </inputParameter>
    <inputParameter>
      <name>beta</name>
      <source>parameters</source>
      <defaultValue>5.475d0</defaultValue>
      <defaultSource>\citep[][average of values in Table~4]{murgia_non-cold_2017}</defaultSource>
      <description>The parameter $\beta$, which controls the shape of the cut-off, appearing in the transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <inputParameter>
      <name>gamma</name>
      <source>parameters</source>
      <defaultValue>-2.0d0</defaultValue>
      <defaultSource>\citep[][average of values in Table~4]{murgia_non-cold_2017}</defaultSource>
      <description>The parameter $\gamma$, which controls the shape of the cut-off, appearing in the transfer function \citep{murgia_non-cold_2017}.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="transferFunction"    name="transferFunctionCDM"  source="parameters"/>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>cosmologyFunctions_%redshiftFromExpansionFactor(cosmologyFunctions_%equalityEpochMatterRadiation(requestTypeExpansionFactor))</defaultValue>
      <description>The redshift of the epoch at which the transfer function is defined.</description>
    </inputParameter>
    !!]
    self=transferFunctionFuzzyDM(transferFunctionCDM,m22,beta,gamma,cosmologyFunctions_%cosmicTime(cosmologyFunctions_%expansionFactorFromRedshift(redshift)),cosmologyParameters_,cosmologyFunctions_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="transferFunctionCDM" />
    !!]
    return
  end function fuzzyDMConstructorParameters

  function fuzzyDMConstructorInternal(transferFunctionCDM,m22,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_) result(self)
    !!{
    Internal constructor for the {\normalfont \ttfamily fuzzyDM} transfer function class.
    !!}
    implicit none
    type            (transferFunctionFuzzyDM)                         :: self
    class           (transferFunctionClass   ), target, intent(in   ) :: transferFunctionCDM
    double precision                                  , intent(in   ) :: m22                 , beta, &
         &                                                               gamma               , time
    class           (cosmologyParametersClass), target, intent(in   ) :: cosmologyParameters_
    class           (cosmologyFunctionsClass ), target, intent(in   ) :: cosmologyFunctions_
    double precision                                                  :: alpha
    !![
    <constructorAssign variables="m22"/>
    !!]

    alpha  =+0.11d0                   &
         &  *self%m22**(-4.0d0/9.0d0) 
    self%transferFunctionMurgia2017=transferFunctionMurgia2017(transferFunctionCDM,alpha,beta,gamma,time,cosmologyParameters_,cosmologyFunctions_)
    return
  end function fuzzyDMConstructorInternal

