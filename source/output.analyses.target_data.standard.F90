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
  Implements the standard :galacticus-class:`outputAnalysisTargetDataClass` class --- a plain struct of axis-labelling and target-dataset fields, with all fields optional at construction.
  !!}

  !![
  <outputAnalysisTargetData name="outputAnalysisTargetDataStandard" docformat="rst">
   <description>
   The standard :galacticus-class:`outputAnalysisTargetDataClass` class --- a simple struct holding axis labels, log-scale flags, and target value/covariance arrays.  All fields are optional at construction; omitted axis labels default to ``'x'``/``'y'``, omitted target/log-scale flags default to empty/false, and the target arrays remain unallocated.
   </description>
  </outputAnalysisTargetData>
  !!]
  type, extends(outputAnalysisTargetDataClass) :: outputAnalysisTargetDataStandard
     !!{RST
     The standard :galacticus-class:`outputAnalysisTargetDataClass` class.
     !!}
     type            (varying_string)                              :: xAxisLabel
     type            (varying_string)                              :: yAxisLabel
     type            (varying_string)                              :: targetLabel
     logical                                                       :: xAxisIsLog
     logical                                                       :: yAxisIsLog
     double precision                , allocatable, dimension(:  ) :: valueTarget     , covarianceTarget1D
     double precision                , allocatable, dimension(:,:) :: covarianceTarget
   contains
     procedure :: hasTarget => standardHasTarget
  end type outputAnalysisTargetDataStandard

  interface outputAnalysisTargetDataStandard
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisTargetDataStandard` output analysis target data class.
     !!}
     module procedure standardConstructorParameters
     module procedure standardConstructorInternal
  end interface outputAnalysisTargetDataStandard

contains

  function standardConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisTargetDataStandard` output analysis target data class which takes a parameter set as input.  Each field has its own ``<inputParameter>``; absent parameters fall back to the same defaults the inline constructor uses.
    !!}
    use :: Error           , only : Error_Report
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisTargetDataStandard)                                :: self
    type            (inputParameters                 ), intent(inout)                 :: parameters
    type            (varying_string                  )                                :: xAxisLabel      , yAxisLabel        , &
         &                                                                               targetLabel
    logical                                                                           :: xAxisIsLog      , yAxisIsLog
    double precision                                  , allocatable  , dimension(:  ) :: valueTarget     , covarianceTarget1D
    double precision                                  , allocatable  , dimension(:,:) :: covarianceTarget

    !![
    <inputParameter docformat="rst">
      <name>xAxisLabel</name>
      <source>parameters</source>
      <variable>xAxisLabel</variable>
      <description>
      Axis label for the property (x-axis) of the output analysis.
      </description>
      <defaultValue>var_str('x')</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yAxisLabel</name>
      <source>parameters</source>
      <variable>yAxisLabel</variable>
      <description>
      Axis label for the function value (y-axis) of the output analysis.
      </description>
      <defaultValue>var_str('y')</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>targetLabel</name>
      <source>parameters</source>
      <variable>targetLabel</variable>
      <description>
      Label identifying the comparison/target dataset, if any.
      </description>
      <defaultValue>var_str('')</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>xAxisIsLog</name>
      <source>parameters</source>
      <variable>xAxisIsLog</variable>
      <description>
      Whether the x-axis should be displayed on a logarithmic scale.
      </description>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yAxisIsLog</name>
      <source>parameters</source>
      <variable>yAxisIsLog</variable>
      <description>
      Whether the y-axis should be displayed on a logarithmic scale.
      </description>
      <defaultValue>.false.</defaultValue>
    </inputParameter>
    !!]
    if (parameters%isPresent('valueTarget')) then
       allocate(valueTarget(parameters%count('valueTarget')))
       !![
       <inputParameter docformat="rst">
         <name>valueTarget</name>
         <source>parameters</source>
         <variable>valueTarget</variable>
         <description>
         Target dataset values to compare against, one per bin of the output analysis.
         </description>
       </inputParameter>
       !!]
       if (parameters%isPresent('covarianceTarget')) then
          if (parameters%count('covarianceTarget') /= parameters%count('valueTarget')**2) call Error_Report('covariance size does not match value size'//{introspection:location})
          allocate(covarianceTarget1D(parameters%count('valueTarget')*parameters%count('valueTarget')))
          allocate(covarianceTarget  (parameters%count('valueTarget'),parameters%count('valueTarget')))
          !![
          <inputParameter docformat="rst">
          <name>covarianceTarget</name>
          <source>parameters</source>
          <variable>covarianceTarget1D</variable>
          <description>
          Target-dataset covariance matrix corresponding to the ``valueTarget`` array.
          </description>
          </inputParameter>
          !!]
          covarianceTarget=reshape(covarianceTarget1D,[parameters%count('valueTarget'),parameters%count('valueTarget')])
       end if
    else
        if (parameters%isPresent('covarianceTarget')) call Error_Report('variance provided but no target values'//{introspection:location})
    end if
    self=outputAnalysisTargetDataStandard(xAxisLabel,yAxisLabel,targetLabel,xAxisIsLog,yAxisIsLog,valueTarget,covarianceTarget)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function standardConstructorParameters

  function standardConstructorInternal(xAxisLabel,yAxisLabel,targetLabel,xAxisIsLog,yAxisIsLog,valueTarget,covarianceTarget) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`outputAnalysisTargetDataStandard` output analysis target data class. All arguments are optional --- omitted axis labels default to ``'x'`` / ``'y'``, omitted log-scale flags default to ``.false.``, and omitted target arrays simply remain unallocated.  These defaults match the per-argument defaults used by the outer 1D function output-analysis constructors when ``targetData_`` itself is absent.
    !!}
    implicit none
    type            (outputAnalysisTargetDataStandard)                                          :: self
    type            (varying_string                  ), intent(in   ), optional                 :: xAxisLabel      , yAxisLabel, &
         &                                                                                         targetLabel
    logical                                           , intent(in   ), optional                 :: xAxisIsLog      , yAxisIsLog
    double precision                                  , intent(in   ), optional, dimension(:  ) :: valueTarget
    double precision                                  , intent(in   ), optional, dimension(:,:) :: covarianceTarget
    !![
    <constructorAssign variables="xAxisLabel='x', yAxisLabel='y', targetLabel, xAxisIsLog=.false., yAxisIsLog=.false., valueTarget, covarianceTarget"/>
    !!]
    ! Maintain a 1D flattening of the covariance for the auto-built descriptor, which serializes
    ! parameter-file fields as 1D lists.  Populated lazily here so every construction path (both
    ! Parameters-from-XML and Internal-from-Fortran) ends up with a consistent `self%covarianceTarget1D`.
    if (allocated(self%covarianceTarget)) &
         & self%covarianceTarget1D=reshape(self%covarianceTarget,[size(self%covarianceTarget)])
    return
  end function standardConstructorInternal

  logical function standardHasTarget(self)
    !!{RST
    Return whether both target arrays are populated.
    !!}
    implicit none
    class(outputAnalysisTargetDataStandard), intent(inout) :: self

    standardHasTarget=allocated(self%valueTarget).and.allocated(self%covarianceTarget)
    return
  end function standardHasTarget
