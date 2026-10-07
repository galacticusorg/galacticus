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

!+    Contributions to this file made by: Andrew Benson, Claude.

  !!{RST
  Implements a synthetic output analysis whose result is a linear function of its parameters.
  !!}

  !![
  <outputAnalysis name="outputAnalysisSyntheticLinear" docformat="rst">
   <description>
   A synthetic output analysis, intended for testing, whose result does not depend on the galaxies of the model. The analysis
   is the straight line :math:`y(x) = a + b x`, where :math:`a` is ``[intercept]`` and :math:`b` is ``[slope]``, evaluated at
   the abscissae ``[x]``. It is compared with target data ``[yTarget]``, with independent Gaussian uncertainties of variance
   ``[varianceTarget]``, so that the log-likelihood,

   .. math::

      \log \mathcal{L} = -{1 \over 2} \sum_i \left[ {(y_i - y(x_i))^2 \over \sigma_i^2} + \log (2 \pi \sigma_i^2) \right],

   is an exactly Gaussian function of :math:`a` and :math:`b`. The results are written in the same form as those of
   function analyses (such as :galacticus-class:`outputAnalysisMeanFunction1D`), so this analysis can stand in for a model
   when testing tools which consume analyses, such as the training of emulators (see :ref:`manual-sec-Emulation`).
   </description>
  </outputAnalysis>
  !!]
  type, extends(outputAnalysisClass) :: outputAnalysisSyntheticLinear
     !!{RST
     A synthetic output analysis whose result is a linear function of its parameters.
     !!}
     private
     type            (varying_string)                            :: label         , comment
     double precision                , allocatable, dimension(:) :: x             , yTarget, &
          &                                                         varianceTarget
     double precision                                            :: intercept     , slope
   contains
     procedure :: analyze       => syntheticLinearAnalyze
     procedure :: finalize      => syntheticLinearFinalize
     procedure :: reduce        => syntheticLinearReduce
     procedure :: logLikelihood => syntheticLinearLogLikelihood
  end type outputAnalysisSyntheticLinear

  interface outputAnalysisSyntheticLinear
     !!{RST
     Constructors for the :galacticus-class:`outputAnalysisSyntheticLinear` output analysis class.
     !!}
     module procedure syntheticLinearConstructorParameters
     module procedure syntheticLinearConstructorInternal
  end interface outputAnalysisSyntheticLinear

contains

  function syntheticLinearConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`outputAnalysisSyntheticLinear` output analysis class which takes a parameter set as
    input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (outputAnalysisSyntheticLinear)                              :: self
    type            (inputParameters              ), intent(inout)               :: parameters
    type            (varying_string               )                              :: label         , comment
    double precision                               , allocatable  , dimension(:) :: x             , yTarget, &
         &                                                                          varianceTarget
    double precision                                                             :: intercept     , slope

    allocate(x             (parameters%count('x')))
    allocate(yTarget       (parameters%count('x')))
    allocate(varianceTarget(parameters%count('x')))
    !![
    <inputParameter docformat="rst">
      <name>label</name>
      <source>parameters</source>
      <defaultValue>var_str('syntheticLinear')</defaultValue>
      <description>
      A label for the analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>comment</name>
      <source>parameters</source>
      <defaultValue>var_str('a synthetic linear function')</defaultValue>
      <description>
      A descriptive comment for the analysis.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>x</name>
      <source>parameters</source>
      <description>
      The abscissae at which the linear function is evaluated.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>yTarget</name>
      <source>parameters</source>
      <description>
      The target data at each abscissa.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>varianceTarget</name>
      <source>parameters</source>
      <description>
      The variance of the target data at each abscissa.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>intercept</name>
      <source>parameters</source>
      <description>
      The intercept, :math:`a`, of the linear function.
      </description>
    </inputParameter>
    <inputParameter docformat="rst">
      <name>slope</name>
      <source>parameters</source>
      <description>
      The slope, :math:`b`, of the linear function.
      </description>
    </inputParameter>
    !!]
    self=outputAnalysisSyntheticLinear(label,comment,x,yTarget,varianceTarget,intercept,slope)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function syntheticLinearConstructorParameters

  function syntheticLinearConstructorInternal(label,comment,x,yTarget,varianceTarget,intercept,slope) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`outputAnalysisSyntheticLinear` output analysis class.
    !!}
    use :: Error, only : Error_Report
    implicit none
    type            (outputAnalysisSyntheticLinear)                              :: self
    type            (varying_string               ), intent(in   )               :: label         , comment
    double precision                               , intent(in   ), dimension(:) :: x             , yTarget, &
         &                                                                          varianceTarget
    double precision                               , intent(in   )               :: intercept     , slope
    !![
    <constructorAssign variables="label, comment, x, yTarget, varianceTarget, intercept, slope"/>
    !!]

    if (size(yTarget) /= size(x) .or. size(varianceTarget) /= size(x)) call Error_Report('`x`, `yTarget`, and `varianceTarget` must have the same size'//{introspection:location})
    if (any(varianceTarget <= 0.0d0)) call Error_Report('`varianceTarget` must be positive'//{introspection:location})
    return
  end function syntheticLinearConstructorInternal

  subroutine syntheticLinearAnalyze(self,node,iOutput)
    !!{RST
    Analyze a node - nothing to do, as this analysis does not depend on the galaxies of the model.
    !!}
    implicit none
    class  (outputAnalysisSyntheticLinear), intent(inout) :: self
    type   (treeNode                     ), intent(inout) :: node
    integer(c_size_t                     ), intent(in   ) :: iOutput
    !$GLC attributes unused :: self, node, iOutput

    return
  end subroutine syntheticLinearAnalyze

  subroutine syntheticLinearReduce(self,reduced)
    !!{RST
    Reduce the analysis - nothing to do, as this analysis accumulates nothing.
    !!}
    implicit none
    class(outputAnalysisSyntheticLinear), intent(inout) :: self
    class(outputAnalysisClass          ), intent(inout) :: reduced
    !$GLC attributes unused :: self, reduced

    return
  end subroutine syntheticLinearReduce

  subroutine syntheticLinearFinalize(self,groupName)
    !!{RST
    Write the analysis to the output file.
    !!}
    use :: Output_HDF5, only : outputFile
    use :: HDF5_Access, only : hdf5Access
    use :: IO_HDF5    , only : hdf5Group
    implicit none
    class           (outputAnalysisSyntheticLinear), intent(inout)                 :: self
    type            (varying_string               ), intent(in   ), optional       :: groupName
    double precision                               , allocatable  , dimension(:,:) :: covariance
    integer                                                                        :: i

    ! The model has no uncertainty; the target data have only diagonal covariance.
    allocate(covariance(size(self%x),size(self%x)))
    covariance=0.0d0
    !$ call hdf5Access%set()
    block
      type(hdf5Group) :: analysesGroup, subGroup, analysisGroup
      analysesGroup=outputFile%openGroup('analyses')
      if (present(groupName)) then
         subGroup     =analysesGroup%openGroup(char(groupName ))
         analysisGroup=subGroup     %openGroup(char(self%label),char(self%comment))
      else
         analysisGroup=analysesGroup%openGroup(char(self%label),char(self%comment))
      end if
      call analysisGroup%writeAttribute(char(self%comment)  ,'description'      )
      call analysisGroup%writeAttribute('function1D'        ,'type'             )
      call analysisGroup%writeAttribute('x'                 ,'xDataset'         )
      call analysisGroup%writeAttribute('y'                 ,'yDataset'         )
      call analysisGroup%writeAttribute('yCovariance'       ,'yCovariance'      )
      call analysisGroup%writeAttribute('yTarget'           ,'yDatasetTarget'   )
      call analysisGroup%writeAttribute('yCovarianceTarget' ,'yCovarianceTarget')
      call analysisGroup%writeAttribute(self%logLikelihood(),'logLikelihood'    )
      call analysisGroup%writeDataset  (self%x                      ,'x'          ,'The abscissae.'                      )
      call analysisGroup%writeDataset  (self%intercept+self%slope*self%x,'y'      ,'The linear function.'                )
      call analysisGroup%writeDataset  (covariance                  ,'yCovariance','The covariance of the linear function.')
      call analysisGroup%writeDataset  (self%yTarget                ,'yTarget'    ,'The target data.'                    )
      do i=1,size(self%x)
         covariance(i,i)=self%varianceTarget(i)
      end do
      call analysisGroup%writeDataset  (covariance                  ,'yCovarianceTarget','The covariance of the target data.')
    end block
    !$ call hdf5Access%unset()
    return
  end subroutine syntheticLinearFinalize

  double precision function syntheticLinearLogLikelihood(self)
    !!{RST
    Return the log-likelihood of the target data given the linear function.
    !!}
    use :: Numerical_Constants_Math, only : Pi
    implicit none
    class(outputAnalysisSyntheticLinear), intent(inout) :: self

    syntheticLinearLogLikelihood=-0.5d0                                                                 &
         &                       *sum(                                                                  &
         &                            +(self%yTarget-self%intercept-self%slope*self%x)**2/self%varianceTarget &
         &                            +log(2.0d0*Pi*self%varianceTarget)                                &
         &                           )
    return
  end function syntheticLinearLogLikelihood
