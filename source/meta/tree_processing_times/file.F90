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
Implements a merger tree processing time estimator using a polynomial relation read from file.
!!}

  !![
  <metaTreeProcessingTime name="metaTreeProcessingTimeFile" docformat="rst">
   <description>
   A merger tree processing time class which estimates processing times using a polynomial relation read from file. Specifically, the time taken to process a tree is estimate to be

   .. math::

      \log_{10} [ \tau_\mathrm{tree}(M)] = \sum_{i=0}^2 C_i (\log_{10} M)^i,

   where :math:`M` is the root mass of the tree and the coefficients :math:`C_i` are read from a file, the name of which is specified via the ``[fileName]`` parameter. This file should be an XML document with the structure:

   .. code-block:: none

      &lt;timing&gt;
       &lt;fit&gt;
         &lt;coefficient&gt;-0.73&lt;/coefficient&gt;
         &lt;coefficient&gt;-0.20&lt;/coefficient&gt;
         &lt;coefficient&gt;0.03&lt;/coefficient&gt;
       &lt;/fit&gt;
      &lt;/timing&gt;

   where the array of coefficients give the values :math:`C_0`, :math:`C_1` and :math:`C_2`.

   Alternatively, if the ``[fileName]`` has an ``.hdf5`` or ``.h5`` extension it is assumed to be a Galacticus output file from a previous run of the same (or a similar) model in which the :galacticus-class:`mergerTreeOperatorTreeProcessingTimer` operator was active. In that case the fit coefficients are read directly from the ``metaData/treeTiming/fitCoefficientMass`` dataset, so that ``&lt;metaTreeProcessingTime value="file"&gt;&lt;fileName value="previousRun.hdf5"/&gt;&lt;/metaTreeProcessingTime&gt;`` &quot;just works&quot; without any manual fitting step.
   </description>
  </metaTreeProcessingTime>
  !!]
  type, extends(metaTreeProcessingTimeClass) :: metaTreeProcessingTimeFile
     !!{RST
     A merger tree processing time estimator using a polynomial relation read from file.
     !!}
     private
     double precision                , dimension(0:2) :: fitCoefficient          , fitCoefficientCountNodes
     logical                                          :: haveCountNodesFit=.false.
     type            (varying_string)                 :: fileName
   contains
     procedure :: time             => fileTime
     procedure :: timeByCountNodes => fileTimeByCountNodes
  end type metaTreeProcessingTimeFile

  interface metaTreeProcessingTimeFile
     !!{RST
     Constructors for the :galacticus-class:`metaTreeProcessingTimeFile` merger tree processing time estimator.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface metaTreeProcessingTimeFile

contains

  function fileConstructorParameters(parameters) result(self)
    !!{RST
    Constructor for the :galacticus-class:`metaTreeProcessingTimeFile` merger tree processing time estimator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(metaTreeProcessingTimeFile)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    type(varying_string            )                :: fileName

    !![
    <inputParameter docformat="rst">
      <name>fileName</name>
      <description>
      The name of the file which contains fit coefficients for the time per tree fitting function.
      </description>
      <source>parameters</source>
    </inputParameter>
    !!]
    self=metaTreeProcessingTimeFile(fileName)
    !![
    <inputParametersValidate source="parameters"/>
    !!]
    return
  end function fileConstructorParameters

  function fileConstructorInternal(fileName) result(self)
    !!{RST
    Internal constructor for the :galacticus-class:`metaTreeProcessingTimeFile` merger tree processing time estimator class.
    !!}
    use :: FoX_DOM           , only : node
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5File             , hdf5Group
    use :: IO_XML            , only : XML_Array_Read_Static, XML_Get_First_Element_By_Tag_Name, XML_Parse
    use :: ISO_Varying_String, only : varying_string       , char                             , operator(//)
    implicit none
    type     (metaTreeProcessingTimeFile)                            :: self
    type     (varying_string            ), intent(in   )             :: fileName
    type     (node                      ), pointer                   :: doc              , fit
    type     (hdf5File                  )                            :: file
    type     (hdf5Group                 )                            :: metaDataGroup    , timingDataGroup
    double precision                     , allocatable, dimension(:) :: coefficients
    integer                                                          :: ioStatus         , lengthName
    character(len=1024                  )                            :: fileNameCharacter
    !![
    <constructorAssign variables="fileName"/>
    !!]

    fileNameCharacter=char    (fileName         )
    lengthName       =len_trim(fileNameCharacter)
    if     (                                                                                      &
         &   (lengthName >= 5 .and. fileNameCharacter(max(1,lengthName-4):lengthName) == ".hdf5") &
         &  .or.                                                                                  &
         &   (lengthName >= 3 .and. fileNameCharacter(max(1,lengthName-2):lengthName) == ".h5"  ) &
         & ) then
       ! The file is a Galacticus HDF5 output file - read the fit coefficients written by the tree processing timer operator.
       !$ call hdf5Access%set()
       file=hdf5File(fileName,readOnly=.true.)
       if (.not.file%hasGroup('metaData')) call Error_Report('tree timing file "'//char(fileName)//'" contains no "metaData" group'//{introspection:location})
       metaDataGroup=file%openGroup('metaData')
       if (.not.metaDataGroup%hasGroup('treeTiming')) call Error_Report('tree timing file "'//char(fileName)//'" contains no "metaData/treeTiming" group'//{introspection:location})
       timingDataGroup=metaDataGroup%openGroup('treeTiming')
       if (.not.timingDataGroup%hasDataset('fitCoefficientMass')) call Error_Report('tree timing file "'//char(fileName)//'" contains no "metaData/treeTiming/fitCoefficientMass" dataset - was the "mergerTreeOperatorTreeProcessingTimer" operator active in the run that produced it?'//{introspection:location})
       call timingDataGroup%readDataset('fitCoefficientMass',coefficients)
       if (size(coefficients) /= 3) call Error_Report('expected 3 fit coefficients in tree timing file "'//char(fileName)//'"'//{introspection:location})
       self%fitCoefficient=coefficients
       ! Also read the node-count-based fit if present (used for the read path, where node counts are known but masses are not).
       if (timingDataGroup%hasDataset('fitCoefficientCountNodes')) then
          call timingDataGroup%readDataset('fitCoefficientCountNodes',coefficients)
          if (size(coefficients) == 3) then
             self%fitCoefficientCountNodes=coefficients
             self%haveCountNodesFit       =.true.
          end if
       end if
       !$ call hdf5Access%unset()
    else
       ! Parse the fit file as an XML document.
       !$omp critical (FoX_DOM_Access)
       doc => XML_Parse(fileName,iostat=ioStatus)
       if (ioStatus /= 0) call Error_Report('Unable to find or parse tree timing file'//{introspection:location})
       fit => XML_Get_First_Element_By_Tag_Name(doc,"fit")
       call XML_Array_Read_Static(fit,"coefficient",self%fitCoefficient)
       !$omp end critical (FoX_DOM_Access)
    end if
    return
  end function fileConstructorInternal

  double precision function fileTime(self,massTree)
    !!{RST
    Return the units of the file property in the SI system.
    !!}
    implicit none
    class           (metaTreeProcessingTimeFile), intent(inout) :: self
    double precision                            , intent(in   ) :: massTree
    integer                                                     :: i
    double precision                                            :: massTreeLogarithmic

    massTreeLogarithmic=log10(massTree)
    fileTime=0.0d0
    do i=0,2
       fileTime=fileTime+self%fitCoefficient(i)*massTreeLogarithmic**i
    end do
    fileTime=10.0d0**fileTime
    return
  end function fileTime

  double precision function fileTimeByCountNodes(self,countNodes)
    !!{RST
    Return an estimate of the time needed to process a tree with the given number of nodes, using the node-count-based fit read from
    file. Returns a negative value if no node-count-based fit is available (e.g. the file was supplied in the legacy XML format,
    which contains only the mass-based fit).
    !!}
    use, intrinsic :: ISO_C_Binding, only : c_size_t
    implicit none
    class           (metaTreeProcessingTimeFile), intent(inout) :: self
    integer         (c_size_t                  ), intent(in   ) :: countNodes
    integer                                                     :: i
    double precision                                            :: countNodesLogarithmic

    if (.not.self%haveCountNodesFit .or. countNodes <= 0_c_size_t) then
       fileTimeByCountNodes=-1.0d0
       return
    end if
    countNodesLogarithmic=log10(dble(countNodes))
    fileTimeByCountNodes=0.0d0
    do i=0,2
       fileTimeByCountNodes=fileTimeByCountNodes+self%fitCoefficientCountNodes(i)*countNodesLogarithmic**i
    end do
    fileTimeByCountNodes=10.0d0**fileTimeByCountNodes
    return
  end function fileTimeByCountNodes

