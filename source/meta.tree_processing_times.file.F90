!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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
Implements a merger tree processing time estimator using a polynomial relation read from file.
!!}

  !![
  <metaTreeProcessingTime name="metaTreeProcessingTimeFile">
   <description>
    A merger tree processing time class which estimates processing times using a polynomial relation read from
    file. Specifically, the time taken to process a tree is estimate to be
    \begin{equation}
     \log_{10} [ \tau_\mathrm{tree}(M)] = \sum_{i=0}^2 C_i (\log_{10} M)^i,
    \end{equation}
    where $M$ is the root mass of the tree and the coefficients $C_i$ are read from a file, the name of which is specified via
    the {\normalfont \ttfamily [fileName]} parameter. This file should be an XML document with the structure:
    \begin{verbatim}
    &lt;timing>
     &lt;fit>
       &lt;coefficient>-0.73&lt;/coefficient>
       &lt;coefficient>-0.20&lt;/coefficient>
       &lt;coefficient>0.03&lt;/coefficient>
     &lt;/fit>
    &lt;/timing>
    \end{verbatim}
    where the array of coefficients give the values $C_0$, $C_1$ and $C_2$.
   </description>
  </metaTreeProcessingTime>
  !!]
  type, extends(metaTreeProcessingTimeClass) :: metaTreeProcessingTimeFile
     !!{
     A merger tree processing time estimator using a polynomial relation read from file.
     !!}
     private
     double precision                , dimension(0:2) :: fitCoefficient
     type            (varying_string)                 :: fileName
   contains
     procedure :: time => fileTime
  end type metaTreeProcessingTimeFile

  interface metaTreeProcessingTimeFile
     !!{
     Constructors for the ``file'' merger tree processing time estimator.
     !!}
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface metaTreeProcessingTimeFile

contains

  function fileConstructorParameters(parameters) result(self)
    !!{
    Constructor for the ``file'' merger tree processing time estimator class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type(metaTreeProcessingTimeFile)                :: self
    type(inputParameters           ), intent(inout) :: parameters
    type(varying_string            )                :: fileName

    !![
    <inputParameter>
      <name>fileName</name>
      <description>The name of the file which contains fit coefficients for the time per tree fitting function.</description>
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
    !!{
    Internal constructor for the ``file'' merger tree processing time estimator class.
    !!}
    use :: FoX_DOM           , only : node                 , parseFile
    use :: Error             , only : Error_Report
    use :: IO_XML            , only : XML_Array_Read_Static, XML_Get_First_Element_By_Tag_Name
    use :: ISO_Varying_String, only : varying_string       , char
    implicit none
    type   (metaTreeProcessingTimeFile)                :: self
    type   (varying_string            ), intent(in   ) :: fileName
    type   (node                      ), pointer       :: doc     , fit
    integer                                            :: ioStatus
    !![
    <constructorAssign variables="fileName"/>
    !!]
    
    ! Parse the fit file.
    !$omp critical (FoX_DOM_Access)
    doc => parseFile(char(fileName),iostat=ioStatus)
    if (ioStatus /= 0) call Error_Report('Unable to find or parse tree timing file'//{introspection:location})
    fit => XML_Get_First_Element_By_Tag_Name(doc,"fit")
    call XML_Array_Read_Static(fit,"coefficient",self%fitCoefficient)
    !$omp end critical (FoX_DOM_Access)
    return
  end function fileConstructorInternal

  double precision function fileTime(self,massTree)
    !!{
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

