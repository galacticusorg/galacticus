!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016
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

!% Implements a file-based transfer function class.

  use Tables
  
  !# <transferFunction name="transferFunctionFile" defaultThreadPrivate="yes">
  !#  <description>
  !# Provides a transfer function from a tabulation given in an XML file.The XML file format for transfer functions looks like:
  !# \begin{verbatim}
  !#  <data>
  !#   <column>k [Mpc^{-1}] - wavenumber</column>
  !#   <column>T(k) - transfer function</column>
  !#   <datum>1.111614e-05 0.218866E+08</datum>
  !#   <datum>1.228521e-05 0.218866E+08</datum>
  !#   <datum>1.357727e-05 0.218866E+08</datum>
  !#   <datum>1.50052e-05 0.218866E+08</datum>
  !#   <datum>1.658335e-05 0.218866E+08</datum>
  !#   <datum>1.83274e-05 0.218865E+08</datum>
  !#   .
  !#   .
  !#   .
  !#   <description>Cold dark matter power spectrum created by CAMB.</description>
  !#   <fileFormat>1</fileFormat>
  !#   <parameter>
  !#     <name>Omega_b</name>
  !#     <value>0.0450</value>
  !#   </parameter>
  !#   <parameter>
  !#     <name>Omega_Matter</name>
  !#     <value>0.250</value>
  !#   </parameter>
  !#   <parameter>
  !#     <name>Omega_DE</name>
  !#     <value>0.750</value>
  !#   </parameter>
  !#   <parameter>
  !#     <name>H_0</name>
  !#     <value>70.0</value>
  !#   </parameter>
  !#   <parameter>
  !#     <name>T_CMB</name>
  !#     <value>2.780</value>
  !#   </parameter>
  !#   <parameter>
  !#     <name>Y_He</name>
  !#     <value>0.24</value>
  !#   </parameter>
  !#   <extrapolation>
  !#     <wavenumber>
  !#       <limit>low</limit>
  !#       <method>extrapolate</method>
  !#     </wavenumber>
  !#     <wavenumber>
  !#       <limit>high</limit>
  !#       <method>extrapolate</method>
  !#     </wavenumber>
  !#   </extrapolation>
  !# </data>
  !# \end{verbatim}
  !# The {\normalfont \ttfamily datum} elements give wavenumber (in Mpc$^{-1}$) and transfer function pairs. The {\normalfont
  !# \ttfamily extrapolation} element defines how the tabulated function should be extrapolated to lower and higher
  !# wavenumbers. The two options for the {\normalfont \ttfamily method} are ``fixed'', in which case the transfer function is
  !# extrapolated assuming that it remains constant, and ``extrapolate'' in which case the extrapolation is performed (typically
  !# this extrapolation is based on extending the cubic spline interpolation used to interpolate the transfer function in a
  !# log-log space). The {\normalfont \ttfamily column}, {\normalfont \ttfamily description} and {\normalfont \ttfamily parameter}
  !# elements are optional, but are encouraged to make the file easier to understand. Finally, the {\normalfont \ttfamily
  !# fileFormat} element should currently always contain the value $1$---this may change in future if the format of this file is
  !# modified.
  !# </description>
  !# </transferFunction>
  type, extends(transferFunctionClass) :: transferFunctionFile
     !% A transfer function class which interpolates a transfer function given in a file.
     private
     type(table1DGeneric) :: transfer
   contains
     !@ <objectMethods>
     !@   <object>transferFunctionFile</object>
     !@   <objectMethod>
     !@     <method>readFile</method>
     !@     <type>void</type>
     !@     <arguments>\textcolor{red}{\textless char(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Read the named transfer function file.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                          fileDestructor
     procedure :: readFile              => fileReadFile
     procedure :: value                 => fileValue
     procedure :: logarithmicDerivative => fileLogarithmicDerivative
     procedure :: halfModeMass          => fileHalfModeMass
  end type transferFunctionFile

  interface transferFunctionFile
     !% Constructors for the file transfer function class.
     module procedure fileConstructorParameters
     module procedure fileConstructorInternal
  end interface transferFunctionFile

  ! Current file format version for transfer function files.
  integer, parameter :: fileFormatVersionCurrent=1

contains

  function fileConstructorParameters(parameters)
    !% Constructor for the file transfer function class which takes a parameter set as input.
    implicit none
    type(transferFunctionFile)                :: fileConstructorParameters
    type(inputParameters     ), intent(inout) :: parameters
    type(varying_string      )                :: fileName
    !# <inputParameterList label="allowedParameterNames" />

    call parameters%checkParameters(allowedParameterNames)    
    !# <inputParameter>
    !#   <name>fileName</name>
    !#   <source>parameters</source>
    !#   <description>The name of the file from which to read a tabulated transfer function.</description>
    !#   <type>string</type>
    !#   <cardinality>0..1</cardinality>
    !# </inputParameter>
    fileConstructorParameters=fileConstructorInternal(char(fileName))
    return
  end function fileConstructorParameters
  
  function fileConstructorInternal(fileName)
    !% Internal constructor for the file transfer function class.
    use Input_Parameters2
    use Cosmology_Parameters
    use FoX_DOM
    use IO_XML
    use Numerical_Comparison
    use Array_Utilities
    use Galacticus_Error
    use Galacticus_Display
    implicit none
    type     (transferFunctionFile)                :: fileConstructorInternal
    character(len=*               ), intent(in   ) :: fileName

    call fileConstructorInternal%readFile(fileName)
    return
  end function fileConstructorInternal
  
  subroutine fileReadFile(self,fileName)
    !% Internal constructor for the file transfer function class.
    use Input_Parameters2
    use Cosmology_Parameters
    use FoX_DOM
    use IO_XML
    use Numerical_Comparison
    use Array_Utilities
    use Galacticus_Error
    use Galacticus_Display
    use Table_Labels
    implicit none
    class           (transferFunctionFile    ), intent(inout)             :: self
    character       (len=*                   ), intent(in   )             :: fileName
    type            (Node                    ), pointer                   :: doc                             , extrapolation              , &
         &                                                                   extrapolationElement            , formatElement              , &
         &                                                                   thisParameter                   , parameters
    type            (NodeList                ), pointer                   :: wavenumberExtrapolationList
    double precision                          , allocatable, dimension(:) :: transfer                        , wavenumber                 , &
         &                                                                   transferLogarithmic             , wavenumberLogarithmic
    class           (cosmologyParametersClass), pointer                   :: cosmologyParameters_            , cosmologyParametersFile
    double precision                          , parameter                 :: toleranceUniformity      =1.0d-6
    type            (inputParameters         )                            :: transferFunctionCosmology
    integer                                                               :: extrapolationMethod             , versionNumber              , &
         &                                                                   iExtrapolation                  , ioError                    , &
         &                                                                   extrapolateWavenumberLow        , extrapolateWavenumberHigh
    character       (len=32                  )                            :: limitType

    ! Get the default cosmology.
    cosmologyParameters_ => cosmologyParameters()
    ! Open and parse the data file.
    !$omp critical (FoX_DOM_Access)
    doc => parseFile(fileName,iostat=ioError)
    if (ioError /= 0) call Galacticus_Error_Report('fileReadFile','Unable to find transfer function file')
    ! Check that the file has the correct format version number.
    formatElement => XML_Get_First_Element_By_Tag_Name(doc,"fileFormat")
    call extractDataContent(formatElement,versionNumber)
    if (versionNumber /= fileFormatVersionCurrent) call Galacticus_Error_Report('fileReadFile','file has the incorrect version number')
    ! Check that parameters match if any are present.
    parameters => XML_Get_First_Element_By_Tag_Name(doc,"parameters")
    !$omp end critical (FoX_DOM_Access)
    transferFunctionCosmology=inputParameters(parameters)
    cosmologyParametersFile => cosmologyParameters(transferFunctionCosmology)
    if (Values_Differ(cosmologyParametersFile%OmegaBaryon    (),cosmologyParameters_%OmegaBaryon    (),absTol=1.0d-3)) &
         & call Galacticus_Display_Message('OmegaBaryon from transfer function file does not match internal value'    )
    if (Values_Differ(cosmologyParametersFile%OmegaMatter    (),cosmologyParameters_%OmegaMatter    (),absTol=1.0d-3)) &
         & call Galacticus_Display_Message('OmegaMatter from transfer function file does not match internal value'    )
    if (Values_Differ(cosmologyParametersFile%OmegaDarkEnergy(),cosmologyParameters_%OmegaDarkEnergy(),absTol=1.0d-3)) &
         & call Galacticus_Display_Message('OmegaDarkEnergy from transfer function file does not match internal value')
    if (Values_Differ(cosmologyParametersFile%HubbleConstant (),cosmologyParameters_%HubbleConstant (),relTol=1.0d-3)) &
         & call Galacticus_Display_Message('HubbleConstant from transfer function file does not match internal value' )
    if (Values_Differ(cosmologyParametersFile%temperatureCMB (),cosmologyParameters_%temperatureCMB (),relTol=1.0d-3)) &
         & call Galacticus_Display_Message('temperatureCMB from transfer function file does not match internal value' )
    ! Get extrapolation methods.
    !$omp critical (FoX_DOM_Access)
    extrapolationElement        => XML_Get_First_Element_By_Tag_Name(doc                 ,"extrapolation")
    wavenumberExtrapolationList => getElementsByTagname             (extrapolationElement,"wavenumber"   )
    extrapolateWavenumberLow    =  extrapolationTypeExtrapolate
    extrapolateWavenumberHigh   =  extrapolationTypeExtrapolate
    do iExtrapolation=0,getLength(wavenumberExtrapolationList)-1
       extrapolation => item(wavenumberExtrapolationList,iExtrapolation)
       call XML_Extrapolation_Element_Decode(                                                   &
            &                                extrapolation                                    , &
            &                                limitType                                        , &
            &                                extrapolationMethod                              , &
            &                                allowedMethods     =[                              &
            &                                                     extrapolationTypeFix        , &
            &                                                     extrapolationTypeExtrapolate  &
            &                                                    ]                              &
            &                               )
       select case (trim(limitType))
       case ('low' )
          extrapolateWavenumberLow =extrapolationMethod
       case ('high')
          extrapolateWavenumberHigh=extrapolationMethod
       case default
          call Galacticus_Error_Report('fileReadFile','unrecognized extrapolation limit')
       end select
    end do
    ! Read the transfer function from file.
    call XML_Array_Read(doc,"datum",wavenumber,transfer)
    ! Destroy the document.
    call destroy(doc)
    !$omp end critical (FoX_DOM_Access)
    ! Construct the tabulated transfer function.
    call self%transfer%destroy()
    wavenumberLogarithmic=log(wavenumber)
    transferLogarithmic  =log(transfer  )   
    ! Create the table.
    call self%transfer%create  (                                                    &
         &                      wavenumberLogarithmic                             , &
         &                      extrapolationType=[                                 &
         &                                         extrapolateWavenumberLow       , &
         &                                         extrapolateWavenumberHigh        &
         &                                        ]                               , &
         &                      interpolationType= FGSL_Interp_cSpline              &
         &                     )
    call self%transfer%populate(                                                    &
         &                      transferLogarithmic                                 &
         &                     )
    return
  end subroutine fileReadFile

  subroutine fileDestructor(self)
    !% Destructor for the file transfer function class.
    implicit none
    type(transferFunctionFile), intent(inout) :: self

    call self%transfer%destroy()
    return
  end subroutine fileDestructor

  double precision function fileValue(self,wavenumber)
    !% Return the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionFile), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber

    fileValue=exp(self%transfer%interpolate(log(wavenumber)))
    return
  end function fileValue

  double precision function fileLogarithmicDerivative(self,wavenumber)
    !% Return the logarithmic derivative of the transfer function at the given wavenumber.
    implicit none
    class           (transferFunctionFile), intent(inout) :: self
    double precision                      , intent(in   ) :: wavenumber

    fileLogarithmicDerivative=+self%transfer%interpolateGradient(log(wavenumber))
    return
  end function fileLogarithmicDerivative
  
  double precision function fileHalfModeMass(self)
    !% Compute the mass corresponding to the wavenumber at which the transfer function is
    !% suppressed by a factor of two relative to a \gls{cdm} transfer function. Not supported in
    !% this implementation.
    use Galacticus_Error
    implicit none
    class(transferFunctionFile), intent(inout) :: self
    
    call Galacticus_Error_Report('fileHalfModeMass','not supported by this implementation')
    return
  end function fileHalfModeMass
