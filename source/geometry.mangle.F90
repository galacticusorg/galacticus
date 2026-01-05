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

!!{
Contains a module which implements functions utilizing \gls{mangle} survey geometry definitions.
!!}

module Geometry_Mangle
  !!{
  Implements functions utilizing \gls{mangle} survey geometry definitions.
  !!}
  use, intrinsic :: ISO_C_Binding, only : c_size_t
  private
  public :: geometryMangleSolidAngle, geometryMangleAngularPower, &
       &    geometryMangleBuild     , window

  type :: cap
     !!{
     A class to hold \gls{mangle} caps.
     !!}
     double precision, dimension(3) :: x
     double precision               :: c
   contains
     !![
     <methods>
       <method description="Return true if the given point lives inside the {\normalfont \scshape mangle} cap." method="pointIncluded" />
     </methods>
     !!]
    procedure :: pointIncluded => capPointIncluded
  end type cap

  type :: polygon
     !!{
     A class to hold \gls{mangle} polygons.
     !!}
     integer                                          :: capCount
     double precision                                 :: weight  , solidAngle
     type            (cap), dimension(:), allocatable :: caps
   contains
     !![
     <methods>
       <method description="Return true if the given point lives inside the {\normalfont \scshape mangle} polygon." method="pointIncluded" />
     </methods>
     !!]
     procedure :: pointIncluded => polygonPointIncluded
  end type polygon

  type :: window
     !!{
     A class to hold \gls{mangle} windows.
     !!}
     integer                                           :: polygonCount
     type   (polygon      ), dimension(:), allocatable :: polygons
     integer(kind=c_size_t), dimension(:), allocatable :: solidAngleIndex
   contains
     !![
     <methods>
       <method description="Read the specified {\normalfont \scshape mangle} polygon file."                        method="read"         />
       <method description="Return true if the given point lives inside the {\normalfont \scshape mangle} window." method="pointIncluded"/>
     </methods>
     !!]
     procedure :: read          => windowRead
     procedure :: pointIncluded => windowPointIncluded
  end type window

contains

  subroutine windowRead(self,fileName)
    !!{
    Read a \gls{mangle} window definition from file.
    !!}
    use :: Display           , only : displayCounter    , displayCounterClear, displayIndent, displayMessage, &
          &                           displayUnindent
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : assignment(=)     , operator(//)       , operator(==) , varying_string
    use :: Sorting           , only : sortIndex
    use :: String_Handling   , only : String_Split_Words, operator(//)
    implicit none
    class           (window        ), intent(inout)              :: self
    character       (len=*         ), intent(in   )              :: fileName
    double precision                , allocatable, dimension(: ) :: solidAngle
    type            (varying_string)             , dimension(11) :: words
    type            (varying_string)                             :: message
    character       (len=1024      )                             :: line
    integer                                                      :: fileUnit  , i, &
         &                                                          ioStatus  , j

    ! Open and parse the file.
    message='Reading mangle window from: '//trim(fileName)
    call displayIndent(message)
    open(newUnit=fileUnit,file=fileName,status='old',form='formatted')
    ! Retrieve count of number of polygons.
    self%polygonCount=-1
    do while (self%polygonCount < 0)
       read (fileUnit,'(a)') line
       if (line(1:1) /= "#") read (line,*) self%polygonCount
    end do
    message='found '
    message=message//self%polygonCount//' polygons'
    call displayMessage(message)
     allocate(self%polygons       (self%polygonCount))
    allocate(solidAngle          (self%polygonCount))
    allocate(self%solidAngleIndex(self%polygonCount))
    ! Read each polygon.
    do i=1,self%polygonCount
       call displayCounter(int(100.0d0*dble(i-1)/dble(self%polygonCount)),isNew=(i==1))
       ioStatus=0
       do while (ioStatus == 0)
          read (fileUnit,'(a)',iostat=ioStatus) line
          if (ioStatus /= 0)  call Error_Report('end of file reached'//{introspection:location})
          call String_Split_Words(words,line)
          if (words(1) == "polygon") exit
       end do
       line=words(4)
       read (line,*) self%polygons(i)%capCount
       line=words(6)
       read (line,*) self%polygons(i)%weight
       line=words(8)
       read (line,*) self%polygons(i)%solidAngle
       solidAngle(i)=self%polygons(i)%solidAngle
       allocate(self%polygons(i)%caps(self%polygons(i)%capCount))
       do j=1,self%polygons(i)%capCount
          read (fileUnit,*) self%polygons(i)%caps(j)%x, &
               &            self%polygons(i)%caps(j)%c
       end do
    end do
    call displayCounterClear()
    close(fileUnit)
    ! Get a sorted index of polygon solid angles.
    self%solidAngleIndex=sortIndex(solidAngle)
    deallocate(solidAngle)
    call displayUnindent('done')
    return
  end subroutine windowRead

  logical function windowPointIncluded(self,point)
    !!{
    Return true if the given Cartesian point lies inside a \gls{mangle} window, i.e. if it lies within any polygon of the window.
    !!}
    implicit none
    class           (window), intent(inout)               :: self
    double precision        , intent(in   ), dimension(3) :: point
    integer                                               :: i

    ! Check each polygon in turn, in order of decreasing solid angle (since, at least for randomly distributed points, we're most
    ! likely to find the point in a polygon with larger solid angle). Exit immediately if the point is found to like within a
    ! polygon.
    windowPointIncluded=.false.
    do i=self%polygonCount,1,-1
       windowPointIncluded=self%polygons(self%solidAngleIndex(i))%pointIncluded(point)
       if (windowPointIncluded) exit
    end do
    return
  end function windowPointIncluded

  logical function polygonPointIncluded(self,point)
    !!{
    Return true if a given Cartesian point lies within a \gls{mangle} polygon, i.e. lies within \emph{all} of the polygons caps.
    !!}
    implicit none
    class           (polygon), intent(inout)               :: self
    double precision         , intent(in   ), dimension(3) :: point
    integer                                                :: i

    polygonPointIncluded=.true.
    do i=1,self%capCount
       polygonPointIncluded=self%caps(i)%pointIncluded(point)
       if (.not.polygonPointIncluded) exit
    end do
    return
  end function polygonPointIncluded

  logical function capPointIncluded(self,point)
    !!{
    Return true if a given Cartesian point lies within a \gls{mangle} cap.
    !!}
    use :: Vectors, only : Vector_Magnitude
    implicit none
    class           (cap), intent(inout)               :: self
    double precision     , intent(in   ), dimension(3) :: point
    double precision                                   :: cosinePolarAngle

    ! Find the cosine of the angle between the north pole of the cap and the given point.
    cosinePolarAngle=sum(point*self%x)/Vector_Magnitude(point)
    ! Determine if the point lies within the cap.
    if (self%c > 0.0d0) then
       capPointIncluded=1.0d0-cosinePolarAngle < abs(self%c)
    else
       capPointIncluded=1.0d0-cosinePolarAngle > abs(self%c)
    end if
    return
  end function capPointIncluded

  subroutine geometryMangleBuild(manglePath,mangleVersion,static)
    !!{
    Download and build the \textsc{mangle} code.
    !!}
    use :: Dependencies      , only : dependencyVersion
    use :: Display           , only : displayMessage   , verbosityLevelWorking
    use :: File_Utilities    , only : Directory_Make   , File_Exists          , File_Lock      , File_Unlock  , &
          &                           lockDescriptor
    use :: Error             , only : Error_Report
    use :: Input_Paths       , only : inputPath        , pathTypeDataDynamic
    use :: ISO_Varying_String, only : operator(//)     , varying_string       , char           , assignment(=)
    use :: String_Handling   , only : stringSubstitute
    use :: System_Download   , only : download
    use :: System_Command    , only : System_Command_Do
    use :: System_Compilers  , only : compiler         , compilerOptions      , languageFortran, languageC
    implicit none
    type   (varying_string), intent(  out)           :: manglePath, mangleVersion
    logical                , intent(in   ), optional :: static
    integer                                          :: status
    type   (lockDescriptor)                          :: fileLock
    type   (varying_string)                          :: command   , staticOptions
    !![
    <optionalArgument name="static" defaultsTo=".false." />
    !!]

    ! Get the version to use.
    mangleVersion=dependencyVersion("mangle")
    manglePath   =inputPath(pathTypeDataDynamic)//"mangle-"//mangleVersion//"/"
    ! Build the mangle code.
    if (.not.File_Exists(manglePath//"bin/harmonize")) then
       call Directory_Make(     manglePath                                         )
       call File_Lock     (char(manglePath//"mangle"),fileLock,lockIsShared=.false.)
       ! Unpack the code.
       if (.not.File_Exists(manglePath//"src/Makefile.in")) then
          ! Download mangle if necessary.
          if (.not.File_Exists(manglePath//".tar.gz")) then
             call displayMessage("downloading mangle code....",verbosityLevelWorking)
             call download("https://github.com/galacticusorg/mangle/archive/refs/tags/v"//char(mangleVersion)//".tar.gz",char(manglePath//".tar.gz"),status=status,retries=5,retryWait=60)
             if (status /= 0 .or. .not.File_Exists(manglePath//".tar.gz")) call Error_Report("unable to download mangle"//{introspection:location})
          end if
          call displayMessage("unpacking mangle code....",verbosityLevelWorking)
          call System_Command_Do("tar -x -v -z -C "//inputPath(pathTypeDataDynamic)//" -f "//manglePath//".tar.gz",status)
          if (status /= 0 .or. .not.File_Exists(manglePath//"src/Makefile.in")) call Error_Report('failed to unpack mangle code'//{introspection:location})
       end if
       staticOptions=""
       if (static_) staticOptions=" -Wl,--whole-archive -lpthread -ldl -Wl,--no-whole-archive"
       command=         'cd '//manglePath//'src; ./configure; '
       command=command//'sed -E -i~ s/"^F77[[:space:]]*=[[:space:]]*[a-zA-Z0-9]+"/"F77 = '//                 compiler       (languageFortran)                         //'"/ Makefile; '
       command=command//'sed -E -i~ s/"^CC[[:space:]]*=[[:space:]]*[a-zA-Z0-9]+"/"CC = '  //                 compiler       (languageC      )                         //'"/ Makefile; '
       command=command//'sed -E -i~ s/"^FFLAGS[[:space:]]*:=(.*)"/"FFLAGS:=\1 '           //stringSubstitute(compilerOptions(languageFortran),"/","\/")//staticOptions//'"/ Makefile; '
       command=command//'sed -E -i~ s/"^CFLAGS[[:space:]]*:=(.*)"/"CFLAGS:=\1 '           //stringSubstitute(compilerOptions(languageC      ),"/","\/")//staticOptions//'"/ Makefile; '
       command=command//'make cleanest; make'
       if (static_) command=command//' static'
       call System_Command_Do(char(command),status);
       if (status /= 0 .or. .not.File_Exists(manglePath//"bin/harmonize")) call Error_Report("failed to build mangle executables"//{introspection:location})
       call File_Unlock(fileLock)
    end if
    return
  end subroutine geometryMangleBuild

 function geometryMangleSolidAngle(fileNames,solidAngleFileName)
    !!{
    Compute the solid angle of a \textsc{mangle} geometry.
    !!}
    use :: File_Utilities          , only : File_Exists       , File_Name_Temporary
    use :: Error                   , only : Error_Report
    use :: HDF5_Access             , only : hdf5Access
    use :: IO_HDF5                 , only : hdf5Object
    use :: ISO_Varying_String      , only : char              , extract            , len               , operator(//), &
          &                                 operator(==)      , varying_string
    use :: Numerical_Constants_Math, only : Pi
    use :: String_Handling         , only : String_Count_Words, String_Join        , String_Split_Words, char
    use :: System_Command          , only : System_Command_Do
    implicit none
    type            (varying_string), intent(in   ), dimension(             : ) :: fileNames
    character       (len=*         ), intent(in   ), optional                   :: solidAngleFileName
    double precision                               , dimension(size(fileNames)) :: geometryMangleSolidAngle
    type            (varying_string), allocatable  , dimension(             : ) :: subFiles
    integer                                                                     :: i                       , j            , &
         &                                                                         status                  , wlmFile
    type            (varying_string)                                            :: fileName                , fileNameTmp  , &
         &                                                                         manglePath              , mangleVersion
    double precision                                                            :: multiplier              , subSolidAngle, &
         &                                                                         w00
    type            (hdf5Object    )                                            :: solidAngleFile

    ! Check for pre-existing calculation.
    if (present(solidAngleFileName).and.File_Exists(solidAngleFileName)) then
       !$ call hdf5Access%set  ()
       call solidAngleFile%openFile         (solidAngleFileName,overWrite=.false.                 )
       call solidAngleFile%readDatasetStatic('solidAngle'      ,          geometryMangleSolidAngle)
       call solidAngleFile%close            (                                                     )
       !$ call hdf5Access%unset()
       return
    end if
    ! Ensure mangle is available.
    call geometryMangleBuild(manglePath,mangleVersion)
    ! Determine the solid angle of each file.
    do i=1,size(fileNames)
       allocate(subFiles(String_Count_Words(char(fileNames(i)),":")))
       call String_Split_Words(subFiles,char(fileNames(i)),":")
       geometryMangleSolidAngle(i)=0.0d0
       do j=1,size(subFiles)
          fileName  =subFiles(j)
          multiplier=+1.0d0
          if (extract(fileName,1,1) == "+") then
             fileName  =extract(fileName,2,len(fileName))
             multiplier=+1.0d0
          else if (extract(fileName,1,1) == "-") then
             fileName  =extract(fileName,2,len(fileName))
             multiplier=-1.0d0
          end if
          fileNameTmp=File_Name_Temporary('geometryMangleSolidAngle')
          call System_Command_Do(manglePath//"bin/harmonize "//fileName//" "//fileNameTmp,status)
          if (status /= 0) call Error_Report('failed to run mangle harmonize'//{introspection:location})
          open(newUnit=wlmFile,file=char(fileNameTmp),status="old",form="formatted")
          read (wlmFile,*)
          read (wlmFile,*) w00
          close(wlmFile)
          subSolidAngle      =2.0d0*sqrt(Pi)*w00
          geometryMangleSolidAngle(i)=geometryMangleSolidAngle(i)+multiplier*subSolidAngle
       end do
       deallocate(subFiles)
    end do
    ! Store the solid angle to file.
    if (present(solidAngleFileName)) then
       !$ call hdf5Access%set  ()
       call solidAngleFile%openFile      (            solidAngleFileName           ,overWrite=.true.      )
       call solidAngleFile%writeAttribute(String_Join(fileNames               ,":"),          'files'     )
       call solidAngleFile%writeDataset  (            geometryMangleSolidAngle     ,          'solidAngle')
       call solidAngleFile%flush         (                                                                )
       call solidAngleFile%close         (                                                                )
       !$ call hdf5Access%unset()
    end if
    return
  end function geometryMangleSolidAngle

  function geometryMangleAngularPower(fileNames,degreeMaximum,angularPowerFileName)
    !!{
    Compute the angular power spectra of a \textsc{mangle} geometry.
    !!}
    use :: File_Utilities    , only : File_Exists       , File_Name_Temporary
    use :: Error             , only : Error_Report
    use :: HDF5_Access       , only : hdf5Access
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char              , extract            , len               , operator(//), &
          &                           operator(==)      , var_str            , varying_string
    use :: String_Handling   , only : String_Count_Words, String_Join        , String_Split_Words, operator(//)
    use :: System_Command    , only : System_Command_Do
    implicit none
    type            (varying_string), intent(in   ), dimension(             :                                                               ) :: fileNames
    integer                         , intent(in   )                                                                                           :: degreeMaximum
    character       (len=*         ), intent(in   ), optional                                                                                 :: angularPowerFileName
    double precision                               , dimension(size(fileNames)*(size(fileNames)+1)/2  , degreeMaximum+1                     ) :: geometryMangleAngularPower
    double precision                               , dimension(                                      2                                      ) :: coefficient
    double precision                               , dimension(size(fileNames)                      ,2,(degreeMaximum+1)*(degreeMaximum+2)/2) :: coefficients
    type            (varying_string), allocatable  , dimension(             :                                                               ) :: subFiles
    type            (varying_string)                                                                                                          :: fileName                  , fileNameTmp  , &
         &                                                                                                                                       manglePath                , mangleVersion
    integer                                                                                                                                   :: degree                    , order        , &
         &                                                                                                                                       status                    , wlmFile      , &
         &                                                                                                                                       i                         , j            , &
         &                                                                                                                                       k                         , l            , &
         &                                                                                                                                       p                         , q
    double precision                                                                                                                          :: multiplier                , weight
    type            (hdf5Object    )                                                                                                          :: angularPowerFile

    ! Read the angular power from file if possible.
    if (present(angularPowerFileName).and.File_Exists(angularPowerFileName)) then
       !$ call hdf5Access%set  ()
       call angularPowerFile      %openFile         (angularPowerFileName                                           )
       l=0
       do p=1,size(fileNames)
          do q=p,size(fileNames)
             l=l+1
             call angularPowerFile%readDatasetStatic(char(var_str('Cl_')//(p-1)//'_'//(q-1)),geometryMangleAngularPower(l,:))
          end do
       end do
       call angularPowerFile      %flush            (                                                               )
       call angularPowerFile      %close            (                                                               )
       !$ call hdf5Access%unset()
       return
    end if
    ! Ensure mangle is available.
    call geometryMangleBuild(manglePath,mangleVersion)
    ! Determine the spherical harmonic coefficients of each file.
    coefficients=0.0d0
    do i=1,size(fileNames)
       allocate(subFiles(String_Count_Words(char(fileNames(i)),":")))
       call String_Split_Words(subFiles,char(fileNames(i)),":")
       do j=1,size(subFiles)
          fileName  =subFiles(j)
          multiplier=+1.0d0
          if (extract(fileName,1,1) == "+") then
             fileName  =extract(fileName,2,len(fileName))
             multiplier=+1.0d0
          else if (extract(fileName,1,1) == "-") then
             fileName  =extract(fileName,2,len(fileName))
             multiplier=-1.0d0
          end if
          fileNameTmp=File_Name_Temporary('geometryMangleAngularPower')
          call System_Command_Do(manglePath//"bin/harmonize -l "//degreeMaximum//" "//fileName//" "//fileNameTmp,status)
          if (status /= 0) call Error_Report('failed to run mangle harmonize'//{introspection:location})
          open(newUnit=wlmFile,file=char(fileNameTmp),status="old",form="formatted")
          read (wlmFile,*)
          k=0
          do degree=0,degreeMaximum
             do order=0,degree
                k=k+1
                read (wlmFile,*) coefficient
                coefficients(i,:,k)=coefficients(i,:,k)+multiplier*coefficient
             end do
          end do
          close(wlmFile)
       end do
       deallocate(subFiles)
    end do
    ! Compute power and cross-power spectra.
    geometryMangleAngularPower=0.0d0
    k=0
    do degree=0,degreeMaximum
       do order=0,degree
          k=k+1
          ! Compute weight, which accounts for the fact that we sum only over positive order (m), as the negative orders are
          ! symmetric (since our window functions are always real).
          if (order == 0) then
             weight=1.0d0
          else
             weight=2.0d0
          end if
          l=0
          do p=1,size(fileNames)
             do q=p,size(fileNames)
                l=l+1
                geometryMangleAngularPower(l,degree+1)=geometryMangleAngularPower(l,degree+1)+weight*sum(coefficients(p,:,k)*coefficients(q,:,k))
             end do
          end do
       end do
       geometryMangleAngularPower(:,degree+1)=geometryMangleAngularPower(:,degree+1)/(2.0d0*float(degree)+1.0d0)
    end do
    ! Store the angular power to file.
    if (present(angularPowerFileName)) then
       !$ call hdf5Access%set  ()
       call angularPowerFile      %openFile      (            angularPowerFileName                ,overWrite=.true.                                 )
       call angularPowerFile      %writeAttribute(String_Join(fileNames                      ,":"),          'files'                                )
       l=0
       do p=1,size(fileNames)
          do q=p,size(fileNames)
             l=l+1
             call angularPowerFile%writeDataset  (            geometryMangleAngularPower(l,:)     ,          char(var_str('Cl_')//(p-1)//'_'//(q-1)))
          end do
       end do
       call angularPowerFile      %flush         (                                                                                                  )
       call angularPowerFile      %close         (                                                                                                  )
       !$ call hdf5Access%unset()
    end if
    return
  end function geometryMangleAngularPower

end module Geometry_Mangle
