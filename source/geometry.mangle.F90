!! Copyright 2009, 2010, 2011, 2012 Andrew Benson <abenson@obs.carnegiescience.edu>
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

!% Contains a module which implements functions utilizing \gls{mangle} survey geometry definitions.

module Geometry_Mangle
  !% Implements functions utilizing \gls{mangle} survey geometry definitions.
  
  type :: cap
     !% A class to hold \gls{mangle} caps.
     double precision, dimension(3) :: x
     double precision               :: c
   contains
     !@ <objectMethods>
     !@   <object>cap</object>
     !@   <objectMethod>
     !@     <method>pointIncluded</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless double(3)\textgreater} point\argin</arguments>
     !@     <description>Return true if the given point lives inside the {\sc mangle} cap.</description>
     !@   </objectMethod>
     !@ </objectMethods>
    procedure :: pointIncluded => capPointIncluded
  end type cap
  
  type :: polygon
     !% A class to hold \gls{mangle} polygons.
     integer                                          :: capCount
     double precision                                 :: weight  , solidAngle
     type            (cap), dimension(:), allocatable :: caps
   contains
     !@ <objectMethods>
     !@   <object>polygon</object>
     !@   <objectMethod>
     !@     <method>pointIncluded</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless double(3)\textgreater} point\argin</arguments>
     !@     <description>Return true if the given point lives inside the {\sc mangle} polygon.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: pointIncluded => polygonPointIncluded
  end type polygon

  type :: window
     !% A class to hold \gls{mangle} windows.
     integer                                     :: polygonCount
     type   (polygon), dimension(:), allocatable :: polygons
   contains
     !@ <objectMethods>
     !@   <object>window</object>
     !@   <objectMethod>
     !@     <method>read</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless character(len=*)\textgreater} fileName\argin</arguments>
     !@     <description>Read the specified {\sc mangle} polygon file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>pointIncluded</method>
     !@     <type>\logicalzero</type>
     !@     <arguments>\textcolor{red}{\textless double(3)\textgreater} point\argin</arguments>
     !@     <description>Return true if the given point lives inside the {\sc mangle} window.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: read          => windowRead
     procedure :: pointIncluded => windowPointIncluded
  end type window
  
contains
  
  subroutine windowRead(self,fileName)
    !% Read a \gls{mangle} window definition from file.
    use Galacticus_Display
    use Galacticus_Error
    use ISO_Varying_String
    use String_Handling
    implicit none
    class    (window        ), intent(inout) :: self
    character(len=*         ), intent(in   ) :: fileName
    type     (varying_string), dimension(11) :: words
    type     (varying_string)                :: message
    character(len=1024      )                :: line
    integer                                  :: fileUnit, i, j

    ! Open and parse the file.
    message='Reading mangle window from: '//trim(fileName)
    call Galacticus_Display_Indent('message')
    open(newUnit=fileUnit,file=fileName,status='old',form='formatted')
    ! Retrieve count of number of polygons.
    read (fileUnit,*) self%polygonCount
    message='found '
    message=message//self%polygonCount//' polygons'
    call Galacticus_Display_Message(message)
    do i=1,3
       read (fileUnit,*) 
    end do
    allocate(self%polygons(self%polygonCount))
    ! Read each polygon.
    do i=1,self%polygonCount
       read (fileUnit,'(a)') line
       call String_Split_Words(words,line)
       if (words(1) /= "polygon") call Galacticus_Error_Report('windowRead','expected polygon')
       line=words(4)       
       read (line,*) self%polygons(i)%capCount
       line=words(6)       
       read (line,*) self%polygons(i)%weight
       line=words(8)       
       read (line,*) self%polygons(i)%solidAngle
       allocate(self%polygons(i)%caps(self%polygons(i)%capCount))
       do j=1,self%polygons(i)%capCount
          read (fileUnit,*) self%polygons(i)%caps(j)%x, &
               &            self%polygons(i)%caps(j)%c
       end do
    end do
    close(fileUnit)
    call Galacticus_Display_Unindent('done')
    return
  end subroutine windowRead
  
  logical function windowPointIncluded(self,point)
    !% Return true if the given Cartesian point lies inside a \gls{mangle} window, i.e. if it lies within any polygon of the window.
    implicit none
    class           (window), intent(inout)               :: self
    double precision        , intent(in   ), dimension(3) :: point
    integer                                               :: i
    
    do i=1,self%polygonCount
       windowPointIncluded=self%polygons(i)%pointIncluded(point)
       if (windowPointIncluded) exit
    end do
    return
  end function windowPointIncluded
  
  logical function polygonPointIncluded(self,point)
    !% Return true if a given Cartesian point lies within a \gls{mangle} polygon, i.e. lies within \emph{all} of the polygons caps.
    use Vectors
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
    !% Return true if a given Cartesian point lies within a \gls{mangle} cap.
    use Vectors
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

end module Geometry_Mangle
