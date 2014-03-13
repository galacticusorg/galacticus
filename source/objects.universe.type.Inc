!! Copyright 2009, 2010, 2011, 2012, 2013, 2014 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% Defines an object class for the simulation universe.

  type, public :: mergerTreeList
     !% A class used for building linked lists of merger trees.
     type(mergerTreeList), pointer :: next
     type(mergerTree    ), pointer :: tree
  end type mergerTreeList

  type, public :: universe
     !% The universe object class.
     type   (mergerTreeList), pointer         :: trees         => null()
     logical                                  :: allTreesBuilt =  .false.
     type   (universeEvent ), pointer, public :: event
   contains
     !@ <objectMethods>
     !@   <object>universe</object>
     !@   <objectMethod>
     !@    <method>createEvent</method>
     !@    <description>Create a {\tt treeEvent} object in this universe.</description>
     !@     <type>\textcolor{red}{\textless *type(universeEvent)\textgreater}</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@    <method>removeEvent</method>
     !@    <description>Remove a {\tt treeEvent} from this universe.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(universeEvent)\textgreater} event\argin</arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@    <method>popTree</method>
     !@    <description>Pop a {\tt mergerTree} from this universe.</description>
     !@     <type>\textcolor{red}{\textless *type(mergerTree)\textgreater}</type>
     !@     <arguments></arguments>
     !@   </objectMethod>
     !@   <objectMethod>
     !@    <method>pushTree</method>
     !@    <description>Pop a {\tt mergerTree} from this universe.</description>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless *type(mergerTree)\textgreater} thisTree\argin</arguments>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: createEvent => universeCreateEvent
     procedure :: removeEvent => universeRemoveEvent
     procedure :: popTree     => universePopTree
     procedure :: pushTree    => universePushTree
  end type universe

  type, public :: universeEvent
     !% Type for events attached to universes.
     private
     integer         (kind=kind_int8   )         , public :: ID   
     type            (universe         ), pointer, public :: universe 
     double precision                            , public :: time 
     type            (universeEvent    ), pointer, public :: next 
     procedure       (universeEventTask), pointer, public :: task 
  end type universeEvent
  
  ! Interface for universe event tasks.
  abstract interface
     logical function universeEventTask(thisEvent,thisUniverse)
       import universeEvent, universe
       class  (universeEvent), intent(in   ) :: thisEvent      
       type   (universe     ), intent(inout) :: thisUniverse       
     end function universeEventTask
  end interface
  
