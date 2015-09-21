!! Copyright 2009, 2010, 2011, 2012, 2013 Andrew Benson <abenson@obs.carnegiescience.edu>
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


  type, public, extends(stateHistory) :: stateCorrelation
     !% Defines a ``correlation'' state class.
     private
     integer                                       :: storedStateCount               , stepComputePrevious, &
          &                                           convergedCorrelationLengthCount
     double precision, allocatable, dimension(:,:) :: states
     integer         , allocatable, dimension(:  ) :: correlationLengths
   contains
     !@ <objectMethods>
     !@   <object>stateCorrelation</object>
     !@   <objectMethod>
     !@     <method>correlationLength</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the current correlation length in the chains.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>postConvergenceCorrelationCount</method>
     !@     <type>\intzero</type>
     !@     <arguments></arguments>
     !@     <description>Return the number of post-convergence correlation lengths that have accrued.</description>
     !@   </objectMethod>
     !@ </objectMethods>
     final     ::                                    stateCorrelationDestroy
     procedure :: update                          => stateCorrelationUpdate
     procedure :: reset                           => stateCorrelationReset
     procedure :: postConvergenceCorrelationCount => stateCorrelationPostConvergenceCorrelationCount
     procedure :: correlationLength               => stateCorrelationCorrelationLength
  end type stateCorrelation
  
  interface stateCorrelation
     !% Constructors for the ``correlation'' state class.
     module procedure stateCorrelationConstructor
  end interface stateCorrelation
