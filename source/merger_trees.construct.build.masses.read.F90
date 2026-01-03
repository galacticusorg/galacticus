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
  Implementation of a merger tree masses class which reads masses from a file.
  !!}

  !![
  <mergerTreeBuildMasses name="mergerTreeBuildMassesRead" abstract="yes">
   <description>A merger tree masses class which samples masses from a distribution.</description>
   <runTimeFileDependencies paths="fileName"/>
  </mergerTreeBuildMasses>
  !!]
  type, abstract, extends(mergerTreeBuildMassesClass) :: mergerTreeBuildMassesRead
     !!{
     Implementation of a merger tree masses class which reads masses from a file.
     !!}
     private
     type            (varying_string) :: fileName
     double precision                 :: massIntervalFractional
   contains
     !![
     <methods>
       <method description="Read the halo masses, and, optionally, weights, from file." method="read"/>
     </methods>
     !!]
     procedure                             :: construct => readConstruct
     procedure(readReadTemplate), deferred :: read
  end type mergerTreeBuildMassesRead

  abstract interface
     subroutine readReadTemplate(self,mass,weight)
       import mergerTreeBuildMassesRead
       class           (mergerTreeBuildMassesRead), intent(inout)                            :: self
       double precision                           , intent(  out), allocatable, dimension(:) :: mass, weight
     end subroutine readReadTemplate
  end interface

contains

  subroutine readConstruct(self,time,mass,massMinimum,massMaximum,weight)
    !!{
    Construct a set of merger tree masses by reading from a file.
    !!}
    use :: Sorting, only : sort
    implicit none
    class           (mergerTreeBuildMassesRead), intent(inout)                            :: self
    double precision                           , intent(in   )                            :: time
    double precision                           , intent(  out), allocatable, dimension(:) :: mass          , weight      , &
         &                                                                                   massMinimum   , massMaximum
    double precision                                                                      :: massPrevious
    integer                                                                               :: iStart        , iEnd        , &
         &                                                                                   iStartPrevious, iEndPrevious, &
         &                                                                                   i
    !$GLC attributes unused :: time

    call self%read(mass,weight)
    if (allocated(weight)) then
       call sort(mass,weight)
    else
       call sort(mass       )
       allocate(massMinimum,mold=mass)
       allocate(massMaximum,mold=mass)
       massPrevious  =-huge(0.0d0)
       iStart        =-huge(0    )
       iEnd          =-huge(0    )
       iStartPrevious=-huge(0    )
       iEndPrevious  =-huge(0    )
       i             =+     0
       do while (i <= size(mass))
          i=i+1
          if (i == size(mass)+1 .or. mass(i) /= massPrevious) then
             ! Process a block of trees of the same mass.
             if (massPrevious > 0.0d0) then
                massMinimum(iStart:iEnd)=+mass(iStart)/sqrt(1.0d0+self%massIntervalFractional)
                massMaximum(iStart:iEnd)=+mass(iStart)*sqrt(1.0d0+self%massIntervalFractional)
                if (iStart > 1 .and. massMinimum(iStart) < massMaximum(iStartPrevious)) then
                   massMinimum(iStart        :iEnd        )=+sqrt(                             &
                        &                                         +mass       (iStartPrevious) &
                        &                                         *mass       (iStart        ) &
                        &                                        )
                   massMaximum(iStartPrevious:iEndPrevious)=+      massMinimum(iStart        )
                end if
             end if
             ! Update to the start of the next block.
             iStartPrevious=     iStart
             iEndPrevious  =     iEnd
             iStart        =     i
             massPrevious  =mass(i     )
          end if
          iEnd=i
       end do
    end if
    return
  end subroutine readConstruct
