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

!+    Contributions to this file made by: Claude.

!!{RST
Contains a module which provides utilities shared by spherical mass distributions.
!!}

module Mass_Distributions_Spherical_Utilities
  !!{RST
  Provides utilities shared by spherical mass distributions.
  !!}
  implicit none
  private
  public :: Mass_Distribution_Time_Freefall_Tabulate

  abstract interface
     double precision function potentialDifferenceScaleFreeTemplate(radius1,radius2)
       !!{RST
       Interface for functions returning the potential difference between two radii in a scale-free mass distribution.
       !!}
       double precision, intent(in   ) :: radius1, radius2
     end function potentialDifferenceScaleFreeTemplate
  end interface

contains

  subroutine Mass_Distribution_Time_Freefall_Tabulate(timeFreefallScaleFree_,timeScaleFree,potentialDifferenceScaleFree,toleranceRelative)
    !!{RST
    Extend the tabulation of the freefall time in a scale-free mass distribution, ``timeFreefallScaleFree_``, until it brackets
    ``timeScaleFree``. The freefall time from each radius is found by integrating over the potential difference given by
    ``potentialDifferenceScaleFree``, to relative tolerance ``toleranceRelative``.
    !!}
    use :: Numerical_Integration, only : integrator
    use :: Tabulations_Inverse  , only : tabulationInverse
    implicit none
    type            (tabulationInverse                   ), intent(inout)               :: timeFreefallScaleFree_
    double precision                                      , intent(in   )               :: timeScaleFree               , toleranceRelative
    procedure       (potentialDifferenceScaleFreeTemplate)                              :: potentialDifferenceScaleFree
    double precision                                      , allocatable  , dimension(:) :: radii
    double precision                                                                    :: radiusStart
    integer                                                                             :: i
    type            (integrator                          )                              :: integrator_

    ! Each point is an independent quadrature from the center out to its own radius, so a point carried over by an extension
    ! is precisely the value which would be computed afresh.
    if (.not.timeFreefallScaleFree_%brackets(timeScaleFree)) then
       integrator_=integrator(timeFreeFallIntegrand,toleranceRelative=toleranceRelative)
       do while (.not.timeFreefallScaleFree_%brackets(timeScaleFree))
          call timeFreefallScaleFree_%expand(timeScaleFree)
          radii=timeFreefallScaleFree_%abscissae()
          do i=1,size(radii)
             if (timeFreefallScaleFree_%isComputed(i)) cycle
             call timeFreefallScaleFree_%set(i,timeFreefallScaleFree(radii(i)))
          end do
          call timeFreefallScaleFree_%build()
       end do
    end if
    return

  contains

    double precision function timeFreefallScaleFree(radius)
      !!{RST
      Evaluate the freefall time from a given radius in the scale-free mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius

      radiusStart          =                            radius
      timeFreefallScaleFree=integrator_%integrate(0.0d0,radius)
      return
    end function timeFreefallScaleFree

    double precision function timeFreeFallIntegrand(radius)
      !!{RST
      Integrand used to find the freefall time in the scale-free mass distribution.
      !!}
      implicit none
      double precision, intent(in   ) :: radius
      double precision                :: potentialDifference

      if (radius == 0.0d0) then
         timeFreeFallIntegrand=+0.0d0
      else
         potentialDifference=+potentialDifferenceScaleFree(radiusStart,radius)
         if (potentialDifference > 0.0d0) then
            timeFreeFallIntegrand=+1.0d0                     &
                 &                /sqrt(                     &
                 &                      +2.0d0               &
                 &                      *potentialDifference &
                 &                     )
         else
            timeFreeFallIntegrand=+0.0d0
         end if
      end if
      return
    end function timeFreeFallIntegrand

  end subroutine Mass_Distribution_Time_Freefall_Tabulate

end module Mass_Distributions_Spherical_Utilities
