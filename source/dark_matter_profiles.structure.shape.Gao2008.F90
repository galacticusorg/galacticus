!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015 Andrew Benson <abenson@obs.carnegiescience.edu>
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

  !% An implementation of dark matter halo profile shapes  using the \cite{gao_redshift_2008} algorithm.

  !# <darkMatterProfileShape name="darkMatterProfileShapeGao2008">
  !#  <description>Dark matter halo shape parameters are computed using the algorithm of \cite{gao_redshift_2008}.</description>
  !# </darkMatterProfileShape>

  type, extends(darkMatterProfileShapeClass) :: darkMatterProfileShapeGao2008
     !% A dark matter halo profile shape parameter class implementing the algorithm of \cite{gao_redshift_2008}.
     private
   contains
     procedure :: shape => gao2008Shape
  end type darkMatterProfileShapeGao2008

  interface darkMatterProfileShapeGao2008
     !% Constructors for the {\normalfont \ttfamily gao2008} dark matter halo profile shape parameter class.
     module procedure gao2008DefaultConstructor
  end interface darkMatterProfileShapeGao2008

contains

  function gao2008DefaultConstructor()
    !% Default constructor for the {\normalfont \ttfamily gao2008} dark matter halo profile shape parameter class.
    use Input_Parameters
    implicit none
    type(darkMatterProfileShapeGao2008), target :: gao2008DefaultConstructor

    return
  end function gao2008DefaultConstructor
  
  double precision function gao2008Shape(self,node)
    !% Return the Einasto profile shape parameter of the dark matter halo profile of {\normalfont \ttfamily node} using the
    !% \cite{gao_redshift_2008} algorithm. More specifically, the parameter is given by:
    !% \begin{equation}
    !% \alpha = \left\{ \begin{array}{ll} 0.155 + 0.0095\nu^2 & \hbox{ if } \nu < 3.907 \\ 0.3 & \hbox{ if } \nu \ge 3.907, \end{array} \right.
    !% \end{equation}
    !% where $\nu=\delta_{\mathrm c}(t)/\sigma(M)$ is the peak height of the halo.
    use Galacticus_Nodes
    use Power_Spectra
    use Critical_Overdensities
    implicit none
    class           (darkMatterProfileShapeGao2008), intent(inout)          :: self
    type            (treeNode                     ), intent(inout), pointer :: node
    double precision                               , parameter              :: nuMaximum=3.907d0
    class           (criticalOverdensityClass     )               , pointer :: criticalOverdensity_
    class           (nodeComponentBasic           )               , pointer :: basic
    double precision                                                        :: nu
    
    ! Get default objects.
    criticalOverdensity_ => criticalOverdensity()
    ! Get the basic component.
    basic => node%basic()
    ! Compute the shape parameter.
    nu     =+criticalOverdensity_%value     (time=basic%time(),mass=basic%mass()) &
         &  /Cosmological_Mass_Root_Variance(                       basic%mass())
    if (nu < nuMaximum) then
       gao2008Shape=0.155d0+0.0095d0*nu**2
    else
       gao2008Shape=0.300d0
    end if    
    return
  end function gao2008Shape
