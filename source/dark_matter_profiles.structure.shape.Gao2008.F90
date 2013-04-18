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

!% Contains a module which implements the \cite{gao_redshift_2008} Einasto halo shape algorithm.

module Dark_Matter_Profiles_Shapes_Gao2008
  !% Implements the \cite{gao_redshift_2008} Einasto halo shape algorithm.
  implicit none
  private
  public :: Dark_Matter_Shapes_Gao2008_Initialize

contains

  !# <darkMatterShapeMethod>
  !#  <unitName>Dark_Matter_Shapes_Gao2008_Initialize</unitName>
  !# </darkMatterShapeMethod>
  subroutine Dark_Matter_Shapes_Gao2008_Initialize(darkMatterShapeMethod,Dark_Matter_Profile_Shape_Get)
    !% Initializes the ``Gao2008'' halo shape module.
    use ISO_Varying_String
    implicit none
    type(varying_string),                 intent(in)    :: darkMatterShapeMethod
    procedure(Dark_Matter_Profile_Shape_Gao2008), pointer, intent(inout) :: Dark_Matter_Profile_Shape_Get
    
    if (darkMatterShapeMethod == 'Gao2008') Dark_Matter_Profile_Shape_Get => Dark_Matter_Profile_Shape_Gao2008
  
    return
  end subroutine Dark_Matter_Shapes_Gao2008_Initialize

  double precision function Dark_Matter_Profile_Shape_Gao2008(thisNode)
    !% Returns the Einasto shape parameter, $alpha$, of the dark matter profile of {\tt thisNode} using the method of
    !% \cite{gao_redshift_2008}. More specifically, the parameter is given by:
    !% \begin{equation}
    !% \alpha = \left\{ \begin{array}{ll} 0.155 + 0.0095\nu^2 & \hbox{ if } \nu < 3.907 \\ 0.3 & \hbox{ if } \nu \ge 3.907, \end{array} \right.
    !% \end{equation}
    !% where $\nu=\delta_{\rm c}(t)/\sigma(M)$ is the peak height of the halo.
    use Galacticus_Nodes
    use Power_Spectra
    use Critical_Overdensity
    implicit none
    type (treeNode          ), intent(inout), pointer :: thisNode
    double precision         , parameter              :: nuMaximum=3.907d0
    class(nodeComponentBasic),                pointer :: thisBasicComponent
    double precision                                  :: nu
    
    ! Get the basic component.
    thisBasicComponent => thisNode%basic()
    ! Compute the shape parameter.
    nu=Critical_Overdensity_for_Collapse(time=thisBasicComponent%time(),mass=thisBasicComponent%mass())/Cosmological_Mass_Root_Variance(thisBasicComponent%mass())
    if (nu < nuMaximum) then
       Dark_Matter_Profile_Shape_Gao2008=0.155d0+0.0095d0*nu**2
    else
       Dark_Matter_Profile_Shape_Gao2008=0.3d0
    end if
    return
  end function Dark_Matter_Profile_Shape_Gao2008
  
end module Dark_Matter_Profiles_Shapes_Gao2008
