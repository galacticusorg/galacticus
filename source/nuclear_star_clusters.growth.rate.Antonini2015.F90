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

  !+    Contributions to this file made by: Matías Liempi

  !!{
  Implementation of the \cite{antonini_coevolution_2015} star formation rate law for galactic \glspl{nsc}.
  !!}
  
  use :: Star_Formation_Rates_Spheroids, only : starFormationRateSpheroidsClass

  !![
  <nuclearStarClusterGrowthRates name="nuclearStarClusterGrowthRatesAntonini2015">
   <description>
    A gas inflow rate implementing the model of \citep{antonini_coevolution_2015} for galactic \glspl{nsc}.
    \begin{equation}
     \dot{M}_\mathrm{gas}^\mathrm{NSC} = A_\mathrm{res}\dot{M}_\star^\mathrm{spheroid},
    \end{equation}    
    where $A_\mathrm{res}=${\normalfont \ttfamily [efficiency]} is a free parameter, and $\dot{M}_\star^\mathrm{spheroid}$ is the
    star formation rate of the spheroid component.
   </description>
  </nuclearStarClusterGrowthRates>
  !!]
  type, extends(nuclearStarClusterGrowthRatesClass) :: nuclearStarClusterGrowthRatesAntonini2015
     !!{
     Implementation of the \cite{antonini_coevolution_2015} gas inflow rate for galactic \glspl{nsc}.
     !!}
     private
     class           (starFormationRateSpheroidsClass), pointer :: starFormationRateSpheroids_ => null()
     double precision                                           :: efficiency
     contains
     final     ::          antonini2015Destructor
     procedure :: rate  => antonini2015Rate
  end type nuclearStarClusterGrowthRatesAntonini2015

  interface nuclearStarClusterGrowthRatesAntonini2015
     !!{
     Constructors for the \refClass{nuclearStarClusterGrowthRatesAntonini2015} gas inflow rate in \glspl{nsc} class.
     !!}
     module procedure antonini2015ConstructorParameters
     module procedure antonini2015ConstructorInternal
  end interface nuclearStarClusterGrowthRatesAntonini2015
    
contains

  function antonini2015ConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{nuclearStarClusterGrowthRatesAntonini2015} gas inflow rate in \glspl{nsc} class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (nuclearStarClusterGrowthRatesAntonini2015)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (starFormationRateSpheroidsClass          ), pointer       :: starFormationRateSpheroids_
    double precision                                                           :: efficiency       

    !![
    <inputParameter>
      <name>efficiency</name>
      <defaultSource>\citep{antonini_coevolution_2015}</defaultSource>
      <defaultValue>1.0d-2</defaultValue>
      <description>Parameter controlling the rate of the gas inflow onto the nuclear star cluster.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="starFormationRateSpheroids" name="starFormationRateSpheroids_" source="parameters"/>
    !!]
    self=nuclearStarClusterGrowthRatesAntonini2015(efficiency,starFormationRateSpheroids_)
    !![
    <inputParametersValidate source="parameters"        />
    <objectDestructor name="starFormationRateSpheroids_"/>
    !!]
    return
  end function antonini2015ConstructorParameters

  function antonini2015ConstructorInternal(efficiency,starFormationRateSpheroids_) result(self)
    !!{
    Internal constructor for the \refClass{nuclearStarClusterGrowthRatesAntonini2015} gas inflow rate from NSCs class.
    !!}
    implicit none
    type            (nuclearStarClusterGrowthRatesAntonini2015)                        :: self
    class           (starFormationRateSpheroidsClass          ), intent(in   ), target :: starFormationRateSpheroids_
    double precision                                           , intent(in   )         :: efficiency
    !![
    <constructorAssign variables="efficiency, *starFormationRateSpheroids_"/>
    !!]
    return
  end function antonini2015ConstructorInternal

  subroutine antonini2015Destructor(self)
    !!{
    Destructor for the \refClass{nuclearStarClusterGrowthRatesAntonini2015} class
    !!}
    implicit none
    type(nuclearStarClusterGrowthRatesAntonini2015), intent(inout) :: self
    !![
    <objectDestructor name="self%starFormationRateSpheroids_"/>
    !!]
    return
  end subroutine antonini2015Destructor

  double precision function antonini2015Rate(self,node) result(rate)
    !!{
    Returns the gas inflow rate (in $M_\odot$ Gyr$^{-1}$) onto the galactic \gls{nsc} of {\normalfont \ttfamily
    node}. The \gls{nsc} is assumed to obey the \cite{antonini_coevolution_2015} gas inflow rate model.
    !!}
    use :: Galacticus_Nodes, only : nodeComponentNSC, nodeComponentSpheroid
    implicit none
    class           (nuclearStarClusterGrowthRatesAntonini2015), intent(inout), target  :: self
    type            (treeNode                                 ), intent(inout)          :: node
    class           (nodeComponentSpheroid                    ),                pointer :: spheroid
    double precision                                                                    :: rateStarFormationSpheroid

    ! Get the spheroid component.
    spheroid                  => node                            %spheroid(    )
    ! Get the star formation rate of the spheroid component.
    rateStarFormationSpheroid =  self%starFormationRateSpheroids_%rate    (node)  
    ! Find the rate of mass accretion onto the nuclear star cluster.
    if (rateStarFormationSpheroid <= 0.0d0) then
      rate    =+0.0d0
    else
      ! Gas accretion rate model from F. Antonini, E. Barausse & J. Silk (2015; https://ui.adsabs.harvard.edu/abs/2015ApJ...812...72A).
       rate   =+self%efficiency                &
            &  *     rateStarFormationSpheroid
    end if 
    return
  end function antonini2015Rate

  
