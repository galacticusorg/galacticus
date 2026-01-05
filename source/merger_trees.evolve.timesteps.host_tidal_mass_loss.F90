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

  use :: Satellite_Tidal_Stripping, only : satelliteTidalStrippingClass

  !![
  <mergerTreeEvolveTimestep name="mergerTreeEvolveTimestepHostTidalMassLoss">
   <description>
    A merger tree evolution timestepping class that enforces
    \begin{eqnarray}
    \Delta t &amp;\le&amp; \epsilon_\mathrm{hostTidalMassLoss} (M_\mathrm{host}/\dot{M}_\mathrm{host}|),
    \end{eqnarray}
    where $\epsilon_\mathrm{hostTidalMassLoss}=${\normalfont \ttfamily [timeStepRelative]}, and $M_\mathrm{host}$ is the
    bound mass of the host satellite. This criterion is intended to prevent any satellite evolving over an excessively
    large time in one step ahead of its host satellite.
   </description>
  </mergerTreeEvolveTimestep>
  !!]
  type, extends(mergerTreeEvolveTimestepClass) :: mergerTreeEvolveTimestepHostTidalMassLoss
     !!{
     Implementation of a merger tree evolution timestepping class which is based on the mass loss rate of the satellite.
     !!}
     private
     class           (satelliteTidalStrippingClass), pointer :: satelliteTidalStripping_ => null()
     double precision                                        :: timeStepRelative                  , fractionTimestepMinimum
     logical                                                 :: refuseToEvolve_
   contains
     final     ::                   hostTidalMassLossDestructor
     procedure :: timeEvolveTo   => hostTidalMassLossTimeEvolveTo
     procedure :: refuseToEvolve => hostTidalMassLossRefuseToEvolve
  end type mergerTreeEvolveTimestepHostTidalMassLoss

  interface mergerTreeEvolveTimestepHostTidalMassLoss
     !!{
     Constructors for the \refClass{mergerTreeEvolveTimestepHostTidalMassLoss} merger tree evolution timestep class.
     !!}
     module procedure hostTidalMassLossConstructorParameters
     module procedure hostTidalMassLossConstructorInternal
  end interface mergerTreeEvolveTimestepHostTidalMassLoss

contains

  function hostTidalMassLossConstructorParameters(parameters) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepHostTidalMassLoss} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (mergerTreeEvolveTimestepHostTidalMassLoss)                :: self
    type            (inputParameters                          ), intent(inout) :: parameters
    class           (satelliteTidalStrippingClass             ), pointer       :: satelliteTidalStripping_
    double precision                                                           :: timeStepRelative        , fractionTimestepMinimum

    !![
    <inputParameter>
      <name>timeStepRelative</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The maximum allowed relative change in time for a single step in the evolution of a node.</description>
      <source>parameters</source>
    </inputParameter>
    <inputParameter>
      <name>fractionTimestepMinimum</name>
      <defaultValue>0.1d0</defaultValue>
      <description>The minimum fraction of the timestep imposed by this timestepper to evolve over. If the timestep allowed is smaller than this fraction, the actual timestep will be reduced to zero. This avoids forcing satellites to take a large number of very small timesteps, and instead defers evolving a satellite until a large timestep can be taken.</description>
      <source>parameters</source>
    </inputParameter>
    <objectBuilder class="satelliteTidalStripping" name="satelliteTidalStripping_" source="parameters"/>
    !!]
    self=mergerTreeEvolveTimestepHostTidalMassLoss(timeStepRelative,fractionTimestepMinimum,satelliteTidalStripping_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="satelliteTidalStripping_"/>
    !!]
    return
  end function hostTidalMassLossConstructorParameters

  function hostTidalMassLossConstructorInternal(timeStepRelative,fractionTimestepMinimum,satelliteTidalStripping_) result(self)
    !!{
    Constructor for the \refClass{mergerTreeEvolveTimestepHostTidalMassLoss} merger tree evolution timestep class which takes a parameter set as input.
    !!}
    implicit none
    type            (mergerTreeEvolveTimestepHostTidalMassLoss)                        :: self
    double precision                                           , intent(in   )         :: timeStepRelative        , fractionTimestepMinimum
    class           (satelliteTidalStrippingClass             ), intent(in   ), target :: satelliteTidalStripping_
    !![
    <constructorAssign variables="timeStepRelative, fractionTimestepMinimum, *satelliteTidalStripping_"/>
    !!]

    return
  end function hostTidalMassLossConstructorInternal

  subroutine hostTidalMassLossDestructor(self)
    !!{
    Destructor for the \refClass{mergerTreeEvolveTimestepHostTidalMassLoss} merger tree evolution timestep class.
    !!}
    implicit none
    type(mergerTreeEvolveTimestepHostTidalMassLoss), intent(inout) :: self

    !![
    <objectDestructor name="self%satelliteTidalStripping_"/>
    !!]
    return
  end subroutine hostTidalMassLossDestructor

  double precision function hostTidalMassLossTimeEvolveTo(self,timeEnd,node,task,taskSelf,report,lockNode,lockType)
    !!{
    Determine a suitable timestep for {\normalfont \ttfamily node} using the {\normalfont \ttfamily hostTidalMassLoss} method.
    This sets the time step size to {\normalfont \ttfamily timeStepRelative}$M_\mathrm{host}/|\dot{M}_\mathrm{host}|$.
    !!}
    use :: Evolve_To_Time_Reports, only : Evolve_To_Time_Report
    use :: Galacticus_Nodes      , only : nodeComponentBasic   , nodeComponentSatellite, treeNode
    use :: ISO_Varying_String    , only : varying_string
    implicit none
    class           (mergerTreeEvolveTimestepHostTidalMassLoss), intent(inout), target            :: self
    double precision                                           , intent(in   )                    :: timeEnd
    type            (treeNode                                 ), intent(inout), target            :: node
    procedure       (timestepTask                             ), intent(  out), pointer           :: task
    class           (*                                        ), intent(  out), pointer           :: taskSelf
    logical                                                    , intent(in   )                    :: report
    type            (treeNode                                 ), intent(  out), pointer, optional :: lockNode
    type            (varying_string                           ), intent(  out)         , optional :: lockType
    type            (treeNode                                 ), pointer                          :: nodeHost
    class           (nodeComponentBasic                       )               , pointer           :: basic            , basicHost
    class           (nodeComponentSatellite                   )               , pointer           :: satellite
    double precision                                                                              :: boundMass        , massLossRate, &
         &                                                                                           timescaleMassLoss
    !$GLC attributes unused :: timeEnd

    self%refuseToEvolve_=.false.
    hostTidalMassLossTimeEvolveTo=huge(0.0d0)
    if (node%isSatellite()) then
       nodeHost  => node    %parent
       basicHost => nodeHost%basic ()
       basic     => node%basic     ()
       if (nodeHost%isSatellite()) then
          ! The host halo is currently a subhalo. Find the current mass loss timescale of the host.
          satellite    =>     nodeHost                          %satellite   (        )
          boundMass    =      satellite                         %boundMass   (        )
          massLossRate =  abs(self     %satelliteTidalStripping_%massLossRate(nodeHost))
          if (self%timeStepRelative > 0.0d0 .and. boundMass > 0.0d0 .and. massLossRate > 0.0d0) then
             timescaleMassLoss            =self%timeStepRelative*boundMass/massLossRate
             hostTidalMassLossTimeEvolveTo=timescaleMassLoss+basicHost%time()
             ! If the allowed timestep is too small, we will refuse to evolve, allowing the host halo to evolve instead.
             if (hostTidalMassLossTimeEvolveTo-basic%time() < self%fractionTimestepMinimum*timescaleMassLoss) self%refuseToEvolve_=.true.
          end if
       else
          ! The host halo is currently not a subhalo. If it will ever become a subhalo, then limit the evolution time of this
          ! node to the time at which the host becomes a subhalo. This prevents this subhalo from evolving in advance of the
          ! host after it may experience tidal stripping.
          do while (associated(nodeHost))
             if (.not.nodeHost%isPrimaryProgenitor() .and. associated(nodeHost%parent)) then
                ! Limit the time to that of the halo where the current host will become a subhalo. Include a small offset here to
                ! avoid deadlocking the tree.
                nodeHost                      => nodeHost %parent
                basicHost                     => nodeHost %basic ()
                hostTidalMassLossTimeEvolveTo =  basicHost%time  ()+1.0d-6
                exit
             end if
             nodeHost => nodeHost%parent
          end do
       end if
       ! Prevent the target time being earlier than the current time of the node.
       hostTidalMassLossTimeEvolveTo =  max(hostTidalMassLossTimeEvolveTo,basic%time())
    end if
    task                               => null()
    taskSelf                           => null()
    if (present(lockNode)) lockNode    => node
    if (present(lockType)) lockType    =  "hostTidalMassLoss"
    if (        report   ) call Evolve_To_Time_Report("hostTidalMassLoss: ",hostTidalMassLossTimeEvolveTo)
    return
  end function hostTidalMassLossTimeEvolveTo

  logical function hostTidalMassLossRefuseToEvolve(self,node)
    !!{
    Refuse to evolve if the timestep is too small.
    !!}
    implicit none
    class(mergerTreeEvolveTimestepHostTidalMassLoss), intent(inout) :: self
    type (treeNode                                 ), intent(inout) :: node
    !$GLC attributes unused :: node

    hostTidalMassLossRefuseToEvolve=self%refuseToEvolve_
    return
  end function hostTidalMassLossRefuseToEvolve
