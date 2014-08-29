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

  !% An implementation of virial orbits which assumes fixed orbital parameters.

  !# <virialOrbit name="virialOrbitFixed">
  !#  <description>Virial orbits assuming fixed orbital parameters.</description>
  !# </virialOrbit>

  type, extends(virialOrbitClass) :: virialOrbitFixed
     !% A virial orbit class that assumes fixed orbital parameters.
     private
     double precision :: radialVelocity, tangentialVelocity
   contains
     procedure :: orbit                     => fixedOrbit
     procedure :: densityContrastDefinition => fixedDensityContrastDefinition
  end type virialOrbitFixed

  interface virialOrbitFixed
     !% Constructors for the {\tt fixed} virial orbit class.
     module procedure fixedDefaultConstructor
     module procedure fixedConstructor
  end interface virialOrbitFixed

  ! Initialization state.
  logical          :: fixedInitialized   =.false.

  ! Default orbital parameters.
  double precision :: fixedRadialVelocity        , fixedTangentialVelocity

contains

  function fixedDefaultConstructor()
    !% Default constructor for the {\tt fixed} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type (virialOrbitFixed), target  :: fixedDefaultConstructor
    
    if (.not.fixedInitialized) then
       !$omp critical(virialOrbitFixedInitialize)
       if (.not.fixedInitialized) then
          ! Get the parameters to use.
          !@ <inputParameter>
          !@   <name>virialOrbitsFixedRadialVelocity</name>
          !@   <defaultValue>0.90</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The radial velocity (in units of the host virial velocity) to used for the fixed virial orbits distribution. Default value matches approximate peak in the distribution of \cite{benson_orbital_2005}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('virialOrbitsFixedRadialVelocity'    ,fixedRadialVelocity    ,defaultValue=-0.90d0)
          !@ <inputParameter>
          !@   <name>virialOrbitsFixedTangentialVelocity</name>
          !@   <defaultValue>0.75</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The radial velocity (in units of the host virial velocity) to used for the fixed virial orbits distribution. Default value matches approximate peak in the distribution of \cite{benson_orbital_2005}.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('virialOrbitsFixedTangentialVelocity',fixedTangentialVelocity,defaultValue= 0.75d0)
          ! Record initialization.
          fixedInitialized=.true.
       end if
       !$omp end critical(virialOrbitFixedInitialize)
    end if
    fixedDefaultConstructor=fixedConstructor(fixedRadialVelocity,fixedTangentialVelocity)
    return
  end function fixedDefaultConstructor

  function fixedConstructor(radialVelocity,tangentialVelocity)
    !% Generic constructor for the {\tt fixed} virial orbits class.
    use Input_Parameters
    implicit none
    type            (virialOrbitFixed), target        :: fixedConstructor
    double precision                  , intent(in   ) :: radialVelocity  , tangentialVelocity

    fixedConstructor%radialVelocity    =radialVelocity
    fixedConstructor%tangentialVelocity=tangentialVelocity
    return
  end function fixedConstructor

  function fixedOrbit(self,node,host,acceptUnboundOrbits)
    !% Return fixed orbital parameters for a satellite.
    use Dark_Matter_Halo_Scales
    implicit none
    type            (keplerOrbit             )                         :: fixedOrbit
    class           (virialOrbitFixed        ), intent(inout)          :: self
    type            (treeNode                ), intent(inout), pointer :: host                , node
    logical                                   , intent(in   )          :: acceptUnboundOrbits
    class           (nodeComponentBasic      )               , pointer :: hostBasic           , basic
    type            (darkMatterHaloScaleClass)               , pointer :: darkMatterHaloScale_
    double precision                                                   :: velocityScale

    ! Reset the orbit.
    call fixedOrbit%reset()
    ! Get required objects.
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Set masses and radius of the orbit.
    basic     => node%basic()
    hostBasic => host%basic()
    call fixedOrbit%massesSet(basic%mass(),hostBasic%mass())
    call fixedOrbit%radiusSet(darkMatterHaloScale_%virialRadius(host))
    velocityScale=darkMatterHaloScale_%virialVelocity(host)
    call fixedOrbit%velocityRadialSet    (self%radialVelocity    *velocityScale)
    call fixedOrbit%velocityTangentialSet(self%tangentialVelocity*velocityScale)
    return
  end function fixedOrbit

  function fixedDensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of fixed virial orbits.
    implicit none
    class(virialDensityContrastClass), pointer       :: fixedDensityContrastDefinition
    class(virialOrbitFixed          ), intent(inout) :: self
    
    fixedDensityContrastDefinition => virialDensityContrast()
    return
  end function fixedDensityContrastDefinition
  
