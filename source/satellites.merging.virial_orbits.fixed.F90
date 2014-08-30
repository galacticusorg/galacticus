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
     double precision                 :: radialVelocity       , tangentialVelocity
     type            (varying_string) :: virialDensityContrast
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
  double precision                 :: fixedRadialVelocity       , fixedTangentialVelocity
  type            (varying_string) :: fixedVirialDensityContrast

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
          !@ <inputParameter>
          !@   <name>virialOrbitsVirialDensityContrast</name>
          !@   <defaultValue>sphericalCollapseMatterLambda</defaultValue>
          !@   <attachedTo>module</attachedTo>
          !@   <description>
          !@     The virial density contrast to assume for the fixed virial orbits class.
          !@   </description>
          !@   <type>real</type>
          !@   <cardinality>1</cardinality>
          !@ </inputParameter>
          call Get_Input_Parameter('virialOrbitsVirialDensityContrast',fixedVirialDensityContrast,defaultValue='sphericalCollapseMatterLambda')
          ! Record initialization.
          fixedInitialized=.true.
       end if
       !$omp end critical(virialOrbitFixedInitialize)
    end if
    fixedDefaultConstructor=fixedConstructor(fixedRadialVelocity,fixedTangentialVelocity,fixedVirialDensityContrast)
    return
  end function fixedDefaultConstructor

  function fixedConstructor(radialVelocity,tangentialVelocity,virialDensityContrast)
    !% Generic constructor for the {\tt fixed} virial orbits class.
    use Input_Parameters
    implicit none
    type            (virialOrbitFixed), target        :: fixedConstructor
    double precision                  , intent(in   ) :: radialVelocity       , tangentialVelocity
    type            (varying_string  ), intent(in   ) :: virialDensityContrast

    fixedConstructor%radialVelocity       =radialVelocity
    fixedConstructor%tangentialVelocity   =tangentialVelocity
    fixedConstructor%virialDensityContrast=virialDensityContrast
    return
  end function fixedConstructor

  function fixedOrbit(self,node,host,acceptUnboundOrbits)
    !% Return fixed orbital parameters for a satellite.
    use Galacticus_Error
    use Cosmology_Parameters
    use Galactic_Structure_Enclosed_Masses
    use Galactic_Structure_Options
    use Numerical_Constants_Physical
    use Dark_Matter_Halo_Scales
    implicit none
    type            (keplerOrbit               )                         :: fixedOrbit
    class           (virialOrbitFixed          ), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: host                      , node
    logical                                     , intent(in   )          :: acceptUnboundOrbits
    class           (nodeComponentBasic        )               , pointer :: hostBasic                 , basic
    class           (darkMatterHaloScaleClass  )               , pointer :: darkMatterHaloScale_
    class           (virialDensityContrastClass), pointer                :: virialDensityContrast_
    class           (cosmologyParametersClass  ), pointer                :: cosmologyParameters_
    double precision                                                     :: velocityHost              , massHost          , &
         &                                                                  radiusHost                , darkMatterFraction, &
         &                                                                  fixedVirialDensityContrast

    ! Reset the orbit.
    call fixedOrbit%reset()
    ! Get required objects.
    cosmologyParameters_ => cosmologyParameters()
    darkMatterHaloScale_ => darkMatterHaloScale()
    ! Find virial density contrast under Benson (2005) definition.
    darkMatterFraction         =  (1.0d0-cosmologyParameters_%omegaBaryon()/cosmologyParameters_%omegaMatter())
    virialDensityContrast_     => self                  %densityContrastDefinition(                )
    fixedVirialDensityContrast =  virialDensityContrast_%densityContrast          (hostBasic%time())*darkMatterFraction
    deallocate(virialDensityContrast_)
    ! Find radius in the host corresponding to the Benson (2005) virial density contrast definition.
    radiusHost  =Galactic_Structure_Radius_Enclosing_Density(host,densityContrast=fixedVirialDensityContrast,massType=massTypeDark,haloLoaded=.false.)
    massHost    =Galactic_Structure_Enclosed_Mass           (host,radiusHost                                ,massType=massTypeDark,haloLoaded=.false.)/darkMatterFraction
    velocityHost=sqrt(gravitationalConstantGalacticus*massHost/radiusHost)
    ! Set basic properties of the orbit.
    basic     => node%basic()
    hostBasic => host%basic()
    call fixedOrbit%massesSet(basic%mass(),hostBasic%mass())
    call fixedOrbit%radiusSet(radiusHost)
    call fixedOrbit%velocityRadialSet    (self%radialVelocity    *velocityHost)
    call fixedOrbit%velocityTangentialSet(self%tangentialVelocity*velocityHost)
    ! Propagate the orbit to the virial radius under the default density contrast definition.
    radiusHost=darkMatterHaloScale_%virialRadius(host)
    if (fixedOrbit%radiusApocenter() >= radiusHost) then
       call fixedOrbit%propagate(radiusHost,infalling=.true.)
    else
       call Galacticus_Error_Report('fixedOrbit','orbit does not reach halo radius')
    end if
    return
  end function fixedOrbit

  function fixedDensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of fixed virial orbits.
    use Galacticus_Error
    implicit none
    class(virialDensityContrastClass), pointer       :: fixedDensityContrastDefinition
    class(virialOrbitFixed          ), intent(inout) :: self
    
    select case (char(self%virialDensityContrast))
    case ('sphericalCollapseMatterLambda')
       allocate(virialDensityContrastSphericalCollapseMatterLambda :: fixedDensityContrastDefinition)
       select type (fixedDensityContrastDefinition)
       type is (virialDensityContrastSphericalCollapseMatterLambda)
          fixedDensityContrastDefinition=virialDensityContrastSphericalCollapseMatterLambda()
       end select
    case default
       call Galacticus_Error_Report('fixedDensityContrastDefinition','only "sphericalCollapseMatterLambda" supported')
    end select
    return
  end function fixedDensityContrastDefinition
  
