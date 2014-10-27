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
     double precision                                      :: radialVelocity                      , tangentialVelocity, &
          &                                                   virialDensityContrastValue
     integer                                               :: virialDensityContrast
     class           (virialDensityContrastClass), pointer :: densityContrast
     logical                                               :: densityContrastInitialized
   contains
     final     ::                              fixedDestructor
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

  ! Virial density contrast types.
  integer, parameter :: fixedDensityContrastSphericalCollapseMatterLambda=0
  integer, parameter :: fixedDensityContrastFixedCritical                =1
  integer, parameter :: fixedDensityContrastFixedMean                    =2
  integer, parameter :: fixedDensityContrastBryanNorman1998              =3

contains

  function fixedDefaultConstructor()
    !% Default constructor for the {\tt fixed} dark matter halo virial density contrast class.
    use Input_Parameters
    implicit none
    type(virialOrbitFixed), target :: fixedDefaultConstructor
    
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
    use Galacticus_Error
    implicit none
    type            (virialOrbitFixed), target        :: fixedConstructor
    double precision                  , intent(in   ) :: radialVelocity       , tangentialVelocity
    type            (varying_string  ), intent(in   ) :: virialDensityContrast
    character       (len=32          )                :: label

    fixedConstructor%radialVelocity            =radialVelocity
    fixedConstructor%tangentialVelocity        =tangentialVelocity
    fixedConstructor%densityContrastInitialized=.false.
    if      (        virialDensityContrast       == 'sphericalCollapseMatterLambda') then
       fixedConstructor%virialDensityContrast=fixedDensityContrastSphericalCollapseMatterLambda
    else if (extract(virialDensityContrast,1, 9) == "fixedMean"                    ) then
       fixedConstructor%virialDensityContrast=fixedDensityContrastFixedMean
       label=extract(virialDensityContrast,10,len(virialDensityContrast))
       read (label,*) fixedConstructor%virialDensityContrastValue
    else if (extract(virialDensityContrast,1,13) == "fixedCritical"                ) then
       fixedConstructor%virialDensityContrast=fixedDensityContrastFixedCritical
       label=extract(virialDensityContrast,14,len(virialDensityContrast))
       read (label,*) fixedConstructor%virialDensityContrastValue
    else if (        virialDensityContrast       == "bryanNorman1998"              ) then
       fixedConstructor%virialDensityContrast=fixedDensityContrastBryanNorman1998
    else
       call Galacticus_Error_Report('fixedDensityContrastDefinition','only "sphericalCollapseMatterLambda", "bryanNorman1998", "fixedMeanXXX", and "fixedCriticalXXX" supported')
    end if
    return
  end function fixedConstructor

  subroutine fixedDestructor(self)
    !% Destructor for the {\tt fixed} virial orbits class.
    implicit none
    type(virialOrbitFixed), intent(inout) :: self

    if (self%densityContrastInitialized.and.self%densityContrast%isFinalizable()) deallocate(self%densityContrast)
    return
  end subroutine fixedDestructor

  function fixedOrbit(self,node,host,acceptUnboundOrbits)
    !% Return fixed orbital parameters for a satellite.
    use Galacticus_Error
    use Dark_Matter_Halo_Scales
    use Dark_Matter_Profile_Mass_Definitions
    implicit none
    type            (keplerOrbit               )                         :: fixedOrbit
    class           (virialOrbitFixed          ), intent(inout)          :: self
    type            (treeNode                  ), intent(inout), pointer :: host                      , node
    logical                                     , intent(in   )          :: acceptUnboundOrbits
    class           (nodeComponentBasic        )               , pointer :: hostBasic                 , basic
    class           (darkMatterHaloScaleClass  )               , pointer :: darkMatterHaloScale_
    class           (virialDensityContrastClass), pointer                :: virialDensityContrast_
    double precision                                                     :: velocityHost              , massHost          , &
         &                                                                  radiusHost                , massSatellite

    ! Reset the orbit.
    call fixedOrbit%reset()
    ! Get required objects.
    darkMatterHaloScale_ => darkMatterHaloScale()
    basic                => node%basic()
    hostBasic            => host%basic()
    ! Find virial density contrast under our definition.
    virialDensityContrast_ => self%densityContrastDefinition()
    ! Find mass, radius, and velocity in the host corresponding to the our virial density contrast definition.
    massHost     =Dark_Matter_Profile_Mass_Definition(host,virialDensityContrast_%densityContrast(hostBasic%mass(),hostBasic%time()),radiusHost,velocityHost)
    massSatellite=Dark_Matter_Profile_Mass_Definition(node,virialDensityContrast_%densityContrast(    basic%mass(),    basic%time())                        )
    if (virialDensityContrast_%isFinalizable()) deallocate(virialDensityContrast_)
    ! Set basic properties of the orbit.
    call fixedOrbit%massesSet(massSatellite,massHost)
    call fixedOrbit%radiusSet(radiusHost)
    call fixedOrbit%velocityRadialSet    (self%radialVelocity    *velocityHost)
    call fixedOrbit%velocityTangentialSet(self%tangentialVelocity*velocityHost)
    ! Propagate the orbit to the virial radius under the default density contrast definition.
    radiusHost=darkMatterHaloScale_%virialRadius(host)
    if (fixedOrbit%radiusApocenter() >= radiusHost .and. fixedOrbit%radiusPericenter() <= radiusHost) then
       call fixedOrbit%propagate(radiusHost  ,infalling=.true.)
       call fixedOrbit%massesSet(basic%mass(),hostBasic%mass())
    else
       call Galacticus_Error_Report('fixedOrbit','orbit does not reach halo radius')
    end if
    return
  end function fixedOrbit

  function fixedDensityContrastDefinition(self)
    !% Return a virial density contrast object defining that used in the definition of fixed virial orbits.
    use Galacticus_Error
    implicit none
    class  (virialDensityContrastClass), pointer       :: fixedDensityContrastDefinition
    class  (virialOrbitFixed          ), intent(inout) :: self
     
    if (.not.self%densityContrastInitialized) then
       select case (self%virialDensityContrast)
       case (fixedDensityContrastSphericalCollapseMatterLambda)
          allocate(virialDensityContrastSphericalCollapseMatterLambda :: self%densityContrast)
          select type (d => self%densityContrast)
          type is (virialDensityContrastSphericalCollapseMatterLambda)
             d=virialDensityContrastSphericalCollapseMatterLambda()
          end select
       case (fixedDensityContrastFixedCritical)
          allocate(virialDensityContrastFixed :: self%densityContrast)
          select type (d => self%densityContrast)
          type is (virialDensityContrastFixed)
             d=virialDensityContrastFixed(self%virialDensityContrastValue,virialDensityContrastFixedDensityTypeCritical)
          end select
       case (fixedDensityContrastFixedMean)
          allocate(virialDensityContrastFixed :: self%densityContrast)
          select type (d => self%densityContrast)
          type is (virialDensityContrastFixed)
             d=virialDensityContrastFixed(self%virialDensityContrastValue,virialDensityContrastFixedDensityTypeMean    )
          end select
       case (fixedDensityContrastBryanNorman1998)
          allocate(virialDensityContrastBryanNorman1998 :: self%densityContrast)
          select type (d => self%densityContrast)
          type is (virialDensityContrastBryanNorman1998)
             d=virialDensityContrastBryanNorman1998()
          end select
       end select
       call self%densityContrast%makeIndestructible()
       self%densityContrastInitialized=.true.
    end if
    fixedDensityContrastDefinition => self%densityContrast
    return
  end function fixedDensityContrastDefinition

