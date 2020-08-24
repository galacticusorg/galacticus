!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020
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

!% Contains a module which provides IO in \gls{irate} format.

module IO_IRATE
  !% Provides IO in \gls{irate} format.
  use :: Cosmology_Functions , only : cosmologyFunctionsClass
  use :: Cosmology_Parameters, only : cosmologyParametersClass
  use :: ISO_Varying_String  , only : varying_string
  private
  public :: irate

  type :: irate
     !% A class for interacting with \gls{irate} format files.
     type (varying_string          )          :: fileName
     class(cosmologyFunctionsClass ), pointer :: cosmologyFunctions_
     class(cosmologyParametersClass), pointer :: cosmologyParameters_
   contains
     !@ <objectMethods>
     !@   <object>irate</object>
     !@   <objectMethod>
     !@     <method>readHalos</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ snapshot\argin, \doublezero\ [redshift]\argout, \doubletwo\ [center]\argout, \doubletwo\ [velocity]\argout, \doubleone\ [mass]\argout</arguments>
     !@     <description>Read a snapshot from the \gls{irate} format file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>readSimulation</method>
     !@     <type>\void</type>
     !@     <arguments>\doublezero\ [boxSize]\argout</arguments>
     !@     <description>Read the request properties of the simulation from an \gls{irate} format file.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>copySimulation</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(irate)\textgreater} targetFile\arginout</arguments>
     !@     <description>Copy ``{\normalfont \ttfamily SimulationProperties}'' group from one \gls{irate} file to another.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>copyCosmology</method>
     !@     <type>\void</type>
     !@     <arguments>\textcolor{red}{\textless type(irate)\textgreater} targetFile\arginout</arguments>
     !@     <description>Copy ``{\normalfont \ttfamily Cosmology}'' group from one \gls{irate} file to another.</description>
     !@   </objectMethod>
     !@   <objectMethod>
     !@     <method>writeHalos</method>
     !@     <type>\void</type>
     !@     <arguments>\intzero\ snapshot\argin, \doublezero\ redshift\argin, \doubletwo\ [center]\argin, \doubletwo\ [velocity]\argin, \doubleone\ [mass]\argin</arguments>
     !@     <description></description>
     !@   </objectMethod>
     !@ </objectMethods>
     procedure :: readHalos      => irateReadHalos
     procedure :: readSimulation => irateReadSimulation
     procedure :: copySimulation => irateCopySimulation
     procedure :: copyCosmology  => irateCopyCosmology
     procedure :: writeHalos     => irateWriteHalos
  end type irate

  interface irate
     module procedure irateConstructor
  end interface irate

contains

  function irateConstructor(fileName,cosmologyParameters_,cosmologyFunctions_)
    !% Constructor for \gls{irate} file interface class.
    use :: ISO_Varying_String, only : assignment(=)
    implicit none
    type     (irate                   )                        :: irateConstructor
    character(len=*                   ), intent(in   )         :: fileName
    class    (cosmologyFunctionsClass ), intent(in   ), target :: cosmologyFunctions_
    class    (cosmologyParametersClass), intent(in   ), target :: cosmologyParameters_

    irateConstructor%cosmologyParameters_ => cosmologyParameters_
    irateConstructor%cosmologyFunctions_  => cosmologyFunctions_
    irateConstructor%fileName             =  trim(fileName)
    return
  end function irateConstructor

  subroutine irateReadHalos(self,snapshot,redshift,center,velocity,mass,IDs)
    !% Read requested properties of halos from an \gls{irate} file.
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: Cosmology_Parameters            , only : hubbleUnitsLittleH
    use            :: IO_HDF5                         , only : hdf5Object
    use            :: ISO_Varying_String              , only : char              , trim
    use            :: Numerical_Constants_Astronomical, only : massSolar         , megaParsec
    use            :: Numerical_Constants_Prefixes    , only : hecto             , kilo
    implicit none
    class           (irate     ), intent(inout)                                        :: self
    integer                     , intent(in   )                                        :: snapshot
    double precision            , intent(  out)                             , optional :: redshift
    double precision            , intent(  out), allocatable, dimension(:,:), optional :: center          , velocity
    double precision            , intent(  out), allocatable, dimension(  :), optional :: mass
    integer         (c_size_t  ), intent(  out), allocatable, dimension(  :), optional :: IDs
    double precision                           , allocatable, dimension(  :)           :: unitsInCGS
    type            (hdf5Object)                                                       :: irateFile       , snapshotGroup, &
         &                                                                                halosGroup      , thisDataset
    character       (len=13    )                                                       :: snapshotLabel
    double precision                                                                   :: redshiftInternal, expansionFactor

    ! Read data from file.
    write (snapshotLabel,'(a,i5.5)') 'Snapshot',snapshot
    call irateFile%openFile(char(self%fileName),readOnly=.true.)
    snapshotGroup=irateFile    %openGroup(snapshotLabel)
    halosGroup   =snapshotGroup%openGroup('HaloCatalog')
    call snapshotGroup%readAttribute("Redshift",redshiftInternal,allowPseudoScalar=.true.)
    expansionFactor=self%cosmologyFunctions_%expansionFactorFromRedshift(redshiftInternal)
    if (present(redshift)) redshift=redshiftInternal
    if (present(IDs     )) then
       thisDataset=halosGroup%openDataset("HaloID"  )
       call thisDataset%readDataset(datasetValue=IDs    )
       call thisDataset%close()
    end if
    if (present(mass    )) then
       thisDataset=halosGroup%openDataset("Mass"    )
       call thisDataset%readDataset(datasetValue=mass   )
       call thisDataset%readAttribute('unitscgs',unitsInCGS)
       call thisDataset%close()
       mass    =mass    *(unitsInCGS(1)/kilo /massSolar )*self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**unitsInCGS(2)
    end if
    if (present(center )) then
       thisDataset=halosGroup%openDataset("Center"  )
       call thisDataset%readDataset(datasetValue=center  )
       call thisDataset%readAttribute('unitscgs',unitsInCGS)
       call thisDataset%close()
       center  =center  *(unitsInCGS(1)/hecto/megaParsec)*self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**unitsInCGS(2)
    end if
    if (present(velocity)) then
       thisDataset=halosGroup%openDataset("Velocity")
       call thisDataset%readDataset(datasetValue=velocity)
       call thisDataset%readAttribute('unitscgs',unitsInCGS)
       call thisDataset%close()
       velocity=velocity*(unitsInCGS(1)/hecto/kilo      )*self%cosmologyParameters_%hubbleConstant(hubbleUnitsLittleH)**unitsInCGS(2)
    end if
    call halosGroup   %close()
    call snapshotGroup%close()
    call irateFile    %close()
    return
  end subroutine irateReadHalos

  subroutine irateReadSimulation(self,boxSize)
    !% Read requested properties of the simulation from an \gls{irate} file.
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char
    implicit none
    class           (irate     ), intent(inout)           :: self
    double precision            , intent(  out), optional :: boxSize
    type            (hdf5Object)                          :: irateFile, simulationGroup

    call irateFile%openFile(char(self%fileName),readOnly=.true.)
    simulationGroup=irateFile%openGroup('SimulationProperties')
    if (present(boxSize)) call simulationGroup%readAttribute("boxSize",boxSize,allowPseudoScalar=.true.)
    call simulationGroup%close()
    call irateFile      %close()
    return
  end subroutine irateReadSimulation

  subroutine irateCopySimulation(self,targetFile)
    !% Copy ``{\normalfont \ttfamily SimulationProperties}'' group from one \gls{irate} file to another.
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char
    implicit none
    class(irate     ), intent(inout) :: self
    type (irate     ), intent(inout) :: targetFile
    type (hdf5Object)                :: selfIRATEFile, targetIRATEFile

    call selfIRATEFile  %openFile(char(self      %fileName),readOnly=.true. )
    call targetIRATEFile%openFile(char(targetFile%fileName),readOnly=.false.)
    call selfIRATEFile  %copy    ("SimulationProperties"   ,targetIRATEFile )
    call selfIRATEFile  %close   (                                          )
    call targetIRATEFile%close   (                                          )
    return
  end subroutine irateCopySimulation

  subroutine irateCopyCosmology(self,targetFile)
    !% Copy ``{\normalfont \ttfamily Cosmology}'' group from one \gls{irate} file to another.
    use :: IO_HDF5           , only : hdf5Object
    use :: ISO_Varying_String, only : char
    implicit none
    class(irate     ), intent(inout) :: self
    type (irate     ), intent(inout) :: targetFile
    type (hdf5Object)                :: selfIRATEFile, targetIRATEFile

    call selfIRATEFile  %openFile(char(self      %fileName),readOnly=.true. )
    call targetIRATEFile%openFile(char(targetFile%fileName),readOnly=.false.)
    call selfIRATEFile  %copy    ("Cosmology"              ,targetIRATEFile )
    call selfIRATEFile  %close   (                                          )
    call targetIRATEFile%close   (                                          )
    return
  end subroutine irateCopyCosmology

  subroutine irateWriteHalos(self,snapshot,redshift,center,velocity,mass,IDs,overwrite,objectsOverwritable)
    !% Read requested properties of halos to an \gls{irate} file.
    use, intrinsic :: ISO_C_Binding                   , only : c_size_t
    use            :: IO_HDF5                         , only : hdf5Object
    use            :: ISO_Varying_String              , only : char
    use            :: Numerical_Constants_Astronomical, only : massSolar , megaParsec
    use            :: Numerical_Constants_Prefixes    , only : hecto     , kilo
    implicit none
    class           (irate     ), intent(inout)                           :: self
    integer                     , intent(in   )                           :: snapshot
    double precision            , intent(in   )                           :: redshift
    double precision            , intent(in   ), dimension(:,:), optional :: center       , velocity
    double precision            , intent(in   ), dimension(  :), optional :: mass
    integer         (c_size_t  ), intent(in   ), dimension(  :), optional :: IDs
    logical                     , intent(in   )                , optional :: overwrite    , objectsOverwritable
    type            (hdf5Object)                                          :: irateFile    , snapshotGroup      , &
         &                                                                   halosGroup   , thisDataset
    character       (len=13    )                                          :: snapshotLabel
    !# <optionalArgument name="overwrite"           defaultsTo=".false."/>
    !# <optionalArgument name="objectsOverwritable" defaultsTo=".false."/>
    
    ! Write data to file.
    write (snapshotLabel,'(a,i5.5)') 'Snapshot',snapshot
    call irateFile%openFile(char(self%fileName),readOnly=.false.,overWrite=overWrite_,objectsOverwritable=objectsOverwritable_)
    snapshotGroup=irateFile    %openGroup(snapshotLabel)
    halosGroup   =snapshotGroup%openGroup('HaloCatalog')
    call snapshotGroup%writeAttribute(redshift,"Redshift")
    if (present(IDs     )) then
       call halosGroup%writeDataset   (IDs                              ,'HaloID'  ,'Halo IDs'                                          )
    end if
    if (present(mass    )) then
       call halosGroup%writeDataset   (mass                             ,'Mass'    ,'Halo masses'           ,datasetReturned=thisDataset)
       call thisDataset%writeAttribute('Msolar'                         ,'unitname'                                                     )
       call thisDataset%writeAttribute([kilo*massSolar   , 0.0d0, 0.0d0],'unitscgs'                                                     )
       call thisDataset%close()
    end if
    if (present(center  )) then
       call halosGroup%writeDataset   (center                           ,'Center'  ,'Halo center positions' ,datasetReturned=thisDataset)
       call thisDataset%writeAttribute('Mpc'                            ,'unitname'                                                     )
       call thisDataset%writeAttribute([hecto*megaparsec , 0.0d0,-1.0d0],'unitscgs'                                                     )
       call thisDataset%close()
    end if
    if (present(velocity)) then
       call halosGroup%writeDataset   (velocity                         ,'Velocity','Halo center velocities',datasetReturned=thisDataset)
       call thisDataset%writeAttribute('Mpc'                           ,'unitname'                                                     )
       call thisDataset%writeAttribute([kilo*hecto      , 0.0d0, 0.0d0],'unitscgs'                                                     )
       call thisDataset%close()
    end if
    call halosGroup   %close()
    call snapshotGroup%close()
    call irateFile    %close()
    return
  end subroutine irateWriteHalos

end module IO_IRATE
