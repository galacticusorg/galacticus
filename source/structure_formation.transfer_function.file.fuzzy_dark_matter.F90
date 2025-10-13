!! Copyright 2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018,
!!           2019, 2020, 2021, 2022, 2023, 2024, 2025
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

!+    Contributions to this file made by: Andrew Benson, Xiaolong Du.

  !!{
  Implements a file-based transfer function class for fuzzy dark matter.
  !!}

  use :: Dark_Matter_Particles, only : darkMatterParticleClass, darkMatterParticleFuzzyDarkMatter

  !![
  <transferFunction name="transferFunctionFileFuzzyDarkMatter">
   <description>
  Provides a fuzzy dark matter transfer function from a tabulation given in an HDF5 file with the following structure:
  \begin{verbatim}
  HDF5 "transferFunction.hdf5" {
  GROUP "/" {
     ATTRIBUTE "description" {
        DATATYPE  H5T_STRING {
           STRSIZE 71;
           STRPAD H5T_STR_NULLTERM;
           CSET H5T_CSET_ASCII;
           CTYPE H5T_C_S1;
        }
        DATASPACE  SCALAR
     }
     ATTRIBUTE "fileFormat" {
        DATATYPE  H5T_STD_I32LE
        DATASPACE  SCALAR
     }
     ATTRIBUTE "redshift" {
        DATATYPE  H5T_STD_I32LE
        DATASPACE  SCALAR
     }
     GROUP "extrapolation" {
        GROUP "wavenumber" {
           ATTRIBUTE "high" {
              DATATYPE  H5T_STRING {
                 STRSIZE 11;
                 STRPAD H5T_STR_NULLTERM;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
              DATASPACE  SCALAR
           }
           ATTRIBUTE "low" {
              DATATYPE  H5T_STRING {
                 STRSIZE 3;
                 STRPAD H5T_STR_NULLTERM;
                 CSET H5T_CSET_ASCII;
                 CTYPE H5T_C_S1;
              }
              DATASPACE  SCALAR
           }
        }
     }
     GROUP "parameters" {
        ATTRIBUTE "HubbleConstant" {
           DATATYPE  H5T_STRING {
              STRSIZE 4;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "OmegaBaryon" {
           DATATYPE  H5T_STRING {
              STRSIZE 6;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "OmegaDarkEnergy" {
           DATATYPE  H5T_STRING {
              STRSIZE 5;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "OmegaMatter" {
           DATATYPE  H5T_STRING {
              STRSIZE 5;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "fuzzyDMMass" {
           DATATYPE  H5T_STRING {
              STRSIZE 6;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
        ATTRIBUTE "fuzzyDMDensityFraction" {
           DATATYPE  H5T_STRING {
              STRSIZE 6;
              STRPAD H5T_STR_NULLTERM;
              CSET H5T_CSET_ASCII;
              CTYPE H5T_C_S1;
           }
           DATASPACE  SCALAR
        }
     }
     DATASET "transferFunction" {
        DATATYPE  H5T_IEEE_F64LE
        DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
     }
     DATASET "wavenumber" {
        DATATYPE  H5T_IEEE_F64LE
        DATASPACE  SIMPLE { ( 1000 ) / ( 1000 ) }
     }
  }
  }
  \end{verbatim}
   </description>
  </transferFunction>
  !!]
  type, extends(transferFunctionFile) :: transferFunctionFileFuzzyDarkMatter
     !!{
     A transfer function class which interpolates a fuzzy dark matter transfer function given in a file.
     !!}
     private
     class(darkMatterParticleClass), pointer :: darkMatterParticle_ => null()
   contains
     procedure :: readFile        => fileFuzzyDarkMatterReadFile
     procedure :: halfModeMass    => fileFuzzyDarkMatterHalfModeMass
     procedure :: quarterModeMass => fileFuzzyDarkMatterQuarterModeMass
  end type transferFunctionFileFuzzyDarkMatter

  interface transferFunctionFileFuzzyDarkMatter
     !!{
     Constructors for the fileFuzzyDarkMatter transfer function class.
     !!}
     module procedure fileFuzzyDarkMatterConstructorParameters
     module procedure fileFuzzyDarkMatterConstructorInternal
  end interface transferFunctionFileFuzzyDarkMatter

contains

  function fileFuzzyDarkMatterConstructorParameters(parameters) result(self)
    !!{
    Constructor for the fileFuzzyDarkMatter transfer function class which takes a parameter set as input.
    !!}
    use :: Input_Parameters, only : inputParameter, inputParameters
    implicit none
    type            (transferFunctionFileFuzzyDarkMatter)                :: self
    type            (inputParameters                    ), intent(inout) :: parameters
    class           (cosmologyParametersClass           ), pointer       :: cosmologyParameters_
    class           (cosmologyFunctionsClass            ), pointer       :: cosmologyFunctions_
    class           (darkMatterParticleClass            ), pointer       :: darkMatterParticle_
    type            (varying_string                     )                :: fileName
    double precision                                                     :: redshift

    !![
    <inputParameter>
      <name>fileName</name>
      <source>parameters</source>
      <description>The name of the file from which to read a tabulated transfer function.</description>
    </inputParameter>
    <inputParameter>
      <name>redshift</name>
      <source>parameters</source>
      <defaultValue>0.0d0</defaultValue>
      <description>The redshift of the transfer function to read.</description>
    </inputParameter>
    <objectBuilder class="cosmologyParameters" name="cosmologyParameters_" source="parameters"/>
    <objectBuilder class="cosmologyFunctions"  name="cosmologyFunctions_"  source="parameters"/>
    <objectBuilder class="darkMatterParticle"  name="darkMatterParticle_"  source="parameters"/>
    !!]
    self=transferFunctionFileFuzzyDarkMatter(char(fileName),redshift,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_)
    !![
    <inputParametersValidate source="parameters"/>
    <objectDestructor name="cosmologyParameters_"/>
    <objectDestructor name="cosmologyFunctions_" />
    <objectDestructor name="darkMatterParticle_" />
    !!]
    return
  end function fileFuzzyDarkMatterConstructorParameters

  function fileFuzzyDarkMatterConstructorInternal(fileName,redshift,cosmologyParameters_,cosmologyFunctions_,darkMatterParticle_) result(self)
    !!{
    Internal constructor for the fileFuzzyDarkMatter transfer function class.
    !!}
    implicit none
    type            (transferFunctionFileFuzzyDarkMatter)                        :: self
    character       (len=*                              ), intent(in   )         :: fileName
    double precision                                     , intent(in   )         :: redshift
    class           (cosmologyParametersClass           ), intent(in   ), target :: cosmologyParameters_
    class           (cosmologyFunctionsClass            ), intent(in   ), target :: cosmologyFunctions_
    class           (darkMatterParticleClass            ), intent(in   ), target :: darkMatterParticle_
    !![
    <constructorAssign variables="fileName, redshift, *cosmologyParameters_, *cosmologyFunctions_, *darkMatterParticle_"/>
    !!]

    self%time=self%cosmologyFunctions_%cosmicTime(self%cosmologyFunctions_%expansionFactorFromRedshift(redshift))
    call self%readFile(fileName)
    return
  end function fileFuzzyDarkMatterConstructorInternal

  subroutine fileFuzzyDarkMatterReadFile(self,fileName)
    !!{
    Read in the transfer function data from a file.
    !!}
    use :: Display             , only : displayMessage
    use :: Error               , only : Error_Report
    use :: HDF5_Access         , only : hdf5Access
    use :: IO_HDF5             , only : hdf5Object
    use :: Numerical_Comparison, only : Values_Differ
    implicit none
    class           (transferFunctionFileFuzzyDarkMatter), intent(inout) :: self
    character       (len=*                              ), intent(in   ) :: fileName
    double precision                                                     :: fuzzyDMMass, fuzzyDMDensityFraction
    type            (hdf5Object                         )                :: fileObject , parametersObject

    ! Open and read the HDF5 data file.
    !$ call hdf5Access%set()
    fileObject=hdf5Object(fileName,readOnly=.true.)
    ! Check that the fuzzy dark matter parameters match.
    parametersObject=fileObject%openGroup('parameters')
    call parametersObject%readAttribute('fuzzyDMMass'           ,fuzzyDMMass           )
    call parametersObject%readAttribute('fuzzyDMDensityFraction',fuzzyDMDensityFraction)
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       if (Values_Differ(fuzzyDMMass           ,darkMatterParticle_%mass           (),relTol=1.0d-3)) &
            & call displayMessage('fuzzyDMMass from transfer function file does not match internal value'           )
       if (Values_Differ(fuzzyDMDensityFraction,darkMatterParticle_%densityFraction(),absTol=1.0d-3)) &
            & call displayMessage('fuzzyDMDensityFraction from transfer function file does not match internal value')
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    !$ call hdf5Access%unset()
    ! Read the data file.
    call self%transferFunctionFile%readFile(fileName)
    return
  end subroutine fileFuzzyDarkMatterReadFile

  double precision function fileFuzzyDarkMatterHalfModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of two relative to a \gls{cdm} transfer function. Here the
    fitting function from \cite{hu_fuzzy_2000} has been used. Note that:
    (1) In \cite{hu_fuzzy_2000}, the half-mode wavenumber is defined as the scale at
        which the matter power spectrum, instead of the transfer function, is suppressed
        by a factor of two. A correction factor has been added in the calculations to be
        consistent with the common definition.
    (2) For a mixed \gls{cdm} and \gls{fdm} model, the half-mode wavenumber
        is not well defined.
    !!}
    use :: Error                       , only : errorStatusSuccess
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (transferFunctionFileFuzzyDarkMatter), intent(inout), target   :: self
    integer                                              , intent(  out), optional :: status
    double precision                                                               :: matterDensity, wavenumberHalfMode
    double precision                                                               :: m22

    fileFuzzyDarkMatterHalfModeMass=0.0d0
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       if (darkMatterParticle_%densityFraction() == 1.0d0) then
          matterDensity=+self%cosmologyParameters_%OmegaMatter    () &
               &        *self%cosmologyParameters_%densityCritical()
          ! Particle mass in units of $10^{-22}$~eV.
          m22                            =+darkMatterParticle_%mass() &
               &                          *kilo                       &
               &                          /1.0d-22
          wavenumberHalfMode             =+1.108d0              &
               &                          *4.5d0                &
               &                          *m22**(4.0d0/9.0d0)
          ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
          ! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].
          fileFuzzyDarkMatterHalfModeMass=+4.0d0                &
               &                          *Pi                   &
               &                          /3.0d0                &
               &                          *matterDensity        &
               &                          *(                    &
               &                            +Pi                 &
               &                            /wavenumberHalfMode &
               &                           )**3
       else
          call Error_Report('half-mode mass is not well defined for a mixed CDM and fuzzy dark matter model'//{introspection:location})
       end if
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    if (present(status)) status=errorStatusSuccess
    return
  end function fileFuzzyDarkMatterHalfModeMass

  double precision function fileFuzzyDarkMatterQuarterModeMass(self,status)
    !!{
    Compute the mass corresponding to the wavenumber at which the transfer function is
    suppressed by a factor of four relative to a \gls{cdm} transfer function.
    !!}
    use :: Error                       , only : errorStatusSuccess
    use :: Numerical_Constants_Math    , only : Pi
    use :: Numerical_Constants_Prefixes, only : kilo
    implicit none
    class           (transferFunctionFileFuzzyDarkMatter), intent(inout), target   :: self
    integer                                              , intent(  out), optional :: status
    double precision                                                               :: matterDensity, wavenumberQuarterMode
    double precision                                                               :: m22

    fileFuzzyDarkMatterQuarterModeMass=0.0d0
    select type (darkMatterParticle_ => self%darkMatterParticle_)
    class is (darkMatterParticleFuzzyDarkMatter)
       if (darkMatterParticle_%densityFraction() == 1.0d0) then
          matterDensity=+self%cosmologyParameters_%OmegaMatter    () &
               &        *self%cosmologyParameters_%densityCritical()
          ! Particle mass in units of $10^{-22}$~eV.
          m22                               =+darkMatterParticle_%mass() &
               &                             *kilo                       &
               &                             /1.0d-22
          wavenumberQuarterMode             =+1.230d0                 &
               &                             *4.5d0                   &
               &                             *m22**(4.0d0/9.0d0)
          ! Compute corresponding mass scale. As a default choice, the wavenumber is converted to a length scale assuming
          ! R = λ/2 = π/k [see Eq.(9) of Schneider et al. (2012; http://adsabs.harvard.edu/abs/2012MNRAS.424..684S)].
          fileFuzzyDarkMatterQuarterModeMass=+4.0d0                   &
               &                             *Pi                      &
               &                             /3.0d0                   &
               &                             *matterDensity           &
               &                             *(                       &
               &                               +Pi                    &
               &                               /wavenumberQuarterMode &
               &                              )**3
       else
          call Error_Report('quarter-mode mass is not well defined for a mixed CDM and fuzzy dark matter model'//{introspection:location})
       end if
    class default
       call Error_Report('transfer function expects a fuzzy dark matter particle'//{introspection:location})
    end select
    if (present(status)) status=errorStatusSuccess
    return
  end function fileFuzzyDarkMatterQuarterModeMass
