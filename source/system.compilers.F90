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

!!{
Contains a module which provides compiler commands.
!!}

module System_Compilers
  !!{
  Provides compiler commands.
  !!}
  implicit none
  private
  public :: compiler, compilerOptions

  !![
  <enumeration>
   <name>language</name>
   <description>Enumerates languages for which compilers are available.</description>
   <visibility>public</visibility>
   <entry label="fortran"  />
   <entry label="c"        />
   <entry label="cplusplus"/>
  </enumeration>
  !!]

contains

  function compiler(language)
    !!{
    Return the name of the compiler to use for a given language.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , var_str, varying_string
    implicit none
    type   (varying_string         )                :: compiler
    type   (enumerationLanguageType), intent(in   ) :: language
    type   (varying_string         )                :: compilerEnvironmentVariable
    integer                                         :: status                     , compilerLength

    select case (language%ID)
    case (languageFortran  %ID)
       compilerEnvironmentVariable=var_str('FCCOMPILER' )
       compiler                   =var_str('gfortran'   )
    case (languageC        %ID)
       compilerEnvironmentVariable=var_str('CCOMPILER'  )
       compiler                   =var_str('gcc'        )
    case (languageCPlusPlus%ID)
       compilerEnvironmentVariable=var_str('CPPCOMPILER')
       compiler                   =var_str('g++'        )
    case default
       call Error_Report('unknown language'//{introspection:location})
    end select
    call Get_Environment_Variable(char(compilerEnvironmentVariable),length=compilerLength,status=status)
    if (status == 0) compiler=compilerRetrieve(char(compilerEnvironmentVariable),compilerLength)
    return
  end function compiler

  function compilerOptions(language)
    !!{
    Return compiler options to use for a given language.
    !!}
    use :: Error             , only : Error_Report
    use :: ISO_Varying_String, only : char        , var_str, varying_string
    implicit none
    type   (varying_string         )                :: compilerOptions
    type   (enumerationLanguageType), intent(in   ) :: language
    type   (varying_string         )                :: compilerEnvironmentVariable
    integer                                         :: status                     , compilerLength

    compilerOptions=var_str('')
    select case (language%ID)
    case (languageFortran  %ID)
       compilerEnvironmentVariable=var_str('GALACTICUS_FCFLAGS' )
    case (languageC        %ID)
       compilerEnvironmentVariable=var_str('GALACTICUS_CFLAGS'  )
    case (languageCPlusPlus%ID)
       compilerEnvironmentVariable=var_str('GALACTICUS_CPPFLAGS')
    case default
       call Error_Report('unknown language'//{introspection:location})
    end select
    call Get_Environment_Variable(char(compilerEnvironmentVariable),length=compilerLength,status=status)
    if (status == 0) compilerOptions=compilerRetrieve(char(compilerEnvironmentVariable),compilerLength)
    return
  end function compilerOptions

  function compilerRetrieve(compilerEnvironmentVariable,compilerLength)
    !!{
    Retrieve the compiler command from an environment variable.
    !!}
    use :: ISO_Varying_String, only : assignment(=), varying_string
    implicit none
    type     (varying_string    )                :: compilerRetrieve
    integer                      , intent(in   ) :: compilerLength
    character(len=*             ), intent(in   ) :: compilerEnvironmentVariable
    character(len=compilerLength)                :: compilerCommand

    call Get_Environment_Variable(compilerEnvironmentVariable,value=compilerCommand)
    compilerRetrieve=compilerCommand
    return
  end function compilerRetrieve

end module System_Compilers
