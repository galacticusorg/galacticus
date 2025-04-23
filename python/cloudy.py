# Provides functions for working with Cloudy
# Andrew Benson (22-April-2025)
import os
import re
import sys
import subprocess
import urllib.request

def initialize(options):
    # Initialize Cloudy by downloading and compiling.
    # Determine Cloudy version.
    if os.environ.get('GALACTICUS_EXEC_PATH') is not None and not 'version' in options:
        dependencies = open(os.environ.get('GALACTICUS_EXEC_PATH')+"/aux/dependencies.yml")
        for line in dependencies:
            match = re.match(r'^cloudy:\s+([\d\.]+)',line)
            if match:
                options['version'] = match.group(1)
        dependencies.close()
    cloudyVersion      = options['version'] if 'version' in options else "23.01"
    cloudyVersionMajor = cloudyVersion[0:cloudyVersion.index(".")]
    # Specify Cloudy path.
    cloudyPath = os.environ.get('GALACTICUS_DATA_PATH')+"/dynamic/c"+cloudyVersion
    os.makedirs(cloudyPath, exist_ok=True)
    # Download the code.
    if not os.path.isfile(cloudyPath+".tar.gz"):
        print("downloading Cloudy code...")
        try:
            urllib.request.urlretrieve("http://data.nublado.org/cloudy_releases/c"+cloudyVersionMajor+"/c"+cloudyVersion+".tar.gz", cloudyPath+".tar.gz")
        except HTTPError as e:
            print(f"FATAL: failed to download Cloudy: HTTP Error: {e.code} - {e.reason}")
            sys.exit()
        except Exception as e:
            print(f"FATAL: failed to download Cloudy: {e}") 
            sys.exit()
    # Unpack the code.
    if not os.path.isfile(cloudyPath+"/source/Makefile"):
        print("unpacking Cloudy code...")
        unpack = subprocess.run(['tar', '-x', '-v', '-z', '-C', os.environ.get('GALACTICUS_DATA_PATH')+'/dynamic', '-f', cloudyPath+'.tar.gz'])
        if unpack.returncode != 0:
            print("FATAL: failed to unpack Cloudy")
            sys.exit()
    if not os.path.isfile(cloudyPath+"/source/Makefile"):
        print("FATAL: failed to unpack Cloudy code.")
        sys.exit()
    # Build the code.
    if not os.path.isfile(cloudyPath+"/source/cloudy.exe"):
        print("compiling Cloudy code...")
        buildCommand      = "cd "+cloudyPath+"/source; chmod u=wrx configure.sh capabilities.pl; "
        buildCommand     += "" if os.environ.get('CLOUDY_COMPILER_PATH') is None else "export PATH="+os.environ.get('CLOUDY_COMPILER_PATH')+":$PATH; "
        buildCommand     += " make -f Makefile_modified";
        extra             = "EXTRA = -static\n" if os.environ.get('CLOUDY_STATIC_BUILD') == "yes" else "EXTRA = \n"
        makeFile          = open(cloudyPath+"/source/Makefile"         ,'r')
        makeFileModified  = open(cloudyPath+"/source/Makefile_modified",'w')
        for line in makeFile:
            lineModified = re.sub(r'^EXTRA\s*=.*', extra, line)
            makeFileModified.write(lineModified)
        makeFile        .close()
        makeFileModified.close()
        build = subprocess.run(buildCommand, shell=True)
        if build.returncode != 0:
            print("FATAL: failed to build Cloudy code")
            sys.exit()
    # Return path and version.
    return (cloudyPath+"/", cloudyVersion)
