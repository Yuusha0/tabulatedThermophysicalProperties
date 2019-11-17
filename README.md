# tabulatedThermophysicalProperties
Tabulated thermophysical properties for OpenFOAM 7

## REQUIREMENTS

### OpenFOAM 7
* https://openfoam.org/

### Python 3
* https://www.python.org/downloads/
* Tested with version 3.5.3 and higher.
  May work with other python 3 versions.
  
## INSTALLATION

* Set the location of the directory in /etc/bashrc file (default $WM_PROJECT_USER_DIR).
* Source this etc/bashrc file in your user bashrc.  
* Compile all sources with ./Allwmake in the src/ directory.


## USAGE

### python script

* Pyhton script can import CSV like tables.
* It can transpose OpenFOAM style 2D tables.
* Use python3 --help to view syntax and options.

### Tabulated thermophysical model

* Compile your solver with:
  * libTabularThermophysicalModels
  * libuserspecie
  * libtabularReactionThermophysicalModels (for multi-species simulation)

* Or use dynamic linking by adding these libraries to your controlDict

* In your constant/thermophysicalProperties file:
  * Change transport, thermo, and equationOfState to the desired values.
  * Use hePsiThermo, heRhoThermo or heTabularThermo for type value.
  * Create subdictionaries for each property using tabulated data.
    * Add fileName value.
    * Add outOfBounds value:
      * ERROR which exits with a Fatal Error
      * CLAMP clamp the values (default)
      * WARN which issues warning and extrapolates value
      * EXTRAPOLATE which extrapolates value
    * Be careful when using WARN and log file, log may be very huge.
    * Add searchMethod value:
      * simple which looking all over the table
      *	uniform which is usefull for reguar tables only
      * bisect which uses a bijection algorithm to find the right value
    * Be careful, uniform does not work on non-uniform table. It is at least 2.5 times faster thant simple for a 238x1 table.
    * Bisection is 1.95 times faster than simple for a 238x1 table.

* If your thermophysical model is not defined:
  * Add it in tabularThermos.C (for compressibility-based solvers) or rhoTabularThermos.C (for density-based solvers)
  * Add it in tabularReactionThermos.C and makeTabularChemistryReaders.C (for multi-species)

## COMPATIBILITY

* OpenFOAM 5.x and OpenFOAM 6 users may have compatibility issues.
* For OpenFOAM 5.x use v2.0.3.
* For OpenFOAM 6 use v4.x
* v2.x is no longer supported.
* Only critical bug fixes are backported to v4.x . So change your OpenFOAM version to 7 if you can.

## MISCELLANEOUS

* Multi-species works but it is very slow for large tables due to OpenFOAM code design.
* Only 1 non reactive model of multi-species is currently implemented.

## LICENSE AND COPYRIGHT

See [license](LICENSE) file.
