# tabulatedThermophysicalProperties
Tabulated thermophysical properties for OpenFOAM 5.x.

## REQUIREMENTS

### OpenFOAM 5.x
* https://openfoam.org/

### Python 3
* https://www.python.org/downloads/
* Tested with version 3.5.3 and higher.
  May work with other python 3 versions.
  
## INSTALLATION

* Copy the source files in your user OpenFOAM tree (often /home/user/OpenFOAM/user-5.x/)
* Compile all sources:
  * Use `wmake libso` in each directory where a Make folder is present
  * Use `wmakeLnInclude -update .` in src/OpenFOAM directory

## USAGE

### python script

* Only table transposition works.
* Importation from CSV file is planned but not implemented.
* To transpose a table:
  * Copy the table in your python script repository.
  * Open thermophysicalTable.py file.
  * Change the values in thermo.read("your_file") and thermo.write("your_file").
  * Use `python3 thermophysicalTable.py`

### Tabulated thermophysical model

* Compile your solver with:
  * libTabularThermophysicalModels
  * libuserspecie
  * libtabularReactionThermophysicalModels (for multi-species simulation)

* In your constant/thermophysicalProperties file:
  * Change transport, thermo, and equationOfState to the desired values.
  * Use hePsiThermo or heRhoThermo for type value.
    heTabularThermo currently doesn't work.
  * Create subdictionaries for each property using tabulated data.
    * Add fileName value.
    * Add outOfBounds value:
      * ERROR which exits with a Fatal Error
      * WARN which issues warning and extrapolates value (default)
      * EXTRAPOLATE which extrapolates value
    Be careful when using WARN and log file, log may be very huge.

* If your thermophysical model is not defined:
  * Add it in tabularThermos.C (for single specie)
  * Add it in tabularReactionThermos.C and makeTabularChemistryReaders.C (for multi-species)

## MISCELLANEOUS

* heTabularThermo currently doesn't work.
* Multi-species works but it is very slow for large tables due to OpenFOAM code design.
* Only 1 non reactive model of multi-species is currently implemented.

## LICENSE AND COPYRIGHT

See [license](LICENSE) file.
