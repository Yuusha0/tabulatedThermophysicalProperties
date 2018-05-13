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
  
## LICENSE AND COPYRIGHT

See [license](LICENSE) file.
