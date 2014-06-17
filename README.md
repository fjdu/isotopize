## Download

Go to any directory you want, then
```
    git clone https://github.com/fjdu/isotopize isotopize
```
and all the files will be downloaded to a subdirectory named ```isotopize```.

## Usage

Compile the Fortran code with ```make``` in the ```src``` directory, which will
create an executable named ```isotopize```, then run
```
    isotopize config.dat
```

See below for the format of the file ```config.dat```.

The converted file will be put in the same directory with "isotopized"
attached to the original file name.

## File format conversion

To convert the herbst format to the format understandable by my code, I have to

1. Delete all the lines beginning with ```#```, or replace ```#``` with ```!```.  The won't appear in the output file anyway.
1. Delete the leading row numbers (with no blank space left at the beginning of each row).
1. Shift the columns so that they follow the format shown below.  If an entry is unapplicable, fill with blank space.
1. Replace ```E``` with ```E-```.  Take care of the spaces and make sure the file is still correctly aligned.
1. Replace ```GRAIN``` with ```Grain```.

All the above can be done within ```vi```, or with a ```python``` script.

To convert the isotopized file back to the original format, follow the reverse of the above procedure.

##Format of the input network file

```fortran
!         R1          R2          R3          P1          P2          P3          P4        A        B        C    T1    T2 iT q cT s
!23456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789123456789123456789123456123456123121231212
```

See also example .dat files in the ```src``` directory.

##Format of the config file

```fortran
&PhysicalParameters
  elementOld = 'C' ! Element to be isotopized.
  elementNew = 'Y' ! Symbol to be used for the isotope.
  nDeutDegree = 1  ! Degree of isotopization.
  nOtherDeutMax = 20 ! Can be ignored as far as it is an integer not too small.
  noDMaxMetal = 20  ! ...
  noDEleAbundance = 0D0  ! Set to zero.
  !
  inputFormat = 'Herbst'
  commentChar = '#'
  inputGrainEleName = '(gr)'  ! Symbol for grain surface species in your input file.
  inputGrainName = 'GRAIN'  ! Symbol for grain itself as a species in your input file.
  outputFormat = 'Herbst'
  outputGrainEleName = '(gr)'
/
&Paths
  path = "./"  ! Folder containing the input network.
  fReactions = "rreacs_herb0308_isotopized.dat"  ! Input network file.
/
```
