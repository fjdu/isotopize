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

The converted file will be put in the same directory as the input network with
the word "isotopized" attached to the original file name.

## File format conversion

_This is not needed anymore because the code now has the option to read directly the "Herbst" format (so to call it)._

To convert the Herbst format to the format understandable by my code, I have to

1. Delete all the lines beginning with ```#```, or replace ```#``` with ```!```.  The won't appear in the output file anyway.
1. Delete the leading row numbers (with no blank space left at the beginning of each row).
1. Shift the columns so that they follow the format shown below.  If an entry is unapplicable, fill with blank space.
1. Replace ```E``` with ```E-```.  Take care of the spaces and make sure the file is still correctly aligned.
1. Replace ```GRAIN``` with ```Grain```.

All the above can be done within ```vi```, or with a ```python``` script.

To convert the isotopized file back to the original format, follow the reverse of the above procedure.

##Format of the input network file

```fortran
# The "Herbst" format
#   N R1         | R2         | P1         | P2         | P3         | P4         |   A      | B       | C       |   T
#23456123456789ABCD123456789ABCD123456789ABCD123456789ABCD123456789ABCD123456789ABCDE123456789A123456789A1234567891234
 6328 COOCH4+      GRAIN-       CH3          CO2          H            GRAIN0         3.14E-10  6.00E+01  5.00E-01  23
```

```fortran
! The "nonHerbst" format
!R1        |R2         |R3         |P1         |P2         |P3         |P4         |A       |B       |C       |T1   |T2   | iT q cT s
!23456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789123456789123456789123456123456123121231212
C5H4N+      Grain-                  C5N         H2          H2          Grain0       3.14e-10 7.80e+01 5.00e-01             23
```

See also example .dat files in the ```src``` directory.

##Format of the config file

```fortran
&PhysicalParameters
  elementOld = 'O' ! Element to be isotopized.
  elementNew = 'Z' ! Symbol to be used for the isotope.
  nDeutDegree = 20 ! Degree of isotopization.
  nOtherDeutMax = 99 ! If the number of metal atoms in a reaction is larger than this number, it will only be isotopized once.
  noDMaxMetal = 99  ! If the number of metal atoms in a reaction is larger than this number, it will not be isotopized.
  noDEleAbundance = 0D0  ! Set to zero.
  !
  inputFormat = 'Herbst'
  commentChar = '#'
  inputGrainEleName = '(gr)'
  inputGrainName = 'GRAIN'
  outputFormat = 'Herbst'
  outputGrainEleName = '(gr)'
  !
  grain_special = .true.
  copy_rates = .true. ! Whether to simply copy the rate coefficients to the newly generated reactions
/
&Paths
  path = "./"  ! Folder containing the input network.
  fReactions = "rreacs_herb0308_test.dat"  ! Network file.
/
```
