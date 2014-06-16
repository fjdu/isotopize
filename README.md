## Usage

Compile the Fortran code with ```make``` in the ```src``` directory, which will
create an executable named ```isotopize```, then run
```
    isotopize config_new.dat
```

The converted file will be produced in the same directory with "isotopized"
attached to the original file name.

## File format conversion

To convert the herbst format to the format understandable by my code, I have to

1. Delete all lines beginning with ```#```, or replace ```#``` with ```!```.  The won't appear in the output file anyway.
1. Shift the columns so that they follow the following format:
    ```!23456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789ABC123456789123456789123456789123456123456123121231212```
1. Replace ```E ``` with ```E-```.  Note the space, otherwise the file will become unaligned.

To convert the isotopized file back to the original format, follow the reverse of the above procedure.

##Format of the config file

```fortran
&PhysicalParameters
  elementOld = 'C' ! Element to be isotopized.
  elementNew = 'X' ! Symbol to be used for the isotope.
  nDeutDegree = 1  ! Degree of isotopization.
  nOtherDeutMax = 20 ! Can be ignored as far as it is an integer not too small.
  noDMaxMetal = 20  ! ...
  noDEleAbundance = 0D0  ! Set to zero.
/
&Paths
  path = "./"  ! Folder containing the input network.
  fReactions = "rate12_umist_reformatted.dat"  ! Network file.
/
