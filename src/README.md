## Usage

Compile the Fortran code with make, which will create an executable named isotopize, then run
```
    isotopize config_new.dat
```

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
