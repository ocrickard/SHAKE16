# SHAKE16

This is a modernized version of SHAKE91 to run with modern fortran compilers. Here's what we changed:

1. Removed the archaic FORTRAN capitalization in most of the files.
2. Converted to modernized `do` loop syntax for core subroutines and improved indentation.
3. Fixed several issues that prevented compilation with modern gfortran.
4. Added a makefile.

## Building

To build SHAKE on a Unix-based system, first install gfortran. If you are on a mac, you can use homebrew:

```
brew install gfortran
```

Then clone this repo to your machine, and run:

```
make install
```

This will compile the source files using gfortran, and copy linked binary to `/usr/local/bin/SHAKE16`.

## Usage

To invoke SHAKE16, ensure that the binary is in your PATH, then simply invoke `SHAKE16` from the command line.

```
$ SHAKE16
```

You will then be asked for the name of the input file. You may provide the relative path to the CWD to the input file.

```
Name of Input File =
INP.DAT
```

Next, the program will ask for the output file names.
```
Name of Output File #1 (input, peak values .. etc) =
> output1.txt
Name of Output File #2 (time histories .. etc) =
> output2.txt
```

Hit enter after finishing inputting each file name. If all is successful the program will perform the computations, and write the output data to the files you specified.

[Example of input file](https://github.com/ocrickard/SHAKE16/blob/master/Input/INP.DAT)

[Example of output file 1](https://github.com/ocrickard/SHAKE16/blob/master/Input/output1.txt)

[Example of output file 2](https://github.com/ocrickard/SHAKE16/blob/master/Input/output2.txt)

## Original Readme

Program: SHAKE-91

Title: Equivalent Linear Seismic Response Analysis of Horizontally 
Layered Soil Deposits

Developer: P. B. Schnabel, J. Lysmer, and H. B. Seed, Department of Civil
Engineering, University of California, Berkeley 1972.

Modified: I. M. Idriss and J. I. Sun, Department of Civil & Environmental
Engineering, University of California, Davis 1992.

Category: Geotechnical

Platform: PC DOS 6, MS PowerStation Fortran, v. 1.0

Reference: Idriss, I.M., and J.I. Sun, "User's Manual for SHAKE91," 
Department of Civil & Environmental Engineering, University of California,
Davis, California, November 1992.

Schnabel, P.B., J. Lysmer, and H.B. Seed, "SHAKE - A Computer Program for
Earthquake Response Analysis of Horizontally Layered Sites," Earthquake 
Engineering Research Center, Report No. UCB/EERC-72/12. University of 
California, Berkeley, December 1972.

Summary: The SHAKE program has been by far the most widely used program 
for computing the seismic response of horizontally layered soil deposits.
The program computes the response of a semi-infinite horizontally layered
soil deposit overlying a uniform half-space subjected to vertically 
propagating shear waves. The analysis is done in the frequency domain, 
and, therefore, for any set of properties, it is a linear analysis. An 
iterative procedure is used to account for the nonlinear behavior of the 
soils. The object motion (i.e., the motion that is considered to be known)
can be specified at the top of any sublayer within the soil profile or at
the corresponding outcrop. 

The main modifications incorporated in SHAKE91 include the following: 
The number of sublayers was increased from 20 to 50; this should permit a
more accurate representation of deeper and/or softer soil deposits. All 
built-in modulus reduction and damping relationships were removed. These 
relationships are now specified by the user. The maximum shear velocity 
or the maximum modulus are now specified for each sublayer; again these 
are part of the input and therefore the program no longer calculates 
modulus values as a function of either confining pressure or shear 
strength. Object motion is now read from a separate file. Other clean-up 
includes: renumbering of options, elimination of infrequently used options,
user specified periods for calculating spectral ordinates.

[SHAKE91 Manual](https://github.com/ocrickard/SHAKE16/raw/master/SHAKE91%20User%20Manual.pdf)
