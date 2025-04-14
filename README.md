SMRAI 2.1 in Fortran

An emulator of REPPU MHD model based on SMRAI2 (Kataoka et al., 2024)


[Instruction]

1. A random number generator is required. Download the Mersenne Twister
in Fortran (mtfort90.f) via 

https://www.math.sci.hiroshima-u.ac.jp/m-mat/MT/VERSIONS/FORTRAN/fortran.html

and rename it as "mtfort90.f90". Remove the main program from this Fortran 
code to extract the Mersenne Twister module. 


2. A data file of 5-minute solar wind data is also necessary.
It should contain year, day of year, hour, minute, solar-wind dynamic
pressure, density, speed, temperature, and interplanetary magnetic field
Bx, By, and Bz in this order. The file "swdata20170327.dat" is a sample 
solar wind data. 


3. It is necessary to generate the poles for a radial basis function network
in advance as follows. 

> python3 genpoles_spiral.py


4. Build the emulator code as follows. 

> make reppu_emulator


5. Run the emulator.

> reppu_emulator [SOLAR WIND DATA FILE]
