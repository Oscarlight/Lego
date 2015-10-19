# Leno

Leno Beta v1.0
by Mingda Oscar Li (ml888@cornell.edu)

1. Introduction

	LEGO is the open source electical device simulator.

	The current version Beta v1.0 is tailed for Thin-TFET exclusively, it include 

		1. Electrostatics: 1D Poisson for 3D and 2D semiconductor. 
		2. Tranport: 2D to 2D interlayer quantum tunneling current simulation.

2. Dependency

	If you would like to recompile it, it depends on several external libaries:

		1. Blas
		2. Lapack
		3. Python2.7

	Luckly, most linux OS has these three pre-installed.

	More compile and link info, see Makefile.

3. How to use it

	*** Current Release is only for Linus ***

        1. Open terminal, go the same folder with Leno_beta1.0
	2. Modify the device file and/or material file as you wish
	3. Run Leno_beta1.0 -d yourDeviceFileName -m yourMaterialFileName
	4. ( Good ) ? Stay_cool : Email ml888@cornell.edu for help

4. Class Tree:

	MathFunct
		*Params
			**Material
				***Device1D
					(Device2D)
			Poisson1D
				(Poisson2D)
			Transport
				Tunneling
			ExtractData
			InOut

5. Reference:
	>"Single particle transport in two-dimensional heterojunction interlayer tunneling field effect transistor" http://dx.doi.org/10.1063/1.4866076
	>"Two-Dimensional Heterojunction Interlayer Tunneling Field Effect Transistors (Thin-TFETs)" DOI:10.1109/JEDS.2015.23906
