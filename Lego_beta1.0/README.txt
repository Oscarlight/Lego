LEGO Beta v1.0
Copyright (c) 2015-2016 Mingda Oscar Li (ml888@cornell.edu)

1. Introduction
	LEGO is the open source electical device simulator. In version Beta v1.0, it include 
		1) Electrostatics: 1D/2D Poisson for 3D and 2D semiconductor. 
		2) Tranport: 2D to 2D interlayer quantum tunneling current simulation.

2. Dependency
	It depends on several external libaries:
		1) Armadillo C++ Linear Algebra Libary 6.100.0
		2) cubature-1.0.2
		3) OpenBLAS-0.2.14 (optional)
		4) SuperLU_4.3

3. Class dependency:
	MathFunct -> Params -> Material -> Device1D -> Device2D

		   	    -> Poisson1D -> Poisson2D

		            -> Current  -> Tunneling


