# crawler
a mini MD engine written in Julia for educational purposes


crawler.jl is a small educational MD engine built from scratch that can energy-minimize and simulate the dynamics
of generic polymers using periodic boundary conditions. It supports bond, angle, dihedral, coulomb, and VDW (LJ) forces.
It also supports bond length constraints with the SHAKE algorithm and Langevin dynamics (beta) with or without constraints
crawler.jl takes in and outputs xyz files for both structures and trajectories and requires a topology file (see below)
The goal is to learn more about MD and learn Julia.

>>> source files <<<
structs.jl: definitions of complex datatypes such as Topology and Energies
IO.jl: read/write-related stuff
integrators.jl: minimization, velocity Verlet, in the future also velocity Verlet with SHAKE
forces.jl: energies and gradients (forces)
ewald.jl: short- and long-range Coulomb energy and forces
constraints.jl: SHAKE algorithm for bond length constraints
util.jl: utility functions used throughout (periodic boundary conditions, etc.)

>>> Sources <<<<
- steepest descent minimization algorithm
	https://en.wikipedia.org/wiki/Gradient_descent, Barzilai–Borwein method for choosing step size
- velocity Verlet algorithm, Langevin dynamics, matrix-SHAKE algorithm
	Mark Tuckerman, Statistical Mechanics: Theory and Molecular Simulation
- angle and dihedral Cartesian derivatives
	https://grigoryanlab.org/docs/dynamics_derivatives.pdf
	https://salilab.org/modeller/9v6/manual/node436.html

Topology format (must be in this order!)

<atomTypes>
<end>

<masses>
1 atomType1
2 ...
<end>

<masses>
1 mass1
2 ...
<end>

<bonds>
atom1Idx atom2Idx length k
<end>

<angles>
atom1 atom2 atom3 angle k
<end>

<dihedrals>
atom1 atom2 atom3 atom4 angle k
<end>

<charges>
1 q1
2 ...
<end>

<vdw>
atom1 atom2 sigma epsilon
<end>

<box>
x y z
<end>

this program uses atomic units everywhere!
	charge: elementary charges
	length: bohr radii (bohr)
	energy: Hartrees
	mass: electron masses
	time: h-bar / E_h
	angles: degrees (in topology), radians (internally)

	bond spring constants are in units of (Hartree / bohr^2), corresponding to N/m in SI
	angle and dihedral spring constants are in units of (Hartree / (bohr.rad))

The force field expression is simple and based on harmonic potentials for most terms:
V(x) = bond + angle + dihedral + coulomb + lennard-jones (VDW), or
V(x) = Σ(1/2)*k*x^2 + Σ(1/2)*kAngle*(θ - θ0)^2 + Σ(1/2)*kDihedral*(α-α0)^2 + ΣΣ(qi*qj)/rij + Σ4ϵ((σ / r)^12 - (σ / r)^6)
