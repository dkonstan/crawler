include("integrators.jl")
include("IO.jl")
using ArgParse
using LinearAlgebra
using SpecialFunctions

#=
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
=#

function main()

	argTable = ArgParseSettings()

	@add_arg_table argTable begin
	    "--topology", "-t"
	        help = "system topology file"
	        arg_type = String
	        required = true
	    "--coordinates", "-c"
	        help = "xyz file with initial coordinates"
	        arg_type = String
	        required = true
	    "--integrator", "-i"
	    	help = "min or vv"
			arg_type = String
			required = true
		"--deltat", "-d"
			help = "timestep size"
			arg_type = Float64
			default = nothing
	    "--steps", "-s"
	        default= 0
	        help = "number of steps"
	        arg_type = Int64
	    "--trajectory", "-x"
	    	help = "trajectory (xyz file)"
	    	default = "traj.xyz"
	    	arg_type = String
	    "--logfile", "-l"
	    	help = "text file to log progress"
	    	default = "log.csv"
	    	arg_type = String
	    "--outfile", "-o"
	        default = "min.xyz"
	    	help = "outfile for minimizations"
	    	arg_type = String
	    "--constrainbonds", "--b"
	    	help = "turn on SHAKE algorithm to constrain bond lengths"
	    	action = :store_true
	    "--loginterval", "-r"
	    	help = "how often to print to logfile (number of steps), default 1"
	    	default = 1
	    	arg_type = Int64
	    "--gamma", "-g"
	    	default = 0.0
	    	arg_type = Float64
	    	help = "collision frequency (1 / time) for Langevin dynamics"
	    "--temperature", "-T"
	    	default = 0.0
	    	arg_type = Float64
	    	help = "desired (starting) temperature"
      	end

    args = parse_args(argTable)

    topologyfile = args["topology"]
    coords = args["coordinates"]
    traj = args["trajectory"]
    dt = args["deltat"]
    integrator = args["integrator"]
    loginterval = args["loginterval"]
    logfile = args["logfile"]
    constrainbonds = args["constrainbonds"]
    steps = args["steps"]
    outfile = args["outfile"]
    temp = args["temperature"]
    collisionFreq = args["gamma"]
	kb = 3.167e-6  # atomic units

	top = readTopology(topologyfile)
	xyz = readCoordinates(coords)

	if args["integrator"] == "min"
		println("finding local energy minimum...")
		minimizeEnergy!(top, xyz, steps, outfile)
	elseif args["integrator"] == "vv"
		println("setting velocities to temperature...")
		vel = zeros((top.nAtoms, 3))
		# https://physics.stackexchange.com/questions/159674/velocity-maxwell-boltzmann-distribution-for-dummies
		for i in 1:size(vel, 1)
			mi = top.masses[i]
			# sample from Maxwell-Boltzmann distribution
			vel[i, 1] = sqrt((kb * temp) / mi) * randn()
			vel[i, 2] = sqrt((kb * temp) / mi) * randn()
			vel[i, 3] = sqrt((kb * temp) / mi) * randn()
		end


		if dt == nothing
			ks = top.bondKs
			ms = top.masses
			dt = (1 / (sqrt(maximum(ks) / minimum(ms)))) * 0.1
			println("setting deltat to $(dt) based on bond spring constants and masses...")
		end

		initiateReport(traj, logfile)
		if constrainbonds == true
			println("running NVE dynamics with constrained bonds...")
			# for now, supports only a single bond length for all bonds, trivial to extend later
			bondLength = norm(xyz[top.bondIdx[1, :], :] .- xyz[top.bondIdx[2, :], :])
			nBonds = size(top.bondIdx, 1)
			λs = zeros(nBonds) # initial guess for Lagrange multipliers
			for frameNumber=1:steps
				if frameNumber % loginterval == 0
					println("step $(frameNumber) of $(steps)")
					report = true
				else
					report = false
				end
				λs = integrateVelocityVerletSHAKE!(top, xyz, vel, dt, frameNumber, report, traj, logfile, bondLength, λs)
			end
		else
			println("running NVE dynamics...")
			for frameNumber=1:steps
				if frameNumber % loginterval == 0
					println("step $(frameNumber) of $(steps)")
					report = true
				else
					report = false
				end
				integrateVelocityVerlet!(top, xyz, vel, dt, frameNumber, report, traj, logfile)
			end
		end

	elseif args["integrator"] == "langevin"

		if args["gamma"] == 0.0
			println("choose a non-zero γ (collision frequency) for Langevin dynamics!")
			exit()
		end

		if args["steps"] == 0
			println("set number of steps [--steps or -s]")
			exit()
		end
		if dt == nothing
			ks = top.bondKs
			ms = top.masses
			dt = (1 / (sqrt(maximum(ks) / minimum(ms)))) * 0.1
			println("setting deltat to $(dt) based on bond spring constants and masses...")
		end

		println("setting velocities to temperature...")
		vel = zeros((top.nAtoms, 3))
		# https://physics.stackexchange.com/questions/159674/velocity-maxwell-boltzmann-distribution-for-dummies
		for i in 1:size(vel, 1)
			mi = top.masses[i]
			# sample from Maxwell-Boltzmann distribution
			vel[i, 1] = sqrt((kb * temp) / (2 * mi)) * randn()
			vel[i, 2] = sqrt((kb * temp) / (2 * mi)) * randn()
			vel[i, 3] = sqrt((kb * temp) / (2 * mi)) * randn()
		end
		initiateReport(traj, logfile)
		if constrainbonds == true
			println("running Langevin dynamics with constrained bonds...")
			bondLength = norm(xyz[top.bondIdx[1, :], :] .- xyz[top.bondIdx[2, :], :])
			nBonds = size(top.bondIdx, 1)
			λs = zeros(nBonds) # initial guess for Lagrange multipliers
			for frameNumber=1:steps
				if frameNumber % loginterval == 0
					println("step $(frameNumber) of $(steps)")
					report = true
				else
					report = false
				end
				λs = integrateLangevinSHAKE!(top, xyz, vel, dt, temp, collisionFreq, frameNumber, report, traj, logfile, bondLength, λs)
			end
		else
			println("running Langevin dynamics...")
			for frameNumber=1:steps
				if frameNumber % loginterval == 0
					println("step $(frameNumber) of $(steps)")
					report = true
				else
					report = false
				end
				integrateLangevin!(top, xyz, vel, dt, temp, collisionFreq, frameNumber, report, traj, logfile)
			end
		end
	else
		error("pick 'min' (steepest descent), 'vv' (velocity Verlet) integrator, or 'langevin' (Langevin dynamics (beta))")
	end
	
	println("Crawler has crawled. Have a nice day.")
end


main()
