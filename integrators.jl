include("io.jl")
include("forces.jl")
include("structs.jl")
include("constraints.jl")

using LinearAlgebra
using Printf

function minimizeEnergy!(top::Topology, xyz::Array{Float64}, nSteps::Int64, outFile::String)
	#=
	gradient descent, Barzilai–Borwein method of choosing step size
	https://en.wikipedia.org/wiki/Gradient_descent
	nSteps = 0 means run until stepSize is NaN (no difference between old and new forces)
	nSteps = something means run for a MAX of something steps (may need fewer)
	=#
	forces = getForces(top, xyz)

	if nSteps == 0
		nSteps = typemax(Int64)  # 9223372036854775807 should be enough for ANY system :-)
	end

	initialStep = 0.001  # just to get things going
	stepSize = initialStep
	for i in 1:nSteps
		maxForce = norm(forces[1, :])
		maxForceIdx = 1
		for j in 2:top.nAtoms
			newForceMagnitude = norm(forces[j, :])
			if newForceMagnitude > maxForce
				maxForce = newForceMagnitude
				maxForceIdx = j
			end
		end

		energies::Energies = getEnergies(top, xyz)
		totalEnergy = energies.bond + energies.angle + energies.dihedral # + energies.coulomb

		@printf("step %d/%d, step size: %.2g, F_max: %.2g on atom %.2g, V(x): %.2g\n",
			    i, nSteps, stepSize, maxForce, maxForceIdx, totalEnergy)
		newCoords = xyz .+ (stepSize .* forces)
		newForces = getForces(top, newCoords)

		# matrix multiplication in numerator
		stepSize = norm(transpose(newCoords .- xyz) * (newForces .- forces)) / (norm(newForces .- forces)^2)
		if isnan(stepSize)
			break
		else
			forces = newForces
			xyz = newCoords
		end
	end

	println("minimization complete, writing coordinates to file $(outFile)...")
	dumpCoordinates(top, xyz, outFile)
end

function integrateVelocityVerlet!(top::Topology, xyz::Array{Float64}, vel::Array{Float64}, dt::Float64, frameNumber::Int64, report::Bool, trajectoryFile::String, reportFile::String)
	forces = getForces(top, xyz)

	for i in 1:top.nAtoms
		xyz[i, :] = xyz[i, :] .+ dt * vel[i, :] .+ ((dt^2) / (2 * top.masses[i])) .* forces[i, :]
	end

	newForces = getForces(top, xyz)

	for i in 1:top.nAtoms
		vel[i, :] = vel[i, :] .+ ((1 / top.masses[i]) .* (forces[i, :] .+ newForces[i, :]) .* dt) / 2
	end

	if report == true
		energies::Energies = getEnergies(top, xyz)
		kineticEnergy = getKineticEnergy(top, vel)
		reportProgress(top, xyz, vel, energies, kineticEnergy, frameNumber, trajectoryFile, reportFile)
	end
end

function integrateVelocityVerletSHAKE!(top::Topology, xyz::Array{Float64}, vel::Array{Float64}, dt::Float64,
									   frameNumber::Int64, report::Bool, trajectoryFile::String, reportFile::String,
									   bondLength::Float64, λs::Array{Float64})
	forces = getForces(top, xyz)

	bondLength = norm(xyz[top.bondIdx[1, 1], :] - xyz[top.bondIdx[1, 2], :])
	# regular velocity Verlet position update
	for i in 1:top.nAtoms
		xyz[i, :] = xyz[i, :] .+ dt * vel[i, :] .+ ((dt^2) / (2 * top.masses[i])) .* forces[i, :]
		vel[i, :] = vel[i, :] .+ ((1 / top.masses[i]) .* forces[i, :] * (dt / 2))  # normal but dt / 2 velocity update
	end

	# adjust positions to fit bond length constraints
	λs = shakePositions!(top, xyz, bondLength, λs)  # updates positions and λs

	newForces = getForces(top, xyz)

	# adjust velocities to fit bond length constraints, step 2
	shakeVelocities!(top, xyz, vel, λs, dt, newForces)
	if report == true
		energies::Energies = getEnergies(top, xyz)
		kineticEnergy = getKineticEnergy(top, vel)
		reportProgress(top, xyz, vel, energies, kineticEnergy, frameNumber, trajectoryFile, reportFile)
	end

	return λs
end

function integrateLangevin!(top::Topology, xyz::Array{Float64}, vel::Array{Float64}, dt::Float64, temp::Float64, collisionFreq::Float64,
						   frameNumber::Int64, report::Bool, trajectoryFile::String, reportFile::String)
	forces = getForces(top, xyz)
	γ = collisionFreq
	kb = 3.167e-6  # Boltzmann's constant in Hartree / K (https://physics.stackexchange.com/questions/635767/boltzmann-constant-in-atomic-units)

	θt = zeros((top.nAtoms, 3))
	ξt = zeros((top.nAtoms, 3))
	At = zeros((top.nAtoms, 3))
	for i in 1:top.nAtoms
		mi = top.masses[i]
		σi = sqrt(2 * kb * temp * γ / mi)
		ξt[i, :] = [randn(), randn(), randn()]
		θt[i, :] = [randn(), randn(), randn()]
		At[i, :] = (1 / 2) * dt^2 * (forces[i, :] .- γ * vel[i, :]) .+ σi * dt^(3 / 2) * ((1 / 2) * ξt[i, :] .+ (1 / (2 * sqrt(3))) * θt[i, :])
		xyz[i, :] = xyz[i, :] .+ dt * vel[i, :] .+ At[i, :]
	end

	newForces = getForces(top, xyz)

	for i in 1:top.nAtoms
		mi = top.masses[i]
		σi = sqrt(2 * kb * temp * γ * mi)
		vel[i, :] = vel[i, :] .+ ((1 / 2) * dt * (newForces[i, :] .+ forces[i, :])) .- (dt * γ * vel[i, :]) .+ (σi * sqrt(dt) * ξt[i, :]) .- (γ * At[i, :])
	end

	for i in 1:top.nAtoms
		mi = top.masses[i]

	end

	if report == true
		energies::Energies = getEnergies(top, xyz)
		kineticEnergy = getKineticEnergy(top, vel)
		reportProgress(top, xyz, vel, energies, kineticEnergy, frameNumber, trajectoryFile, reportFile)
	end

end


function integrateLangevinSHAKE!(top::Topology, xyz::Array{Float64}, vel::Array{Float64}, dt::Float64, temp::Float64, collisionFreq::Float64,
						        frameNumber::Int64, report::Bool, trajectoryFile::String, reportFile::String, bondLength::Float64, λs::Array{Float64})


	forces = getForces(top, xyz)

	γ = collisionFreq
	kb = 3.167e-6  # Boltzmann's constant in Hartree / K (https://physics.stackexchange.com/questions/635767/boltzmann-constant-in-atomic-units)

	θt = zeros((top.nAtoms, 3))
	ξt = zeros((top.nAtoms, 3))
	At = zeros((top.nAtoms, 3))
	for i in 1:top.nAtoms
		mi = top.masses[i]
		σi = sqrt(2 * kb * temp * γ / mi)
		ξt[i, :] = [randn(), randn(), randn()]
		θt[i, :] = [randn(), randn(), randn()]
		At[i, :] = (1 / 2) * dt^2 * (forces[i, :] .- γ * vel[i, :]) .+ σi * dt^(3 / 2) * ((1 / 2) * ξt[i, :] .+ (1 / (2 * sqrt(3))) * θt[i, :])

		# initial position update by full delta-t before Lagrange correction
		xyz[i, :] = xyz[i, :] .+ dt * vel[i, :] .+ At[i, :]

		# update velocities halfway by delta-t / 2 before Lagrange correction
		At[i, :] = (1 / 2) * (dt / 2)^2 * (forces[i, :] .- γ * vel[i, :]) .+ σi * (dt / 2)^(3 / 2) * ((1 / 2) * ξt[i, :] .+ (1 / (2 * sqrt(3))) * θt[i, :])
		# adapted Langevin adjustment to be delta-t / 2
		vel[i, :] = vel[i, :] .+ ((dt / 2) * forces[i, :]) .- ((dt / 2) * γ * vel[i, :]) .+ (σi * sqrt(dt / 2) * ξt[i, :]) .- (γ * At[i, :])
	end

	# adjust positions to fit bond length constraints (Lagrange correction)
	λs = shakePositions!(top, xyz, bondLength, λs)  # updates positions and λs

	newForces = getForces(top, xyz)

	σis = zeros(top.nAtoms)
	for i in 1:top.nAtoms
		mi = top.masses[i]
		σis[i] = sqrt(2 * kb * temp * γ / mi)
		ξt[i, :] = [randn(), randn(), randn()]
		θt[i, :] = [randn(), randn(), randn()]

		# At for the second half of the velocity adjustment, delta-t / 2
		At[i, :] = (1 / 2) * (dt / 2)^2 * (newForces[i, :] .- γ * vel[i, :]) .+ σis[i] * (dt / 2)^(3 / 2) * ((1 / 2) * ξt[i, :] .+ (1 / (2 * sqrt(3))) * θt[i, :])
	end

	# adjust velocities to fit bond length constraints, step 2, and do second half of Langevin adjustment within the function
	shakeVelocitiesLangevin!(top, xyz, vel, λs, dt, newForces, At, γ, ξt, σis)

	if report == true
		energies::Energies = getEnergies(top, xyz)
		kineticEnergy = getKineticEnergy(top, vel)
		reportProgress(top, xyz, vel, energies, kineticEnergy, frameNumber, trajectoryFile, reportFile)
	end

	return λs
end
