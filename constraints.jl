include("structs.jl")
include("util.jl")

using LinearAlgebra

function sqNorm(a::Array{Float64}, b::Array{Float64})::Float64
	sqNorm = 0.0
	for i in 1:3
		sqNorm += (a[i] - b[i])^2
	end

	return sqNorm
end

function shakePositions!(top::Topology, xyz::Array{Float64}, bondLength::Float64, λs::Array{Float64})

	# Tuckerman,
	idx = top.bondIdx
	m = top.masses
	d = bondLength
	N = top.nAtoms
	nC = size(idx, 1)  # nConstraints

	# calculate constraints at original position
	σ = zeros((nC))
	for k in 1:nC
		i = idx[k, 1]
		j = idx[k, 2]
		pos1 = xyz[i, :]
		pos2 = xyz[j, :]
		σ[k] = sqNorm(pos1, pos2) - d^2
	end

	# calculate gradient of constraints at original position
	∇σ = zeros((nC, N, 3))  # gradients of constraints
	∇σMatrix = zeros((nC, nC))
	for k in 1:nC
		i = idx[k, 1]
		j = idx[k, 2]
		ri = xyz[i, :]
		rj = xyz[j, :]
		∇σ[k, i, :] = 2 * (ri .- rj)  # gradient of constraint wrt i
		∇σ[k, j, :] = 2 * (rj .- ri)  # gradient of constraint wrt j
	end

	# calculate initial adjustment based on guesses for Lagrange multipliers
	for i in 1:N
		adjustment = zeros(3)
		for k in 1:nC
			adjustment .+= (1 / m[i]) * (λs[k] * ∇σ[k, i, :])
		end
		xyz[i, :] = xyz[i, :] .+ adjustment
	end

	# now time to calculate the remaining correction to the position

	tolerance = 1e-10
	convergenceIndicator = 1.0  # just a random large value to get started
	while convergenceIndicator > tolerance
		# calculate constraints at adjusted position
		σ = zeros((nC))
		for k in 1:nC
			pos1 = xyz[idx[k, 1], :]
			pos2 = xyz[idx[k, 2], :]
			σ[k] = sqNorm(pos1, pos2) - d^2
		end

		# calculate gradient of constraints at adjusted position
		∇σ = zeros((nC, N, 3))  # gradients of constraints
		∇σMatrix = zeros((nC, nC))
		for k in 1:nC
			i = idx[k, 1]
			j = idx[k, 2]
			ri = xyz[i, :]
			rj = xyz[j, :]
			∇σ[k, i, :] = 2 * (ri .- rj)  # gradient of constraint wrt i
			∇σ[k, j, :] = 2 * (rj .- ri)  # gradient of constraint wrt j
		end

		for l in 1:nC
			for k in 1:nC
				for i in 1:N
				 	∇σMatrix[l, k] += (1 / m[i]) * dot(∇σ[l, i, :], ∇σ[k, i, :])
				end
			end
		end

		# solve matrix equation for all dλs at once (M-SHAKE algorithm)
		dλs = inv(∇σMatrix) * -σ

		λs += dλs
		# calculate adjustment based on dλs
		adjustmentMagnitudes = zeros(N)
		for i in 1:N
			adjustment = zeros(3)
			for k in 1:nC
				adjustment .+= (1 / m[i]) * (dλs[k] * ∇σ[k, i, :])
			end
			xyz[i, :] = xyz[i, :] .+ adjustment
			adjustmentMagnitudes = norm(adjustment)
		end
		convergenceIndicator = maximum(adjustmentMagnitudes)
	end

	return λs
end


function shakeVelocities!(top::Topology, xyz::Array{Float64}, vel::Array{Float64}, λs::Array{Float64}, dt::Float64, newForces::Array{Float64})
	idx = top.bondIdx
	m = top.masses
	N = top.nAtoms
	nC = size(idx, 1)  # nConstraints

	∇σ = zeros((nC, N, 3))  # gradients of constraints at time t0 + delta-t
	for k in 1:nC
		i = idx[k, 1]
		j = idx[k, 2]
		ri = xyz[i, :]
		rj = xyz[j, :]
		∇σ[k, i, :] = 2 * (ri .- rj)  # gradient of constraint wrt i
		∇σ[k, j, :] = 2 * (rj .- ri)  # gradient of constraint wrt j
	end

	for i in 1:N
		adjustment1 = zeros(3)
		for k in 1:nC
			adjustment1 .+= (1 / m[i]) * (1 / dt) * (λs[k] * ∇σ[k, i, :])
		end
		adjustment1 += ((dt / 2) * (1 / m[i]) * newForces[i, :])  # second dt / 2 adjustment
		vel[i, :] = vel[i, :] .+ adjustment1  # v'i
	end

	∇σMatrix = zeros((nC, nC))
	∇σVector = zeros(nC)
	for k in 1:nC
		for i in 1:N
			∇σVector[k] += dot(∇σ[k, i, :],  vel[i, :])
		end
	end
	for k in 1:nC
		for l in 1:nC
			for i in 1:N
				∇σMatrix[k, l] += dot(∇σ[k, i, :], (1 / m[i]) * ∇σ[l, i, :])
			end
		end
	end

	μs = inv(∇σMatrix) * -∇σVector

	for i in 1:N
		adjustment2 = zeros(3)
		for k in 1:nC
			adjustment2 .+= (1 / m[i]) * (μs[k] * ∇σ[k, i, :])  # Lagrange adjustment
		end
		vel[i, :] = vel[i, :] .+ adjustment2
	end
end


function shakeVelocitiesLangevin!(top::Topology, xyz::Array{Float64}, vel::Array{Float64}, λs::Array{Float64},
								  dt::Float64, newForces::Array{Float64},
								  At::Array{Float64}, γ::Float64, ξt::Array{Float64}, 
								  σis::Array{Float64})
	idx = top.bondIdx
	m = top.masses
	N = top.nAtoms
	nC = size(idx, 1)  # nConstraints

	∇σ = zeros((nC, N, 3))  # gradients of constraints at time t0 + delta-t
	for k in 1:nC
		i = idx[k, 1]
		j = idx[k, 2]
		ri = xyz[i, :]
		rj = xyz[j, :]
		∇σ[k, i, :] = 2 * (ri .- rj)  # gradient of constraint wrt i
		∇σ[k, j, :] = 2 * (rj .- ri)  # gradient of constraint wrt j
	end

	for i in 1:N
		adjustment1 = zeros(3)
		for k in 1:nC
			adjustment1 .+= (1 / m[i]) * (1 / dt) * (λs[k] * ∇σ[k, i, :])
		end
		# update velocities halfway by delta-t / 2 including Langevin adjustment
		adjustment1 += ((dt / 2) * (1 / m[i]) * newForces[i, :]) - ((dt / 2) * γ * vel[i, :]) .+ (σis[i] * sqrt(dt / 2) * ξt[i, :]) .- (γ * At[i, :])

		vel[i, :] = vel[i, :] .+ adjustment1  # v'i
	end

	∇σMatrix = zeros((nC, nC))
	∇σVector = zeros(nC)
	for k in 1:nC
		for i in 1:N
			∇σVector[k] += dot(∇σ[k, i, :],  vel[i, :])
		end
	end
	for k in 1:nC
		for l in 1:nC
			for i in 1:N
				∇σMatrix[k, l] += dot(∇σ[k, i, :], (1 / m[i]) * ∇σ[l, i, :])
			end
		end
	end

	μs = inv(∇σMatrix) * -∇σVector

	for i in 1:N
		adjustment2 = zeros(3)
		for k in 1:nC
			adjustment2 .+= (1 / m[i]) * (μs[k] * ∇σ[k, i, :])  # Lagrange adjustment is done last (of course)
		end
		vel[i, :] = vel[i, :] .+ adjustment2
	end
end
