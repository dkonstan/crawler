include("IO.jl")
include("util.jl")
include("structs.jl")
using LinearAlgebra
using SpecialFunctions  # erf, erfc

function getEwaldEnergy(top::Topology, xyz)::Float64
	# cubic box only for now
	tryUntil = 2
	L = top.box[1]
	rc = L / 2.0
	α = 3.5 / rc
	σ = 1 / (sqrt(2) * α)

	λ = 1.0
	ns = getCombinations(-tryUntil, tryUntil, 3)

	selfEnergy = 0.0
	@inbounds for i in 1:top.nAtoms
		selfEnergy += (1 / (sqrt(2 * pi) * σ)) * (top.charges[i] * top.charges[i])
	end

	shortRangeEnergy = 0.0
	@inbounds for l in 1:size(ns, 1)
		for i in 1:top.nAtoms
			for j in i:top.nAtoms
				ci = top.charges[i]
				cj = top.charges[j]

				idx = top.bondIdx
				bonded  = false

				# check if i/j pair is bonded
				for k in 1:size(idx, 1)
					if (i == idx[k, 1] && j == idx[k, 2]) || (j == idx[k, 1] && i == idx[k, 2])
						bonded = true
						break
					end
				end

				@views if bonded == false && !(i == j && ns[l, :] == [0.0, 0.0, 0.0])
					@views r = norm(xyz[i, :] .- xyz[j, :] .+ L .* ns[l, :])
					if r < rc - λ
						shortRangeEnergy += ((ci * cj) / r) * erfc(r / (sqrt(2) * σ))
					elseif r < rc && r > rc - λ
						# switching function turns on
						shortRangeEnergy += ((ci * cj) * (rc - r)^2 * (-2 * rc + 2 * r + 3 * λ) * erfc(r / (sqrt(2) * σ))) / (2 * r * λ^3)
					end
				end
			end
		end
	end


	tryUntil = 4
	ns = getCombinations(-tryUntil, tryUntil, 3)
	ms = ns

	# long-range energy
	longRangeEnergy = 0.0
	volume::Float64 = L^3

	@inbounds for l in 1:size(ms)[1]  # no k = [0, 0, 0]
		@views if ms[l, :] == [0.0, 0.0, 0.0]
			continue
		end

		# define reciprocal space vectors
		kVector = (2 * pi) * (ms[l, :] / L)
		k = norm(kVector)

		# check this equal 1 good
		# println(exp(-im * dot(kVector, L * [ms[l, 1], ms[l, 2], ms[l, 3]])))
		strucFactor = 0.0 + 0.0 * im
		for i in 1:top.nAtoms
			@views strucFactor += top.charges[i] * complexExp(im * dot(kVector, xyz[i, :]))
		end

		longRangeEnergy += ((4 * pi) / (volume)) * ((exp(-σ^2 * k^2 / 2) / k^2) * (abs2(strucFactor))) # (abs(strucFactor)^2) cab be replaced with abs2(strucFactor)
	end


	# this roughly cancels out long range energy due to bonded interactions (Tuckerman pg. 663)
	bondEnergy::Float64 = 0.0
	idx = top.bondIdx
	@inbounds for b::Int64 in 1:size(idx, 1)
		for i in 1:top.nAtoms
			for j in 1:top.nAtoms
				@views if ((idx[b, 1] == i && idx[b, 2] == j) || (idx[b, 1] == j && idx[b, 2] == i)) # && !(i in isdone)
					# println("firing j: ", j)
					ci = top.charges[i]
					cj = top.charges[j]
					@views rij = norm(xyz[i, :] .- xyz[j, :])
					bondEnergy += (ci * cj * erf(rij / (sqrt(2) * σ))) / rij
				end
			end
		end
	end

	ewaldEnergy = shortRangeEnergy + longRangeEnergy - selfEnergy - bondEnergy

	return ewaldEnergy
end


function getEwaldForces(top::Topology, xyz)::Array{Float64}
	forces = zeros((top.nAtoms, 3))
	bondForces = zeros((top.nAtoms, 3))
	shortRangeForces = zeros((top.nAtoms, 3))
	longRangeForces = zeros((top.nAtoms, 3))
	tryUntil = 2 # 5 box sizes brings the energy to within 5 decimal places
	L = top.box[1]
	rc = L / 2
	ns = getCombinations(-tryUntil, tryUntil, 3)
	λ = 1.0

	α = 3.5 / rc
	σ = 1 / (sqrt(2) * α)

	# no force due to self-energy correction

	# short-range forces
	@inbounds for l in 1:size(ns, 1)
		for i in 1:top.nAtoms
			for j in 1:top.nAtoms

				ci = top.charges[i]
				cj = top.charges[j]
				idx = top.bondIdx
				bonded = false

				# check if i/j pair is bonded
				for k in 1:size(idx, 1)
					@views if (i == idx[k, 1] && j == idx[k, 2]) || (j == idx[k, 1] && i == idx[k, 2])
						bonded = true
						break
					end
				end

				@views if bonded == false && !(i == j && ns[l, :] == [0.0, 0.0, 0.0])
					@views ξ = norm(xyz[i, :] .- xyz[j, :] .+ L .* ns[l, :])
					rx = xyz[i, 1] .- xyz[j, 1] + L .* ns[l, 1]
					ry = xyz[i, 2] .- xyz[j, 2] + L .* ns[l, 2]
					rz = xyz[i, 3] .- xyz[j, 3] + L .* ns[l, 3]
					Ξ = ξ^2
					if ξ < rc - λ
						shortRangeForces[i, :] .+= ((exp(-Ξ/(2 * σ^2)) * (sqrt(2 / pi)) * ci * cj * [rx, ry, rz]) / (σ * Ξ) + ((ci * cj * [rx, ry, rz] 
															* erfc(ξ / (sqrt(2) * σ))) / (Ξ^(3 / 2))))
					elseif ξ < rc && ξ > rc - λ
						shortRangeForces[i, :] .+= (ci * cj * [rx, ry, rz] * (
														(
															(exp(-ξ / (2 * σ^20)) * sqrt(2 / pi) * (2 * rc - 3 * λ - 2 * Ξ) * (rc - Ξ)^2) / (ξ * σ)
														) +
												  		(
												  			((2 * rc - 3 * λ - 2 * Ξ) * (rc - Ξ)^2 * erfc(Ξ / (sqrt(2) * σ))) / (ξ^(3 / 2))
												  		) +
												  		(
												  			(6 * (-rc * Ξ) * (-rc * λ + Ξ) * erfc(Ξ / (sqrt(2) * σ))) / (Ξ^2)
												  		)
												  					)
													) / (2 * λ^3)
					end
				end
			end
		end
	end


	volume = L^3

	tryUntil = 4
	ns = getCombinations(-tryUntil, tryUntil, 3)
	ms = ns

	@inbounds for l in 1:size(ms, 1)
		@views if ms[l, :] == [0.0, 0.0, 0.0]
			continue
		end

		# kVector::Array{Float64} = ms[l, 1] * b1 + ms[l, 2] * b2 + ms[l, 3] * b3
		@views kVector = (2 * pi) .* (ms[l, :] ./ L)
		k = norm(kVector)

		# strucFactorS = conj(strucFactor)
		strucFactor = 0.0 + 0.0 * im
		for i in 1:top.nAtoms
			@views strucFactor += top.charges[i] * complexExp(im * dot(kVector, xyz[i, :]))
		end

		strucFactorS = conj(strucFactor)

		for i in 1:top.nAtoms
			# dSdr = im * kVector * top.charges[i] * exp(im * dot(kVector, xyz[i, :]))
			@views dSdr = im * kVector * top.charges[i] * complexExp(im * dot(kVector, xyz[i, :]))
			dSsdr = conj(dSdr)
			longRangeForces[i, :] .+= -((4 * pi) / volume) * ((exp(-σ^2 * k^2 / 2)) / k^2) .* (strucFactorS .* dSdr + strucFactor .* dSsdr)
		end
	end

	# this roughly cancels out long range forces due to bonded interactions (Tuckerman pg. 663)
	# need a faster way to get all of i's bonded neighbors than this but this works for now
	# isdone = Int64[]
	idx = top.bondIdx
	@inbounds for b in 1:size(idx, 1)
		for i in 1:top.nAtoms
			# println("i: ", i)
			for j in 1:top.nAtoms
				@views if ((idx[b, 1] == i && idx[b, 2] == j) || (idx[b, 1] == j && idx[b, 2] == i)) # && !(i in isdone)
					# println("firing j: ", j)
					ci = top.charges[i]
					cj = top.charges[j]
					rx = xyz[i, 1] - xyz[j, 1]
					ry = xyz[i, 2] - xyz[j, 2]
					rz = xyz[i, 3] - xyz[j, 3]
					@views rij = norm(xyz[i, :] .- xyz[j, :])
					bondForces[i, :] .+= -((2 * exp(-rij^2 * α^2) * ci * cj * [rx, ry, rz] * α) / (sqrt(pi) * rij^2)) + ((ci * cj * [rx, ry, rz] * erf(rij * α)) / (rij^3))
				end
			end
		end
	end

	return shortRangeForces .+ longRangeForces .- bondForces
end

# top = readTopology("top.txt")
# xyz = readCoordinates("crd.xyz")

# println(typeof(xyz))
# exit()
# energy = getEwaldEnergy(top, xyz)
# @time energy = getEwaldEnergy(top, xyz)
# println("energy", " ", energy)

# forces = getEwaldForces(top, xyz)
# @time forces = getEwaldForces(top, xyz)
# for i in 1:top.nAtoms
# 	println(forces[i, :])
# end
