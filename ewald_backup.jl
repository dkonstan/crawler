function getEwaldForces(top::Topology, xyz::Array{Float64})::Array{Float64}
	forces::Array{Float64} = zeros((top.nAtoms, 3))
	bondForces::Array{Float64} = zeros((top.nAtoms, 3))
	shortRangeForces::Array{Float64} = zeros((top.nAtoms, 3))
	longRangeForces::Array{Float64} = zeros((top.nAtoms, 3))
	tryUntil::Int64 = 4 # 5 box sizes brings the energy to within 5 decimal places
	tolerance::Float64 = 1e-24
	L::Float64 = top.box[1]
	rc::Float64 = L / 2
	ns::Array{Int64} = getCombinations(-tryUntil, tryUntil, 3)
	λ::Float64 = 1.0

	α::Float64 = 3.5 / (top.box[1] / 2)
	σ::Float64 = 1 / (sqrt(2) * α)

	# no force due to self-energy correction

	shortRangeForcesBefore = copy(shortRangeForces)
	# short-range forces
	for l::Int64 in 1:size(ns)[1]

		for i in 1:top.nAtoms
			for j in 1:top.nAtoms
				ci::Float64 = top.charges[i]
				cj::Float64 = top.charges[j]
				idx::Array{Array{Int64}} = top.bondIdx
				bonded::Bool = false

				# check if i/j pair is bonded
				for k::Int64 in 1:length(idx)
					if (i == idx[k][1] && j == idx[k][2]) || (j == idx[k][1] && i == idx[k][2])
						bonded = true
						break
					end
				end

				if bonded == false && !(i == j && ns[l, :] == [0.0, 0.0, 0.0])
					ξ = norm(xyz[i, :] - xyz[j, :] + L * ns[l, :])
					rx = xyz[i, 1] - xyz[j, 1] + L * ns[l, 1]
					ry = xyz[i, 2] - xyz[j, 2] + L * ns[l, 2]
					rz = xyz[i, 3] - xyz[j, 3] + L * ns[l, 3]
					Ξ = ξ^2
					if ξ < rc - λ
						shortRangeForces[i, :] += 0.5 * ((exp(-Ξ/(2 * σ^2)) * (sqrt(2 / pi)) * ci * cj * [rx, ry, rz]) / (σ * Ξ) + ((ci * cj * [rx, ry, rz] 
															* erfc(ξ / (sqrt(2) * σ))) / (Ξ^(3 / 2))))
					elseif ξ < rc && ξ > rc - λ
						shortRangeForces[i, :] += (ci * cj * [rx, ry, rz] * (
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


	volume::Float64 = L^3
	ms = ns

	longRangeForcesBefore = copy(longRangeForces)
	for l in 1:size(ms)[1]
		if ms[l, :] == [0.0, 0.0, 0.0]
			continue
		end

		# kVector::Array{Float64} = ms[l, 1] * b1 + ms[l, 2] * b2 + ms[l, 3] * b3
		kVector::Array{Float64} = (2 * pi) * (ms[l, :] ./ L)
		# println(kVector)
		k = norm(kVector)

		# strucFactorS = conj(strucFactor)
		strucFactor::ComplexF64 = 0.0 + 0.0 * im
		for i in 1:top.nAtoms
			strucFactor += top.charges[i] * exp(im * dot(kVector, xyz[i, :]))
		end

		strucFactorS = conj(strucFactor)

		for i in 1:top.nAtoms
			dSdr::Array{ComplexF64} = im * kVector * top.charges[i] * exp(im * dot(kVector, xyz[i, :]))
			dSsdr::Array{ComplexF64} = conj(dSdr)
			longRangeForces[i, :] += -((4 * pi) / volume) * ((exp(-σ^2 * k^2 / 2)) / k^2) .* (strucFactorS .* dSdr + strucFactor .* dSsdr)
		end
	end

	# this roughly cancels out long range forces due to bonded interactions (Tuckerman pg. 663)
	# need a faster way to get all of i's bonded neighbors than this but this works for now
	# isdone = Int64[]
	idx = top.bondIdx
	for b::Int64 in 1:size(idx)[1]
		for i in 1:top.nAtoms
			# println("i: ", i)
			for j in 1:top.nAtoms
				if ((idx[b][1] == i && idx[b][2] == j) || (idx[b][1] == j && idx[b][2] == i)) # && !(i in isdone)
					# println("firing j: ", j)
					ci::Float64 = top.charges[i]
					cj::Float64 = top.charges[j]
					rx::Float64 = xyz[i, 1] - xyz[j, 1]
					ry::Float64 = xyz[i, 2] - xyz[j, 2]
					rz::Float64 = xyz[i, 3] - xyz[j, 3]
					rij::Float64 = norm(xyz[i, :] - xyz[j, :])
					bondForces[i, :] += -((2 * exp(-rij^2 * α^2) * ci * cj * [rx, ry, rz] * α) / (sqrt(pi) * rij^2)) + ((ci * cj * [rx, ry, rz] * erf(rij * α)) / (rij^3))
				end
			end
		end
	end

	return shortRangeForces + longRangeForces - bondForces
end
