include("IO.jl")
include("structs.jl")
include("util.jl")
include("ewald.jl")
using LinearAlgebra


function getPotentialEnergy(energies::Energies)
	return energies.bond + energies.angle + energies.dihedral + energies.coulomb
end

function getBondEnergy(top::Topology, xyz)::Float64
	bondEnergy = 0.0
	for i in 1:size(top.bondIdx, 1)
		atom1Idx = top.bondIdx[i, 1]
		atom2Idx = top.bondIdx[i, 2]
		k = top.bondKs[i]
		bondLength = top.bondLengths[i]
		@views r1r2 = xyz[atom1Idx, :] - xyz[atom2Idx, :]
		r1r2 = pbcAdjust(top, r1r2)

		distance = norm(r1r2)
		stretch = distance - bondLength
		bondEnergy += 0.5 * k * stretch^2
	end

	return bondEnergy
end


function getBondForces(top::Topology, xyz)::Array{Float64}
	bondForces = zeros((top.nAtoms, 3))

	for i in 1:size(top.bondIdx, 1)
		atom1Idx = top.bondIdx[i, 1]
		atom2Idx = top.bondIdx[i, 2]
		k = top.bondKs[i]
		bondLength = top.bondLengths[i]

		@views r1r2 = xyz[atom1Idx, :] .- xyz[atom2Idx, :]
		r1r2 = pbcAdjust(top, r1r2)
		distance = norm(r1r2)
		stretch = distance - bondLength
		unitVector = r1r2 / distance

		bondForces[atom1Idx, :] += -k .* stretch .* unitVector
		bondForces[atom2Idx, :] += -k .* stretch .* -(unitVector)
	end

	return bondForces
end


function getAngleEnergy(top::Topology, xyz)::Float64
	angleEnergy = 0.0

	for i in 1:size(top.angleIdx, 1)
		@views atom1 = xyz[top.angleIdx[i, 1], :]
		@views atom2 = xyz[top.angleIdx[i, 2], :]  # central atom
		@views atom3 = xyz[top.angleIdx[i, 3], :]
		atom1Idx = top.angleIdx[i, 1]
		atom2Idx = top.angleIdx[i, 2]
		atom3Idx = top.angleIdx[i, 3]
		kTheta = top.angleThetas[i]

		vec12 = atom1 .- atom2
		vec32 = atom3 .- atom2

		vec12 = pbcAdjust(top, vec12)
		vec32 = pbcAdjust(top, vec32)

		vec12 /= norm(vec12)
		vec32 /= norm(vec32)

		dot1232 = dot(vec12, vec32)
		# correct numerical errors
		if dot1232 < -1.0
			dot1232 = -1.0
		elseif dot1232 > 1.0
			dot1232 = 1.0
		end

		angle = acos(dot1232)
		angle0 = top.angleValues[i]

		angleEnergy += 0.5 * kTheta * (angle - angle0)^2
	end

	return angleEnergy
end

function getAngleForces(top::Topology, xyz)::Array{Float64}
	# https://salilab.org/modeller/9v6/manual/node436.html
	angleForces = zeros((top.nAtoms, 3))
	limitAngle = 0.001

	for i in 1:size(top.angleIdx, 1)
		@views atom1 = xyz[top.angleIdx[i, 1], :]
		@views atom2 = xyz[top.angleIdx[i, 2], :]  # central atom
		@views atom3 = xyz[top.angleIdx[i, 3], :]
		atom1Idx = top.angleIdx[i, 1]
		atom2Idx = top.angleIdx[i, 2]
		atom3Idx = top.angleIdx[i, 3]
		kTheta = top.angleThetas[i]

		vec12 = pbcAdjust(top, atom1 .- atom2)
		vec32 = pbcAdjust(top, atom3 .- atom2)

		vec12 /= norm(vec12)
		vec32 /= norm(vec32)

		dot1232 = dot(vec12, vec32)

		# correct numerical errors
		if dot1232 < -1.0
			dot1232 = -1.0
		elseif dot1232 > 1.0
			dot1232 = 1.0
		end

		angle = acos(dot1232)
		angle0 = top.angleValues[i]

		if abs(angle) < limitAngle || abs(pi - angle) < limitAngle
			forcei = 0.0
			forcek = 0.0
		else
			prefactor = (1 / sqrt(1 - cos(angle)^2))
			rijVec = pbcAdjust(top, atom1 .- atom2)
			rkjVec = pbcAdjust(top, atom3 .- atom2)
			rij = norm(rijVec)
			rkj = norm(rkjVec)
			forcei = prefactor * (1 / rij) * ((rijVec / rij) * cos(angle) .- (rkjVec / rkj))
			forcek = prefactor * (1 / rkj) * ((rkjVec / rkj) * cos(angle) .- (rijVec / rij))
		end

		force1 = -kTheta * (angle - angle0) .* forcei
		force3 = -kTheta * (angle - angle0) .* forcek
		force2 = -(force1 .+ force3)

		angleForces[atom1Idx, :] .+= force1
		angleForces[atom2Idx, :] .+= force2
		angleForces[atom3Idx, :] .+= force3
	end

	return angleForces
end


function getDihedralEnergy(top::Topology, xyz)::Float64
	dihedralEnergy = 0.0

	for i in 1:size(top.dihedralIdx, 1)
		@views atom1 = xyz[top.dihedralIdx[i, 1], :]
		@views atom2 = xyz[top.dihedralIdx[i, 2], :]
		@views atom3 = xyz[top.dihedralIdx[i, 3], :]
		@views atom4 = xyz[top.dihedralIdx[i, 4], :]
		atom1Idx = top.dihedralIdx[i, 1]
		atom2Idx = top.dihedralIdx[i, 2]
		atom3Idx = top.dihedralIdx[i, 3]
		atom4Idx = top.dihedralIdx[i, 4]
		kTheta = top.dihedralThetas[i]

		vec12 = pbcAdjust(top, atom1 .- atom2)
		vec32 = pbcAdjust(top, atom3 .- atom2)
		vec23 = -vec32
		vec43 = pbcAdjust(top, atom4 .- atom3)

		normal123 = cross(vec12, vec32)
		normal234 = cross(vec43, vec23)

		normal123 /= norm(normal123)
		normal234 /= norm(normal234)

		dott = dot(normal123, normal234)

		if dott > 1.0
			dott = 1.0
		elseif dott < -1.0
			dott = -1.0
		end

		dihedralAngle = acos(dott)
		dihedralAngle0 = top.dihedralValues[i]

		dihedralEnergy += 0.5 * kTheta * (dihedralAngle - dihedralAngle0)^2
	end

	return dihedralEnergy
end

function getDihedralForces(top::Topology, xyz)::Array{Float64}
	# https://salilab.org/modeller/9v6/manual/node436.html
	dihedralForces = zeros((top.nAtoms, 3))

	for i in 1:size(top.dihedralIdx, 1)

		@views atom1 = xyz[top.dihedralIdx[i, 1], :]
		@views atom2 = xyz[top.dihedralIdx[i, 2], :]
		@views atom3 = xyz[top.dihedralIdx[i, 3], :]
		@views atom4 = xyz[top.dihedralIdx[i, 4], :]
		atom1Idx = top.dihedralIdx[i, 1]
		atom2Idx = top.dihedralIdx[i, 2]
		atom3Idx = top.dihedralIdx[i, 3]
		atom4Idx = top.dihedralIdx[i, 4]
		kTheta = top.dihedralThetas[i]

		vec12 = pbcAdjust(top, atom1 .- atom2)
		vec32 = pbcAdjust(top, atom3 .- atom2)
		vec23 = -vec32
		vec43 = pbcAdjust(top, atom4 .- atom3)

		normal123 = cross(vec12, vec32)
		normal234 = cross(vec43, vec23)

		norm123 = norm(normal123)
		norm234 = norm(normal234)
		# println(norm123, " ", norm234)
		if norm123 < 0.00001 || norm234 < 0.00001
			# happens if the polymer is straight, avoid a singularity
			force1 = [0.0, 0.0, 0.0]
			force2 = [0.0, 0.0, 0.0]
			force3 = [0.0, 0.0, 0.0]
			force4 = [0.0, 0.0, 0.0]
		else
			dotprod = dot(normal123 / norm123, normal234 / norm234)
			if dotprod > 1.0
				dotprod = 1.0
			elseif dotprod < -1.0
				dotprod = -1.0
			end
			dihedralAngle = acos(dotprod)
			dihedralAngle0 = top.dihedralValues[i]

			rijVec = pbcAdjust(top, atom2 .- atom1)
			rkjVec = pbcAdjust(top, atom2 .- atom3)
			rklVec = pbcAdjust(top, atom4 .- atom3)
			rmjVec = cross(rijVec, rkjVec)
			rnkVec = cross(rkjVec, rklVec)

			rij = norm(rijVec)
			rkj = norm(rkjVec)
			rmj = norm(rmjVec)
			rnk = norm(rnkVec)
			forcei = (rkj / rmj^2) * rmjVec
			forcel = -(rkj / rnk^2) * rnkVec

			forcej = ((dot(rijVec, rkjVec) / rkj^2) - 1.0) * forcei - (dot(rklVec, rkjVec) / rkj^2) * forcel
			forcek = ((dot(rklVec, rkjVec) / rkj^2) - 1.0) * forcel - (dot(rijVec, rkjVec) / rkj^2) * forcei

			force1 = -kTheta * (dihedralAngle - dihedralAngle0) .* forcei
			force2 = -kTheta * (dihedralAngle - dihedralAngle0) .* forcej
			force3 = -kTheta * (dihedralAngle - dihedralAngle0) .* forcek
			force4 = -kTheta * (dihedralAngle - dihedralAngle0) .* forcel
		end

		dihedralForces[atom1Idx, :] .+= force1
		dihedralForces[atom2Idx, :] .+= force2
		dihedralForces[atom3Idx, :] .+= force3
		dihedralForces[atom4Idx, :] .+= force4
	end

	return dihedralForces
end


function getVDWEnergy(top::Topology, xyz::Array{Float64})::Float64
	vdwEnergy = 0.0
	for i in 1:size(top.vdwIdx, 1)
		atom1Idx = top.vdwIdx[i, 1]
		atom2Idx = top.vdwIdx[i, 2]

		σ = top.vdwSigmas[i]
		ϵ = top.vdwEpsilons[i]
		@views r1r2 = xyz[atom1Idx, :] - xyz[atom2Idx, :]
		r1r2 = pbcAdjust(top, r1r2)
		r = norm(r1r2)

		vdwEnergy += 4 * ϵ * ((σ / r)^12 - (σ / r)^6)
	end
	return vdwEnergy
end

function getVDWForces(top::Topology, xyz::Array{Float64})
	vdwForces = zeros((top.nAtoms, 3))

	for i in 1:size(top.vdwIdx, 1)
		atom1Idx = top.vdwIdx[i, 1]
		atom2Idx = top.vdwIdx[i, 2]
		σ = top.vdwSigmas[i]
		ϵ = top.vdwEpsilons[i]

		@views r1r2 = xyz[atom1Idx, :] .- xyz[atom2Idx, :]
		r1r2 = pbcAdjust(top, r1r2)
		r = norm(r1r2)
		unitVector = r1r2 / r
		vdwForces[atom1Idx, :] += -4 * ϵ * (((6 * σ^6) / r^7) - ((12 * σ^12) / r^13)) * unitVector
		vdwForces[atom2Idx, :] += -4 * ϵ * (((6 * σ^6) / r^7) - ((12 * σ^12) / r^13)) * -unitVector
	end

	return vdwForces
end


function getForces(top::Topology, xyz)::Array{Float64}
	bondForces = getBondForces(top, xyz)
	angleForces = getAngleForces(top, xyz)
	dihedralForces = getDihedralForces(top, xyz)
	vdwForces = getVDWForces(top, xyz)
	# coulombForces = getEwaldForces(top, xyz)
	# vdwForces = zeros((top.nAtoms, 3))
	coulombForces = zeros((top.nAtoms, 3))
	return bondForces + angleForces .+ dihedralForces .+ coulombForces .+ vdwForces
end

function getEnergies(top::Topology, xyz)::Energies
	bondEnergy = getBondEnergy(top, xyz)
	angleEnergy = getAngleEnergy(top, xyz)
	dihedralEnergy = getDihedralEnergy(top, xyz)
	vdwEnergy = getVDWEnergy(top, xyz)
	# coulombEnergy = getEwaldEnergy(top, xyz)
	coulombEnergy = 0.0
	energies::Energies = Energies(bondEnergy, angleEnergy, dihedralEnergy, coulombEnergy, vdwEnergy)
	return energies
end


function getKineticEnergy(top::Topology, vel::Array{Float64})::Float64
	kineticEnergy = 0.0
	for i in 1:top.nAtoms
		@views kineticEnergy += 0.5 * top.masses[i] * dot(vel[i, :], vel[i, :])
	end

	return kineticEnergy
end
