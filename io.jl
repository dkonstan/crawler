include("structs.jl")
using StaticArrays

function readTopology(filename)

	masses = []
	bondIdx = []
	bondLengths = []
	bondKs = []
	angleIdx = []
	angleValues = []
	angleThetas = []
	dihedralIdx = []
	dihedralValues = []
	dihedralThetas = []
	charges = []
	vdwIdx = []
	vdwSigmas = []
	vdwEpsilons = []
	atomTypes = []
	box = [0.0, 0.0, 0.0]

	try
		open(filename, "r") do top
			contents = readlines(top)

			i = 1
			while i <= length(contents)
				if occursin("<atomTypes>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						# printlns(contents[i])
						atomTypeContents = split(contents[i], " ")
						atomType = atomTypeContents[2]
						push!(atomTypes, atomType)
						i += 1
					end
				elseif occursin("<masses>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						massContents = split(contents[i], " ")
						mass = parse(Float64, massContents[2])
						append!(masses, mass)
						i += 1
					end
				elseif occursin("<bonds>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						bondContents = split(contents[i], " ")
						bond1 = parse(Int64, bondContents[1])
						bond2 = parse(Int64, bondContents[2])
						bondLength = parse(Float64, bondContents[3])
						bondK = parse(Float64, bondContents[4])
						push!(bondIdx, [bond1, bond2])
						append!(bondLengths, bondLength)
						append!(bondKs, bondK)
						i += 1
					end
				elseif occursin("<angles>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						angleContents = split(contents[i], " ")
						angle1 = parse(Int64, angleContents[1])
						angle2 = parse(Int64, angleContents[2])
						angle3 = parse(Int64, angleContents[3])
						angleValue = parse(Float64, angleContents[4])
						angleTheta = parse(Float64, angleContents[5])
						push!(angleIdx, [angle1, angle2, angle3])
						append!(angleValues, angleValue)
						append!(angleThetas, angleTheta)
						i += 1
					end
				elseif occursin("<dihedrals>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						dihedralContents = split(contents[i], " ")
						dihedral1 = parse(Int64, dihedralContents[1])
						dihedral2 = parse(Int64, dihedralContents[2])
						dihedral3 = parse(Int64, dihedralContents[3])
						dihedral4 = parse(Int64, dihedralContents[4])
						dihedralValue = parse(Float64, dihedralContents[5])
						dihedralTheta = parse(Float64, dihedralContents[6])
						push!(dihedralIdx, [dihedral1, dihedral2, dihedral3, dihedral4])
						append!(dihedralValues, dihedralValue)
						append!(dihedralThetas, dihedralTheta)
						i += 1
					end
				elseif occursin("<charges>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						chargeContents = split(contents[i], " ")
						charge = parse(Float64, chargeContents[2])
						append!(charges, charge)
						i += 1
					end
				elseif occursin("<vdw>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						vdwContents = split(contents[i], " ")
						vdw1 = parse(Int64, vdwContents[1])
						vdw2 = parse(Int64, vdwContents[2])
						vdwSigma = parse(Float64, vdwContents[3])
						vdwEpsilon = parse(Float64, vdwContents[4])
						push!(vdwIdx, [vdw1, vdw2])
						append!(vdwSigmas, vdwSigma)
						append!(vdwEpsilons, vdwEpsilon)
						i += 1
					end
				elseif occursin("<box>", contents[i])
					i += 1
					while !occursin("<end>", contents[i])
						boxContents = split(contents[i], " ")
						boxX = parse(Float64, boxContents[1])
						boxY = parse(Float64, boxContents[2])
						boxZ = parse(Float64, boxContents[3])
						box = [boxX, boxY, boxZ]
						i += 1
					end
				end
				i += 1
			end
		end
	catch
		println("fix topology!")
		exit()
	end

	nAtoms = length(masses)

	bondIdx = Array{Int64}(hcat(bondIdx...)')
	angleIdx = Array{Int64}(hcat(angleIdx...)')
	dihedralIdx = Array{Int64}(hcat(dihedralIdx...)')
	vdwIdx = Array{Int64}(hcat(vdwIdx...)')

	topology = Topology(masses,
						bondIdx, bondKs, bondLengths,
						angleIdx, angleValues * (pi / 180), angleThetas,
						dihedralIdx, dihedralValues * (pi / 180), dihedralThetas,
						charges,
						vdwIdx, vdwSigmas, vdwEpsilons,
						nAtoms, atomTypes, box)

	if checkTopology(topology::Topology) == false
		println("fix topology!")
		exit()
	end

	return topology
end


function checkTopology(top)
	#=
	masses = []
	bondIdx = []
	bondLengths = []
	bondKs = []
	angleIdx = []
	angleValues = []
	angleThetas = []
	dihedralIdx = []
	dihedralValues = []
	dihedralThetas = []
	charges = []
	vdwSigma1s = []
	vdwSigma2s = []
	vdwEpsilons = []
	atomTypes = []
	box = [0.0, 0.0, 0.0]
	=#
	if length(top.masses) != length(top.charges)
		return false
	else
		return true
	end
end

function readCoordinates(filename)
	# xyz format
	xyz = []
	open(filename, "r") do top
		contents = readlines(top)

		i = 3  # first line is the number of atoms, then a comment line
		while i <= length(contents)
			try
				if !isletter(contents[i][1])
					break
				end
				xyzContents = split(contents[i], " ")
				x = parse(Float64, xyzContents[2])
				y = parse(Float64, xyzContents[3])
				z = parse(Float64, xyzContents[4])
				push!(xyz, [x, y, z])
				i += 1
			catch
				break
			end
		end
	end

	xyz = Array{Float64}(hcat(xyz...)')
	# xyz = MArray{Tuple{size(xyz, 1), 3}}(xyz)
	# println(xyz)
	# exit()
	return xyz # turn array of arrays into pure 2D array
end


function initiateReport(trajectoryFile, reportFile)
	open(trajectoryFile, "w+") do traj
	end

	open(reportFile, "w+") do report
		write(report, "frameNumber,total_energy,kinetic_energy,bond_energy,angle_energy,dihedral_energy,coulomb_energy,vdw_energy,vx\n")
	end
end

function reportProgress(top, xyz, vel, energies, kineticEnergy, frameNumber, trajectoryFile, reportFile)
	open(trajectoryFile, "a+") do traj
		write(traj, "$(top.nAtoms)\n")
		write(traj, "$(frameNumber)\n")
		for i in 1:top.nAtoms
			write(traj, "$(top.atomTypes[i]) $(xyz[i, 1]) $(xyz[i, 2]) $(xyz[i, 3])\n")
		end
	end

	total_energy = kineticEnergy + energies.bond + energies.angle + energies.dihedral + energies.coulomb + energies.vdw
	open(reportFile, "a+") do report
		write(report, "$(frameNumber),$(total_energy),$(kineticEnergy),$(energies.bond),$(energies.angle),$(energies.dihedral),$(energies.coulomb),$(energies.vdw),$(vel[1, 1])\n")
	end
end


function dumpCoordinates(top, xyz, reportFile)
	# xyz format
	open(reportFile, "w+") do report
		write(report, "$(top.nAtoms)\n")
		write(report, "comment line\n")
		for i in 1:top.nAtoms
			write(report, "$(top.atomTypes[i]) $(xyz[i, 1]) $(xyz[i, 2]) $(xyz[i, 3])\n")
		end
	end
end


