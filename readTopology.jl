using LinearAlgebra

function readTopology(filename)

	masses::Array{Float64} = []
	xyz::Array{Array{Float64}} = []
	bondIdx::Array{Array{Int64}} = []
	bondLengths::Array{Float64} = []
	bondKs::Array{Float64} = []
	angleIdx::Array{Array{Int64}} = []
	angleValues::Array{Float64} = []
	angleThetas::Array{Float64} = []
	charges::Array{Float64} = []
	vdwSigma1s::Array{Float64} = []
	vdwSigma2s::Array{Float64} = []
	vdwEpsilons::Array{Float64} = []

	open(filename, "r") do top
		contents = readlines(top)

		i = 1
		while i <= length(contents)
			if occursin("<masses>", contents[i])
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
					vdwSigma1 = parse(Float64, vdwContents[2])
					vdwSigma2 = parse(Float64, vdwContents[3])
					vdwEpsilon = parse(Float64, vdwContents[4])
					append!(vdwSigma1s, vdwSigma1)
					append!(vdwSigma2s, vdwSigma2)
					append!(vdwEpsilons, vdwEpsilon)
					i += 1
				end
			end
			i += 1
		end
	end

	vel::Array{Float64} = zeros((length(xyz), 3))

	topology = Dict("xyz" => xyz,
					"vel" => vel,
					"masses" => masses,
					"bondIdx" => bondIdx,
					"bondLengths" => bondLengths,
					"bondKs" => bondKs,
					"angleIdx" => angleIdx,
					"angleValues" => angleValues * (pi / 180),  # use ° in topology but radians internally
					"angleThetas" => angleThetas * (pi / 180),
					"charges" => charges,
					"vdwSigma1s" => vdwSigma1s,
					"vdwSigma2s" => vdwSigma2s,
					"vdwEpsilons" => vdwEpsilons,
					"nAtoms" => length(xyz))

	return topology
end


