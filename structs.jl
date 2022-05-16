struct Topology
	masses::Array{Float64}
	bondIdx::Array{Int64}
	bondKs::Array{Float64}
	bondLengths::Array{Float64}
	angleIdx::Array{Int64}
	angleValues::Array{Float64}
	angleThetas::Array{Float64}
	dihedralIdx::Array{Int64}
	dihedralValues::Array{Float64}
	dihedralThetas::Array{Float64}
	charges::Array{Float64}
	vdwIdx::Array{Int64}
	vdwSigmas::Array{Float64}
	vdwEpsilons::Array{Float64}
	nAtoms::Int64
	atomTypes::Array{String}
	box::Array{Float64}
end

mutable struct Energies
	bond::Float64
	angle::Float64
	dihedral::Float64
	coulomb::Float64
	vdw::Float64
end