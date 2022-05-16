function complexExp(complexNum::ComplexF64)
	# faster than the built in exp(im*)
	return cos(complexNum) + im * sin(complexNum)
end

function pbcAdjust(top::Topology, crd)
	for k in 1:3
		if crd[k] > top.box[k] / 2
			crd[k] -= top.box[k]
		elseif crd[k] <= -top.box[k] / 2
			crd[k] += top.box[k]
		end
	end

	return crd
end

function pbcAdjust(top::Topology, crdk::Float64, k::Int64)
	if crdk > top.box[k] / 2
		crdk -= top.box[k]
	elseif crdk <= -top.box[k] / 2
		crdk += top.box[k]
	end

	return crdk
end

function getCombinations(start::Int64, finish::Int64, nDim::Int64)::Array{Int64}
	#= produces an integer array like [[0, 0], [0, 1], [1, 0], [1, 1]] etc.
	   start > finish and start < finish are acceptable =#

	n = abs(finish - start) + 1
	if finish - start < 0
		sign = -1
	else
		sign = 1
	end	

	N = n^nDim
	comb = zeros(Int64, (N, nDim))
	@inbounds for i in 1:N
		for j in 1:nDim
			comb[i, j] = start + sign * ((i - 1) ÷ (N / n^j)) % n
		end
	end


	return comb
end

function testgetCombinations()
	@time comb::Array{Int64} = getCombinations(2, -1, 3)

	for i::Int64 in 1:size(comb)[1]
		println(comb[i, :])
	end

	@time comb = getCombinations(-1, 2, 3)

	for i::Int64 in 1:size(comb)[1]
		println(comb[i, :])
	end
end
# testgetCombinations()

function getBondedNeighbors(bondIdx::Array{Array{Int64}}, nAtoms::Int64)::Array{Array{Int64}}
	network::Array{Array{Int64}} = [Int64[] for i in 1:nAtoms] 
	for i::Int64 in 1:nAtoms
		for b in 1:size(bondIdx)[1]
			if bondIdx[b][1] == i
				if !(bondIdx[b][2] in network[i])
					push!(network[i], bondIdx[b][2])
				end
			elseif bondIdx[b][2] == i
				if !(bondIdx[b][1] in network[i])
					push!(network[i], bondIdx[b][1])
				end
			end
		end
	end

	return network

end
