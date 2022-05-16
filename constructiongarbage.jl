	# x1 = atom1[1]
	# y1 = atom1[2]
	# z1 = atom1[3]

	# x2 = atom2[1]
	# y2 = atom2[2]
	# z2 = atom2[3]

	# x3 = atom3[1]
	# y3 = atom3[2]
	# z3 = atom3[3]

	# x4 = atom4[1]
	# y4 = atom4[2]
	# z4 = atom4[3]
	# for i in 1:size(top.dihedralIdx, 1)


	# atom1 = xyz[top.dihedralIdx[i][1], :]
	# 	atom2 = xyz[top.dihedralIdx[i][2], :]
	# 	atom3 = xyz[top.dihedralIdx[i][3], :]
	# 	atom4 = xyz[top.dihedralIdx[i][4], :]
	# 	atom1Idx = top.dihedralIdx[i][1]
	# 	atom2Idx = top.dihedralIdx[i][2]
	# 	atom3Idx = top.dihedralIdx[i][3]
	# 	atom4Idx = top.dihedralIdx[i][4]
	# 	kTheta = top.dihedralThetas[i]

	# 	vec12 = atom1 .- atom2
	# 	vec32 = atom3 .- atom2
	# 	vec12 = pbcAdjust(top, vec12)
	# 	vec32 = pbcAdjust(top, vec32)
	# 	vec23 = -vec32
	# 	vec43 = atom4 .- atom3
	# 	vec43 = pbcAdjust(top, vec43)

	# 	normal123 = cross(vec12, vec32)
	# 	normal234 = cross(vec43, vec32)

	# 	# ?
	# 	# normal123 /= norm(normal123)
	# 	# normal234 /= norm(normal234)

	# 	dotprod = dot(normal123 / norm(normal123), normal234 / norm(normal234))
	# 	if dotprod > 1.0
	# 		dotprod = 1.0
	# 	elseif dotprod < -1.0
	# 		dotprod = -1.0
	# 	end
	# 	dihedralAngle = acos(dotprod)
	# 	# println(dihedralAngle * (180 / pi))
	# 	dihedralAngle0 = top.dihedralValues[i]

	# 	force1 = x
	# 	dihedralForces[atom1Idx, :] += force1
	# 	dihedralForces[atom2Idx, :] += force2
	# 	dihedralForces[atom3Idx, :] += force3
	# 	dihedralForces[atom4Idx, :] += force4
	# end
	# for i in 1:size(top.dihedralIdx, 1)
	# 	atom1 = xyz[top.dihedralIdx[i][1], :]
	# 	atom2 = xyz[top.dihedralIdx[i][2], :]
	# 	atom3 = xyz[top.dihedralIdx[i][3], :]
	# 	atom4 = xyz[top.dihedralIdx[i][4], :]
	# 	atom1Idx = top.dihedralIdx[i][1]
	# 	atom2Idx = top.dihedralIdx[i][2]
	# 	atom3Idx = top.dihedralIdx[i][3]
	# 	atom4Idx = top.dihedralIdx[i][4]
	# 	kTheta = top.dihedralThetas[i]

	# 	vec12 = atom1 .- atom2
	# 	vec32 = atom3 .- atom2
	# 	vec12 = pbcAdjust(top, vec12)
	# 	vec32 = pbcAdjust(top, vec32)
	# 	vec23 = -vec32
	# 	vec43 = atom4 .- atom3
	# 	vec43 = pbcAdjust(top, vec43)

	# 	normal123 = cross(vec12, vec32)
	# 	normal234 = cross(vec43, vec32)

	# 	# ?
	# 	# normal123 /= norm(normal123)
	# 	# normal234 /= norm(normal234)

	# 	dotprod = dot(normal123 / norm(normal123), normal234 / norm(normal234))
	# 	if dotprod > 1.0
	# 		dotprod = 1.0
	# 	elseif dotprod < -1.0
	# 		dotprod = -1.0
	# 	end
	# 	dihedralAngle = acos(dotprod)
	# 	# println(dihedralAngle * (180 / pi))
	# 	dihedralAngle0 = top.dihedralValues[i]

	# 	x1 = normal123[1]
	# 	y1 = normal123[2]
	# 	z1 = normal123[3]
	# 	x2 = 0.0
	# 	y2 = 0.0
	# 	z2 = 0.0
	# 	x3 = normal234[1]
	# 	y3 = normal234[2]
	# 	z3 = normal234[3]
	# 	r1::Array{Float64} = [x1, y1, z1] .- [x2, y2, z2]
	# 	r2::Array{Float64} = [x3, y3, z3] .- [x2, y2, z2]

	# 	r1 = pbcAdjust(top, r1)
	# 	r2 = pbcAdjust(top, r2)

	# 	L1::Float64 = norm(r1)
	# 	L2::Float64 = norm(r2)

	# 	limitAngle::Float64 = 0.001

	# 	if abs(dihedralAngle) < limitAngle || abs(pi - dihedralAngle) < limitAngle
	# 		# println("yes")
	# 		if dihedralAngle < limitAngle
	# 			signInFront = 1  # if close to 0
	# 		else
	# 			signInFront = -1  # if close to pi
	# 		end

	# 		if dot(r2 - r1, [1.0, 0.0, 0.0]) > 0.0
	# 			signX = 1
	# 		else
	# 			signX = -1
	# 		end

	# 		if dot(r2 - r1, [0.0, 1.0, 0.0]) > 0.0
	# 			signY = 1
	# 		else
	# 			signY = -1
	# 		end

	# 		if dot(r2 - r1, [0.0, 0.0, 1.0]) > 0.0
	# 			signZ = 1
	# 		else
	# 			signZ = -1
	# 		end

	# 		dA0dN1x = signInFront * signX * (1 / L2) * sin(acos(r2[1] / L2))
	# 		dA0dN1y = signInFront * signY * (1 / L2) * sin(acos(r2[2] / L2))
	# 		dA0dN1z = signInFront * signZ * (1 / L2) * sin(acos(r2[3] / L2))
	# 		dA0dN2x = signInFront * -signX * (1 / L1) * sin(acos(r1[1] / L1))
	# 		dA0dN2y = signInFront * -signY * (1 / L1) * sin(acos(r1[2] / L1))
	# 		dA0dN2z = signInFront * -signZ * (1 / L1) * sin(acos(r1[3] / L1))

	# 		# dA0dN1x = 0.0
	# 		# dA0dN1y = 0.0
	# 		# dA0dN1z = 0.0
	# 		# dA0dN2x = 0.0
	# 		# dA0dN2y = 0.0
	# 		# dA0dN2z = 0.0
	# 	else
	# 		p::Float64 = dot(r1, r2)
	# 		d::Float64 = p / (L1 * L2)
	# 		LL::Float64 = L1^2 * L2^2
	# 		prefactor = (-(1 / sqrt(1 - d^2))) * (1 / LL)

	# 		X1 = r1[1]
	# 		Y1 = r1[2]
	# 		Z1 = r1[3]
	# 		X2 = 0.0
	# 		Y2 = 0.0
	# 		Z2 = 0.0
	# 		X3 = r2[1]
	# 		Y3 = r2[2]
	# 		Z3 = r2[3]
	# 		# dA0dN1x = prefactor * (X2 * L1 * L2 - p * (L2 / L1) * X1)
	# 		# dA0dN2x = prefactor * (X1 * L1 * L2 - p * (L1 / L1) * X2)

	# 		# dA0dN1y = prefactor * (Y2 * L1 * L2 - p * (L2 / L1) * Y1)
	# 		# dA0dN2y = prefactor * (X1 * L1 * L2 - p * (L1 / L2) * Y2)

	# 		# dA0dN1z = prefactor * (Z2 * L1 * L2 - p * (L2 / L1) * Z1)
	# 		# dA0dN2z = prefactor * (Z1 * L1 * L2 - p * (L1 / L2) * Z2)

	# 		dA0dN1x = prefactor * (pbcAdjust(top, X3 - X2, 1) * L1 * L2 - p * ((L2 / L1) * pbcAdjust(top, X1 - X2, 1)))
	# 		dA0dN1y = prefactor * (pbcAdjust(top, Y3 - Y2, 2) * L1 * L2 - p * ((L2 / L1) * pbcAdjust(top, Y1 - Y2, 2)))
	# 		dA0dN1z = prefactor * (pbcAdjust(top, Z3 - Z2, 3) * L1 * L2 - p * ((L2 / L1) * pbcAdjust(top, Z1 - Z2, 3)))

	# 		dA0dN2x = prefactor * (pbcAdjust(top, X1 - X2, 1) * L1 * L2 - p * ((L1 / L2) * pbcAdjust(top, X3 - X2, 1)))
	# 		dA0dN2y = prefactor * (pbcAdjust(top, Y1 - Y2, 2) * L1 * L2 - p * ((L1 / L2) * pbcAdjust(top, Y3 - Y2, 2)))
	# 		dA0dN2z = prefactor * (pbcAdjust(top, Z1 - Z2, 3) * L1 * L2 - p * ((L1 / L2) * pbcAdjust(top, Z3 - Z2, 3)))
	# 	end

	# 	if -dot(cross(vec12, vec32), vec43) > 0.0
	# 		sigma = 1
	# 	else
	# 		sigma = -1
	# 	end

	# 	x1 = atom1[1]
	# 	y1 = atom1[2]
	# 	z1 = atom1[3]

	# 	x2 = atom2[1]
	# 	y2 = atom2[2]
	# 	z2 = atom2[3]

	# 	x3 = atom3[1]
	# 	y3 = atom3[2]
	# 	z3 = atom3[3]

	# 	x4 = atom4[1]
	# 	y4 = atom4[2]
	# 	z4 = atom4[3]

	# 	forceX1 = sigma * (dA0dN1y * pbcAdjust(top, z2 - z3, 3) + dA0dN1z * pbcAdjust(top, y3 - y2, 2))
	# 	forceY1 = sigma * (dA0dN1x * pbcAdjust(top, z3 - z2, 3) + dA0dN1z * pbcAdjust(top, x2 - x3, 1))
	# 	forceZ1 = sigma * (dA0dN1x * pbcAdjust(top, y2 - y3, 2) + dA0dN1y * pbcAdjust(top, x3 - x2, 1))

	# 	forceX2 = sigma * (dA0dN1y * pbcAdjust(top, z3 - z1, 3) + dA0dN1z * pbcAdjust(top, y1 - y3, 2) + dA0dN2y * pbcAdjust(top, z3 - z4, 3) + dA0dN2z * pbcAdjust(top, y4 - y3, 2))
	# 	forceY2 = sigma * (dA0dN1x * pbcAdjust(top, z1 - z3, 3) + dA0dN1z * pbcAdjust(top, x3 - x1, 1) + dA0dN2x * pbcAdjust(top, z4 - z3, 3) + dA0dN2z * pbcAdjust(top, x3 - x4, 1))
	# 	forceZ2 = sigma * (dA0dN1x * pbcAdjust(top, y3 - y1, 2) + dA0dN1y * pbcAdjust(top, x1 - x3, 1) + dA0dN2x * pbcAdjust(top, y3 - y4, 2) + dA0dN2y * pbcAdjust(top, x4 - x3, 1))
	# 	forceX3 = sigma * (dA0dN1y * pbcAdjust(top, z1 - z2, 3) + dA0dN1z * pbcAdjust(top, y2 - y1, 2) + dA0dN2y * pbcAdjust(top, z4 - z2, 3) + dA0dN2z * pbcAdjust(top, y2 - y4, 2))
	# 	forceY3 = sigma * (dA0dN1x * pbcAdjust(top, z2 - z1, 3) + dA0dN1z * pbcAdjust(top, x1 - x2, 1) + dA0dN2x * pbcAdjust(top, z2 - z4, 3) + dA0dN2z * pbcAdjust(top, x4 - x2, 1))
	# 	forceZ3 = sigma * (dA0dN1x * pbcAdjust(top, y1 - y2, 2) + dA0dN1y * pbcAdjust(top, x2 - x1, 1) + dA0dN2x * pbcAdjust(top, y4 - y2, 2) + dA0dN2y * pbcAdjust(top, x2 - x4, 1))

	# 	forceX4 = sigma * (dA0dN2y * pbcAdjust(top, z2 - z3, 3) + dA0dN2z * pbcAdjust(top, y3 - y2, 2))
	# 	forceY4 = sigma * (dA0dN2x * pbcAdjust(top, z3 - z2, 3) + dA0dN2z * pbcAdjust(top, x2 - x3, 1))
	# 	forceZ4 = sigma * (dA0dN2x * pbcAdjust(top, y2 - y3, 2) + dA0dN2y * pbcAdjust(top, x3 - x2, 1))
		
	# 	force1 = -kTheta * (dihedralAngle - dihedralAngle0) .* [forceX1, forceY1, forceZ1]
	# 	force2 = -kTheta * (dihedralAngle - dihedralAngle0) .* [forceX2, forceY2, forceZ2]
	# 	force3 = -kTheta * (dihedralAngle - dihedralAngle0) .* [forceX3, forceY3, forceZ3]
	# 	force4 = -kTheta * (dihedralAngle - dihedralAngle0) .* [forceX4, forceY4, forceZ4]

	# 	dihedralForces[atom1Idx, :] += force1
	# 	dihedralForces[atom2Idx, :] += force2
	# 	dihedralForces[atom3Idx, :] += force3
	# 	dihedralForces[atom4Idx, :] += force4

	# end

	# return dihedralForces









			# if angle < limitAngle
			# 	signInFront = 1  # if close to 0
			# else
			# 	signInFront = -1  # if close to pi
			# end

			# vane::Array{Float64} = r2 - r1
			# if dot(vane, [1.0, 0.0, 0.0]) > 0.0
			# 	signX::Float64 = 1
			# else
			# 	signX = -1
			# end

			# if dot(vane, [0.0, 1.0, 0.0]) > 0.0
			# 	signY::Float64 = 1
			# else
			# 	signY = -1
			# end

			# if dot(vane, [0.0, 0.0, 1.0]) > 0.0
			# 	signZ::Float64 = 1
			# else
			# 	signZ = -1
			# end

			# forceX1 = signInFront * signX * (1 / L2) * sin(acos(r2[1] / L2))
			# forceY1 = signInFront * signY * (1 / L2) * sin(acos(r2[2] / L2))
			# forceZ1 = signInFront * signZ * (1 / L2) * sin(acos(r2[3] / L2))

			# forceX3 = signInFront * -signX * (1 / L1) * sin(acos(r1[1] / L1))
			# forceY3 = signInFront * -signY * (1 / L1) * sin(acos(r1[2] / L1))
			# forceZ3 = signInFront * -signZ * (1 / L1) * sin(acos(r1[3] / L1))


			# x1::Float64 = atom1[1]
			# y1::Float64 = atom1[2]
			# z1::Float64 = atom1[3]
			# x2::Float64 = atom2[1]
			# y2::Float64 = atom2[2]
			# z2::Float64 = atom2[3]
			# x3::Float64 = atom3[1]
			# y3::Float64 = atom3[2]
			# z3::Float64 = atom3[3]
			# r1::Array{Float64}  = atom1 .- atom2
			# r2::Array{Float64} = atom3 .- atom2

			# r1 = pbcAdjust(top, r1)
			# r2 = pbcAdjust(top, r2)

			# L1::Float64 = norm(r1)
			# L2::Float64 = norm(r2)

			# p = dot(r1, r2)
			# d = p / (L1 * L2)
			# LL = L1^2 * L2^2
			# prefactor = (-(1 / sqrt(1 - d^2))) * (1 / LL)

			# forceX1 = prefactor * (pbcAdjust(top, x3 - x2, 1) * L1 * L2 - p * ((L2 / L1) * pbcAdjust(top, x1 - x2, 1)))
			# forceY1 = prefactor * (pbcAdjust(top, y3 - y2, 2) * L1 * L2 - p * ((L2 / L1) * pbcAdjust(top, y1 - y2, 2)))
			# forceZ1 = prefactor * (pbcAdjust(top, z3 - z2, 3) * L1 * L2 - p * ((L2 / L1) * pbcAdjust(top, z1 - z2, 3)))

			# forceX3 = prefactor * (pbcAdjust(top, x1 - x2, 1) * L1 * L2 - p * ((L1 / L2) * pbcAdjust(top, x3 - x2, 1)))
			# forceY3 = prefactor * (pbcAdjust(top, y1 - y2, 2) * L1 * L2 - p * ((L1 / L2) * pbcAdjust(top, y3 - y2, 2)))
			# forceZ3 = prefactor * (pbcAdjust(top, z1 - z2, 3) * L1 * L2 - p * ((L1 / L2) * pbcAdjust(top, z3 - z2, 3)))


		# force1 = -kTheta * (angle - angle0) .* [forceX1, forceY1, forceZ3]
		# force3 = -kTheta * (angle - angle0) .* [forceX3, forceY3, forceZ3]
