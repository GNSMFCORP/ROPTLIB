
# JTestPhaseRetrieval.jl
# Apply Riemannian method for Phase Retrieval problem
#
# --- by Wen Huang

function JPhaseRetrieval_Exp2_Recovery()
	ls = [3]#[2, 3, 4, 5, 6, 7, 8, 9, 10]
	num = 100
	data = zeros(3, length(ls), num)
	for i in 1 : length(ls)
		for j in 1 : num
			(data[1, i, j], data[2, i, j], data[3, i, j]) = JPhaseRetrievalLRBFGSvsFlow(ls[i], j)
		end
	end
	println(pwd())
	writedlm("data.txt", data)
	for i in 1 : 3
		for j in 1 : length(ls)
			for k in 1 : num
				data[i, j, k] = (data[i, j, k] < 1e-5) ? 1 : 0
			end
		end
	end

	tab = zeros(3, length(ls))
	for i = 1 : length(ls)
		tab[1, i] = sum(data[1, i, :]) / num
		tab[2, i] = sum(data[2, i, :]) / num
		tab[3, i] = sum(data[3, i, :]) / num
	end

	open("./Paper_LRBFGS_WF_Recovery.txt", "w") do f
		for i in 1 : size(data, 1)
			for j in 1 : size(data, 2)
				datastr = outputfloat(tab[i, j])
				if(j < size(data, 2))
					write(f, "\$" * datastr * "\$ & ")
				else
					write(f, "\$" * datastr * "\$ ")
				end
			end
			write(f, "\\\\\n")
		end
	end
end

function outputfloat(x)
	if(x <= 0)
		return "0"
	end
	
	p = log(x) / log(10)
	p = -ceil(-p)
	x = round(x * 10^(-p) * 100)
	x = x / 100
	# Formatting need be installed first by command Pkg.add("Formatting")
	# using Formatting
	strx = sprintf1("%3.2f", x)
	strp = sprintf1("%d", p)
	if(p != 0)
        str = strx * "_{" * strp * "}"
    else
        str = strx
    end
	return str
end

function JPhaseRetrievalLRBFGSvsFlow(ll, randseed)
	global n1 = 128
	global n2 = 1
	global n = n1 * n2
	global l = ll
	global m = l * n
	global p = 1
	global delta = 0.95
	global kappa = 1 / n
	global tol = 1e-10

	# problem related parameters
	global seed = floor(rand() * 1000000)#randseed
	println("seed:", seed)
	srand(Int64(seed))
	global masks = complex(randn(n1, n2, l), randn(n1, n2, l))
#	global masks = octanaryrandom(n1, n2, l)
	global phase = complex(randn(n1, n2), randn(n1, n2))
	phase = phase / norm(phase)
	global vphase = reshape(phase, n, 1)
	global b = zeros(m, 1)
	for i in 1 : l
		zi = reshape(fft(phase .* masks[:, :, i]) / sqrt(n), n, 1)
		b[(i - 1) * n + 1 : i * n] = conj(zi) .* zi;
	end
	# initial iterate for Riemannian approach.
	initialX = InitialIterate(p, 50)
	
	# initial step size at the first iterate
	(Xopt, nfftrank1) = JTestPhaseRetrieval(initialX, tol, 50 * 2 * l * p)
	R1 = vecnorm(vphase - Xopt * (Xopt'*vphase) / norm(Xopt'*vphase)) / vecnorm(phase)
	println("p=1: dist(sol, xinit) / |sol|:", R1)

	p = 2
	# initial iterate for Riemannian approach.
	initialX = InitialIterate(p, 50)#Int64(floor(50 / p)))

	# initial step size at the first iterate
	(Xopt, nfftrank2) = JTestPhaseRetrieval(initialX, tol, 50 * 2 * l * p)
	R2 = vecnorm(vphase - Xopt * (Xopt'*vphase) / norm(Xopt'*vphase)) / vecnorm(phase)
	println("p=2: dist(sol, xinit) / |sol|:", R2)

	p = 4
	# initial iterate for Riemannian approach.
	initialX = InitialIterate(p, 50)#Int64(floor(50 / p)))

	# initial step size at the first iterate
	(Xopt, nfftrank2) = JTestPhaseRetrieval(initialX, tol, 50 * 2 * l * p)
	R4 = vecnorm(vphase - Xopt * (Xopt'*vphase) / norm(Xopt'*vphase)) / vecnorm(phase)
	println("p=4: dist(sol, xinit) / |sol|:", R4)

return (R1, R2, R4)
	# initial iterate for Wirtinger Flow.
	p = 1
	initialX = InitialIterate(p, Int64(floor(50 / p)))
	tic()
	(Xopt, nfftWF) = WFlow(initialX, tol, Int64(floor(50 / p)) * 2 * l)
	toc()
	WF = vecnorm(vphase - Xopt * (Xopt'*vphase) / norm(Xopt'*vphase)) / vecnorm(phase)
	println("WF: dist(sol, xinit) / |sol|:", WF)
	return (WF, R1, R2, R4)
end

function octanaryrandom(n1, n2, l)
	masks = complex(zeros(n1, n2, l), zeros(n1, n2, l))
	for i in 1 : n1
		for j in 1 : n2
			for k in 1 : l
				tmp = rand()
				if(tmp < 0.25)
					b1 = 1
				elseif(tmp < 0.5)
					b1 = -1
				elseif(tmp < 0.75)
					b1 = -im
				else
					b1 = im
				end
				tmp = rand()
				if(tmp < 0.8)
					b2 = sqrt(2.0)/2
				else
					b2 = sqrt(3.0)
				end
				masks[i, j, k] = b1 * b2
			end
		end
	end
	return masks
end

function InitialIterate(k, niter)

	srand(Int64(seed))
	tic()
	lambda = sqrt(sum(b) / real(sum(conj(masks) .* masks)) * n)

	initialX = complex(randn(n, k), randn(n, k))
	initialX = initialX / norm(initialX)
	for i in 1 : niter
		initialX = Ax(initialX, k)
	end
	if(niter >= 1)
		initialX = initialX / norm(initialX)
	end
#	initialX = qr(complex(randn(n, k), randn(n, k)))[1]
#	for i in 1 : niter
#		initialX = qr(Ax(initialX, k))[1]
#	end
	toc()
	return lambda * initialX
end

function WFlow(initialX, tol, initnfft)
	gradf = egf(initialX)
	ngf = norm(gradf)
	ngf0 = ngf
	t = 0
	deno = norm(initialX)^2
	X = initialX
	while(ngf > tol && t < 2500)
		# this stepsize works well for n = 2^2 to 256^2
		if(l == 3)
			stepsize = min(1.0 - exp(-(1.0 + t) / 330.0), 0.2) * (20) / deno
		elseif(l == 4)
			stepsize = min(1.0 - exp(-(1.0 + t) / 330.0), 0.2) * (1) / deno
		end
		gradf = egf(X)
		if(mod(t, 1) == 0)
			ngf = norm(gradf)
			println("iter", t, ",ngf:", ngf, ", ngf/ngf0:", ngf / ngf0, ", stepsize:", stepsize)
		end
		X -= stepsize * gradf
		t+=1
	end
	ngf = norm(gradf)
	println("final iter:", t, ",ngf:", ngf, ", ngf/ngf0:", ngf / ngf0, ", nfft:", initnfft + (t+1)*2*l)
	return (X, initnfft + (t+1)*2*l)
end

function Ax(x, k)
	sqrtn = sqrt(n)
	ZY = complex(zeros(m, k), zeros(m, k))
	for i in 1 : l
		for j in 1 : k
			temp = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
			ZY[(i - 1) * n + 1 : i * n, j] = temp
		end
	end
	DZY = complex(zeros(m, k), zeros(m, k))
	for i in 1 : k
		DZY[:, i] = b .* ZY[:, i]
	end

	Y = complex(zeros(n, k), zeros(n, k))
	for i in 1 : l
		temp = complex(zeros(n, k), zeros(n, k))
		for j in 1 : k
			temp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
			temp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* temp[:, j]
		end
		Y += temp
	end
	return Y
end

#--------------------main test function for Riemannian method-----------------------------------------
function JTestPhaseRetrieval(initialX, tol, initnfft)
	# set domain manifold to be the complex noncompact Stiefel manifold quotient the unitary group
	mani1 = "CpxNStQOrth"; ManArr = [pointer(mani1)]
	UseDefaultArr = [-1] # -1 means that the default value in C++ is used.
	numofmani = [1]
	# The size is n by p
	ns = [n]
	ps = [p]
	Mparams = ManiParams(0, length(ManArr), pointer(ManArr), pointer(numofmani), pointer(UseDefaultArr), pointer(UseDefaultArr), pointer(ns), pointer(ps))

	# set function handles
	fname = "func"
	gfname = "gfunc"
	hfname = "hfunc"
	isstopped = "stopfunc"
	LinesearchInput = "" # or "LSfunc" to use the function defined below
	Handles = FunHandles(pointer(fname), pointer(gfname), pointer(hfname), pointer(isstopped), pointer(LinesearchInput))

	# set solvers by modifying the default one.
	method = "LRBFGS"
	Sparams.name = pointer(method)
	Sparams.OutputGap = 10
	Sparams.LengthSY = 2
	Sparams.Min_Iteration = 10
	Sparams.Max_Iteration = 2500
	Sparams.DEBUG = 1
	Sparams.Tolerance = tol
	Sparams.Stop_Criterion = 2
	Sparams.Accuracy = 1e-6
	Sparams.Finalstepsize = 1
	Sparams.IsCheckGradHess = 0
	Sparams.IsCheckParams = 0
	Sparams.isconvex = 1
	#Sparams.Initstepsize = stepsize0

	# use locking condition or not
	HasHHR = 0::Int64
	FinalIterate = []; FF = []; GG = []; TT = []; time = 0; nf = 0; ng = 0; nR = 0; nH = 0; nV = 0; nVp = 0; nfft = initnfft; nn = 0;
	while true
		global p
		println("rank is ", p)
		ps = [p]
		Mparams.p = pointer(ps)
		# Call the solver and get results. See the user manual for details about the outputs.
		(FinalIterate, fv, gfv, gfgf0, iter, nfi, ngi, nRi, nVi, nVpi, nHi, ComTime, funs, grads, times) = DriverJuliaOPT(Handles, Sparams, Mparams, HasHHR, initialX)

		if(length(TT) > 0)
			FF = [FF; funs]
			GG = [GG; grads]
			TT = [TT; times + TT[end]]
		else
			FF = funs
			GG = grads
			TT = times
		end
		time += ComTime; nf += nfi; ng += ngi; nR += nRi; nH += nHi; nV += nVi; nVp += nVpi; 
		nfft += (nfi + ngi) * l * p; nn += 4 * ngi * p^2 + 4 * (length(FF) - 1) * p^2 + nVi * 5 * p^2;
		if(p == 1)
			break;
		end
		if(isnan(sum(FinalIterate)) || isinf(sum(FinalIterate)))
			return (zeros(n, 1), 0)
		end
		(U, S) = svd(FinalIterate)
		normSdsqrk = norm(S) / sqrt(p)
		for i = 1 : p
			if(S[i] / normSdsqrk < delta)
				initialX = U[:, 1 : i - 1] * diagm(S[1 : i - 1])
				break;
			end
		end
		p = size(initialX, 2)
		if(p == 1)
			println("dist(sol, xinit by rank reduce) / |sol|:", vecnorm(vphase - initialX * (initialX'*vphase) / norm(initialX'*vphase)) / vecnorm(vphase))
		end
	end
	println("Num of total iter.:", length(FF), ", total time:", time, ", nf:", nf, ", ng:", ng, ", nH:", nH, ", nV:", nV, ", nR:", nR,
			", nfft:", nfft)
	return (FinalIterate, nfft)
end


# Define function handles
# The function names are assigned to the "FunHandles" struct.
# See lines 17-21
function func(xreal, inTmpreal)
	x = real2complex(xreal)
	x = reshape(x, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	sqrtn = sqrt(n)
	outTmp = complex(zeros(m, p + 1), zeros(m, p + 1)) # the first k columns are for ZY and the last column is for D.
	for i in 1 : l
		sZqi = zeros(n, 1)
		for j in 1 : p
			temp = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
			outTmp[(i - 1) * n + 1 : i * n, j] = temp
			sZqi += conj(temp) .* temp
		end
		outTmp[(i - 1) * n + 1 : i * n, p + 1] = sZqi - b[(i - 1) * n + 1 : i * n]
	end
	output = norm(outTmp[:, p + 1])^2 / norm(b)^2
	if(p != 1)
		output += kappa * vecnorm(x)^2
	end
	return (output, complex2real(outTmp)) # The temparary data "outTmp" will replace the "inTmp"
end

function gfunc(xreal, inTmpreal) # The inTmp is the temparary data computed in "func".
	x = real2complex(xreal)
	inTmp = real2complex(inTmpreal)
	x = reshape(x, n, p) # All the input argument is a vector. One has to reshape it to have a proper size
	inTmp = reshape(inTmp, m, p + 1) # All the input argument is a vector. One has to reshape it to have a proper size
	sqrtn = sqrt(n)
	DZY = complex(zeros(m, p), zeros(m, p))
	for i = 1 : p
		DZY[:, i] = inTmp[:, end] .* inTmp[:, i]
	end

	Y = complex(zeros(n, p), zeros(n, p))
	for i in 1 : l
		temp = complex(zeros(n, p), zeros(n, p))
		for j in 1 : p
			temp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
			temp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* temp[:, j]
		end
		Y += temp
	end
	Y = 4.0 * Y / norm(b)^2
	if(p != 1)
		Y += kappa * 2 * x
	end
	gf = Y / (x' * x)

	return (complex2real(gf), []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

# action of the Hessian is not used. So use identity to avoid error.
function hfunc(xreal, inTmpreal, etareal)
	return (etareal, []) # If one does not want to change the temparary data, then let the outTmp be an empty array.
end

# Users can define their own stopping criterion by passing the name 
# of the function to the "isstopped" field in the object of structure FunHandles
function stopfunc(xreal, gfreal, fx, ngfx, ngfx0)
# x: the current iterate
# gf: the gradient at x
# gx: the function value at x
# ngfx: the norm of gradient at x
# ngfx0: the norm of gradient at the initial iterate
	if(isnan(sum(xreal)) || isinf(sum(xreal)))
		return true
	end
	x = real2complex(xreal)
	x = reshape(x, n, p)
	if(p > 1)
		S = svdvals(x)
		normSdsqrk = norm(S) / sqrt(p)
		if(isnan(normSdsqrk) || isinf(normSdsqrk))
			return true
		end
		for i = 1 : p
			if(S[i] / normSdsqrk < delta)
				return true
			end
		end
	end
	return (ngfx < tol)
end


function egf(x) # the Euclidean gradient for Wirtinger flow method
	sqrtn = sqrt(n)
	ZY = complex(zeros(m, p), zeros(m, p))
	D = zeros(m, 1)
	for i in 1 : l
		sZqi = zeros(n, 1)
		for j in 1 : p
			temp = reshape(fft(reshape(x[:, j], n1, n2) / sqrtn .* masks[:, :, i]), n, 1)
			ZY[(i - 1) * n + 1 : i * n, j] = temp
			conjtemp = conj(temp);
			sZqi += conj(temp) .* temp
		end
		D[(i - 1) * n + 1 : i * n] = sZqi - b[(i - 1) * n + 1 : i * n]
	end

	DZY = complex(zeros(m, p), zeros(m, p))
	for i = 1 : p
		DZY[:, i] = D .* ZY[:, i]
	end
	Y = complex(zeros(n, p), zeros(n, p))
	for i in 1 : l
		temp = complex(zeros(n, p), zeros(n, p))
		for j in 1 : p
			temp[:, j] = reshape(ifft(reshape(DZY[(i - 1) * n + 1 : i * n, j] * sqrtn, n1, n2)), n, 1)
			temp[:, j] = conj(reshape(masks[:, :, i], n, 1)) .* temp[:, j]
		end
		Y += temp
	end
	Y = 4.0 * Y / norm(b)^2
	return Y
end

