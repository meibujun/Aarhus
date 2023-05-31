#Calculate MAF
function _counts(s)
    cc = zeros(Int, (4, size(s, 2)))
    m, n = size(s)
    @inbounds for j in 1:n
                for i in 1:m                    
                    cc[s[i, j] + 1, j] += 1
                end
            end
    return cc
end
function maf(s)
    cc = _counts(s)
    col = size(s, 2)
    out = Array{Float64, 1}(undef, col)
    @inbounds for j in 1:col
        out[j] = (cc[3, j] + 2cc[4, j]) / 2(cc[1, j] + cc[3, j] + cc[4, j])
        (out[j] > 0.5) && (out[j] = 1 - out[j])
    end
    out
end
#Function to calculate allele frequencies
function getaf(genmatrix, snpdim=2)
  af = zeros(Float64, size(genmatrix)[snpdim])
  for c in 1:size(genmatrix)[snpdim]
    if snpdim == 2
      af[c] = mean(genmatrix[:,c])/2
    else
      af[c] = mean(genmatrix[c,:])/2
    end
  end
  return af
end
#Compute genotype frequencies
function gf(af1,af2)
	nqtl = size(af1)[1]
	gfreq = zeros(Float64, nqtl, 3)
	
	for i in 1:nqtl
		pAA = af1[i] * af2[i]
		pAa = (af1[i] * (1-af2[i])) + ((1-af1[i]) * af2[i])
		paa = (1-af1[i]) * (1-af2[i])
	
		gfreq[i,1] = pAA
		gfreq[i,2] = pAa
		gfreq[i,3] = paa
	end
	return gfreq
end
#Function to get statistical avarage, dominance, and epistatic effects from functional effects
function epi(gfreq, af, epicoef, a, epimodel, epicomb)
	nqtl = size(af)[1]
	ncombs = size(epicomb)[1]	
	epi_effects = zeros(Float64, ncombs, 9)
	epi_vars = zeros(Float64, ncombs, 9)
	if epimodel == "complementary"
		epimat = [1 1 0; 1 1 0; 0 0 0]
	elseif  epimodel == "aa"
		epimat = [2 1 0; 1 1 1; 0 1 2]
	elseif epimodel == "dd"
		epimat = [0 1 0; 1 1 1; 0 1 0]
	elseif epimodel == "interaction"
		epimat = [2 1 0]' * [2 1 0]
	else
		epimat = [0 0 0; 0 0 0; 0 0 0]
	end	
	#scale and center the epimat
	if epimodel !== "none"
		epimat = (epimat .- mean(epimat)) ./ std(epimat)
	end	
	for c in 1:ncombs		
		#get locus numbers
		i = epicomb[c,1]
		j = epicomb[c,2]		
		#allele frequencies
		af_A = af[i]
		af_B = af[j]	
		#create y, the effect of epistatic coefficients on the phenotype
		y2 = epimat .* (epicoef[c] * sqrt(abs(a[i])) * sqrt(abs(a[j])))
		y = vec(y2')		
		#Get genototype frequencies to create D and W
		pAA = gfreq[i,1]
		pAa = gfreq[i,2]
		paa = gfreq[i,3]
		pBB = gfreq[j,1]
		pBb = gfreq[j,2]
		pbb = gfreq[j,3]
		#D matrix with genotype probabilities
		DA = [pAA 0 0; 0 pAa 0 ; 0 0 paa]
		DB = [pBB 0 0; 0 pBb 0 ; 0 0 pbb]
		D = kron(DB, DA)
		if det(D) == 0
			epi_vars[c,:] = [af_A af_B 0 0 0 0 0 0 0]
			continue
		end  
		#create W, following NOIA parameterization
		WA = zeros(Float64,3,3)
		WB = zeros(Float64,3,3)
		WA[:,1] .= 1
		sc_a = pAA+paa-((pAA-paa)^2)
		sc_b = pBB+pbb-((pBB-pbb)^2)			
		WA[1,2] = -(-pAa - (2*paa))
		WA[2,2] = -(1-pAa - (2*paa))
		WA[3,2] = -(2-pAa - (2*paa))
		if sc_a == 0.0   #locus is fixed
			WA[:,3] .= 0.0
		else
			WA[1,3] = -2*pAa*paa/sc_a
			WA[2,3] = 4*pAA*paa/sc_a
			WA[3,3] = -2*pAA*pAa/sc_a
		end
		WB[:,1] .= 1
		WB[1,2] = -(-pBb - (2*pbb))
		WB[2,2] = -(1-pBb - (2*pbb))
		WB[3,2] = -(2-pBb - (2*pbb))

		if sc_b == 0.0  #locus is fixed
			WB[:,3] .= 0.0
		else
			WB[1,3] = -2*pBb*pbb/sc_b
			WB[2,3] = 4*pBB*pbb/sc_b
			WB[3,3] = -2*pBB*pBb/sc_b
		end
		W = kron(WB, WA)		
		b = inv(W' * D * W) * W' * D * y
		B = Diagonal(b)
		V = B' * W' * D * W * B
		VA = V[2,2]+V[4,4]
		VD = V[3,3]+V[7,7]
		VAA = V[5,5]
		VAD = V[6,6]+V[8,8]
		VDD = V[9,9]
		epi_vars[c,:] = [af_A af_B VA VD VAA VAD VDD VA+VD+VAA+VAD+VDD epicoef[c]]		
		epi_effects[c,:] = b'
	end
	return epi_vars, epi_effects
end
#Compute alphas and all other effects following NOIA
function get_effects(a,d,p,epi_comb,epi_effects)
	println("Computing orthogonal alphas and ds")
	nqtl = size(a)[1]
	alphas = zeros(Float64, nqtl)
	ds = zeros(Float64, nqtl)	
	for i in 1:nqtl
		a1 = a[i] + (1 - (2*p[i])) *d[i]		
		flt1 = epi_comb[:,1] .== i
		ea1 = epi_effects[flt1,2]
		ed1 = epi_effects[flt1,3]		
		flt2 = epi_comb[:,2] .== i
		ea2 = epi_effects[flt2,4]
		ed2 = epi_effects[flt2,7]		
		alphas[i] = a1 + sum(ea1) + sum(ea2)
		ds[i] = d[i] + sum(ed1) + sum(ed2)		
	end
	return alphas, ds
end
#Make NOIA parmterization
function make_noia(geno)	
	println("Creating NOIA matrices..")	
	noia = convert(Array{Float32,2}, geno)
	noia_dom = copy(noia)	
	n = size(noia)[1]
	for i in 1:size(noia)[2]
		pAA = size(noia[noia[:,i] .== 2, :])[1] / n
		pAa = size(noia[noia[:,i] .== 1, :])[1] / n
		paa = size(noia[noia[:,i] .== 0, :])[1] / n
		sc = pAA+paa-((pAA-paa)^2)		
		if sc == 0.0
			noia_dom[:,i] .= 0.0
		else
			noia_dom[geno[:,i] .== 2,i] .= -(2*pAa*paa)/sc
			noia_dom[geno[:,i] .== 1,i] .= 4*pAA*paa/sc
			noia_dom[geno[:,i] .== 0,i] .= -(2*pAA*pAa)/sc
		end
		noia[geno[:,i] .== 2, i] .= -(-pAa - (2*paa))
		noia[geno[:,i] .== 1, i] .= -(1-pAa - (2*paa))
		noia[geno[:,i] .== 0, i] .= -(2-pAa - (2*paa))
	end	
	println("Done creating NOIA matrices!")	
	return noia, noia_dom
end
function Multiply_adjacent_columns(M)
    col = size(M,2)
    row = size(M,1)
    Mnew =  Array{Int64,2}(undef, row, Int(col/2))
    k = 1
    for i = 1: Int(col/2)
        Mnew[:,i] = M[:,k] .* M[:,k+1]
        k = k + 2
    end
    return Mnew
end
function epistaticLocMean(epiLocM)
    row = size(epiLocM,1)
    #row = size(M,1)
    Mnew =  Array{Float64,1}(undef, Int(row/2))
    k = 1
    for i = 1: Int(row/2)
        #Mnew[i] = Int.((M[k] .+ M[k+1])/2)
        Mnew[i] = median([epiLocM[k], epiLocM[k+1]])
        k = k + 2
    end
    return Mnew
end
function winTest(W,y,Va,Ve,wSize=100)
    W = W .- mean(W,dims=1)
    V = W*W'*Va + I*Ve    
    Vi = inv(V)
    n,k = size(W)
    X = [ones(n) zeros(n)]  
    nWin = ceil(Int64,k/wSize)
    testStat = zeros(nWin)
    wStartV  = Array{Int64,1}(undef,nWin)
    wEndV    = Array{Int64,1}(undef,nWin)
    for i=1:nWin
        wStart = (i-1)*wSize + 1
        wEnd   = wStart + wSize - 1
        wEnd   = wEnd <= k ? wEnd : k
        X = [ones(n) W[:,collect(wStart:wEnd)]]
        lhsi   = inv(X'Vi*X + I*0.001)
        rhs = X'Vi*y
        betaHat = lhsi*rhs
        vBetaH = lhsi - I*Va
        testStat[i] = betaHat[2:end]'inv(vBetaH[2:end,2:end])*betaHat[2:end]
        wStartV[i] = wStart
        wEndV[i] = wEnd
    end
    return wStartV,wEndV,testStat
end
function TruePosEpi_com(epi_combinations,epistaticQTLPos, QTLPos)
	QeInQtlPos = [1:size(epistaticQTLPos)[1] findall(x->x in epistaticQTLPos, QTLPos)]
	epi_combinationsR, epi_combinationsC = size(epi_combinations)
	QeInQtlPosR, QeInQtlPosC = size(QeInQtlPos)
	for i in 1:epi_combinationsR
		for j in 1:epi_combinationsC
			for k in 1:QeInQtlPosR            
				if epi_combinations[i,j] == QeInQtlPos[k,1]
					epi_combinations[i,j] = QeInQtlPos[k,2]
				end            
			end
		end
	end
	return epi_combinations
end
function TruePosDom(SortQTLDompos,QTLPosIndex)
	QTLPosIndexR, QTLPosIndexC = size(QTLPosIndex)
	sortDominantQTLPosR = size(SortQTLDompos)[1]
	for i in 1:sortDominantQTLPosR    
		for k in 1:QTLPosIndexR            
			if SortQTLDompos[i] == QTLPosIndex[k,1]
				SortQTLDompos[i] = QTLPosIndex[k,2]
			end            
		end    
	end
	return SortQTLDompos
end
function TruePosEpis(epistaticQTLPos,QTLPosIndex)
	sortepistaticQTLPos = sort(epistaticQTLPos)
	QTLPosIndexR, QTLPosIndexC = size(QTLPosIndex)
	sortepistaticQTLPosR = size(sortepistaticQTLPos)[1]
	for i in 1:sortepistaticQTLPosR    
		for k in 1:QTLPosIndexR            
			if sortepistaticQTLPos[i] == QTLPosIndex[k,1]
				sortepistaticQTLPos[i] = QTLPosIndex[k,2]
			end            
		end    
	end
	return sortepistaticQTLPos
end
function TruePos(epi_combinations,QTLPosIndex)
	QTLPosIndexNew = [1:size(QTLPosIndex)[1] QTLPosIndex[:,2]]
	QTLPosIndexNewR, QTLPosIndexNewC = size(QTLPosIndexNew)
	epi_combinationsR,epi_combinationsC = size(epi_combinations)
	for i in 1:epi_combinationsR
		for j in 1:epi_combinationsC
			for k in 1:QTLPosIndexNewR            
				if epi_combinations[i,j] == QTLPosIndexNew[k,1]
					epi_combinations[i,j] = QTLPosIndexNew[k,2]
				end
			end
		end    
	end
	Mnew =  Array{Float64,1}(undef,epi_combinationsR)
	for i in 1:epi_combinationsR
		Mnew[i] = median([epi_combinations[i,1], epi_combinations[i,2]])
	end
	return Mnew
end
function epi_combination(SNPpairs, intperloc=1)
	breaker = 0
	while breaker == 0
		intmatrix = zeros(Int8, 2*SNPpairs, 2*SNPpairs)
		loclist = collect(1:1:2*SNPpairs)
		#global loclist = epistaticQTLPos 
		nint = zeros(Int64, 2*SNPpairs)   #track number of interactions per qtl
		for i in loclist  
			if size(nint[nint .< intperloc])[1] == 0  #break when all is 5
				#println("Done")
				breaker = 1
				break
			end  
			ll = shuffle(loclist[nint .< intperloc])   #randomize the list of qtl
			ll = ll[ll .!= i]    
			if nint[i] < intperloc    #if the qtl has less than 5 interactions, sample some new ones
				nnew = intperloc - nint[i]      
				if nnew > 0    #there may be too few loci left
					if nnew > size(ll)[1]
						println("Too few loci left")
						@goto label1
					end        
					smp = ll[1:nnew]
					intmatrix[smp,i] .= 1
					nint[smp] .+= 1        
				end
			end    
			nint[i] = intperloc    
		end
		@label label1
	end
	#check the sums
	rowsums = convert(Array, sum(intmatrix,dims=1))
	colsums = convert(Array, sum(intmatrix,dims=2))
	allsums = rowsums .+ colsums'
	size(allsums[allsums .== intperloc])
	#find the combinations
	indx = findall(x-> x == 1, intmatrix)
	indx1 = [ x[1] for x in indx]
	indx2 = [ x[2] for x in indx]
	#randomize the numbers
	locrand = shuffle(loclist)
	indx1 = locrand[indx1]
	indx2 = locrand[indx2]
	epi_combinations = [indx1 indx2]
	return epi_combinations
end
function Dominant(QTLPos,QTLd,α,domval = "medium")
	DominantQTLPos = sample(QTLPos,QTLd,replace=false)
	domval = "medium"
	if domval == "none"
		dommean = 0.0
		domsd = 0.0
	elseif domval == "low"
		dommean = 0.2
		domsd = 0.3
	elseif domval == "medium"
		dommean = 0.2
		domsd = 0.7
	else
		dommean = 0.2
		domsd = 1.5
	end
	#dominance
	if domsd .!= 0.0
		domcoeff = rand(Normal(dommean,domsd), QTLd)		
	else
		domcoeff = zeros(Float64, QTLd)
	end
	SortQTLDompos = sort(DominantQTLPos)
	SortQTLpos = sort(QTLPos)
	d = [abs(x) for x in α[findall(x->x in SortQTLDompos, QTLPos)]] .* domcoeff   # .* sample([-1,1],QTLd)
	return d, SortQTLDompos, SortQTLpos
end
function epistatic(QTLPos, SNPpairs, intperloc=1, epival = "medium")
	if epival == "none"
		epimean = 0.0
		episd = 0.0
	elseif epival == "low"
		epimean = 0.0
		episd = sqrt( (0.2^2 + 0.3^2)/intperloc )
	elseif epival == "medium"
		epimean = 0.0
		episd = sqrt( (0.2^2 + 0.7^2)/intperloc )
	else
		epimean = 0.0
		episd = sqrt( (0.2^2 + 1.5^2)/intperloc )
	end	
	ncombs = Int(ceil(2*SNPpairs*intperloc/2))
	if episd != 0.0
		epi_coef = rand(Normal(epimean, episd), ncombs)		
	else
		epi_coef = zeros(Float64, ncombs)
	end
	epistaticQTLPos = sample(QTLPos,2*SNPpairs,replace=false)		
	breaker = 0	
	while breaker == 0		
		global intmatrix = zeros(Int8, 2*SNPpairs, 2*SNPpairs)
		global loclist = collect(1:1:2*SNPpairs)		 
		nint = zeros(Int64, 2*SNPpairs)   #track number of interactions per qtl
		for i in loclist  
			if size(nint[nint .< intperloc])[1] == 0  #break when all is 5				
				breaker = 1
				break
			end  
			ll = shuffle(loclist[nint .< intperloc])   #randomize the list of qtl
			ll = ll[ll .!= i]    
			if nint[i] < intperloc    #if the qtl has less than 5 interactions, sample some new ones
				nnew = intperloc - nint[i]      
				if nnew > 0    #there may be too few loci left
					if nnew > size(ll)[1]
						println("Too few loci left")
						@goto label1
					end        
					smp = ll[1:nnew]
					intmatrix[smp,i] .= 1
					nint[smp] .+= 1        
				end
			end    
			nint[i] = intperloc    
		end
		@label label1
	end
	#check the sums
	rowsums = convert(Array, sum(intmatrix,dims=1))
	colsums = convert(Array, sum(intmatrix,dims=2))
	allsums = rowsums .+ colsums'
	size(allsums[allsums .== intperloc])
	#find the combinations
	indx = findall(x-> x == 1, intmatrix)
	indx1 = [ x[1] for x in indx]
	indx2 = [ x[2] for x in indx]
	#randomize the numbers
	locrand = shuffle(loclist)
	indx1 = locrand[indx1]
	indx2 = locrand[indx2]
	epi_combinations = [indx1 indx2]	
	return epistaticQTLPos, epi_coef, epi_combinations  #, epi_combinations, epi_coef
end
function epistaticCom(noia, epi_combinations, noia_dom, nObs, epimodel = "none")
	if epimodel == "none"
		aa = zeros(Float32, nObs)
		ad = zeros(Float32, nObs)
		da = zeros(Float32, nObs)
		dd = zeros(Float32, nObs)
	else
		a1 = noia[:,epi_combinations[:,1]]
		a2 = noia[:,epi_combinations[:,2]]
		d1 = noia_dom[:,epi_combinations[:,1]]
		d2 = noia_dom[:,epi_combinations[:,2]]
		aa = a1 .* a2 * epi_effects[:,5]
		ad = a1 .* d2 * epi_effects[:,8]
		da = d1 .* a2 * epi_effects[:,6]
		dd = d1 .* d2 * epi_effects[:,9]
	end
	return aa    #, ad, da, dd
end
function Zmatrix(M, p)
	Row, Col = size(M)
	Z = Array{Float64,2}(undef, Row, Col)
	W = Array{Float64,2}(undef, Row, Col)
	for i = 1:Col
		for j = 1:Row
			if M[j,i] .== 2
				Z[j,i] = 2 .- 2 .* p[i]
				W[j,i] = 0 .- 2 .* (1 .- p[i])^2
			elseif M[j,i] .== 1
				Z[j,i] = 1 .- 2 .* p[i]
				W[j,i] = 2 .* p[i] .* (1 .- p[i])
			elseif M[j,i] .== 0
				Z[j,i] = 0 .- 2 .* p[i]
				W[j,i] = 0 .- 2 .* p[i]^2
			else
				Z[j,i] = 9999
				W[j,i] = 9999
			end
		end	
	end
	return Z, W
end
function sdCor2(A,B,n_covariate=0)
	ss = size(A, 1)
	A2 = A.^2
	B2 = B.^2
	Exij = A' * B ./ss 
	Exij2 = A2' * B2 ./ ss
	if n_covariate != 0
		cc = Array{Float64,2}(undef, size(A, 2), size(A, 2))
		for i = 1:size(n_covariate, 2)
			tmp = A' * (B .* n_covariate[:,i])./(ss .- 1)
			cc = cc .+ tmp .^ 2
		end
		VARxij = (ss ./ (ss.-1)) .* (Exij2 .- Exij .^ 2)
		SDxij = sqrt.(VARxij .- cc)
	else
		SDxij = sqrt.(ss./(ss .- 1)) .* sqrt.(Exij2 .- Exij .^ 2)
	end
	return SDxij
end
function innerCor2(A,B,nY)
	ss = size(A, 1)
	if size(nY, 2) == 1
		r = A' * (B .* nY[:,1]) ./ (ss .- 1)
	else
		r = A' * (B .* nY) ./ (ss .- 1)
	end
	return r
end
function matrixCor(nX,nY)
	return nX' * nY ./ (size(nX, 1) - 1)
end
function matrixPval(r,df)
	t = sqrt.(df) .* abs.(r) ./ sqrt.(abs.(1 .- r .^2))
	dist = Distributions.TDist(df)	
	return 2 .* (1 .- Distributions.cdf(dist, t))
end
function standardizes_norm(x)
	norm = (x .- mean(x))./ sd.(x)
	replace!(norm, NaN => 0)
	return norm
end
function matrixPar(r12,r13,r23)
	return (r12 .- r13 .* r23) ./ sqrt.(abs.(1 .- r13 .^2) .* abs.(1 .- r23 .^2))
end
function removeNAN(x, v=0)
	return replace!(x, NaN => v)
end
function as_dummy(x::AbstractVector)
	u = unique(x)
	m = Matrix{Bool}(undef, length(x), length(u))
	for i in eachindex(u)
		@. m[:, i] = x .== u[i]
	end
    return m
end
#p2c converts the given p value to the absolute value of corresponding correlation coefficient
function p2c(pval,df)
	dist = Distributions.TDist(df)
	t = 1 - Distributions.cdf(dist, pval / 2)
	return sqrt(t ^ 2 / (t ^ 2 + df))
end
function smoothQQplot(x, y)
	#import Gadfly.ElementOrFunction
	#Gadfly.plot(qq::QQPair, elements::ElementOrFunction...) = Gadfly.plot(x=qq.qx, y=qq.qy, Geom.point, Theme(highlight_width=0px), elements...)
	#qqplot(x, y, elements::ElementOrFunction...) = Gadfly.plot(qqbuild(x, y), elements...)
	if typeof(x) == Array{Float64, 2} x = UpperTriangular(x) end
	if typeof(y) == Array{Float64, 2} y = UpperTriangular(y) end
	replace!(x, NaN => 0)
	replace!(y, NaN => 0)
	x = -log10.(x)
	y = -log10.(y)	
	l = minimum([x; y])
	r = ceil(maximum([x; y]))
	#RowX = size(x, 1)
	#RowY = size(y, 1)
	#"deepskyblue"
	myplot = Gadfly.plot(layer(x=x,y=y, Geom.point, Theme(point_size = 8pt, default_color=colorant"blue",alphas=[0.5], highlight_width = 1pt), order=1), Guide.xlabel("-log10 Pvalue"), Guide.ylabel("-log10 Pvalue"),
	                     layer(x=[l,r], y=[l,r], Geom.line, Theme(default_color=colorant"red",style(line_width=0.5mm)), order=2),
						 layer(x=quantile(x, LinRange(0, 1, 1000)),y=quantile(y, LinRange(0, 1, 1000)),Geom.point,Theme(default_color=colorant"red",point_size = 5pt, highlight_width = 0pt), order=3),
						 Coord.cartesian(xmin=l, xmax=r, ymin=l, ymax=r)						 
						 )
	draw(PNG("smoothQQplot.png", 8inch, 8inch), myplot)
end
function name_SNP(X)
    name = collect(1:size(X,2))
    SNP_name = Matrix{String}(undef,size(X,2),1)
    for i = 1:size(X,2)
        SNP_name[i] = string("SNP",name[i])
    end
    return SNP_name
end
function SNP_pairs_method(snpA, snpB, trait, covariate=false)
    trait = Float64.(trait)
    n_covariate = 0
    if covariate==false
        n_trait = standardize(ZScoreTransform,trait,dims=1)
        n_snpA = standardize(ZScoreTransform,snpA,dims=2)
        n_snpB = standardize(ZScoreTransform,snpB,dims=2)
        df1 = size(snpA, 1)  - 2  #df for epistasis
        df2 = size(snpA, 1) - 4  #df for additive and epistasis
    else
        pc_covariate = svd(covariate).U
        n_covariate = standardize(ZScoreTransform,pc_covariate,dims=2)
        dat = DataFrame(trait = trait,n_covariate = n_covariate)
        model = lm(@formula(trait ~ 1 + n_covariate), dat)
        res_trait = dispersion(model.model)
        n_trait = standardize(ZScoreTransform,res_trait,dims=1)
        #snpA
        h_snpA = snpA - mean(snpA[:,2])
        res_snpA = h_snpA - n_covariate *(n_covariate' * h_snpA)/(size(snpA, 1)-1)
        n_snpA = standardize(ZScoreTransform,res_snpA,dims=2)
        #snpB
        h_snpB = snpB - mean(snpB[:,2])
        res_snpB = h_snpB - n_covariate *(n_covariate' * h_snpB)/(size(snpB, 1)-1)
        n_snpB = standardize(ZScoreTransform,res_snpB,dims=2)
        df1 = size(snpA, 1) - size(covariate, 2) - 2  #df for only epistasis with adjustment of covariates
        df2 = size(snpA, 1) - size(covariate, 2) - 4  #df for additive and epistasis with adjustment of covariates
    end    
    SDxij = sdCor2(snpA,snpB,n_covariate)
    r_xy_t =  innerCor2(snpA,snpB,n_trait) ./ SDxij
    r_xy_y = innerCor2(snpA,snpB,n_snpB) ./ SDxij
    r_xy_x = transpose(innerCor2(snpB,snpA,n_snpA)) ./ SDxij
    r_x_t = matrixCor(n_snpA,n_trait)
    #r_x_t =   reshape(matrixCor(n_snpA,n_trait) , size(snpA,2) , size(snpB,2))
    r_y_t = transpose(matrixCor(n_snpB,n_trait))
    #r_y_t = transpose(matrix(matrixCor(n_snpB,n_trait), size(snpB,2) , size(snpA,2)))
    #r_y.t = matrix( matrixCor(n_snpB,n_trait) , nrow=ncol(snpA) , ncol=ncol(snpB) , byrow=TRUE )
    r_x_y = matrixCor(n_snpA,n_snpB)
    r_xy_t_x = matrixPar(r_xy_t, r_xy_x, r_x_t)
    r_xy_y_x = matrixPar(r_xy_y, r_xy_x, r_x_y)
    r_y_t_x = matrixPar(r_y_t, r_x_t, r_x_y)
    r_xy_t_x_y = matrixPar(r_xy_t_x, r_xy_y_x, r_y_t_x)
    r_xy_t = removeNAN(r_xy_t)
    r_xy_t_x_y = removeNAN(r_xy_t_x_y)    
    return r_xy_t, df1 
end
function Single_point_GWAS(trait,Markers,covariate=false)
    M_col = size(Markers,2)
    trait = Float64.(trait)
    Markers = Float64.(Markers)
    #data = DataFrame(Markers = Markers,trait =trait)
    n_covariate = 0
    cc = Array{Float64,2}(undef, M_col, 5)
    n_Markers = standardize(ZScoreTransform,Markers,dims=2)
    if covariate==false
        n_trait = standardize(ZScoreTransform,trait,dims=1)                
        for i = 1 : M_col
            data = DataFrame(trait =n_trait,Markers = n_Markers[:,i])
            ols = lm(@formula(trait ~ 1 + Markers), data)            
            cc[i,1] = i
            cc[i,2] = coeftable(ols).cols[1][2]
            cc[i,3] = coeftable(ols).cols[2][2]
            cc[i,4] = coeftable(ols).cols[3][2]
            cc[i,5] = coeftable(ols).cols[4][2]
        end       
    else
        data1 = DataFrame(trait = trait,covariate = covariate)
        model = lm(@formula(trait ~ 1 + covariate), data1)
        res_trait = dispersion(model.model)
        n_trait = standardize(ZScoreTransform,res_trait,dims=1)        
        for i = 1 : M_col
            data = DataFrame(trait =n_trait,Markers = n_Markers[:,i])
            ols = lm(@formula(n_trait ~ 1 + n_Markers(:,i)), data)            
            cc[i,1] = i
            cc[i,2] = coeftable(ols).cols[1][2]
            cc[i,3] = coeftable(ols).cols[2][2]
            cc[i,4] = coeftable(ols).cols[3][2]
            cc[i,5] = coeftable(ols).cols[4][2]
        end
    end    
    replace!(cc, NaN => 0)
    replace!(cc, Inf => 0)
    cc = DataFrame(SNP=Int64.(cc[:,1]),Est=cc[:,2],SE=cc[:,3],t=cc[:,4],Pvalu=cc[:,5])
    return cc
end
function SNP_SNP_GWAS1(trait,Markers1,Markers2,covariate=false)
    M_col1 = size(Markers1,2)
    M_col2 = size(Markers2,2)
    trait = Float64.(trait)
    Markers1 = Float64.(Markers1)
    Markers2 = Float64.(Markers2)    
    n_covariate = 0
    cc = Array{Float64,2}(undef, M_col1*M_col2, 6)
    n_Markers1 = standardize(ZScoreTransform,Markers1,dims=2)
    n_Markers2 = standardize(ZScoreTransform,Markers2,dims=2)
    n = 0
    if covariate==false
        n_trait = standardize(ZScoreTransform,trait,dims=1)
        for i = 1 : M_col1
            for j = 1 : M_col2
                data = DataFrame(trait =n_trait,Markers1 = n_Markers1[:,i],Markers2 = n_Markers2[:,j])
                ols = lm(@formula(trait ~ 1 + Markers1 * Markers2), data)
                n = n + 1
                cc[n,1] = i
                cc[n,2] = j
                cc[n,3] = coeftable(ols).cols[1][4]
                cc[n,4] = coeftable(ols).cols[2][4]
                cc[n,5] = coeftable(ols).cols[3][4]
                cc[n,6] = coeftable(ols).cols[4][4]
            end
        end
    else
        data1 = DataFrame(trait = trait,covariate = covariate)
        model = lm(@formula(trait ~ 1 + covariate), data1)
        res_trait = dispersion(model.model)
        n_trait = standardize(ZScoreTransform,res_trait,dims=1) 
        for i = 1 : M_col1
            for j = 1 : M_col2
                data = DataFrame(trait =n_trait,Markers1 = n_Markers1[:,i],Markers2 = n_Markers2[:,j])
                ols = lm(@formula(trait ~ 1 + Markers1 * Markers2), data)
                n = n + 1
                cc[n,1] = i
                cc[n,2] = j
                cc[n,3] = coeftable(ols).cols[1][4]
                cc[n,4] = coeftable(ols).cols[2][4]
                cc[n,5] = coeftable(ols).cols[3][4]
                cc[n,6] = coeftable(ols).cols[4][4]
            end
        end
    end
    replace!(cc, NaN => 0)
    replace!(cc, Inf => 0)
    cc = DataFrame(SNP1=Int64.(cc[:,1]),SNP2=Int64.(cc[:,2]),Est=cc[:,3],SE=cc[:,4],t=cc[:,5],Pvalu=cc[:,6])
    return cc
end
function SNP_SNP_GWAS2(trait,Markers,covariate=false)
    M_col = size(Markers,2)    
    trait = Float64.(trait)
    Markers = Float64.(Markers)    
    n_covariate = 0
    cc = Array{Float64,2}(undef, Int64.(((1 + M_col) * M_col)/2), 6)
    n_Markers = standardize(ZScoreTransform,Markers,dims=2)    
    n = 0
    if covariate==false
        n_trait = standardize(ZScoreTransform,trait,dims=1)
        for i = 1 : M_col
            for j = 1 : i
                data = DataFrame(trait =n_trait,Markers1 = n_Markers[:,i],Markers2 = n_Markers[:,j])
                ols = lm(@formula(trait ~ 1 + Markers1 * Markers2), data)
                n = n + 1
                cc[n,1] = i
                cc[n,2] = j
                cc[n,3] = coeftable(ols).cols[1][4]
                cc[n,4] = coeftable(ols).cols[2][4]
                cc[n,5] = coeftable(ols).cols[3][4]
                cc[n,6] = coeftable(ols).cols[4][4]
            end
        end
    else
        data1 = DataFrame(trait = trait,covariate = covariate)
        model = lm(@formula(trait ~ 1 + covariate), data1)
        res_trait = dispersion(model.model)
        n_trait = standardize(ZScoreTransform,res_trait,dims=1) 
        for i = 1 : M_col
            for j = 1 : i
                data = DataFrame(trait =n_trait,Markers1 = n_Markers[:,i],Markers2 = n_Markers[:,j])
                ols = lm(@formula(trait ~ 1 + Markers1 * Markers2), data)
                n = n + 1
                cc[n,1] = i
                cc[n,2] = j
                cc[n,3] = coeftable(ols).cols[1][4]
                cc[n,4] = coeftable(ols).cols[2][4]
                cc[n,5] = coeftable(ols).cols[3][4]
                cc[n,6] = coeftable(ols).cols[4][4]
            end
        end
    end
    replace!(cc, NaN => 0)
    replace!(cc, Inf => 0)
    cc = DataFrame(SNP1=Int64.(cc[:,1]),SNP2=Int64.(cc[:,2]),Est=cc[:,3],SE=cc[:,4],t=cc[:,5],Pvalu=cc[:,6])
    return cc
end