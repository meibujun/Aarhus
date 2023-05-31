#EXECUTE 
using DataFrames, SparseArrays, LinearAlgebra, Random, JWAS, CSV, Distributions
using Printf

mutable struct ModelTerm
    trmString::AbstractString
    iModel::Int64                           # the model (trait) this term belongs to
    startPos::Int64
    endPos::Int64
    randomType::String
end

# Composite type for including SNP effects in the model
mutable struct Genotypes                    # ZM where Z is a sparse incidence matrix and M is a dense matrix of covariates           
    trmString::AbstractString               # "Ztrm"
    Z::SparseMatrixCSC{Float64,Int64}       # incidence matrix Z to map rows of M to phenotypes according to "Ind" 
    M::Array{Float64,2}                     # dense matrix of covariates
    MArray::Array{Any,1}                    # column views of M
    MPMArray::Array{Float64,1}              # diagonals of M'M
    idGeno::Array{String}                   # row names of M (Ind must be a subset of idGeno)
    α::Array{Float64,2}                     # nCovariates x nTraits matrix of effects
    αMean::Array{Float64,2}                 # mean of sampled values  
    yAdjArray::Array{Any,1}                 # Array of Arrays that will be used to access adjusted y-vector
    νpr::Float64
    Spr::Array{Float64,2}
    Vi::Array{Float64,2}
    meanV::Array{Float64,2}
    estimate::Bool
    outSamples::Bool
    outFreq::Int64
    outStream   
end 

mutable struct RandomEffect
    modelTermVec::Array{ModelTerm,1}      
    νpr::Float64
    Spr::Array{Float64,2}
    Vi::Array{Float64,2}
    meanV::Array{Float64,2}
    Ai
    estimate::Bool
    outSamples::Bool
    outFreq::Int64
    outStream
end

mutable struct MME
    modelEquations::Array{AbstractString,1} #"[y1 = A + B + A*B ;y2 = A + B"]
    covVec::Array{String}                   #["age"]
    modelTermVec::Array{ModelTerm,1}        #[modelTerm("A") , modelTerm("A*B")]
    modelTermDict::Dict{AbstractString,ModelTerm}
    genotypes                               # dense component of model
    randomEffectsVec::Array{RandomEffect,1} #[RandomEffect("Ind"), RandomEffect("Mat")]
    ped                                     # false or a pedigee structure defined in JWAS.PedModule
    depVarVec::Array{AbstractString,1}      #["y1","y2"]
    X::SparseMatrixCSC{Float64,Int64}       #design matrix
    mmeLhs::SparseMatrixCSC{Float64,Int64}  #left-hand side for MME
    mmeRhs::Array{Float64,1}                #right-hand side for MME
    mmeSpl::Array{Float64,1}                #sample of location parameters
    meanEffects::Array{Float64,1}           #mean of location parameters
    varRes::Array{Float64,2}                #residual variance 
    νRes::Float64                           #prior degrees of freedom for varRes
    SRes::Array{Float64,2}                  #prior scale parameter for resVar
    y::Array{Float64,1}                     #dependent variable
    yAdj::Array{Float64,1}
    meanVarRes::Array{Float64,2}                     
end 

function ModelTerm(str::AbstractString,iModel)
    ModelTerm(str,iModel,0,0,"fixed")
end

#EXECUTE 
function initMME(modelEquations::AbstractString,varRes;νRes=5.0)  # "y1 = A + A*B; y2 = A + B + !(A)M
    if modelEquations==""
        error("modelEquations is empty\n")
    end
    modelVec   = [strip(i) for i in split(modelEquations,[';','\n'],keepempty=false)]
    nModels    = size(modelVec,1)
    depVarVec  = Array{AbstractString,1}()
    modelTermVec  = Array{ModelTerm,1}()                 # modelTerms across all the models
    modelTermDict = Dict{AbstractString,ModelTerm}()  
    for (m,model) = enumerate(modelVec)
        leftRight  = split(model,"=")                 # "y", "A+A*B"
        depVarVec  = [depVarVec; strip(leftRight[1])] # ["y1","y2"]
        modelParts = strip(leftRight[2])              # "A+A*B" 
        termsVec   = split(modelParts,"+")            # "A","A*B"
        modelTermVeci = [ModelTerm(string(m)*":"*strip(trmStr),m) for trmStr in termsVec]
        modelTermVec  = [modelTermVec;modelTermVeci]    
    end
    for i in modelTermVec
        modelTermDict[i.trmString] = i
    end
    genotypes = false
    randomEffectsVec = []
    ped = false
    covVec = []
    X = spzeros(0,0)
    mmeLhs = spzeros(0,0)
    mmeRhs = []
    mmeSpl = []
    meanEffects = []
    y      = Array{Float64,1}(undef,0)
    yAdj   = Array{Float64,1}(undef,0)
    m = size(varRes,1)
    varRes = m==1 ? fill(varRes,1,1) : varRes
    SRes   = varRes*(νRes - m - 1)
    meanVarRes = zero(varRes)
    return MME(modelVec,covVec,modelTermVec,modelTermDict,genotypes,randomEffectsVec,ped,
               depVarVec,X,mmeLhs,mmeRhs,mmeSpl,meanEffects,varRes,νRes,SRes,y,yAdj,meanVarRes)
end

#EXECUTE
function addGenotypes!(mme::MME,df::DataFrame,ZTrm::AbstractString,M::Array{Float64,2},
            idGeno::Array{String,1},V;
            estimate::Bool=false,
            νPrior::Float64=4.1,
            outSamples::Bool=false,
            outFreq::Int64=100
        )
        m = size(V,1)        
        if outSamples==true
            fileName = "G.samples"
            outStream = open(fileName,"w")
            for i=1:m,j=i:m
                if i==j==1
                    @printf(outStream,"%6s", "V[j]")
                else
                    @printf(outStream," %6s", "V[j]")    
                end        
            end
            @printf(outStream,"\n")        
        else
            outStream = nothing
        end        
        var = m==1 ? fill(V,1,1) : V
        d = Dict()
        for (i,s) = enumerate(idGeno)
            d[s] = i
        end
        n = size(df,1) 
        ii = collect(1:n) 
        jj = Array{Int64,1}(undef,n)
        idPheno = string.(df[!,Symbol(ZTrm)])    
        for (i,v) = enumerate(idPheno)
            jj[i] = d[v]
        end
        Z = sparse(ii,jj,1.0,n,size(idGeno,1))
        α = zeros(0,0)
        αMean = zeros(0,0)    
        Spr = (νPrior - m - 1)*var 
        yAdjArray = Array{Array{Float64,1},1}(undef,0)
        MArray = Array{Any,1}(undef,0)
        MPMArray = Array{Any,1}(undef,0)  
        mme.genotypes = Genotypes(ZTrm,Z,M,MArray,MPMArray,idGeno,α,αMean,yAdjArray,νPrior,Spr,
            inv(var),zero(var),estimate,outSamples,outFreq,outStream)
        return
end

#EXECUTE 
# This function returns returns a dictionary with the names in the pedigree as the keys and their 
# sequential numbers as the associated values
# It also returns vector with the keys in sequential order. 

function mkDict(ped::JWAS.PedModule.Pedigree)
    d = Dict()
    names = Array{String}(undef,length(ped.idMap))
    for i in ped.idMap    
        d[i.first] = i.second.seqID
        names[i.second.seqID] = i.first
    end
    return d,names
end

# This function returns 
# a dictionary with the unique values in the vector "a" as the keys and their 
# sequential numbers as the associated values
# It also returns vector with the keys in sequential order. 
function mkDict(a)
    d = Dict()
    aUnique = unique(a)
    names = Array{String}(undef,size(aUnique,1))
    for (i,s) in enumerate(aUnique)
    names[i] = s
    d[s] = i
    end
    return d,names
end

function getX(mme,modelTerm::ModelTerm,covariables,df)
    n = size(df,1)
    nModels = size(mme.modelEquations,1)
    trmString = split(modelTerm.trmString,':')[2]
    if trmString == "intercept"
        ii = (modelTerm.iModel-1)*n .+ (1:n) # row numbers
        jj = ones(n)
        val = 1.0
        p = 1
        X  = sparse(ii,jj,val,n*nModels,p)
        colNames = modelTerm.trmString
        return X,colNames
    end
    factors = strip.(split(trmString,"*"))
    covs = [i in covariables for i in factors]
    
    if covs[1] == false
        str = string.(df[:,Symbol(factors[1])])
        val = 1.0
    else
        str = fill(factors[1],n) 
        val = df[:,Symbol(factors[1])]    
    end       

    for i in 2:length(factors)
        if covs[i] == false
            str = str .*" x ".*string.(df[:,Symbol(factors[i])])
            val = val .* 1.0 
        else
            str = str .*" x ".*fill(factors[i],n) 
            val = val .* df[:,Symbol(factors[i])]    
        end 
    end
    dict,colNames   = modelTerm.randomType != "A" ? mkDict(str) : mkDict(mme.ped)
    str = modelTerm.randomType != "A" ? str : string.(df[:, Symbol(factors[factors .!= mme.covVec][1])])
    ii = (modelTerm.iModel-1)*n .+ (1:n)      # row numbers 
    jj = [dict[i] for i in str]  # column numbers
    p = length(colNames)        
    X  = sparse(ii,jj,val,n*nModels,p)
    return X, strip(modelTerm.trmString)*": ".*colNames   
end

function setRandom!(mme::MME,trmStrings::String,V,Ai;
                estimate::Bool=false,
                νPrior::Float64=4.1,
                outSamples::Bool=false,
                outFreq::Int64=100
            )
    m = size(V,1)        
    if outSamples==true
        fileName = "V($trmStrings).samples"
        outStream = open(fileName,"w")
        for i=1:m,j=i:m
            if i==j==1
                @printf(outStream,"%6s", "V[j]")
            else
                @printf(outStream," %6s", "V[j]")    
            end        
        end
        @printf(outStream,"\n")        
    else
        outStream = nothing
    end        
    var = m==1 ? fill(V,1,1) : V        
    trmVec = strip.(split(trmStrings,[',';' '],keepempty=false))
    modelTerms = []
    for modelTerm in mme.modelTermVec
        trmString = split(modelTerm.trmString,':')[2]
        if trmString in trmVec
            push!(modelTerms,modelTerm)
        end                
    end
    if length(modelTerms) != size(var,1)
        println("size of var: (length(modelTerms))")
        return
    end 
	[trm.randomType = "Kmatrix" for trm in modelTerms]

    Spr = (νPrior - m - 1)*var    
    randomEffect = RandomEffect(modelTerms,νPrior,Spr,inv(var),zero(var),Ai,
                        estimate,outSamples,outFreq,outStream)
    push!(mme.randomEffectsVec,randomEffect)   
end

function addGiMats!(mme::MME)
    for randomEffect in mme.randomEffectsVec
        for (i,modelTermi) in enumerate(randomEffect.modelTermVec), (j,modelTermj) in enumerate(randomEffect.modelTermVec)
            starti = modelTermi.startPos
            startj = modelTermj.startPos
            endi = modelTermi.endPos
            endj = modelTermj.endPos
            mme.mmeLhs[starti:endi,startj:endj] = mme.mmeLhs[starti:endi,startj:endj] + pinv(randomEffect.Ai)*randomEffect.Vi[i,j]
        end
    end
end


#EXECUTE
function getLhsRhs!(mme,df)
    n = size(df,1)
    X,colNames = getX(mme,mme.modelTermVec[1],mme.covVec,df)
    mme.modelTermVec[1].startPos = 1
    mme.modelTermVec[1].endPos  = mme.modelTermVec[1].startPos + size(X,2) - 1
    for i = 2:size(mme.modelTermVec,1)
        Xi,namesi = getX(mme,mme.modelTermVec[i],mme.covVec,df)
        X = [X Xi]
        mme.modelTermVec[i].startPos = mme.modelTermVec[i-1].endPos + 1
        mme.modelTermVec[i].endPos   = mme.modelTermVec[i].startPos + size(Xi,2) - 1
        colNames = [colNames; namesi]
    end
    y = Array{Float64,1}()
    for v in mme.depVarVec
        y = [y;df[:,Symbol(v)]]
    end
    Ri = kron(inv(mme.varRes),sparse(I, n, n))
    mme.X = X
    mme.mmeLhs = X'Ri*X
    mme.mmeRhs = X'Ri*y
    mme.mmeSpl = zero(mme.mmeRhs)
    mme.meanEffects = zero(mme.mmeRhs)
    mme.y      = y
    mme.yAdj = copy(mme.y)
    addGiMats!(mme)
    if mme.genotypes != false
        mme.genotypes.M = mme.genotypes.M .- mean(mme.genotypes.M,dims=1)
        mme.genotypes.M = mme.genotypes.Z * mme.genotypes.M
        nModels  = size(mme.varRes,1)
        nMarkers = size(mme.genotypes.M,2)
        yAdjArray = Array{Any,1}(undef,nModels)
        for i=1:nModels
            startPos = (i-1)*n + 1
            endPos = startPos + n - 1
            yAdjArray[i] = @view mme.yAdj[startPos:endPos]
        end
        MArray   = Array{Any,1}(undef,nMarkers)
        MPMArray = Array{Float64,1}(undef,nMarkers)
        for i=1:nMarkers
            m = @view mme.genotypes.M[:,i]
            MArray[i] = m
            MPMArray[i] = m'm
        end
        mme.genotypes.yAdjArray = yAdjArray
        mme.genotypes.MArray = MArray
        mme.genotypes.MPMArray = MPMArray
        mme.genotypes.α = zeros(nMarkers,nModels)
        mme.genotypes.αMean = zeros(nMarkers,nModels)
    end
    return mme.mmeLhs,mme.mmeRhs,colNames
end
function updateLhsRhs!(mme)
    m = size(mme.varRes,1)
    n = Int(size(mme.y,1)/m)
    Ri = kron(inv(mme.varRes),sparse(I, n, n))
    mme.mmeLhs = mme.X'Ri*mme.X
    mme.mmeRhs = mme.X'Ri*mme.yAdj
    addGiMats!(mme)
end
