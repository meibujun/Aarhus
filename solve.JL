#The method of successive over-relaxation
function sor(mat::Array{Float64},sol::Array{Float64},rhs::Array{Float64},
             wgt::Float64=1.1,
             maxit::Int=5000,tol::Float64=1.0e-8) 
    nbrows, nbcols = size(mat)
    result = deepcopy(sol)
    deltax = zeros(nbrows, 1)
    numit = 0; nrmdx = 1;
    while numit < maxit
        numit = numit + 1
        for i=1:nbrows
            deltax[i] = wgt*rhs[i]
            for j=1:nbcols
                deltax[i] = deltax[i] - mat[i,j]*wgt*result[j]
            end
            deltax[i] = deltax[i]/mat[i,i]
            result[i] = result[i] + deltax[i]
        end
        nrmdx = norm(deltax)
        strdx = @sprintf("%.2e", nrmdx)
        #println("||dx|| = $strdx")
        if norm(deltax) <= tol
            return result  #(result, numit, nrmdx, false)
        end
    end
    return result  #(result, numit, nrmdx, true)
end

#The Conjugate Gradient Method
function CGM(A::Array{Float64,2},x0::Array{Float64,1},b::Array{Float64,1},             
             maxit::Int64=5000,tol::Float64=1.0e-8,
             verbose=false)
    sol = deepcopy(x0)
    r = b - A*sol
    p = deepcopy(r)
    if verbose
        println(" norm(r)    alpha       beta")
    end
    for i=1:maxit
        res = norm(r)
        if verbose
            sres = @sprintf("%.2e", res)
            print("$sres")
        end 
        if res < tol
            if verbose
                println("  succeeded after ", i, " steps")
            end
            return sol  #(sol, res, i, false)
        end
        alpha = (transpose(r)*r)/(transpose(p)*A*p)
        if verbose
            salpha = @sprintf("%.4e", alpha)
            print("  $salpha")
        end
        sol = sol + alpha*p
        r1 = r - alpha*A*p
        beta = (transpose(r1)*r1)/(transpose(r)*r)
        if verbose
            sbeta = @sprintf("%.4e", beta)
            println("  $sbeta")
        end
        p = r1 + beta*p
        r = r1
    end
    return sol  #(sol, norm(r), maxit, true)
end

#Jacobi, Gauss-Seidel, and Gibbs (general,lambda,one iteration)
function Jacobi(A,x,b,p=0.7;tolerance=0.000001,printout_frequency=10,maxiter=1000)
    n       = size(A,1)   #number of linear equations
    D       = diag(A)
    error   = b - A*x
    diff    = sum(error.^2)/n

    iter    = 0
    while (diff > tolerance) && (iter < maxiter)
        iter   += 1
        error   = b - A*x
        x_temp  = error./D + x
        x       = p*x_temp + (1-p)*x
        diff    = sum(error.^2)/n

        if iter%printout_frequency == 0
            println("at iteration ",iter,": ",diff)
        end
    end
    return x
end

function GaussSeidel(A,x,b;tolerance=0.000001,printout_frequency=10,maxiter=1000)
    n = size(A,1)
    for i = 1:n
        x[i] = (b[i] - A[:,i]'x)/A[i,i] + x[i]
    end
    error = b - A*x
    diff  = sum(error.^2)/n

    iter  = 0
    while (diff > tolerance) & (iter < maxiter)
        iter += 1
        for i = 1:n
            x[i] = (b[i] - A[:,i]'x)/A[i,i] + x[i]
        end

        error = b - A*x
        diff  = sum(error.^2)/n
        if iter%printout_frequency == 0
            println("at iteration ",iter,": ",diff)
        end
    end
    return x
end

function gauss_seidel(mat::Array{Float64},rhs::Array{Float64},
                      sol::Array{Float64},
                      maxit::Int=1000,tol::Float64=1.0e-8) 
    nbrows, nbcols = size(mat)
    result = deepcopy(sol)
    deltax = zeros(nbrows, 1)
    numit = 0; nrmdx = 1;
    while numit < maxit
        numit = numit + 1
        for i=1:nbrows
            deltax[i] = rhs[i]
            for j=1:nbcols
                deltax[i] = deltax[i] - mat[i,j]*result[j]
            end
            deltax[i] = deltax[i]/mat[i,i]
            result[i] = result[i] + deltax[i]
        end
        nrmdx = norm(deltax)
        strdx = @sprintf("%.2e", nrmdx)
        #println("||dx|| = $strdx")
        if norm(deltax) <= tol
            return result  #(result, numit, nrmdx, false)
        end
    end
    return result  #(result, numit, nrmdx, true)
end

function gauss_seidel2(mat::Array{Float64},sol::Array{Float64},rhs::Array{Float64},                      
                      maxit::Int=1000,tol::Float64=1.0e-8) 
    nbrows, nbcols = size(mat)
    result = deepcopy(sol)
    deltax = zeros(nbrows, 1)
    numit = 0; nrmdx = 1;
    while numit < maxit
        numit = numit + 1
        for i=1:nbrows
            deltax[i] = rhs[i]
            for j=1:nbcols
                deltax[i] = deltax[i] - mat[i,j]*result[j]
            end
            deltax[i] = deltax[i]/mat[i,i]
            result[i] = result[i] + deltax[i]
        end
        nrmdx = norm(deltax)
        strdx = @sprintf("%.2e", nrmdx)
        #println("||dx|| = $strdx")
        if norm(deltax) <= tol
            return result  #(result, numit, nrmdx, false)
        end
    end
    return result  #(result, numit, nrmdx, true)
end

#Gibbs for lambda version of MME (single-trait)
function Gibbs(A,x,b,vare,niter::Int64;printout_frequency=100)
    n = size(x,1)
    xmean = zeros(n)
    for iter = 1:niter
        if iter%printout_frequency==0
            println("at iteration: ",iter)
        end
        for i=1:n
            invlhs = 1.0/A[i,i]
            μ      = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]   = randn()*sqrt(invlhs*vare) + μ
        end
        xmean += (x - xmean)/iter
    end
    return xmean
end

#General Gibbs (multi-trait)
function Gibbs(A,x,b,niter::Int64;printout_frequency=100)
    n = size(x,1)
    xmean = zeros(n)
    for iter = 1:niter
        if iter%printout_frequency==0
            println("at iteration: ",iter)
        end
        for i=1:n
            if A[i,i] != 0.0 #issue70, zero diagonals in MME
                invlhs = 1.0/A[i,i]
                μ      = invlhs*(b[i] - A[:,i]'x) + x[i]
                x[i]   = randn()*sqrt(invlhs) + μ
            end
        end
        xmean += (x - xmean)/iter
    end
    return xmean
end

#One iteration of Gibbs sampler for lambda version of MME (single-trait)
function Gibbs(A,x,b,vare)
    for i = 1:size(x,1)
        if A[i,i] != 0.0 #issue70, zero diagonals in MME
            invlhs  = 1.0/A[i,i]
            μ       = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn()*sqrt(invlhs*vare) + μ
        end
    end
end

#One iteration of Gibbs sampler for general version of MME (multi-trait)
function Gibbs(A,x,b)
    for i = 1:size(x,1)
        if A[i,i] != 0.0 #issue70, zero diagonals in MME
            invlhs  = 1.0/A[i,i]
            μ       = invlhs*(b[i] - A[:,i]'x) + x[i]
            x[i]    = randn()*sqrt(invlhs) + μ
        end
    end
end
