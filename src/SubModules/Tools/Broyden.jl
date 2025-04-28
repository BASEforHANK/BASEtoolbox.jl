@doc raw"""
    broyden(f, x0, critF, critX, maxiter, tau, tauS)

Solve for the root of a function `f` using the "good" Broyden algorithm, an iterative method
    to approximate the root of a nonlinear system of equations. This method is an
    approximation to Newton's method that uses updates to the Jacobian, improving efficiency
    in solving for the root.

# Arguments
- `f::Function`: The function for which the root is being found. It should return a vector
  of values.
- `x0::AbstractVector`: The starting guess for the root, provided as a column vector.
- `critF::Real`: The precision required for the function values to converge.
- `critX::Real`: The precision required for the change in `x` to converge.
- `maxiter::Int`: The maximum number of iterations allowed for the algorithm.
- `tau::Real`: A scaling factor used to adjust the update step.
- `tauS::Real`: A scaling factor for adjusting the step size dynamically with each
  iteration.

# Returns
- `xstar::AbstractVector`: The estimated root of the function.
- `fval::AbstractVector`: The value of the function at the root.
- `iter::Int`: The number of iterations performed.
- `distF::AbstractVector`: The distance from zero at each iteration, representing the
  convergence of the function.
"""
function broyden(f, x0, critF, critX, maxiter, tau, tauS)
    distF = zeros(maxiter)
    distF[1] = 1
    distX = 9999
    iter = 0
    xnow = x0[:] # x needs to be a column vector
    Fnow = f(xnow)
    Fnow = Fnow[:]  # F needs to be a column vector
    Bnow = I + zeros(length(xnow), length(xnow))

    while distF[iter + 1] > critF && distX > critX && iter < maxiter - 1
        iter = iter + 1 # count number of iterations
        tauF = min(tau, iter * tauS)

        Fold = copy(Fnow) # Store last function values
        xold = copy(xnow) # Store last x guess
        xnow = xnow - tauF .* Bnow * Fnow # Update x guess
        Fnow = f(xnow) # Evaluate the function
        Fnow = Fnow[:] # If Function is matrix valued then vectorize
        Dx = xnow - xold # Change in x
        Dy = Fnow - Fold # Change in f(x)

        # update inverse Jacobian
        Bnow .= Bnow .+ (Dx - Bnow * Dy) * (Dx' * Bnow) / (Dx' * Bnow * Dy)
        distF[iter + 1] = maximum(abs.(Fnow)) # Calculate the distance from zero
        distX = maximum(abs.(Dx)) # Calculate the change in variable
    end
    fval = Fnow
    xstar = xnow
    if distF[iter + 1] < critF
        display("Broyden has found a root of function f")
    else
        display("Failed to find a root of function f")
    end
    return xstar, fval, iter, distF
end
