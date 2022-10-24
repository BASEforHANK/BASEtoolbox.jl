function  broyden(f,x0,critF,critX,maxiter,tau,tauS)
  # # The function uses the "good" broyden algorithm to solve for the root of a function F.
  # # X0 is the starting guess. CRITX the precision in variable X, CRITF the precision
  # # to which the root needs to be found. MAXITER is the maximum number of iterations.
  distF = zeros(maxiter);
  distF[1] = 1;
  distX = 9999;
  iter = 0;
  xnow = x0[:]; # x needs to be a column vector
  Fnow = f(xnow);
  Fnow = Fnow[:];  # F needs to be a column vector
  Bnow = I + zeros(length(xnow),length(xnow));

  while distF[iter+1]>critF && distX>critX && iter<maxiter-1
    iter = iter+1; # count number of iterations

    #  tauF = 2.0 .* tau ./ (1.0 .+ exp(iter.*tauS)) .+ (1.0 .- 2.0 ./ (1.0 .+ exp(iter.*tauS)))
     tauF = min(tau,iter*tauS)


     Fold = copy(Fnow); # Store last function values
     xold = copy(xnow); # Store last x guess
     xnow = xnow - tauF.* Bnow*Fnow; # Update x guess
     Fnow = f(xnow); # Evaluate the function
     Fnow = Fnow[:]; # If Function is matrix valued then vectorize
     Dx = xnow - xold; # Change in x
     Dy = Fnow - Fold; # Change in f(x)
    # update inverse Jacobian
     Bnow .= Bnow .+ (Dx - Bnow*Dy)*(Dx'*Bnow)/(Dx'*Bnow*Dy);
     distF[iter+1] = maximum(abs.(Fnow)); # Calculate the distance from zero
     distX = maximum(abs.(Dx)); # Calculate the change in variable


    #  println( distF[iter+1])
  end
  fval=Fnow;
  xstar=xnow;
  if distF[iter+1]<critF;
    display("Broyden has found a root of function f")
  else
    display("Failed to find a root of function f")
  end
return xstar, fval, iter, distF
end
