function bs = mathieu_b(Ne, q)
  % This uses a finite-difference approximation to
  % the Mathieu equation to create an eigenvalue
  % problem.  The solution to the eigenvalue problem
  % gives the odd Mathieu eigenvalues.  This fcn returns
  % the first Ne eigenvalues.  The return is a row vector
  % of ascending order.
  
  % Number of sample points.  I should make this depend
  % upon how many eigenvalues are requested.
  N = 251;
  
  % My playing field -- fcn domain.
  v = linspace(-pi, pi, N)';  % Set domain including zero end points.

  % Preallocate space for eigenvalues.
  bs = zeros(1,Ne);
  
  % Create finite difference matrix of Mathieu ODE.
  % This matrix doesn't include the zero end points.
  Ae = make_matrix_o(N-2,q,v);
  
  % Compute eigenvalues
  opts = struct();
  opts.maxit = 2000;
  %opts.disp = 1;
  opts.p = 50;
  opts.tol = 1e-10;
  %[S,D,flag] = eigs(Ae,2*Ne,'largestreal',opts);
  [S,D,flag] = eigs(Ae,2*Ne,'sm',opts);  
  DD = -diag(D);

  % Must sort eigenvalues.
  [DD,idx] = sort(DD);
  
  % Store away eigenvalues.  Must select only ever other one.
  for j=2:2:length(DD)
    bs(j/2) = DD(j);
  end

  
end
