function as = mathieu_a(Ne, q)
  % This uses a finite-difference approximation to
  % the Mathieu equation to create an eigenvalue
  % problem.  The solution to the eigenvalue problem
  % gives the even Mathieu eigenvalues.  This fcn returns
  % the first Ne eigenvalues.  The return is a row vector
  % of ascending order.
  
  % Number of sample points.  I should make this depend
  % upon how many eigenvalues are requested.
  N = 251;
  
  % My playing field -- fcn domain.
  v = linspace(-pi, pi, N)';

  % Preallocate space for eigenvalues.
  as = zeros(1,Ne);
  
  % Create finite difference matrix of Mathieu ODE.
  Ae = make_matrix_e(N,q,v);
  
  % Compute eigenvalues
  opts = struct();
  opts.maxit = 2000;
  %opts.disp = 1;
  opts.p = 50;
  opts.tol = 1e-10;
  [S,D,flag] = eigs(Ae,2*Ne,'largestreal',opts);
  DD = -diag(D);

  % Store away eigenvalues.  Must select only ever other one.
  for j=1:2:length(DD)
    as((j+1)/2) = DD(j);
  end

  
end
