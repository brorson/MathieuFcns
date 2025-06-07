function A = make_matrix_o(N, q, vn)
  % This creates the finite-difference formula for
  % the Mathieu equation.
  % N = row/col length of matrix.  Doesn't include zero ends.
  % q = frequency parameter
  % vn = domain.  Includes zero ends.
  
  h = vn(2)-vn(1);
  
  c = 2*q*cos(2*vn(2:end-1));
  hm2 = ones(N,1)./(h*h);
  
  A = spdiags([hm2, -2*hm2-c, hm2],[-1,0,1], N, N);

  % Boundary conds: None here.  Zeros at ends are
  % set in caller.

end
