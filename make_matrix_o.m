function A = make_matrix_o(N, q, vn)
  % This creates the finite-difference formula for
  % the Mathieu equation.
  
  h = vn(2)-vn(1);
  
  c = 2*q*cos(2*vn);
  hm2 = ones(N,1)./(h*h);
  
  A = spdiags([hm2, -2*hm2-c, hm2],[-1,0,1], N, N);

  % Boundary conds: None here.  Zeros at ends are
  % set in caller.

end
