function A = make_matrix_e(N, q, vn)
  % This creates the finite-difference formula for
  % the Mathieu equation.
  
  h = vn(2)-vn(1);
  
  c = 2*q*cos(2*vn);
  hm2 = ones(N,1)./(h*h);
  
  A = spdiags([hm2, -2*hm2-c, hm2],[-1,0,1], N, N);
  % Boundary conditions -- cyclic
  %A(1,end) = hm2(1); % A(1,2);
  %A(end,1) = hm2(1); % A(end,end-1);

  % Boundary conditions -- e
  % These are zero-slope BCs.
  A(1,2) = 2*hm2(1);
  A(end,end-1) = 2*hm2(1);

end
