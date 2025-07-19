function se = mathieu_se(Ne,q,N)
  % This uses a finite-difference approximation to
  % the Mathieu equation to create an eigenvalue
  % problem.  The Mathieu functions appear as sampled
  % functions in the col vectors of S.
  
  % My playing field -- fcn domain.  I include all points
  % in the domain.
  v = linspace(-pi, pi, N)';
  h = v(2)-v(1);

  % N must be >= 2*Ne+2. so check that here.
  if (N < (2*Ne+2))
    error(' N must be >= 2*Ne+2\n')
  end

  %----------------------------------------------------------
  % Compute se using collocation method.

  % Make matrix.  Ask for matrix with row/col lengths two less than
  % full input domain.  
  A = make_matrix_o(N-2,q,v);
  % fprintf('cond(A) = %e\n', cond(full(A)))
  
  % Compute eigenvalues & vectors
  opts = struct();
  opts.maxit = 2000;
  %opts.disp = 1;
  opts.p = 75;  % 200
  %opts.tol = 1e-15;
  opts.v0 = exp(-(1:(N-2))');  % May not be helpful...

  % Use this when running Matlab
  %[S,D,flag] = eigs(A,2*Ne,'largestreal',opts);
  
  % Use this when running Octave
  [S,D] = eigs(A,2*Ne,'sm', opts);

  % I need to check for convergence since I get bad
  % results for some values of q,
  if (flag ~= 0)
    error('eigs did not converge!')
  end
  
  % Since position of eigenvalues & vectors from eigs jumps around,
  % sort them prior to use.  Also extract odd fcns.
  [DD, idx] = sort(diag(-D),'ascend');
  %fprintf('DD(1) = %18.15e, DD(2) = %18.15e \n', DD(1), DD(2))
  diff = DD(2)-DD(1);
  fprintf('diff = %e\n', diff)
  if (abs(diff)<1e-8)
    error('Eigenvalues too close to resolve!')
  end
  
  S = S(:,idx);

  % Tack zeros to top and bottom of the eigenvector matrix.
  S = [zeros(1, 2*Ne); S; zeros(1, 2*Ne)];  

  % Now extract the even cols and put them into se
  se = S(:,2:2:2*Ne);
  
  % By defintion, all slopes are positive for v = 0.  Modify fcns
  % to obey this definition.
  zidx = floor(N/2);  % Location of v=0
  for j=1:Ne
    if ( (se(zidx+1,j)-se(zidx,j)) < 0 )
      % Slope was negative.
      %fprintf('Flipping sign, j = %d\n', j)
      se(:,j) = -se(:,j);
    end
  end
  
  % Now normalize fcns.  Use DLMF normalization.
  % cf.  DLMF 28.4.13.
  for j=1:Ne
    %fprintf('---------------------\n')
    C2 = h*trapz(se(:,j).^2)/pi;
    C = sqrt(C2);
    se(:,j) = se(:,j)/C;
    %fprintf('j = %d, C = %f\n', j, C)
  end
  
end
