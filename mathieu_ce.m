function ce = mathieu_ce(Ne,q,N)
  % This uses a finite-difference approximation to
  % the Mathieu equation to create an eigenvalue
  % problem.  The Mathieu functions appear as sampled
  % functions in the col vectors of S.
  
  % My playing field -- fcn domain.
  v = linspace(-pi, pi, N)';
  h = v(2)-v(1);
  % Find location of v = 0.  Used to flip sign of some
  % eigenvectors.
  %zidx = find( abs(v) < (v(end)-v(1))/N );
  %zidx = zidx(1);
  zidx = floor(N/2);
  
  %----------------------------------------------------------

  A = make_matrix_e(N,q,v);
  %fprintf('cond(A) = %e\n', cond(full(A)))
  
  % Compute eigenvalues & vectors
  opts = struct();
  opts.maxit = 2000;
  %opts.disp = 1;
  %opts.p = 500;  % Was 75
  opts.p = 75;
  %opts.tol = 1e-8*eps();
  opts.v0 = sin(v).^2;

  % Use this when running Matlab
  %[S,D,flag] = eigs(A,2*Ne,'largestreal',opts);
  
  % Use this when running Octave
  %[S,D,flag] = eigs(A+abs(q)*eye(size(A)),3*Ne,'sm', opts);
  [S,D,flag] = eigs(A,2*Ne,'sm', opts);  

  % I need to check for convergence since I get bad
  % results for some values of q,
  if (flag ~= 0)
    error('eigs did not converge!')
  end
  
  % Try Lapack -- it takes forever to create all GVs
  %[S,D] = eig(A+abs(q)*eye(size(A)),'nobalance');  
  %[S,D] = eig(A,'nobalance');    
 
  % Since position of eigenvalues & vectors from eigs jumps around,
  % sort them prior to use.  Also extract even fcns.
  [DD, idx] = sort(diag(-D),'ascend');
  %fprintf('DD(1) = %18.15e, DD(2) = %18.15e \n', DD(1), DD(2))
  diff = DD(2)-DD(1);
  fprintf('diff = %e\n', diff)
  if (abs(diff)<1e-8)
    error('Eigenvalues too close to resolve!')
  end
  
  S = S(:,idx);

  % Now extract the odd cols and put them into ce
  % (the even cols correspond to solns we don't want).
  ce = S(:,1:2:2*Ne);
  
  % Correct sign of fcns.  By definition, all fcns are
  % positive at v = 0.
  for j=1:Ne
    if (ce(zidx,j) < 0)
      % Must flip sign.
      %fprintf('Flipping sign for j=%d\n', j)
      ce(:,j) = -ce(:,j);
    end
  end
  
  % Now normalize fcns.  Use DLMF normalization.
  % cf.  DLMF 28.4.13.
  for j=1:Ne
    %fprintf('---------------------\n')
    C2 = h*trapz(ce(:,j).^2)/pi;
    C = sqrt(C2);
    ce(:,j) = ce(:,j)/C;
    %fprintf('j = %d, C = %f\n', j, C)
  end
  
end
