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
  % fprintf('cond(A) = %e\n', cond(full(A)))
  
  % Compute eigenvalues & vectors
  opts = struct();
  opts.maxit = 2000;
  %opts.disp = 1;
  %opts.p = 500;  % Was 75
  opts.p = 75;  % Was 75  
  %opts.tol = 1e-15;
  opts.v0 = exp(-(1:N)');

  % Use this when running Matlab
  %[S,D,flag] = eigs(A,2*Ne,'largestreal',opts);
  
  % Use this when running Octave
  [S,D] = eigs(A,2*Ne,'sm', opts);
 
  % Since position of eigenvalues & vectors from eigs jumps around,
  % sort them prior to use.  Also extract even fcns.
  [DD, idx] = sort(diag(-D),'ascend');
  S = S(:,idx);

  % Now extract the odd cols and put them into ce
  % (the even cols correspond to solns we don't want).
  ce = S(:,1:2:2*Ne);
  
  % Correct sign of fcns.  By definition, all fcns are
  % positive at v = 0.
  for j=1:Ne
    if (ce(zidx+5,j) < 0)
      % Must flip sign.
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
