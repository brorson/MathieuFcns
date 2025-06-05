function se = mathieu_se(Ne,q,v)
  % This uses a finite-difference approximation to
  % the Mathieu equation to create an eigenvalue
  % problem.  The Mathieu functions appear as sampled
  % functions in the col vectors of S.
  
  % Number of sample points
  N = length(v);
  
  % My playing field -- fcn domain.
  v = linspace(-pi, pi, N)';
  h = v(2)-v(1);
  % Find location of v = 0.  Used to flip sign of some
  % eigenvectors.
  zidx = find( abs(v) < (v(end)-v(1))/N );
  zidx = zidx(1);

  %----------------------------------------------------------
  % Compute se using collocation method.
  
  A = make_matrix_o(N-2,q,v(2:end-1));
  % fprintf('cond(A) = %e\n', cond(full(A)))
  
  % Compute eigenvalues & vectors
  opts = struct();
  opts.maxit = 2000;
  %opts.disp = 1;
  opts.p = 100;
  %opts.tol = 1e-15;
  [S,D,flag] = eigs(A,2*Ne,'largestreal',opts);
  
  % Since position of eigenvalues & vectors from eigs jumps around,
  % sort them prior to use.  Also extract odd fcns.
  [DD, idx] = sort(diag(-D),'ascend');
  % Tack zeros to top and bottom of the eigenvector matrix.
  S = [zeros(1, 2*Ne); S(:,idx); zeros(1, 2*Ne)];

  % Now extract the odd fcns and put them into se
  se = S(:,2:2:2*Ne);
  
  % Correct sign of fcns.  By definition, all fcns are
  % positive for v > 0.
  for j=1:Ne
    if (se(zidx+2,j) < -1e-3)
      % Must flip sign.
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
