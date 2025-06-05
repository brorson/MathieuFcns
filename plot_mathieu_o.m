function plot_mathieu_o()
  % This uses a finite-difference approximation to
  % the Mathieu equation to create an eigenvalue
  % problem.  The solution to the eigenvalue problem
  % gives the Mathieu eigenvalues.  The Mathieu
  % functions appear as sampled functions in the eigen-
  % vectors.
  
  % Number of sample points
  N = 151;
  
  % Number of eigenvalues to track
  Ne = 5;

  % My playing field -- fcn domain.
  v = linspace(-pi, pi, N)';
  h = v(2)-v(1);
  % Find location of v = 0.  Used to flip sign of some
  % eigenvectors.
  zidx = find( abs(v) < (v(end)-v(1))/N );
  zidx = zidx(1);

  %----------------------------------------------------------
  % First plot Mathieu eigs vs. q to reproduce plot
  % on https://dlmf.nist.gov/28.2

  % Domain of q values to examine
  qs = linspace(0,10,N)';
  
  % Preallocate a vector to store values.
  bs = zeros(length(qs), Ne);
  
  % Loop over qs.
  tic
  for i = 1:length(qs)
    q = qs(i);  % Get this value of q.
    
    % Get eigenvalues.  mathieu_a returns eigenvalues up
    % to Ne as a row vector.  Here I just stack up the row
    % vectors.
    bs(i,:) = mathieu_b(Ne, q);
  end
  toc

  fprintf('Odd eigs close to q=0:')
  disp(bs(1,:))
  
  % Make plot of eigenvalues vs. q
  figure(1)
  c = {};
  for j=1:Ne
    hold on
    plot(qs,bs(:,j),'-')
    c = [c, ['b',num2str(j-1)]];
  end
  ylim([-5,20]);
  title('First Mathieu odd eigenvalues vs. q')
  xlabel('q')
  ylabel('eigenvalue')
  legend(c)

%disp(DD)
%return
  
  %----------------------------------------------------------
  % Now compute and plot Mathieu functions for fixed q
  Ne = 8;
  q = 1.0;

  se = mathieu_se(Ne,q,v);

  % Now make plots
  see_leg = {};
  seo_leg = {};  
  for j=1:Ne
    if (mod(j,2) == 0)
      fprintf('j = %d -- even se\n', j)
      figure(2)
      plot(v,se(:,j),'-')
      title('se functions -- even j')
      see_leg = [see_leg, ['se',num2str(j)] ];
    else
      fprintf('j = %d -- odd se\n', j)
      figure(3)
      plot(v,se(:,j),'-')
      title('se functions -- odd j')
      seo_leg = [seo_leg, ['se',num2str(j)] ];    
    end
    xlim([0,pi/2])
    hold on
  end
  
  for j=1:Ne  
    if (mod(j,2) == 0)
      figure(2)
      legend(see_leg);
    else
      figure(3)
      legend(seo_leg);
    end
  end
  
end
