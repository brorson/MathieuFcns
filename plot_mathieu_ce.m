function plot_mathieu_ce()
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

  % Domain of q values to examine (for plotting)
  qs = linspace(0,10,N)';
  
  % Preallocate a vector to store values.
  as = zeros(length(qs), Ne);
  
  % Loop over qs.
  tic
  for i = 1:length(qs)
    q = qs(i);  % Get this value of q.
    
    % Get eigenvalues.  mathieu_a returns eigenvalues up
    % to Ne as a row vector.  Here I just stack up the row
    % vectors.
    as(i,:) = mathieu_a(Ne, q);
  end
  toc

  fprintf('Even eigs close to q=0:')
  disp(as(1,:))
  
  % Make plot of eigenvalues vs. q
  figure(1)
  c = {};
  for j=1:Ne
    hold on
    plot(qs,as(:,j),'-')
    c = [c, ['a',num2str(j-1)]];
  end
  ylim([-5,20]);
  title('First Mathieu even eigenvalues vs. q')
  xlabel('q')
  ylabel('eigenvalue')
  legend(c)

%disp(DD)
%return
  
  %----------------------------------------------------------
  % Now compute and plot Mathieu functions for fixed q
  Ne = 8;
  q = 1.0;

  ce = mathieu_ce(Ne,q,v);

  % Now make plots
  cee_leg = {};
  ceo_leg = {};  
  for j=1:Ne
    if (mod(j-1,2) == 0)
      fprintf('j-1 = %d -- even ce\n', j-1)
      figure(2)
      plot(v,ce(:,j),'-')
      title('ce functions -- even j')
      cee_leg = [cee_leg, ['ce',num2str((j-1))]];    
    else
      fprintf('j-1 = %d -- odd ce\n', j-1)
      figure(3)
      plot(v,ce(:,j),'-')
      title('ce functions -- odd j')
      ceo_leg = [ceo_leg, ['ce',num2str((j-1))]];    
    end
    xlim([0,pi/2])
    hold on
  end
  
  if (mod(j-1,2) == 0)
    figure(2)
    legend(cee_leg);
  else
    figure(3)    
    legend(ceo_leg);
  end
  
  
end
