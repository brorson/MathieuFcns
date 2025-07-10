function plot_eigs()
  % Make plot of Mathieu eigs vs. q to reproduce plot
  % on https://dlmf.nist.gov/28.2
  
  % Number of sample points
  N = 151;
  
  % My playing field -- fcn domain.
  v = linspace(-pi, pi, N)';
  h = v(2)-v(1);

  % Domain of q values to examine (for plotting)
  qs = linspace(0,10,N)';
  

  %---------------------------------------------
  % First plot ce eigs
  % Number of even eigenvalues to track
  Ne = 5;

  % Preallocate a vector to store values.
  as = zeros(length(qs), Ne);
  
  % Loop over qs.
  fprintf('Calculating a eigenvalues ... ')
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


  %---------------------------------------------
  % Next plot se eigs
  
  % Number of sample points. For some reason I need
  % more points to get almost the same accuracy as for
  % the even eigenvalues.
  N = 251;

  % Number of odd eigenvalues to track
  Ne = 5;

  % Preallocate a vector to store values.
  bs = zeros(length(qs), Ne);
  
  % Loop over qs.
  fprintf('Calculating b eigenvalues ... ')  
  tic
  for i = 1:length(qs)
    q = qs(i);  % Get this value of q.
    
    % Get eigenvalues.  mathieu_b returns eigenvalues up
    % to Ne as a row vector.  Here I just stack up the row
    % vectors.
    bs(i,:) = mathieu_b(Ne, q);
  end
  toc

  fprintf('Odd eigs close to q=0:')
  disp(bs(1,:))

  % Make plot of eigenvalues vs. q
  %figure(1)
  for j=1:Ne
    hold on
    plot(qs,bs(:,j),'--')
    c = [c, ['b',num2str(j)]];
  end
  ylim([-5,20]);
  title('First Mathieu eigenvalues vs. q')
  xlabel('q')
  ylabel('eigenvalue')
  legend(c)

  
end