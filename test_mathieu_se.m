function test_mathieu_se()
  % This computes the residual of my computed Mathiue fcn
  % when inserted into a finite diff version of the original ODE.
   
  N = 250;   % Number of grid points.
  q = 1;
  v = linspace(-pi, pi, N)';  % Domain of interest.
  h = v(2) - v(1);
  Ne = 4;    % Top order of fcn to request.
  
  % First ask for Mathieu eigs and fcns.  Ask for first four orders.
  bs = mathieu_b(Ne, q);
  Ss = mathieu_se(Ne,q,N);
  
  % Now test each fcn individually
  for j=1:Ne
    b = bs(j);
    S = Ss(:,j);
    
    % Finite diff version of ODE
    r = zeros(size(v));
    %r(1) = (S(2) - 2*S(1) + S(end))/(h*h) + (b-2*q*cos(2*v(1))).*S(1);
    r(2:end-1) = (S(3:end) - 2*S(2:end-1) + S(1:end-2))/(h*h) + (b-2*q*cos(2*v(2:end-1))).*S(2:end-1);
    %r(end) = (S(1) - 2*S(end) + S(end-1))/(h*h) + (b-2*q*cos(2*v(end))).*S(end);
  
    figure(1)
    plot(v,r)
    hold on
    
  end
  title('Residual of se vs. v')
  legend('se1', 'se2', 'se3', 'se4')
  xlabel('v')
  ylabel('se_n')
  xlim([0, pi/2])
  
end

  
  