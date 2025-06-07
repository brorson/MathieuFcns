function test_eigs_b()
  % This requests the value of b(q) for small q and then
  % tests the return against the power series given in the
  % DLMF: https://dlmf.nist.gov/28.6#i
    
  % Currently, the test is to simply plot my values against the
  % computed series value.  Later I can turn this into a real test.

  N = 30;
  
  Ne = 4;  % Number of eigenvalues to test.
  
  % The DLMF approximations seem to work over a limited domain.
  qs = linspace(0.001, 5.0, N)';
  
  q2 = qs.*qs;
  q4 = q2.*q2;
  q6 = q2.*q4;
  q8 = q4.*q4;

  % Expansions from the DLMF.
  b1 = 1 - qs - (1/8)*q2 + (1/64)*qs.*q2 - (1/1536)*q4 - (11/36864)*qs.* q4 ...
       + (49/589824)*q6 -(55/9437184)*qs.*q6 - (83/35389440)*q8;

  b2 = 4 - (1/12)*q2 + (5/13824)*q4 - (289/79626240)*q6 + (21391/458647142400)*q8;

  b3 = 9 + (1/16)*q2 - (1/64)*qs.*q2 + (13/20480)*q4 + (5/16384)*qs.*q4 ...
       - (1961/23592960)*q6 + (609/104857600)*qs.*q6;

  b4 = 16 + (1/30)*q2 - (317/86400)*q4 + (10049/2721600000)*q6;
  
    
  % Fill up calculated eigenvalue matrices.  Eigenvalue orders are
  % across rows, q values are down cols.
  b1_calc = zeros(N,1);
  b2_calc = zeros(N,1);
  b3_calc = zeros(N,1);   
  b4_calc = zeros(N,1);   
  for i=1:N
    q = qs(i);
    E = mathieu_b(Ne,q);  % 1 row, Ne cols
    %disp(size(E))
    b1_calc(i) = E(1,1);
    b2_calc(i) = E(1,2);
    b3_calc(i) = E(1,3);    
    b4_calc(i) = E(1,4);        
  end

  figure(1)
  plot(qs, b1,'-.')
  hold on
  plot(qs, b1_calc,'-')

  plot(qs, b2,'-.')
  plot(qs, b2_calc,'-')

  plot(qs, b3,'-.')
  plot(qs, b3_calc,'-')
  
  plot(qs, b4,'-.')
  plot(qs, b4_calc,'-')
  
  legend('Series b1','Calculated b1', 'Series b2', 'Calculated b2', ...
  	 'Series b3', 'Calculated b3', 'Series b4', 'Calculated b4', ...
  	 'Location','SouthWest')
  title('se eigenvalues')

  ylabel('Eigenvalues')
  xlabel('q')

end
