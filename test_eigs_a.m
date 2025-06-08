function test_eigs_a()
  % This requests the value of a(q) for small q and then
  % tests the return against the power series given in the
  % DLMF: https://dlmf.nist.gov/28.6#i
    
  % Currently, the test is to simply plot my values against the
  % computed series value.  Later I can turn this into a real test.

  N = 30;
  
  Ne = 4;  % Number of eigenvalues to test.
  
  % The DLMF approximations seem to work over a limited domain.
  qs = linspace(0.001, 2.5, N)';
  
  q2 = qs.*qs;
  q4 = q2.*q2;
  q6 = q2.*q4;
  q8 = q4.*q4;

  % Expansions from the DLMF.
  a0 = -(1/2)*q2 + (7/128)*q4 - (29/2304)*q6 + (68687/18874368)*q8;

  a1 = (1 + qs - (1/8)*q2 - (1/64)*qs.*q2 - (1/1536)*q4 + (11/36864)*q4.*qs ...
         + (49/589824)*q6 + (55/9437187)*qs.*q6 - (83/35389440)*q8);

  a2 = (4 + (5/12)*q2 - (763/13824)*q4 + (1002401/79626240)*q6 ...
	 - (1669068401/458647142400)*q8);

  a3 = 9 + (1/16)*q2 + (1/64)*qs.*q2 + (13/20480)*q4 - (5/16382)*qs.*q4 ...
       - (1961/23592960)*q6 - (609/104857600)*qs.*q6;
  
    
  % Fill up calculated eigenvalue matrices.  Eigenvalue orders are
  % across rows, q values are down cols.
  a0_calc = zeros(N,1);
  a1_calc = zeros(N,1);
  a2_calc = zeros(N,1);   
  a3_calc = zeros(N,1);   
  for i=1:N
    q = qs(i);
    E = mathieu_a(Ne,q);  % 1 row, Ne cols
    %disp(size(E))
    a0_calc(i) = E(1,1);
    a1_calc(i) = E(1,2);
    a2_calc(i) = E(1,3);    
    a3_calc(i) = E(1,4);        
  end

  figure(1)
  plot(qs, a0,'-.')
  hold on
  plot(qs, a0_calc,'-')

  plot(qs, a1,'-.')
  plot(qs, a1_calc,'-')

  plot(qs, a2,'-.')
  plot(qs, a2_calc,'-')
  
  plot(qs, a3,'-.')
  plot(qs, a3_calc,'-')
  
  legend('Series a0','Calculated a0', 'Series a1', 'Calculated a1', ...
	 'Series a2', 'Calculated a2', 'Series a3', 'Calculated a3', ...
	 'Location','NorthWest')
  title('Computed ce eigenvalues compared to power series')

  ylabel('Eigenvalues')
  xlabel('q')

  %------------------------------------------------------------
  % Test asymptotic expansion from section 3.3.2 in Bricombe's paper
  % https://arxiv.org/pdf/2008.01812v2.
  % SDB says: These don't work at all.  Why?
  return
  
  qs = linspace(20, 100, N)';
  h = sqrt(qs);
  a = zeros(N,4);  % Put eigenvalues into cols
  for m=0:3
    s = 2*m+1;
    a(:,m+1) = -2*h.*h + 2*s*h - (s*s+1)/8 - (1./(2^7*h))*(s^3+3*s) ...
	     - (1./(2^12*h.*h))*(5*s^4 + 34*s^2 + 9);
  end
  
  figure(2)
  plot(qs, a(:,1),'-.')
  hold on
  plot(qs, a0_calc,'-')

  plot(qs, a(:,2),'-.')
  plot(qs, a1_calc,'-')

  plot(qs, a(:,3),'-.')
  plot(qs, a2_calc,'-')
  
  plot(qs, a(:,4),'-.')
  plot(qs, a3_calc,'-')
  
  legend('Ince a0', 'Calculated a0', 'Ince a1', 'Calculated a1', ...
	 'Ince a2', 'Calculated a2', 'Ince a3', 'Calculated a3', ...
	 'Location','SouthWest')
  title('Computed ce eigenvalues compared to Ince')

  ylabel('Eigenvalues')
  xlabel('q')
  
  
end
