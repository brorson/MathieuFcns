function write_mathieu_ce_gvs()
  % This creates a file with golden values in columns used
  % to test other impls of the Mathieu ce fcns.
   
  N = 1000;   % Number of v values
  q = 1;
  % mathieu_ce only operates over the domain [-pi, pi]
  v = linspace(-pi,pi,N);
  Ne = 12;    % Top order of fcn to request.

  % Compute fcn values
  Ss = mathieu_ce(Ne,q,N);  % GVs for different orders are
                            % arranged in columns.
  
  % Make plots to check the fcns.
  leg = {};
  for i=1:Ne
    ce = Ss(:,i);   % Extract one order.
    plot(v,ce)
    hold on
    leg = [leg, num2str(i-1)];
  end
  xlabel('v')
  ylabel('ce')
  legend(leg)
  title('Golden values')
  
  % Write GVs to a file along with the v value.
  fh = fopen('mathieu_ce_gvs.csv','w');
  fmt = ['%f, ',repmat('%f, ',[1,Ne-1]),'%f \n'];
  for i=1:length(v)
    fprintf(fh, fmt, v(i), Ss(i,:));
  end
  fclose(fh);
  
end

  
  