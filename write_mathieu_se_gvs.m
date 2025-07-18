function write_mathieu_se_gvs()
  % This creates a file with golden values in columns used
  % to test other impls of the Mathieu se fcns.
   
  fid = stdin();
  q = fscanf(fid,'%f');
  
  fprintf('q = %f\n', q)
    
    
  N = 2500;   % Number of v values
  % mathieu_ce only operates over the domain [-pi, pi]
  v = linspace(-pi,pi,N);
  Ne = 35;    % Top order of fcn to request.

  % Compute fcn values
  Ss = mathieu_se(Ne,q,N);  % GVs for different orders are
                            % arranged in columns.
  
  % Make plots to check the fcns.
  if 0
    leg = {};
    for i=1:Ne
      se = Ss(:,i);   % Extract one order.
      plot(v,se)
      hold on
      leg = [leg, num2str(i-1)];
    end
    xlabel('v')
    ylabel('se')
    legend(leg)
    title('Golden values')
  end
  
  % Write GVs to a file along with the v value.
  filename = ['mathieu_se_gvs_q',num2str(q),'.csv'];
  fh = fopen(filename,'w');
  % First write q value to file.
  fprintf(fh,'%f\n',q)
  % Then write fcn values.
  fmt = ['%f, ',repmat('%f, ',[1,Ne-1]),'%f \n'];
  for i=1:length(v)
    fprintf(fh, fmt, v(i), Ss(i,:));
  end
  fclose(fh);
  
end

  
