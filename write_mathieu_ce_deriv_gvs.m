function write_mathieu_ce_deriv_gvs()
  % This creates a file with golden values in columns used
  % to test other impls of the Mathieu ce derivs.
   
  fid = stdin();
  q = fscanf(fid,'%f');
  
  fprintf('q = %f\n', q)
 
  N = 2500;   % Number of v values
  % mathieu_ce only operates over the domain [-pi, pi]
  v = linspace(-pi,pi,N);
  dv = v(2) - v(1);
  Ne = 5;    % Top order of fcn to request.

  % Compute fcn values
  Ss = mathieu_ce(Ne,q,N);  % GVs for different orders are
                            % arranged in columns.

  % Make plots to check the fcns.
  if 0
    leg = {};
    for i=1:Ne
      ce = Ss(:,i);   % Extract one order.
      plot(v,ce)
      hold on
      leg = [leg, ['ce ',num2str(i-1)]];      
    end
  end

  % Compute center difference derivs.  Since the fcn has cyclic
  % BCs I create a matrix with the same number of elements, and
  % and handle the domain ends separately.
  Ssd = zeros(size(Ss));
  Ssd(1,:) = (Ss(2,:)-Ss(end,:))/(2*dv);
  Ssd(2:end-1,:) = (Ss(3:end,:)-Ss(1:end-2,:))/(2*dv);
  Ssd(end,:) = (Ss(1,:)-Ss(end-1,:))/(2*dv);

  % Make plots to check the fcns.
  if 0
    for i=1:Ne
      ced = Ssd(:,i);   % Extract one order.
      plot(v,ced)
      hold on
      leg = [leg, ['ced ',num2str(i-1)]];
    end
    xlabel('v')
    ylabel('ce deriv')
    legend(leg)
    title('Golden values')
  end

  
  % Write GVs to a file along with the v value.
  filename = ['mathieu_ce_deriv_gvs_q',num2str(q),'.csv'];
  fh = fopen(filename,'w');
  % First write q value to file.
  fprintf(fh,'%f\n',q)
  % Then write fcn values.
  fmt = ['%f, ',repmat('%f, ',[1,Ne-1]),'%f \n'];
  for i=1:length(v)
    fprintf(fh, fmt, v(i), Ssd(i,:));
  end
  fclose(fh);
  
end

  
  