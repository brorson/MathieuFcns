function write_mathieu_se_deriv_gvs()
  % This creates a file with golden values in columns used
  % to test other impls of the Mathieu se derivs.

  fid = stdin();
  q = fscanf(fid,'%f');
  
  fprintf('q = %f\n', q)
   
  N = 2500;   % Number of v values

  % mathieu_se only operates over the domain [-pi, pi]
  v = linspace(-pi,pi,N);
  dv = v(2) - v(1);
  Ne = 5;    % Top order of fcn to request.

  % Compute fcn values
  Ss = mathieu_se(Ne,q,N);  % GVs for different orders are
                            % arranged in columns.

  % This impl of mathieu_se is zero on both ends.  Remove zero
  % at end of vector so I get the right answer when computing
  % the deriv.
  Ss = Ss(1:end-1,:);
  v = v(1:end-1);
			    
  % Make plots to check the fcns.
  if 0
    leg = {};
    for i=1:Ne
      se = Ss(:,i);   % Extract one order.
      plot(v,se)
      hold on
      leg = [leg, ['se ',num2str(i-1)]];      
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
      sed = Ssd(:,i);   % Extract one order.
      plot(v,sed)
      hold on
      leg = [leg, ['sed ',num2str(i-1)]];
    end
    xlabel('v')
    ylabel('se deriv')
    legend(leg)
    title('Golden values')
  end

			    
  % Write GVs to a file along with the v value.
  filename = ['mathieu_se_deriv_gvs_q',num2str(q),'.csv'];
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
