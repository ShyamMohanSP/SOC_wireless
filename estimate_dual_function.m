function res = estimate_dual_function(parameters,N1,N2,N3,M,u)
    
    a_opt = parameters.A0;
    p_pre = parameters.p_data(1) - (parameters.p_data(2)-parameters.p_data(1))/(parameters.re_timesteps(2)-parameters.re_timesteps(1))*(parameters.re_timesteps(1)+parameters.deltaCalib);
    pdot_pre = parameters.pdot_data(1) - (parameters.pdot_data(2)-parameters.pdot_data(1))/(parameters.re_timesteps(2)-parameters.re_timesteps(1))*(parameters.re_timesteps(1)+parameters.deltaCalib);

    % fun1 = @(x,y) interp3(0:1/N1:1,0:1/N2:1,0:1/N3:1,reshape(u(:,1),[N1+1,N2+1,N3+1]),a_opt,x,(y-parameters.xi_u)/(parameters.xi_bar-parameters.xi_u),'makima')*normpdf(x,p_pre+pdot_pre*parameters.deltaCalib,sqrt(2*parameters.alpha*parameters.theta0*p_pre*(1-p_pre)*parameters.deltaCalib))*gampdf(y-parameters.xi_u,parameters.sigma,1);
    % int_fun = @(x,y) arrayfun(@(a,b) fun1(a,b),x,y);
    % 
    % xmin = 0;
    % xmax = 1;
    % ymin = parameters.xi_u;
    % ymax = Inf;
    % 
    % res = integral2(int_fun,xmin,xmax,ymin,ymax);
    dual_cost = zeros(1,M);
    for m = 1:1:M        
        r_path = p_pre + pdot_pre*parameters.deltaCalib + sqrt(2*parameters.alpha*parameters.theta0*p_pre*(1-p_pre))*sqrt(parameters.deltaCalib)*randn;
        if r_path < 0
            r_path = abs(r_path);
        end
        xi_path = gamrnd(parameters.sigma,1)/(parameters.xi_bar-parameters.xi_u);
        dual_cost(m) = interp3(0:1/N1:1,0:1/N2:1,0:1/N3:1,reshape(u(:,1),[N1+1,N2+1,N3+1]),a_opt,r_path,xi_path,'makima');
    end

    %res = mean(dual_cost);

    % uncertainty bars

    nBootstrap = 1000;
    nSamples = 0.9*M;
    bootstrapMeans = zeros(nBootstrap,1);

    for i=1:nBootstrap
        bootstrapSample = datasample(dual_cost,nSamples,2);
        bootstrapMeans(i) = mean(bootstrapSample,2);
    end

    confInterval = prctile(bootstrapMeans,[2.5 97.5]);
    res = [mean(bootstrapMeans),confInterval];

end