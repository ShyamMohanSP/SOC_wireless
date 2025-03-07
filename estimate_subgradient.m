function [res,mc_prob] = estimate_subgradient(parameters,Nt,N_tau,N1,N2,N3,M,delu)
    
    sg = zeros(parameters.lambda_ell,M);
    dt = 1/Nt;
    tgrid = 0:dt:parameters.horizon;
    
    a_opt = zeros(M,length(tgrid));
    a_opt(:,1) = parameters.A0*ones(M,1);
    r_path = zeros(M,length(tgrid));
    %obtain initial distribution
    p_pre = parameters.p_data(1) - (parameters.p_data(2)-parameters.p_data(1))/(parameters.re_timesteps(2)-parameters.re_timesteps(1))*(parameters.re_timesteps(1)+parameters.deltaCalib);
    pdot_pre = parameters.pdot_data(1) - (parameters.pdot_data(2)-parameters.pdot_data(1))/(parameters.re_timesteps(2)-parameters.re_timesteps(1))*(parameters.re_timesteps(1)+parameters.deltaCalib);

    for m=1:1:M
        r_path(m,1) = p_pre + pdot_pre*parameters.deltaCalib + sqrt(2*parameters.alpha*parameters.theta0*p_pre*(1-p_pre))*sqrt(parameters.deltaCalib)*randn;
        if r_path(m,1) < 0
            r_path(m,1) = abs(r_path(m,1));
        end
    end
    
    xi_path = zeros(M,length(tgrid));
    xi_path(:,1) = gamrnd(parameters.sigma,1,M,1)/(parameters.xi_bar-parameters.xi_u);
%    xi_path(:,1) = (gamrnd(parameters.sigma,1,M,1)-parameters.xi_u)/(parameters.xi_bar-parameters.xi_u);

    dWr = sqrt(dt)*randn(M,length(tgrid)-1);
    dWxi = sqrt(dt)*randn(M,length(tgrid)-1);
    outage = zeros(M,length(tgrid)-1);
    
    for m = 1:1:M

        control_path = zeros(4,length(tgrid)-1);
        
        for n=1:1:length(tgrid)-1
            
            r_path(m,n+1) = r_path(m,n) + parameters.F2_fn(tgrid(n),r_path(m,n))*dt + sqrt(2*parameters.G1_fn(r_path(m,n)))*dWr(m,n);

            if r_path(m,n+1) < 0
                r_path(m,n+1) = -r_path(m,n+1);
            elseif r_path(m,n+1) > 1
                r_path(m,n+1) = 2 - r_path(m,n+1);
            end

            xi_path(m,n+1) = xi_path(m,n) + parameters.F3_fn(xi_path(m,n))*dt + sqrt(2*parameters.G2_fn(xi_path(m,n)))*dWxi(m,n);

            if xi_path(m,n+1) < 0
                xi_path(m,n+1) = abs(xi_path(m,n+1));
            end

            delu_a_now = zeros(size(delu,1),1);

            for i=1:1:size(delu,1)
                delu_a_now(i) = interp1(0:parameters.horizon/N_tau:parameters.horizon,delu(i,:),tgrid(n));
            end
            
            delu_a_now_here = interp3(0:1/N1:1,0:1/N2:1,0:1/N3:1,reshape(delu_a_now,[N1+1,N2+1,N3+1]),a_opt(m,n),r_path(m,n),xi_path(m,n),'makima');

            controls = compute_controls(parameters,tgrid(n),[a_opt(m,n),r_path(m,n),xi_path(m,n)],delu_a_now_here);            
            control_path(:,n) = controls;
            
            outage(m,n) = 1 - parameters.user_cdf(xi_path(m,n),controls);
            a_opt(m,n+1) = a_opt(m,n) + parameters.F1_fn(tgrid(n),r_path(m,n),controls)*dt;

        end
        
        for l=1:1:parameters.lambda_ell
            for n=(1+((length(tgrid)-1)*(l-1)/parameters.lambda_ell)):1:((length(tgrid)-1)*l/parameters.lambda_ell)
                sg(l,m) = sg(l,m) + dt*parameters.subgrad_integrand(tgrid(n),xi_path(m,n),control_path(:,n));
            end
        end

    end

    res = mean(sg,2);

    if nargout > 1
        mc_prob = mean(heaviside(outage - parameters.phi_th - eps),1); 
    end
    
end

% figure
% plot(24*tgrid,pars.Abar*mean(a_opt),'LineWidth',1);
% grid on;
% hold on;
% plot(24*tgrid,pars.Abar*(mean(a_opt) + 1.96*std(a_opt)),'--');
% hold on;
% plot(24*tgrid,pars.Abar*(mean(a_opt) - 1.96*std(a_opt)),'--');
% axis([0 24 0 pars.Abar]);
% xlabel('Time');
% ylabel('A (in W)');
% title('Charge in battery (optimal path)');
% legend('Expectation','95% CI')
% 
% figure
% plot(24*pars.re_timesteps,pars.Pbar_R*pars.p_data,'LineWidth',1);
% grid on;
% hold on;
% plot(24*tgrid,pars.Pbar_R*mean(r_path),'LineWidth',1);
% hold on;
% plot(24*tgrid,pars.Pbar_R*mean(r_path)+1.96*pars.Pbar_R*std(r_path),'--');
% hold on;
% plot(24*tgrid,pars.Pbar_R*mean(r_path)-1.96*pars.Pbar_R*std(r_path),'--');
% axis([0 24 0 pars.Pbar_R]);
% xlabel('Time');
% ylabel('R (in W)');
% legend('Forecast','Expectation','95% CI');
% title('Wind power');
% 
% figure
% plot(24*tgrid,pars.xi_bar*mean(xi_path),'LineWidth',3);
% grid on;
% hold on;
% for m=1:1:M
%     plot(24*tgrid,pars.xi_bar*xi_path(m,:),'--');
%     hold on;
% end
% axis([0 24 0 pars.xi_bar]);
% xlabel('Time');
% ylabel('\xi');
% legend('Expectation','Sample paths')
% title('Wireless channel fading');
% 
% figure
% plot(24*tgrid,pars.Nu(tgrid),'LineWidth',1);
% grid on;
% axis([0 24 0 pars.Nu_max]);
% xlabel('Time');
% ylabel('N_U');
% title('Number of cellular users');
% 
% figure
% plot(24*tgrid(1:end-1),control_path(1,:),'LineWidth',1);
% hold on;
% plot(24*tgrid(1:end-1),pars.C_scal*pars.Nu(tgrid(1:end-1)).*control_path(2,:) + pars.C_offset,'LineWidth',1);
% hold on;
% plot(24*tgrid(1:end-1),control_path(3,:),'LineWidth',1);
% hold on;
% plot(24*tgrid(1:end-1),control_path(4,:),'LineWidth',1);
% grid on;
% axis([0 24 0 200]);
% xlabel('Time');
% ylabel('Power (in W)');
% legend('FFS','Net consumption','Sold','Battery');
% title('Optimal controls');
% 
% % figure
% % plot(0,0,'x','LineWidth',5);
% % hold on;
% % grid on;
% % for n=1:1:Nt+1
% %     user_locs = -0.5 + rand(2,pars.Nu(t));
% %     plot(user_locs(1,:),user_locs(2,:),'o');
% %     axis([-0.5 0.5 -0.5 0.5])
% %     title('User distribution');
% %     F(n) = getframe;
% % end
% % legend('Base Station','Users');
% 
% % figure
% % plot(24*tgrid(1:end-1),1 - min(1,pi*kappa/A_omega/sigma0/(10^(SNR_th/10))*control_path(2,:).*xi_path(M,1:end-1)),'LineWidth',1);
% % grid on;
% % hold on;
% % plot(24*tgrid(1:end-1),phi_th*ones(1,200),'--','LineWidth',1);
% % axis([0 24 0 1]);
% % xlabel('Time');
% % ylabel('Users (in percentage)');
% % title('Proportion of users in outage');
% 
% % u_grid = reshape(u(:,1),[6,6,6]);
% 
% % [X,Y] = meshgrid(Tau_a,Tau_r);
% % 
% % % figure
% % % surf(X,Y,u_grid(:,:,end));
% % % xlabel('a');
% % % ylabel('r');
% % % zlabel('u');
% % % title('Value function at \xi = 1 at T=0');
% % 
% % for n=1:1:Nt+1
% %     u_grid = reshape(u(:,n),[6,6,6]);
% %     surf(X,Y,u_grid(:,:,2));
% % %    axis([0 1 0 1 0 1])
% %     xlabel('a');
% %     ylabel('r');
% %     zlabel('u');
% %     axis([0 1 0 1 min(u,[],"all") max(u,[],"all")]);
% %     title('Value function at \xi = 0.2');
% %     F(n) = getframe;
% % end
% % 
% % for n=1:1:Nt+1
% %     PF_grid = reshape(PF(:,n),[6,6,6]);
% %     surf(X,Y,PF_grid(:,:,1));
% % %    axis([0 1 0 1 0 1])
% %     xlabel('a');
% %     ylabel('r');
% %     zlabel('u');
% %     axis([0 1 0 1 min(PF,[],"all") max(PF,[],"all")]);
% %     title('FFS Power bought at \xi = 0');
% %     F(n) = getframe;
% % end
% % 
% % for n=1:1:Nt+1
% %     Ptx_grid = reshape(Ptx(:,n),[6,6,6]);
% %     surf(X,Y,Ptx_grid(:,:,1));
% % %    axis([0 1 0 1 0 1])
% %     xlabel('a');
% %     ylabel('r');
% %     zlabel('u');
% %     axis([0 1 0 1 min(Ptx,[],"all") max(Ptx,[],"all")]);
% %     title('Transmitted power at \xi = 0');
% %     F(n) = getframe;
% % end
% % 
% % for n=1:1:Nt+1
% %     PS_grid = reshape(PS(:,n),[6,6,6]);
% %     surf(X,Y,PS_grid(:,:,1));
% % %    axis([0 1 0 1 0 1])
% %     xlabel('a');
% %     ylabel('r');
% %     zlabel('u');
% %     axis([0 1 0 1 min(PS,[],"all") max(PS,[],"all")]);
% %     title('Sold power at \xi = 0');
% %     F(n) = getframe;
% % end
% % 
% % for n=1:1:Nt+1
% %     PA_grid = reshape(PA(:,n),[6,6,6]);
% %     surf(X,Y,PA_grid(:,:,1));
% % %    axis([0 1 0 1 0 2])
% %     xlabel('a');
% %     ylabel('r');
% %     zlabel('u');
% %     axis([0 1 0 1 min(PA,[],"all") max(PA,[],"all")]);
% %     title('Used RE power at \xi = 0');
% %     F(n) = getframe;
% % end
% % Nt = 100;
% % 
% % u_plot_ta = zeros(N1+1,Nt+1);
% % for i=1:1:N1+1
% %     u_plot_ta(i,:) = value_fxn(convert_coord(i,N2+1,N3+1,N1,N2,N3),:);
% % end
% % 
% % u_plot_tr = zeros(N2+1,Nt+1);
% % for i=1:1:N2+1
% %     u_plot_tr(i,:) = value_fxn(convert_coord(N1+1,i,N3+1,N1,N2,N3),:);
% % end
% % 
% % u_plot_txi = zeros(N3+1,Nt+1);
% % for i=1:1:N3+1
% %     u_plot_txi(i,:) = value_fxn(convert_coord(N1+1,N2+1,i,N1,N2,N3),:);
% % end
% % 
% % [X,Y] = meshgrid(0:T/Nt:1,0:1/N1:1);
% % figure
% % surf(X,Y,u_plot_ta);
% % xlabel('t');
% % ylabel('a');
% % zlabel('u');
% % 
% % [X,Y] = meshgrid(0:T/Nt:1,0:1/N2:1);
% % figure
% % surf(X,Y,u_plot_tr);
% % xlabel('t');
% % ylabel('r');
% % zlabel('u');
% % 
% % [X,Y] = meshgrid(0:T/Nt:1,0:1/N3:1);
% % figure
% % surf(X,Y,u_plot_txi);
% % xlabel('t');
% % ylabel('\xi');
% % zlabel('u');