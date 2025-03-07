classdef MyParameters
    properties
        Pareto
        z_BS
        kappa
        sigma0
        SNR_th
        eta
        C1
        C2
        phi_th
        epsilon
        Pbar_R
        alpha
        theta0
        deltaCalib
        sigma
        B
        Pbar_tx
        A_omega
        Nu_max
        std_users
        A0
        Abar
        Pu_A
        Pbar_A
        C_scal
        C_offset
        xi_u
        xi_bar
        re_timesteps
        p_data
        pdot_data
        kb_timesteps
        k_data
        lambda_ell
        lambda_discrete
        tau_discrete
        tau_grid
        tau_threshold
        user_dist
        cell_price
        P_k
        horizon
    end
    methods 
        function res = max_flux(obj,x)
            res = min(obj.Pbar_A*x*5/obj.Abar,obj.Pbar_A);
        end
        function res = min_flux(obj,x)
            res = max(obj.Pu_A*(x-obj.Abar)*5/obj.Abar,-obj.Pu_A);
        end
        function res = pdot_fn(obj,t)
            res = interp1(obj.re_timesteps(1:end-1),obj.pdot_data,t,'linear','extrap');
        end
        function res = p_fn(obj,t)
            res = interp1(obj.re_timesteps,obj.p_data,t);
        end
        function res = theta_fn(obj,t)
            res = max(obj.theta0,(obj.alpha*obj.theta0 + abs(obj.pdot_fn(t)))./min(obj.p_fn(t),1-obj.p_fn(t)));
        end
        function res = Nu(obj,t)
            res = max(100,obj.Nu_max*0.125*(1+sin(pi*t*24/6 + pi))^3);
 %           res = 1000;
        end
        function res = K_b(obj,t)
            p = polyfit(obj.kb_timesteps,obj.k_data,7);
            res = 1e-4*(p(1)*t.^7 + p(2)*t.^6 + p(3)*t.^5 + p(4)*t.^4 + p(5)*t.^3 + p(6)*t.^2 + p(7)*t + p(8));
        end
        function res = K_s(obj,t)
            p = polyfit(obj.kb_timesteps,obj.k_data,7);
            res = 1e-4*(p(1)*t.^7 + p(2)*t.^6 + p(3)*t.^5 + p(4)*t.^4 + p(5)*t.^3 + p(6)*t.^2 + p(7)*t + p(8));
        end
        function res = pi_fn(obj,t)
            res = obj.cell_price;
        end
        function res = F1_fn(obj,t,x,phi)
            res = (obj.Pbar_R*x - phi(3) - phi(4))/obj.Abar;
        end
        function res = F2_fn(obj,t,x)
            res = obj.pdot_fn(t) - obj.theta_fn(t)*(x - obj.p_fn(t));
        end
        function res = F3_fn(obj,x)
            res = obj.B*(obj.sigma/(obj.xi_bar-obj.xi_u) - x);
%            res = -obj.B*(x - (obj.sigma - obj.xi_u)/(obj.xi_bar-obj.xi_u));
        end
        function res = G1_fn(obj,x)
            res = obj.alpha*obj.theta0*x*(1-x);
        end
        function res = G2_fn(obj,x)
            res = obj.B*x/(obj.xi_bar-obj.xi_u);
%            res = obj.B*(obj.xi_u+x*(obj.xi_bar-obj.xi_u))/(obj.xi_bar-obj.xi_u)/(obj.xi_bar-obj.xi_u);
        end
        function res = lambda(obj,t)
            lambda_grid = 0:(obj.horizon/obj.lambda_ell):obj.horizon;
            for l = 1:1:obj.lambda_ell
                if t>=lambda_grid(l) && t<=lambda_grid(l+1)
                    res = obj.lambda_discrete(l);
                end
            end
        end
        % function res = tau_star(obj,t)
        %     res = interp1(obj.tau_grid,obj.tau_discrete,t);
        % end
        function res = user_cdf(obj,x,phi)
            x_scal = obj.xi_u + x*(obj.xi_bar-obj.xi_u);
            if obj.user_dist == 0
                res = pi*phi(2)*x_scal*obj.kappa/obj.A_omega/obj.sigma0/(10^(obj.SNR_th/10));
            elseif obj.user_dist == 1
                res = 1 - exp(-phi(2)*x_scal*obj.kappa/2/(obj.std_users^2)/obj.sigma0/(10^(obj.SNR_th/10)));
            end
        end
        function res = subgrad_integrand(obj,t,x,phi)
%            res = smooth_indicator(1-obj.user_cdf(x,phi),obj.phi_th) - obj.epsilon;
            res = heaviside(1-obj.user_cdf(x,phi)-obj.phi_th-eps) - obj.epsilon;
%            res = exp(log(obj.tau_star(t))+(1-obj.user_cdf(x,phi)-obj.phi_th)/obj.tau_star(t)) - obj.tau_star(t)*obj.epsilon;
        end
        function res = H_fn(obj,t,x,phi)
            res = obj.Pareto*obj.K_b(t)*phi(1) - obj.Pareto*obj.K_s(t)*phi(3) - obj.Pareto*obj.pi_fn(t)*obj.Nu(t)*obj.user_cdf(x,phi) ...
    + (1-obj.Pareto)*(obj.C1*phi(1) + obj.C2*phi(1)^2) + obj.lambda(t)*obj.subgrad_integrand(t,x,phi);
        end
    end
end