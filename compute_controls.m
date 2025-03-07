function P = compute_controls(parameters,t,x,delu)
    x3_scal = parameters.xi_u + x(3)*(parameters.xi_bar-parameters.xi_u);

    %first minimization problem
    d5 = -2*(parameters.std_users^2)*parameters.sigma0*10^(0.1*parameters.SNR_th)*log(parameters.phi_th)/x3_scal/parameters.kappa;
%    d5 = 0;
    if parameters.user_dist == 0
        a1 = delu/parameters.Abar + parameters.Pareto*parameters.K_b(t) + (1-parameters.Pareto)*parameters.C1;
        a2 = (1-parameters.Pareto)*parameters.C2;
        b1 = -(delu/parameters.Abar) - parameters.Pareto*parameters.K_s(t);
        c1 = -(parameters.C_scal*parameters.Nu(t)*delu/parameters.Abar) - parameters.Pareto*parameters.pi_fn(t)*parameters.Nu(t)*pi*parameters.kappa*x3_scal/parameters.A_omega/parameters.sigma0/(10^(parameters.SNR_th/10));
        c2 = 0;
        c3 = 0;
        d1 = parameters.C_scal*parameters.Nu(t);
        d2 = -parameters.max_flux(parameters.Abar*x(1)) - parameters.Pbar_R*x(2) + parameters.C_offset;
        d3 = parameters.Pbar_tx/parameters.Nu(t);
        d4 = parameters.min_flux(parameters.Abar*x(1)) + parameters.Pbar_R*x(2) - parameters.C_offset;
    elseif parameters.user_dist == 1
        a1 = delu/parameters.Abar + parameters.Pareto*parameters.K_b(t) + (1-parameters.Pareto)*parameters.C1;
        a2 = (1-parameters.Pareto)*parameters.C2;
        b1 = -(delu/parameters.Abar) - parameters.Pareto*parameters.K_s(t);
        c1 = -parameters.C_scal*parameters.Nu(t)*delu/parameters.Abar;
        c2 = parameters.Pareto*parameters.pi_fn(t)*parameters.Nu(t);
        c3 = x3_scal*parameters.kappa/2/(parameters.std_users^2)/parameters.sigma0/(10^(parameters.SNR_th/10));
        d1 = parameters.C_scal*parameters.Nu(t);
        d2 = -parameters.max_flux(parameters.Abar*x(1)) - parameters.Pbar_R*x(2) + parameters.C_offset;
        d3 = parameters.Pbar_tx/parameters.Nu(t);
        d4 = parameters.min_flux(parameters.Abar*x(1)) + parameters.Pbar_R*x(2) - parameters.C_offset;
    end
    
    fun1 = @(z) c1 - c2*c3*exp(-c3*z) + a1*d1 + 2*a2*d1*d1*z + 2*a2*d1*parameters.C_offset;
    fun2 = @(z) c1 - c2*c3*exp(-c3*z) + a1*d1 + 2*a2*d1*d1*z + 2*a2*d1*d2;
    fun3 = @(z) c1 - c2*c3*exp(-c3*z) + d1*(a1+2*a2*(d1*z-d4));

    if a1 <=0 && b1 >=0 && c1 <= c2*c3*exp(-c3*d5) && c1 >= c2*c3*exp(-c3*d3) && d2 <= -a1/2/a2 + d1*log(c1/c2/c3)/c3 && d4 <= a1/2/a2 - d1*log(c1/c2/c3)/c3 && a1/2/a2 >= d1*log(c1/c2/c3)/c3 - parameters.C_offset
        PF = -a1/2/a2;  %case1
        Ptx = -log(c1/c2/c3)/c3;
        PS = 0;

    elseif a1 >= 0 && b1 >=0 && c1 <= c2*c3*exp(-c3*d5) && c1 >= c2*c3*exp(-c3*d3) && d2 <= d1*log(c1/c2/c3)/c3 && d4 <= -d1*log(c1/c2/c3)/c3 && parameters.C_offset >= d1*log(c1/c2/c3)/c3
        PF = 0;%case3
        Ptx = -log(c1/c2/c3)/c3;
        PS = 0;

    elseif a1 <= 0 && b1 >= 0 && c1 >= c2*c3*exp(-c3*d5) && d2 <= -a1/2/a2 - d1*d5 && d4 <= a1/2/a2 + d1*d5 && parameters.C_offset >= -a1/2/a2 - d1*d5
        PF = -a1/2/a2;  %case4
        PS = 0;
        Ptx = d5;

    elseif -a1/2/a2 >= d1*d5 + parameters.C_offset && b1 >= 0 && c1 >= c2*c3*exp(-c3*d5) - d1*(a1+2*a2*(d1*d5+parameters.C_offset)) && d2 <= parameters.C_offset && d4 <= -parameters.C_offset
        PF = d1*d5 + parameters.C_offset;   %case5
        Ptx = d5;
        PS = 0;

    elseif a1 >= 0 && b1 >= 0 && c1 >= c2*c3*exp(-c3*d5) && d2 <= -d1*d5 && d4 <= d1*d5
        PF = 0;     %case6
        Ptx = d5;
        PS = 0;

    elseif a1 <= 0 && b1 >= 0 && c1 <= c2*c3*exp(-c3*d3) && d2 <= -a1/2/a2 - d1*d3 && d4 <= a1/2/a2 + d1*d3 && -a1/2/a2 <= d1*d3 + parameters.C_offset
        PF = -a1/2/a2;%case7
        Ptx = d3;
        PS = 0;

    elseif -a1/2/a2 >= d1*d3 + parameters.C_offset && b1 >= 0 && c1 <= c2*c3*exp(-c3*d3) - d1*(a1 + 2*a2*(d1*d3+parameters.C_offset)) && d2 <= parameters.C_offset && d4 <= -parameters.C_offset
        PF = d1*d3 + parameters.C_offset;   %case8
        Ptx = d3;
        PS = 0;

    elseif a1 >= 0 && b1 >= 0 && c1 <= c2*c3*exp(-c3*d3) && d2 <= -d1*d3 && d4 <= d1*d3 
        PF = 0;     %case9
        PS = 0;
        Ptx = d3;

    elseif b1 <= 0 && c1 <= c2*c3*exp(-c3*d5) + b1*d1 && c1 >= b1*d1 + c2*c3*exp(-c3*d3) && d2 <= d1*log((c1-b1*d1)/c2/c3)/c3
        PF = 0;     %case10
        Ptx = -log((c1-b1*d1)/c2/c3)/c3;
        PS = d1*log((c1-b1*d1)/c2/c3)/c3 - d2;

    elseif d2 <= -d1*d5 && b1 <=0 && c1>=c2*c3*exp(-c3*d5)+b1*d1
        PF = 0;     %case12
        Ptx = d5;
        PS = -d2 - d1*d5;

    elseif d2 <= -d1*d5 && c1 <= c2*c3*exp(c3*d2/d1) && c1 <= c2*c3*exp(c3*d2/d1) + b1*d1 && c1 >= c2*c3*exp(c3*d2/d1) - a1*d1 && d2 <= parameters.C_offset && d2 >= -d1*d3 
        PF = 0;     %case13
        PS = 0;
        Ptx = -d2/d1;

    elseif d2 >= -d1*d5 && -a1/2/a2 <= d2 + d1*d5 && c1 >= c2*c3*exp(-c3*d5) - d1*(a1+2*a2*(d2+d1*d5)) && d2 <= parameters.C_offset
        PF = d2 + d1*d5;    %case14
        PS = 0;
        Ptx = d5;

    elseif d2 >= -d1*d3 && -a1/2/a2 <= d1*d3+d2 && c1 <= c2*c3*exp(-c3*d3) - d1*(a1+2*a2*(d1*d3+d2)) && d2 <= parameters.C_offset
        PF = d1*d3+d2;      %case15
        PS = 0;
        Ptx = d3;

    elseif d2 <= -d1*d3 && b1 <= 0 && c1 <= c2*c3*exp(-c3*d3) + b1*d1
        PF = 0;           %case16
        Ptx = d3;
        PS = -d2-d1*d3;

    elseif b1 >= 0 && c1 <= c2*c3*exp(-c3*d5) + b1*d1 && c1 >= c2*c3*exp(-c3*d3) + b1*d1 && d4 >= -d1*log((c1-b1*d1)/c2/c3)/c3
        PF = 0;       %case17
        PS = d4 + d1*log((c1-b1*d1)/c2/c3)/c3;
        Ptx = -log((c1-b1*d1)/c2/c3)/c3;

    elseif d4 >= d1*d5 && b1 >= 0 && c1 >= c2*c3*exp(-c3*d5) + b1*d1
        PF = 0;     %case19
        Ptx = d5;
        PS = d4 - d1*d5;

    elseif d4 >= d1*d5 && d4 <= d1*d3 && c1 >= c2*c3*exp(-c3*d4/d1) && c1 >= c2*c3*exp(-c3*d4/d1) - a1*d1 && c1 <= b1*d1 + c2*c3*exp(-c3*d4/d1) && d4 >= -parameters.C_offset
        PF = 0;   %case20
        PS = 0;
        Ptx = d4/d1;

    elseif d4 <= d1*d5 && -a1/2/a2 >= d1*d5-d4 && c1 >= c2*c3*exp(-c3*d5) - d1*(a1+2*a2*(d1*d5-d4)) && d4 >= -parameters.C_offset
        PF = d1*d5-d4;       %case21
        Ptx = d5;
        PS = 0;

    elseif d4 >= d1*d3 && b1 >= 0 && c1 <= b1*d1 + c2*c3*exp(-c3*d3)
        PF = 0;       %case22
        Ptx = d3;
        PS = d4-d1*d3;

    elseif d4 <= d1*d3 && -a1/2/a2 >= d1*d3-d4 && c1 <= c2*c3*exp(-c3*d3)-d1*(a1+2*a2*(d1*d3-d4)) && d4 >= -parameters.C_offset
        Ptx = d3;     %case23
        PS = 0;
        PF = d1*d3-d4;

    else
        x1 = min([d3,-(a1/2/a2+parameters.C_offset)/d1]);
        x2 = max([d5,-d2/d1,-(a1/2/a2+d2)/d1]);
        x3 = max([d5,d4/d1]);
        x4 = min([d3,(-a1/2/a2+d4)/d1]);
        if b1 >= 0 && d2 <= parameters.C_offset && d4 <= -parameters.C_offset && x1 >= d5 && fun1(d5)*fun1(x1) < 0
            bounds = bisection(fun1,d5,x1,1e-10);
            Ptx = fzero(fun1,bounds);   %case2
            PF = d1*Ptx+parameters.C_offset;
            PS = 0;
        elseif d2 <= parameters.C_offset && x2 <= d3 && fun2(x2)*fun2(d3) < 0
            bounds = bisection(fun2,x2,d3,1e-10);
            Ptx = fzero(fun2,bounds);       %case11
            PF = d1*Ptx+d2;
            PS = 0;
        elseif d4 >= -parameters.C_offset && x3 <= x4 && fun3(x3)*fun3(x4) < 0
            bounds = bisection(fun3,x3,x4,1e-10);
            Ptx = fzero(fun3,bounds);       %case18
            PF = d1*Ptx-d4;
            PS = 0;
        end
    end
 
    PA = parameters.C_scal*parameters.Nu(t)*Ptx + parameters.C_offset - PF;
    fval1 = a1*PF + a2*PF*PF + b1*PS + c1*Ptx + c2*exp(-c3*Ptx);
    P1 = [PF,Ptx,PS,PA];

    %second minimization problem
    if parameters.user_dist == 0
        a1 = delu/parameters.Abar + parameters.Pareto*parameters.K_b(t) + (1-parameters.Pareto)*parameters.C1;
        a2 = (1-parameters.Pareto)*parameters.C2;
        b1 = -(delu/parameters.Abar) - parameters.Pareto*parameters.K_s(t);
        c1 = -(parameters.C_scal*parameters.Nu(t)*delu/parameters.Abar) - parameters.Pareto*parameters.pi_fn(t)*parameters.Nu(t)*pi*parameters.kappa*x3_scal/parameters.A_omega/parameters.sigma0/(10^(parameters.SNR_th/10));
        c2 = 0;
        c3 = 0;
        d1 = parameters.C_scal*parameters.Nu(t);
        d2 = -parameters.max_flux(parameters.Abar*x(1)) - parameters.Pbar_R*x(2) + parameters.C_offset;
        d3 = d5;
        d4 = parameters.min_flux(parameters.Abar*x(1)) + parameters.Pbar_R*x(2) - parameters.C_offset;
    elseif parameters.user_dist == 1
        a1 = delu/parameters.Abar + parameters.Pareto*parameters.K_b(t) + (1-parameters.Pareto)*parameters.C1;
        a2 = (1-parameters.Pareto)*parameters.C2;
        b1 = -(delu/parameters.Abar) - parameters.Pareto*parameters.K_s(t);
        c1 = -parameters.C_scal*parameters.Nu(t)*delu/parameters.Abar;
        c2 = parameters.Pareto*parameters.pi_fn(t)*parameters.Nu(t);
        c3 = x3_scal*parameters.kappa/2/(parameters.std_users^2)/parameters.sigma0/(10^(parameters.SNR_th/10));
        d1 = parameters.C_scal*parameters.Nu(t);
        d2 = -parameters.max_flux(parameters.Abar*x(1)) - parameters.Pbar_R*x(2) + parameters.C_offset;
        d3 = d5;
        d4 = parameters.min_flux(parameters.Abar*x(1)) + parameters.Pbar_R*x(2) - parameters.C_offset;
    end
    
    fun1 = @(z) c1 - c2*c3*exp(-c3*z) + a1*d1 + 2*a2*d1*d1*z + 2*a2*d1*parameters.C_offset;
    fun2 = @(z) c1 - c2*c3*exp(-c3*z) + a1*d1 + 2*a2*d1*d1*z + 2*a2*d1*d2;
    fun3 = @(z) c1 - c2*c3*exp(-c3*z) + d1*(a1+2*a2*(d1*z-d4));

    if a1 <=0 && b1 >= 0 && c1 <= c2*c3 && c1 >= c2*c3*exp(-c3*d3) && d2 <= -a1/2/a2 + d1*log(c1/c2/c3)/c3 && d4 <= a1/2/a2 - d1*log(c1/c2/c3)/c3 && a1/2/a2 >= d1*log(c1/c2/c3)/c3 - parameters.C_offset
        PF = -a1/2/a2;
        Ptx = -log(c1/c2/c3)/c3;
        PS = 0;

    elseif a1 >= 0 && b1 >=0 && c1 <= c2*c3 && c1 >= c2*c3*exp(-c3*d3) && d2 <= d1*log(c1/c2/c3)/c3 && d4 <= -d1*log(c1/c2/c3)/c3 && parameters.C_offset >= d1*log(c1/c2/c3)/c3
        PF = 0;
        Ptx = -log(c1/c2/c3)/c3;
        PS = 0;

    elseif a1 <=0 && b1 >= 0 && c1 >= c2*c3 && d2 <= -a1/2/a2 && d4 <= a1/2/a2 && parameters.C_offset >= -a1/2/a2
        PF = -a1/2/a2;
        PS = 0;
        Ptx = 0;

    elseif parameters.C_offset <= -a1/2/a2 && b1 >= 0 && c1 >= c2*c3 - d1*(a1+2*a2*parameters.C_offset) && d2 <= parameters.C_offset && d4 <= -parameters.C_offset
        PF = parameters.C_offset;
        Ptx = 0;
        PS = 0;

    elseif a1 >= 0 && b1 >= 0 && c1 >= c2*c3 && d2 <= 0 && d4 <= 0
        PF = 0;
        Ptx = 0;
        PS = 0;

    elseif a1 <= 0 && b1 >= 0 && c1 <= c2*c3*exp(-c3*d3) && d2 <= -a1/2/a2 - d1*d3 && d4 <= a1/2/a2 + d1*d3 && -a1/2/a2 <= d1*d3 + parameters.C_offset
        PF = -a1/2/a2;
        Ptx = d3;
        PS = 0;

    elseif -a1/2/a2 >= d1*d3 + parameters.C_offset && b1 >= 0 && c1 <= c2*c3*exp(-c3*d3) - d1*(a1 + 2*a2*(d1*d3+parameters.C_offset)) && d2 <= parameters.C_offset && d4 <= -parameters.C_offset
        PF = d1*d3 + parameters.C_offset;
        Ptx = d3;
        PS = 0;

    elseif a1 >= 0 && b1 >= 0 && c1 <= c2*c3*exp(-c3*d3) && d2 <= -d1*d3 && d4 <= d1*d3 
        PF = 0;
        PS = 0;
        Ptx = d3;

    elseif b1 <= 0 && c1 <= c2*c3 + b1*d1 && c1 >= b1*d1 + c2*c3*exp(-c3*d3) && d2 <= d1*log((c1-b1*d1)/c2/c3)/c3
        PF = 0;
        Ptx = -log((c1-b1*d1)/c2/c3)/c3;
        PS = d1*log((c1-b1*d1)/c2/c3)/c3 - d2;

    elseif d2 <= 0 && b1 <= 0 && c1 >= c2*c3+b1*d1
        PF = 0;
        Ptx = 0;
        PS = -d2;

    elseif d2 <= 0 && c1 <= c2*c3*exp(c3*d2/d1) && c1 <= c2*c3*exp(c3*d2/d1)+b1*d1 && c1 >= c2*c3*exp(c3*d2/d1)-a1*d1 && d2 <= parameters.C_offset && d2 >= -d1*d3 
        PF = 0;
        PS = 0;
        Ptx = -d2/d1;

    elseif d2 >= 0 && -a1/2/a2 <= d2 && c1 >= c2*c3-d1*(a1+2*a2*d2) && d2 <= parameters.C_offset
        PF = d2;
        PS = 0;
        Ptx = 0;

    elseif d2 >= -d1*d3 && -a1/2/a2 <= d1*d3+d2 && c1 <= c2*c3*exp(-c3*d3) - d1*(a1+2*a2*(d1*d3+d2)) && d2 <= parameters.C_offset
        PF = d1*d3+d2;
        PS = 0;
        Ptx = d3;

    elseif d2 <= -d1*d3 && b1 <= 0 && c1 <= c2*c3*exp(-c3*d3)+b1*d1
        PF = 0;
        Ptx = d3;
        PS = -d2-d1*d3;

    elseif b1 >= 0 && c1 <= c2*c3+b1*d1 && c1 >= c2*c3*exp(-c3*d3)+b1*d1 && d4 >= -d1*log((c1-b1*d1)/c2/c3)/c3
        PF = 0;
        PS = d4 + d1*log((c1-b1*d1)/c2/c3)/c3;
        Ptx = -log((c1-b1*d1)/c2/c3)/c3;

    elseif d4 >= 0 && b1 >= 0 && c1 >= c2*c3 + b1*d1
        PF = 0;
        Ptx = 0;
        PS = d4;

    elseif d4 >= 0 && d4 <= d1*d3 && c1 >= c2*c3*exp(-c3*d4/d1) && c1 >= c2*c3*exp(-c3*d4/d1)-a1*d1 && c1 <= b1*d1 + c2*c3*exp(-c3*d4/d1) && d4 >= -parameters.C_offset
        PF = 0;
        PS = 0;
        Ptx = d4/d1;

    elseif d4 <= 0 && d4 >= a1/2/a2 && c1 >= c2*c3 + 2*a2*d1*d4 - a1*d1 && d4 >= -parameters.C_offset
        PF = -d4;
        Ptx = 0;
        PS = 0;

    elseif d4 >= d1*d3 && b1 >= 0 && c1 <= b1*d1 + c2*c3*exp(-c3*d3)
        PF = 0;
        Ptx = d3;
        PS = d4 - d1*d3;

    elseif d4 <= d1*d3 && -a1/2/a2 >= d1*d3-d4 && c1 <= c2*c3*exp(-c3*d3) - d1*(a1+2*a2*(d1*d3-d4)) && d4 >= -parameters.C_offset
        Ptx = d3;
        PS = 0;
        PF = d1*d3-d4;

    else
        x1 = min([d3,-(a1/2/a2+parameters.C_offset)/d1]);
        x2 = max([0,-d2/d1,-(a1/2/a2+d2)/d1]);
        x3 = max([0,d4/d1]);
        x4 = min([d3,(-a1/2/a2+d4)/d1]);
        if b1 >= 0 && d2 <= parameters.C_offset && d4 <= -parameters.C_offset && x1 >= 0 && fun1(0)*fun1(x1) < 0
            bounds = bisection(fun1,0,x1,1e-10);
            Ptx = fzero(fun1,bounds);
            PF = d1*Ptx+parameters.C_offset;
            PS = 0;
        elseif d2 <= parameters.C_offset && x2 <= d3 && fun2(x2)*fun2(d3) < 0
            bounds = bisection(fun2,x2,d3,1e-10);
            Ptx = fzero(fun2,bounds);
            PF = d1*Ptx+d2;
            PS = 0;
        elseif d4 >= -parameters.C_offset && x3 <= x4 && fun3(x3)*fun3(x4) < 0
            bounds = bisection(fun3,x3,x4,1e-10);
            Ptx = fzero(fun3,bounds);
            PF = d1*Ptx-d4;
            PS = 0;
        end
    end
    
    PA = parameters.C_scal*parameters.Nu(t)*Ptx + parameters.C_offset - PF;
    fval2 = a1*PF + a2*PF*PF + b1*PS + c1*Ptx + c2*exp(-c3*Ptx) + parameters.lambda(t)*(1-parameters.epsilon);
    P2 = [PF,Ptx,PS,PA];

    if fval1<fval2
        P = P1;
    else
        P = P2;
    end

end