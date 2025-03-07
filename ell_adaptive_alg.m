%% initialization

global pars N1 N2 N3 Nt Nt1 M_SG M_dual;

pars = MyParameters;

pars.horizon = 2;       %number of days horizon for algorithm
pars.Pareto = 0.5;        %pareto opt constant - equal weightage to profits and clean energy;

pars.z_BS = [0,0];   %location of base station (center of populated area)
pars.user_dist = 1;     %distribution of users(0 if uniform, 1 if gaussian)

pars.kappa = 1;      %transmission loss constant
pars.sigma0 = 3.1623e-8;    %transmission noise level(in W)
pars.SNR_th = 15;            %signal to noise threshold (in dB)
pars.eta = 2;                %transmission loss exponent for free-space path loss

pars.C1 = 4e-4;         %environmental penalty coefficients
pars.C2 = 1e-4;

pars.phi_th = 0.001;       %outage proportion threshold

pars.epsilon = 0.1;       %95p.c. confidence of satisfying QOS

pars.Pbar_R = 10000;       %maximum RE capacity (in W)

pars.sigma = 3;          %wireless fading SDE coefficients
pars.B = 1;

pars.Pbar_tx = 5000;     %maximum transmitted power(W)
pars.A_omega = 4.9e5;     %area of region in m^2 (700 m) squared
pars.Nu_max = 2000;      %maximum number of users in region
pars.std_users = 300;    %standard deviation of users around the base station (in m)
pars.cell_price = 1;       %unit price of cellular service

pars.Abar = 10000;       %maximum battery charge (in Wh)
pars.A0 = 0.5;          %percentage of initial battery charge
pars.Pu_A = 7500;         %maximum battery influx (in W)       
pars.Pbar_A = 30000;       %maximum battery outflux (in W)

pars.C_scal = 7.84;      %cellular BS power model coefficients
pars.C_offset = 71.5;    %(in W)

% pars.xi_u = 0.3295;     %99p.c. confidence interval foor gamma distribution
% pars.xi_bar = 9.229;
% fun = @(x) x*exp(1-x) - pars.epsilon;
% pars.tau_threshold = 0.5*pars.phi_th/fzero(fun,10);
pars.xi_u = -2*pars.Nu_max*pars.std_users^2*pars.sigma0*10^(0.1*pars.SNR_th)*log(pars.phi_th)/pars.Pbar_tx/pars.kappa;
fun = @(x) gamcdf(x-pars.xi_u,pars.sigma,1) - 0.95;   %95p.c. quantile of invariant distribution
pars.xi_bar = fzero(fun,pars.sigma);         %xi space boundary 
pars.xi_bar = ceil(pars.xi_bar);

%% data processing

pars.alpha = 0.34063;        %RE SDE coefficients calibrated from 2023 data
pars.theta0 = 2.3948;
pars.deltaCalib = 0.054;

raw_data = csvread('50Hertz/2024_wind_data_50Hertz.csv');
day_num = 100;            %choose out of 238 days

pars.re_timesteps = pars.horizon - 1 + raw_data(day_num+pars.horizon-1,1:97);
pars.p_data = raw_data(day_num+pars.horizon-1,98:194);
if pars.horizon > 1
    for n=pars.horizon-1:-1:1
        pars.re_timesteps = [n - 1 + raw_data(day_num+n-1,1:96),pars.re_timesteps];
        pars.p_data = [raw_data(day_num+n-1,98:193),pars.p_data];
    end
end

pars.pdot_data = zeros(1,length(pars.p_data)-1);
pars.pdot_data(1) = (pars.p_data(2)-pars.p_data(1))/(pars.re_timesteps(2)-pars.re_timesteps(1));
for i=2:1:length(pars.p_data)-1
    pars.pdot_data(i) = (pars.p_data(i+1)-pars.p_data(i-1))/(pars.re_timesteps(i+1)-pars.re_timesteps(i-1));
end

clear raw_data;

raw_data = csvread('50Hertz/dayahead_price_50Hertz_2024.csv');

pars.kb_timesteps = pars.horizon - 1 + raw_data(day_num+pars.horizon-1,1:25);
pars.k_data = raw_data(day_num+pars.horizon-1,26:50);
if pars.horizon > 1
    for n=pars.horizon-1:-1:1
        pars.kb_timesteps = [n - 1 + raw_data(day_num+n-1,1:24),pars.kb_timesteps];
        pars.k_data = [raw_data(day_num+n-1,26:49),pars.k_data];
    end
end

clear raw_data;

%determine price of selling battery energy at the end of the horizon
num = 0;
denom = 0;
temp = linspace(pars.horizon-1,pars.horizon,97);
for i = 1:1:length(temp)
    num = num + pars.p_fn(temp(i))*pars.K_b(temp(i));
    denom = denom + pars.p_fn(temp(i));
end
pars.P_k = num/denom;

T = pars.horizon;          %final time
Nx = 10;        %uniform discretization in all space dimensions
N1 = Nx;         %discretization in a-space
N2 = Nx;         %discretization in r-space
N3 = Nx;         %discretization in xi-space
Nt = pars.horizon*ceil(max([4*Nx*Nx*pars.B,4*Nx*Nx*pars.alpha*pars.theta0,8*pars.Pbar_A*Nx/pars.Abar]));     %time steps
Nt1 = 2^6;      %SDE discretization

%initialize lambda and tau_star as constant one function
pars.lambda_ell = 1;
pars.lambda_discrete = ones(pars.lambda_ell,1);

%% initialization algorithm

TOL_naive = 1;       %relative subgradient tolerance
tol_outer = TOL_naive*pars.epsilon;     %absolute tolerance with respect to confidence of chance constraint
max_iter_outer = 10;     %maximum iterations of lambda search     
i = 1;
betaF = 5;              %factor of increase in lambda search
M_SG = 10;                 %number of samples used to estimate subgradient
M_dual = 10;               %number of samples used to evaluate dual function

[value_fxn,its_derivative] = hjb_solver(pars,T,Nt,N1,N2,N3);
temp_dual = estimate_dual_function(pars,N1,N2,N3,M_dual,value_fxn);
dual = temp_dual(1);
sg = estimate_subgradient(pars,Nt1,Nt,N1,N2,N3,M_SG,its_derivative);
init_array = [i,pars.lambda_discrete',temp_dual,sg',norm(sg)];

while norm(sg) > sqrt(pars.lambda_ell)*tol_outer && i <= max_iter_outer && any(sg > 0)
    
    pars.lambda_discrete = betaF.*pars.lambda_discrete;
    [value_fxn,its_derivative] = hjb_solver(pars,T,Nt,N1,N2,N3);
    temp_dual = estimate_dual_function(pars,N1,N2,N3,M_dual,value_fxn);
    dual = temp_dual(1);
    sg = estimate_subgradient(pars,Nt1,Nt,N1,N2,N3,M_SG,its_derivative);
    i = i + 1;
    init_array = [init_array;[i,pars.lambda_discrete',temp_dual,sg',norm(sg)]];
    
end

%% deterministic LMBM

rpar = [0.1,-1,0,1e-2,1e-2,0.5,0.2,10];
ipar = [0, 50, 50, 5, -1, 0, 0];
M_SG = 10;                 %number of samples used to estimate subgradient
M_dual = 10;               %number of samples used to evaluate dual function

[x,fval,niter,nfeval,term,time] = lmbm_driver('my_lmbm_tfunc',pars.lambda_discrete,pars.lambda_ell,1,0,3,5,3,rpar,ipar);

[ans1,ans2] = my_lmbm_tfunc(x);

%% adaptive dual optimizer

max_iter = 2;           %maximum outer loop iterations
TOL_A = 0.1;           %relative subgradient tolerance of double loop algorithm
tol = pars.epsilon*TOL_A;
data_array = [];
lambda_array = [];
violation_data = [];
tgrid = 0:1/Nt1:pars.horizon;
warmup_iter = 10;

for j=1:1:6
    
    pars.lambda_ell = 2*pars.lambda_ell;
    if isscalar(pars.lambda_discrete)
        pars.lambda_discrete = repelem(pars.lambda_discrete,2)';
    else
        pars.lambda_discrete = repelem(pars.lambda_discrete,2);
    end
    [value_fxn,its_derivative] = hjb_solver(pars,T,Nt,N1,N2,N3);
    temp_dual = estimate_dual_function(pars,N1,N2,N3,M_dual,value_fxn);
    dual = temp_dual(1);
    [sg,mc_prob] = estimate_subgradient(pars,Nt1,Nt,N1,N2,N3,M_SG,its_derivative);
%    const = pars.lambda_discrete/pars.epsilon/TOL_naive;
    const = pars.lambda_discrete;
    k = 1;
    
    while k <= max_iter && norm(sg) > sqrt(pars.lambda_ell)*tol
        if k <= warmup_iter
            alpha_k = const/norm(sg);
        else
            alpha_k = const/norm(sg)/sqrt(k+1);
        end
        pars.lambda_discrete = pars.lambda_discrete + alpha_k.*sg;
        [value_fxn,its_derivative] = hjb_solver(pars,T,Nt,N1,N2,N3);
        temp_dual = estimate_dual_function(pars,N1,N2,N3,M_dual,value_fxn);
        dual = temp_dual(1);
        [sg,mc_prob] = estimate_subgradient(pars,Nt1,Nt,N1,N2,N3,M_SG,its_derivative);
        k = k + 1;
    end
    
    data_array = [data_array;[j,temp_dual,norm(sg)/sqrt(pars.lambda_ell)]];
    violation_data = [violation_data;mc_prob];
    
    temp_array = zeros(1,length(tgrid));
    for z = 1:1:length(tgrid)
        temp_array(z) = pars.lambda(tgrid(z));
    end
    lambda_array = [lambda_array; temp_array];

end