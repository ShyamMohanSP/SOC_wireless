function [f, g] = my_lmbm_tfunc(lambda)
    
    global pars Nt N1 N2 N3 M_dual Nt1 M_SG;
    pars.lambda_ell = length(lambda);
    pars.lambda_discrete = lambda;
    T = pars.horizon;
    [value_fxn,its_derivative] = hjb_solver(pars,T,Nt,N1,N2,N3);
    dual = estimate_dual_function(pars,N1,N2,N3,M_dual,value_fxn);
    f = -1.0*dual(1);
    sg = estimate_subgradient(pars,Nt1,Nt,N1,N2,N3,M_SG,its_derivative);
    g = -1.0*sg(1:pars.lambda_ell);

end