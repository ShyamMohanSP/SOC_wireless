function [u,delu_a] = hjb_solver(parameters,T,Nt,N1,N2,N3)

    dt = T/Nt;
    Tau_t = 0:dt:T;
    da = 1/N1;
    Tau_a = 0:da:1;
    dr = 1/N2;
    Tau_r = 0:dr:1;
    dxi = 1/N3;
    Tau_xi = 0:dxi:1;

    % solution arrays
    u = zeros((N1+1)*(N2+1)*(N3+1),Nt+1);   
    delu_a = zeros((N1+1)*(N2+1)*(N3+1),Nt+1);

    % control arrays
    PF = zeros((N1+1)*(N2+1)*(N3+1),Nt+1);      %FFS bought
    Ptx = zeros((N1+1)*(N2+1)*(N3+1),Nt+1);     %transmitted power
    PS = zeros((N1+1)*(N2+1)*(N3+1),Nt+1);      %sold RE power
    PA = zeros((N1+1)*(N2+1)*(N3+1),Nt+1);      %used RE power

    %final time condition
    u(:,end) = zeros((N1+1)*(N2+1)*(N3+1),1);   
    delu_a(:,end) = zeros((N1+1)*(N2+1)*(N3+1),1);
    for i=1:1:N1+1
        for j=1:1:N2+1
            for k=1:1:N3+1
                u(convert_coord(i,j,k,N1,N2,N3),end) = -parameters.P_k*parameters.Abar*Tau_a(i);
                delu_a(convert_coord(i,j,k,N1,N2,N3),end) = -parameters.P_k*parameters.Abar;
            end
        end
    end
    
    for n = Nt+1:-1:2

        F1save = zeros((N1+1)*(N2+1)*(N3+1),1);

        %update rule for value function at time n-1
        for i=1:1:N1+1
            for j=1:1:N2+1
                for k=1:1:N3+1

                    controls = compute_controls(parameters,Tau_t(n),[Tau_a(i),Tau_r(j),Tau_xi(k)],delu_a(convert_coord(i,j,k,N1,N2,N3),n));
                    PF(convert_coord(i,j,k,N1,N2,N3),n) = controls(1);
                    Ptx(convert_coord(i,j,k,N1,N2,N3),n) = controls(2);
                    PS(convert_coord(i,j,k,N1,N2,N3),n) = controls(3);
                    PA(convert_coord(i,j,k,N1,N2,N3),n) = controls(4);

                    F1 = parameters.F1_fn(Tau_t(n),Tau_r(j),controls);
                    F2 = parameters.F2_fn(Tau_t(n),Tau_r(j));
                    F3 = parameters.F3_fn(Tau_xi(k));
                    G1 = parameters.G1_fn(Tau_r(j));
                    G2 = parameters.G2_fn(Tau_xi(k));
                    H = parameters.H_fn(Tau_t(n),Tau_xi(k),controls);

                    if i == 1
                        ucoord2 = 0;
                    else
                        ucoord2 = u(convert_coord(i-1,j,k,N1,N2,N3),n);
                    end

                    if i == N1+1
                        ucoord1 = 0;
                    else
                        ucoord1 = u(convert_coord(i+1,j,k,N1,N2,N3),n);
                    end

                    if j == 1
                        ucoord4 = 0;
                    else
                        ucoord4 = u(convert_coord(i,j-1,k,N1,N2,N3),n);
                    end

                    if j== N2+1
                        ucoord3 = 0;
                    else
                        ucoord3 = u(convert_coord(i,j+1,k,N1,N2,N3),n);
                    end

                    if k == 1
                        ucoord6 = 0;
                    else
                        ucoord6 = u(convert_coord(i,j,k-1,N1,N2,N3),n);
                    end

                    if k == N3+1
                        u(convert_coord(i,j,k,N1,N2,N3),n-1) = u(convert_coord(i,j,k,N1,N2,N3),n)*(1 - 2*G1*dt/dr/dr - abs(F1)*dt/da - abs(F2)*dt/dr + F3*dt/dxi) ...
                        + ucoord1*fplus(F1)*dt/da + ucoord2*dt*fminus(F1)/da + ucoord3*(G1*dt/dr/dr + fplus(F2)*dt/dr) + ucoord4*(G1*dt/dr/dr + fminus(F2)*dt/dr) ...
                        + ucoord6*(-F3*dt/dxi) + dt*H;
                    else
                       u(convert_coord(i,j,k,N1,N2,N3),n-1) = u(convert_coord(i,j,k,N1,N2,N3),n)*(1 - 2*G1*dt/dr/dr - 2*G2*dt/dxi/dxi - abs(F1)*dt/da - abs(F2)*dt/dr - abs(F3)*dt/dxi) ...
                        + ucoord1*fplus(F1)*dt/da + ucoord2*dt*fminus(F1)/da + ucoord3*(G1*dt/dr/dr + fplus(F2)*dt/dr) + ucoord4*(G1*dt/dr/dr + fminus(F2)*dt/dr) ...
                        + u(convert_coord(i,j,k+1,N1,N2,N3),n)*(G2*dt/dxi/dxi + fplus(F3)*dt/dxi) + ucoord6*(G2*dt/dxi/dxi + fminus(F3)*dt/dxi) + dt*H; 
                    end

                    if (1 - 2*G1*dt/dr/dr - 2*G2*dt/dxi/dxi - abs(F1)*dt/da - abs(F2)*dt/dr - abs(F3)*dt/dxi) < 0
                        disp('Unstable scheme. Check discretization parameters');
                        return
                    end

                    F1save(convert_coord(i,j,k,N1,N2,N3),1) = F1;

                end
            end
        end

        %save gradient information at time step n-1
        for i=1:1:N1+1
            for j=1:1:N2+1
                for k=1:1:N3+1
                    if i == 1
                        delu_a(convert_coord(i,j,k,N1,N2,N3),n-1) = (u(convert_coord(i+1,j,k,N1,N2,N3),n-1) - u(convert_coord(i,j,k,N1,N2,N3),n-1))/da;
                    elseif i == N1+1
                        delu_a(convert_coord(i,j,k,N1,N2,N3),n-1) = (u(convert_coord(i,j,k,N1,N2,N3),n-1) - u(convert_coord(i-1,j,k,N1,N2,N3),n-1))/da;
                    else
                        if F1save(convert_coord(i,j,k,N1,N2,N3),1) >= 0
                            delu_a(convert_coord(i,j,k,N1,N2,N3),n-1) = (u(convert_coord(i+1,j,k,N1,N2,N3),n-1) - u(convert_coord(i,j,k,N1,N2,N3),n-1))/da;
                        else
                            delu_a(convert_coord(i,j,k,N1,N2,N3),n-1) = (u(convert_coord(i,j,k,N1,N2,N3),n-1) - u(convert_coord(i-1,j,k,N1,N2,N3),n-1))/da;
                        end
                    end

                end
            end
        end

    end

end