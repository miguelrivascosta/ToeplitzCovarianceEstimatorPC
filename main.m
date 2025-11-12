clear, close all

n = 10;
nrf = 2;
Q = n/nrf;



rng default

r_vec = [randn(1);randn(n-1,1)+1j*randn(n-1,1)];
r_ext = zeros(2*n-1,1);
r_ext(1) = r_vec(1);
r_ext(2:2:end) = real(r_vec(2:end));
r_ext(3:2:end) = imag(r_vec(2:end));
% psi_u = 2*pi/n*u;
% psi_v = 2*pi/n*v;
% lambda = exp(-1j*psi_u);
% mu = exp(1j*psi_v);
% lambda = 0.3;
% mu = 0.8;
% i = 0;
% j = 1;
% li_lambda = eval_li(n,lambda,Q,i);
% lj_mu = eval_li(n,mu,Q,j);
% R = toeplitz(r,r');
% S_true = li_lambda.'*R*lj_mu;
% S_ut = eval_Sij(r,lambda,mu,i,j,Q);

R = toeplitz(r_vec,r_vec');

i = 0;
j = 0;

u = 0;
v = 1;

li_u = eval_liu(n,u,Q,i);
lj_v = eval_liu(n,v,Q,j);
Sij_uv = li_u'*R*lj_v;

%%% comprobacion para N=4 y Nrf=2 si la matriz de beamforming extendida
%%% tiene suficiente rango.
% B0 = [eval_liu(4,0,2,0),eval_liu(4,0,2,1)];
% B1 = [eval_liu(4,1,2,0),eval_liu(4,1,2,1)];
% B2 = [eval_liu(4,0,2,0),eval_liu(4,1,2,1)];
% B3 = [eval_liu(4,1,2,0),eval_liu(4,0,2,1)];
% B_ext = [kron(B0.',conj(B0));kron(B1.',conj(B1))];
% rank(B_ext)
% F = dftmtx(n);
% rank(kron(conj(F),F))
%%% Check eq. 8
Sij_uv_ut = 0;
for r = i*Q:(i+1)*Q-1
    for s = j*Q:(j+1)*Q-1
        Sij_uv_ut = Sij_uv_ut + exp(-1j*2*pi/Q*u*r)*exp(1j*2*pi/Q*v*s)*R(r+1,s+1);
    end
end
error = norm(Sij_uv-Sij_uv_ut); 
fprintf("Error in eq.8: %d\n",error);


%%% Check eq. 9
Sij_uv_ut = 0;
for q = -(Q-1):Q-1
    idx_r = (i-j)*Q+q;
    if idx_r>0
        r_internal = r_vec(idx_r+1);
    else
        r_internal = conj(r_vec(abs(idx_r)+1));
    end
    
    aux = 0;
    for s = max(0,q):min(Q-1,Q-1+q)
        aux = aux + exp(1j*2*pi/Q*(v-u)*s);
    end

    Sij_uv_ut = Sij_uv_ut + r_internal*exp(-1j*2*pi/Q*v*q)*aux;
end
error = norm(Sij_uv-Sij_uv_ut); 
fprintf("Error in eq.9: %d\n",error);


%%% Check eq. 9.1 
Sij_uv_ut = 0;
for q = -(Q-1):Q-1
    idx_r = (i-j)*Q+q;
    if idx_r>0
        r_internal = r_vec(idx_r+1);
    else
        r_internal = conj(r_vec(abs(idx_r)+1));
    end
    
    aux = 0;
    if u==v
        Sij_uv_ut = Sij_uv_ut + r_internal*(Q-abs(q))*exp(-1j*2*pi/Q*u*q);
    else
        if q<0
            Sij_uv_ut = Sij_uv_ut + r_internal*exp(-1j*2*pi/Q*v*q)*(1-exp(1j*2*pi/Q*(v-u)*q))/(1-exp(1j*2*pi/Q*(v-u)));
        else
            Sij_uv_ut = Sij_uv_ut - r_internal*exp(-1j*2*pi/Q*v*q)*(1-exp(1j*2*pi/Q*(v-u)*q))/(1-exp(1j*2*pi/Q*(v-u)));
        end       
    end
end
error = norm(Sij_uv-Sij_uv_ut); 
fprintf("Error in eq.9.1: %d\n",error);


%%% Check eq. 10

Sij_uv_ut = 0;
q_vec = (-(Q-1):Q-1)';
rd = zeros(2*Q-1,1);
for q = -(Q-1):Q-1
    idx_r = (i-j)*Q+q;
    if idx_r>0
        rd(q+Q)= r_vec(idx_r+1);
    else
        rd(q+Q) = conj(r_vec(abs(idx_r)+1));
    end
end

r_neg = zeros(Q,1);
r_pos = zeros(Q,1);

for q = 0:Q-1
    idx_r_neg = (i-j)*Q-q;
    if idx_r_neg>=0
        r_neg(q+1) = r_vec(idx_r_neg+1);
    else
        r_neg(q+1) = conj(r_vec(abs(idx_r_neg)+1));
    end

    idx_r_pos = (i-j)*Q+q;
    if idx_r_pos>=0
        r_pos(q+1) = r_vec(idx_r_pos+1);
    else
        r_pos(q+1) = conj(r_vec(abs(idx_r_pos)+1));
    end
end

weight_u_neg = [1/2;exp(1j*2*pi/Q*u*(1:Q-1).')];
weight_u_pos = [1/2;exp(-1j*2*pi/Q*u*(1:Q-1).')];
weight_v_neg = [1/2;exp(1j*2*pi/Q*v*(1:Q-1).')];
weight_v_pos = [1/2;exp(-1j*2*pi/Q*v*(1:Q-1).')];

exp_u = exp(-1j*2*pi/Q*u*q_vec);
exp_v = exp(-1j*2*pi/Q*v*q_vec);

Sij_uv_ut = 0;
if u == v
    % Su = sum(rd.*exp_u);
    % Su_tilde = sum(abs(q_vec).*rd.*exp_u);
    % Sij_uv_ut = Q*Su-Su_tilde;
    
    Su_pos = sum(r_pos.*weight_u_pos);
    Su_neg = sum(r_neg.*weight_u_neg);

    dot_Su_pos = sum(r_pos.*(0:Q-1).'.*weight_u_pos);
    dot_Su_neg = sum(r_neg.*(0:Q-1).'.*weight_u_neg);

    % Sij_uv_ut = Q*(Su_pos+Su_neg) - (dot_Su_neg+dot_Su_pos);
    Sij_uv_ut = (Q*Su_pos-dot_Su_pos) + (Q*Su_neg-dot_Su_neg);
else
    % Su = sum(rd.*exp_u);
    % Sv = sum(rd.*exp_v);
    % Sij_uv_ut = -1/(1-exp(1j*2*pi/Q*(v-u)))*(Su-Sv);
    
    Su_pos = sum(r_pos.*weight_u_pos);
    Su_neg = sum(r_neg.*weight_u_neg);
    Sv_pos = sum(r_pos.*weight_v_pos);
    Sv_neg = sum(r_neg.*weight_v_neg);

    Sij_uv_ut = 1/(1-exp(1j*2*pi/Q*(v-u)))*((Sv_neg-Sv_pos)-(Su_neg-Su_pos));
    % for q = 0:Q-1
    %     idx_r0 = (i-j)*Q-q;
    %     idx_r1 = (i-j)*Q+q;
    % 
    %     if idx_r0>=0
    %         r_internal0 = r_vec(idx_r0+1);
    %     else
    %         r_internal0 = conj(r_vec(abs(idx_r0)+1));
    %     end
    % 
    %     if idx_r1>=0
    %         r_internal1 = r_vec(idx_r1+1);
    %     else
    %         r_internal1 = conj(r_vec(abs(idx_r1)+1));
    %     end
    % 
    %     % Sij_uv_ut = Sij_uv_ut + r_internal0*(exp(1j*2*pi/Q*v*q)-exp(1j*2*pi/Q*u*q))...
    %     % -r_internal1*(exp(-1j*2*pi/Q*v*q)-exp(-1j*2*pi/Q*u*q));
    % 
    %     Sij_uv_ut = Sij_uv_ut + r_internal0*exp(1j*2*pi/Q*v*q)-r_internal1*exp(-1j*2*pi/Q*v*q)...
    %         - r_internal0*exp(1j*2*pi/Q*u*q)+r_internal1*exp(-1j*2*pi/Q*u*q);
    % end
    % Sij_uv_ut = 1/(1-exp(1j*2*pi/Q*(v-u)))*Sij_uv_ut;
end

error = norm(Sij_uv-Sij_uv_ut); 
fprintf("Error in eq.10: %d\n",error);


%% Check equation xx
aq = eval_aq(Q);
bq = eval_bq(u,v,Q);

Sij_uv_ut = 0;
if u == v
    Sij_uv_ut = sum(aq.*r_neg.*weight_u_neg +aq.*r_pos.*weight_u_pos);
else
    Sij_uv_ut = 1/(1-exp(1j*2*pi/Q*(v-u)))*sum(-conj(bq).*r_neg.*weight_u_neg+bq.*r_pos.*weight_u_pos);
end
error = norm(Sij_uv-Sij_uv_ut); 
fprintf("Error in eq.xx: %d\n",error);


%%%%%%%%%%%%%%%%%%%%%%%%%%

% i = 0;
% j = 0;
% u = 1;
% v = 1;
% 
% J0 = zeros(2*Q-1,2*n-1);
% J0(1,2) = 1;
% J0(1,3) = -1j;
% J0(2,1)=1;
% J0(3,2) = 1;
% J0(3,3) = 1j;
% Sij_uv = eval_liu(n,u,Q,i)'*R*eval_liu(n,u,Q,j);
% 
% q_vec = (-(Q-1):Q-1).';
% rd = [r_neg(end:2);r_pos];
% exponentials = exp(-1j*2*pi/Q*q_vec*u);
% aq_ext = Q-abs(q_vec);
% 
% Sij_uv_ut = sum(aq_ext.*rd.*exponentials);
% 
% error = norm(Sij_uv-Sij_uv_ut); 
% fprintf("Error in eq.xx: %d\n",error);
% 
% 
% lu = gen_lluu(u,Q);
% 
% Sij_uv_ut = lu.'*J0*r_ext;
% error = norm(Sij_uv-Sij_uv_ut); 
% fprintf("Error in eq.xx: %d\n",error);
% 
% L = [gen_lluu(0,Q).';gen_lluu(1,Q).'];


%% comprobacion de

%%
i = 0;
j = 0;
u = 0;
v = 1;

psi_u = 2*pi/Q*u;
psi_u_tilde = 2*pi/Q*(u+1/2);

q_vec = (-(Q-1):Q-1).';

aq = Q-abs(q_vec);
bq = sin(pi/Q*abs(q_vec));
exp_psi = exp(-1j*psi_u*q_vec);
exp_psi_tilde = exp(-1j*psi_u_tilde*q_vec);

Sij_uv = eval_liu(n,u,Q,i)'*R*eval_liu(n,v,Q,j);

if u==v
    Sij_uv_ut = sum(aq.*rd.*exp_psi);
else
    Sij_uv_ut = 2j/(1-exp(1j*2*pi/Q))*sum(bq.*rd.*exp_psi_tilde);
end
error = norm(Sij_uv-Sij_uv_ut); 
fprintf("Error in eq.xx: %d\n",error);







J = gen_J(n);
q_vec = -(Q-1):Q-1;
i = 1;
j = 0;
% u = 0;
% v = 1;
d = abs(i-j);
Jd = J(n+d*Q+q_vec,:);
rd = Jd*r_ext;
%%%% check relationship of Sij with the 2DFT of the matrix:
FQ = exp(1j*2*pi/Q*(0:Q-1).'.*(0:Q-1));

Rd = toeplitz(rd(Q:end),rd(Q:-1:1));
Sd = FQ'*Rd*FQ;

Sij_mat = zeros(Q,Q);
for u = 0 : Q-1
    for v = 0 : Q-1
        Sij_mat(u+1,v+1) = eval_Sijuv(i,j,u,v,Q,R);
    end
end
norm(Sij_mat-Sd)

% lambda = exp(-1j*2*pi/Q*u);
% mu = exp(1j*2*pi/Q*v);
% li = eval_li(n,lambda,Q,i);
% lj = eval_li(n,mu,Q,j);
% S = li.'*R*lj;
% eval_Sij_uv(r,u,v,i,j,Q)
% eval_Sij(r,lambda,mu,i,j,Q)

function y = eval_Sijuv(i,j,u,v,Q,R)
    n = size(R,1);
    y = eval_liu(n,u,Q,i)'*R*eval_liu(n,v,Q,j);
end


function J = gen_J(n)
    %%%% Diseño de J
    J = zeros(2*n-1);
    
    % superior part
    ini = (2*n-2:2*n-1);
    for i = 1:n-1
       J(i,ini-(i-1)*2) = [1,-1j];
    end
    % middle
    ini = 2:3;
    J(n,1) = 1;
    % inferior part
    for i = 1:n-1
        J(n+i,ini+(i-1)*2) = [1,1j];
    end
end

function y = gen_lluu(u,Q)
    q_vec = -(Q-1):(Q-1);
    exponentials = exp(-1j*2*pi/Q*q_vec*u);
    ak = (Q-abs(q_vec));
    y = ak.'.*exponentials.';
end

function y = eval_aq(Q)
    y = Q-(0:Q-1)';
end
function y = eval_bq(u,v,Q)
    y = 1-exp(-1j*2*pi/Q*(v-u)*(0:Q-1)');
end

function y = eval_liu(n,u,Q,i)
    aux = power(exp(1j*2*pi/Q*u),0:n-1);
    y = zeros(n,1);
    idx = i*Q+1:(i+1)*Q;
    y(idx) = aux(idx);
end

function y = eval_li(n,lambda,Q,i)
    aux = power(lambda,(0:n-1)');
    y = zeros(n,1);
    idx = i*Q:(i+1)*Q-1;
    y(idx+1) = aux(idx+1);
end

function Sij = eval_Sij_uv(r,u,v,i,j,Q)
    Sij = 0;
    d = i-j;
    psi_u = 2*pi/Q*u;
    psi_v = 2*pi/Q*v;

    mu = exp(1j*psi_v);
    lambda = exp(-1j*psi_u);
    if u~=v
        % for s = -(Q-1):Q-1
        %     idx_r = (i-j)*Q+s;
        %     if idx_r>=0
        %         r_idx = r(idx_r+1);
        %     else
        %         r_idx = conj(r(abs(idx_r)+1));
        %     end
        % 
        % 
        %     % if s>=0
        %     %     Sij = Sij + r_idx*( exp(-1j*2*pi/Q*u*s)-exp(-1j*2*pi/Q*v*s)  )/(1-exp(1j*2*pi/Q*(v-u)));
        %     % 
        %     % else
        %     %     Sij = Sij + r_idx*( exp(-1j*2*pi/Q*v*s)-exp(-1j*2*pi/Q*u*s)  )/(1-exp(1j*2*pi/Q*(v-u)));
        %     % end
        % 
        %     if s>=0
        %         Sij = Sij + r_idx*( exp(-1j*2*pi/Q*u*s)-exp(-1j*2*pi/Q*v*s)  )/(1-exp(1j*2*pi/Q*(v-u)));
        %     else
        %         Sij = Sij + r_idx*( exp(-1j*2*pi/Q*v*s)-exp(-1j*2*pi/Q*u*s)  )/(1-exp(1j*2*pi/Q*(v-u)));
        %     end
        % end

        for s = 1 : Q-1
            idx_r0 = (i-j)*Q+s;
            idx_r1 = (i-j)*Q-s;
            if idx_r0>=0
                r0 = r(idx_r0+1);
            else
                r0 = conj(r(abs(idx_r0)+1));
            end

            if idx_r1>=0
                r1 = r(idx_r1+1);
            else
                r1 = conj(r(abs(idx_r1)+1));
            end

            Sij = Sij + r0*(exp(-1j*2*pi/Q*u*s)-exp(-1j*2*pi/Q*v*s))/(1-exp(1j*2*pi/Q*(v-u))) + ...
                  r1*(exp(1j*2*pi/Q*v*s)-exp(1j*2*pi/Q*u))/(1-exp(1j*2*pi/Q*(v-u)));
        end
    else
    end
end

function Sij = eval_Sij(r,lambda,mu,i,j,Q) 
    % Sij = 0;
    % for p = 0 : Q-1
    %     for q = 0 : Q-1
    %         idx = (i-j)*Q+(p-q);
    % 
    %         if idx>=0
    %             Sij = Sij + r(idx+1)*power(lambda,i*Q+p)*power(mu,j*Q+q);
    %         else
    %             Sij = Sij + conj(r(abs(idx)+1))*power(lambda,i*Q+p)*power(mu,j*Q+q);
    %         end
    % 
    %     end
    % end
    % Sij = 0;
    % for s = -(Q-1):Q-1
    %     idx_r = (i-j)*Q+s;
    % 
    %     pmin = max(0,s);
    %     pmax = min(Q-1,Q-1+s);
    % 
    %     aux = sum(power(lambda*mu,pmin:pmax));
    %     if idx_r>=0
    %         Sij = Sij + r(idx_r+1)*mu^(-s)*aux;
    %     else
    %         Sij = Sij + conj(r(abs(idx_r)+1))*mu^(-s)*aux;
    %     end
    % end
    % Sij = lambda^(i*Q)*mu^(j*Q)*Sij;

    Sij = 0;
    eta = mu*lambda;
    for s = -(Q-1):Q-1
        idx_r = (i-j)*Q+s;
        if idx_r>=0
            r_idx = r(idx_r+1);
        else
            r_idx = conj(r(abs(idx_r)+1));
        end
        
        if norm(eta-1)>1e-12
            if s>=0
                Sij = Sij + r_idx*power(mu,-s)*(power(eta,s) - power(eta,Q) )/(1-eta);
            else
                Sij = Sij + r_idx*power(mu,-s)*(1 - power(eta,Q+s) )/(1-eta);
            end
        else
            if s>=0
                Sij = Sij + r_idx*power(mu,-s)*(Q-s);
            else
                Sij = Sij + r_idx*power(mu,-s)*(Q+s);
            end
        end
    end
    Sij = lambda^(i*Q)*mu^(j*Q)*Sij;



end