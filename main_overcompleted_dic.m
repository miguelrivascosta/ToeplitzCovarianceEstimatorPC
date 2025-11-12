clear, close all

rng default
n = 8;
nrf = 4;
Q  = n/nrf;
bq = 2*Q-1;
angles = 2*pi/bq*((0:bq-1)-(bq-1)/bq);
M = 2*Q-1;
codebook = zeros(n,nrf,M);
for i = 1:M
    for j = 0 : nrf-1
        idx = (j*Q:(j+1)*Q-1)+1;
        codebook(idx,j+1,i) = 1/sqrt(Q)*exp(1j*angles(i)*(0:Q-1).');
    end
end


rk = [randn(1);randn(n-1,1)+1j*randn(n-1,1)];
R = toeplitz(rk,conj(rk));

r = [R(1,end:-1:1,:).';R(2:end,1)];




r_real = zeros(2*n-1,1);
r_real(1) = rk(1);
r_real(2:2:end) = real(rk(2:end));
r_real(3:2:end) = imag(rk(2:end));

q = -(nrf-1);

aux = (-(n-1):(n-1)).';
P = [];
% P2 = [];
% for q = -(nrf-1):(nrf-1)
%     Pq = [zeros(2*Q-1,(q+nrf-1)*Q),eye(2*Q-1),zeros(2*Q-1,2*n-(q+nrf+1)*Q)];
%     % (q+nrf-1)*Q
%     2*n-(q+nrf+1)*Q
%     P = [P;Pq];
% end

for q = 0:2*nrf-2
    Puv = [zeros(2*Q-1,q*Q),eye(2*Q-1),zeros(2*Q-1,2*n-(q+2)*Q)];
    P = [P;Puv];
end




Sigma = zeros(nrf*nrf,2*nrf-1);
Sigma(:,nrf) = reshape(eye(nrf),[nrf*nrf,1]);
for i = 1 : nrf-1
    Qi = [zeros(nrf-i,i), eye(nrf-i);zeros(i,nrf-i), zeros(i,i)];
    Sigma(:,nrf-i) = reshape(Qi,[nrf*nrf,1]);
    Sigma(:,nrf+i) = reshape(Qi.',[nrf*nrf,1]); 
end




pinv_Sigma = (Sigma.'*Sigma)\Sigma.';


Sm_dic = pagemtimes(pagemtimes(pagectranspose(codebook),R),codebook);
sm_tilde = zeros(2*nrf-1,M);
sm_vec = zeros(nrf*nrf,M);
for m = 1 : M
    Sm = Sm_dic(:,:,m);
    sm_vec(:,m) = Sm(:);
    sm_tilde(:,m) = pinv_Sigma*Sm(:);
end
s_vec = sm_vec(:);

toeplitz_mat = zeros(2*n-1,2*n-1);

q_vec = -(Q-1):(Q-1);
aq = 1 - abs(q_vec)/Q;

u = 3;
Psi = [];
for u = 1 : M
    phi_u = angles(u);
    
    
    Psi_u = zeros(2*nrf-1,2*n-1);
    Psi_u(1,1) = aq(Q);
    Psi_u(1,2:2:2*Q-1) = 2*aq(Q+1:2*Q-1).*cos(phi_u*q_vec(Q+1:2*Q-1));
    Psi_u(1,3:2:2*Q-1) = 2*aq(Q+1:2*Q-1).*sin(phi_u*q_vec(Q+1:2*Q-1));
    
    idx_real_ini = 1:2:2*(2*Q-1);
    idx_imag_ini = 2:2:2*(2*Q-1);
    for i = 1:nrf-1
        idx_real = (idx_real_ini + (i-1)*2*Q)+1;
        idx_imag = (idx_imag_ini + (i-1)*2*Q)+1;
    
        Psi_u(2*i,idx_real) = aq.*exp(-1j*phi_u*q_vec);
        Psi_u(2*i,idx_imag) = 1j*aq.*exp(-1j*phi_u*q_vec);
        Psi_u(2*i+1,idx_real) = aq.*exp(1j*phi_u*q_vec);
        Psi_u(2*i+1,idx_imag) = -1j*aq.*exp(1j*phi_u*q_vec);
    end
    Psi = [Psi;Psi_u];
end
Sigma_ext = kron(eye(M),Sigma);
Gamma = Sigma_ext*Psi;

pinv_Gamma = (Gamma'*Gamma)\Gamma';
norm(r_real-pinv_Gamma*s_vec)



% ext = [aux,zeros(1,2*n-1-(2*Q-1))];
% aux*r_real
% sm_tilde(1,u)