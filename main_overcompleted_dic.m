clear

% rng default
n = 8;
nrf = 4
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


% rk = [randn(1);randn(n-1,1)+1j*randn(n-1,1)];
% R = toeplitz(rk,conj(rk));
% 
% r = [R(1,end:-1:1,:).';R(2:end,1)];


Sigma = zeros(nrf*nrf,2*nrf-1);
Sigma(:,nrf) = reshape(eye(nrf),[nrf*nrf,1]);
for i = 1 : nrf-1
    Qi = [zeros(nrf-i,i), eye(nrf-i);zeros(i,nrf-i), zeros(i,i)];
    Sigma(:,nrf-i) = reshape(Qi,[nrf*nrf,1]);
    Sigma(:,nrf+i) = reshape(Qi.',[nrf*nrf,1]); 
end
% pinv_Sigma = (Sigma.'*Sigma)\Sigma.';
% Sm_dic = pagemtimes(pagemtimes(pagectranspose(codebook),R),codebook);
% sm_tilde = zeros(2*nrf-1,M);
% sm_vec = zeros(nrf*nrf,M);
% for m = 1 : M
%     Sm = Sm_dic(:,:,m);
%     sm_vec(:,m) = Sm(:);
%     sm_tilde(:,m) = pinv_Sigma*Sm(:);
% end
% s_vec = sm_vec(:);


% r_real = zeros(2*n-1,1);
% r_real(1) = rk(1);
% r_real(2:2:end) = real(rk(2:end));
% r_real(3:2:end) = imag(rk(2:end));




P = [];
for q = 0:2*nrf-2
    Puv = [zeros(2*Q-1,q*Q),eye(2*Q-1),zeros(2*Q-1,2*n-(q+2)*Q)];
    P = [P;Puv];
end

q_vec = (-(Q-1):Q-1).';
a = 1-abs(q_vec)/Q;
selected_frequencies = [];
Omega = a.'.*exp(-1j*angles.'.*q_vec.');
I_bnrf = eye(2*nrf-1);
for m = 0 : M-1
    phi_m = 2*pi/(2*Q-1)*(m-(2*Q-2)/(2*Q-1));
    fm = exp(-1j*phi_m*q_vec);
    Omega_m = kron(I_bnrf,(a.*fm).');
    em = zeros(2*Q-1,1);
    em(m+1) = 1;
    
    sel = kron(eye(2*nrf-1),em.'*Omega);
    selected_frequencies = [selected_frequencies;sel];
end

% norm(s_vec-kron(eye(M),Sigma)*Omega_ext*P*r)

Psi = kron(eye(M),Sigma)*selected_frequencies*P;

% norm(s_vec-Psi*r)



MC = 1e3;
snrdb = linspace(-10,40,30);
snr = power(10,snrdb/10);
theta_deg = [-45,45];
theta = theta_deg*pi/180;
n_vec = (0:n-1)';
L = length(theta);

K = 5*315;
Km = K/M

sigman2 = 1;
nrf2 = nrf*nrf;
se_gls = zeros(MC,length(snrdb));
se_dig = zeros(MC,length(snrdb));

crb_hyb = ToolsObj.eval_crb_theta_per_snr(1/2,1,theta_deg,snr,codebook,Km);
rcrb_hyb_deg = sqrt(mean(crb_hyb,1))*180/pi;

crb_dig = ToolsObj.eval_crb_theta_per_snr(1/2,1,theta_deg,snr,eye(n),K);
rcrb_dig_deg = sqrt(mean(crb_dig,1))*180/pi;

Rs = eye(L);

for j = 1 : length(snr)
    A = exp(1j*pi*sin(theta).*n_vec);
    sigman2 = 1/snr(j);
    for i = 1 : MC
        s = sqrt(Rs/2)*(randn(L,K)+1j*randn(L,K));
        noise = sqrt(sigman2/2)*(randn(n,K)+1j*randn(n,K));
        
        x = A*s+noise;
        xm_dic = reshape(x,[n,Km,M]);
        
        R_epsilon = zeros(M*nrf2,M*nrf2);
        hat_sm_dic = zeros(nrf2,M);
        for m = 1 : M
            Ym = codebook(:,:,m)'*xm_dic(:,:,m);
            hat_Sm = 1/Km*(Ym*Ym');
            hat_sm_dic(:,m) = hat_Sm(:);
            R_epsilon((m-1)*nrf2+1:m*nrf2,(m-1)*nrf2+1:m*nrf2) = 1/Km*kron(hat_Sm.',hat_Sm);
        end
        hat_sm = hat_sm_dic(:);
        hat_r_ext = (Psi'/R_epsilon*Psi)\Psi'/R_epsilon*hat_sm;
        hat_r = hat_r_ext(n:end);
        hat_r(1) = real(hat_r(1));
        hat_R = toeplitz(hat_r,conj(hat_r));

        hat_theta_deg = rootmusic(hat_R,L);
        se_gls(i,j) = mean(power((theta_deg-hat_theta_deg),2));

        
        % hat_R_dig = (x*x')/K;
        % hat_theta_deg_dig = rootmusic(hat_R_dig,L);
        % se_dig(i,j) = power((theta_deg-hat_theta_deg_dig),2);
    end
    fprintf("j: %d\n",j)
end

rmse_gls_deg = sqrt(mean(se_gls));
% rmse_dig = sqrt(mean(se_dig));

hold on
semilogy(snrdb,rmse_gls_deg,'-or')
semilogy(snrdb,rcrb_hyb_deg,'-b')
% semilogy(snrdb,rcrb_dig_deg,'-k')
hold off
set(gca,"YScale","log")

results = struct();
results.n = n;
results.nrf = nrf;
results.M = M;
results.K = K;
results.snrdb = snrdb;
results.theta_deg = theta_deg;
results.MC = MC;
results.rmse_gls_deg = rmse_gls_deg;
results.rcrb_hyb_deg = rcrb_hyb_deg;
results.rcrb_dig_deg = rcrb_dig_deg;
name = "results_nrf4.mat"
save(name, '-struct', 'results');

