classdef ToolsObj
    methods(Static)
        function a = a_ula(N,psi)
            a = exp(1j *(0:N-1).' .* psi);
        end
		
		function [CRB,J] = crb_unc_uncorrelated_psi(psi,Rs,sigman2,codebook,Km)
			[N,Nrf,M] = size(codebook);
			L = size(psi,2);
            A = ToolsObj.a_ula(N,psi);
    		R = A*Rs*A'+sigman2*eye(N);
			J_psi_psi_m = zeros(L,L,M);
			J_psi_s_m = zeros(L,L,M);
			J_psi_n_m = zeros(L,1,M);
			J_s_s_m = zeros(L,L,M);
			J_s_n_m = zeros(L,1,M);
			J_n_n_m = zeros(1,1,M);
    		for m = 1 : M
        		Bm = codebook(:,:,m);
        		Rm = Bm'*R*Bm;
        		inv_Rm = eye(Nrf)/Rm;
        		inv_Rm2 = inv_Rm*inv_Rm;
        		
        		inv_Rm = Bm*eye(Nrf)/Rm*Bm';
        		inv_Rm2 = Bm*inv_Rm2*Bm';
    		
        		D = 1j*diag(0:N-1)*A;
        		J_psi_psi_m(:,:,m) = 2*Km*real((Rs*A'*inv_Rm*D).*(Rs*A'*inv_Rm*D).' + (Rs*A'*inv_Rm*A*Rs).*(D'*inv_Rm*D).');
        		J_psi_s_m(:,:,m) = 2*Km*real((Rs*A'*inv_Rm*A).*(A'*inv_Rm*D).');
        		J_psi_n_m(:,:,m) = 2*Km*real(diag(real(Rs*A'*inv_Rm2*D)));
        		J_s_s_m(:,:,m) = Km*real((A'*inv_Rm*A).*(A'*inv_Rm*A).');
        		J_s_n_m(:,:,m) = Km*real(diag(A'*inv_Rm2*A));
        		J_n_n_m(:,:,m) = Km*real(trace(inv_Rm2));
    		end
    		
    		J_psi_psi = sum(J_psi_psi_m,3);
    		J_psi_s = sum(J_psi_s_m,3);
    		J_psi_n = sum(J_psi_n_m,3);
    		J_s_s = sum(J_s_s_m,3);
    		J_s_n = sum(J_s_n_m,3);
    		J_n_n = sum(J_n_n_m,3);
    		J = [J_psi_psi, J_psi_s, J_psi_n;
         		J_psi_s', J_s_s, J_s_n;
         		J_psi_n',J_s_n',J_n_n];
    		CRB = eye(2*L+1)/J;
		end
        
        function [CRB,J] = crb_unc_uncorrelated_theta(d,lambda,theta,Rs,sigman2,codebook,Km)
			[N,Nrf,M] = size(codebook);
            psi = 2*pi*d/lambda*sind(theta);
            A = ToolsObj.a_ula(N,psi);
			L = size(A,2);
    		R = A*Rs*A'+sigman2*eye(N);
			J_psi_psi_m = zeros(L,L,M);
			J_psi_s_m = zeros(L,L,M);
			J_psi_n_m = zeros(L,1,M);
			J_s_s_m = zeros(L,L,M);
			J_s_n_m = zeros(L,1,M);
			J_n_n_m = zeros(1,1,M);
    		for m = 1 : M
        		Bm = codebook(:,:,m);
        		Rm = Bm'*R*Bm;
        		inv_Rm = eye(Nrf)/Rm;
        		inv_Rm2 = inv_Rm*inv_Rm;

        		inv_Rm = Bm*eye(Nrf)/Rm*Bm';
        		inv_Rm2 = Bm*inv_Rm2*Bm';

        		D = 1j*pi*cosd(theta).*(diag(0:N-1)*A);
        		J_psi_psi_m(:,:,m) = 2*Km*real((Rs*A'*inv_Rm*D).*(Rs*A'*inv_Rm*D).' + (Rs*A'*inv_Rm*A*Rs).*(D'*inv_Rm*D).');
        		J_psi_s_m(:,:,m) = 2*Km*real((Rs*A'*inv_Rm*A).*(A'*inv_Rm*D).');
        		J_psi_n_m(:,:,m) = 2*Km*real(diag(real(Rs*A'*inv_Rm2*D)));
        		J_s_s_m(:,:,m) = Km*real((A'*inv_Rm*A).*(A'*inv_Rm*A).');
        		J_s_n_m(:,:,m) = Km*real(diag(A'*inv_Rm2*A));
        		J_n_n_m(:,:,m) = Km*real(trace(inv_Rm2));
    		end

    		J_psi_psi = sum(J_psi_psi_m,3);
    		J_psi_s = sum(J_psi_s_m,3);
    		J_psi_n = sum(J_psi_n_m,3);
    		J_s_s = sum(J_s_s_m,3);
    		J_s_n = sum(J_s_n_m,3);
    		J_n_n = sum(J_n_n_m,3);
    		J = [J_psi_psi, J_psi_s, J_psi_n;
         		J_psi_s', J_s_s, J_s_n;
         		J_psi_n',J_s_n',J_n_n];
    		CRB = eye(2*L+1)/J;
		end
        
        function CRB = eval_crb_theta_per_snr(d,lambda,theta,snr,codebook,Km)
            L = length(theta);
            CRB = zeros(L,length(snr));
            Rs = eye(L);

            for i = 1 : length(snr)
               sigman2 = 1/snr(i);
               aux = ToolsObj.crb_unc_uncorrelated_theta(d,lambda,theta,Rs,sigman2,codebook,Km);
               CRB(:,i) = diag(aux(1:L,1:L));
            end
		end

        function CRB = eval_crb_theta_per_Kmvec(d,lambda,theta,snr,codebook,Km_vec)
            L = length(theta);
            CRB = zeros(L,length(Km_vec));
            Rs = eye(L);
            sigman2 = 1/snr;
            for i = 1 : length(Km_vec)
               aux = ToolsObj.crb_unc_uncorrelated_theta(d,lambda,theta,Rs,sigman2,codebook,Km_vec(i));
               CRB(:,i) = diag(aux(1:L,1:L));
            end
		end

        function CRB = eval_crb_theta_per_theta(d,lambda,theta,snr,codebook,Km)
            CRB = zeros(1,length(theta));
            Rs = eye(1);
            sigman2 = 1/snr;
            for i = 1 : length(theta)
               aux = ToolsObj.crb_unc_uncorrelated_theta(d,lambda,theta(i),Rs,sigman2,codebook,Km);
               CRB(:,i) = aux(1,1);
            end
        end
		
        function [spec,theta_grid_deg] = eval_diagram(B,Q)
			n = size(B,1);
			theta_grid_deg = linspace(-85,85,Q);
			psi_grid = pi*sind(theta_grid_deg);
			a_grid = ToolsObj.a_ula(n,psi_grid);
			spec = abs(B'*a_grid).^2;
            figure
            plot(theta_grid_deg,spec)
        end
        
        function [X] = generate_uncorrelated_gaussian_channel(A,Rs,sigman2,K)
            [n,L] = size(A);
            s = sqrt(Rs/2)*(randn(L,K)+1j*randn(L,K));
            noise = sqrt(sigman2/2)*(randn(n,K)+1j*randn(n,K));
            X = A*s+noise;
        end

        function doas = get_doas_rootmusic(R,L)
            doas = asind(sort(rootmusic(R,L),'ascend')/pi).';
        end
    end
end
