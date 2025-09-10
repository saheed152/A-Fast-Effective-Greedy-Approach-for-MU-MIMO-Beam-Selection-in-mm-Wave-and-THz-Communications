clear,clc
close all
addpath(pwd);
%cd manopt;
addpath(genpath(pwd));
%cd ..;
Ns = 2; %number of data streams to be transmitted
NRF =2;  %This is Nrf for one user. actual NRF of the system = NRF*user
User=5;
SNR_dB = -15:5:25;
SNR = 10.^(SNR_dB./10);
smax = length(SNR);% enable the parallel
realiza=1;
In_PE=0;
%WRF=1
power = .1;
Nr=36;
Nt=144;
% n = wgn(Nr,1,power); %AWGN Noise
% n1 = wgn(Ns,1,power);

load("awgn_noise.mat");


%load('../../datasets/all_u5_ns2_1.mat'); %Load the channnel matrix

tEnd = 0;
r_sel = 4;
Comb_Fopt = zeros(Nt,User*r_sel);
Comb_Wopt = zeros(Nr,User*r_sel);
Fopt_corr = zeros(Nt,User*r_sel);
Comb_sv = zeros(User*r_sel);
R_sumo = zeros(length(SNR),1);
Max_R = zeros(length(SNR),1);
load('all_u5_ns2.mat');
tStart = tic;
for tl = 1:10

    for i=1:User
        [U,S,V] = svd(H(:,:,i));
        Comb_Fopt(:,(i-1)*r_sel+1:i*r_sel) = V(:,1:r_sel);
        Comb_Wopt(:,(i-1)*r_sel+1:i*r_sel) = U(:,1:r_sel);
        temp_F = V*S';
        Fopt_corr(:,(i-1)*r_sel+1:i*r_sel) = temp_F([1:Nt],[1:r_sel]); % Getting V multiplied by sigmas
        temp_sv = diag(S);
        Comb_sv((i-1)*r_sel+1:i*r_sel) = temp_sv(1:r_sel);
    end

    big_Fopt=Fopt_corr'*Fopt_corr;
    u = User;
    r = r_sel;
    Ns_red = 1;
    A = nchoosek([1:r],1);
    num_comb = size(A,Ns_red);  % This is the number of combinations of one user's vectors (rCNs)

    allPerms = permn([1:num_comb], User);  % Generate all permutations. Function 'permn' needs to be available in your path.

    ind = zeros(size(allPerms,1), User*Ns_red);
    for i = 1:size(allPerms,1)
        ind(i,:) = reshape(A(allPerms(i,:),:)',1, User*Ns_red);
    end

    % The following loop sets indices for the combined matrix
    for i = 1:User-1
        ind(:, Ns_red*i+1:Ns_red*(i+1)) = ind(:, Ns_red*i+1:Ns_red*(i+1)) + i*r;
    end
    ind_tot = [];
    for k = 1:Ns


        Min_corr = 10^9;
        max_sv = sum(Comb_sv(ind(1,:)));
        size_big=r_sel*User;
        kount1=0;
        for j = 1:size(ind,1)

            if (sum(Comb_sv(ind(j,:)))>0.65*max_sv)
                kount1 = kount1 + 1;
                temp_corr_sel = sel_matrix(ind(j,:),size_big).*big_Fopt;
                temp_diag = temp_corr_sel - diag(diag(temp_corr_sel));
                Corr_F = norm(temp_diag,'fro');

                if ((Corr_F < Min_corr))
                    Min_corr = Corr_F;
                    ind_corr_min = j;
                end
            end

        end
        ind_tot = [ind_tot,ind(ind_corr_min,:)];

        ind = removeRowsWithValues(ind, ind(ind_corr_min,:));

    end


    ind_tot = sort(ind_tot);

    Fopt = Comb_Fopt(:,ind_tot);

    
end


tEnd = tEnd + toc(tStart);
Corr_F = norm(abs(Fopt'*Fopt),'fro')^2 - trace((Fopt'*Fopt).^2);

Fopt = reshape(Fopt,[Nt,Ns,User]);


Wopt = Comb_Wopt(:,ind_tot);
Wopt = reshape(Wopt,[Nr,Ns,User]);


for s = 1:smax
    R_o=0;
    SNR(s)
    for U=1:User
        R_opt=(SNR(s)/User)*pinv(Wopt(:,:,U)) * H(:,:,U) * Fopt(:,:,U) * Fopt(:,:,U)' * H(:,:,U)' * Wopt(:,:,U);
        In_opt=0;
        for p=1:User
            if p~=U
                In_opt=In_opt+(SNR(s)/User)*pinv(Wopt(:,:,U)) * H(:,:,U) * Fopt(:,:,p) * Fopt(:,:,p)' * H(:,:,U)' * Wopt(:,:,U)+Wopt(:,:,U)'* n*n'*Wopt(:,:,U);
            end
        end
        R_o=R_o+ log2(det(eye(Ns) +R_opt/In_opt));
    end
    R_sumo(s)=R_o; %Suboptimal
end

plot(SNR_dB,R_sumo,'k-->','Marker','>','LineWidth',1);

xlabel('SNR [dB]')
ylabel('Spectral efficiency (bits/s/Hz)')
legend('IOSVB_greedy','Location','NW')
