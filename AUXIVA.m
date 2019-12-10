%% var_opt == 0 : Spherical laplacian distribution
%% var_opt == 1 : Multivariate gaussian distribution with Time-varying Variance
function [y] = AUXIVA(x, nfft, var_opt, Maxiter)
nol = fix(3*nfft/4);
nshift = nfft-nol;

%% STFT
X_tmp = STFT( x, nfft, nshift);
X = X_tmp;

[K, Nframe, M] = size(X);
% tic

%% intialization
R = zeros(1,Nframe,M);
V = zeros(M, M, K); % weighted covariance matrix
Y = zeros(K,Nframe,M);
W = zeros(M, M, K); % demixing matrix

for k = 1 : K
    W(:,:,k) = eye(M); 
    Y(k,:,:) = (W(:,:,k)*squeeze(X(k,:,:)).').'; 
end

P = max(abs(Y).^2,eps); % ����� �Ŀ��� P�� �Ҵ�
Xp = permute(X,[3,2,1]); % Nmic x Nframe x Nfreq

%% AUXIVA processing
fprintf('Iteration:    ');
for iter = 1 : Maxiter
    fprintf('\b\b\b\b%4d', iter);    
    for m  = 1 : M
        XpHt = conj( permute( Xp, [2,1,3] ) ); % Nframe x Nmic x Nfreq (matrix-wise Hermitian transpose)
        if var_opt == 0
            R(:,:,m) = sqrt(sum(P(:,:,m),1));
        else
            R(:,:,m) = sum(P(:,:,m),1)/K;
        end
        Rp = permute(R,[3,2,1]); % Nmic x Nframe
        eleinvRp =repmat(1./Rp(m,:),M,1,1);
        for k = 1 : K
           phi = eleinvRp(:,:,1);
           V(:,:,k) = (Xp(:,:,k).*phi)*XpHt(:,:,k)/Nframe; % Nmic x Nmic x Nfreq
           invWD(:,:,k) = inv(W(:,:,k)*V(:,:,k));
           invWDE(:,k) = squeeze(invWD(:,m,k));
        end
        %% normalizing (37)               
        normCoef = zeros(1,K);
        for k = 1 : K
            normCoef(1,k) = invWDE(:,k)'*V(:,:,k)*invWDE(:,k);
        end
        w = (invWDE./max(sqrt(normCoef),eps)).'; 
        
        W(m, :,:) = w';
        %% obtain output Y = W * X
        for k = 1:K
            Y(k,:,m) = W(m,:,k)*Xp(:,:,k);
        end         
    end    
    P = max(abs(Y).^2,eps);
end

%% MDP
for k =  1:K    
    Wmdp(:,:,k) = diag(diag(inv(W(:,:,k))))*W(:,:,k);
    Y(k,:,:) = (Wmdp(:,:,k)*Xp(:,:,k)).';
end

Y(isnan(Y)) = eps ;
y = ISTFT( Y, nfft, nshift, size(x,1) );

end

