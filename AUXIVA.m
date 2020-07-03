%% var_opt == 0 : Spherical laplacian distribution
%% var_opt == 1 : Multivariate gaussian distribution with Time-varying Variance
function [y] = AUXIVA(x, nfft, var_opt, Maxiter)
nol = fix(3*nfft/4);
nshift = nfft-nol;

%% STFT
X = STFT( x, nfft, nshift);
[K, N, M] = size(X);

%% intialization
R = zeros(1,N,M);
V = zeros(M, M, K); % weighted covariance matrix
Y = zeros(K,N,M);
W = zeros(M, M, K); % demixing matrix
normCoef = zeros(1,K);

for k = 1 : K
    W(:,:,k) = eye(M); 
    Y(k,:,:) = (W(:,:,k)*squeeze(X(k,:,:)).').'; 
end
P = max(abs(Y).^2,eps);
Xp = permute(X,[3,2,1]);

%% AUXIVA processing
fprintf('Iteration:    ');
for iter = 1 : Maxiter
    fprintf('\b\b\b\b%4d', iter);    
    for m  = 1 : M
        %% obtain output Y = W * X
        for k = 1:K
            Y(k,:,m) = W(m,:,k)*Xp(:,:,k);
        end
        P(:,:,m) = max(abs(Y(:,:,m)).^2,eps);
        XpHt = conj( permute( Xp, [2,1,3] ) );
        if var_opt == 0
            R(:,:,m) = sqrt(sum(P(:,:,m),1));
        else
            R(:,:,m) = sum(P(:,:,m),1)/K;
        end
        Rp = permute(R,[3,2,1]);
        for k = 1 : K
           V(:,:,k) = (Xp(:,:,k).*max(Rp(m,:),10^-6))*XpHt(:,:,k)/N; 
           invWD(:,:,k) = inv(W(:,:,k)*V(:,:,k));
           invWDE(:,k) = squeeze(invWD(:,m,k));
           normCoef(1,k) = invWDE(:,k)'*V(:,:,k)*invWDE(:,k);
        end
        W(m, :,:) = conj(invWDE./max(sqrt(normCoef),eps));     
    end    
end

%% MDP
for k =  1:K
    Y(k,:,:) = (diag(diag(inv(W(:,:,k))))*W(:,:,k)*Xp(:,:,k)).';
end
Y(isnan(Y)) = eps ;
y = ISTFT( Y, nfft, nshift, size(x,1) );
end

