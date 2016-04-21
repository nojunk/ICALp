clear all;
M_Lp_ICA_G = 0;
M_Lp_ICA_F = 0;
M_PCA_L1 = 0;
M_f = 0;

Arr_LG = zeros(100,5);
Arr_LF = zeros(100,5);
Arr_LE = zeros(100,5);
Arr_L1 = zeros(100,5);
Arr_f = zeros(100,5);
Arr_ef = zeros(100,5);
Arr_1f = zeros(100,5);

Time_LG = zeros(100,5);
Time_LF = zeros(100,5);
Time_LE = zeros(100,5);
Time_L1 = zeros(100,5);
Time_f = zeros(100,5);
Time_ef = zeros(100,5);
Time_1f = zeros(100,5);

rf = 0;
%num of sources
D_sub = 3;
D_sup = 3;
%D = D_sub+D_sup;
%num of samples
N = 500;

for setting=1:6
    
    %{
    switch setting
    %number of samples
    case 1
        N = 100;
        D_sub = 3;
        D_sup = 3;
    case 2
        N = 300;
        D_sub = 3;
        D_sup = 3;
    case 3
        N = 500;
        D_sub = 3;
        D_sup = 3;
    case 4
        N = 700;
        D_sub = 3;
        D_sup = 3;
    case 5
        N = 900;
        D_sub = 3;
        D_sup = 3;
    end
    %}
    %{
    %number of sources
    switch setting
    case 1
        N = 500;
        D_sub = 1;
        D_sup = 1;
    case 2
        N = 500;
        D_sub = 3;
        D_sup = 3;
    case 3
        N = 500;
        D_sub = 5;
        D_sup = 5;
    case 4
        N = 500;
        D_sub = 7;
        D_sup = 7;
    case 5
        N = 500;
        D_sub = 9;
        D_sup = 9;
    %gaussian noise
    end
    %}
    D = D_sub+D_sup;
    %{
    case 11
        N = 500;
        D_sub = 3;
        D_sup = 3;
    case 12
        N = 500;
        D_sub = 3;
        D_sup = 3;
    case 13
        N = 500;
        D_sub = 3;
        D_sup = 3;
    %uniform noise    
    case 14
        N = 500;
        D_sub = 3;
        D_sup = 3;
    case 15
        N = 500;
        D_sub = 3;
        D_sup = 3;
    case 16
        N = 500;
        D_sub = 3;
        D_sup = 3;
end
    %}
    %{
    switch setting
        case 1
            N_gauss = N*0.1;
        case 2
            N_gauss = N*0.2;
        case 3
            N_gauss = N*0.3;
        case 4
            N_gauss = N*0.4;
        case 5
            N_gauss = N*0.5;
    end
    %}
    %D = D_sub +D_sup;
    N_gauss = 0;
for r=1:100
    %% Data generation
    if D_sub == 0
        S0 = laprnd(D_sup,N-N_gauss);  % super-G: Laplacian, p = 1
    elseif D_sup == 0
        S0 = (rand(D_sub,N-N_gauss)-0.5)*sqrt(3)*2;    % sub-G : uniform, p = infty
    else
        S1 = (rand(D_sub,N-N_gauss)-0.5)*sqrt(3)*2;    % sub-G : uniform, p = infty
        S2 = laprnd(D_sup,N-N_gauss);  % super-G: Laplacian, p = 1
        S0 = [S1;S2];
    end
    
    %noise generation
    
    switch setting
        case 1
            noise = randn(D,N)*0.001;
        case 2
            noise = randn(D,N)*0.01;
        case 3
            noise = randn(D,N)*0.1;
        case 4
            noise = (rand(D,N)-0.5)*0.001;
        case 5
            noise = (rand(D,N)-0.5)*0.01;
        case 6
            noise = (rand(D,N)-0.5)*0.1;
    end
    
    %gaussian source ratio
     Ga = randn(D,N_gauss);
     S0 = [S0,Ga];
    
    S = zeros(size(S0));
    % making zerom mean unit variance
    for i=1:D
        S(i,:) = (S0(i,:)-mean(S0(i,:))*ones(1,N))/std(S0(i,:));
    end
    S = S+noise;
    
    % mixing
    A0 = rand(D,D)-0.5;
    X0 = A0*S;
    %{
    %additive noise if needed
    if setting < 11
        X0 = A0*S;% + noise;
    else
        X0 = A0*S + noise;
    end
    %}
    %% Whitening
    [U,L,V] = svd(X0);
    sv = diag(L);
    invL = diag(1./sv)*sqrt(N);
    
    X = invL*U'*X0;
    A = invL*U'*A0;

    %% Lp_ICA_G
    tic;
    W_Lp_ICA_G = Lp_ICA_G(X,D,N);
    elapsedTime = toc;
    Time_LG(r,setting) = elapsedTime;
    
    %% Lp_ICA_exact (p is also estimated)
    tic;
    W_Lp_ICA_exact = Lp_ICA_exact(X,D);
    elapsedTime = toc;
    Time_LE(r,setting) = elapsedTime;
    
    %% Lp_ICA_F
    tic;
    W_Lp_ICA_F = Lp_ICA_F(X,D_sub,D_sup,N);
    elapsedTime = toc;
    Time_LF(r,setting) = elapsedTime;
    
    %% PCA-L1
    if D_sup == 0
        tic;
        W_PCA_L1 = ICA_PCA_L1(X,D,N);
        elapsedTime = toc;
        Time_L1(r,setting) = elapsedTime;
    end
    
    %% FastICA    
    tic;
    [W_fast, A_fast] = fastica(X,'maxNumIterations',2000, 'epsilon', 0.001, 'verbose', 'off', 'displayMode', 'off' );
    elapsedtime = toc; 
    Time_f(r,setting) = elapsedtime;
    
    %% EFICA    
    tic;
    W_ef = efica(X);
    W_ef = W_ef';
    elapsedtime = toc; 
    Time_ef(r,setting) = elapsedtime;
    
    %% 1FICA    
    tic;
    W_1f = FICA1(X);
    W_1f = W_1f';
    elapsedtime = toc; 
    Time_1f(r,setting) = elapsedtime;
    
    
    %% measure      
    
    %%%%%%%%%%%%%%%%%% for Lp-ICA-G
    % ideally, B should be a permutation matrix
    B = W_Lp_ICA_G'*A;
    
    P = pol_n_permute(B);
    
    % estimated source with an appropriate order
    Y = P'*W_Lp_ICA_G'*X;
    
    % mean squared error
    Arr_LG(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
    %M_Lp_ICA_G = M_Lp_ICA_G+mse_Lp_ICA_G;
    
    disp('Lp_ICA gradient descent');
    disp([Arr_LG(r,setting), Time_LG(r, setting)]);
    
    %%%%%%%%%%%%%%%%%% for Lp-ICA-exact
    % ideally, B should be a permutation matrix
    B = W_Lp_ICA_exact'*A;
    
    P = pol_n_permute(B);
    
    % estimated source with an appropriate order
    Y = P'*W_Lp_ICA_exact'*X;
    
    % mean squared error
    Arr_LE(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
    %M_Lp_ICA_G = M_Lp_ICA_G+mse_Lp_ICA_G;
    
    disp('Lp_ICA exact gradient descent');
    disp([Arr_LE(r,setting), Time_LE(r, setting)]);
    
    
    %%%%%%%%%%%%%%%%%% for Lp-ICA-F
    % ideally, B should be a permutation matrix
    B = W_Lp_ICA_F'*A;
    
    P = pol_n_permute(B);
    % estimated source with an appropriate order
    Y = P'*W_Lp_ICA_F'*X;
    
    % mean squared error
    Arr_LF(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
    %M_Lp_ICA_F = M_Lp_ICA_F+mse_Lp_ICA_F;
    
    disp('Lp_ICA fast');
    disp( [Arr_LF(r,setting), Time_LF(r, setting)]);
    
    %%%%%%%%%%%%%%%%%% for PCA-L1
    if D_sup == 0
        % ideally, B should be a permutation matrix
        B = W_PCA_L1'*A;
        
        P = pol_n_permute(B);
        % estimated source with an appropriate order
        Y = P'*W_PCA_L1'*X;
        
        % mean squared error
        Arr_L1(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
        %M_PCA_L1 = M_PCA_L1+mse_PCA_L1;
        
        disp('PCA-L1');
        disp([Arr_L1(r,setting), Time_L1(r, setting)]);
    end
    
    %%%%%%%%%%%%%%%%%%% for FASTICA
    [rows,cols] = size(W_fast);
    
    if rows == cols
        B = W_fast'*A;
        P = pol_n_permute(B);
        Y = P'*W_fast'*X;
        Arr_f(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
        %M_f = M_f+mse_fast;
        %rf = rf+1;
        disp('fastICA');
        disp([Arr_f(r,setting),Time_f(r, setting)]);
    else
        Arr_f(r,setting) = NaN;
    end
    
     %%%%%%%%%%%%%%%%%%% for EFICA
    [rows,cols] = size(W_ef);
    
    if rows == cols
        B = W_ef'*A;
        P = pol_n_permute(B);
        Y = P'*W_ef'*X;
        Arr_ef(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
        %M_f = M_f+mse_fast;
        %rf = rf+1;
        disp('EFICA');
        disp([Arr_ef(r,setting),Time_ef(r, setting)]);
    else
        Arr_ef(r,setting) = NaN;
    end
    
      %%%%%%%%%%%%%%%%%%% for 1FICA
    [rows,cols] = size(W_1f);
    
    if rows == cols
        B = W_1f'*A;
        P = pol_n_permute(B);
        Y = P'*W_1f'*X;
        Arr_1f(r,setting) = norm(S-Y,'fro')/sqrt(D*N);
        %M_f = M_f+mse_fast;
        %rf = rf+1;
        disp('1FICA');
        disp([Arr_1f(r,setting),Time_1f(r, setting)]);
    else
        Arr_1f(r,setting) = NaN;
    end
    %% end of measurement
    disp([num2str(setting),' : ',num2str(r),' trial done']);
end
%{
Arr_LG(setting) = M_Lp_ICA_G/r;
Arr_LF(setting) = M_Lp_ICA_F/r;
Arr_f(setting) = M_f/rf;

M_Lp_ICA_G = 0;
M_Lp_ICA_F = 0;
M_f = 0;
rf = 0;
%}
end

save('gauss_additive_noise.mat', ...
    'Arr_1f', 'Time_1f', ...
    'Arr_ef', 'Time_ef', ...
    'Arr_f', 'Time_f', ...
    'Arr_LF', 'Time_LF',...
    'Arr_LE', 'Time_LE', ...
    'Arr_LG', 'Time_LG', ...
    'Arr_L1', 'Time_L1');