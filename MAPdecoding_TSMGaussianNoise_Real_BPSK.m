function symbols_proba = MAPdecoding_TSMGaussianNoise_Real_BPSK(y,modulation_alphabet,mod_alphabet_apriori,trans_mat,init_state_prob,GaussPDFmean,GaussPDFvariance)

% [1] "Power Line Commucations" (p. 277-278) Ferreira, Lampe et al., Wiley 2010

% [2] "Channel Coding in Communication Networks. From Theory to Turbocodes", A. Glavieux et al., ISTE 2007

    N=size(trans_mat,1); % number of states
    
    L=length(y); M=length(modulation_alphabet);

    symbols_proba=zeros(M,L);

% Reshaping

    init_state_prob=reshape(init_state_prob,N,1);

    GaussPDFmean=reshape(GaussPDFmean,N,1);

    GaussPDFvariance=reshape(GaussPDFvariance,N,1);

    modulation_alphabet=reshape(modulation_alphabet,M,1);

    mod_alphabet_apriori=reshape(mod_alphabet_apriori,M,1);

% Forward-backward filters initialization

    alpha=zeros(N,L); alpha(:,1)=init_state_prob;

    beta=zeros(N,L+1); % only beta(:,2:L+1) elements are used in the algo

    beta(:,end)=1;

  % initial alpha and beta normalization in order to avoid recursive computation imprecisions [2]
    
    %alpha(:,1)=alpha(:,1)/sum(alpha(:,1));
    
    %beta(:,end)=beta(:,end)/sum(beta(:,end));
    
  % end normalization

% Forward-backward filters recursive computation

    temp_ones=ones(M,1);

    temp_GaussPDFvariance=kron(GaussPDFvariance,temp_ones);

    temp_GaussPDFmean=kron(GaussPDFmean,temp_ones);

    temp_modulation_alphabet=repmat(modulation_alphabet,N,1); % see if necessary

    temp_mod_alphabet_apriori=repmat(mod_alphabet_apriori,N,1);

    temp_previous_alpha=kron(alpha(:,1),temp_ones);

    next_beta=beta(:,end);

for k=1:L-1
        
        p_nk_alpha=(1./sqrt(2*pi*temp_GaussPDFvariance)).*exp(-((y(k)-temp_modulation_alphabet-temp_GaussPDFmean).^2)./(2*temp_GaussPDFvariance));
    
    for  s=1:N
    
        % Forward filter: alpha
        
        % Computation of F_k 
    
        % noise probability given each state and for each symbol,

        % Remark: The first M elements of p_nk_alpha correspond to the noise probability p_nk(n_k=y_k-c_k|s_k=s1) for each symbol c_k given the first state of the Markov chain 
        
        ps_alpha=trans_mat(:,s); % p(s_k+1|s_k) for each state s_k
        
        temp_ps_alpha=kron(ps_alpha,temp_ones);
        
        F_k_alpha=temp_ps_alpha.*p_nk_alpha; % the reading of F_k_alpha is similar to p_nk_alpha (see the above remark)
        
        temp_alpha=F_k_alpha.*temp_mod_alphabet_apriori.*temp_previous_alpha; % SEE IF KRON CAN BE USED TO COMPUTE temp_alpha LIKE FOR THE BACKWARD FILTER
        
        alpha(s,k+1)=sum(temp_alpha);
        
        % Backward filter: beta. SEE IF INTERESTING TO EXPLICITELY EXPRESS F_k
        
        p_nk_beta=(1/sqrt(2*pi*GaussPDFvariance(s)))*exp(-((y(L-k+1)-modulation_alphabet-GaussPDFmean(s)).^2)/(2*GaussPDFvariance(s)));
        
        temp1=p_nk_beta.*mod_alphabet_apriori;
                    
        ps_beta=trans_mat(s,:)';  % p(s_k+1|s_k) for each state s_k+1
        
        temp2=ps_beta.*next_beta;               
        
        temp_beta=kron(temp2,temp1);
        
        beta(s,L-k+1)=sum(temp_beta);      
        
    end
    
    % alpha and beta normalization in order to avoid recursive computation imprecisions [2]
    
    %alpha(:,k+1)=alpha(:,k+1)/sum(alpha(:,k+1));
    
    %beta(:,L-k+1)=beta(:,L-k+1)/sum(beta(:,L-k+1));
    
    sum_alpha=sum(alpha(:,k+1)); sum_beta=sum(beta(:,L-k+1));
    
    if sum_alpha==0 % situation where alpha(s,k+1)=0 for all states "s" due to numerical imprecisions (mainly from exponential computation in p_nk)
        
        alpha(:,k+1)=1/N;
        
    else
        
        alpha(:,k+1)=alpha(:,k+1)/sum_alpha; 
        
    end
      
    if sum_beta==0
        
        beta(:,L-k+1)=1/N;
        
    else
        
        beta(:,L-k+1)=beta(:,L-k+1)/sum_beta;
        
    end
    
    % end normalization
    
    temp_previous_alpha=kron(alpha(:,k+1),temp_ones);
    
    next_beta=beta(:,L-k+1);

    end
    
for k=1:L
    alpha_replicate=repmat(alpha(:,k),1,N);
        
    beta_replicate=repmat(beta(:,k+1)',N,1);
        
    for m=1:M
        
        p_nk=(1./sqrt(2*pi*GaussPDFvariance)).*exp(-((y(k)-modulation_alphabet(m)-GaussPDFmean).^2)./(2*GaussPDFvariance));
        
        p_nk_replicate=repmat(p_nk,1,N);
        
        F_k=p_nk_replicate.*trans_mat;
        
        temp=F_k.*alpha_replicate.*beta_replicate;
        
        symbols_proba(m,k)=mod_alphabet_apriori(m)*sum(sum(temp)); % represents p(c_k,y_k), which is proportional to p(c_k|y_k)
        
        % Note: sum(p(c_k=+/-1,y_k)) is not equal to 1, whereas sum(p(c_k=+/-1|y_k))=1
    
    end   
    sum_symb_proba=sum(symbols_proba(:,k));
    
    if sum_symb_proba==0 % means symbols_proba(:,k) is a null vector. This can arise from numerical imprecision of the exponential computation in p_nk
        
        symbols_proba(:,k)=1/M; 
        
    end
    
end

return;
