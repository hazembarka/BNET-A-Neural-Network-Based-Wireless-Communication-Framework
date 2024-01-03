% Communication chain
% Communication chain

clear all; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555555

% Parameters
snr_db     =[6]; snr_lin=10.^(snr_db./10); snr_size=length(snr_lin);
%snr_db     =[8]; snr_lin=10.^(snr_db./10); snr_size=length(snr_lin);
R_list     =[50];
gamma_list =[100];
p_B_list = [0.1];

frameNb = 256; % number of transmitted frames. the length of a frame is given by the number of encoder input bits (e.g. 32400 bits for the DVB-S2 LDPC encoder)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555555

load('all_MAP','enc','dec');
MAPdetectorInputLength=64800; % number of input symbols of the detector. should be a factor of L
N_sample=frameNb*MAPdetectorInputLength;

% Mapping parameters
modulation_alphabet=[1 -1]; 
mod_alphabet_apriori=[.5 .5]; % QPSK Mapping. mod_alphabet_apriori: a priori probabilities of the mapping elements


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555555


% BER evaluation




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555555
for P_B_id  = 1:length(p_B_list)
    p_B = p_B_list(P_B_id)

for gamma_id  = 1:length(gamma_list)
    gamma = gamma_list(gamma_id)
    
for R_id  = 1:length(R_list)
    R = R_list(R_id)



    %Two state Markov-Gaussian derived parameters
    p_G=1-p_B;
    init_state_prob=[p_G p_B]; %(p_GB+p_BG)= 0.1;
    p_BG=(1/gamma)*p_G; p_GB=(1/gamma)*p_B;
    p_BB=1-p_BG;

    trans_mat=[1-p_GB p_GB; p_BG 1-p_BG];
    GaussPDFmean=[0 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55555555
trans_signal_power=1;


for snr=1:snr_size

 chSNR=snr_lin(snr);
 snr_db_value = num2str(snr_db(snr));
 snr_db_value
noise_power_input=trans_signal_power/chSNR;
s_0=sqrt(noise_power_input);
s_1=sqrt((R)*noise_power_input);

all_snr = [];
all_R = [];
all_p_B = [];
all_gamma = [];
all_y_SD = [];
all_Optimum_LLR = [];
all_MAP = [];
all_bits = [];
all_decoded_original = [];
T = table(all_snr,all_y_SD,all_MAP);
T_original_bits = table(all_snr,all_bits,all_decoded_original);


    for nb=1:frameNb

    % Coding
    len = 32400*2;
    source_bits=randi([0 1],32400,1); % enc.NumInfoBits size is 32400 for the default LDPC encoder
    encodedBits = enc(source_bits); % LDPC encode
    x = 1-2*encodedBits; %BPSK mapping
    x = sqrt(trans_signal_power)*x; % the power of emit_signal is trans_signal_power for BPSK
    modulation_alphabet_current=sqrt(trans_signal_power)*modulation_alphabet;
    x = reshape(x,[len,1]);


    %%Generating Noise samples
    corrNoise_SD = TSMG(p_B,gamma,s_0,s_1,length(x));
    corrNoise_SD = reshape(corrNoise_SD,[len,1]);



% Generate Rayleigh bloc fading Channel, with Tc = 32 Ts
    Tc = 32;
    h = sqrt(0.5)*(randn(len/Tc,1)+ 1j*randn(len/Tc,1));
    h = repelem(h,Tc);

    %%The receiver Signal
    y    = h.*x  + corrNoise_SD;

    % equalization of received data by channel information at the receiver
    y_SD = real(y./h);   
    
       
    % BER for a MAP detection
    
    MAPdetectorInputLength = Tc;

    sigma_A=noise_power_input; sigma_B=noise_power_input*(R);

    
    kk_max=length(x)/MAPdetectorInputLength-1; 
    metric_ratio_MAP_SD=zeros(1,length(x));
    metric_ratio_iid_SD=zeros(1,length(x));
    % Other noise derived parameters

for kk=0:kk_max
    
            vect_temp=kk*MAPdetectorInputLength+1:(kk+1)*MAPdetectorInputLength;
            mean_h2 = mean(abs(h(vect_temp)).^2);
            GaussPDFvariance =[sigma_A/mean_h2 sigma_B/mean_h2];
            y_SD_detect=y_SD(vect_temp);
            symbols_proba_SD_MAP = MAPdecoding_TSMGaussianNoise_Real_BPSK(y_SD_detect,modulation_alphabet_current,mod_alphabet_apriori,trans_mat,init_state_prob,GaussPDFmean,GaussPDFvariance);
            metric_ratio_MAP_SD(vect_temp)=(log(symbols_proba_SD_MAP(1,:)./symbols_proba_SD_MAP(2,:)));
end
            
    metric_ratio_MAP_SD = reshape(metric_ratio_MAP_SD,[len,1]);
    all_y_SD = round(y_SD,3);
    all_MAP = round(metric_ratio_MAP_SD,3);
    dec_bits_di_MAP = dec(metric_ratio_MAP_SD);
    decoded = dec(metric_ratio_MAP_SD);
    
    all_bits =  source_bits;
    all_decoded_original =  decoded;
    all_snr = ones(64800,1)*snr_db(snr);
    
    T__ = table(all_snr,all_y_SD,all_MAP);
    T = [T;T__];
    all_snr = ones(32400,1)*snr_db(snr);
    T__original_bits = table(all_snr,all_bits,all_decoded_original);
    T_original_bits = [T_original_bits;T__original_bits];
    end
    writetable(T_original_bits,['test datasets/G',num2str(gamma),'_R',num2str(R),'/',num2str(snr_db_value),' original bits.csv']);
    writetable(T,['test datasets/G',num2str(gamma),'_R',num2str(R),'/',num2str(snr_db_value),'received.csv']);
 
end
end
end
end
% 