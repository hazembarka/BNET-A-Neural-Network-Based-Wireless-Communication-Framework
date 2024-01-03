clear all; close all;

% Parameters
snr_db = 0:2:10; snr_lin=10.^(snr_db./10); snr_size=length(snr_lin);

frameNb = 2; % number of transmitted frames. the length of a frame is given by the number of encoder input bits (e.g. 32400 bits for the DVB-S2 LDPC encoder)

all_y_SD = [];
all_MAP = [];

enc = comm.LDPCEncoder;
dec = comm.LDPCDecoder;

dec.DecisionMethod = 'Hard decision';

dec.OutputValue = 'Information part';

dec.MaximumIterationCount = 50;

% Stop if all parity-checks are satisfied

dec.FinalParityChecksOutputPort  = true;


MAPdetectorInputLength=1000000; % number of input symbols of the detector. should be a factor of L

N_sample=frameNb*MAPdetectorInputLength;

% Mapping parameters

modulation_alphabet=[1 -1]; 

mod_alphabet_apriori=[.5 .5]; % QPSK Mapping. mod_alphabet_apriori: a priori probabilities of the mapping elements

% Two-state Markov-Gaussian global parameters

p_B=0.1; gamma = 100 ; R = 50; p_B_Gaussian=0;% values taken from fertonani paper
%Two state Markov-Gaussian derived parameters

p_G=1-p_B;

init_state_prob=[p_G p_B]; %(p_GB+p_BG)= 0.1;

p_BG=(1/gamma)*p_G; p_GB=(1/gamma)*p_B;

p_BB=1-p_BG;

trans_mat=[1-p_GB p_GB; p_BG 1-p_BG];

GaussPDFmean=[0 0];

% BER evaluation

trans_signal_power=1;


error_Direct_MAP=zeros(1,snr_size); 
error_Direct_MAP_h=zeros(1,snr_size);
error_Direct_blank_fixed=zeros(1,snr_size);
error_Direct_clipped_fixed=zeros(1,snr_size);
error_Direct_clip_blank_fixed=zeros(1,snr_size);
error_Direct_Gaussian=zeros(1,snr_size);
error_Direct_AWGN=zeros(1,snr_size);


for snr=1:snr_size
    
snr_db(snr)
    
chSNR=snr_lin(snr);    
noise_power_input=trans_signal_power/chSNR;
s_0=sqrt(noise_power_input);
s_1=sqrt(R*noise_power_input);

for nb=1:frameNb
    
% Coding
len = 2000000;
source_bits=randi([0 1],len,1); % enc.NumInfoBits size is 32400 for the default LDPC encoder
encodedBits = source_bits;%enc(source_bits); % LDPC encode
x = 1-2*encodedBits; %BPSK mapping
x = sqrt(trans_signal_power)*x; % the power of emit_signal is trans_signal_power for BPSK
modulation_alphabet_current=sqrt(trans_signal_power)*modulation_alphabet;
x = reshape(x,[len,1]);


%%Generating Noise samples
corrNoise_SD = TSMG(p_B,gamma,s_0,s_1,length(x));
corrNoise_SD = reshape(corrNoise_SD,[len,1]);

AWGN_SD      = TSMG(p_B_Gaussian,gamma,s_0,s_0,length(x));
AWGN_SD      = reshape(AWGN_SD,[len,1]);


% Generate Rayleigh bloc fading Channel, with Tc = 32 Ts
Tc = 64800;
h = sqrt(0.5)*(randn(len,1)+ 1j*randn(len,1));
%h = repelem(1,64800);
%h = reshape(h,[len,1]);
%%The receiver Signal
y          = x  + corrNoise_SD;

y_AWGN     = x   + AWGN_SD;



tic
gamma = 100;
R = 50;
p_B = 0.1;
p_BG=(1/gamma)*p_G; p_GB=(1/gamma)*p_B;
p_BB=1-p_BG;
trans_mat=[1-p_GB p_GB; p_BG 1-p_BG];
sigma_A=noise_power_input; sigma_B= noise_power_input*(R);

x_rec = real(y./h);   

%% BER for a MAP detection
MAPdetectorInputLength = 2000000;
kk_max=length(x)/MAPdetectorInputLength-1; 
metric_ratio_MAP_SD=zeros(1,length(x));
x_rec_GR=zeros(1,length(x));
trans_mat = ones(19)/19;
init_state_prob =  ones(1,19)/19;
GaussPDFmean =  ones(1,19)/19;
GaussPDFvariance =  ones(1,19)/19;

for kk=0:kk_max
    
            %%MAP decoding Part
            vect_temp=kk*MAPdetectorInputLength+1:(kk+1)*MAPdetectorInputLength;
            mean_h2 = mean(abs(h(vect_temp)).^2);
            x_rec_detect=x_rec(vect_temp);

symbols_proba_SD_MAP = MAPdecoding_TSMGaussianNoise_Real_BPSK(x_rec_detect,modulation_alphabet_current,mod_alphabet_apriori,trans_mat,init_state_prob,GaussPDFmean,GaussPDFvariance);
metric_ratio_MAP_SD(vect_temp)=(log(symbols_proba_SD_MAP(1,:)./symbols_proba_SD_MAP(2,:)));
end
toc


% equalization of received data by channel information at the receiver

%ZF
%x_rec = real(y./h);   
%x_rec_AWGN= real(y_AWGN./h);


% % LMMSE
% Pt = trans_signal_power; 
% h_2 = abs(h).^2;
% x_rec = real((Pt.*conj(h)./(Pt.*h_2 + s_0^2)).*y);%




%noise derived parameters
gamma = 100;
R = 50;
p_B = 0.1;
p_BG=(1/gamma)*p_G; p_GB=(1/gamma)*p_B;
p_BB=1-p_BG;
trans_mat=[1-p_GB p_GB; p_BG 1-p_BG];
sigma_A=noise_power_input; sigma_B= noise_power_input*(R);




%% Clipping and Blanking performance with fixed threshold
x_rec_abs=abs(x_rec);
T_clip_fixed=0.9;
T_blank_fixed=3.4;


%Blanking performance with fixed threshold
x_rec_blanked_fixed=zeros(1,length(x));
for t=1:length(x);
    if x_rec_abs(t)<=T_blank_fixed
        x_rec_blanked_fixed(t)=x_rec(t);
    elseif x_rec_abs(t)>T_blank_fixed
        x_rec_blanked_fixed(t)=0;
    end
end
x_rec_blanked_fixed = reshape(x_rec_blanked_fixed,[len,1]);


%%Clipping performance with fixed threshold
x_rec_clipped_fixed=zeros(1,length(x));
for t=1:length(x);
    if x_rec_abs(t)<=T_clip_fixed
        x_rec_clipped_fixed(t)=x_rec(t);
    elseif x_rec_abs(t)>T_clip_fixed
        x_rec_clipped_fixed(t)=T_clip_fixed*sign(x_rec(t));
    end
end
x_rec_clipped_fixed = reshape(x_rec_clipped_fixed,[len,1]);



%%Combined clipping&Blanking performance with fixed threshold
T_com_clip=T_clip_fixed;
T_com_blank=T_blank_fixed;
x_rec_clip_blank_fixed=zeros(1,length(x));
for t=1:length(x)
    
    if x_rec_abs(t)<=T_com_clip
       x_rec_clip_blank_fixed(t)=x_rec(t);
    elseif(x_rec_abs(t)>T_com_clip) && (x_rec_abs(t)<=T_com_blank)
         x_rec_clip_blank_fixed(t)=T_com_clip*sign(x_rec(t));
    elseif x_rec_abs(t)>T_com_blank
        x_rec_clip_blank_fixed(t)=0;
    end
end
x_rec_clip_blank_fixed = reshape(x_rec_clip_blank_fixed,[len,1]);



tic

%% BER for a MAP detection
MAPdetectorInputLength = 1000000;
kk_max=length(x)/MAPdetectorInputLength-1; 
metric_ratio_MAP_SD=zeros(1,length(x));
x_rec_GR=zeros(1,length(x));


for kk=0:kk_max
    
            %%MAP decoding Part
            vect_temp=kk*MAPdetectorInputLength+1:(kk+1)*MAPdetectorInputLength;
            mean_h2 = mean(abs(h(vect_temp)).^2);
            GaussPDFvariance =[sigma_A/mean_h2 sigma_B/mean_h2];
            x_rec_detect=x_rec(vect_temp);

symbols_proba_SD_MAP = MAPdecoding_TSMGaussianNoise_Real_BPSK(x_rec_detect,modulation_alphabet_current,mod_alphabet_apriori,trans_mat,init_state_prob,GaussPDFmean,GaussPDFvariance);
metric_ratio_MAP_SD(vect_temp)=(log(symbols_proba_SD_MAP(1,:)./symbols_proba_SD_MAP(2,:)));

toc

%%Rayleigh + AWGN Channel Decoding
% x_rec_AWGN(vect_temp) = (2*sqrt(trans_signal_power)*mean_h2/noise_power_input)*x_rec_AWGN(vect_temp);
x_rec_AWGN(vect_temp) = (2*sqrt(trans_signal_power)*mean_h2/noise_power_input)*x_rec_AWGN(vect_temp);

%%Rayleigh + TSMG channel Using Gaussian Receiver Decoding
x_rec_GR(vect_temp) = (2*sqrt(trans_signal_power)*mean_h2/noise_power_input)*x_rec(vect_temp);

%%Rayleigh + TSMG channel Using blanking
x_rec_blanked_fixed(vect_temp) = (2*sqrt(trans_signal_power)*mean_h2/noise_power_input)*x_rec_blanked_fixed(vect_temp);

%%Rayleigh + TSMG channel Using Clipping
x_rec_clipped_fixed(vect_temp) = (2*sqrt(trans_signal_power)*mean_h2/noise_power_input)*x_rec_clipped_fixed(vect_temp);

%%Rayleigh + TSMG channel Using Blanking and Clipping
x_rec_clip_blank_fixed(vect_temp) = (2*sqrt(trans_signal_power)*mean_h2/noise_power_input)*x_rec_clip_blank_fixed(vect_temp);


end

metric_ratio_MAP_SD = reshape(metric_ratio_MAP_SD,[len,1]);
dec_bits_di_MAP=dec(metric_ratio_MAP_SD);
error_Direct_MAP(snr)=error_Direct_MAP(snr)+length(find(dec_bits_di_MAP~=source_bits));



% 
% 
% 
% performance of TSMG channel with AWGN receiver
x_rec_GR = reshape(x_rec_GR,[len,1]);
dec_bits_Gaussian=dec(x_rec_GR);
error_Direct_Gaussian(snr)=error_Direct_Gaussian(snr)+length(find(dec_bits_Gaussian~=source_bits));

% performance over AWGN channel
dec_bits_AWGN=dec(x_rec_AWGN);
error_Direct_AWGN(snr)=error_Direct_AWGN(snr)+length(find(dec_bits_AWGN~=source_bits));

% performance using blanking
dec_bits_di_blank_fixed=dec(x_rec_blanked_fixed);
error_Direct_blank_fixed(snr)=error_Direct_blank_fixed(snr)+length(find(dec_bits_di_blank_fixed~=source_bits));

% performance using clipping
dec_bits_di_clipped_fixed=dec(x_rec_clipped_fixed);
error_Direct_clipped_fixed(snr)=error_Direct_clipped_fixed(snr)+length(find(dec_bits_di_clipped_fixed~=source_bits));

% performance using blanking and clipping
dec_bits_di_clip_blank_fixed=dec(x_rec_clip_blank_fixed);
error_Direct_clip_blank_fixed(snr)=error_Direct_clip_blank_fixed(snr)+length(find(dec_bits_di_clip_blank_fixed~=source_bits));
















end
end


%
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
 
ber_Direct_MAP=error_Direct_MAP./(N_sample/2);

ber_Direct_blanked_fixed=error_Direct_blank_fixed./(N_sample/2);

ber_Direct_clipped_fixed=error_Direct_clipped_fixed./(N_sample/2);

ber_Direct_clip_blank_fixed=error_Direct_clip_blank_fixed./(N_sample/2);

ber_Direct_Gaussian=error_Direct_Gaussian./(N_sample/2);

ber_Direct_AWGN=error_Direct_AWGN./(N_sample/2);

semilogy(snr_db,ber_Direct_MAP,'-*r',snr_db,ber_Direct_blanked_fixed,'-db',snr_db,ber_Direct_clipped_fixed,'-dr',snr_db,ber_Direct_clip_blank_fixed,'-dg',snr_db,ber_Direct_Gaussian,'-b',snr_db,ber_Direct_AWGN,'-r');

legend('MAP','blanked fixed T','clipped fixed T','clip blank fixed T','AWGN Receiver','AWGN Noise');


ylim([10^(-5) 1]);
l = size(snr_db);
l = l(2);
snr_db = reshape(snr_db,[l,1]);
BCJR = reshape(ber_Direct_MAP,[l,1]);
AWGN_Channel = reshape(ber_Direct_AWGN,[l,1]);
AWGN_Receiver = reshape(ber_Direct_Gaussian,[l,1]);
blanking = reshape(ber_Direct_blanked_fixed,[l,1]);
clipping = reshape(ber_Direct_clipped_fixed,[l,1]);
cliping_blanking = reshape(ber_Direct_clip_blank_fixed,[l,1]);

T = table(snr_db,AWGN_Channel,AWGN_Receiver,BCJR,blanking,clipping,cliping_blanking);
writetable(T,'conventional methods R = 50, G = 100.csv');



%save SimulationData.mat;
