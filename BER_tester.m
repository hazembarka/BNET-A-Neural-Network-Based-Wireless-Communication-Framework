% Communication chain
% Communication chain

clear all; close all;

frameNb_test = 256;


tic;

% Parameters


R     = 50;
GAMMA = 100;
P_B = 0.1;
backbone_string = ' + LSTM';


SEQUENCE_LENGTH = 32;






load('all_MAP','enc','dec');


snr_db  = [4:0.5:6];
snr_lin=10.^(snr_db./10); snr_size=length(snr_lin);


BER_BCJRs = [];
BER_BNETs = [];
BER_Nothing_s = [];
snr_s = [];


for snr=1:snr_size

    chSNR=snr_lin(snr);
    noise_power_input=1/chSNR;
    snr_db_value = snr_db(snr)
    
    
    source_bits = readtable(['test datasets/G',num2str(GAMMA),'_R',num2str(R),'/',num2str(snr_db_value),' original bits.csv']);
    sub_LLR_table  = readtable(['test datasets/G',num2str(GAMMA),'_R',num2str(R),'/',num2str(snr_db_value),' predictions BNET',backbone_string,'.csv']);

   
    y_SD = sub_LLR_table.('Received');
    y_SD = ((2/noise_power_input)* y_SD);
    
    MAP_BNET = sub_LLR_table.('MAP_BNET');
    MAP_BCJR_decoded =source_bits.('all_decoded_original');
    SOURCE_BITS = source_bits.('all_bits');

    %%Decoding
    temp_BER_BCJR =[];
    temp_BER_Nothing =[];
    temp_BER_BNET  =[];
    
    for l=0:frameNb_test-1
        Chunk_y_SD = y_SD(l*64800+1:(l+1)*64800);
        Chunk_SOURCE_BITS = SOURCE_BITS(l*32400+1:(l+1)*32400);
        Chunk_MAP_BNET = MAP_BNET(l*64800+1:(l+1)*64800);
        
        decoded_Nothing = dec(Chunk_y_SD);
        decoded_BCJR =  MAP_BCJR_decoded(l*32400+1:(l+1)*32400);
        decoded_BNET = dec(Chunk_MAP_BNET);

        BER_BCJR = length(find(decoded_BCJR~=Chunk_SOURCE_BITS))/(32400);
        BER_BNET  = length(find(decoded_BNET~=Chunk_SOURCE_BITS))/(32400);
        BER_Nothing  = length(find(decoded_Nothing~=Chunk_SOURCE_BITS))/(32400);
        
        temp_BER_BCJR =[temp_BER_BCJR;BER_BCJR];
        temp_BER_Nothing =[temp_BER_Nothing;BER_Nothing];
        temp_BER_BNET  =[temp_BER_BNET;BER_BNET];

    end
    BER_BCJR = mean(temp_BER_BCJR)
    BER_BNET  = mean(temp_BER_BNET)
    BER_Nothing  = mean(temp_BER_Nothing)
        
    
    BER_Nothing_s = [BER_Nothing_s;BER_Nothing];
    BER_BCJRs = [BER_BCJRs;BER_BCJR];
    BER_BNETs = [BER_BNETs;BER_BNET];
    
    snr_s = [snr_s;snr_db_value];

end

T = table(snr_s,BER_BCJRs,BER_BNETs,BER_Nothing_s);
writetable(T,['results 3/BNET',backbone_string,', gamma = ',num2str(GAMMA),', R = ',num2str(R),', p_B = ',num2str(P_B),'.csv']);


toc;


