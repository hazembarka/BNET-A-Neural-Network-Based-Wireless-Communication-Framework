function [corrNoise] = TwoStateMarkovGaussianSamples_BPSK_realnoise(p_B,gamma,s_0,s_1,N_sample)

%N=10e6;     %number of samples to genrate
p_G=1-p_B;

p_BG=(1/gamma)*p_G; p_GB=(1/gamma)*p_B;

p_BB=1-p_BG;
 
corrNoise=sqrt(0.5)*(randn(1,N_sample)+1i*randn(1,N_sample));       % Sample vector

usample=rand(1,N_sample);          % Dice vector
 
%p01=Np/(N-Np);      % Probablilty to create an impulse
%p11=1-(1/td);       % Probability to stay in an impulse
c0=0;               % Transition to Background counter
c1=0;               % Transition to impulse counter
state=0;            %initial state
 
count0=0;           % Totla number of samples in state 0
count1=0;           % Totla number of samples in state 1
 
for q=1:N_sample
    
    if(state==0 && usample(q)<p_GB)
        state=1;
        c1=c1+1;
    elseif (state==1 && usample(q)>p_BB)
        state=0; 
        c0=c0+1;
    end
    
    if state==0
        corrNoise(q)=s_0.*corrNoise(q);
        count0=count0+1;
    elseif state==1
        corrNoise(q)=s_1.*corrNoise(q);
        count1=count1+1;
    end
end

 