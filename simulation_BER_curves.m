%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------

close all
clear all

%% INPUT DATA

M = 64;                                         % Size of signal constellation
k = log2(M);                                    % Number of bits per symbol    
N = 1024;                                       % Number of total carriers
usedN = 600;                                    % Number of data carriers
unusedN = N-usedN;                              % Number of guard carriers

nSymbOFDM = 1000;                               % Number of OFDM Symbols input
n = usedN*k*nSymbOFDM;                          % Number of bits

CP = N/8;                                       % Cyclic Prefix length samples
ZT = N/8;                                       % Zero Tail length samples

nSymbEst = 2;                                   % Number of channel estimation OFDM symbols

t = 0:1:120;                                    % Time vector to represent multi-path channel
BW = 20e6;  % Hz                                % Bandwidth of the system
Ts = (1/BW)*1e9; %ns                            % Sampling period

% EbN0 vector to generate the BER curve
EbN0_dB = 0:1:20;                               % Bit energy to noise ratio of the simulation

% Intialization
HdB = -inf.*ones(1,length(t));                  % Multi-path channel
SER = ones(1,length(EbN0_dB)).*inf;             % Symbol Error Rate
BER = ones(1,length(EbN0_dB)).*inf;             % Bit Error Rate OFDM
BER_ZT = ones(1,length(EbN0_dB)).*inf;          % Bit Error Rate OFDM-Zero Tail 
BER_QAM = ones(1,length(EbN0_dB)).*inf;         % Bit Error Rate QAM

% Different channel models
%   Uncoment the choosen one and coment the rest to avoid overwrite
%--------------------------------------------------------------------------

% Single-path channel
HdB(1) = 0;

% % [EPA] Extended Pedestrian A model 
% HdB(1) = 0;                                                              
% HdB(ceil(51/Ts)) = -1;
% HdB(ceil(71/Ts)) = -2;
% HdB(ceil(91/Ts)) = -3;
% HdB(ceil(111/Ts)) = -8;
% HdB(ceil(191/Ts)) = -17.2;
% HdB(ceil(411/Ts)) = -20.8;

% % [EVA] Extended Vehicular A model
% HdB(1) = 0;
% HdB(ceil(51/Ts)) = -1.5;
% HdB(ceil(151/Ts)) = -1.4;
% HdB(ceil(311/Ts)) = -3.6;
% HdB(ceil(371/Ts)) = -0.6;
% HdB(ceil(711/Ts)) = -9.1;
% HdB(ceil(1091/Ts)) = -7;
% HdB(ceil(1731/Ts)) = -12;
% HdB(ceil(2511/Ts)) = -16.9;

% % [ETU] Extended Typical Urban model
% HdB(1) = -1;
% HdB(ceil(51/Ts)) = -1;
% HdB(ceil(121/Ts)) = -1;
% HdB(ceil(201/Ts)) = 0;
% HdB(ceil(231/Ts)) = 0;
% HdB(ceil(501/Ts)) = 0;
% HdB(ceil(1601/Ts)) = -3;
% HdB(ceil(2301/Ts)) = -5;
% HdB(ceil(5001/Ts)) = -7;

%--------------------------------------------------------------------------

H = 10.^(HdB/10);                               % Convert channel taps to natural values

figure
stem(t,HdB)
title('CHANNEL')
xlabel('Time Samples')
ylabel('Signal realtion (dB)')

% Input bits
dataIn = randi([0 1],n,1);                      % Generate vector of random binary data

%--------------------------------------------------------------------------
%   All the simulations have the same input data to compare between them. 
%   Uncoment the whole section (%% OFDM______) to simulate
%--------------------------------------------------------------------------

%% OFDM____________________________________________________________________ 

% TX ----------------------------------------------------------------------
[ofdm , dataMod] = TX_OFDM(dataIn,M,N,usedN,CP);

% NOISE -------------------------------------------------------------------


for z=1:length(EbN0_dB)
    
    [channelCorrection] = CHANNEL_ESTIMATION...
        (H,nSymbEst,EbN0_dB(z),k,N,usedN,CP);

    [ofdmChannel] = CHANNEL_OFDM(ofdm,H);
    
    [ofdmAWGN] = AWGN_OFDM(EbN0_dB(z),ofdmChannel,k,N,usedN,CP);
    
    [dataInRx , dataModRx] = RX_CHANNEL_OFDM...
        (ofdmAWGN,M,N,usedN,CP,channelCorrection);
    
    dataSymbolsIn = qamdemod(dataMod,M,0,'gray');
    dataSymbolsInRx = qamdemod(dataModRx,M,0,'gray');
    
    [SER(z)] = sum(dataSymbolsIn ~= dataSymbolsInRx)./...
        length(dataSymbolsIn);
    [~, BER(z)] = biterr(dataIn,dataInRx);
    
end

% %% OFDM ZT_________________________________________________________________
% 
% % TX ----------------------------------------------------------------------
% [ofdmZT , dataModZT] = TX_OFDM_ZEROTAIL(dataIn,M,N,usedN,ZT);
% 
% % NOISE -------------------------------------------------------------------
% 
% for z=1:length(EbN0_dB)
%     
%     [channelCorrectionZT] = CHANNEL_ESTIMATION_ZEROTAIL...
%         (H,nSymbEst,EbN0_dB(z),k,N,usedN,ZT);
% 
%     [ofdmChannelZT] = CHANNEL_OFDM(ofdmZT,H);
%     
%     [ofdmAWGNZT] = AWGN_OFDM(EbN0_dB(z),ofdmChannelZT,k,N,usedN,ZT);
%     
%     [dataInRxZT , dataModRxZT] = RX_CHANNEL_OFDM_ZEROTAIL...
%         (ofdmAWGNZT,M,N,usedN,ZT,channelCorrectionZT);
%     
%     dataSymbolsZT = qamdemod(dataModZT,M,0,'gray');
%     dataSymbolsRxZT = qamdemod(dataModRxZT,M,0,'gray');
%     
%     [SER_ZT(z)] = sum(dataSymbolsZT ~= dataSymbolsRxZT)./...
%         length(dataSymbolsZT);
%     [~, BER_ZT(z)] = biterr(dataIn,dataInRxZT);
%     
% end

%% QAM_____________________________________________________________________

% MODUALTION---------------------------------------------------------------
% Reshape data into binary 4-tuples
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);                          
% Convert to integers
dataSymbolsIn = bi2de(dataInMatrix);                                        
% Modulate the data with Gray code
dataMod = qammod(dataSymbolsIn,M,0,'gray');                                 

% NOISE -------------------------------------------------------------------

for z=1:length(EbN0_dB)

    % Signal to noise ratio
    snrdB = EbN0_dB(z) + 10*log10(k);                                       
    % AWGN channel
    dataModNoise = awgn(dataMod,snrdB,'measured');                          
    
    % Demodulate with Gray code
    dataSymbolsInRx = qamdemod(dataModNoise,M,0,'gray');                    
    % Convert to binary data
    dataInMatrixRx = de2bi(dataSymbolsInRx);                                
    % Reshape into a binary vector
    dataInRx = reshape (dataInMatrixRx,size(dataInMatrixRx,1)...            
        *size(dataInMatrixRx,2),1);
    
    [SER(z)] = sum(dataSymbolsIn ~= dataSymbolsInRx)./...
        length(dataSymbolsIn);
    [~, BER_QAM(z)] = biterr(dataIn,dataInRx);
    
end

%% CURVES__________________________________________________________________

% Theoretical BER curve
EbN0 = 10.^(EbN0_dB/10);
SER_MQAM = 2*erfc(sqrt((3*k*EbN0)/(2*(M-1))));
BER_MQAM = SER_MQAM./k;                         % Bit Error Rate Theoretical

figure
semilogy(EbN0_dB,[BER;BER_ZT;BER_QAM;BER_MQAM],'LineWidth',2)
grid on;
legend('Simulation OFDM','Simulation OFDM ZT','Simulation QAM','Theoretical');
xlabel('EbN0 (dB)');ylabel('Bit Error Rate');
title('BER OFDM')
