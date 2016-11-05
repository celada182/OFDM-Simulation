%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------

clc
clear all
close all

%% INPUT DATA

M = 16;                                         % Size of signal constellation
k = log2(M);                                    % Number of bits per symbol
N = 200;                                       % Number of total carriers
usedN = 150;                                    % Number of data carriers
unusedN = N-usedN;                              % Number of guard carriers

nSymbOFDM = 2;                                % Number of OFDM Symbols input
n = usedN*k*nSymbOFDM;                          % Number of bits

CP = N/8;                                       % Cyclic Prefix length samples
ZT = N/8;                                       % Zero Tail length samples

nSymbEst = 2;                                   % Number of channel estimation OFDM symbols

EbN0_dB = inf;                                   % Bit energy to noise ratio of the simulation
disp('EbN0 dB FOR ALL TRANSMISSIONS:');
disp(EbN0_dB);

t = 0:1:120;                                    % Time vector to represent multi-path channel
BW = 20e6;  % Hz                                % Bandwidth of the system
Ts = (1/BW)*1e9; %ns                            % Sampling period

% Intialization
HdB = -inf.*ones(1,length(t));                  % Multi-path channel

% Different channel models
%   Uncoment the choosen one and coment the rest to avoid overwrite
%--------------------------------------------------------------------------

% % Single-path Channel
% HdB(1) = 0;

% [EPA] Extended Pedestrian A model 
HdB(1) = 0;                                                              
HdB(ceil(51/Ts)) = -1;
HdB(ceil(71/Ts)) = -2;
HdB(ceil(91/Ts)) = -3;
HdB(ceil(111/Ts)) = -8;
HdB(ceil(191/Ts)) = -17.2;
HdB(ceil(411/Ts)) = -20.8;

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

% Theoretical Bit Erro Rate
EbN0 = 10.^(EbN0_dB/10);
SER_MQAM = 2*erfc(sqrt((3*k*EbN0)/(2*(M-1))));
BER_MQAM = SER_MQAM./k;

disp('THEORETICAL BER');
disp(BER_MQAM);

% Input bits
dataIn = randi([0 1],n,1);                      % Generate vector of random binary data

%--------------------------------------------------------------------------
%   All the simulations have the same input data to compare between them. 
%   Uncoment the whole section (%% OFDM______) to simulate
%--------------------------------------------------------------------------

%% OFDM CHANNEL FIXED SIMULATION___________________________________________

% CHANNEL ESTIMATION ------------------------------------------------------
[channelCorrection] = CHANNEL_ESTIMATION(H,nSymbEst,EbN0_dB,k,N,usedN,CP);  

% TX ----------------------------------------------------------------------
[ofdm , dataMod] = TX_OFDM(dataIn,M,N,usedN,CP);

% CHANNEL -----------------------------------------------------------------
[ofdmChannel] = CHANNEL_OFDM(ofdm, H);

% NOISE -------------------------------------------------------------------
[ofdmAWGN] = AWGN_OFDM(EbN0_dB,ofdmChannel,k,N,usedN,CP);

% RX ----------------------------------------------------------------------
[dataInRx , dataModRxFixed] = RX_CHANNEL_OFDM...
    (ofdmAWGN,M,N,usedN,CP,channelCorrection);

% BER----------------------------------------------------------------------
[~, BER] = biterr(dataIn,dataInRx);             % Calculate BER comparing output/input
disp('BER OFDM CHANNEL FIXED');
disp(BER);

% Constellation -----------------------------------------------------------
sPlotFig = scatterplot(dataModRxFixed,1,0,'g.');
grid on
hold on
scatterplot(dataMod,1,0,'k*',sPlotFig)
title('OFDM CHANNEL FIXED CONSTELLATION')

keyboard

%% OFDM ZEROTAIL CHANNEL FIXED SIMULATION__________________________________

% CHANNEL ESTIMATION ------------------------------------------------------
[channelCorrection] = CHANNEL_ESTIMATION_ZEROTAIL...
    (H,nSymbEst,EbN0_dB,k,N,usedN,ZT);

% TX ----------------------------------------------------------------------
[ofdmZT , dataMod] = TX_OFDM_ZEROTAIL(dataIn,M,N,usedN,ZT);

% CHANNEL -----------------------------------------------------------------
[ofdmChannel] = CHANNEL_OFDM(ofdmZT, H);

% NOISE -------------------------------------------------------------------
[ofdmAWGN] = AWGN_OFDM(EbN0_dB,ofdmChannel,k,N,usedN,ZT);

% RX ----------------------------------------------------------------------
[dataInRx , dataModRxFixed] = RX_CHANNEL_OFDM_ZEROTAIL...
    (ofdmAWGN,M,N,usedN,ZT,channelCorrection);

% BER----------------------------------------------------------------------
[~, BER] = biterr(dataIn,dataInRx);             % Calculate BER comparing output/input
disp('BER OFDM ZEROTAIL CHANNEL FIXED');
disp(BER);

% Constellation -----------------------------------------------------------
sPlotFig = scatterplot(dataModRxFixed,1,0,'g.');
hold on
grid on
scatterplot(dataMod,1,0,'k*',sPlotFig)
title('OFDM ZEROTAIL CHANNEL FIXED CONSTELLATION')

keyboard

%% QAM_____________________________________________________________________

% MODUALTION---------------------------------------------------------------
% Reshape data into binary 4-tuples
dataInMatrix = reshape(dataIn,length(dataIn)/k,k);                          
% Convert to integers
dataSymbolsIn = bi2de(dataInMatrix);                                        
% Modulate the data with Gray code
dataMod = qammod(dataSymbolsIn,M,0,'gray');                                 

% NOISE -------------------------------------------------------------------
% Signal to noise ratio
snrdB = EbN0_dB + 10*log10(k);                                              
% AWGN channel
dataModNoise = awgn(dataMod,snrdB,'measured');                              

% Demodulate with Gray code
dataSymbolsInRx = qamdemod(dataModNoise,M,0,'gray');                        
% Convert to binary data
dataInMatrixRx = de2bi(dataSymbolsInRx);                                    
% Reshape into a binary vector
dataInRx = reshape (dataInMatrixRx,size(dataInMatrixRx,1)...                
    *size(dataInMatrixRx,2),1);

% BER----------------------------------------------------------------------
[~, BER] = biterr(dataIn,dataInRx);                 % Calculate BER comparing output/input
disp('BER QAM');
disp(BER);

% Constellation -----------------------------------------------------------
sPlotFig = scatterplot(dataModNoise,1,0,'g.');
hold on
grid on
scatterplot(dataMod,1,0,'k*',sPlotFig)
title('QAM CONSTELLATION')

