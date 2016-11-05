%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------

clear all
close all

%% INPUT DATA

M = 64;                                             % Size of signal constellation
k = log2(M);                                        % Number of bits per symbol
N = 2048;                                           % Number of carriers
usedN = 1200;                                       % Number of data carriers
unusedN = N-usedN;                                  % number of guard carriers
nSymbOFDM = 100;                                    % Number of OFDM Symbols input

n = usedN*k*nSymbOFDM;                              % Input number of bits

CP = N/8;                                           % Cyclic Prefix length
ZT = N/8;                                           % Zero tail length

dataIn = randi([0 1],n,1);                          % Generate vector of binary data

%% OFDM ZERO TAIL DFT-spread ----------------------------------------------

% TX ----------------------------------------------------------------------
[ofdmZTDFT , ~] = TX_OFDM_ZEROTAIL_DFT(dataIn,M,N,usedN,ZT);

%% OFDM ZEROTAIL ----------------------------------------------------------

% TX ----------------------------------------------------------------------
[ofdmZT , ~] = TX_OFDM_ZEROTAIL(dataIn,M,N,usedN,ZT);

%% OFDM -------------------------------------------------------------------

% TX ----------------------------------------------------------------------
[ofdm , ~] = TX_OFDM(dataIn,M,N,usedN,CP);

%% Waveforms --------------------------------------------------------------
figure
subplot(3,1,1)
plot(real(ofdm(1:N+300)))
title('OFDM')
subplot(3,1,2)
plot(real(ofdmZT(1:N+300)),'r')
title('OFDM ZT')
subplot(3,1,3)
plot(real(ofdmZTDFT(1:N+300)),'g')
title('OFDM ZT DFTs')

%% Spectrum ---------------------------------------------------------------
figure
[pxx,f] = periodogram(ofdmZT);
plot(f/pi,10*log10(pxx),'r')
hold on
[pxx,f] = periodogram(ofdm);
plot(f/pi,10*log10(pxx))
hold on
[pxx,f] = periodogram(ofdmZTDFT);
plot(f/pi,10*log10(pxx),'g')
ylabel('Power/frecuency (dB/rad/sample)')
xlabel('Normalized Frecuency (xpi rad/sample)')
title('OFDM ZEROTAIL vs DFTs')
legend('OFDM ZT', 'OFDM','OFDM ZT DFTs')
