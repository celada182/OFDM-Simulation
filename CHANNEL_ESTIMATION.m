%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   CHANNEL ESTIMATION FOR OFDM
%   Estimation of the channel used in the simulation to obtain the
%   correction values for the OFDM carriers to avoid the efect of the
%   multipath channel.
%   
%   [ channelCorrection ] = CHANNEL_ESTIMATION( H , nSymbEst, EbN0_dB, k , N , usedN, CP )
%
%   Funcitions returns channelCorrection (Nx1) and needs this parameters input:
%
%   H->         Channel
%   nSymbEst->  Number of OFDM test symbols to do the estimation
%   EbN0_dB->   Bit energy to noise ratio in dB
%   k->         Bits per QAM symbol
%   N->         Number of OFDM carriers
%   usedN->     Number of data carriers
%   CP->        Length of cyclic prefix

function [ channelCorrection ] = CHANNEL_ESTIMATION...
    ( H , nSymbEst, EbN0_dB, k , N ,usedN, CP)

unusedN = N - usedN;

%TX
dataModUsedN = ones(usedN,nSymbEst);                % Data test to send                                              
dataModN = zeros(N,nSymbEst);

ofdm = zeros(size(dataModN,2)*N,1);
ofdmSymbol = zeros(N,size(dataModN,2));
ofdmSymbolCP = zeros(N+CP,size(dataModN,2));

for j=1:size(dataModN,2)
    dataModN(:,j) = vertcat(zeros(unusedN/2,1),dataModUsedN(:,j),zeros(unusedN/2,1));
    % Create Symbol OFDM
    ofdmSymbol(:,j) = ifft(dataModN(:,j),N);
    % Add the cyclic prefix at the beginning of the OFDM Symbol
    ofdmSymbolCP(:,j) = vertcat(ofdmSymbol(N-CP+1:N,j),ofdmSymbol(:,j));    
    % Create the vector with all the OFDM symbols                                                               
    ofdm((j-1)*length(ofdmSymbolCP)+1:j*length(ofdmSymbolCP)) = ofdmSymbolCP(:,j);
end

%CHANNEL
ofdmChannel = filter (H,1,ofdm);                    % Channel Multipath  

%AWGN
EbN0 = 10^(EbN0_dB/10);
snr = (N/(N+CP))*(usedN/N)*EbN0*k;
snrdB = 10*log10(snr);                              % SNR from EbN0
ofdmAWGN = awgn(ofdmChannel,snrdB,'measured');      % Channel AWGN

ofdmSymbolCPRx = reshape (ofdmAWGN,N+CP,length(ofdmAWGN)/(N+CP));
ofdmSymbolRx = ofdmSymbolCPRx(CP+1:end,:);

dataModNRx = zeros(N,size(ofdmSymbolRx,2));
channelCorrectionMatrix = zeros(N,size(ofdmSymbolRx,2));

for j=1:size(ofdmSymbolRx,2)
    % Info of the CP carriers is ignored
    % dataModNRx = H
    dataModNRx(:,j) = fft(ofdmSymbolRx(:,j),N);                             
    % X / Y = 1 / H
    channelCorrectionMatrix(:,j)=dataModN(:,j)./dataModNRx(:,j);            
end

channelCorrection = mean(channelCorrectionMatrix,2);    % Estimation nSymb

% DEMO 2
% Channel estimation figure
figure
subplot(4,1,1)
stem(real(dataModN(:,1)))
title('Test Data Sent')
subplot(4,1,2)
stem(real(dataModNRx(:,1)))
title('Test Data Received')
subplot(4,1,3)
stem(real(channelCorrection))
title('Channel Correction')
subplot(4,1,4)
stem(real(channelCorrection.*dataModNRx(:,1)))
title('Test Data Received X Channel Correction')

keyboard
end

