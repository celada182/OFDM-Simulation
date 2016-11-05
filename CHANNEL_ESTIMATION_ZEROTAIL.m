%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   CHANNEL ESTIMATION FOR OFDM ZEROTAIL 
%   Estimation of the channel used in the simulation to obtain the
%   correction values for the OFDM carriers to avoid the efect of the
%   multipath channel.
%   
%   [ channelCorrection ] = CHANNEL_ESTIMATION_ZEROTAIL...
%                   ( H , nSymbEst, EbN0_dB, k , N , usedN, ZT )
%
%   Funcitions returns channelCorrection (Nx1) and needs this parameters input:
%
%   H->         Channel
%   nSymbEst->  Number of OFDM test symbols to do the estimation
%   EbN0_dB->   Bit energy to noise ratio
%   k->         Bits per QAM symbol
%   N->         Number of OFDM carriers
%   usedN->     Numer of data carriers
%   ZT->        Length of the zero tail

function [ channelCorrection ] = CHANNEL_ESTIMATION_ZEROTAIL...
    ( H , nSymbEst, EbN0_dB, k , N ,usedN, ZT )

unusedN = N - usedN;

%TX
dataModUsedN = ones(usedN,nSymbEst);            % Data test to send 
dataModN = zeros(N,nSymbEst);                   % through the channel

ofdm = zeros(size(dataModN,2)*N,1);
ofdmSymbol = zeros(N,size(dataModN,2));
ofdmSymbolZT = zeros(N+ZT,size(dataModN,2));
zerotail = zeros (ZT,size(dataModN,2));

for j=1:size(dataModN,2)
    dataModN(:,j) = vertcat(zeros(unusedN/2,1),dataModUsedN(:,j),zeros(unusedN/2,1));
    % Create Symbol OFDM
    ofdmSymbol(:,j) = ifft(dataModN(:,j),N);                                
    % Add the ending zeros
    ofdmSymbolZT(:,j) = vertcat(ofdmSymbol(:,j),zerotail(:,j));             
    % Create the vector with all the OFDM symbols  
    ofdm((j-1)*length(ofdmSymbolZT)+1:j*length(ofdmSymbolZT)) = ofdmSymbolZT(:,j);
end

%CHANNEL
ofdmChannel = filter (H,1,ofdm);                % Channel Multipath

%AWGN
EbN0 = 10^(EbN0_dB/10);
snr = (N/(N+ZT))*(usedN/N)*EbN0*k;
snrdB = 10*log10(snr);                          % SNR from EbN0
ofdmAWGN = awgn(ofdmChannel,snrdB,'measured');  % Channel AWGN

%RX
ofdmSymbolZTRx = reshape (ofdmAWGN,N+ZT,length(ofdmAWGN)/(N+ZT));
dataModNRx = zeros(N,size(ofdmSymbolZTRx,2));
channelCorrectionMatrix = zeros(N,size(ofdmSymbolZTRx,2));

for j=1:size(ofdmSymbolZTRx,2)
% Info of the zero carriers is added to the frist ZT ones                                                                     
    ofdmSymbolZTRx(1:ZT,j) = ofdmSymbolZTRx(1:ZT,j) + ofdmSymbolZTRx(N+1:end,j);
    ofdmSymbolRx = ofdmSymbolZTRx(1:N,:);
    % dataModNRx = H
    dataModNRx(:,j) = fft(ofdmSymbolRx(:,j),N);                             
    % X / Y = 1 / H
    channelCorrectionMatrix(:,j)=dataModN(:,j)./dataModNRx(:,j);            
end

channelCorrection = mean(channelCorrectionMatrix,2);                        

% % Channel estimation figure
% figure
% subplot(4,1,1)
% stem(real(dataModN(:,1)))
% title('Test Data Sent')
% subplot(4,1,2)
% stem(real(dataModNRx(:,1)))
% title('Test Data Received')
% subplot(4,1,3)
% stem(real(channelCorrection))
% title('Channel Correction')
% subplot(4,1,4)
% stem(real(channelCorrection.*dataModNRx(:,1)))
% title('Test Data Received X Channel Correction')

end