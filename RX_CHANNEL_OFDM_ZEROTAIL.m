%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   RECEPTION OF THE OFDM ZEROTAIL SIGNAL WITH MULTIPATH CHANNEL CORRECTION
%   Simulation of the extraction of the OFDM symbols, their useful
%   information and the final chain of QAM symbols to demodulate.
%
%   [ dataInRx,dataModRxFixed ] = ...
%       RX_CHANNEL_OFDM_ZEROTAIL( ofdm,M,N,usedN,ZT,channelCorrection )
%
%   MUST HAVE channelCorrection vector simulated in the function:
%   CHANNEL_ESTIMATION_ZEROTAIL.m
%
%   OUTPUT:
%   dataInRx -> Bit vector obtain from the ofdm signal
%   dataModRxFixed -> QAM symbols inside the carriers OFDM with the
%                       channel correction necessary.
%
%   INPUT:
%   ofdm -> Chain of OFDM symbols
%   M -> Number of symbols of the QAM constellation
%   N -> Number of OFDM carriers
%   usedN -> Number of data OFDM carriers
%   ZT -> Length of the zero tail
%   channelCorrection -> 1xN vector equivalent to 1/H, the compensation of
%                           the multipath channel effect.

function [ dataInRx,dataModRxFixed ] = RX_CHANNEL_OFDM_ZEROTAIL( ofdm,M,N,usedN,ZT,channelCorrection )

unusedN = N - usedN;

ofdmSymbolZTRx = reshape (ofdm,N+ZT,length(ofdm)/(N+ZT));

dataModNRx = zeros(N,size(ofdmSymbolZTRx,2));
dataModUsedNRx = zeros(usedN,size(ofdmSymbolZTRx,2));
dataModNRxFixed = zeros(N,size(ofdmSymbolZTRx,2));

for j=1:size(ofdmSymbolZTRx,2)
    % Info of the zero carriers is added to the frist ZT ones
    ofdmSymbolZTRx(1:ZT,j) = ofdmSymbolZTRx(1:ZT,j) + ofdmSymbolZTRx(N+1:end,j);
    ofdmSymbolRx = ofdmSymbolZTRx(1:N,j);
    % Carries with QAM symbols
    dataModNRx(:,j) = fft(ofdmSymbolRx,N);
    % Carriers channel fixed
    dataModNRxFixed (:,j) = dataModNRx(:,j).*channelCorrection;
    % Recover data carriers
    dataModUsedNRx(:,j) = dataModNRxFixed(unusedN/2+1:N-unusedN/2,j);
end

% Info of the OFDM carriers to a chain of QAM simbols to get data vector 
dataModRxFixed = reshape...
    (dataModUsedNRx,size(dataModUsedNRx,1)*size(dataModUsedNRx,2),1);
dataSymbolsInRx = qamdemod(dataModRxFixed,M,0,'gray');
dataInMatrixRx = de2bi(dataSymbolsInRx);
dataInRx = reshape...
    (dataInMatrixRx,size(dataInMatrixRx,1)*size(dataInMatrixRx,2),1);

end
