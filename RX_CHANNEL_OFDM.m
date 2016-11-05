%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   RECEPTION OF THE OFDM SIGNAL WITH CHANNEL CORRECTION
%   Simulation of the extraction of the OFDM symbols, their useful
%   information and the final chain of QAM symbols to demodulate.
%
%   [ dataInRx,dataModRxFixed ] = RX_CHANNEL_OFDM( ofdm,M,N,usedN,CP,channelCorrection )
%
%   MUST HAVE channelCorrection vector simulated in the function:
%   CHANNEL_ESTIMATION.m
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
%   usedN -> Number of OFDM data carriers
%   CP -> Length of the cyclic preix
%   channelCorrection -> 1xN vector equivalent to 1/H, the compensation of
%                                               the multipath channel effect.

function [ dataInRx,dataModRxFixed ] = RX_CHANNEL_OFDM( ofdm,M,N,usedN,CP,channelCorrection )

unusedN = N - usedN;

ofdmSymbolCPRx = reshape (ofdm,N+CP,length(ofdm)/(N+CP));
ofdmSymbolRx = ofdmSymbolCPRx(CP+1:end,:);

dataModNRx = zeros(N,size(ofdmSymbolRx,2));
dataModUsedNRx = zeros(usedN,size(ofdmSymbolRx,2));
dataModNRxFixed = zeros(N,size(ofdmSymbolRx,2));

for j=1:size(ofdmSymbolRx,2)
    % Info of the CP carriers is ignored
    % Carries with QAM symbols
    dataModNRx(:,j) = fft(ofdmSymbolRx(:,j),N);                             
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

