%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   TRANSMISSION OF THE OFDM DFT-SPREAD ZEROTAIL SIGNAL
%   Simulation of the creation of the OFDM signal from data bits input.
%   Distribution of the bits to create the QAM simbols and place them in the
%   OFDM carriers to form independent OFDM symbols with the zero tail.
%   
%   [ ofdm,dataMod ] = TX_OFDM_ZEROTAIL_DFT( dataIn,M,N,usedN,ZT )
%
%   OUTPUT:
%   ofdm -> Chain of OFDM symbols
%   dataMod -> QAM simbols inside the OFDM carriers
%
%   INPUT:
%   dataIn -> Vector of bits
%   M -> Number of simbols of the QAM constellation
%   N -> Number of carriers OFDM
%   usedN -> Number of data carriers OFDM
%   ZT -> Length of the zero tail

function [ ofdm,dataMod ] = TX_OFDM_ZEROTAIL_DFT( dataIn,M,N,usedN,ZT )

k = log2(M);
unusedN = N - usedN;
ZTh = 10;                                               % Zero header
ZTt = ZT - ZTh;                                         % Zero tail    

dataInMatrix = reshape(dataIn,length(dataIn)/k,k);      % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                    % Convert to integers

dataMod = qammod(dataSymbolsIn,M,0,'gray');

% Used carriers now are usedN - ZT because of the ifft size
% The dataModUsedN matrix discards the last QAM symbols
reshapeUsedN = floor(length(dataMod)/(usedN-ZT));       % Number of OFDM symbols to send
disp('Discarded last QAM symbols in ZT-DFTs');
disp(length(dataMod) - reshapeUsedN*(usedN-ZT));
dataMod = dataMod(1:reshapeUsedN*(usedN-ZT));

dataModUsedN = reshape(dataMod,usedN-ZT,reshapeUsedN);  % Data modulated for the used carriers
dataModN = zeros (N,size(dataModUsedN,2));              % Data added guard carriers

ofdm = zeros(size(dataModN,2)*(N+ZT),1);
ofdmSymbol = zeros(usedN-ZTh-ZTt,size(dataModN,2));
ofdmSymbolZT = zeros(usedN,size(dataModN,2));
ofdmSymbolDFT = zeros(N,size(dataModN,2));

for j=1:size(dataModN,2)
    ofdmSymbol(:,j) = ifft(dataModUsedN(:,j),usedN-ZT);
    % Insert zero head and tail in the ofdmSymbol
    ofdmSymbolZT(:,j) = vertcat(zeros(ZTh,1),ofdmSymbol(:,j),zeros(ZTt,1));
    % DFT-s
    dataModUsedNZT = fft(ofdmSymbolZT,usedN);
    % Insert guard carriers keeping the data modulated in the middle
    dataModN(:,j) = ...
        vertcat(zeros(unusedN/2,1),dataModUsedNZT(:,j),zeros(unusedN/2,1));
    ofdmSymbolDFT(:,j) = ifft(dataModN(:,j),N);
    ofdm((j-1)*size(ofdmSymbolDFT)+1:j*size(ofdmSymbolDFT)) = ...
        ofdmSymbolDFT(:,j);
end

end

