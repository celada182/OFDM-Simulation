%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   TRANSMISSION OF THE OFDM ZEROTAIL SIGNAL
%   Simulation of the creation of the OFDM signal from data bits input.
%   Distribution of the bits to create the QAM simbols and place them in the
%   OFDM carriers to form independent OFDM symbols with the zero tail.
%
%   [ ofdm,dataMod ] = TX_OFDM_ZEROTAIL( dataIn,M,N,usedN,ZT )
%
%   OUTPUT:
%   ofdm -> Chain of OFDM symbols
%   dataMod -> QAM simbols inside the OFDM carriers
%
%   INPUT:
%   dataIn -> Vector of bits
%   M -> Number of simbols of the QAM constellation
%   N -> Number of carriers OFDM
%   usedN -> Number of data carriers
%   ZT -> Length of the zero tail

function [ ofdm,dataMod ] = TX_OFDM_ZEROTAIL( dataIn,M,N,usedN,ZT )

k = log2(M);
unusedN = N - usedN;

dataInMatrix = reshape(dataIn,length(dataIn)/k,k);      % Reshape data into binary 4-tuples
dataSymbolsIn = bi2de(dataInMatrix);                    % Convert to integers

dataMod = qammod(dataSymbolsIn,M,0,'gray');

% Data modulated for the used carriers
dataModUsedN = reshape(dataMod,usedN,length(dataMod)/usedN);                
% Data added guard carriers
dataModN = zeros (N,size(dataModUsedN,2));                                  

ofdm = zeros(size(dataModN,2)*N,1);
ofdmSymbol = zeros(N,size(dataModN,2));
ofdmSymbolZT = zeros(N+ZT,size(dataModN,2));
zerotail = zeros (ZT,size(dataModN,2));

for j=1:size(dataModN,2)
    % Insert guard carriers keeping the data modulated in the middle
    dataModN(:,j) = ...
        vertcat(zeros(unusedN/2,1),dataModUsedN(:,j),zeros(unusedN/2,1));
    ofdmSymbol(:,j) = ifft(dataModN(:,j),N);
    % Insert zero tail at the end of the ofdmSymbol
    ofdmSymbolZT(:,j) = vertcat(ofdmSymbol(:,j),zerotail(:,j));
    ofdm((j-1)*size(ofdmSymbolZT)+1:j*size(ofdmSymbolZT)) = ...
        ofdmSymbolZT(:,j);
end

% %DEMO 4
% figure
% stem(real(dataModN(:,1)),'r')
% title('Carriers of first OFDM symbol')
% xlabel('Carriers')
% ylabel('Real part of QAM symbol')
% 
% figure
% plot(real(ofdm),'g')
% title('Complete OFDM signal')
% xlabel('Tiem samples')
% ylabel('Real part')
% 
% figure
% plot(real(ofdm(1:N+ZT)))
% title('First OFDM symbol')
% xlabel('Tiem samples')
% ylabel('Real part')
% 
% keyboard

end

