%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   OFDM THROUGH AWGN CHANNEL
%   Simulation of the effect produced by this type of channel to a OFDM
%   signal.
%
%   [ ofdmAWGN ] = AWGN_OFDM( EbN0_dB , ofdm, k )
%
%   Funcitions returns ofdmNoise (same dimension as ofdm) and needs 
%   this parameters input:
%
%   EbN0_dB ->  Bit energy to noise ratio in dB
%   ofdm ->     Chain of OFDM symbols (or one of them)
%   k->         Bits per QAM symbol
%   N->         Number of OFDM carriers
%   usedN->     Number of data carriers
%   CP->        Cyclic prefix lenght (or Zero tail)

function [ ofdmAWGN ] = AWGN_OFDM( EbN0_dB , ofdm, k , N , usedN , CP )

EbN0 = 10^(EbN0_dB/10);
snr = (N/(N+CP))*(usedN/N)*EbN0*k;
snrdB = 10*log10(snr);                              % SNR from EbN0
ofdmAWGN = awgn(ofdm,snrdB,'measured');             % Channel AWGN

% % AWGN addition figure
% figure
% plot(real(ofdmAWGN),'r')
% hold on
% plot(real(ofdm))
% plot(real(ofdmAWGN)-real(ofdm),'g')
% title('OFDM NOISE')
% legend('OFDM AWGN','OFDM','NOISE')

end

