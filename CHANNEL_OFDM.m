%--------------------------------------------------------------------------
%   Javier Celada Muñoz
%   Universidad Carlos III, Madrid
%   Grado Ingeniería de Sistemas Audiovisuales
%
%   Final Project
%   Simulation of waveforms for 5G systems  
%--------------------------------------------------------------------------
%   OFDM THROUGH MULTIPATH CHANNEL  
%   Simulation of the effect produced by this type of channel to a OFDM
%   signal.
%   
%   [ ofdmChannel ] = CHANNEL_OFDM( ofdm , H )
%
%   Funcitions returns ofdmEchoes (same dimension as ofdm) and needs 
%   this parameters input:
%   ofdm -> Chain of OFDM symbols (or one of them)
%   H-> Channel

function [ ofdmChannel ] = CHANNEL_OFDM( ofdm , H )

ofdmChannel = filter (H,1,ofdm);                    % OFDM through multi-path channel                                              

end