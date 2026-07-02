function [impulse_response, excitation, lag] = pulse_setup(f0, pulse_duration, fb, fs)
%PULSE_SETUP Create excitation and impulse-response waveforms.
%
% This function defines the transmit excitation and Gaussian probe impulse
% response used by Field II for synthetic pulse-echo RF channel-data
% generation. It also estimates the two-way impulse-response lag used when
% assigning USTB transmit delays.
%
% Inputs:
%   f0             : Transducer center frequency [Hz].
%   pulse_duration : Pulse duration [cycles].
%   fb             : Fractional bandwidth of the probe.
%   fs             : RF sampling frequency [Hz].
%
% Outputs:
%   impulse_response : Gaussian impulse response waveform.
%   excitation       : Square-wave transmit excitation waveform.
%   lag              : Two-way impulse-response lag in RF samples.

dt = 1/fs;

fractional_bandwidth = fb;
t0 = (-1/fractional_bandwidth/f0):dt:(1/fractional_bandwidth/f0);
impulse_response = gauspuls(t0, f0, fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response);

te = (-pulse_duration/2/f0):dt:(pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response, excitation);
two_way_ir = conv(one_way_ir, impulse_response);
lag = length(two_way_ir)/2+1;
