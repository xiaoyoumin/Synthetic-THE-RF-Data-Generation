function [Th, Rh] = FIELD_setup(prb, fs, attenfreq, f0, c0, impulse_response, excitation)
%FIELD_SETUP Initialize Field II and create probe aperture handles.
%
% This helper sets the global Field II acoustic parameters, builds either a
% linear or convex array from the probe structure, and assigns the transmit
% excitation and impulse responses required for pulse-echo RF simulation.
%
% Inputs:
%   prb              : Probe structure. Required fields include N,
%                      element_width, element_height, kerf and pitch;
%                      curvilinear probes also include radius.
%   fs               : RF sampling frequency [Hz].
%   attenfreq        : Acoustic attenuation [dB/MHz/cm]. Use 0 or [] to
%                      disable attenuation.
%   f0               : Transducer center frequency [Hz].
%   c0               : Speed of sound [m/s].
%   impulse_response : Probe impulse response waveform.
%   excitation       : Transmit excitation waveform.
%
% Outputs:
%   Th : Field II transmit aperture handle.
%   Rh : Field II receive aperture handle.

lambda = c0/f0;

field_init(-1);

set_field('c', c0);
set_field('fs', fs);
set_field('use_rectangles', 1);

if ~isempty(attenfreq) && attenfreq > 0
    set_field('att', attenfreq*100*f0/1e6);
    set_field('freq_att', attenfreq*100/1e6);
    set_field('att_f0', f0);
    set_field('use_att', 1);
else
    set_field('use_att', 0);
end

noSubAz = round(prb.element_width/(lambda/8));
noSubEl = round(prb.element_height/(lambda/8));

if isfield(prb, 'radius')
    Th = xdc_convex_array(prb.N, prb.element_width, prb.element_height, ...
        prb.kerf, prb.radius, noSubAz, noSubEl, [0 0 Inf]);
    Rh = xdc_convex_array(prb.N, prb.element_width, prb.element_height, ...
        prb.kerf, prb.radius, noSubAz, noSubEl, [0 0 Inf]);
else
    Th = xdc_linear_array(prb.N, prb.element_width, prb.element_height, ...
        prb.kerf, noSubAz, noSubEl, [0 0 Inf]);
    Rh = xdc_linear_array(prb.N, prb.element_width, prb.element_height, ...
        prb.kerf, noSubAz, noSubEl, [0 0 Inf]);
end

xdc_excitation(Th, excitation);
xdc_impulse(Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th, [0 0 0]);
xdc_impulse(Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh, [0 0 0]);
