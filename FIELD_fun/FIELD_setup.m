function [Th, Rh] = FIELD_setup(prb,fs,attenfreq,f0,c0,impulse_response, excitation)

lambda = c0/f0;

field_init(-1);

set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

set_field('att', attenfreq*100*f0/1e6);  %% db/m
set_field('freq_att',attenfreq*100/1e6);  %% dB/m/Hz
% set_field('freq_att', 0);
set_field('att_f0', f0);
set_field('use_att',1);

noSubAz=round(prb.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(prb.element_height/(lambda/8));       % number of subelements in the elevation direction

if isfield(prb,'radius')
    Th = xdc_convex_array (prb.N, prb.element_width, prb.element_height, prb.kerf, prb.radius, noSubAz, noSubEl, [0 0 Inf]); 
    Rh = xdc_convex_array (prb.N, prb.element_width, prb.element_height, prb.kerf, prb.radius, noSubAz, noSubEl, [0 0 Inf]); 
else
    Th = xdc_linear_array (prb.N, prb.element_width, prb.element_height, prb.kerf, noSubAz, noSubEl, [0 0 Inf]);
    Rh = xdc_linear_array (prb.N, prb.element_width, prb.element_height, prb.kerf, noSubAz, noSubEl, [0 0 Inf]);
end

xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);
