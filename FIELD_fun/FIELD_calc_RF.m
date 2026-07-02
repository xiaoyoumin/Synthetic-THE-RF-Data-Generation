function RF = FIELD_calc_RF(prb, sca_mesh, amp, TF_rev, tx_sequence, scene_depth, fs, attenfreq, f0, c0, impulse_response, excitation)
%FIELD_CALC_RF Compute Field II RF data for moving scatterers.
%
% This function is the pulse-echo ultrasound simulation step of the
% framework. For each time-varying scatterer frame and each transmit event,
% it configures the Field II aperture, simulates received RF channel data,
% and resamples the result onto a common RF sample grid.
%
% Inputs:
%   prb              : Probe structure. Required fields include N,
%                      element_width, element_height, kerf and pitch;
%                      curvilinear probes also include radius.
%   sca_mesh         : Time-varying scatterer coordinates in mesh
%                      coordinates, arranged as
%                      [4, num_scatterers, num_frames].
%   amp              : Column vector of scatterer amplitudes.
%   TF_rev           : 4-by-4 transform from mesh coordinates to probe
%                      coordinates.
%   tx_sequence      : Transmit sequence. Linear probes use scalar steering
%                      angles [rad]; curvilinear probes use focal points
%                      arranged as [num_transmits, 3] [m].
%   scene_depth      : Axial imaging depth used to crop RF data [m].
%   fs               : RF sampling frequency [Hz].
%   attenfreq        : Acoustic attenuation [dB/MHz/cm].
%   f0               : Transducer center frequency [Hz].
%   c0               : Speed of sound [m/s].
%   impulse_response : Probe impulse response waveform.
%   excitation       : Transmit excitation waveform.
%
% Outputs:
%   RF : Simulated RF channel data arranged as
%        [samples, elements, transmits, frames].

dt = 1/fs;
F = size(sca_mesh, 3);
if isfield(prb, 'radius')
    Na = size(tx_sequence, 1);
else
    Na = numel(tx_sequence);
end
cropat = round(2*scene_depth/c0/dt);

RF = zeros(cropat, prb.N, Na, F, 'single');

disp('Field II: Computing RF dataset');
parfor f = 1:F
    disp(['Calculating frame ', num2str(f), ' of ', num2str(F)]);

    sca_prb = TF_rev*sca_mesh(:,:,f);
    sca_prb = sca_prb(1:3,:)';

    try
        [Th, Rh] = FIELD_setup(prb, fs, attenfreq, f0, c0, impulse_response, excitation);

        for n = 1:Na
            disp(['Calculating frame ', num2str(f), ' of ', num2str(F), ...
                ', transmit ', num2str(n), ' of ', num2str(Na)]);

            if isfield(prb, 'radius')
                FIELD_TRx_aperture(prb, Th, Rh, c0, tx_sequence(n,:)); %#ok<PFBNS>
            else
                FIELD_TRx_aperture(prb, Th, Rh, c0, tx_sequence(n));
            end

            [v, t] = calc_scat_multi(Th, Rh, sca_prb, amp);
            v_inq = interpolate_FIELD_response(v, t, cropat, prb.N, dt);

            RF(:,:,n,f) = single(v_inq);
        end

        field_end;
    catch ME
        try
            field_end;
        catch
        end
        error('FIELD_calc_RF:FieldIIFailed', ...
            'Field II simulation failed for frame %d. Unfinished RF entries remain zero. Error: %s', ...
            f, ME.message);
    end
end
