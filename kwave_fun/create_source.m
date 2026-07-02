function source = create_source(kgrid, source_freq)
%CREATE_SOURCE Define the harmonic shear-wave source for k-Wave.
%
% The source represents the external mechanical vibrator used in
% time-harmonic elastography. It applies a sinusoidal shear stress on the
% left side of the phantom grid.
%
% Inputs:
%   kgrid       : kWaveGrid object with the simulation time array.
%   source_freq : Vibration frequency of the harmonic source [Hz].
%
% Outputs:
%   source : k-Wave source structure containing the source mask and the
%            sinusoidal shear-stress time series.

Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

source.s_mask = zeros(Nx, Ny, Nz);
source.s_mask(:, 1, :) = 1;

source_mag_1 = 2; % stress [Pa]
angular_freq = 2 * pi * source_freq;

window = 1;
source.syz = repmat(source_mag_1 * window * sin(angular_freq * kgrid.t_array), ...
    sum(source.s_mask, 'all'), 1);
source.syz = source.syz(any(source.syz, 2), :);
source.sxz = source.syz;
