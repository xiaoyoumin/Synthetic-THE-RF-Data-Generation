function source = create_source(kgrid, source_freq)

% Create source on the left side of the mesh and define the sinusoidal
% signal.

Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

% create k-Wave source object
source.s_mask = zeros(Nx, Ny, Nz);
source.s_mask(:, 1, :) = 1; % source on the left

% Define source signal
source_mag_1 = 2; % stress [Pa]

% Calculate the angular frequency from the source frequency
angular_freq = 2 * pi * source_freq;

% shear stress
window = 1;
source.syz = [];
source.syz ...
    =repmat( source_mag_1 * window * sin(angular_freq * kgrid.t_array), sum(source.s_mask, 'all'),1);
source.syz = source.syz(any(source.syz, 2), :);
source.sxz = source.syz;
