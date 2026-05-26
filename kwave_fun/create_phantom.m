function [kgrid, medium] = create_phantom(shear_params)
% function [kgrid, medium] = create_phantom(shear_params)
%
% Create basic phantom model in k-Wave.
%
% Inputs:
%   shear_params : Structure containing parameters for phantom geometry and
%                  material physical properties:
%       - c_shear_bkg                            : background shear wave speed [m/s]
%       - c_shear_incl                           : inclusion shear wave speed [m/s]
%       - source_freq                            : shaker frequency [Hz]
%       - rho0                                   : medium density [kg/m3]
%       - xrange                                 : mesh length [m]
%       - zrange                                 : mesh height [m]
%       - yrange                                 : mesh width  [m]
%       - cx1, cy1, cz1, cr1, cx2, cy2, cz2, cr2 : center position and 
%                                                  radius of the inclusions
%       - alpha_coeff_compression                : compression wave
%                                                  attenuation (option)
%       - alpha_coeff_shear                      : shear wave attenuation
%                                                  (option)
%
% Outputs:
%   kgrid : k-Wave grid object
%   medium : k-Wave medium object
%
% =========================================================================




c_shear1 = shear_params.c_shear_bkg;
c_shear2 = shear_params.c_shear_incl;
source_freq = shear_params.source_freq;
rho0 = shear_params.rho0;
xrange = shear_params.xrange;
zrange = shear_params.zrange;
yrange = shear_params.yrange;
cx1 = shear_params.cx1;
cy1 = shear_params.cy1;
cz1 = shear_params.cz1;
cr1 = shear_params.cr1;
cx2 = shear_params.cx2;
cy2 = shear_params.cy2;
cz2 = shear_params.cz2;
cr2 = shear_params.cr2;



% set shear wave simulation grid
dx_shear = 1/source_freq/10; % use the lowest possible shear wave speed 1m/s to calculate wavelength
Nz_shear = round(zrange/dx_shear); % [grid points]
Nx_shear = round(xrange/dx_shear); % [grid points]
Ny_shear = round(yrange/dx_shear); % [grid points]
kgrid = kWaveGrid(Nx_shear, dx_shear, Ny_shear, dx_shear, Nz_shear, dx_shear);

% define shear medium
v_comp = 20; % compression wave speed used in shear wave simulation
medium.density = rho0 * ones(Nx_shear, Ny_shear, Nz_shear);
medium.sound_speed_compression = v_comp * ones(Nx_shear, Ny_shear, Nz_shear);
medium.sound_speed_shear = c_shear1 .* ones(Nx_shear, Ny_shear, Nz_shear);

% inclusion 1
cn_shear=sqrt((kgrid.x-cx1).^2+(kgrid.y-(cy1)).^2+(kgrid.z-(cz1-zrange/2)).^2)<cr1;
medium.sound_speed_shear(cn_shear) = c_shear2;

% inclusion 2
cn_shear=sqrt((kgrid.x-cx2).^2+(kgrid.y-(cy2)).^2+(kgrid.z-(cz2-zrange/2)).^2)<cr2;
medium.sound_speed_shear(cn_shear) = c_shear2;

% define the absorption properties
if isfield(shear_params, 'alpha_coeff_compression')
    medium.alpha_coeff_compression = shear_params.alpha_coeff_compression; % [dB/(MHz^2 cm)]
end
if isfield(shear_params, 'alpha_coeff_shear')
    medium.alpha_coeff_shear = shear_params.alpha_coeff_compression; % [dB/(MHz^2 cm)]
end

% make simulation time series
t_end_shear =0.1; % [s]
kgrid.makeTime(medium.sound_speed_compression, [], t_end_shear);

