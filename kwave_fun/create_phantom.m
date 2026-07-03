function [kgrid, medium] = create_phantom(shear_params)
%CREATE_PHANTOM Build the k-Wave phantom shear-wave simulation model.
%
% This function creates the regular hexahedral k-Wave grid and assigns
% piecewise material properties for the phantom-mimicking case study. The
% background and two spherical inclusions define the shear-wave speed map
% used to simulate the time-dependent particle velocity field.
%
% Inputs:
%   shear_params : Structure with phantom geometry and mechanical
%                  parameters.
%       c_shear_bkg  : Background shear-wave speed [m/s].
%       c_shear_incl : Inclusion shear-wave speed [m/s].
%       source_freq  : External vibration frequency [Hz].
%       rho0         : Mass density [kg/m3].
%       xrange       : Phantom length in x [m].
%       yrange       : Phantom width in y [m].
%       zrange       : Phantom depth in z [m].
%       cx1, cy1, cz1 : Center of the first spherical inclusion [m].
%       cr1          : Radius of the first spherical inclusion [m].
%       cx2, cy2, cz2 : Center of the second spherical inclusion [m].
%       cr2          : Radius of the second spherical inclusion [m].
%       alpha_coeff_compression : Optional compression-wave attenuation.
%       alpha_coeff_shear       : Optional shear-wave attenuation.
%       reflect                 : Optional flag to add high-density
%                                 reflecting boundaries.
%       reflect_thickness       : Optional reflecting-boundary thickness in
%                                 grid points. Defaults to 5.
%
% Outputs:
%   kgrid  : kWaveGrid object defining the computational grid and time axis.
%   medium : k-Wave medium structure containing density, compressional-wave
%            speed, shear-wave speed and optional attenuation fields.

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

dx_shear = 1/source_freq/5;
Nz_shear = round(zrange/dx_shear);
Nx_shear = round(xrange/dx_shear);
Ny_shear = round(yrange/dx_shear);
kgrid = kWaveGrid(Nx_shear, dx_shear, Ny_shear, dx_shear, Nz_shear, dx_shear);

v_comp = 20;
medium.density = rho0 * ones(Nx_shear, Ny_shear, Nz_shear);
medium.sound_speed_compression = v_comp * ones(Nx_shear, Ny_shear, Nz_shear);
medium.sound_speed_shear = c_shear1 .* ones(Nx_shear, Ny_shear, Nz_shear);

cn_shear = sqrt((kgrid.x-cx1).^2 + (kgrid.y-cy1).^2 + ...
    (kgrid.z-(cz1-zrange/2)).^2) < cr1;
medium.sound_speed_shear(cn_shear) = c_shear2;

cn_shear = sqrt((kgrid.x-cx2).^2 + (kgrid.y-cy2).^2 + ...
    (kgrid.z-(cz2-zrange/2)).^2) < cr2;
medium.sound_speed_shear(cn_shear) = c_shear2;

if isfield(shear_params, 'alpha_coeff_compression')
    medium.alpha_coeff_compression = shear_params.alpha_coeff_compression;
end
if isfield(shear_params, 'alpha_coeff_shear')
    medium.alpha_coeff_shear = shear_params.alpha_coeff_shear;
end

if isfield(shear_params, 'reflect') && shear_params.reflect
    reflect_thickness = 5;
    if isfield(shear_params, 'reflect_thickness')
        reflect_thickness = shear_params.reflect_thickness;
    end
    reflect_thickness = max(0, round(reflect_thickness));

    nx_shell = min(reflect_thickness, floor((Nx_shear-1)/2));
    ny_shell = min(reflect_thickness, max(Ny_shear-1, 0));
    nz_shell = min(reflect_thickness, max(Nz_shear-1, 0));

    if nx_shell > 0
        medium.density(1:nx_shell,:,:) = 7850*10;
        medium.density(end-nx_shell+1:end,:,:) = 7850*10;
    end
    if ny_shell > 0
        medium.density(:,end-ny_shell+1:end,:) = 7850*10;
    end
    if nz_shell > 0
        medium.density(:,:,end-nz_shell+1:end) = 7850*10;
    end
end

t_end_shear = 0.1;
kgrid.makeTime(medium.sound_speed_compression, [], t_end_shear);
