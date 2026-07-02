function [amp, incl] = define_phantom_scatterer_amplitudes(sca_mesh, shear_params)
%DEFINE_PHANTOM_SCATTERER_AMPLITUDES Assign stronger amplitudes near inclusions.
%
% This helper defines scatterer backscattering amplitudes for the phantom
% case study. Scatterers close to the spherical inclusion boundaries are
% assigned a stronger amplitude to make the inclusion interfaces visible in
% the synthetic pulse-echo data.
%
% Inputs:
%   sca_mesh     : Time-varying scatterer coordinates in mesh coordinates,
%                  arranged as [4, num_scatterers, num_frames].
%   shear_params : Phantom geometry structure containing inclusion centers,
%                  radii and phantom depth.
%
% Outputs:
%   amp  : Column vector of Field II scatterer amplitudes.
%   incl : Logical mask identifying scatterers near inclusion boundaries.

background_amp = 1;
inclusion_amp = 2;
boundary_width = 4e-4;

num_scat = size(sca_mesh, 2);
amp = background_amp * ones(num_scat, 1);

dist1 = sqrt((sca_mesh(1,:,1)-shear_params.cx1).^2 + ...
    (sca_mesh(2,:,1)-shear_params.cy1).^2 + ...
    (sca_mesh(3,:,1)-(shear_params.cz1-shear_params.zrange/2)).^2);
dist2 = sqrt((sca_mesh(1,:,1)-shear_params.cx2).^2 + ...
    (sca_mesh(2,:,1)-shear_params.cy2).^2 + ...
    (sca_mesh(3,:,1)-(shear_params.cz2-shear_params.zrange/2)).^2);

incl = (dist1 < shear_params.cr1+boundary_width & ...
    dist1 > shear_params.cr1-boundary_width) | ...
    (dist2 < shear_params.cr2+boundary_width & ...
    dist2 > shear_params.cr2-boundary_width);

amp(incl) = inclusion_amp;
