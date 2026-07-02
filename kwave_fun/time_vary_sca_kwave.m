function sca_mesh = time_vary_sca_kwave(kgrid, sensor, motion, scatterers, TF)
%TIME_VARY_SCA_KWAVE Encode k-Wave motion into moving scatterers.
%
% This function updates the initial scatterer coordinates according to the
% local PRF-sampled particle displacement from the k-Wave shear-wave
% simulation. The resulting time-varying scatterer distribution is the
% coupling representation used by Field II to generate synthetic RF channel
% data.
%
% Inputs:
%   kgrid      : kWaveGrid object used for interpolation coordinates.
%   sensor     : k-Wave sensor structure defining sampled displacement
%                nodes.
%   motion     : PRF-sampled motion structure with fields x, y and z.
%   scatterers : Initial scatterer coordinates in probe coordinates,
%                arranged as [num_scatterers, 3].
%   TF         : 4-by-4 transform from probe coordinates to k-Wave mesh
%                coordinates.
%
% Outputs:
%   sca_mesh : Homogeneous scatterer coordinates in mesh coordinates,
%              arranged as [4, num_scatterers, num_frames].

num_scat = size(scatterers, 1);
num_updates = size(motion.z, 2);
Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

sca_mesh = zeros(4, num_scat, num_updates + 1);
sca_mesh(:,:,1) = TF*[scatterers, ones(num_scat, 1)]';

for i = 1:num_updates
    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = motion.x(:,i);
    Dx_intrp = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = motion.y(:,i);
    Dy_intrp = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = motion.z(:,i);
    Dz_intrp = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    Dx_intrp(isnan(Dx_intrp)) = 0;
    Dy_intrp(isnan(Dy_intrp)) = 0;
    Dz_intrp(isnan(Dz_intrp)) = 0;

    sca_mesh(:,:,i+1) = sca_mesh(:,:,i) + ...
        [Dx_intrp; Dy_intrp; Dz_intrp; zeros(1, num_scat)];
end
