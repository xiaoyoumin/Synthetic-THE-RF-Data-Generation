function scatterers = create_scatterers(prb, scene_depth, pulse_duration, lambda, sca_per_cell)
%CREATE_SCATTERERS Generate the initial acoustic scatterer distribution.
%
% Scatterers are randomly distributed over the ultrasound imaging region to
% provide realistic speckle statistics for the mesh-free Field II
% simulation. The number of scatterers is controlled per approximate
% ultrasound resolution cell. Linear probes produce Cartesian coordinates;
% curvilinear probes produce fan-shaped coordinates in the probe frame.
%
% Inputs:
%   prb            : Probe structure. Required fields are N, pitch,
%                    element_width and element_height; curvilinear probes
%                    also include radius.
%   scene_depth    : Axial depth of the imaging region [m].
%   pulse_duration : Transmit pulse duration [cycles].
%   lambda         : Ultrasound wavelength [m].
%   sca_per_cell   : Number of scatterers per resolution cell.
%
% Outputs:
%   scatterers : Initial scatterer coordinates in probe coordinates,
%                arranged as [num_scatterers, 3].

if isfield(prb, 'radius')
    dtheta = 2*asin(prb.pitch/2/prb.radius);
    l_range = abs((prb.N-1)*dtheta);
    array_length = 2*prb.radius*sin(l_range/2);
    radius = prb.radius;
else
    l_range = (prb.N-1)*prb.pitch;
    array_length = l_range;
    radius = 0;
end

lateral_res = lambda/array_length;
axis_res = pulse_duration*lambda/2;

x_edges = -l_range/2:lateral_res:l_range/2;
z_edges = radius+1e-3:axis_res:radius+scene_depth;
num_scat = numel(x_edges)*numel(z_edges)*sca_per_cell;

xxp = zeros(num_scat, 1);
zzp = zeros(num_scat, 1);

idx = 1;
for x = x_edges
    for z = z_edges
        sca_idx = idx:idx+sca_per_cell-1;
        xxp(sca_idx) = random('unif', x, x+lateral_res, sca_per_cell, 1);
        zzp(sca_idx) = random('unif', z, z+axis_res, sca_per_cell, 1);
        idx = idx + sca_per_cell;
    end
end

if isfield(prb, 'radius')
    xxp_speckle = zzp.*sin(xxp);
    zzp_speckle = zzp.*cos(xxp)-prb.radius;
else
    xxp_speckle = xxp;
    zzp_speckle = zzp;
end

scatterers = [xxp_speckle, zeros(size(xxp_speckle)), zzp_speckle];
