function disp_gt = get_GT_fem(model, motion, vec_T, prb, scene_depth, TF)
%GET_GT_FEM Interpolate FEM motion onto the liver ultrasound image grid.
%
% The liver FEM displacement field is interpolated from the tetrahedral
% mesh onto a curvilinear ultrasound image grid and projected along the
% probe beam direction. This provides the known ground-truth motion field
% associated with the synthetic RF channel data.
%
% Inputs:
%   model       : MATLAB PDE model containing the liver tetrahedral mesh.
%   motion      : PRF-sampled FEM motion structure with fields x, y and z.
%   vec_T       : Time of each sampled ultrasound pulse interval [s].
%   prb         : Curvilinear probe structure with fields pitch, radius and
%                 N.
%   scene_depth : Axial imaging depth measured from the probe surface [m].
%   TF          : 4-by-4 transform from probe coordinates to FEM mesh
%                 coordinates.
%
% Outputs:
%   disp_gt : Ground-truth displacement projected onto the probe direction,
%             arranged as [depth samples, azimuth samples, time frames].
%   gt_x    : Lateral coordinates of the curvilinear image grid [m].
%   gt_z    : Axial coordinates of the curvilinear image grid [m].
%   PrbDir  : Unit vector of the probe beam direction in FEM coordinates.

grid_size = [256, 400]; % [azimuth samples, depth samples]

dtheta = 2*asin(prb.pitch/2/prb.radius);
max_angle = abs((prb.N-1)*dtheta)/2;

[theta, rho] = meshgrid(linspace(-max_angle, max_angle, grid_size(1)), ...
    linspace(prb.radius, prb.radius+scene_depth, grid_size(2)));

gt_x = rho.*sin(theta);
gt_z = rho.*cos(theta)-prb.radius;
gt_scatterers = [gt_x(:), zeros(size(gt_x(:))), gt_z(:)];
gt_mesh = TF*[gt_scatterers, ones(size(gt_x(:)))]';

Ux_sim = createPDEResults(model, motion.x, vec_T, "time-dependent");
Uy_sim = createPDEResults(model, motion.y, vec_T, "time-dependent");
Uz_sim = createPDEResults(model, motion.z, vec_T, "time-dependent");

time_idx = 1:length(vec_T);
Ux_gt = interpolateSolution(Ux_sim, gt_mesh(1:3,:), time_idx);
Uy_gt = interpolateSolution(Uy_sim, gt_mesh(1:3,:), time_idx);
Uz_gt = interpolateSolution(Uz_sim, gt_mesh(1:3,:), time_idx);

PrbDir = TF(1:3,1:3)*[0; 0; 1];
PrbDir = PrbDir(:).'/norm(PrbDir);

disp_gt = Ux_gt.*-PrbDir(1) + Uy_gt.*-PrbDir(2) + ...
    Uz_gt.*-PrbDir(3);
disp_gt = reshape(disp_gt, [size(gt_x), length(vec_T)]);
