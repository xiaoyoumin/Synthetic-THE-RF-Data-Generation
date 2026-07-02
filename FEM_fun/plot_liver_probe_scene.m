function plot_liver_probe_scene(model, sca_mesh, liver_idx)
%PLOT_LIVER_PROBE_SCENE Visualize liver mesh and initial scatterer classes.
%
% This diagnostic plot shows the liver geometry together with scatterers
% classified as inside or outside the liver region. It is useful for
% checking the selected probe position and the initial backscatterer
% distribution before Field II RF simulation.
%
% Inputs:
%   model     : MATLAB PDE model containing the liver tetrahedral mesh.
%   sca_mesh  : Time-varying scatterer coordinates in mesh coordinates,
%               arranged as [4, num_scatterers, num_frames].
%   liver_idx : Logical mask identifying scatterers inside the liver cell.
%
% Outputs:
%   No output arguments. The function creates a figure.

figure
pdegplot(model, 'CellLabels', 'on', 'FaceAlpha', 0.2)
hold on

liver_sca = sca_mesh(:,liver_idx,1);
plot3(liver_sca(1,:), liver_sca(2,:), liver_sca(3,:), "ok", ...
    "MarkerFaceColor", "g", "MarkerSize", 0.1)
hold on

background_sca = sca_mesh(:,~liver_idx,1);
plot3(background_sca(1,:), background_sca(2,:), background_sca(3,:), "ok", ...
    "MarkerFaceColor", "r", "MarkerSize", 0.1)

axis equal
set(gca, 'ZDir', 'reverse')
