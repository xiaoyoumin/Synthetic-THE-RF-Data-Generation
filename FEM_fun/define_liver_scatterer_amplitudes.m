function [amp, liver_idx] = define_liver_scatterer_amplitudes(model, sca_mesh)
%DEFINE_LIVER_SCATTERER_AMPLITUDES Assign stronger amplitudes inside liver.
%
% This helper classifies the initial scatterer positions using the liver
% cell of the tetrahedral FEM geometry. Scatterers inside the liver region
% are assigned stronger backscattering amplitudes than background
% scatterers for the synthetic pulse-echo simulation.
%
% Inputs:
%   model    : MATLAB PDE model containing the liver tetrahedral mesh and
%              cell labels.
%   sca_mesh : Time-varying scatterer coordinates in mesh coordinates,
%              arranged as [4, num_scatterers, num_frames].
%
% Outputs:
%   amp       : Column vector of Field II scatterer amplitudes.
%   liver_idx : Logical mask identifying scatterers located inside the
%               liver cell in the initial frame.

liver_cell_id = 2;
background_amp = 1;
liver_amp = 5;

num_scat = size(sca_mesh, 2);
amp = background_amp * ones(num_scat, 1);

Ef = findElements(model.Mesh, "region", Cell=liver_cell_id);
elements = model.Mesh.Elements(:,Ef);
node_idx = unique(elements(:));
[~, local_elements] = ismember(elements, node_idx);
TR = triangulation(local_elements', model.Mesh.Nodes(:,node_idx)');
SI = pointLocation(TR, sca_mesh(1:3,:,1)');

liver_idx = ~isnan(SI);
amp(liver_idx) = liver_amp;
