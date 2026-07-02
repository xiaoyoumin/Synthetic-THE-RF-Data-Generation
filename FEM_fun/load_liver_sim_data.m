function [model, mesh, Uxyz, TimeVector] = load_liver_sim_data(h5file)
%LOAD_LIVER_SIM_DATA Read liver FEM mesh and motion data from an HDF5 file.
%
% The HDF5 file stores the anatomically realistic liver tetrahedral mesh
% and the time-dependent shear-wave displacement field used by the liver
% case study. This function reconstructs a MATLAB PDE model from the mesh
% so that FEM interpolation can be used in the scatterer-motion and
% ground-truth steps.
%
% Inputs:
%   h5file : Path to the liver simulation HDF5 file.
%
% Outputs:
%   model      : MATLAB PDE model with geometry created from the
%                tetrahedral mesh and region labels.
%   mesh       : Structure containing raw mesh arrays and convenience
%                fields nodes, elements, groupsID and source_freq.
%   Uxyz       : FEM displacement array arranged as
%                [num_nodes, 3 components, num_time_samples].
%   TimeVector : Row vector of FEM time samples [s].

meshNames = {'p','b','t','nv','nbe','nt','labels'};
mesh = struct();
for k = 1:numel(meshNames)
    mesh.(meshNames{k}) = h5read(h5file, sprintf('/mesh/%s', meshNames{k}));
end

mesh.nodes = mesh.p(1:3,:);
mesh.elements = mesh.t(1:4,:);
mesh.groupsID = mesh.t(5,:);

model = createpde;
geometryFromMesh(model, mesh.nodes, mesh.elements, mesh.groupsID);

Uxyz = h5read(h5file, '/Uxyz');
TimeVector = h5read(h5file, '/TimeVector');
TimeVector = TimeVector(:).';

try
    mesh.source_freq = h5read(h5file, '/source_freq');
catch
    mesh.source_freq = [];
end
