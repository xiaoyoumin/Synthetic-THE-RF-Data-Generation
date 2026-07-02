function sca_mesh = time_vary_sca_fem(model, motion, vec_T, scatterers, TF)
%TIME_VARY_SCA_FEM Encode FEM motion into moving scatterers.
%
% This function updates the initial scatterer coordinates according to the
% local FEM displacement interpolated from the tetrahedral liver mesh. The
% resulting time-varying scatterer distribution links the shear-wave FEM
% simulation to the Field II pulse-echo RF simulation.
%
% Inputs:
%   model      : MATLAB PDE model containing the liver tetrahedral mesh.
%   motion     : PRF-sampled FEM motion structure with fields x, y and z.
%   vec_T      : Time of each sampled ultrasound pulse interval [s].
%   scatterers : Initial scatterer coordinates in probe coordinates,
%                arranged as [num_scatterers, 3].
%   TF         : 4-by-4 transform from probe coordinates to FEM mesh
%                coordinates.
%
% Outputs:
%   sca_mesh : Homogeneous scatterer coordinates in mesh coordinates,
%              arranged as [4, num_scatterers, num_frames].

num_scat = size(scatterers, 1);
num_updates = length(vec_T);

motion_sim.x = createPDEResults(model, motion.x, vec_T, "time-dependent");
motion_sim.y = createPDEResults(model, motion.y, vec_T, "time-dependent");
motion_sim.z = createPDEResults(model, motion.z, vec_T, "time-dependent");

sca_mesh = zeros(4, num_scat, num_updates + 1);
sca_mesh(:,:,1) = TF*[scatterers, ones(num_scat, 1)]';

for i = 1:num_updates
    Dx_intrp = interpolateSolution(motion_sim.x, sca_mesh(1:3,:,i), i);
    Dy_intrp = interpolateSolution(motion_sim.y, sca_mesh(1:3,:,i), i);
    Dz_intrp = interpolateSolution(motion_sim.z, sca_mesh(1:3,:,i), i);

    Dx_intrp(isnan(Dx_intrp)) = 0;
    Dy_intrp(isnan(Dy_intrp)) = 0;
    Dz_intrp(isnan(Dz_intrp)) = 0;

    sca_mesh(:,:,i+1) = sca_mesh(:,:,i) + ...
        [Dx_intrp'; Dy_intrp'; Dz_intrp'; zeros(1, num_scat)];
end
