function [TF, TF_rev] = coor_transformation(prb_center, prb_theta)
% function [TF, TF_rev] = coor_transformation(prb_center, prb_theta)
%
% Define the coordinate transformation matrix, between probe coordinate and
% mesh coordinate.
%
% Inputs:
%   prb_center : position of the probe on the mesh, (x, y, z), [mm]
%   prb_theta : rotation of the probe, (yaw, roll), [rad]. Pitch rotation is not
%               allowed. Yaw: Rotate angle around vertical axis; roll:
%               rotation angle around long axis of the probe.
%
% Outputs:
%   TF : transformation matrix from probe coordinate to mesh coordinate.
%   TF_rev : transformation matrix from mesh coordinate to probe
%            coordinate.
%
% =========================================================================



prb_theta_x = prb_theta(1);
prb_theta_z = prb_theta(2);


% Define Transform matrix from probe coordinate to mesh coordinate

% Yaw rotation matrix
R_yaw = [cos(prb_theta_x), -sin(prb_theta_x), 0, 0;
         sin(prb_theta_x), cos(prb_theta_x), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];

% Roll rotation matrix
R_roll = [1, 0, 0, 0;
          0, cos(prb_theta_z), -sin(prb_theta_z), 0;
          0, sin(prb_theta_z), cos(prb_theta_z), 0;
          0,0,0,1];

% Translation matrix
T = [1,0,0,prb_center(1);
    0,1,0,prb_center(2);
    0,0,1,prb_center(3);
    0,0,0,1];

TF = T*R_yaw*R_roll;%*[scatterers, ones(size(xxp_speckle))]';

% Reverse transform matrix from mesh to probe

T_re = [1,0,0,-prb_center(1);
    0,1,0,-prb_center(2);
    0,0,1,-prb_center(3);
    0,0,0,1];

R_roll_re = [1, 0, 0, 0;
          0, cos(-prb_theta_z), -sin(-prb_theta_z), 0;
          0, sin(-prb_theta_z), cos(-prb_theta_z), 0;
          0,0,0,1];

R_yaw_re = [cos(-prb_theta_x), -sin(-prb_theta_x), 0, 0;
         sin(-prb_theta_x), cos(-prb_theta_x), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];

TF_rev = R_roll_re*R_yaw_re*T_re;%*sca_mesh;