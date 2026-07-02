function [TF, TF_rev] = coor_transformation(prb_center, prb_theta)
%COOR_TRANSFORMATION Build probe-to-mesh coordinate transforms.
%
% The ultrasound scatterers and image grid are first defined in probe
% coordinates. This function defines the rigid transform that places the
% probe coordinate system inside the shear-wave simulation mesh, together
% with the inverse transform used before Field II RF simulation.
%
% Inputs:
%   prb_center : Probe origin in mesh coordinates [x, y, z] [m].
%   prb_theta  : Probe rotation [yaw, roll] [rad]. Pitch rotation is not
%                used. Yaw rotates around the vertical axis and roll rotates
%                around the long axis of the probe.
%
% Outputs:
%   TF     : 4-by-4 homogeneous transform from probe coordinates to mesh
%            coordinates.
%   TF_rev : 4-by-4 homogeneous transform from mesh coordinates to probe
%            coordinates.

prb_theta_x = prb_theta(1);
prb_theta_z = prb_theta(2);

R_yaw = [cos(prb_theta_x), -sin(prb_theta_x), 0, 0;
         sin(prb_theta_x), cos(prb_theta_x), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];

R_roll = [1, 0, 0, 0;
          0, cos(prb_theta_z), -sin(prb_theta_z), 0;
          0, sin(prb_theta_z), cos(prb_theta_z), 0;
          0, 0, 0, 1];

T = [1, 0, 0, prb_center(1);
     0, 1, 0, prb_center(2);
     0, 0, 1, prb_center(3);
     0, 0, 0, 1];

TF = T*R_yaw*R_roll;

T_re = [1, 0, 0, -prb_center(1);
        0, 1, 0, -prb_center(2);
        0, 0, 1, -prb_center(3);
        0, 0, 0, 1];

R_roll_re = [1, 0, 0, 0;
             0, cos(-prb_theta_z), -sin(-prb_theta_z), 0;
             0, sin(-prb_theta_z), cos(-prb_theta_z), 0;
             0, 0, 0, 1];

R_yaw_re = [cos(-prb_theta_x), -sin(-prb_theta_x), 0, 0;
            sin(-prb_theta_x), cos(-prb_theta_x), 0, 0;
            0, 0, 1, 0;
            0, 0, 0, 1];

TF_rev = R_roll_re*R_yaw_re*T_re;
