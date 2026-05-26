function disp_gt = get_GT_kwave(kgrid, sensor, array_length, scene_depth, disp_x,disp_y, disp_z, prb_theta, TF)




prb_theta_x = prb_theta(1);
prb_theta_z = prb_theta(2);

% solve Normal of the scene plane
syms x y;
eq1 = x*cos(prb_theta_x)+y*sin(prb_theta_x) == 0;
% eq2 = norm([x,y,cos(pi/2-prb_theta_z)]) == 1;
eq2 = x^2+y^2+cos(pi/2-prb_theta_z)^2 == 1;
eqns = [eq1,eq2];
S = solve(eqns,[x y]);
N_scene = [double(S.x(2)), double(S.y(2)), cos(pi/2-prb_theta_z)];

% Solve the direction vector of the probe
syms x y;
eq1 = x*cos(prb_theta_x)+y*sin(prb_theta_x) == 0;
eq2 = N_scene(1)*x + N_scene(2)*y + N_scene(3) == 0;
eqns = [eq1,eq2];
S = solve(eqns,[x y]);
PrbDir = [double(S.x), double(S.y), 1];
PrbDir = PrbDir/norm(PrbDir);



Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

[X,Z] = meshgrid(-(array_length)/2:kgrid.dx:(array_length)/2, 5e-3:kgrid.dz:scene_depth);
gt_grid = TF*[X(:), zeros(size(X(:))), Z(:), ones(size(X(:)))]';

for i = 1:size(disp_z,2)

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = disp_x(:,i);
    Dx_gt = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = disp_y(:,i);
    Dy_gt = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = disp_z(:,i);
    Dz_gt = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');
    
    % deal with NaN
    Dx_gt(isnan(Dx_gt)) = 0;
    Dy_gt(isnan(Dy_gt)) = 0;
    Dz_gt(isnan(Dz_gt)) = 0;

    disp_gt(:,:,i) = reshape(Dx_gt.*-PrbDir(1) + Dy_gt.*-PrbDir(2) + Dz_gt.*-PrbDir(3), size(X));

end
