function sca_mesh = time_vary_sca_kwave(kgrid, sensor, disp_x, disp_y, disp_z, scatterers, TF)



sca_mesh = TF*[scatterers, ones(size(scatterers,1),1)]';
Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

for i = 1:size(disp_z,2)

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = disp_x(:,i);
    Dx_intrp = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = disp_y(:,i);
    Dy_intrp = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    disp_whole = zeros(Nx, Ny, Nz);
    disp_whole(logical(sensor.mask)) = disp_z(:,i);
    Dz_intrp = interpn(kgrid.x, kgrid.y, kgrid.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');
    
    % deal with NaN
    Dx_intrp(isnan(Dx_intrp)) = 0;
    Dy_intrp(isnan(Dy_intrp)) = 0;
    Dz_intrp(isnan(Dz_intrp)) = 0;

    sca_mesh(:,:,i+1) = sca_mesh(:,:,i) + [Dx_intrp; Dy_intrp; Dz_intrp; zeros(1, length(scatterers))];

end