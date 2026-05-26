function sensor = create_sensor(kgrid, array_length, array_height, scene_depth, TF)


dx = kgrid.dx;
dy = kgrid.dy;
dz = kgrid.dz;
Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

% sensor mask in probe coord
[mask_coor.x, mask_coor.y, mask_coor.z] = meshgrid(-(array_length+2e-3)/2:dx/15:(array_length+2e-3)/2, ...
    -array_height/2:dy/15:array_height/2, ...
    0:dz/15:scene_depth+2e-3);

% transform coord
coor_mesh = TF*[mask_coor.x(:), mask_coor.y(:), mask_coor.z(:), ones(size(mask_coor.x(:)))]';

% find the mesh grids that corresponding to the transformed coord
mask_idx = knnsearch([kgrid.x(:), kgrid.y(:), kgrid.z(:)], coor_mesh(1:3,:)');

sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(mask_idx) = 1;

sensor.record = {'u'};
