function sensor = create_sensor(kgrid, array_length, array_height, scene_depth, TF)
%CREATE_SENSOR Create the k-Wave sensor mask for the ultrasound image plane.
%
% The sensor mask samples the shear-wave field in a rectangular slab that
% corresponds to the pulse-echo imaging region of the linear probe. Points
% are defined in probe coordinates, transformed to the k-Wave mesh, and
% mapped to the nearest grid nodes.
%
% Inputs:
%   kgrid        : kWaveGrid object for the phantom shear-wave simulation.
%   array_length : Lateral width of the ultrasound aperture/scene [m].
%   array_height : Elevation height of the ultrasound probe [m].
%   scene_depth  : Axial depth of the imaging region [m].
%   TF           : 4-by-4 transform from probe coordinates to k-Wave
%                  mesh coordinates.
%
% Outputs:
%   sensor : k-Wave sensor structure. The mask selects the image slab and
%            sensor.record is set to record displacement components.

dx = kgrid.dx;
dy = kgrid.dy;
dz = kgrid.dz;
Nx = kgrid.Nx;
Ny = kgrid.Ny;
Nz = kgrid.Nz;

[mask_coor.x, mask_coor.y, mask_coor.z] = meshgrid(-(array_length+2e-3)/2:dx/15:(array_length+2e-3)/2, ...
    -array_height/2:dy/15:array_height/2, ...
    0:dz/15:scene_depth+2e-3);

coor_mesh = TF*[mask_coor.x(:), mask_coor.y(:), mask_coor.z(:), ones(size(mask_coor.x(:)))]';

mask_idx = knnsearch([kgrid.x(:), kgrid.y(:), kgrid.z(:)], coor_mesh(1:3,:)');

sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(mask_idx) = 1;

sensor.record = {'u'};
