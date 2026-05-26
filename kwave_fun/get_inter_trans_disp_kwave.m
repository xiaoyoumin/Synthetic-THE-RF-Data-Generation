function displace = get_inter_trans_disp_kwave(sensor_data_shear, kgrid, t_current, t_next)




% method to shift the k-wave medium according to the simulated displacement
% field (velocity of each scatterer)

% check variable is vaild
if t_current < kgrid.t_array(1)
    error('t_current out of range')
end

if  t_next > kgrid.t_array(end)
    error('t_next out of range')
end

data_size = 1:size(sensor_data_shear.ux, 1);
%%%%%%%%%%% Tackle with the interpolation in case     %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% t_current at the middle of two time step  %%%%%%%%%%%%%%%%%%%%%

% find the indexes of t_current, interpolate if need
if ~ismember(t_current, kgrid.t_array)

    displacement_current.ux = interpn(data_size, kgrid.t_array, sensor_data_shear.ux, data_size, t_current);
    displacement_current.uy = interpn(data_size, kgrid.t_array, sensor_data_shear.uy, data_size, t_current);
    displacement_current.uz = interpn(data_size, kgrid.t_array, sensor_data_shear.uz, data_size, t_current);
    
    idx_new = ceil(t_current/kgrid.dt);
    
    % find the medium displacement in [m]
    displace.x = (kgrid.t_array(idx_new) - t_current) .* displacement_current.ux;
    displace.y = (kgrid.t_array(idx_new) - t_current) .* displacement_current.uy;
    displace.z = (kgrid.t_array(idx_new) - t_current) .* displacement_current.uz;


else
    idx_new = floor(t_current/kgrid.dt);
    displace.x = zeros(data_size,1);
    displace.y = zeros(data_size,1);
    displace.z = zeros(data_size,1);
end


%%%%%%%%%%%% Estimate the rest time step  %%%%%%%%%%%%%%%%%
for idx = idx_new : floor(t_next/kgrid.dt)

    % find the medium displacement in [m]
    displace.x = displace.x + kgrid.dt .* sensor_data_shear.ux(:, idx);
    displace.y = displace.y + kgrid.dt .* sensor_data_shear.uy(:, idx);
    displace.z = displace.z + kgrid.dt .* sensor_data_shear.uz(:, idx);

end


%%%%%%%%%%% Tackle with the interpolation in case   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% t_next at the middle of two time step   %%%%%%%%%%%%%%%%%%%%%

if ~ismember(t_next, kgrid.t_array)

    displacement_current.ux = interpn(data_size, kgrid.t_array, sensor_data_shear.ux, data_size, t_next);
    displacement_current.uy = interpn(data_size, kgrid.t_array, sensor_data_shear.uy, data_size, t_next);
    displacement_current.uz = interpn(data_size, kgrid.t_array, sensor_data_shear.uz, data_size, t_next);

    idx_floor = floor(t_next/kgrid.dt);
    
    % find the medium displacement in [m]
    displace.x = displace.x + (t_next - kgrid.t_array(idx_floor)) .* displacement_current.ux;
    displace.y = displace.y + (t_next - kgrid.t_array(idx_floor)) .* displacement_current.uy;
    displace.z = displace.z + (t_next - kgrid.t_array(idx_floor)) .* displacement_current.uz;

end


%%%%%%%%%% finer the computation grid of shear wave result %%%%%%%%%%%%%%%
% the grid used for bmode imaging is much finer than the shear wave
% simulation grid. Shear wave result need to be extended to match the size
% of the bmode grid or scene. On the other hand, the displacement usually
% very small (around 10^-7 m), so the medium grid need to finner to this level in
% order to perform the move

% displace.x = interpn(kgrid_shear.x, kgrid_shear.y, displace.x, kgrid.x, kgrid.y, 'linear');
% displace.z = interpn(kgrid_shear.x, kgrid_shear.y, displace.z, kgrid.x, kgrid.y, 'linear');
% 
% deal with NaN
displace.x(isnan(displace.x)) = 0;
displace.y(isnan(displace.y)) = 0;
displace.z(isnan(displace.z)) = 0;