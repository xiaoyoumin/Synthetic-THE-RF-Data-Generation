function disp_interval = get_inter_trans_disp(motion_data, time_vector, ...
    t_current, t_next, sample_scale, endpoint_scale)
%GET_INTER_TRANS_DISP Sample motion over one ultrasound pulse interval.
%
% This helper accumulates the shear-wave motion that occurs between two
% ultrasound transmissions. It follows the same interval accumulation
% pattern used in the original phantom and liver scripts: interpolate the
% motion at the start and end of the pulse interval when needed, and add
% the samples that lie between the two transmit times.
%
% Inputs:
%   motion_data : Motion array arranged as
%                 [num_spatial_samples, 3 components, num_time_samples].
%   time_vector : Time samples corresponding to the third dimension of
%                 motion_data [s].
%   t_current   : Start time of the pulse interval [s].
%   t_next      : End time of the pulse interval [s].
%   sample_scale  : Displacement contribution of one full stored sample.
%                   For velocity samples, use the simulation time step. For
%                   stored displacement increments, use 1. Defaults to the
%                   time step.
%   endpoint_scale: Optional scale for interpolated endpoint samples. Use
%                   [] to weight endpoints by the fractional interval
%                   duration, as in the original phantom code. Use 1 to add
%                   interpolated endpoint samples directly, as in the
%                   original liver code.
%
% Outputs:
%   disp_interval : Displacement accumulated over the pulse interval,
%                   arranged as [num_spatial_samples, 3 components].

if nargin < 5
    sample_scale = [];
end
if nargin < 6
    endpoint_scale = [];
end

time_vector = time_vector(:).';

if t_current < time_vector(1)
    error('t_current out of range')
end

if t_next > time_vector(end)
    error('t_next out of range')
end

delta_t = time_vector(2)-time_vector(1);
if isempty(sample_scale)
    sample_scale = delta_t;
end

if any(abs(diff(time_vector)-delta_t) > 100*eps(max(abs(time_vector))))
    error('time_vector must be uniformly sampled.')
end

disp_interval = accumulate_interval(motion_data, time_vector, ...
    t_current, t_next, sample_scale, endpoint_scale, delta_t);
disp_interval(isnan(disp_interval)) = 0;

end

function disp_interval = accumulate_interval(motion_data, time_vector, ...
    t_current, t_next, sample_scale, endpoint_scale, delta_t)
num_points = size(motion_data, 1);
num_components = size(motion_data, 2);

disp_interval = zeros(num_points, num_components, 'like', motion_data);

idx_current = fractional_time_index(t_current, time_vector);
idx_next = fractional_time_index(t_next, time_vector);

if ~is_time_sample(t_current, time_vector)
    idx_new = ceil(idx_current);
    start_weight = get_endpoint_weight(time_vector(idx_new)-t_current, ...
        sample_scale, endpoint_scale, delta_t);
    disp_interval = start_weight .* ...
        interpolate_motion_sample(motion_data, idx_current, floor(idx_current));
else
    idx_new = find_time_index(t_current, time_vector);
end

if is_time_sample(t_next, time_vector)
    idx_last = find_time_index(t_next, time_vector);
else
    idx_last = floor(idx_next);
end

for idx = idx_new:idx_last
    disp_interval = disp_interval + sample_scale .* motion_data(:,:,idx);
end

if ~is_time_sample(t_next, time_vector)
    idx_floor = floor(idx_next);
    end_weight = get_endpoint_weight(t_next-time_vector(idx_floor), ...
        sample_scale, endpoint_scale, delta_t);
    disp_interval = disp_interval + ...
        end_weight .* interpolate_motion_sample(motion_data, idx_next, idx_floor);
end
end

function weight = get_endpoint_weight(duration, sample_scale, endpoint_scale, delta_t)
if isempty(endpoint_scale)
    weight = duration/delta_t * sample_scale;
else
    weight = endpoint_scale;
end
end

function tf = is_time_sample(t, time_vector)
tol = 10*eps(max(1, max(abs(time_vector))));
tf = any(abs(time_vector-t) <= tol);
end

function idx = find_time_index(t, time_vector)
tol = 10*eps(max(1, max(abs(time_vector))));
idx = find(abs(time_vector-t) <= tol, 1, 'first');
end

function idx = fractional_time_index(t, time_vector)
delta_t = time_vector(2)-time_vector(1);
idx = (t-time_vector(1))/delta_t+1;
end

function data_sample = interpolate_motion_sample(motion_data, idx, idx_floor)
frac = idx-idx_floor;
data_sample = (1-frac) .* motion_data(:,:,idx_floor) ...
    + frac .* motion_data(:,:,idx_floor+1);
end
