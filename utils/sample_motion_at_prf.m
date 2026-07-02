function [motion, vec_T] = sample_motion_at_prf(model_type, sim_data, time_info, prf, t_start)
%SAMPLE_MOTION_AT_PRF Sample simulated motion over ultrasound pulse periods.
%
% The ultrasound RF simulation uses one scatterer distribution per pulse
% transmission. This function converts the shear-wave simulation output into
% PRF-sampled displacement increments for those transmission times. It
% supports both k-Wave particle-velocity output and FEM displacement output.
%
% Inputs:
%   model_type : "kwave" for k-Wave sensor data or "fem" for FEM data.
%   sim_data   : k-Wave sensor_data_shear structure or FEM Uxyz array
%                arranged as [num_nodes, 3 components, num_time_samples].
%   time_info  : kWaveGrid object for k-Wave data or FEM TimeVector [s].
%   prf        : Ultrasound pulse repetition frequency [Hz].
%   t_start    : Time of the first ultrasound pulse [s]. Defaults to 0.
%
% Outputs:
%   motion : Structure with fields x, y and z, each arranged as
%            [num_spatial_samples, num_pulse_intervals].
%   vec_T  : Start time of each sampled pulse interval [s].

if nargin < 5
    t_start = 0;
end

model_type = lower(string(model_type));

switch model_type
    case "kwave"
        motion_data = permute(cat(3, sim_data.ux, sim_data.uy, sim_data.uz), [1 3 2]);
        time_vector = time_info.t_array;
        sample_scale = time_vector(2)-time_vector(1);
        endpoint_scale = [];
    case "fem"
        motion_data = sim_data;
        time_vector = time_info;
        sample_scale = 1;
        endpoint_scale = 1;
    otherwise
        error('model_type must be "kwave" or "fem".')
end

time_vector = time_vector(:).';
t_end = time_vector(end);

t_current = t_start;
t_next = t_current + 1/prf;

motion.x = [];
motion.y = [];
motion.z = [];
vec_T = [];

while t_next < t_end
    disp_interval = get_inter_trans_disp(motion_data, time_vector, ...
        t_current, t_next, sample_scale, endpoint_scale);

    motion.x = cat(2, motion.x, disp_interval(:,1));
    motion.y = cat(2, motion.y, disp_interval(:,2));
    motion.z = cat(2, motion.z, disp_interval(:,3));
    vec_T(end+1) = t_current; %#ok<AGROW>

    t_current = t_next;
    t_next = t_current + 1/prf;
end
