function v_inq = interpolate_FIELD_response(v, t, cropat, num_elements, dt)
%INTERPOLATE_FIELD_RESPONSE Align Field II output onto the RF sample grid.
%
% Field II returns RF samples with an initial time offset for each simulated
% scatterer response. This helper shifts and interpolates the channel data
% onto a fixed sample grid so that all frames and transmit events have the
% same RF data dimensions.
%
% Inputs:
%   v            : Raw Field II RF data arranged as
%                  [raw_samples, num_elements].
%   t            : Initial time returned by Field II [s].
%   cropat       : Number of RF samples to keep.
%   num_elements : Number of receive elements.
%   dt           : RF sampling interval [s].
%
% Outputs:
%   v_inq : RF data interpolated onto the common sample grid, arranged as
%           [cropat, num_elements].

t_vec = 0:dt:(size(v,1)-1)*dt;
t_inq = -t:dt:-t+(cropat-1)*dt;

v_inq = zeros(cropat, num_elements);
for i = 1:num_elements
    v_inq(:,i) = interp1(t_vec, v(:,i), t_inq, 'linear', 0);
end
