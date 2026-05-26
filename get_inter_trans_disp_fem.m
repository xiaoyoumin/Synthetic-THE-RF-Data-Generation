function Disp_temp = get_inter_trans_disp_fem(Uxyz, TimeVector, t_current, t_next)


% method to shift the k-wave medium according to the simulated displacement
% field (velocity of each scatterer)

% check variable is vaild
if t_current < TimeVector(1)
    error('t_current out of range')
end

if  t_next > TimeVector(end)
    error('t_next out of range')
end


%%%%%%%%%%% Tackle with the interpolation in case     %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% t_current at the middle of two time step  %%%%%%%%%%%%%%%%%%%%%

% find the indexes of t_current, interpolate if need
if ~ismember(t_current, TimeVector)
    idx = t_current/DeltaT+1;
    idx_new = ceil(idx);
    idx_floor = floor(idx);
    U_temp = (1-(idx-idx_floor)) .* Uxyz(:,:, idx_floor) ...
        + (idx-idx_floor) .* Uxyz(:,:, idx_floor + 1);
    
    % find the medium displacement in [m]
    % Disp_temp = (TimeVector(idx_new) - t_current) .* U_temp;
    Disp_temp = U_temp;


else
    idx_new = floor(t_current/DeltaT)+1;
    Disp_temp = zeros(size(Uxyz, [1,2]));
end


%%%%%%%%%%%% Estimate the rest time step  %%%%%%%%%%%%%%%%%
for idx = idx_new : floor(t_next/DeltaT)+1

    % find the medium displacement in [m]
    Disp_temp = Disp_temp + Uxyz(:, :, idx);
end


%%%%%%%%%%% Tackle with the interpolation in case   %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% t_next at the middle of two time step   %%%%%%%%%%%%%%%%%%%%%

if ~ismember(t_next, TimeVector)
    idx = t_next/DeltaT+1;
    idx_floor = floor(idx);
    
    U_temp = (1-(idx-idx_floor)) .* Uxyz(:,:, idx_floor) ...
        + (idx-idx_floor) .* Uxyz(:,:, idx_floor + 1);
    
    % find the medium displacement in [m]
    Disp_temp = Disp_temp + U_temp;

end


