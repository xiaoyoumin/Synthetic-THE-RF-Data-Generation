

clear;
% close all;

addpath('/fp/homes01/u01/ec-chaoranh/Documents/MATLAB/ustb')
addpath('/fp/homes01/u01/ec-chaoranh/Documents/MATLAB/field_II/')


%% General setting

DATA_CAST       = 'gpuArray-single';     % set to 'single' or 'gpuArray-single' to speed up computations

prf = 3000;
source_freq = 200; % [Hz]
alpha_max = deg2rad(15);%atan(1/2/F_number);
Na=1;                                      % number of plane waves 

SAVE_DIR = strcat('/fp/projects01/ec35/homes/ec-chaoranh/data/uffc_3D_phantom_sim4_attenu_',num2str(source_freq),'Hz.mat');

c0=1540;     % Speed of sound [m/s]
rho0 = 1079;    % medium density [kg/m3]


%% Transducer definition L7-4v, 128-element linear array transducer

probe = uff.linear_array();
f0                      = 5.1333e+06;      % Transducer center frequency [Hz]
lambda                  = c0/f0;           % Wavelength [m]
probe.element_height    = 5e-3;            % Height of element [m]
probe.pitch             = 0.300e-3;        % probe.pitch [m]
kerf                    = 0.03e-03;        % gap between elements [m]
probe.element_width     = probe.pitch-kerf;% Width of element [m]
lens_el                 = 20e-3;           % position of the elevation focus
probe.N                 = 128;             % Number of elements
pulse_duration          = 2.5;             % pulse duration [cycles]


% probe = uff.linear_array();
% f0                      = 4464430;      % Transducer center frequency [Hz]
% lambda                  = c0/f0;           % Wavelength [m]
% probe.element_height    = 7.5e-3;            % Height of element [m]
% probe.pitch             = 0.298e-3;        % probe.pitch [m]
% probe.element_width     = 0.25e-3;% Width of element [m]
% kerf                    = probe.pitch-probe.element_width;        % gap between elements [m]
% lens_el                 = 25e-3;           % position of the elevation focus
% probe.N                 = 128;             % Number of elements
% pulse_duration          = 2.5;             % pulse duration [cycles]



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Shear wave definition                 %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Setting scene
xrange = 100e-3;
zrange = 95e-3;
yrange = 120e-3;

%% shear wave velocity 
c_shear1 = 2.32; 
c_shear2 = 4.67;

v_comp = 20; % compression wave speed used in shear wave simulation


%% Set the frequency of the source

% Define the magnitude of the source
source_mag_1 = 2; % stress [Pa]


%% set shear wave simulation grid
dx_shear = 1/source_freq/10; % use the lowest possible shear wave speed 1m/s to calculate wavelength
Nz_shear = round(zrange/dx_shear); % [grid points]
Nx_shear = round(xrange/dx_shear); % [grid points]
Ny_shear = round(yrange/dx_shear); % [grid points]

kgrid_shear = kWaveGrid(Nx_shear, dx_shear, Ny_shear, dx_shear, Nz_shear, dx_shear);

%% define shear medium
medium_shear.density = rho0 * ones(Nx_shear, Ny_shear, Nz_shear);
medium_shear.sound_speed_compression = v_comp * ones(Nx_shear, Ny_shear, Nz_shear);
medium_shear.sound_speed_shear = c_shear1 .* ones(Nx_shear, Ny_shear, Nz_shear);

% include hyperechoic cyst
cx1=0; cy1=40e-3; cz1=15e-3; cr1 = 5e-3;
cn_shear=sqrt((kgrid_shear.x-cx1).^2+(kgrid_shear.y-(cy1-yrange/2)).^2+(kgrid_shear.z-(cz1-zrange/2)).^2)<cr1;
medium_shear.sound_speed_shear(cn_shear) = c_shear2;

% include hyperechoic cyst 2
cx2=0; cy2=80e-3; cz2=35e-3; cr2 = 10e-3;
cn_shear=sqrt((kgrid_shear.x-cx2).^2+(kgrid_shear.y-(cy2-yrange/2)).^2+(kgrid_shear.z-(cz2-zrange/2)).^2)<cr2;
medium_shear.sound_speed_shear(cn_shear) = c_shear2;

% define the absorption properties
medium.alpha_coeff_compression = 1; % [dB/(MHz^2 cm)]
medium.alpha_coeff_shear       = 1; % [dB/(MHz^2 cm)]


%% make simulation time series
t_end_shear =0.1;% 0.064; % [s]
kgrid_shear.makeTime(medium_shear.sound_speed_compression, [], t_end_shear);


%% create shear wave source

source_shear.s_mask = zeros(Nx_shear, Ny_shear, Nz_shear);
source_shear.s_mask(:, 1, :) = 1; % source on the left

%% Shear source

% Calculate the angular frequency from the source frequency
angular_freq = 2 * pi * source_freq;

% shear stress
window = 1;%tukeywin(length(source_shear.s_mask(end, :)), 0.75);
% source_shear.sxy(sub2ind(size(source_shear.s_mask), ones(1,kgrid_shear.Nx)*kgrid_shear.Nx, 1:kgrid_shear.Ny), :) ...
%     =repmat( source_mag_1 * window * sin(angular_freq * kgrid_shear.t_array), kgrid_shear.Ny,1);
source_shear.syz = [];
source_shear.syz ...
    =repmat( source_mag_1 * window * sin(angular_freq * kgrid_shear.t_array), sum(source_shear.s_mask, 'all'),1);
% source_shear.sxy(sub2ind(size(source_shear.s_mask), aaaa, bbbb), :) = ...
%     repmat (6*source_mag_1 * window * sin(angular_freq * kgrid_shear.t_array), length(bbbb), 1);
source_shear.syz = source_shear.syz(any(source_shear.syz, 2), :);
source_shear.sxz = source_shear.syz;

% Set the source pressure fields in the x and y directions to zero
% source_shear.sxx = 0*source_shear.syz;
% source_shear.syy = 0*source_shear.syz;


%% define scene

minZ = kgrid_shear.z_vec(1);

% Center of the aperture
prb_center = [0 20e-3 minZ]; % use the lower bigger inclusion

scene_depth = 60e-3;

% Direction of probe
prb_theta_x = deg2rad(90);  % 90 degree
prb_theta_y = pi/2-prb_theta_x;
prb_theta_z = deg2rad(0);

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


% Array geometry
array_length = probe.x(end)-probe.x(1);
array_height = probe.element_height;

dx_prb = array_length/2*cos(prb_theta_x);
dy_prb = array_length/2*sin(prb_theta_x);

% Range of the scene
range1 = prb_center - [dx_prb, dy_prb, 0] + PrbDir*scene_depth;
range2 = prb_center + [dx_prb, dy_prb, 0] + PrbDir*scene_depth;
range3 = prb_center - [dx_prb, dy_prb, 0];
range4 = prb_center + [dx_prb, dy_prb, 0];


%% Visualize the scene of probe
% figure % TODO
% pcolor3(permute(kgrid_shear.x, [2,1,3]), permute(kgrid_shear.y, [2,1,3]), permute(kgrid_shear.z, [2,1,3]), permute(medium_shear.sound_speed_shear, [2,1,3]), ...
%     'alpha', 0.05)
% hold on
% surf([range1(1), range2(1); range3(1), range4(1)], ...
%     [range1(2), range2(2); range3(2), range4(2)], ...
%     [range1(3), range2(3); range3(3), range4(3)], ...
%     'EdgeColor','none','FaceAlpha',0.5,'FaceColor','g')
% % hold on
% % quiver3(prb_center(1),prb_center(2),prb_center(3),N_scene(1),N_scene(2),N_scene(3))
% % hold on
% % quiver3(prb_center(1),prb_center(2),prb_center(3),PrbDir(1),PrbDir(2),PrbDir(3), 'Color', 'k')
% % 
% axis equal
% % xlim([minX,maxX])
% % ylim([minY,maxY])
% % zlim([minZ,maxZ])
% set(gca, 'ZDir','reverse')
% colormap cool

%% Define Transform matrix from probe coordinate to mesh coordinate

% Yaw rotation matrix
R_yaw = [cos(prb_theta_x), -sin(prb_theta_x), 0, 0;
         sin(prb_theta_x), cos(prb_theta_x), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];

% Roll rotation matrix
R_roll = [1, 0, 0, 0;
          0, cos(prb_theta_z), -sin(prb_theta_z), 0;
          0, sin(prb_theta_z), cos(prb_theta_z), 0;
          0,0,0,1];

% Translation matrix
T = [1,0,0,prb_center(1);
    0,1,0,prb_center(2);
    0,0,1,prb_center(3);
    0,0,0,1];

TF = T*R_yaw*R_roll;%*[scatterers, ones(size(xxp_speckle))]';

%% Reverse transform matrix from mesh to probe

T_re = [1,0,0,-prb_center(1);
    0,1,0,-prb_center(2);
    0,0,1,-prb_center(3);
    0,0,0,1];

R_roll_re = [1, 0, 0, 0;
          0, cos(-prb_theta_z), -sin(-prb_theta_z), 0;
          0, sin(-prb_theta_z), cos(-prb_theta_z), 0;
          0,0,0,1];

R_yaw_re = [cos(-prb_theta_x), -sin(-prb_theta_x), 0, 0;
         sin(-prb_theta_x), cos(-prb_theta_x), 0, 0;
         0, 0, 1, 0;
         0, 0, 0, 1];

TF_rev = R_roll_re*R_yaw_re*T_re;%*sca_mesh;



%% create a sensor mask covering the entire computational domain using the

% sensor mask in probe coord
[mask_coor.x, mask_coor.y, mask_coor.z] = meshgrid(-(array_length+2e-3)/2:dx_shear/15:(array_length+2e-3)/2, ...
    -array_height/2:dx_shear/15:array_height/2, ...
    0:dx_shear/15:scene_depth+2e-3);

% transform coord
coor_mesh = TF*[mask_coor.x(:), mask_coor.y(:), mask_coor.z(:), ones(size(mask_coor.x(:)))]';

% find the mesh grids that corresponding to the transformed coord
mask_idx = knnsearch([kgrid_shear.x(:), kgrid_shear.y(:), kgrid_shear.z(:)], coor_mesh(1:3,:)');

sensor_shear.mask = zeros(Nx_shear, Ny_shear, Nz_shear);
sensor_shear.mask(mask_idx) = 1;

sensor_shear.record = {'u'};


%% Shear wave simulation

display_mask = source_shear.s_mask;
input_args_shear = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false, 'DataCast', DATA_CAST};
% Run the simulation TODO
sensor_data_shear = pstdElastic3D(kgrid_shear, medium_shear, source_shear, sensor_shear, input_args_shear{:});

% sensor_data_shear.ux = gather(sensor_data_shear.ux); % TODO
% sensor_data_shear.uy = gather(sensor_data_shear.uy);
% sensor_data_shear.uz = gather(sensor_data_shear.uz);
%% displacement over time

% find the time that shear wave propergate through the entire medium
t_start = sqrt(xrange^2+yrange^2+zrange^2)/c_shear1;
% t_start = xrange/c_shear1;

% set bmode image starting time
t_current = t_start + 50 * kgrid_shear.dt;
t_next = t_current + 1/prf;


disp_x = [];
disp_y = [];
disp_z = [];

while t_next < kgrid_shear.t_array(end)

    est_disp_3D;

    disp_x = cat(2, disp_x, displace.x);
    disp_y = cat(2, disp_y, displace.y);
    disp_z = cat(2, disp_z, displace.z);

    % update pulse emit time
    t_current = t_next;
    t_next = t_current + 1/prf;

end

%% Generate ground truth

[X,Z] = meshgrid(-(array_length)/2:dx_shear:(array_length)/2, 5e-3:dx_shear:scene_depth);
gt_grid = TF*[X(:), zeros(size(X(:))), Z(:), ones(size(X(:)))]';

for i = 1:size(disp_z,2)

    disp_whole = zeros(Nx_shear, Ny_shear, Nz_shear);
    disp_whole(logical(sensor_shear.mask)) = disp_x(:,i);
    Dx_gt = interpn(kgrid_shear.x, kgrid_shear.y, kgrid_shear.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    disp_whole = zeros(Nx_shear, Ny_shear, Nz_shear);
    disp_whole(logical(sensor_shear.mask)) = disp_y(:,i);
    Dy_gt = interpn(kgrid_shear.x, kgrid_shear.y, kgrid_shear.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');

    disp_whole = zeros(Nx_shear, Ny_shear, Nz_shear);
    disp_whole(logical(sensor_shear.mask)) = disp_z(:,i);
    Dz_gt = interpn(kgrid_shear.x, kgrid_shear.y, kgrid_shear.z, disp_whole, ...
        gt_grid(1,:), gt_grid(2,:), gt_grid(3,:), 'spline');
    
    % deal with NaN
    Dx_gt(isnan(Dx_gt)) = 0;
    Dy_gt(isnan(Dy_gt)) = 0;
    Dz_gt(isnan(Dz_gt)) = 0;

    disp_gt(:,:,i) = reshape(Dx_gt.*-PrbDir(1) + Dy_gt.*-PrbDir(2) + Dz_gt.*-PrbDir(3), size(X));

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initial scatterers for ultrasound image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% resolution cell. Ref. Wagner, Robert F., et al. 
lateral_res = lambda/array_length;
axis_res = pulse_duration*lambda/2;

% % Num of scatterers
% lateral_sca_num = round(array_length/lateral_res*20);
% axis_sca_num = round(scene_depth/axis_res*3);
% [X_sca_prb,Z_sca_prb] = meshgrid(linspace(-array_length/2, array_length/2, lateral_sca_num), ...
%     linspace(3e-3, scene_depth, axis_sca_num));   % regular speckle, bad

sca_per_cell = 50;

xxp_speckle = [];
zzp_speckle = [];
for x = -array_length/2:lateral_res:array_length/2
    for z = 5e-3:axis_res:scene_depth
        xxp_speckle = [xxp_speckle; random('unif',x,x+lateral_res,sca_per_cell,1)];
        zzp_speckle = [zzp_speckle; random('unif',z,z+axis_res,sca_per_cell,1)];
    end
end


% Scatterers in probe coordinate
scatterers = [xxp_speckle, zeros(size(xxp_speckle)), zzp_speckle];
sca_mesh = TF*[scatterers, ones(size(xxp_speckle))]';



%% Define scatterers intensity

% Initial amplitude
amp = ones([length(scatterers),1]);

% Find scatterers belong to inclusions
incl=sqrt((sca_mesh(1,:)-cx1).^2+(sca_mesh(2,:)-(cy1-yrange/2)).^2+(sca_mesh(3,:)-(cz1-zrange/2)).^2)<cr1+4e-4 & ...
    sqrt((sca_mesh(1,:)-cx1).^2+(sca_mesh(2,:)-(cy1-yrange/2)).^2+(sca_mesh(3,:)-(cz1-zrange/2)).^2)>cr1-4e-4 | ...
    sqrt((sca_mesh(1,:)-cx2).^2+(sca_mesh(2,:)-(cy2-yrange/2)).^2+(sca_mesh(3,:)-(cz2-zrange/2)).^2)<cr2+4e-4 & ...
    sqrt((sca_mesh(1,:)-cx2).^2+(sca_mesh(2,:)-(cy2-yrange/2)).^2+(sca_mesh(3,:)-(cz2-zrange/2)).^2)>cr2-4e-4;


% Set higher amp for liver scatterers
amp(incl) = 2;


%% Scatterers moving

for i = 1:size(disp_z,2)

    disp_whole = zeros(Nx_shear, Ny_shear, Nz_shear);
    disp_whole(logical(sensor_shear.mask)) = disp_x(:,i);
    Dx_intrp = interpn(kgrid_shear.x, kgrid_shear.y, kgrid_shear.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    disp_whole = zeros(Nx_shear, Ny_shear, Nz_shear);
    disp_whole(logical(sensor_shear.mask)) = disp_y(:,i);
    Dy_intrp = interpn(kgrid_shear.x, kgrid_shear.y, kgrid_shear.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');

    disp_whole = zeros(Nx_shear, Ny_shear, Nz_shear);
    disp_whole(logical(sensor_shear.mask)) = disp_z(:,i);
    Dz_intrp = interpn(kgrid_shear.x, kgrid_shear.y, kgrid_shear.z, disp_whole, ...
        sca_mesh(1,:,i), sca_mesh(2,:,i), sca_mesh(3,:,i), 'spline');
    
    % deal with NaN
    Dx_intrp(isnan(Dx_intrp)) = 0;
    Dy_intrp(isnan(Dy_intrp)) = 0;
    Dz_intrp(isnan(Dz_intrp)) = 0;

    sca_mesh(:,:,i+1) = sca_mesh(:,:,i) + [Dx_intrp; Dy_intrp; Dz_intrp; zeros(1, length(scatterers))];

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%               Ultrasound simulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fs=100e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s] 

attenfreq = 0.55;  % Attenuation [dB/MHz/cm]

%% field II initialization

field_init(0);
set_field('c',c0);              % Speed of sound [m/s]
set_field('fs',fs);             % Sampling frequency [Hz]
set_field('use_rectangles',1);  % use rectangular elements

set_field('att', attenfreq*100*f0/1e6);  %% db/m
set_field('freq_att',attenfreq*100/1e6);  %% dB/m/Hz
% set_field('freq_att', 0);
set_field('att_f0', f0);
set_field('use_att',1);

%% Pulse definition

pulse = uff.pulse();
pulse.fractional_bandwidth = 0.65;        % probe bandwidth [1]
pulse.center_frequency = f0;
t0 = (-1/pulse.fractional_bandwidth/f0): dt : (1/pulse.fractional_bandwidth/f0);
impulse_response = gauspuls(t0, f0, pulse.fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2+1;   

% display the pulse to check that the lag estimation is on place 
% (and that the pulse is symmetric)

% figure;
% plot((0:(length(two_way_ir)-1))*dt -lag*dt,two_way_ir); hold on; grid on; axis tight
% plot((0:(length(two_way_ir)-1))*dt -lag*dt,abs(hilbert(two_way_ir)),'r')
% plot([0 0],[min(two_way_ir) max(two_way_ir)],'g');
% legend('2-ways pulse','Envelope','Estimated lag');
% title('2-ways impulse response Field II');

%% Aperture Objects

noSubAz=round(probe.element_width/(lambda/8));        % number of subelements in the azimuth direction
noSubEl=round(probe.element_height/(lambda/8));       % number of subelements in the elevation direction
Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 
Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]); 

% set the excitation, impulse response and baffle as below:
xdc_excitation (Th, excitation);
xdc_impulse (Th, impulse_response);
xdc_baffle(Th, 0);
xdc_center_focus(Th,[0 0 0]);
xdc_impulse (Rh, impulse_response);
xdc_baffle(Rh, 0);
xdc_center_focus(Rh,[0 0 0]);


%% Define plane wave sequence
% Define the start_angle and number of angles
F_number = 1.7;
F = size(sca_mesh,3);                                        % number of frames
if Na == 1
    alpha=0;%linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
else
    alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
end


%% output data
cropat=round(2*scene_depth/c0/dt);    % maximum time sample, samples after this will be dumped
CPW=zeros(cropat,probe.N,Na,F);  % impulse response channel data


%% Compute CPW signals
time_index=0;
disp('Field II: Computing CPW dataset');
parfor f=1:F
    disp(['Calculating frame ',num2str(f),' of ',num2str(F)]);

    sca_prb = TF_rev*sca_mesh(:,:,f); % reverse the scatterers back to probe coordinate
    sca_prb = sca_prb(1:3,:)';

    field_init(-1);

    set_field('c',c0);              % Speed of sound [m/s]
    set_field('fs',fs);             % Sampling frequency [Hz]
    set_field('use_rectangles',1);  % use rectangular elements
    
    set_field('att', attenfreq*100*f0/1e6);  %% db/m
    set_field('freq_att',attenfreq*100/1e6);  %% dB/m/Hz
    % set_field('freq_att', 0);
    set_field('att_f0', f0);
    set_field('use_att',1);


    Th = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
    Rh = xdc_linear_array (probe.N, probe.element_width, probe.element_height, kerf, noSubAz, noSubEl, [0 0 Inf]);
    xdc_excitation (Th, excitation);
    xdc_impulse (Th, impulse_response);
    xdc_baffle(Th, 0);
    xdc_center_focus(Th,[0 0 0]);
    xdc_impulse (Rh, impulse_response);
    xdc_baffle(Rh, 0);
    xdc_center_focus(Rh,[0 0 0]);
         

    for n=1:Na
        disp(['Calculating frame ',num2str(f),' of ',num2str(F), ', angle ',num2str(n),' of ',num2str(Na)]);
        % transmit aperture
        xdc_apodization(Th,0,ones(1,probe.N));
        xdc_times_focus(Th,0,probe.geometry(:,1)'.*sin(alpha(n))/c0);
        
        % receive aperture
        xdc_apodization(Rh, 0, ones(1,probe.N));
        xdc_focus_times(Rh, 0, zeros(1,probe.N));

        % do calculation
        [v,t]=calc_scat_multi(Th, Rh, sca_prb, amp);

        % absorb t
        t_vec = 0:dt:(size(v,1)-1)*dt;
        t_inq = -t:dt:-t+(cropat-1)*dt;
        v_inq = zeros(cropat,probe.N);
        for i = 1:probe.N
            v_inq(:,i) = interp1(t_vec, v(:,i), t_inq, 'linear', 0);
        end
         
        % build the dataset
        CPW(:,:,n,f)=single(v_inq);
        % if size(v,1)<cropat
        %     CPW(:,:,n)=padarray(v,[cropat-size(v,1) 0],0,'post');    
        % else
        %     CPW(:,:,n)=v(1:cropat,:);
        % end
         
        
    end
    
    field_end;
end

%% Save transmit sequence
for n = 1:Na
    seq(n)=uff.wave();
    seq(n).probe=probe;
    seq(n).source.azimuth=alpha(n);
    seq(n).source.distance=Inf;
    seq(n).sound_speed=c0;
    seq(n).delay = -lag*dt;
end

%% Channel Data
channel_data = uff.channel_data();
channel_data.sampling_frequency = fs;
channel_data.sound_speed = c0;
channel_data.initial_time = 0;
channel_data.pulse = pulse;
channel_data.probe = probe;
channel_data.sequence = seq;
channel_data.PRF = prf;
channel_data.data = CPW;

%%
save(SAVE_DIR, "channel_data", 'disp_gt', 'disp_x', 'disp_y', 'disp_z' , 'source_freq', 'prb_center', 'prb_theta_x', 'prb_theta_z', 'scene_depth', '-v7.3')
% TODO
%% Scan
%
% The scan area is defines as a collection of pixels spanning our region of 
% interest. For our example here, we use the *linear_scan* structure, 
% which is defined with two components: the lateral range and the 
% depth range. *scan* too has a useful *plot* method it can call.

sca=uff.linear_scan('x_axis',linspace(-array_length/2,array_length/2,512).', 'z_axis', linspace(5e-3,scene_depth,1024).');

%% Pipeline
%
% With *channel_data* and a *scan* we have all we need to produce an
% ultrasound image. We now use a USTB structure *pipeline*, that takes an
% *apodization* structure in addition to the *channel_data* and *scan*.

das=midprocess.das();
das.channel_data=channel_data;
das.scan=sca;

% das.dimension = dimension.receive();

das.transmit_apodization.window=uff.window.tukey25;
das.transmit_apodization.f_number=F_number;

das.receive_apodization.window=uff.window.tukey25;
das.receive_apodization.f_number=F_number;
das.code = "matlab";

% beamforming
b_data=das.go();

%% Display images
b_data.plot([], [], 40);
