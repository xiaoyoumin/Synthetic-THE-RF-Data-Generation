



%% ========================================================================
%
%   This script perform ultrasound time harmonic shear wave
%   elastography simulation using phantom geometry setup.
%   
%   Dependencies:
%       - k-Wave (http://www.k-wave.org/)
%       - Field II (https://field-ii.dk/)
%       - USTB (https://ultrasoundtoolbox.com/)
%
% =========================================================================

clear; clc; close all;

% Add dependencies
addpath('/fp/homes01/u01/ec-chaoranh/Documents/MATLAB/ustb')
addpath('/fp/homes01/u01/ec-chaoranh/Documents/MATLAB/field_II/')

addpath("FIELD_fun")
addpath("kwave_fun")
addpath("utils")
%% ---------------- Parameters for shear-wave simulation ------------------

% Define phantom geometry
shear_params.xrange = 100e-3;  % length [m]
shear_params.zrange = 95e-3;   % height [m]
shear_params.yrange = 120e-3;  % width  [m]

% Define stiffer inclusion 1
shear_params.cx1=0; shear_params.cy1=-20e-3; shear_params.cz1=15e-3;  % center coordinate
shear_params.cr1 = 5e-3;                   % radius

% Define stiffer inclusion 2
shear_params.cx2=0; shear_params.cy2=20e-3; shear_params.cz2=35e-3;  % center coordinate
shear_params.cr2 = 10e-3;                  % radius



% Define material physical properties
shear_params.c_shear_bkg = 2.32;   % background shear wave speed [m/s]
shear_params.c_shear_incl = 4.67;   % inclusion shear wave speed [m/s]
shear_params.source_freq = 200; % shaker frequency [Hz]
shear_params.rho0 = 1079;       % medium density [kg/m3]




%% -------- Parameters for ultrasound pulse-echo imaging simulation -------

% Define probe
c0                      = 1540;                   % Speed of sound compression wave [m/s]
f0                      = 5.1333e+06;             % Transducer center frequency [Hz]
lambda                  = c0/f0;                  % Wavelength [m]
prb.element_height      = 5e-3;                   % Height of element [m]
prb.pitch               = 0.300e-3;               % probe.pitch [m]
prb.kerf                = 0.03e-03;               % gap between elements [m]
prb.element_width       = prb.pitch-prb.kerf;     % Width of element [m]
lens_el                 = 20e-3;                  % position of the elevation focus
prb.N                   = 128;                    % Number of elements
pulse_duration          = 2.5;                    % pulse duration [cycles]
prf                     = 3000;                   % pulse repetition frequency [Hz]

fs=100e6;    % Sampling frequency [Hz]
dt=1/fs;     % Sampling step [s] 

attenfreq = 0.55;  % Attenuation [dB/MHz/cm]

% Scene of the image
array_length = (prb.N-1)*prb.pitch;
scene_depth = 50e-3;

%% Define k-Wave objects for phantom model

[kgrid, medium] = create_phantom(shear_params);

%% Define wave source k-Wave object

source = create_source(kgrid, shear_params.source_freq);

%% Define k-Wave sensor object

prb_center = [0 -20e-3 kgrid.z_vec(1)]; % location of the probe in the mesh
prb_theta = [deg2rad(45), deg2rad(0)]; % rotation of the probe

[TF, TF_rev] = coor_transformation(prb_center, prb_theta);
sensor = create_sensor(kgrid, array_length, prb.element_height, scene_depth, TF);

%% render

viewer = viewer3d(BackgroundColor="white",BackgroundGradient="off");

volshow(permute(medium.sound_speed_shear, [2,1,3]),Parent=viewer, Colormap=cool, Alphamap=[0.005,0.1]);
volshow(permute(source.s_mask, [2,1,3]),Parent=viewer, Colormap=hot, Alphamap=[0.00,0.2]);
volshow(permute(sensor.mask, [2,1,3]),Parent=viewer, Colormap=parula, Alphamap=[0.00,0.1]);


%% Shear wave simulation

display_mask = source.s_mask;
input_args_shear = {'DisplayMask', display_mask, 'PMLInside', false, 'PlotPML', false, 'DataCast', DATA_CAST};
% Run the simulation TODO
sensor_data_shear = pstdElastic3D(kgrid, medium, source, sensor, input_args_shear{:});


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

    displace = get_inter_trans_disp_kwave(sensor_data_shear, kgrid, t_current, t_next);

    disp_x = cat(2, disp_x, displace.x);
    disp_y = cat(2, disp_y, displace.y);
    disp_z = cat(2, disp_z, displace.z);

    % update pulse emit time
    t_current = t_next;
    t_next = t_current + 1/prf;

end


%% create initial scatterer map

sca_per_cell = 50;
scatterers = create_scataterers(prb, scene_depth, pulse_duration, lambda, sca_per_cell);

%% time-varying scatterer position

sca_mesh = time_vary_sca_kwave(kgrid, sensor, disp_x, disp_y, disp_z, scatterers, TF);

%% Define scatterers intensity

% Initial amplitude
amp = ones([size(scatterers,1),1]);

% Find scatterers belong to inclusions
incl=sqrt((sca_mesh(1,:,1)-shear_params.cx1).^2+(sca_mesh(2,:,1)-(shear_params.cy1)).^2+(sca_mesh(3,:,1)-(shear_params.cz1-shear_params.zrange/2)).^2)<shear_params.cr1+4e-4 & ...
    sqrt((sca_mesh(1,:,1)-shear_params.cx1).^2+(sca_mesh(2,:,1)-(shear_params.cy1)).^2+(sca_mesh(3,:,1)-(shear_params.cz1-shear_params.zrange/2)).^2)>shear_params.cr1-4e-4 | ...
    sqrt((sca_mesh(1,:,1)-shear_params.cx2).^2+(sca_mesh(2,:,1)-(shear_params.cy2)).^2+(sca_mesh(3,:,1)-(shear_params.cz2-shear_params.zrange/2)).^2)<shear_params.cr2+4e-4 & ...
    sqrt((sca_mesh(1,:,1)-shear_params.cx2).^2+(sca_mesh(2,:,1)-(shear_params.cy2)).^2+(sca_mesh(3,:,1)-(shear_params.cz2-shear_params.zrange/2)).^2)>shear_params.cr2-4e-4;


% Set higher amp for liver scatterers
amp(incl) = 2;


%% plane wave transmit sequence

alpha_max = deg2rad(15);%atan(1/2/F_number);
Na=1;                                      % number of plane waves 

if Na == 1
    alpha=0;
else
    alpha=linspace(-alpha_max,alpha_max,Na);    % vector of angles [rad]
end


%% ultrasound pulse setup

[impulse_response, excitation, lag] = pulse_setup(f0,pulse_duration,0.65,fs);


%% Compute CPW signals

F = size(sca_mesh,3);                                        % number of frames

cropat=round(2*scene_depth/c0/dt);    % maximum time sample, samples after this will be dumped

% output channel data [samples, elements, n transmits, n frames]
CPW=zeros(cropat,prb.N,Na,F);  % impulse response channel data

disp('Field II: Computing CPW dataset');
parfor f=1:F
    disp(['Calculating frame ',num2str(f),' of ',num2str(F)]);

    sca_prb = TF_rev*sca_mesh(:,:,f); % reverse the scatterers back to probe coordinate
    sca_prb = sca_prb(1:3,:)';

    % FIELD initialization
    [Th, Rh] = FIELD_setup(prb,fs,attenfreq,f0,c0,impulse_response, excitation);

    for n=1:Na
        disp(['Calculating frame ',num2str(f),' of ',num2str(F), ', angle ',num2str(n),' of ',num2str(Na)]);
        
        FIELD_TRx_aperture(prb,Th,Rh,c0,alpha(n));


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
