function scatterers = create_scataterers(prb, scene_depth, pulse_duration, lambda, sca_per_cell)



if isfield(prb,'radius')
    dtheta=2*asin(prb.pitch/2/prb.radius); 
    l_range = abs((prb.N-1)*dtheta);
    radius = prb.radius;
else
    l_range = (prb.N-1)*prb.pitch;
    radius = 0;
end


% resolution cell. Ref. Wagner, Robert F., et al. 
lateral_res = lambda/l_range;
axis_res = pulse_duration*lambda/2;



xxp = [];
zzp = [];
for x = -l_range/2:lateral_res:l_range/2
    for z = radius+1e-3:axis_res:radius+scene_depth
        xxp = [xxp; random('unif',x,x+lateral_res,sca_per_cell,1)];
        zzp = [zzp; random('unif',z,z+axis_res,sca_per_cell,1)];
    end
end


if isfield(prb,'radius')
   xxp_speckle = zzp.*sin(xxp);
   zzp_speckle = zzp.*cos(xxp)-prb.radius; % minus radius put the 
else
   xxp_speckle = xxp;
   zzp_speckle = zzp;

end


% Scatterers in probe coordinate
scatterers = [xxp_speckle, zeros(size(xxp_speckle)), zzp_speckle];
