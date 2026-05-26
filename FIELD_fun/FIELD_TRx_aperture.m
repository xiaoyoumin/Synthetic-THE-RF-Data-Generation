function FIELD_TRx_aperture(prb,Th,Rh,c0,alpha)

% transmit aperture
xdc_apodization(Th,0,ones(1,prb.N));

if isfield(prb,'radius')
    focal = [0, 0 -prb.radius];
    xdc_center_focus (Th, [0 0 0]);
    xdc_focus (Th, 0, focal);
else
    x0=(1:prb.N)*prb.pitch;
    x0=x0-mean(x0);
    xdc_times_focus(Th,0,x0.*sin(alpha)/c0);
end

% receive aperture
xdc_apodization(Rh, 0, ones(1,prb.N));
xdc_focus_times(Rh, 0, zeros(1,prb.N));
