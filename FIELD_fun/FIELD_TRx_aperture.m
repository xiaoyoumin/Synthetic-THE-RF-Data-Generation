function FIELD_TRx_aperture(prb, Th, Rh, c0, tx_param)
%FIELD_TRX_APERTURE Apply transmit and receive aperture focusing in Field II.
%
% This helper configures one transmit event and the corresponding receive
% aperture for either the linear-probe phantom case or the curvilinear-probe
% liver case. Linear probes use a plane-wave steering angle, while
% curvilinear probes use a virtual source/focal point for diverging-wave
% transmission.
%
% Inputs:
%   prb      : Probe structure. Linear probes require N and pitch;
%              curvilinear probes also include radius.
%   Th       : Field II transmit aperture handle.
%   Rh       : Field II receive aperture handle.
%   c0       : Speed of sound used by Field II [m/s].
%   tx_param : Transmit parameter. For a linear probe, this must be a
%              scalar plane-wave angle [rad]. For a curvilinear probe, this
%              must be a 1-by-3 focal point [x, y, z] [m].
%
% Outputs:
%   No output arguments. The Field II aperture handles are updated in place.

xdc_apodization(Th, 0, ones(1, prb.N));

if isfield(prb, 'radius')
    if ~isnumeric(tx_param) || ~isequal(size(tx_param), [1, 3])
        error('For curvilinear probes, tx_param must be a 1-by-3 focal point [x, y, z].')
    end
    focal = tx_param;
    xdc_center_focus(Th, [0 0 0]);
    xdc_focus(Th, 0, focal);
else
    if ~isnumeric(tx_param) || ~isscalar(tx_param)
        error('For linear probes, tx_param must be a scalar transmit angle.')
    end
    alpha = tx_param;
    x0 = (1:prb.N)*prb.pitch;
    x0 = x0-mean(x0);
    xdc_times_focus(Th, 0, x0.*sin(alpha)/c0);
end

xdc_apodization(Rh, 0, ones(1, prb.N));
xdc_focus_times(Rh, 0, zeros(1, prb.N));
