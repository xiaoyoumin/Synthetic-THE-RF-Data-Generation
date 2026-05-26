function [impulse_response, excitation, lag] = pulse_setup(f0,pulse_duration,fb,fs)

dt = 1/fs;

fractional_bandwidth = fb;        % probe bandwidth [1]
t0 = (-1/fractional_bandwidth/f0): dt : (1/fractional_bandwidth/f0);
impulse_response = gauspuls(t0, f0, fractional_bandwidth);
impulse_response = impulse_response-mean(impulse_response); % To get rid of DC

te = (-pulse_duration/2/f0): dt : (pulse_duration/2/f0);
excitation = square(2*pi*f0*te+pi/2);
one_way_ir = conv(impulse_response,excitation);
two_way_ir = conv(one_way_ir,impulse_response);
lag = length(two_way_ir)/2+1;   
