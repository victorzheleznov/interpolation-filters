% shifted linear interpolation [1]
% input:
%   t --- time vector;
%   u --- signal samples;
%   fs --- sampling rate [Hz];
%   tau --- shift.
% output:
%   pp --- piecewise linear function;
%   c --- pre-filter output.
% references:
% [1] T. Blu, P. Thevenaz and M. Unser, "Linear interpolation revitalized," 
%     in IEEE Transactions on Image Processing, vol. 13, no. 5, pp. 710-719, 
%     May 2004, doi: 10.1109/TIP.2004.826093.
function [pp,c] = shlinear(t,u,fs,tau)
    % pre-filter
    b = 1/(1-tau);
    a = [1, tau/(1-tau)];
    c0 = u(1);
    zi = -tau/(1-tau)*c0;
    c = filter(b,a,u(2:end),zi);
    c = [c0;c];

    % calculate coeffitients
    coefs = zeros(length(c),2);
    coefs(:,1) = (zeroshift(c,-1,1)-c)*fs;
    coefs(:,2) = c;

    % create piecewise polynomial
    breaks = t + tau/fs;
    pp = mkpp(breaks, coefs(1:end-1,:));
end