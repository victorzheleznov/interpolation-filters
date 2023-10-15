% shifted linear interpolation [1]
% input:
%   t --- time vector (column);
%   u --- input signal (column);
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

%% FUNCTIONS
% shift 2d array values and fill boundary with zeros
% input:
%   in --- input 2d array;
%   l --- integer shift (positive shifts toward the end and negative shifts toward the beginning);
%   d --- dimension (d = 1 to shift rows, d = 2 to shift columns).
% output:
%   out --- shifted array.
function out = zeroshift(in,l,d)
    out = circshift(in,l,d);
    if l > 0
        if d == 1
            out(1:l,:) = 0;
        elseif d == 2
            out(:,1:l) = 0;
        end
    elseif l < 0
        if d == 1
            out(end+l+1:end,:) = 0;
        elseif d == 2
            out(:,end+l+1:end) = 0;
        end
    end
end