% MOMS interpolation [1,2]
% input:
%   t --- time vector (column);
%   u --- input signal (column);
%   fs --- sampling rate [Hz].
% output:
%   pp --- piecewise cubic polynomial approximation of input signal.
% references:
% [1] T. Blu, P. Thévenaz and M. Unser, "High-quality causal interpolation 
%     for online unidimensional signal processing," 2004 12th European 
%     Signal Processing Conference, Vienna, Austria, 2004, pp. 1417-1420.
% [2] Petrinovic, D.. (2009). Continuous time domain properties of causal 
%     cubic splines. Signal Processing. 89. 1941-1958. 
%     10.1016/j.sigpro.2009.03.031.
function pp = moms(t,u,fs)
    % pre-filter
    b = 3;
    a = [2 1];
    w = filter(b,a,u);

    % output filter
    b3 = [1 -3 3 -1]*fs^3/6;
    b2 = [-1 6 -9 4]*fs^2/6;
    b1 = [0 1 4 -5] *fs/6;
    b0 = [0 0 4 2]/6;

    coefs(:,1) = filter(b3,1,w);
    coefs(:,2) = filter(b2,1,w);
    coefs(:,3) = filter(b1,1,w);
    coefs(:,4) = filter(b0,1,w);

    % compensate delay and remove first segment
    % (output FIR filters require one past sample => n = 1 is undefied)
    del = 2;
    coefs = zeroshift(coefs,-del,1);
    coefs(1,:) = 0;
    coefs = coefs(1:end-1,:);
    breaks = t;

    % make piecewise polynomial
    pp = mkpp(breaks, coefs);
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