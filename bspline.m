% noncausal B-spline interpolation [1,2,3,4]
% input:
%   t --- time vector (column);
%   u --- input signal (column);
%   fs --- sampling rate [Hz].
% output:
%   pp --- piecewise cubic polynomial approximation of input signal.
% references:
% [1] Unser, Michael & Aldroubi, Akram & Eden, Murray. (1991). Fast 
%     B-Spline Transforms for Continuous Image Representation and 
%     Interpolation. IEEE Trans. Pattern Anal. Mach. Intell.. 13. 277-285. 
%     10.1109/34.75515.
% [2] M. Unser, A. Aldroubi and M. Eden, "B-spline signal processing. I. 
%     Theory," in IEEE Transactions on Signal Processing, vol. 41, no. 2, 
%     pp. 821-833, Feb. 1993, doi: 10.1109/78.193220.
% [3] M. Unser, A. Aldroubi and M. Eden, "B-spline signal processing. II. 
%     Efficiency design and applications," in IEEE Transactions on Signal 
%     Processing, vol. 41, no. 2, pp. 834-848, Feb. 1993, 
%     doi: 10.1109/78.193221.
% [4] M. Unser, "Splines: a perfect fit for signal and image processing," 
%     in IEEE Signal Processing Magazine, vol. 16, no. 6, pp. 22-38, Nov. 
%     1999, doi: 10.1109/79.799930.
function pp = bspline(t,u,fs)
    % pole of the direct B-spline filter
    alpha = -2 + sqrt(3);
    
    % direct B-spline cascade filters coefficients
    b1 = 1;
    a1 = [1 -alpha];
    b2 = -alpha;
    a2 = [1 -alpha];

    % forward IIR
    k0 = ceil(log(eps) / log(abs(alpha)));
    wp0 = u(k0+1);
    for i = k0:-1:1
        wp0 = u(i) + alpha*wp0;
    end
    zi = -a1(2)*wp0; % filter delay
    wp = filter(b1,a1,u(2:end),zi);
    wp = [wp0;wp];

    % backward IIR
    w0 = (-alpha/(1-alpha^2))*(2*wp(end)-u(end));
    zi = -a2(2)*w0; % filter delay
    wp = flip(wp);  % backward time
    w = filter(b2,a2,wp(2:end),zi);
    w = [w0;w];
    w = flip(w);

    % output filter
    b3 = [1,-3,3,-1]*fs^3;
    b2 = [0,3,-6,3] *fs^2;
    b1 = [0,3,0,-3] *fs;

    coefs(:,1) = filter(b3,1,w);
    coefs(:,2) = filter(b2,1,w);
    coefs(:,3) = filter(b1,1,w);
    coefs(:,4) = zeroshift(u,2,1);

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