% causal cubic spline based on truncation of anticausal IIR filter [1,2]
% (cascade pre-filter and De Boor's output filter)
% input:
%   t --- time vector;
%   u --- signal samples;
%   fs --- sampling rate [Hz];
%   M --- pre-filter truncation.
% output:
%   pp --- piecewise cubic polynomial approximation of input signal.
% references:
% [1] Petrinovic, Davor. (2008). Causal Cubic Splines: Formulations, 
%     Interpolation Properties and Implementations. Signal Processing, 
%     IEEE Transactions on. 56. 5442 - 5453. 10.1109/TSP.2008.929133. 
% [2] Petrinovic, D.. (2009). Continuous time domain properties of causal 
%     cubic splines. Signal Processing. 89. 1941-1958. 
%     10.1016/j.sigpro.2009.03.031. 
function pp = cspline(t,u,fs,M)
    % filters parameters
    alpha = -2 + sqrt(3);
    s = 1;%1/(1-alpha^(M+1));

    % direct casual B-spline cascade filter coefficients
    b = 1;
    a = [1 -alpha];
    b_fir = alpha.^(M:-1:1).';

    %% cascade casual B-spline filter + De Boor's output filter
    % FIR filter
    u_fir = filter(b_fir,1,u);
    u_zM = circshift(u,M);
    u_zM(1:M) = 0;

    % IIR filter
    w_tmp = -alpha*s*(u_fir+u_zM);
    w = filter(b,a,w_tmp);
    
    % output filter (De Boor's formulation)
    b31 = [3 3 -3 -3]*fs^3;
    b32 = [0 -2 2]   *fs^3;
    b21 = [-3 -6 3 6]*fs^2;
    b22 = [0 3 -3]   *fs^2;
    b11 = [0 3 0 -3] *fs;

    coefs = zeros(length(u), 4);
    coefs(:,1) = filter(b31,1,w) + filter(b32,1,u_zM);
    coefs(:,2) = filter(b21,1,w) + filter(b22,1,u_zM);
    coefs(:,3) = filter(b11,1,w);
    coefs(:,4) = circshift(u_zM, 2);
    coefs(1:2,4) = 0;

    % compensate delay and remove first segment
    % (output FIR filters require one past sample => n = 1 is undefied)
    del = M+2;
    coefs = zeroshift(coefs,-del,1);
    coefs(1,:) = 0;
    coefs = coefs(1:end-1,:);
    breaks = t;

    % create piecewise polynomial
    pp = mkpp(breaks, coefs);
end