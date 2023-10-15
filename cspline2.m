% causal cubic spline in Meinsma et al. formulation [1]
% input:
%   t --- time vector (column);
%   u --- input signal (column);
%   fs --- sampling rate [Hz];
%   l --- look-ahead.
% output:
%   pp --- piecewise cubic polynomial approximation of input signal.
% references:
% [1] G. Meinsma and L. Mirkin, "L2 Sampled signal reconstruction with 
%     causality constraints - Part I: Setup and solutions," in IEEE 
%     Transactions on Signal Processing, vol. 60, no. 5, pp. 2260-2272, 
%     May 2012, doi: 10.1109/TSP.2012.2185228.
function pp = cspline2(t,u,fs,l)
    alpha = -2 + sqrt(3);
    
    % causal filter
    b1 = [0 6*alpha];
    a1 = [1 -alpha];
    tmp = filter(b1,a1,u);
    w11 = (4-sqrt(3))*u + tmp;
    w12 = (3-sqrt(3))*u + tmp;
    w1 = [w11.';w12.'];

    % l-causal filter
    b2 = alpha.^(l-1:-1:1);
    tmp = filter(b2,1,u);
    tmp = zeroshift(tmp,-(l-1),1); % compensate FIR filter delay
    tmp = 6*fs^3*[(3*sqrt(3)-3)*tmp.' + (-4+3*sqrt(3))*u.';...
                  -sqrt(3)*tmp.'      + (1-sqrt(3))*u.'];
    tmp = zeroshift(tmp,-1,2);

    % coupling
    X = (fs^3)*[6*sqrt(3)-6, 3-3*sqrt(3); 3-3*sqrt(3), sqrt(3)];
    AX = (alpha^l*fs^3)*[6*sqrt(3), 3-3*sqrt(3); -3-3*sqrt(3), sqrt(3)];
    tmp = tmp - X*w1;
    w2 = tmp + zeroshift(AX*w1,-l,2);

    % output coefficients
    coefs = zeros(length(u),4);
    w2 = w2.';
    if l == 1
        w2(:,2) = 0; % analytical result for lookahead = 1 (strictly causal)
    end
    w1 = w1.';
    coefs(:,1) = -w2(:,1)/6;
    coefs(:,2) = (w2(:,1)+w2(:,2))/(2*fs);
    coefs(:,3) = (w2(:,1)+w2(:,2))/(2*sqrt(3)*fs^2) + w1(:,2)*fs;
    coefs(:,4) = w1(:,1)-w1(:,2);

    % create polynomial
    pp = mkpp(t, coefs(1:end-1,:));
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