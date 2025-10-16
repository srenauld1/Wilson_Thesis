function data = fit_motion_binary_model(fly_data)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
dff = fly_data.dat.ts(1,:); % dat.ts(1,:)(:)
dff = dff(:);
moving_binary = fly_data.daq.motion.moving_or_not(:);

X = [moving_binary ones(size(moving_binary))];
beta = X \ dff;
k = beta(1);
c = beta(2);

% Predicted values for plotting or assessment:
dff_fit = X*beta;

tbl = table(moving_binary, dff);
mdl = fitlm(tbl, 'dff ~ moving_binary');  % Intercept is fitted by default
% Get parameters:
k_lm = mdl.Coefficients.Estimate(2);
c_lm = mdl.Coefficients.Estimate(1);

% Example plot
plot(dff, 'k'); hold on;
plot(dff_fit, 'r');
legend('DFF data', 'Fit');
xlabel('Sample');
ylabel('Signal');
title('Fit: dff = k * movenotmove + c');

% subtract out motion signal if k > 0
if k > 0
    dff_motion_removed = dff - dff_fit;
end


end