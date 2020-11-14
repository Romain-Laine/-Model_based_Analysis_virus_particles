function [ hFigure ] = Review_loc_file( LocInfo )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

% profile on

X = LocInfo (:,1); % in nm
Y = LocInfo (:,2); % in nm
Frame = LocInfo (:,3);
ADC = LocInfo (:,4);

n_frame = max(Frame);
n_loc = max(size(X));    % total number of localizations

hFigure = figure('Color','white','Units','normalized','position',[0.2 0.2 0.55 0.55]);
subplot(2,2,1)
hist(Frame,n_frame)
xlabel 'Frame number'
ylabel 'Number of localization / frame'
axis tight
title(['Total number of localizations: ', num2str(n_loc)])

n_bins = ceil(n_loc/2);
[n,x] = hist(ADC, n_bins);
N = cumsum(n)/sum(n);

% max_display_ADC = interp1(N(N<0.99),x(N<0.99),0.95);
% ADC_10p = interp1(N(N<0.4),x(N<0.4),0.1);
% ADC_90p = interp1(N(N<0.99),x(N<0.99),0.9);

ADC_sorted = sort(ADC);
max_display_ADC = ADC_sorted(round(n_loc*0.95));    % n_loc * 0.95 gives the 95% top value
ADC_10p = ADC_sorted(round(n_loc*0.1));
ADC_90p = ADC_sorted(round(n_loc*0.9));
ADC_99p = ADC_sorted(round(n_loc*0.99));

subplot(2,2,4)
plot(x,100*N)
xlabel 'ADC'
ylabel 'Percentage (%)'
grid on
xlim([0 max_display_ADC])
title(['ADC 10%: ', num2str(ADC_10p,'%6.0f'),' - ADC 90%: ',num2str(ADC_90p,'%6.0f')])


subplot(2,2,2)
hist(ADC,n_bins)
xlabel 'ADC'
ylabel 'Occurences'
axis tight
xlim([0 max_display_ADC])
title(['Min ADC: ',num2str(min(ADC)),' & Max ADC: ',num2str(max(ADC))])


Frame_unique = unique(Frame);
ADC_av = zeros(max(size(Frame_unique)),1);

for i = 1: max(size(Frame_unique))
    ADC_av(i) = mean(ADC(Frame == Frame_unique(i)));
end

ADC_av = smooth(ADC_av,100);

subplot(2,2,3)
plot(Frame, ADC,'.')
xlabel 'Frame number'
ylabel 'Localization amplitude (ADC)'
axis tight
ylim([0 ADC_99p])
title 'ADC as a function of time'
hold on
plot(Frame_unique, ADC_av,'r-')

% profile viewer

end

