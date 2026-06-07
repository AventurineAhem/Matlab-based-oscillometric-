clc; clear; close all;
tic

filename = 'sample_07_06.dat';
data = importdata(filename);

time = data(:,1);
signal1 = 950*(data(:,2)-0.712);   % oscillometric signal
signal2 = data(:,3);               % cuff pressure (mmHg)

target = 130;

buff = abs(signal1 - target);
[minVal, indx] = min(buff);

if minVal > 5
    [~, indx] = max(signal1);
    fprintf("Target not found. Using MAX value.\n");
end

cutoutsignal = signal1(indx:end);
time = time(indx:end);
signal2 = signal2(indx:end);

buff = abs(cutoutsignal - 46);
[~, indx2] = min(buff);

cutoutsignal = cutoutsignal(1:indx2);
time = time(1:indx2);
signal2 = signal2(1:indx2);

%% Filtering
Fs = 1000;
Fc_highpass = 0.5;
Fc_lowpass  = 20;
order = 4;

Wn = [Fc_highpass Fc_lowpass]/(Fs/2);
[b,a] = butter(order, Wn, 'bandpass');

filtered_signal = filtfilt(b, a, cutoutsignal);

figure(1);
plot(time, cutoutsignal,'g','LineWidth',1.5);
grid on; title('Original Signal');

figure(2);
plot(time, filtered_signal,'b','LineWidth',1.5);
grid on; title('Filtered Signal');

%% Minima & Maxima
inv_sig = -filtered_signal;

is_val = [false; (inv_sig(2:end-1) > inv_sig(1:end-2)) & ...
                 (inv_sig(2:end-1) >= inv_sig(3:end)); false];

locs_vl = find(is_val);
pks_vl  = filtered_signal(locs_vl);

thrs = 600;

i=1;
while i<length(locs_vl)
    if (locs_vl(i+1)-locs_vl(i))<thrs
        if pks_vl(i)>pks_vl(i+1)
            locs_vl(i)=[]; pks_vl(i)=[];
        else
            locs_vl(i+1)=[]; pks_vl(i+1)=[];
        end
        i=max(i-1,1);
    else
        i=i+1;
    end
end

locs_pk=[]; pks_pk=[];

for k=1:length(locs_vl)-1
    s1=locs_vl(k);
    s2=locs_vl(k+1);

    seg=filtered_signal(s1:s2);
    [mx,idx]=max(seg);

    locs_pk(end+1)=s1+idx-1;
    pks_pk(end+1)=mx;
end

figure(3);
plot(time, filtered_signal,'b'); hold on;
plot(time(locs_pk),pks_pk,'r+');
plot(time(locs_vl),pks_vl,'w+');
grid on; title('Maxima & Minima');

%% Function through minimas f(x)
f = interp1(time(locs_vl), pks_vl, time, 'spline','extrap');

figure(4);
plot(time, filtered_signal,'b'); hold on;
plot(time, f,'m','LineWidth',2);
grid on; title('F(x)');

%%Area per cycle
area_per_cycle=zeros(1,length(locs_vl)-1);
time_mid=zeros(1,length(locs_vl)-1);

for i=1:length(locs_vl)-1
    idx=locs_vl(i):locs_vl(i+1);

    segment = filtered_signal(idx) - f(idx);
    segment(segment<0)=0;

    area_per_cycle(i)=trapz(time(idx),segment);
    time_mid(i)=mean(time(idx));
end

figure(5);
plot(time_mid,area_per_cycle,'g-o');
grid on; title('Area per Cycle');

time_dense = linspace(min(time_mid),max(time_mid),500);
area_interp = interp1(time_mid,area_per_cycle,time_dense,'linear');
area_interp = smoothdata(area_interp,'gaussian',25);

%%Gaussian fit
ft = fittype('a*exp(-((x-b)^2)/(2*c^2))');

a0=max(area_interp);
b0=time_dense(area_interp==max(area_interp));
c0=10;

curve = fit(time_dense',area_interp',ft,...
    'StartPoint',[a0,b0(1),c0]);

figure(6);
plot(time_dense,area_interp,'g','LineWidth',2); hold on;
plot(time_dense,curve(time_dense),'r--','LineWidth',2);
grid on; title('Gaussian Fit');

gauss_vals = curve(time_dense);

% MAP 
[MA, idx_max] = max(gauss_vals);
t_MAP = time_dense(idx_max);

% Ratio levels
SBPA_level = 0.5 * MA;
DBPA_level = 0.8 * MA;

% SBP
left_side = 1:idx_max;
[~, idx_sbpa_local] = min(abs(gauss_vals(left_side) - SBPA_level));
idx_sbpa = left_side(idx_sbpa_local);
t_SBPA = time_dense(idx_sbpa);

% DBP 
right_side = idx_max:length(gauss_vals);
[~, idx_dbpa_local] = min(abs(gauss_vals(right_side) - DBPA_level));
idx_dbpa = right_side(idx_dbpa_local);
t_DBPA = time_dense(idx_dbpa);

plot(t_MAP, MA, 'ks','MarkerSize',10,'LineWidth',2);
plot(t_SBPA, SBPA_level, 'ro','MarkerSize',10,'LineWidth',2);
plot(t_DBPA, DBPA_level, 'bo','MarkerSize',10,'LineWidth',2);

legend('Envelope','Gaussian','MAP','SBPA','DBPA');

%% FINAL OUTPUT
SBP = interp1(time, cutoutsignal, t_SBPA);
MAP_val = interp1(time, cutoutsignal, t_MAP);
DBP = interp1(time, cutoutsignal, t_DBPA);

fprintf('SBP = %.2f mmHg\n', SBP);
fprintf('MAP = %.2f mmHg\n', MAP_val);
fprintf('DBP = %.2f mmHg\n', DBP);


toc
