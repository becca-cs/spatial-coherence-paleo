%% generate a climate signal :)
clear;
N = 10000; dt = 1;
t = -N:dt:0; 
idx = 5:10:N;
timepoints = t(idx);
ntp = N+1;
s = (1/N:1/N:ntp/2/N)'; % frequency vector (1/kyr)
smax = s(end);

flg = 0; % 0 for power, 1 for red
sigbdot = 1;

if flg==0
    nu = 0.25; % slope
    Pmax = 2*dt*sigbdot^2; smax = s(end);
    PSD = Pmax*(s/smax).^-nu;

% alternatively, try
elseif flg==1
    r = 0.9; 
    R = -log(r); % decorrelation rate
    omega = 2*pi*s./smax;
    PSD = sigbdot.*R./(R.^2+omega.^2);
end

temp = generate_timeseries(PSD,N./dt); % timeseries

%% generate a sedimentation record :)

sed_mean = 50; % this is in cm per kyr! need to transform into yr per cm
sed_std = 100;

[depth_sample,depth_interp,age_model,age_true,sed_true] = age_depth_model(timepoints,1/(sed_mean/1000),1/(sed_std/1000));

sed = interp1(age_true,sed_true,timepoints); 
timepoints_est = interp1(age_true,age_model,timepoints); % modeled ages!
timepoints_est = timepoints_est-min(timepoints_est); % go forward

sed_acc_rate = 1./sed*1000; % now in cm per kyr

% filter out NaN values
idx = ~isnan(sed_acc_rate);

sed_acc_rate = sed_acc_rate(idx);
timepoints = timepoints(idx);
timepoints_est = timepoints_est(idx);

%% add a seasonal cycle!
clim_signal = repmat(temp,[1 12]);

% add a seasonal cycle
mo = 0:11;
seasonal = 2*sin(mo*pi/6);

clim_signal = clim_signal + repmat(seasonal,[N+1 1]);

figure('Position',[10 10 1800 600])

for calib = 1:2
    if calib == 1
        calibration_type = 'Uk37';
        n_samples = Inf;
        meas_bias = 0.1;
    else
        calibration_type = 'MgCa';
        n_samples = 30;
        meas_bias = 0.1;
    end

    sigma_meas = 0.23; sigma_ind = 0.1; 

    % make time go forward
    time_true = age_true-min(age_true);
    time_model = age_model-min(age_model);
    timepoints_forward = timepoints-min(timepoints);
    time_forward = t-min(t);
    sed_acc_rate = sed_acc_rate;

    if calib == 2
        [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
            proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
            proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,time_forward(:),timepoints_forward(:),calibration_type,n_samples,...
            sigma_meas, sigma_ind, meas_bias,sed_acc_rate(:),@ForamGrowthfT,10,'ruber');

    else
        [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
            proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
            proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,time_forward(:),timepoints_forward(:),calibration_type,n_samples,...
            sigma_meas, sigma_ind, meas_bias,sed_acc_rate(:));
    end

    time_out_err = timepoints_est(valid_inds);

    subplot(1,3,1)
    if calib==1
        plot(time_out,clim_signal_ann,'k'); hold on;
        plot(time_out,reconstructed_climate,'r','linewidth',1)
    else
        plot(time_out,reconstructed_climate,'color',rgb('sky blue'),'linewidth',1)
        legend('climate','reconstructed climate (Uk37)','reconstructed climate (MgCa)');
        grid on; set(gca,'fontsize',12);
        ylabel('temperature (^oC)'); xlabel('time (yr)')
    end

    if calib==1
        subplot(1,3,2); title('Uk37')
    else
        subplot(1,3,3); title('MgCa')
    end
    plot(time_forward(1:10:end),mean(proxy_clim_signal(1:10:end,:),2),'k','linewidth',1); hold on;
    plot(time_out,proxy_bt,'color',rgb('pale green'),'linewidth',1);
    plot(time_out,proxy_bt_sb,'color',rgb('light purple'),'linewidth',1);

    if strcmp(calibration_type, ...
            'Uk37')
        plot(time_out,proxy_bt_sb_inf_b,'color',rgb('dark yellow'),'linewidth',1);
        plot(time_out,proxy_bt_sb_inf_b_n,'color',rgb('hazel'),'linewidth',1);
        plot(time_out_err,proxy_bt_sb_inf_b_n,'r','linewidth',1);
        legend('original signal','bioturbation','bioturbation + habitat bias', ...
            'bioturbation + habitat bias + calibration bias', 'bioturbation + habitat bias + calibration bias + noise',...
            'bioturbation + habitat bias + calibration bias + noise + time err','location','southwest');
    else
        plot(time_out,proxy_bt_sb_sampYM,'color',rgb('dark yellow'),'linewidth',1);
        plot(time_out,proxy_bt_sb_sampYM_b,'color',rgb('grey'),'linewidth',1);
        plot(time_out,proxy_bt_sb_sampYM_b_n,'color',rgb('magenta'),'linewidth',1);
        plot(time_out_err,proxy_bt_sb_sampYM_b_n,'color',rgb('sky blue'),'linewidth',1);
        legend('original signal','bioturbation','bioturbation + habitat bias', ...
            'bioturbation + habitat bias + aliasing', ...
            'bioturbation + habitat bias + aliasing + calibration bias', ...
            'bioturbation + habitat bias + aliasing + clibration bias + noise',...
            'bioturbation + habitat bias + aliasing + clibration bias + noise + time err', ...
            'location','southwest');
    end

    if calib==1
        title('Uk37');
    else
        title('MgCa')
    end

    grid on; set(gca,'fontsize',12);
    ylabel('proxy units'); xlabel('time (yr)')
end

%% okay, let's try this like 1000 times

% generate a climate signal :)
clear; 

for nruns = 1:1000

    N = 10000; dt = 1;
    t = -N:dt:0;
    idx = 5:10:N;
    timepoints = t(idx);
    ntp = N+1;
    s = (1/N:1/N:ntp/2/N)'; % frequency vector (1/kyr)
    smax = s(end);

    flg = 0; % 0 for power, 1 for red
    sigbdot = 1;

    if flg==0
        nu = 0.25; % slope
        Pmax = 2*dt*sigbdot^2; smax = s(end);
        PSD = Pmax*(s/smax).^-nu;

        % alternatively, try
    elseif flg==1
        r = 0.9;
        R = -log(r); % decorrelation rate
        omega = 2*pi*s./smax;
        PSD = sigbdot.*R./(R.^2+omega.^2);
    end

    temp = generate_timeseries(PSD,N./dt); % timeseries

    % generate a sedimentation record :)

    sed_mean = 50; % this is in cm per kyr! need to transform into yr per cm
    sed_std = 100;

    [depth_sample,depth_interp,age_model,age_true,sed_true] = age_depth_model(timepoints,1/(sed_mean/1000),1/(sed_std/1000));

    sed = interp1(age_true,sed_true,timepoints);
    timepoints_est = interp1(age_true,age_model,timepoints); % modeled ages!
    timepoints_est = timepoints_est-min(timepoints_est); % go forward

    sed_acc_rate = 1./sed*1000; % now in cm per kyr

    % filter out NaN values
    idx = ~isnan(sed_acc_rate);

    sed_acc_rate = sed_acc_rate(idx);
    timepoints = timepoints(idx);
    timepoints_est = timepoints_est(idx);

    % add a seasonal cycle!
    clim_signal = repmat(temp,[1 12]);

    % add a seasonal cycle
    mo = 0:11;
    seasonal = 2*sin(mo*pi/6);

    clim_signal = clim_signal + repmat(seasonal,[N+1 1]);

    for calib = 1:2
        if calib == 1
            calibration_type = 'Uk37';
            n_samples = Inf;
            meas_bias = 0.1;
        else
            calibration_type = 'MgCa';
            n_samples = 30;
            meas_bias = 0.1;
        end

        sigma_meas = 0.23; sigma_ind = 0.1;

        % make time go forward
        time_true = age_true-min(age_true);
        time_model = age_model-min(age_model);
        timepoints_forward = timepoints-min(timepoints);
        time_forward = t-min(t);
        sed_acc_rate = sed_acc_rate;

        if calib == 2
            [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
                proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
                proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,time_forward(:),timepoints_forward(:),calibration_type,n_samples,...
                sigma_meas, sigma_ind, meas_bias,sed_acc_rate(:),@ForamGrowthfT,10,'ruber');

        else
            [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
                proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
                proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,time_forward(:),timepoints_forward(:),calibration_type,n_samples,...
                sigma_meas, sigma_ind, meas_bias,sed_acc_rate(:));
        end

        time_out_err = timepoints_est(valid_inds);


        % just save original & reconstructed climate
        if calib==1
            reconstructed_Uk37{nruns} = reconstructed_climate;
            reconstructed_time_Uk37{nruns} = time_out_err;
        else
            reconstructed_MgCa{nruns} = reconstructed_climate;
            reconstructed_time_MgCa{nruns} = time_out_err;
        end
    end
end

%% Okay, lets plot power spectra & coherence WOO

% create big matrices
s_Uk37 = NaN*ones(500,1000); 
P_Uk37 = NaN*ones(500,1000);
s_MgCa = NaN*ones(500,1000);
P_MgCa = NaN*ones(500,1000);

for nruns = 1:1000
    [P1,s1] = pmtmLS(detrend(reconstructed_Uk37{nruns}),reconstructed_time_Uk37{nruns},3,0);
    [P2, s2] = pmtmLS(detrend(reconstructed_MgCa{nruns}),reconstructed_time_MgCa{nruns},3,0);

    % save stuff
    s_Uk37(1:length(s1),nruns) = s1;
    P_Uk37(1:length(s1),nruns) = P1;
    s_MgCa(1:length(s2),nruns) = s2;
    P_MgCa(1:length(s2),nruns) = P2;
end

%%
s_interp = 10.^linspace(-5,1,500); %s = repmat(s,[1000 1]);
%s = s';

% there's definitely a faster way to do this but i'm feeling lazy

P_Uk37_interp = NaN*ones(500,1000); P_MgCa_interp = NaN*ones(500,1000);

for n = 1:1000
    idx1 = ~isnan(s_Uk37(:,n)); idx2 = ~isnan(s_MgCa(:,n));
    P1 = interp1(s_Uk37(idx1,n),P_Uk37(idx1,n),s_interp);
    P_Uk37_interp(1:length(P1),n) = P1;
    P2 = interp1(s_MgCa(idx2,n),P_MgCa(idx2,n),s_interp);
    P_MgCa_interp(1:length(P2),n) = P2;
end

p = [0.025 0.5 0.975];

Q_Uk37 = quantile(P_Uk37_interp',p);
Q_MgCa = quantile(P_MgCa_interp',p);

idx = ~isnan(Q_Uk37(1,:)) & ~isnan(Q_MgCa(1,:));

shade(s_interp(idx),Q_Uk37(3,idx),s_interp(idx),Q_Uk37(1,idx),'filltype',[1 2],'fillalpha',0.2, 'fillcolor',...
    'r','color','none'); hold on;
h(1) = plot(s_interp(idx),Q_Uk37(2,idx),'r','linewidth',2);
shade(s_interp(idx),Q_MgCa(3,idx),s_interp(idx),Q_MgCa(1,idx),'filltype',[1 2],'fillalpha',0.2, 'fillcolor',...
    rgb('sky blue'),'color','none'); 
h(2) = plot(s_interp(idx),Q_MgCa(2,idx),'color',rgb('sky blue'),'linewidth',2);
h(3) = plot(s,PSD,'k','linewidth',2);
set(gca,'yscale','log','xscale','log');
set(gca,'fontsize',14); legend(h,'Uk37','MgCa','original climate signal')
%shade(s,Q_MgCa(1,:),s,Q_MgCa(5,:));