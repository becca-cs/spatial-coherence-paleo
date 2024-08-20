% run PSM in MATLAB

addpath(genpath('sedproxy'))
close all; clear;

habitat_wts = readmatrix('seasonality.csv'); habitat_wts = habitat_wts(:,2); habitat_wts = habitat_wts(2:end);
timein = readmatrix('timein.csv'); timein = timein(:,2); timein = timein(2:end);
timepts = readmatrix('timepts.csv'); timepts = timepts(:,2); timepts = timepts(2:end);
sed_acc_rate = readmatrix('sed_rate.csv'); sed_acc_rate = sed_acc_rate(:,2); sed_acc_rate = sed_acc_rate(2:end);
clim_signal = readmatrix('climin.csv'); clim_signal = clim_signal(2:end,2:end);

%t = -max(timein):1:-min(timein); % this needs to be in yrs 

figure('Position',[10 10 1100 400])

for calib = 1:2

    if calib == 1
        calibration_type = 'Uk37';
        n_samples = Inf;
        meas_bias = 0;
    else
        calibration_type = 'MgCa';
        n_samples = 100000;
        meas_bias = 0;
    end

    sigma_meas = 0; sigma_ind = 0;

    % make time go forward
    %time_forward = t-min(t); 

    sed_acc_rate = sed_acc_rate;

    if calib == 2
        [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
            proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
            proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,timein(:),timepts,calibration_type,n_samples,...
            sigma_meas, sigma_ind, meas_bias,sed_acc_rate(:),habitat_wts);

    else
        [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
            proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
            proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,timein,timepts,calibration_type,n_samples,...
            sigma_meas, sigma_ind, meas_bias,sed_acc_rate(:));
    end

    %time_out = time_out-min(time_out); % make go forwardw

    if calib == 1
        Rproxy = readmatrix('simulated_Uk37.csv');
    else
        Rproxy = readmatrix('simulated_MgCa.csv');
    end

    Rproxy = Rproxy(2:end,2:end);

    subplot(1,2,calib)
    plot(time_out/1000,clim_signal_ann,'k','LineWidth',2); hold on;
    plot(time_out/1000,reconstructed_climate,'r','linewidth',1.5); hold on;
    plot(Rproxy(:,1)/1000,Rproxy(:,end),'color',rgb('light blue'),'linewidth',1.5); 

    
    grid on; set(gca,'fontsize',12);
    ylabel('temperature (^oC)'); xlabel('time (kya)')
    if calib==1
        legend('raw temp','reconstructed (Becca)','reconstructed (Thom)','Location','northwest')
        title('Uk37')
    else
        title('MgCa')
    end

    set(gca,'XDir','Reverse')

end