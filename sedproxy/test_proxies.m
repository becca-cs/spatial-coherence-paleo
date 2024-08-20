N = 10000;
clim_signal = randn(N,12);
time_in = 0:N-1;
timepoints = 0:10:N-1;
calibration_type = 'MgCa';
n_samples = 100;
sigma_meas = 0; sigma_ind = 0; meas_bias = 0;
sed_acc_rate = 50*ones(length(timepoints),1);

[time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal,proxy_bt,proxy_bt_sb,...
    proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
    proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate] = ClimToProxyClim(clim_signal,time_in,timepoints,calibration_type,n_samples,...
    sigma_meas, sigma_ind, meas_bias,sed_acc_rate);

%%
sed_mean = mean(sed_acc_rate); % this is in cm per kyr! need to transform into yr per cm
sed_std = 100;

[depth_sample,depth_interp,age_model,age_true,sed_true] = age_depth_model(timepoints,1/(sed_mean/1000),1/(sed_std/1000));

%%
figure('Position',[10 10 1000 400])
subplot(1,2,1)
plot(time_out,clim_signal_ann,'k','linewidth',1); hold on;
plot(time_out,reconstructed_climate,'r','linewidth',1)
legend('climate','reconstructed climate');
set(gca,'fontsize',14)

subplot(1,2,2)
plot(timepoints,mean(proxy_clim_signal(timepoints+1,:),2),'k','linewidth',1); hold on;
plot(time_out,proxy_bt,'r','linewidth',1); 
plot(time_out,proxy_bt_sb,'g','linewidth',1); 
plot(time_out,proxy_bt_sb_inf_b,'y','linewidth',1); 
plot(time_out,proxy_bt_sb_inf_b_n,'b','linewidth',1); 
plot(age_model,proxy_bt_sb_inf_b_n)
legend('input','bioturbation',...
    'bioturbation + habitat','bioturbation + habitat + calib bias','bioturbation + habitat + calib bias + meas error')

set(gca,'fontsize',14);