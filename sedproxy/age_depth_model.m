function[depth_sample,depth_interp,age_model,age_true,sed_true,age_model_timepoints,age_true_timepoints] = age_depth_model(time_in,sed_mean,sed_std)

% sed_mean is in yr/cm
% sed_std is in yr/cm

%sed_mean = sed_mean./1000; sed_std = sed_std./1000;

time_in = time_in(:);
max_time = max(abs(time_in));

x = sed_mean; y = sed_std; %1/sed_mean; y = 1/sed_std; % now in yrs per cm

% Parameters for sediment accumulation!
a = x^2./y^2; % maybe should have these vary more across sections? but probably okay for now...
b = y^2./x;
% Mean value is a/b (yr per cm) - this should be somewhere b/w 10 (100 cm/kyr) & 100 (10 cm/kyr)
% Variance is a/b^2 - for now, just kinda play around with this!

% also can have w be drawn from a distribution
% beta distribution - the correlation of accumulation rates of any two sections of the core separated by ds cm
a_w = 7;
b_w = 3;
ds = 1; %separation in cm
delc = 10; %how often to change the accumulation rate
%high memory

% for low memory, b = 1/25, aw = 2, bw = 28

w = 0.01:0.01:1; % possible memories
fR = betapdf(w.^(ds/delc),a_w,b_w);
fR = fR/sum(fR);
fw = ds*w.^(ds/delc-1)/delc;
fw = fw./sum(fw);
f_w = fR.*fw;
f_w = f_w./sum(f_w);

% build distribution
w_samp = randsample(w, round(max_time), true, f_w);

% Parameters for age-depth model
age_samp = 10; % how many age samples to take
age_err = 40; % yrs - what is the error in each radiocarbon date? should do some research on this distribution
d_depth = 1; % how often to sample SST (cm)

% Produce sedimentation record!
dt = 1;  % assuming an annual resolution
yr = 1:dt:max_time;
alpha = gamrnd(a,b,round(max_time),1);
sed_rate = zeros(1,round(max_time));
sed_rate(1) = a*b; % just take the mean value 
total_depth = sed_rate;

%%

for i = 2:round(max_time) 
  ind = floor((sum(1./sed_rate(1:i-1)))./delc)+1;
  sed_rate(i) = w_samp(ind)*sed_rate(i-1) + (1-w_samp(ind))*alpha(ind);
end

sed_rate_yr = 1./sed_rate;

for i = 2:round(max_time)
  total_depth(i) = sum(sed_rate_yr(2:i))*dt;
end

total_depth(1) = 0; total_depth = fliplr(total_depth);

seg = max(total_depth):-max(total_depth)/age_samp:0;
idx(1) = 1;
  
age = yr - max(yr);

for i = 2:length(seg)
    idx(i) = find(abs(total_depth - seg(i))==min(abs(total_depth - seg(i)))); %which.min(abs(total_depth - seg[i])))
end

ct = 1;
err1 = age(1);
time_err = zeros(1,length(idx)); 

% Add dating uncertainty, normally distributed (bounded)
for i = idx
  if i==1
      time_err(1) = age(i) + age_err*randn();
      err1 = time_err(ct) + 1; ct = ct+1;
  elseif(i<length(total_depth)) 
      flg = 1;
      while flg
          err = age(i) + age_err*randn();
          if err>err1 && err<0
              time_err(ct) = err; err1 = time_err(ct) + 1; ct = ct+1;
              flg = 0;
          end
      end
  else
      time_err(ct) = age(i);
  end
end

% Okay, this should now predict AGES for every SST measurement
depth_interp = max(total_depth):-d_depth:0; % interpolation depths
depth_sample = total_depth(idx); % sampled depths
age_model = interp1(depth_sample,time_err,depth_interp); 
age_true = interp1(total_depth,age,depth_interp);
sed_true = interp1(total_depth,sed_rate,depth_interp);

depth_model_timepoints = interp1(age_model,depth_interp,time_in);
depth_true_timepoints = interp1(age_true,depth_interp,time_in);

% just at timepoints of interest!
age_model_timepoints = interp1(depth_interp,age_model,depth_model_timepoints);
age_true_timepoints = interp1(depth_interp,age_true,depth_true_timepoints);

end
