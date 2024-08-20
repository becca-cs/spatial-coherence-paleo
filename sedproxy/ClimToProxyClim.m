% #'   The following aspects of proxy creation are currently modelled.
%      Based on sedproxy
%      (https://cp.copernicus.org/preprints/cp-2018-13/cp-2018-13.pdf)
% #'
% #'   1. Seasonal bias in the encoding of a proxy due to the interaction between
% #'   climate seasonality and any seasonality in the life cycle of the organism
% #'   encoding the climate signal (e.g. Foraminifera for Mg/Ca ratios, or
% #'   phytoplankton for Alkenone unsaturation indices).
% #'
% #'   2. Bioturbation of the sediment archived proxy. For each requested
% #'   timepoint, the simulated proxy consists of a weighted mean of the climate
% #'   signal over a time window that is determined by the sediment accumulation
% #'   rate \emph{sed.acc.rate} and the bioturbation depth \emph{bio.depth} which defaults
% #'   to 10 cm. The weights are given by the depth solution to an impulse
% #'   response function (Berger and Heath, 1968).
% #'
% #'   3. Aliasing of seasonal and inter-annual climate variation onto to
% #'   bioturbated (smoothed) signal. For proxies measured on a small number of
% #'   discrete particles both seasonal and inter-annual climate variation is
% #'   aliased into the proxy record. For example, Foraminifera have a life-cycle
% #'   of approximately 1 month, so they record something like the mean
% #'   temperature from a single month. If Mg/Ca is measured on e.g.
% #'   \code{n.samples} = 30 individuals, the measured proxy signal is a mean of
% #'   30 distinct monthly mean temperatures and will thus be a stochastic sample
% #'   of the true mean climate.
% #'
% #'   4. Measurement noise/error is added as a pure Gaussian white noise process
% #'   with mean = 0, standard deviation =
% #'    \code{sqrt(sigma.meas^2 + sigma.ind^2/n.samples)}.
% #'
% #'   5. Additionally, a random *bias* can be added to each realisation of a
% #'   proxy record. Bias is simulated as a Gaussian random variable with mean =
% #'   0, standard deviation = \code{meas.bias}. The same randomly generated bias
% #'   value is applied to all timepoints in a simulated proxy record, when
% #'   multiple replicate proxies are generated (\emph{n.replicates} > 1) each
% #'   replicate has a different bias applied.
% #'
% #'   \code{ClimToProxyClim} returns one or more replicates of the final
% #'   simulated proxy as well as several intermediate stages (see section
% #'   **Value** below).
% #'
% #' @param clim.signal The "assumed true" climate signal, e.g. climate model
% #'   output or instrumental record. A \code{\link{ts}} object consisting of a
% #'   years x 12 (months) x n habitats (e.g. depths) matrix of temperatures. The
% #'   time series should be at annual resolution and in reverse, i.e. "most
% #'   recent timepoint first" order.
% #' @param timepoints The timepoints for which the proxy record is to be modelled
% #' @param calibration.type Type of proxy, e.g. Uk'37 or MgCa, to which the
% #'   clim.signal is converted before the archiving and measurement of the proxy
% #'   is simulated. Defaults to "identity" which means no conversion takes place.
% #' @param noise.type Determines whether additive or multiplicative measurement
% #'   noise is added. The appropriate type depends on the units of the proxy.
% #'   Defaults to multiplicative for MgCa, additive for Uk'37 and identity (none)
% #'   calibration types. Can be overidden with a string, "additive" or
% #'   "multiplicative" in the case that pre-converted climate signal and
% #'   measurement noise values are used in combination with an "identity"
% #'   calibration type.
% #' @param plot.sig.res The resolution, in years, of the smoothed (block
% #'   averaged) version of the input climate signal returned for plotting. This
% #'   does not affect what the proxy model uses as input. If set to NA, no
% #'   smoothed climate output is generated, this can speed up some simulations.
% #' @param habitat.weights Production weights for the proxy / proxy-carrier
% #' either as a vector of values with length = ncol(clim.signal), i.e. 1 weight
% #' for each month x habitat combination, a matrix of the same dimensions as the
% #' input climate signal matrix, or a function that produces an index of
% #' productivity as a function of temperature.
% #' Defaults to a vector of length = ncol(clim.signal) of equal weights.
% #' @param habitat.wt.args A named list of parameter values to be passed to
% #' a function named in habitat.weights.
% #' @param bio.depth Depth of the bioturbated layer in cm, defaults to 10 cm.
% #' @param layer.width the width of the sediment layer from which samples were
% #'   taken, e.g. foraminifera were picked or alkenones were extracted, in cm.
% #'   Defaults to 1 cm. If bio.depth and layer.width are both set to zero,
% #'   each timepoint samples from a single year of the clim.signal, equivalent to
% #'   sampling an annually laminated sediment core.
% #' @param sed.acc.rate Sediment accumulation rate in cm per 1000 years. Defaults
% #'   to 50 cm per ka. Either a single value, or vector of same length as
% #'   "timepoints"
% #' @param sigma.meas The standard deviation of the measurement error
% #' added to each simulated proxy value.
% #' @param sigma.ind The standard deviation of error between individuals
% #' (e.g. Forams) not otherwise modelled. This could included "vital effects" or
% #' aliasing of depth habitat variation not modelled via a depth resolved input
% #' climate signal and habitat weights. sigma.ind is scaled by n.samples
% #' before being combined with sigma.meas.
% #' @param n.samples Number of e.g. Foraminifera sampled per timepoint, this can
% #'   be either a single number, or a vector of length = timepoints. Can be set
% #'   to Inf for non-discrete proxies, e.g. for Uk’37.
% #' @param meas.bias The amount of bias to add to each simulated proxy
% #'   time-series. Each replicate proxy time-series has a constant bias added,
% #'   drawn from a normal distribution with mean = 0, sd = meas.bias. Bias
% #'   defaults to zero.
% #' @param scale.noise Scale noise to proxy units. Defaults to TRUE if
% #' calibration.type is not "identity"
% #' @param n.replicates Number of replicate proxy time-series to simulate from
% #'   the climate signal
% #' @param n.bd Number of multiples of the bioturbation width at which to truncate
% #' the bioturbation filter
% #' @param top.of.core The theoretical minimum age at the top of the core, ie.
% #' the year the core was sampled, defaults to the start of clim.in
% #' @inheritParams ProxyConversion
% #' @return \code{ClimToProxyClim} returns an object of class "sedproxy.pfm", a list with three elements:
% #'
% #'   1. a dataframe \code{simulated.proxy}
% #'   2. a dataframe \code{smoothed.signal}
% #'   3. a dataframe \code{everything}
% #'
% #'
% #'   The dataframe \code{simulated.proxy} contains a single realisation of the
% #'   final forward modelled proxy, as well as the intermediate stages and the
% #'   original climate signal at the requested timepoints.
% #'
% #'   The dataframe \code{smoothed.signal} contains a block averaged version the
% #'   input climate signal, defaults to 100 year means but this is set by the
% #'   parameter plot.sig.res. This is useful for plotting against the
% #'   resulting simulated proxy.
% #'
% #'   The dataframe \code{everything} contains all of the above but with multiple
% #'   replicates of the pseudo-proxy records if requested. The data are in
% #'   "long form", with the column "stage" inidcating the proxy stage or input
% #'   climate resolution and column "value" giving the values.
% #'
% #' **Named elements of the returned proxy record:**
% #'
% #' \describe{
% #'    \item{timepoints}{Requested timepoints}
% #'    \item{clim.signal.ann}{Input climate signal at requested timepoints at annual resolution}
% #'    \item{clim.signal.smoothed}{Input climate signal at regular time intervals and resolution = plot.sig.res}
% #'    \item{clim.timepoints.ssr}{Input climate signal at requested timepoints, smoothed to resolution = plot.sig.res}
% #'    \item{proxy.bt}{Climate signal after bioturbation}
% #'    \item{proxy.bt.sb}{Climate signal after bioturbation and habitat bias}
% #'    \item{proxy.bt.sb.inf.b}{Climate signal after bioturbation, habitat bias, and calibration bias}
% #'    \item{proxy.bt.sb.inf.b.n}{Climate signal after bioturbation, habitat bias, and measurement error}
% #'    \item{proxy.bt.sb.sampY}{Climate signal after bioturbation, habitat bias, and aliasing of inter-annual variation}
% #'    \item{proxy.bt.sb.sampYM}{Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats}
% #'    \item{proxy.bt.sb.sampYM.b}{Climate signal after bioturbation, habitat bias, and aliasing of inter-annual and intra-annual variation such as monthly temperatures or depth habitats, and calibration bias}
% #'    \item{proxy.bt.sb.sampYM.b.n}{Climate signal after bioturbation, habitat bias, aliasing, and measurement error}
% #'    \item{simulated.proxy}{Final simulated pseudo-proxy, this will be same as proxy.bt.sb.inf.b.n when n.samples = Inf, and proxy.bt.sb.sampYM.b.n when n.samples is finite}
% #'    \item{observed.proxy}{True observed proxy (when supplied)}
% #' }
% #'
% #' @importFrom dplyr rename
% #' @importFrom rlang .data
% #' @export
% #'
% #'@examples
% #' library(ggplot2)
% #' set.seed(26052017)
% #' clim.in <- ts(N41.t21k.climate[nrow(N41.t21k.climate):1,] - 273.15, start = -39)
% #'
% #' PFM <- ClimToProxyClim(clim.signal = clim.in,
% #'                        timepoints = round(N41.proxy$Published.age),
% #'                        calibration.type = "identity",
% #'                        habitat.weights = N41.G.ruber.seasonality,
% #'                        sed.acc.rate = N41.proxy$Sed.acc.rate.cm.ka,
% #'                        layer.width = 1,
% #'                        sigma.meas = 0.46,
% #'                        sigma.ind = 0, n.samples = Inf,
% #'                        plot.sig.res = 10, meas.bias = 1,
% #'                        n.replicates = 10)
% #'
% #' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "seq") +
% #'   facet_wrap(~stage)
% #'
% #' PlotPFMs(PFM$everything, max.replicates = 1, stage.order = "var")
% #'
% #' PlotPFMs(PFM$everything, stage.order = "var", plot.stages = "all")
% #'
function [time_out,clim_signal_ann,clim_timepoints_ssr,proxy_clim_signal_ann,proxy_bt,proxy_bt_sb,...
    proxy_bt_sb_inf_b,proxy_bt_sb_inf_b_n,proxy_bt_sb_sampY,proxy_bt_sb_sampYM,...
    proxy_bt_sb_sampYM_b,proxy_bt_sb_sampYM_b_n,reconstructed_climate,valid_inds] = ClimToProxyClim(clim_signal,time_in,timepoints,calibration_type,n_samples,...
    sigma_meas, sigma_ind, meas_bias, varargin)
if isempty(varargin)
    sed_acc_rate = 50;
    habitat_weights = 1/size(clim_signal,2)*ones(1,size(clim_signal,2));
    bio_depth = 10;
    habitat_wt_args = [];
    proxy_units = 1;
elseif length(varargin)==1
    sed_acc_rate = varargin{1};
    habitat_weights = 1/size(clim_signal,2)*ones(1,size(clim_signal,2));
    bio_depth = 10;
    habitat_wt_args = [];
    proxy_units = 1;
elseif length(varargin)==2
    sed_acc_rate = varargin{1};
    habitat_weights = varargin{2};
    bio_depth = 10; habitat_wt_args = [];
    proxy_units = 1;
elseif length(varargin)==3
    sed_acc_rate = varargin{1}(:);
    habitat_weights = varargin{2}(:);
    bio_depth = varargin{3}(:);
    habitat_wt_args = [];
    proxy_units = 1;
elseif length(varargin)==4
    sed_acc_rate = varargin{1}(:);
    habitat_weights = varargin{2};
    bio_depth = varargin{3}(:);
    habitat_wt_args = varargin{4}; % args for habitat function
    proxy_units = 1;
else
    sed_acc_rate = varargin{1}(:);
    habitat_weights = varargin{2};
    bio_depth = varargin{3}(:);
    habitat_wt_args = varargin{4}; % args for habitat function
    proxy_units = varargin{5};
end

if isa(habitat_weights,'function_handle') && isempty(habitat_wt_args)
    habitat_species = 'sacculifer'; % pick a species :)
elseif isa(habitat_weights,'function_handle')
    habitat_species = habitat_wt_args; % 
end

slp_int_means = [];
slp_int_vcov = [];
% plot.sig.res = 100,
% habitat.weights = rep(1/ncol(clim.signal),
% ncol(clim.signal)),
% habitat.wt.args = NULL,
% bio.depth = 10,
% sed.acc.rate = 50,
layer_width = 1;
% sigma.meas = 0,
% sigma.ind = 0,
% meas.bias = 0,
% n.samples = Inf,

n_replicates = 1;
n_timepoints = length(timepoints);
top_of_core = NaN;
n_bd = 3;
plot_sig_res = 100;
cal_uncert = 1;

% make sure stuff is the right size!
time_in = time_in(:); 
timepoints = timepoints(:); 

% input some standard stuff

if (length(n_samples) ~= 1 && length(n_samples)~=n_timepoints)
    warning(['n_samples must be either a single value, or a vector the same',...
        ' length as timepoints'])
end

if ~ismember(calibration_type,{'identity','Uk37','MgCa'})
    warning('Please input a calibration type: identity, Uk37, or MgCa');
end

if ismember(calibration_type,'identity')
    calibration = NaN; noise_type = 'additive';
    scale_noise = 0;
elseif ismember(calibration_type,'Uk37')
    calibration = 'Mueller global'; noise_type = 'additive';
    scale_noise = 1;
else
    calibration = 'Ten planktonic species_350-500';
    noise_type = 'multiplicative';
    scale_noise = 1;
end

if strcmp(calibration_type,'Uk37') && sum(isfinite(n_samples))
    n_samples = Inf;
    warning('n_samples should be infitinite; fixed input')
end

if strcmp(calibration_type,'MgCa') && sum(~isfinite(n_samples))
    n_samples = 30;
    warning('n_samples should be fitinite; fixed input')
end


if (length(sigma_meas) ~= 1 && length(sigma_meas) ~= n_timepoints)
    warning(['sigma_meas must be either a single value, or a vector the same',...
        ' length as timepoints'])
end

if (length(sigma_ind) ~= 1 && length(sigma_ind)~=n_timepoints)
    warning(['sigma_ind must be either a single value, or a vector the same',...
        ' length as timepoints'])
end

if (length(sed_acc_rate)~=1 && length(sed_acc_rate) ~= n_timepoints)
    warning(['length of sed_acc_rate must be either a single value, or a vector the same',...
        'length as timepoints']);
end

if (length(habitat_weights) ~= size(clim_signal,2) | size(habitat_weights) ~= size(clim_signal))
    warning(['habitat.weights must be either a vector of weights with length = ncol(clim.signal),',...
        'a matrix of weights of the same dimensions as the input climate signal, or a function.',...
        'Function names should be given unquoted, e.g. dnorm, not \dnorm\']);
end

if isnan(top_of_core)
    top_of_core = time_in(1);
elseif top_of_core < time_in(1)
    warning('top_of_core cannot be younger than the start of clim_signal');
end

% -------- Convert to proxy units if requested --------

if strcmp(calibration_type,"identity")==0
    proxy_clim_signal = ProxyConversion(clim_signal, [], calibration_type,...
        slp_int_means, slp_int_vcov, "point", 1);
else
    proxy_clim_signal = clim_signal;
end

% -------- Calculate timepoint invariant values ------
max_clim_signal = time_in(end); %
min_clim_signal = time_in(1);

% -------- Create smoothed climate signal --------
if isnan(plot_sig_res)
    timepoints_smoothed = NaN;
    clim_signal_smoothed = NaN;
else
    timepoints_smoothed = min_clim_signal-1:plot_sig_res:max_clim_signal-1;
    clim_signal_smoothed = ChunkMatrix(timepoints_smoothed, plot_sig_res,...
        clim_signal, time_in); % could also use movmean
end

% -------- Check whether bioturbation window will extend beyond climate signal for any of the timepoints

% bioturbation window will be focal.timepoint - bio.depth.timesteps - layer.width.years / 2 to
% focal.timepoint + 3*bio.depth.timesteps

if size(layer_width)~=size(sed_acc_rate)
    layer_width = layer_width';
end

bio_depth_timesteps = round(1000 * bio_depth./ sed_acc_rate)';
layer_width_years = ceil(1000 * layer_width./ sed_acc_rate)';

max_windows = ceil(timepoints + n_bd.* bio_depth_timesteps');
min_windows = floor(timepoints - bio_depth_timesteps' - layer_width_years'./2);

max_ind = max_windows >= max_clim_signal;
min_ind = min_windows <  min_clim_signal;

% -------- Use Rapid or Slow version ----------------------

if (length(bio_depth)==1 && length(sed_acc_rate)==1 && ...
        length(layer_width)==1 && length(n_samples)==1 && ...
        isa('habitat_weights','function_handle') == 0 && isvector(habitat_weights))

    % ---- Rapid ------
    disp("Using Rapid version")

    % -------- Find mixed layer points ------
    % keep points in the mixed layer as well as those below but still inside time-signal
    tpts_above_core_top = (timepoints < top_of_core);
    valid_inds = (max_ind == 0 & tpts_above_core_top == 0);

    % identify mixed layer
    mixed_layer_inds =  (min_ind == 1 & tpts_above_core_top == 0);
    mixed_layer_inds = (mixed_layer_inds(valid_inds));

    if sum(max_ind)>0 % warning for bioturbation
        warning(['One or more requested timepoints is too old. Bioturbation window(s) for timepoint(s) ',...
            num2str(timepoints(max_ind)'),...
            ' extend(s) beyond end of input climate signal. Returning pseudo-proxy for valid timepoints.']);
    end

    if sum(min_windows < min_clim_signal)>0 % warning for mixed layer
        warning(['Timepoint(s) ', num2str(timepoints(mixed_layer_inds)'),...
            ' are in the mixed layer']);
    end

    if sum(tpts_above_core_top)>0
        warning(['One or more requested timepoints is too recent. Timepoint(s) ',...
        num2str(timepoints(tpts_above_core_top)),...
            ' are more recent than the top of the core.'])
    end

    % remove too old or young timepoints
    timepoints = timepoints(valid_inds);
    n_timepoints = length(timepoints);

    % adjusted timepoints for the mixed layer
    % in the mixed layer the bioturbation window is centred around the
    % bottom of the mixed layer
    timepoints_adj = timepoints;
    timepoints_adj(mixed_layer_inds) = 1 + bio_depth_timesteps + layer_width_years / 2;

    % Scale sigma.ind by n.samples and create combined error term
    sigma_ind_scl = sigma_ind/sqrt(n_samples); sigma_ind_scl(~isfinite(n_samples)) = 0;
    sigma_meas_ind = sqrt(sigma_meas^2 + sigma_ind_scl^2);

    % Ensure seasonal productivities are weights
    habitat_weights = habitat_weights./ sum(habitat_weights);


    % Get relative bioturbation window ----------
    first_tp = -bio_depth_timesteps - layer_width_years / 2;
    last_tp = n_bd * bio_depth_timesteps;

    bioturb_window = first_tp:last_tp;

    % Get bioturbation weights --------
    bioturb_weights = BioturbationWeights(bioturb_window, 0, layer_width,...
        sed_acc_rate, bio_depth,'time');


    % Get bioturbation X no-seasonality weights matrix ---------
    biot_sig_weights = repmat(bioturb_weights',[1, size(proxy_clim_signal,2)]);
    biot_sig_weights = biot_sig_weights./sum(biot_sig_weights(:));

    % Get bioturbation X seasonality weights matrix ---------
    clim_sig_weights = bioturb_weights.*habitat_weights';
    clim_sig_weights = clim_sig_weights./ sum(clim_sig_weights(:));
    clim_sig_weights = clim_sig_weights';

    % Check weights sum to 1, within tolerance
    weight_err = abs(sum(clim_sig_weights(:)) - 1);
    if (weight_err > 1e-10)
        warning('ah! weight error is big')
    end

    % Do sampling ------
    if (isfinite(n_samples))
        % call sample once for all replicates and timepoints together, then take means of
        % groups of n_samples
        % Get indices not values
        tot_n_samples = n_samples * n_replicates * n_timepoints;

        samp_indices = datasample(1:numel(clim_sig_weights),tot_n_samples, ...
            'weights',clim_sig_weights(:)); % I'm pretty sure about this? but not 100%
    end

    % For each timepoint ------
    for tp = 1:n_timepoints

        % Get portion of clim.signal corresponding to bioturbation window for this timepoint -------
        if size(clim_signal,2)>1
            clim_sig_window = proxy_clim_signal(bioturb_window + timepoints_adj(tp) - min_clim_signal +1,:);
        else
            clim_sig_window = proxy_clim_signal(bioturb_window + timepoints_adj(tp) - min_clim_signal +1);
        end

        % Calculate mean clim.signal -------

        % Just bioturbation
        proxy_bt(tp) = sum(biot_sig_weights(:).* clim_sig_window(:));

        % Bioturbation + seasonal bias
        proxy_bt_sb(tp) = sum(clim_sig_weights(:).* clim_sig_window(:));

        % Bioturbation + seasonal bias + aliasing
        if (~isfinite(n_samples))
            proxy_bt_sb_sampY(:,tp) = NaN*zeros(1,n_replicates);
            proxy_bt_sb_sampYM(:,tp) = NaN*zeros(1,n_replicates);
        else
            % Use previously sampled indices
            % get indices for this timepoint
            samp_indices_tp = samp_indices(((tp-1)*n_samples*n_replicates+1):(tp*n_samples*n_replicates));

            % convert vector to matrix (cheap only attributes changed), then means
            % can be taken across columns to get per replicate means

            samp = reshape(clim_sig_window(samp_indices_tp), n_samples,[]);
            proxy_bt_sb_sampYM(:,tp) = mean(samp); % this is probably wrong?

            %Get without seasonal aliasing (bioturbation aliasing only)
            clim_sig_window_ann = sum(clim_sig_window'.* habitat_weights');
            col_indices = mod(samp_indices_tp - 1, size(clim_sig_window, 1)) + 1;

            samp_bt = reshape(clim_sig_window_ann(col_indices),[],n_samples);
            proxy_bt_sb_sampY(:,tp) = mean(samp_bt,2);
        end
    end

else
    % Slow ----
    % Find mixed layer points ------
    % keep points in the mixed layer as well as those below but still inside time-signal
    tpts_above_core_top = timepoints < top_of_core;

    % identify mixed layer
    % find oldest timepoint in mixed layer

    oldest_in_mix = timepoints(min_ind);

    if length(oldest_in_mix)>0
        mixed_layer_inds = (min_ind == 1 & tpts_above_core_top(:) == 0);
    else
        mixed_layer_inds = zeros(1,n_timepoints);
    end

    valid_inds = (max_ind == 0 & tpts_above_core_top(:) == 0);

    if sum(max_ind)>0
        warning(['One or more requested timepoints is too old. Bioturbation window(s) for timepoint(s) ',...
            num2str(timepoints(max_ind)'),...
            ' extend(s) beyond end of input climate signal. Returning pseudo-proxy for valid timepoints.'])
    end
    if sum(mixed_layer_inds)>0
        warning(['Timepoint(s) ',...
            num2str(timepoints(mixed_layer_inds)'), ...
            ' are in the mixed layer']);
    end

    if sum(tpts_above_core_top)>0
        warning(['One or more requested timepoints is too recent. Timepoint(s) ',...
            num2str(timepoints(tpts_above_core_top)),...
            ' are more recent than the top of the core.'])
    end

    timepoints = timepoints(valid_inds);
    n_timepoints = length(timepoints);
    mixed_layer_inds = mixed_layer_inds(valid_inds);

    % Scale sigma.ind by n.samples and create combined error term
    if isfinite(n_samples)
        sigma_ind_scl = sigma_ind / sqrt(n_samples);
    else
        sigma_ind_scl = 0;
    end

    sigma_meas_ind = sqrt(sigma_meas^2 + sigma_ind_scl^2);

    % Create vectors from "scalar" inputs
    if (length(sed_acc_rate) == 1)
        sed_acc_rate = sed_acc_rate * ones(n_timepoints,1);
    end

    if (length(layer_width) == 1)
        layer_width = layer_width*ones(n_timepoints,1);
    end

    if (length(n_samples) == 1)
        n_samples = n_samples*ones(n_timepoints,1);
    end

    if (length(sigma_meas) == 1)
        sigma_meas = sigma_meas*ones(n_timepoints,1);
    end

    if (length(sigma_ind) == 1)
        sigma_ind = sigma_ind*ones(n_timepoints,1);
    end

    % Remove timepoint specific parameters that exceed clim.signal ------

    if (length(sed_acc_rate) > n_timepoints)
        sed_acc_rate = sed_acc_rate(valid_inds);
    end
    if (length(layer_width) > n_timepoints)
        layer_width = layer_width(valid_inds);
    end

    if (length(n_samples) > n_timepoints)
        n_samples = n_samples(valid_inds);
    end

    if (length(sigma_meas) > n_timepoints)
        sigma_meas = sigma_meas(valid_inds);
    end

    if (length(sigma_ind) > n_timepoints)
        sigma_ind = sigma_ind(valid_inds);
    end


    % Generate productivity weights from function if supplied
    % COME BACK AND FIX THIS
    if isa(habitat_weights,'function_handle')
        habitat_function = habitat_weights;
        habitat_weights = habitat_function(habitat_species,clim_signal+ 273.15);
    end

    % If vector ensure habitat.weights are weights and matrix
    if (isvector(habitat_weights))
        habitat_weights = habitat_weights./ sum(habitat_weights);
        habitat_weights = repmat(habitat_weights,[size(clim_signal,1) 1]);
    end

    timepoints_adj = timepoints;

    max_windows = max_windows(valid_inds);
    min_windows = min_windows(valid_inds);

    % set mixed layer sed.acc.rate to the lowest
    if sum(mixed_layer_inds)>0
        sed_acc_rate(mixed_layer_inds) = min(sed_acc_rate(mixed_layer_inds));

        % reset mixed window for mixed layer points
        bio_depth_timesteps = round(1000 * bio_depth./ sed_acc_rate);
        layer_width_years = ceil(1000 * layer_width./ sed_acc_rate);

        % adjusted timepoints for the mixed layer
        % in the mixed layer the bioturbation window is centred around the
        % bottom of the mixed layer

        timepoints_adj(mixed_layer_inds) = 1 + ...
            bio_depth_timesteps(mixed_layer_inds(:)) + ...
            layer_width_years(mixed_layer_inds) / 2;


        max_windows(mixed_layer_inds) = ceil((n_bd+1) * bio_depth_timesteps(mixed_layer_inds) +...
            layer_width_years(mixed_layer_inds) / 2);
        min_windows(mixed_layer_inds) = top_of_core;

    end

    % For each timepoint ------

    for tp = 1:n_timepoints

        % Get bioturbation window ----------
        first_tp = min_windows(tp);
        last_tp = max_windows(tp);
        bioturb_window = first_tp:last_tp;

        % Get bioturbation weights --------
        bioturb_weights = BioturbationWeights(bioturb_window, timepoints_adj(tp),...
            layer_width(tp), sed_acc_rate(tp),bio_depth,'time');

        if size(proxy_clim_signal,2)>1
            clim_sig_window =  proxy_clim_signal(time_in >= first_tp & time_in <= last_tp,:);
        else
            clim_sig_window =  proxy_clim_signal(time_in >= first_tp & time_in <= last_tp);
        end

        % Get bioturbation X no-seasonality weights matrix ---------
        biot_sig_weights = bioturb_weights'.*ones(1, size(proxy_clim_signal,2));
        biot_sig_weights = biot_sig_weights./ sum(biot_sig_weights(:));



        % Get bioturbation X seasonality weights matrix ---------
        if size(habitat_weights,2)>1
            habitat_weights_mat = habitat_weights(time_in >= first_tp & time_in <= last_tp,:);
        else
            habitat_weights_mat = habitat_weights(time_in >= first_tp & time_in <= last_tp);
        end

        clim_sig_weights = repmat(bioturb_weights',[1, size(clim_signal,2)]).* habitat_weights_mat;
        clim_sig_weights = clim_sig_weights./sum(clim_sig_weights(:));

        % Check weights sum to 1, within tolerance
        weight_err = abs(sum(clim_sig_weights(:)) - 1);
        if (weight_err > 1e-10)
            warning('weight error is big!')
        end


        % Calculate mean clim.signal -------

        % Just bioturbation
        proxy_bt(tp) = sum(biot_sig_weights(:).* clim_sig_window(:));

        % Bioturbation + seasonal bias
        proxy_bt_sb(tp) = sum(clim_sig_weights(:).* clim_sig_window(:));

        % Bioturbation + seasonal bias + aliasing
        if (~isfinite(n_samples(tp)))
            proxy_bt_sb_sampY(tp) = NaN; %*ones(1,n_replicates);
            proxy_bt_sb_sampYM(tp) = NaN; %*ones(1,n_replicates);
        else

            % call sample once for all replicates together, then take means of
            % groups of n.samples
            % Get indices not values
            % account for function being weird - if no habitat weights
            % available, just weight everything equally...
            if sum(clim_sig_weights(:)>0) == 0 || sum(~isnan(clim_sig_weights(:))) == 0
                clim_sig_weights = repmat(bioturb_weights',[1, size(clim_signal,2)]);
                clim_sig_weights = clim_sig_weights./sum(clim_sig_weights(:));
            end
           
            samp_indices =  datasample(1:numel(clim_sig_window),n_samples(tp)*n_replicates, ...
                'weights',clim_sig_weights(:));

            %# convert vector to matrix (cheap only attributes changed), then means
            %# can be taken across columns to get per replicate means
            samp = reshape(clim_sig_window(samp_indices), n_samples(tp),[]);
            proxy_bt_sb_sampYM(tp) = mean(samp);

            % Get without seasonal aliasing (bioturbation aliasing only)
            % habitat weights rows need to sum to 1
            habitat_weights_r1 = habitat_weights_mat./ sum(habitat_weights_mat,2);

            clim_sig_window_ann = sum(clim_sig_window.* habitat_weights_r1,2);
            row_indices = mod(samp_indices-1, size(clim_sig_window,1)) + 1;
            samp_bt = clim_sig_window_ann(row_indices, :);
            samp_bt = reshape(samp_bt,n_samples(tp),[]);
            proxy_bt_sb_sampY(tp) = mean(samp_bt);
        end
    end
end

% not sure if the formatting for the next paragraph is right?
if (n_replicates == 1)
    proxy_bt_sb_sampYM = reshape(proxy_bt_sb_sampYM,1,[]);
end
%proxy_bt_sb_sampYM = proxy_bt_sb_sampYM';

if (n_replicates == 1)
    proxy_bt_sb_sampY = reshape(proxy_bt_sb_sampY,1,[]);
end
%proxy_bt_sb_sampY = proxy_bt_sb_sampY';

% Add bias and noise --------
% Rescale noise if using a calibration -----
if (scale_noise == 1)
    % scale.noise will either be TRUE because a non-identity calibration is being used
    % or it will be a string to identify the correct calibration if the input time-series
    % has already been converted.
    disp('Rescaling noise')

    % If cal type is identity re-scaling still required
    if strcmp(calibration_type,'identity')
        pct = scale_noise;
    else
        pct = calibration_type;
    end

    % mean temperature in temperature units at each timepoint - use bioturbated signal
    mean_temperature =  ProxyConversion([],proxy_bt,...
        pct, slp_int_means, slp_int_vcov,'point',1); 

    sigma_meas_ind = ProxyConversion(mean_temperature + sigma_meas_ind,[],...
        pct,slp_int_means, slp_int_vcov,'sample',1) - ...
        ProxyConversion(mean_temperature,[],...
        pct, slp_int_means, slp_int_vcov,'sample',1);

    if strcmp(noise_type,'multiplicative')
        % noise SD needs to be divided by the mean temperature in proxy units in
        % order to maintain a consistent SD in temperature units.
        sigma_meas_ind = sigma_meas_ind'./ proxy_bt;
    end
end

if strcmp(noise_type,'additive')
    noise = 0 + sigma_meas_ind.*randn(n_replicates * n_timepoints,1);

    if (meas_bias ~= 0)
        bias = 0 + meas_bias.*randn(n_replicates,1);
    else
        bias = zeros(n_replicates,1);
    end

    % Add bias and noise to infinite sample --------

    proxy_bt_sb_inf_b = proxy_bt_sb + bias;
    proxy_bt_sb_inf_b_n = proxy_bt_sb_inf_b + noise';

    if (sum(isfinite(n_samples)))==length(n_samples)
        proxy_bt_sb_inf_b = NaN;
        proxy_bt_sb_inf_b_n = NaN;
    end

    % Add bias and noise to finite sample --------

    proxy_bt_sb_sampYM_b = proxy_bt_sb_sampYM + bias;
    proxy_bt_sb_sampYM_b_n = proxy_bt_sb_sampYM_b + noise';

    % set intermediate bias stages to NA if no bias modelled
    if (meas_bias == 0)
        proxy_bt_sb_inf_b = NaN.*proxy_bt_sb;
        proxy_bt_b_sampYM_b = NaN.*proxy_bt_sb;
    end

elseif strcmp(noise_type,'multiplicative')
    noise = exp(sigma_meas_ind'.*randn(n_replicates * n_timepoints,1));

        if (meas_bias ~= 0)
            bias = exp(0+ meas_bias.*randn(n_replicates,1));
        else
            bias = ones(1, n_replicates);
        end

        % Add bias and noise to infinite sample --------

        proxy_bt_sb_inf_b = proxy_bt_sb.*bias;
        proxy_bt_sb_inf_b_n = proxy_bt_sb_inf_b.* noise';

        if sum(isfinite(n_samples))==length(n_samples)
            proxy_bt_sb_inf_b= proxy_bt_sb.*NaN;
            proxy_bt_sb_inf_b_n = proxy_bt_sb.*NaN;
        end

        % Add bias and noise to finite sample --------
        proxy_bt_sb_sampYM_b = proxy_bt_sb_sampYM.* bias;
        proxy_bt_sb_sampYM_b_n = proxy_bt_sb_sampYM_b.* noise';

        % set intermediate bias stages to NA if no bias modelled
        if (meas_bias == 0)
            proxy_bt_sb_inf_b = NaN.*proxy_bt_sb;
            proxy_bt_sb_sampYM_b = NaN.*proxy_bt_sb;
        end
end


% Create smoothed climate signal -----
if (isnan(plot_sig_res))
    clim_timepoints_ssr = NaN;
else
    clim_timepoints_ssr = ChunkMatrix(timepoints, plot_sig_res, clim_signal,time_in);
end

% Add items to output list -----------
for i = 1:length(timepoints)
    idx(i) = find(timepoints(i)==time_in);
end
clim_signal_ann = sum(clim_signal(idx,:),2)./size(clim_signal,2); % needs to be in
proxy_clim_signal_ann = sum(proxy_clim_signal(idx,:),2)./size(clim_signal,2);

if (sum(isfinite(n_samples))==length(n_samples))
    simulated_proxy = proxy_bt_sb_sampYM_b_n;
else
    simulated_proxy = proxy_bt_sb_inf_b_n;
end

% Add calibration uncertainty -------
% If n.replicates > 1
% First convert back to temperature units with fixed parameters
% Then re-convert to proxy units with random parameters
if ~strcmp(calibration_type,'identity')

    if n_replicates > 1
        simulated_proxy_cal_err = ProxyConversion([],simulated_proxy,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);

        simulated_proxy_cal_err = ProxyConversion(simulated_proxy_cal_err,...
            calibration_type,slp_int_means,slp_int_vcov,"sample", n_replicates);
    else
        out_simulated_proxy_cal_err = simulated_proxy;
    end

    % Do this in all cases, not just if n.replicates == 1
    reconstructed_climate = ProxyConversion([],simulated_proxy,...
        calibration_type, slp_int_means,slp_int_vcov,"point", 1);
else
    out_simulated_proxy_cal_err = simulated_proxy;
    reconstructed_climate = simulated_proxy; 
end

time_out = timepoints;

if proxy_units==0

    proxy_bt = ProxyConversion([],proxy_bt,...
        calibration_type,slp_int_means,slp_int_vcov,"point", 1);
    proxy_bt_sb = ProxyConversion([],proxy_bt_sb,...
        calibration_type,slp_int_means,slp_int_vcov,"point", 1);

    if ismember(calibration_type,{'MgCa'})
        proxy_bt_sb_sampY = ProxyConversion([],proxy_bt_sb_sampY,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);
        proxy_bt_sb_sampYM = ProxyConversion([],proxy_bt_sb_sampYM,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);
        proxy_bt_sb_sampYM_b = ProxyConversion([],proxy_bt_sb_sampYM_b,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);
        proxy_bt_sb_sampYM_b_n = ProxyConversion([],proxy_bt_sb_sampYM_b_n,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);
    elseif ismember(calibration_type,{'Uk37'})
        proxy_bt_sb_inf_b = ProxyConversion([],proxy_bt_sb_inf_b,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);
        proxy_bt_sb_inf_b_n = ProxyConversion([],proxy_bt_sb_inf_b_n,...
            calibration_type,slp_int_means,slp_int_vcov,"point", 1);
    end

end

% slp_int_means =
% if (isnan(slp_int_means) && calibration_type ~= "identity")
%     if (calibration_type == "MgCa" & isnan(calibration))
%         calibration = "Ten planktonic species_350-500";
%     else
%         calibration
%     end
%
%             cp <- data.frame(calibration.parameters)
%             cfs.vcov <- cp[cp$calibration.type == calibration.type &
%                 cp$calibration == calibration,]
%             matrix(c(cfs.vcov$slope, cfs.vcov$intercept),
%             ncol = 2,
%             byrow = TRUE)
%         else
%             slp.int.means
%         end
%
%         calibration.pars <- list(calibration.type = calibration.type,
%         calibration = calibration,
%         slp.int.means = slp.int.means)
%
%         attr(simulated.proxy, "calibration.pars") <-  calibration.pars
%         attr(everything, "calibration.pars") <-  calibration.pars


end