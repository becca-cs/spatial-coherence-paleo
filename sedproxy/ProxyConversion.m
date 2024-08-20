function [out] = ProxyConversion(temperature, proxy_value,...
                            calibration_type,...
                            slp_int_means, slp_int_vcov,...
                            point_or_sample, n)

%CALIB_FUNCTIONS Summary of this function goes here
%   Detailed explanation goes here
% #' Convert between Temperature in Degrees C and Proxy Units
% #'
% #' @description A wrapper function for accessing proxy - temperature conversion functions
% #'
% #' @param temperature Temperature in degrees C
% #' @param proxy.value Temperature in proxy units
% #' @param calibration.type Type of proxy, e.g. Uk37 or MgCa
% #' @param point.or.sample Use the "best estimate" calibration parameters, or
% #'   parameters sampled from the fitted calibration model
% #' @param n the number of replicate conversions to make in the case of sampled
% #'   calibration parameters
% #' @param calibration The name of a specific calibration for which calibration parameters
% #'   are provided by sedproxy. Currently applies only to calibration.type MgCa.
% #' @param slp.int.means Optional user supplied vector of values for the slope
% #'   and intercept of the calibration function. Overides the defaults.
% #' @param slp.int.vcov Optional user supplied variance covariance matrix
% #'   calibration parameters. Overides the defaults.
% #' @details Valid entries for calibration are: "Ten planktonic species_350-500", "G. aequilateralis_350-500", "G.
% #'   aequilateralis_500-1000", "G. conglobatus_350-500", "G. hirsuta_350-500",
% #'   "G. inflata_350-500", "G. ruber pink_250-350", "G. ruber pink_350-500", "G.
% #'   ruber white_250-350", "G. ruber white_350-500", "G. sacculifer with
% #'   sac_350-500", "G. sacculifer without sac_350-500", "G.
% #'   truncatulinoides_350-500", "G. truncatulinoides_500-1000", "N.
% #'   dutertrei_350-500", "O. univesa_350-500", "P. obliquiloculata_350-500"
% #' @return a vector of temperatures or proxy values
% #' @export
% #' @importFrom mvtnorm rmvnorm
% #' @family calib
% #' @examples
% #' # From temperature to UK'37
% #' ## With fixed calibration
% #' ProxyConversion(temperature = c(10, 20), point.or.sample = "point",
% #'                 calibration.type = "Uk37")
% #' 
% #' ## With random calibration, 5 replicates
% #' ProxyConversion(temperature = c(1, 2), n = 5, point.or.sample = "sample",
% #'                 calibration.type = "Uk37")
% #' 
% #' 
% #' ## Back-transformation with same calibration
% #' ProxyConversion(
% #'   proxy.value = as.vector(
% #'     ProxyConversion(
% #'       temperature = c(21, 22),
% #'       calibration.type = "Uk37",
% #'       point.or.sample = "point"
% #'     )
% #'   ),
% #'   point.or.sample = "point",
% #'   calibration.type = "Uk37"
% #' )
% #' 
% #' ## Back-transformation with random calibration
% #' ProxyConversion(
% #'   proxy.value = as.vector(
% #'     ProxyConversion(
% #'       temperature = c(21, 22),
% #'       calibration.type = "Uk37",
% #'      point.or.sample = "point"
% #'     )
% #'   )
% #'   ,
% #'   n = 5,
% #'   point.or.sample = "sample",
% #'   calibration.type = "Uk37"
% #' )
% #' 
% #' ## Incompatible arguments
% #' \dontrun{
% #' ProxyConversion(temperature = 1, proxy.value = 1)
% #' }

if strcmp(calibration_type,'identity')
    calibration = NaN;
elseif strcmp(calibration_type,'MgCa')
    calibration = 'Ten planktonic species_350-500';
elseif strcmp(calibration_type,'Uk37')
    calibration = 'Mueller global';
end

if (isempty(temperature) && isempty(proxy_value)) || ...
      (~isempty(temperature) & ~isempty(proxy_value))
    disp('One and only one of temperature or proxy.value must be supplied')
end

% ## Check dimensions if matrix
if (min(size(temperature))>1 || min(size(proxy_value)>1))
    if (point_or_sample == "sample" && max(size(temperature,2), size(proxy_value,2)) ~= n)
        stop('If input is matrix and sample, n must equal ncol(input)')
    end
end

if (strcmp(point_or_sample,"point") && n > 1)
    stop('Multiple replicates only returned if sample');
end

T = readtable('sedproxy/calibration_parameters.csv');

% ## Get calibration parameters
if ~strcmp(calibration_type,"identity")

    idx = find(strcmp(T.calibration_type,calibration_type) & strcmp(T.calibration,calibration));

    if ~isempty(slp_int_means)
        cfs = slp_int_means;
    else
        cfs = [T.slope(idx) T.intercept(idx)];
    end

    if ~isempty(slp_int_vcov)
        vcov = slp_int_vcov;
    else
        V = split(T.vcov(idx),'|');
        vcov(1,1) = str2num(cell2mat(V(1))); vcov(2,2) = str2num(cell2mat(V(4)));
        vcov(1,2) = str2num(cell2mat(V(2))); vcov(2,1) = str2num(cell2mat(V(3)));
    end

    if strcmp(point_or_sample,'sample')
        if ~isempty(slp_int_means) & isempty(slp_int_vcov)
        warning(['Sampling calibration parameters using user supplied values' ,...
                ' for the mean slope and intercept but the variance covariance matrix for the ',...
                'default or named calibration.'])

        cfs = mvnrnd(cfs, vcov, n);
        end
    end

    % Do converstion

  %## check if vector input and convert to 1 column matrix
  is_vec = sum([isvector(temperature), isvector(proxy_value)])>0;

  if (is_vec) % this might be wrong, idk?
    if (isempty(temperature))
      proxy_value = repmat(proxy_value,n);
    elseif (isempty(proxy_value))
      temperature = repmat(temperature,n);
    end
  end

  if strcmp(calibration_type,'identity')
      if (isempty(temperature))
          out = proxy_value;
      else
          out = temperature;
      end
  elseif strcmp(calibration_type,'MgCa')
      cfs(:,2) = exp(cfs(:,2));

      % convert from temperature to MgCa
      if (isempty(proxy_value))
          out = (cfs(:,2).*exp(temperature'*cfs(:,1)));
      else % convert from MgCa to temperature
          out = log(proxy_value' ./ cfs(:,2)) ./ cfs(:,1);
      end

  elseif strcmp(calibration_type,'Uk37')
      % convert from temperature to UK'37
      if (isempty(proxy_value))
          out = cfs(:, 2) + temperature' * cfs(:, 1);
      elseif (isempty(temperature))
          % convert from UK'37 to temperature
          out = (proxy_value' - cfs(:, 2)) / cfs(:, 1);
      end
  end
end

if isempty(temperature)
    if size(out)~=size(proxy_value)
        out = out';
    end
else
    if size(out)~=size(temperature)
        out = out';
    end
end


% convert back to vector if vector input and n == 1
if (is_vec && n == 1)
    out = out(:);
end

% # # Visual checks of calibration data and parameters
% #
% # ## Uk37
% #
% # df <- tibble(x = 1:30,
% #              y = ProxyConversion(temperature = x, point.or.sample = "point",
% #                                  calibration.type = "Uk37"))
% #
% # dat <- climproxycalibration::mueller.uk37.sst
% #
% # dat %>%
% #   ggplot(aes(x = `SST (1-12) [°C]`, y = `UK'37`)) +
% #   geom_point( ) +
% #   geom_line(data = df, aes(x=x, y =y))
% #
% # t1 <- 1:30
% #
% # p1 <- ProxyConversion(temperature = t1, point.or.sample = "point",
% #                       calibration.type = "Uk37")
% #
% # t2 <- ProxyConversion(proxy.value = p1, point.or.sample = "point",
% #                       calibration.type = "Uk37")
% #
% # all.equal(t1, t2)
% #
% #
% # ## MgCa
% #
% # df <- tibble(x = 1:30,
% #              y = ProxyConversion(temperature = x, point.or.sample = "point",
% #                                  calibration.type = "MgCa"))
% #
% # dat <- climproxycalibration::anand.2003.hydro
% #
% # dat %>%
% #   ggplot(aes(x = calc.t, y = `Mg/Ca (mmol/mol)`)) +
% #   geom_point( ) +
% #   geom_line(data = df, aes(x=x, y =y))
% #
% # t1 <- 1:30
% #
% # p1 <- ProxyConversion(temperature = t1, point.or.sample = "point",
% #                       calibration.type = "MgCa")
% #
% # t2 <- ProxyConversion(proxy.value = p1, point.or.sample = "point",
% #                       calibration.type = "MgCa")
% #
% # all.equal(t1, t2)
end

