function [out] = ForamGrowthfT(foram,temperature_K, varargin)
% #' @title Foraminifer Growth Rate Function from Lombard et al. (2009)
% #' @description Implements the function for foraminifer growth rate as a function of
% #' temperature from Lombard et al. (2009) with parametrization from FAME 1.0
% #' (Roche et al, 2018).
% #'
% #' @param foram Name of foram.
% #' @param temperature_K Temperature in Kelvin
% #' @param norm Optional normalizing factor
% #' @param min.growth.thresh Sets a lower cutoff for growth as a proportion of
% #'   the maximum growth rate for that taxon. For example in Roche et al (2018) a cutoff of
% #'   0.1 was used, meaning all growth rates less than 0.1*max were set to zero.
% #' @references
% #'
% #' Lombard, F., Labeyrie, L., Michel, E., Spero, H. J., and Lea, D. W.:
% #' Modelling the temperature dependent growth rates of planktic foraminifera,
% #' Marine Micropaleontology, 70, 1–7,
% #' https://doi.org/10.1016/j.marmicro.2008.09.004, 2009.
% #'
% #' Roche, D. M., Waelbroeck, C., Metcalfe, B. and Caley, T.: FAME
% #'   (v1.0): a simple module to simulate the effect of planktonic foraminifer
% #'   species-specific habitat on their oxygen isotopic content, Geosci. Model
% #'   Dev. Discuss., 2017, 1–22, doi:10.5194/gmd-2017-251, 2017.
% #'
% #' @return A numerical vector or matrix with the same dimensions as the object
% #'   passed to temperature_K. Units are daily growth rate, unless norm == TRUE.
% #' @export
% #'
% #' @examples
% #' ForamGrowthfT(foram = 'ruber', temperature_K = (c(280, 290)), norm = 1)

if isempty(varargin)
    norm = 0; min_growth_thresh = 0;
elseif length(varargin)<2
    norm = varargin{1};
    min_growth_thresh = 0;
else
    norm = varargin{1};
    min_growth_thresh = varargin{2};
end

if ~ismember(foram,{'sacculifer', 'bulloides', 'pachy_d',...
        'siphonifera', 'universa', 'pachy_s',...
        'dutertrei', 'ruber'})
    disp('warning! try another value for foram')
end

idx = find(strcmp(foram,{'sacculifer', 'bulloides', 'pachy_d',...
        'siphonifera', 'universa', 'pachy_s',...
        'dutertrei', 'ruber'})==1);

T = readtable('sedproxy/l09_cnsts_dic');
pars = table2array(T(:,idx));
muT1 = pars(1); TA = pars(2); TL = pars(3); TH = pars(4); TAL = pars(5); TAH = pars(6);

out = muT1 * exp(TA ./ 293 - TA ./ temperature_K) ./ ...
        (1 + exp(TAL ./ temperature_K - TAL ./ TL) + exp(TAH ./ TH - TAH ./ temperature_K));

T = readtable('sedproxy/l09_maxgrowth_dic.csv'); T = table2array(T);
max_func = T(idx);

out(out < min_growth_thresh.* max_func) = 0;

if (norm == 1)  
    out = out ./ max_func;
end

out(temperature_K < -2 + 273.15) = 0;

end

