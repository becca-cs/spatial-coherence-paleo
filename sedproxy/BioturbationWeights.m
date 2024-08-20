function [fz] = BioturbationWeights(z, focal_z, layer_width, sed_acc_rate, bio_depth, scale)
%BIOTURBATIONWEIGHTS Summary of this function goes here
%   scale should be "depth" or "time
% layer_width should be 1

if ~isempty(layer_width)
    layer_width = 1;
end

sed_acc_rate = sed_acc_rate/1000; % in cm / yr

if strcmp(scale,"depth")
    z = z / sed_acc_rate;
    focal_z = focal_z / sed_acc_rate;
end


lwy = ceil(layer_width / sed_acc_rate);
mdy = ceil(bio_depth / sed_acc_rate);

if (lwy == 0 && mdy == 0) 
    lwy = 1;
end

C = lwy/2;
lam = 1/mdy;

z = z - focal_z + mdy;

if (mdy <= 1 && lwy > 0)
    fz = unifpdf(z,-C,C);
elseif (lwy == 0)
    fz = exppdf(x,mdy);
else
    fz = (z < -C) * 0 + ...
      (z >= -C & z <= C).* (lam*(1/lam-exp(-lam*C-lam*z)/lam))/(2*C)  + ...
      (z > C).* (lam*(exp(lam*C-lam*z)/lam-exp(-lam*C-lam*z)/lam))/(2*C);
end
  
if (sum(fz) == 0)
    fz = fz;
else
    fz = fz / sum(fz);
end

end

