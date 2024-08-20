function [mm] = ChunkMatrix(timepoints, width, climate_matrix,time_in)

  rel_wind = 1:width -ceil(width/2);

  strt = time_in(1);
  n_row = size(climate_matrix,1);

  for i = 1:length(timepoints)
    inds = rel_wind + timepoints(i) - strt + 1;
    inds = inds(inds>0 & inds < n_row);
    if size(climate_matrix,2) >1
        m = climate_matrix(inds,:);
    else
    m = climate_matrix(inds);
    end
    mm(i) = mean(m(:));
  end
end

