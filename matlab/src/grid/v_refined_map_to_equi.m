function [X,V] = v_refined_map_to_equi(X,V,original_indizes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function gets a unstructured refined grid that was
% computed with refine_velocity_grid and returns
% the equidistant base grid
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X = X(original_indizes, :);
V = V(original_indizes, :);

end