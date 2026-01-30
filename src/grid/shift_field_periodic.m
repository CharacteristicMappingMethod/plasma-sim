function [f_shifted] = shift_field_periodic(f, shift)
% Shift a 2D data field periodically in x direction
% f: 2D array (Nv x Nx), shift: integer shift amount
% Positive shift = right, negative shift = left

    shift = round(shift);
    [~, Nx] = size(f);
    shift = mod(shift, Nx);
    f_shifted = circshift(f, [0, shift]);
    
end
