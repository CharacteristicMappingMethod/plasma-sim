%% Compilation script for B-spline MEX functions
% This script compiles the C++ MEX files for B-spline interpolation

fprintf('Compiling B-spline MEX functions...\n');

% Get the current directory
current_dir = pwd;
mex_dir = fileparts(mfilename('fullpath'));

% Change to the MEX directory
cd(mex_dir);

try
    % Compile the B-spline evaluation MEX file
    mex('bspline_periodic_eval_2d_mex.cpp', '-output', 'bspline_periodic_eval_2d_mex_cpp');
    fprintf('✓ Successfully compiled bspline_periodic_eval_2d_mex_cpp\n');

catch ME
    fprintf('✗ Compilation failed: %s\n', ME.message);
    fprintf('Make sure you have a C++ compiler configured for MATLAB MEX.\n');
    fprintf('You can configure a compiler using: mex -setup C++\n');
end

% Return to original directory
cd(current_dir);

fprintf('Compilation complete.\n');