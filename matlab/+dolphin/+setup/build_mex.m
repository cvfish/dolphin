function build_mex(filename, deps, varargin)
%BUILD_MEX Builds mex files in Dolphin Library
%
%   BUILD_MEX(filename);
%
%       Build a file in default setting
%
%   BUILD_MEX(filename, deps);
%
%       Build a file with dependent libraries. Here, deps is a cell
%       array of dependency names (e.g. 'svml', 'mwblas', etc)
%
%   BUILD_MEX(filename, deps, ...);
%
%       Additionally specifies other options. If a file has no extra
%       library dependency, deps can be set to an empty array.
%

%% verify arguments

if ~ischar(filename)
    error('build_mex:invalidarg', 'filename must be a string.');
end

if nargin < 2 || isempty(deps)
    has_deps = false;
else
    if ~iscellstr(deps)
        error('build_mex:invalidarg', ...
            'deps should be a cell array of strings.');
    end
    has_deps = true;
end

%% main

lmat_inc = getenv('LIGHT_MATRIX_HOME');
if isempty(lmat_inc)
    error('build_mex:fatal_error', ...
        'The environment variable LIGHT_MATRIX_HOME is not set.');
end

dolp_inc = fileparts(fileparts(fileparts(fileparts(mfilename('fullpath')))));

outdir = fileparts(filename);
if isempty(outdir)
    outdir = '.';
end

if has_deps
    dep_opts = cellfun(@(s) ['-l', s], deps, 'UniformOutput', false);
else
    dep_opts = {};
end

opts = [{['-I' lmat_inc], ['-I', dolp_inc], '-outdir', outdir, '-DLMAT_USE_INTEL_SVML'}, ...
    dep_opts, varargin];

mex(filename, opts{:});
    