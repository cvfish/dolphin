function build_all(varargin)
%BUILD_ALL Builds all mex files (or those in a specific list)
%
%   BUILD_ALL
%
%       Builds all mex files in Dolphin
%
%   BUILD_ALL name1 name2 ...
%
%       Builds a specific list of targets
%
%   BUILD_ALL clean
%   BUILD_ALL clean name1 name2 ...
%
%       Clears all or a specific list of targets
%

%% extract list

targets = dolphin.setup.load_build_list();
do_clean = false;

if ~isempty(varargin)  
    if ~iscellstr(varargin)
        error('build_all:invalidarg', ...
            'The arguments should be a list of strings.');
    end
    
    if strcmpi(varargin{1}, 'clean')
        do_clean = true;
        tarlist = varargin(2:end);
    else
        tarlist = varargin;
    end
else
    tarlist = [];
end
   
if ~isempty(tarlist)
    
    [tf, loc] = ismember(tarlist, {targets.name});
    if ~all(tf)
        i = find(~tf, 1);
        error('build_all:invalidarg', ...
            'Unrecognized target: %s', tarlist{i});
    end
    targets = targets(loc);
end

%% building

rootdir = fileparts(fileparts(mfilename('fullpath')));

for i = 1 : length(targets)
    t = targets(i);
    srcpath = fullfile(rootdir, [t.name '.cpp']);
    mexpath = fullfile(rootdir, [t.name '.' mexext]); 
    
    if ~do_clean % Build
           
        if ~exist(srcpath, 'file')
            error('build_all:filenotfound', ...
                'The source file %s is not found.', srcpath);
        end
        
        % determine whether the mex is updated
        
        updated = false;
        
        if exist(mexpath, 'file')
            info_src = dir(srcpath);
            info_mex = dir(mexpath);
            if time_lt(info_src.date, info_mex.date)
                updated = true;
            end
        end
        
        if updated
            fprintf('[TARGET] %s is updated (last modified: %s)\n', ...
                t.name, info_mex.date);
            
            continue;
        end
        
        % build mex
        
        fprintf('[TARGET] %s building ...\n', t.name);
        dolphin.setup.build_mex(srcpath, t.deps, '-O');
    
    else % Clean
       
        if exist(mexpath, 'file')
            delete(mexpath);
        end                    
    end
end


%% auxiliary function

function b = time_lt(a, b)

a = datevec(a);
b = datevec(b);

assert(length(a) == 6 && length(b) == 6)

for i = 1 : 6
    if a(i) < b(i)
        b = true;
        return
    elseif a(i) > b(i)
        b = false;
        return
    end    
end

b = false;






