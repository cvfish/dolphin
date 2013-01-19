function targets = load_build_list()
%LOAD_BUILD_LIST Loads the mex-build list for Dolphin
%
%   targets = LOAD_BUILD_LIST;
%

fdir = fileparts(mfilename('fullpath'));
listfile = fullfile(fdir, 'build_list.txt');

fid = fopen(listfile, 'r');
if fid < 0
    error('load_build_list:ioerror', 'Failed to open file: %s', listfile);
end
lines = textscan(fid, '%s', 'Delimiter', '\r\n');
fclose(fid);
lines = lines{1};

targets = cell(length(lines), 1);
n = 0;

for i = 1 : length(lines)
    line = strtrim(lines{i});
    if isempty(line) || line(1) == '#'
        continue
    end
    
    % parse line
    
    icolon = strfind(line, ':');
    if isempty(icolon)
        name = line;
        deps = [];
    else
        name = strtrim(line(1:icolon-1));
        deps = strtrim(line(icolon+1:end));
        
        if ~isempty(deps)
            deps = regexp(deps, '\s*,\s*', 'split');            
        else
            deps = [];
        end
    end
        
    n = n + 1;
    tar = [];
    tar.name = name;
    tar.deps = deps;
    targets{n} = tar;
    
end

targets = vertcat(targets{1:n});
