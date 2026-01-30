function folderList = listFolders(directory, pattern)
    % listFolders - Returns a cell array of all folders in the given directory that match a regex pattern
    %
    % Syntax:
    %   folderList = listFolders(directory, pattern)
    %
    % Input:
    %   directory - The path to the directory to search
    %   pattern   - Regular expression pattern to match folder names
    %
    % Output:
    %   folderList - Cell array containing the relative paths of matching subfolders
    %
    % Example:
    %   folders = listFolders('C:\Users\Username\Documents', '^data.*')
    
    if nargin < 1
        directory = pwd; % Default to current directory if none specified
    end
    if nargin < 2
        pattern = '.*'; % Default to match all folders
    end
    
    contents = dir(directory);
    isFolder = [contents.isdir]; % Logical array indicating directories
    folderNames = {contents(isFolder).name};
    
    % Remove '.' and '..' which represent current and parent directories
    folderNames = folderNames(~ismember(folderNames, {'.', '..'}));
    
    % Apply regex filter
    matches = cellfun(@(name) ~isempty(regexp(name, pattern, 'once')), folderNames);
    folderNames = folderNames(matches);
    
    % Construct relative paths
    folderList = fullfile(directory, folderNames);
end
