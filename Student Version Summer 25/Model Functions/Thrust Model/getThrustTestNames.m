function [TestInfo] = getThrustTestNames(fileLoc)
%% Get Thrust Test Names Summary:
% This function is meant to return a table of file paths for each file in
% the test directory and the volume of water corresponding to that test.
% This function is not meant to be modified by the students and is simply
% meant to organize the data for their usage

list = dir([fileLoc, '*Test*']); % This lists all files in fileLoc with 'Test" in the file name
numFiles = length(list);
waterVol = zeros(numFiles, 1);
for i = 1:numFiles
    fileNames{i} = [fileLoc, list(i).name]; % This makes a string of the complete file name with the path in front of it
    waterVol(i) = str2num(char(extractBetween(list(i).name, '_W', '_B'))); % Finds the water volume in the name and make it a usable number
end
fileNames = char(fileNames); % Convert to char array from cell array
TestInfo = table(waterVol, fileNames); % make a table :)
end

