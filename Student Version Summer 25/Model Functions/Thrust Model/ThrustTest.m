clear; clc; close all;

%% Thrust Summary
% This funciton will take in the file location of the two test setups and
% using the file names in those directories, will pull out all of the
% avialable tests, process them, and fit the data into a standard
% formatting for output

%% Outputs:
% ThrustCurves:
%   A table containing 0.5 seconds of thrust data for each of the cases
%   available, this data will have formatting such that there are 501
%   evenly spaced thrust data points (rows) for each test (columns). The
%   ordering of the columns will go from max to min water volume in the 2L
%   bottle and then max to min in the 1.25L bottle
%
% Time:
%   A 1D array corresponding to the times of the thrust data points in the
%   ThrustCurves table
%
% <User defined variables for statistics>

%% Define data locations
fileLoc_2L = 'Static Test Stand Data/2000mL Bottle/Variable Volume/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = 'Static Test Stand Data/1250mL Bottle/'; % path to the data files, be sure to include a trailing slash

%% Read in all of the avilable data and find what data there is
testInfo_2L = getThrustTestNames(fileLoc_2L);
configs_2L = unique(testInfo_2L.waterVol);
numConfigs_2L = length(configs_2L);

testInfo_1pt25L = getThrustTestNames(fileLoc_1pt25L);
configs_1pt25L = unique(testInfo_1pt25L.waterVol);
numConfigs_1pt25L = length(configs_1pt25L);

numConfigs = numConfigs_2L + numConfigs_1pt25L;

% Set known sampling frequency
f= 1652; % [Hz]

%% Preallocate variables of interest
Time = 0:0.001:0.5; % just go ahead and define this, note that it will be 501 long
ThrustCurves = zeros(length(Time),numConfigs);

ThrustCurvesNames = {};

%% Loop over all of the configurations
for N = 1:numConfigs % use upper case N to distiguish that it is counting something different from the aerodynamic modeling loops

    %% Dertemine what configuration to use for this iteration in the loop
    if N <=  numConfigs_2L % determine if we should be reading 2L or 1.25L data
        bottleSize = '2000'; % [ml]
        waterSize = configs_2L(N);
        testIndexes = find(testInfo_2L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_2L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    else
        bottleSize = '1250'; % [ml]
        waterSize = configs_1pt25L(N-numConfigs_2L);
        testIndexes = find(testInfo_1pt25L.waterVol == waterSize); % finds the index of the relavant tests
        numTests = length(testIndexes); % finds the number of tests performed
        testNames = testInfo_1pt25L.fileNames(testIndexes, :); % pulls all of the test names of interest, weird indexing is due to string arrays
    end

    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    numDiscard = 0; % make a counter to count how mant data sets we throw away for a given configuration
    configMean = zeros(1, 750); % re-zero for each configuration
    close all
    for i = 1:numTests
        %% Load data
        fileName = testNames(i, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
        data = readmatrix(fileName); % load the data
        data = data(:,3)*4.448; % take only the third column and converting from lbf to N

        % get the sampling frequency
        % f = readlines(fileName);
        % f = f(1);
        % f = str2double(extractBetween(f,"Sampling Frequency ", " kHz"))*1000;  % [Hz]
        %
        % figure(1)
        % plot(data)
        % hold on

        %% Data Conditioning
        % Find ending offset
        offset = data(mean(end-100:end));
        if offset < 50 && offset > 0 % Check that the offset is reasonable, if not just skip the data set by not entering the loop
            % Force negative data to be 0 since that is unphysical
            data(data <= 0) = 0;

            data = data(100:end); % cut off the front to avaid any startup errors
            % data(isnan(data)) = 0; % also set nan to 0
            % data(data > 500) = 0; % also remove super high values

            % Find the index of max thrust
            [~,iMax] = max(data);

            % Find start of test - defined by 5 samples before last sample below 10N prior to max
            iStart = find(data(1:iMax) < 10, 1, 'last');
            data = data(iStart-5:end); % Truncate data to begin at the start of the test - note that this means our max index is now wrong

            % Trim the data to about 0.5 seconds
            data = data(1:750);

            % Find ending offset
            offset = data(mean(end-100:end));

            % Find end of test - defined by last time the absolute value of the
            % derivative is greater than 1 (arbitrary choice)
            iEnd = find(abs(diff(data)) >= 1, 1, 'last');
            data = data(1:iEnd);

            % Find total test time
            testTime = linspace(0,length(data)/f,length(data))';

            % Apply linear offset that implies constant mass loss
            data = data - offset*linspace(0, 1, length(data))';

            % Configure for averaging
            thrustOut = zeros(1, 750);
            tMean = 0:1/f:1/f*749;
            thrustOut(1:length(data)) = data;
            configMean = configMean + thrustOut;

            figure(10)
            plot(tMean, thrustOut)
            hold on

        else
            numDiscard = numDiscard + 1;
        end
    end
    configMean = configMean/(numTests-numDiscard); % This is the average of all of the tests

    figure(10)
    plot(tMean, configMean, 'k')
    hold on

    % Find the max of the average
    [Tmax,iMax] = max(configMean);

    % Find end of the mean - again defined by last time the absolute value of the
    % derivative is greater than 1 (arbitrary choice)
    iEnd = find(abs(diff(configMean)) >= 1, 1, 'last');
    configMean = configMean(1:iEnd);
    tMean = tMean(1:iEnd);


    %% Data Fitting
    % Fit thrust rise (start --> max)
    coeff1 = polyfit(tMean(1:iMax),configMean(1:iMax),3);
    timeFit1 = 0:0.001:ceil(tMean(iMax)*1000)/1000;  %time vector to sample onto
    dataFit1 = polyval(coeff1,timeFit1);
    dataFit1(dataFit1 < 0) = 0; % force data to be non-nagative

    % Fit thrust fall off (max --> end),
    coeff2 = polyfit(tMean(iMax:end),configMean(iMax:end),3);
    timeFit2 = timeFit1(end)+0.001:0.001:ceil(tMean(end)*1000)/1000; % time vector; each element is space by 1 ms
    dataFit2 = polyval(coeff2,timeFit2);
    dataFit2(dataFit2 < 0) = 0; % force data to be non- negative


    %% Sample onto the standard output array format
    % Make an array to make thrust = 0 for the rest of the time until 0.5 s
    dataEnd = zeros(1, 500 - ceil(tMean(end)*1000));

    figure(5)
    plot(tMean, configMean, 'k')
    hold on
    % Make the fitted trust array
    thrustOut = [dataFit1, dataFit2, dataEnd];

    % thrustOut = interp1(testTime, data, Time, 'linear', 0);
    figure(5)
    plot(Time, thrustOut)

% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
    %% Convert to table for output
    % It is very important that the data is 501 elements long corresponding
    % to 0-0.5 seconds of time at this point!!!
    ThrustCurves(:, N) = thrustOut;
    % Header naming convention of <bottle size (in ml)>_<water volume (in ml)>
    ThrustCurvesNames{N} = [bottleSize, '_', num2str(waterSize)];
end
ThrustCurves = array2table(ThrustCurves);
ThrustCurves.Properties.VariableNames = ThrustCurvesNames;

