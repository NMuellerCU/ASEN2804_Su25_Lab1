function [ThrustCurves, Time, ThrustStruct] = Thrust(Plot_Thrust_Data)
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
% statsStruct
%   A structure for holding statistics data. This will be broken into 2L
%   and 1.25L, with cases undeaneath and finally values of interest at the
%   bottom

%% Define data locations
fileLoc_2L = 'Model Functions/Thrust Model/Static Test Stand Data/2000mL Bottle/'; % path to the data files, be sure to include a trailing slash
fileLoc_1pt25L = 'Model Functions/Thrust Model/Static Test Stand Data/1250mL Bottle/'; % path to the data files, be sure to include a trailing slash

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
    configSTD = zeros(1, 750); % re-zero for each configuration

    TmaxMean = 0;
    TmaxSTD = 0;
    timeMean = 0;
    timeSTD = 0;
    impulseMean = 0;
    impulseSTD = 0;

    for i = 1:numTests
        %% Load data
        fileName = testNames(i, :); % again weird indexing is due to string arrays, we have to ask for all the characters in a row
        data = readmatrix(fileName); % load the data
        data = data(:,3)*4.448; % take only the third column and converting from lbf to N

        %% Data Conditioning
        % Find ending offset
        offset = data(mean(end-100:end));
        if offset < 50 && offset > 0 % Check that the offset is reasonable, if not just skip the data set by not entering the loop
            % Force negative data to be 0 since that is unphysical
            data(data <= 0) = 0;

            data = data(50:end); % cut off the front to avoid any startup errors
            % data(isnan(data)) = 0; % also set nan to 0

            % Find the index of max thrust
            [~,iMax] = max(data);

            % Find start of test - defined by 5 samples before last sample below 10N prior to max
            iStart = find(data(1:iMax) < 10, 1, 'last');
            data = data(iStart-5:end); % Truncate data to begin at the start of the test - note that this means our max index is now wrong

            % Trim the data to about 0.5 seconds
            data = data(1:750);

            % Find ending offset
            offset = mean(data(end-100:end));

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

            %% Statistics and averaging of specific quantities
            configMean = configMean + thrustOut;
            configSTD = configSTD + thrustOut.^2;
            % We need to do this for individual vaiables becuase max(mean)
            % is NOT equal to mean(max) in this case (it would be if all of
            % the max's happened at the same point in time)
            Tmax = max(thrustOut);
            TmaxMean = TmaxMean + Tmax; 
            TmaxSTD = TmaxSTD + Tmax^2;
            timeMean = timeMean + testTime(end);
            timeSTD = timeSTD + testTime(end)^2;
            impulse = trapz(tMean, thrustOut);
            impulseMean = impulseMean + impulse;
            impulseSTD = impulseSTD + impulse^2;

        else
            numDiscard = numDiscard + 1;
        end
    end
    
    configMean = configMean/(numTests-numDiscard); % This is the average of all of the tests
    configSTD = configSTD/(numTests-numDiscard); % Need the mean of the squares of the values for...
    configSTD = sqrt(configSTD - configMean.^2); % std = sqrt[mean(x^2) - mean(x)^2] where x is some data

    TmaxMean = TmaxMean/(numTests-numDiscard);
    TmaxSTD = TmaxSTD/(numTests-numDiscard);
    TmaxSTD = sqrt(TmaxSTD - TmaxMean^2);
    timeMean = timeMean/(numTests-numDiscard);
    timeSTD = timeSTD/(numTests-numDiscard);
    timeSTD = sqrt(timeSTD - timeMean^2);
    impulseMean = impulseMean/(numTests-numDiscard);
    impulseSTD = impulseSTD/(numTests-numDiscard);
    impulseSTD = sqrt(impulseSTD - impulseMean^2);

    %% Pack up statistics of interest
    BottleName = ['B', bottleSize];
    WaterName = ['W', num2str(waterSize)];
    ThrustStruct.(BottleName).(WaterName).Wwater = waterSize/1000*9.8;
    ThrustStruct.(BottleName).(WaterName).MaxThrust = TmaxMean;
    ThrustStruct.(BottleName).(WaterName).stdMaxThrust = TmaxSTD;
    ThrustStruct.(BottleName).(WaterName).TestTime = timeMean;
    ThrustStruct.(BottleName).(WaterName).stdTestTime = timeSTD;
    ThrustStruct.(BottleName).(WaterName).Impulse = impulseMean;
    ThrustStruct.(BottleName).(WaterName).stdImpulse = impulseSTD;
    ThrustStruct.(BottleName).(WaterName).Isp = impulseMean/(waterSize/1000*9.8);
    ThrustStruct.(BottleName).(WaterName).MeanProfile = configMean;
    ThrustStruct.(BottleName).(WaterName).stdMeanProfile = configSTD;

    %% Finding quantities of interest on the mean data
    % Find the max of the average
    [~,iMax] = max(configMean);

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
    dataFit2(dataFit2 < 0) = 0; % force data to be non-negative

    %% Sample onto the standard output array format
    % Make an array to make thrust = 0 for the rest of the time until 0.5 s
    dataEnd = zeros(1, 500 - ceil(tMean(end)*1000));

    % Make the fitted trust array
    thrustOut = [dataFit1, dataFit2, dataEnd];

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

%% Plots for this function (Figure 800 - 899)
if Plot_Thrust_Data == 1

    % Static Test Thrust Profle Model Plots
    cmap = colormap(lines(12));
    set(0,'DefaultAxesColorOrder',cmap)
    set(gca,'ColorOrder',cmap);
    
    %2L Thrust Profiles
    figure(800)
    for i = 1:numConfigs_2L
        plot(Time, ThrustCurves{:, i}, DisplayName = [ThrustCurves.Properties.VariableNames{i}(6:end), ' mL'])
        if i == 1
            hold on
        end
    end
    xlabel('Test Time [s]');
    ylabel('Fitted Thrust [N]');
    title('2L Fitted Thrust Profiles');
    legend();
    grid on
    hold off
    
    %1.25L Thrust Profiles
    figure(801)
    for i = (numConfigs_2L+1):width(ThrustCurves)
        plot(Time, ThrustCurves{:, i}, DisplayName = [ThrustCurves.Properties.VariableNames{i}(6:end), ' mL'])
        if i == (numConfigs_2L+1)
            hold on
        end
    end
    xlabel('Test Time [s]');
    ylabel('Fitted Thrust [N]');
    title('1.25L Fitted Thrust Profiles');
    legend();
    grid on
    hold off

    %Reset default color order
    set(0,'DefaultAxesColorOrder','default')

end

end
