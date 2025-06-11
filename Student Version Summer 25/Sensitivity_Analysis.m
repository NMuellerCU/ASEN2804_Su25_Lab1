%% Clean Workspace and Housekeeping
clear
% clearvars
close all

% removes warnings for table variable names for a cleaner output
warning('OFF', 'MATLAB:table:ModifiedAndSavedVarnames')

%Add folder and subfolder path for standard Design Input Files, Model
%Functions, and Static Test Stand Data
addpath(genpath('Design Input Files'));
addpath(genpath('Model Functions'));

%% Import and Read Aircraft Design File
Design_Inputs = {
    % [File, Column Being Changed, Output Name, Title, xlabel, ylabel]
    % ["Design Input File_V25-00-SUMMER.xlsx", "Dia_f", "bestGlide", "Diameter Impact on Glide Performance", "Diameter [m]", "Glide Range [m]"];
};

len = length(Design_Inputs);
for i = 1:len
    current_input = Design_Inputs{i};
    [Design_Input, GlideData] = ASEN2804_Aerospace_Design_Model_Main( current_input(1) );

    eval( sprintf("param = Design_Input.%s;", current_input(2)) );
    eval( sprintf("output = GlideData.%s;", current_input(3)) );

    figure();
    hold on;
    plot(param, GlideData.bestGlide);
    title( current_input(4) )
    xlabel( current_input(5) );
    ylabel( current_input(6) );
end


