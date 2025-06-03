function [apogee, hApogee, stateStruct] = BoostAscent(Design_Input, ATMOS, Parasite_Drag_Data, Weight_Data, ThrustCurves, Time, Count, g,Plot_Boost_Data)
%% BoostAscent Summary:
% The main purpose of this function is to set up for the call to ODE 45 to
% propigate our eqations of motion, and then to format the output of ODE 45
% to be passed back into main for the glide function and for plotting.

%% Outputs:
%
% apogee:
%   A vector of values of the maximum height achieved in [m] for each case 
%   input
%
% hApogee:
%   A table of inertial heading vectors in cartisian coordinates for each
%   case input
%
% stateStruct:
%   A structure containing the full output data from ODE45 for each case
%   input. This structure's second level will correspond to the different
%   input cases, with the data being containied within being the table of
%   states for all time steps. This is meant only for in-depth analysis and
%   is not meant to be munipulated by the students

%% Preallocate variables of interest
apogee = zeros(Count,1); % Height of apogee above the ground (this will be postitive despite the coordinate system) [m]
hApogee = zeros(Count,3); % Body heading/pointing unit vector at apogee

%% Preallocate ODE varaibles
consts = zeros(10, 1);
S0 = zeros(7, 1);

%% Some useful constants
rho_w = 1000; % Density of water [kg/m^3]
mu_k = 0.2; % Launch rail coefficient of dynamic friction []
A_exit = 0.021; % Area of bottle outlet [m^2]

%% Loop through different configurations
for n = 1:Count
    %% Pick the correct thrust curve
    waterSize = Design_Input.Water_Vol(n);
    bottleSize = Design_Input.Bottle_Vol(n)*1000; % Convert to [ml]

    thrustCurveName = [num2str(bottleSize), '_', num2str(waterSize)];

    thrustVec = ThrustCurves.(thrustCurveName); % use parenthasis to index into a table with a variable
    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    %% Make variables for the constants vector
    m_empty = Weight_Data.Wo(n)/g; % [kg]
    m0 = m_empty + waterSize/1000; % Note that the input water volume should be in ml which approx = grams
    % /////////////////////////////////////////////////////////////////////////
    % MODIFY THIS SECTION
    % /////////////////////////////////////////////////////////////////////////
    %% Pack Constants Vector
    % Basic Properties
    consts(1) = g; % Accelatation due to gravity [m/s^2]
    consts(2) = rho_w; % Density of water [kg/m^3]
    consts(3) = ATMOS.rho(n); % Density of air [kg/m^3]
    consts(4) = mu_k; % Launch rail coefficient of dynamic friction []
    % Vehicle info
    consts(5) = A_exit; % Area of bottle outlet [m^2]
    consts(6) = Parasite_Drag_Data.CDo(n); % C_D of the vehicle (assume zero lift) []
    consts(7) = Design_Input.Sref_w(n); % Wing reference area [m^2]
    consts(8) = m_empty; % Weight of the rocet with no water [kg]
    % % Launch Direction
    consts(9) = Design_Input.Launch_El(n); % Launch elevation [degrees]
    consts(10) = Design_Input.Launch_Az(n); % Launch Azimuth [degrees], measured CW from north when looking down on the map; also known as compass heading

    %% Make an initial condition state vector (S0)
    S0(1) = 0; % inertial velocity in x-direction [m/s]
    S0(2) = 0; % inertial velocity in y-direction [m/s]
    S0(3) = 0; % inertial velocity in z-direction [m/s]
    S0(4) = 0; % position in x (inertial) [m]
    S0(5) = 0; % position in y (inertial) [m]
    S0(6) = 0; % position in z (inertial) [m]
    S0(7) = m0; % current total mass [kg]

    %% Final setup
    % set and event to stop the integration 
    opts = odeset('Events', @stoppingPoint, 'RelTol',1e-8);

    % Set a start and max time to integrate over (pick one of the below)
    %intSpan = [0, 8]; % [s] - Lets matlab pick its timestep (reccomended)
    intSpan = 0:0.001:5; % [s] - forces output to have 1ms spacing

    %% Call to ODE 45
    % The @ part is telling ODE45 what our independent vaiable is (t) and
    % what the dependent variables to propigate in time are (S)
    %
    % Then, we pass it the odefun we made that can take in the current
    % time, state, and whatever else it requires and retuns the derivative
    % of the state variables at the current time
    %
    % Finally, we also give ODE45 a span of time to integrate over
    % (intSpan), an initial conditions vector (S0), and any options we set
    % (opts), which in this case holds the function to stop the integration
    [t, S] = ode45(@(t, S) BoostAscent_odefun(t, S, consts, thrustVec, Time), intSpan, S0, opts); %S = State Vector (Vx,Vy,Vz,X,Y,Z,Mass); t = Time Vector (sec)

    %% Find outputs of interest
    [apogee(n), iApogee] = max(abs(S(:, 6)));
    hApogee(n, :) = S(iApogee, 4:6)/norm(S(iApogee, 4:6));

    %% Store full output
    configName = ['Config_' num2str(n)];
    stateStruct.(configName).time = t; %Struct Column Variables: time in sec
    stateStruct.(configName).data = S; %Struct Columns (1-7): Vx,Vy,Vz,X,Y,Z,mass


end
%% Convert to tables for output
dirNames = {'x', 'y', 'z'}; %Create Inertial Axis Labels for table
hApogee = array2table(hApogee); % Convert to table
hApogee.Properties.VariableNames = dirNames;

%% Plots for this function (Figure 900 - 999)
if Plot_Boost_Data == 1

    % Boost_Ascent Flight Profile Plots
    % Some setup to make plots more readable in color, look up the
    % documentation for 'cmap' for other color map options
    cmap = colormap(lines(Count));
    set(0,'DefaultAxesColorOrder',cmap)
    set(gca(),'ColorOrder',cmap);
    
    %2D Altitude vs Total Downrange Ballistic Trajectory
    fields = fieldnames(stateStruct);
    figure(900)
    for n = 1:Count
        %[~, iApogee] = max(stateStruct.(fields{n}).data(:, 6)); % find location of max z
        distBoost = vecnorm([stateStruct.(fields{n}).data(:, 4), stateStruct.(fields{n}).data(:, 5)], 2, 2);
        plot(distBoost,...
                -stateStruct.(fields{n}).data(:, 6), ...
                DisplayName=Design_Input.Properties.RowNames{n}, Color=cmap(n, :))
        if n == 1
            hold on
        end
    end
    xlabel('Total Distance Traveled [m]');
    ylabel('Altitude Achieved [m]');
    title('Boost 2D Altitude vs Total Range Trajectory');
    legend();
    grid on
    hold off
    
    %Reset default color order
    set(0,'DefaultAxesColorOrder','default')

end


   