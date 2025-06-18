function [SensitivityData] =...
    SensitivityModeling(Design_Input, Airfoil, SensVar, SensVarRange, n, Material_Data, Parasite_Drag_Data, GlideData, Plot_Sensitivity_Data)

%% Sensitivity Data Summary:
% takes in all data and a variable that we will be performing sensitivity
% analysis on, and then shows plots of different variables vs that
% sensitivity variable, where the sensitivity variable varies by +/- 20,
% 40, 60, 80 100% of its original value.
%% outputs

% Parasite_Drag_Data_Mod = cell(n,1);
% FF_Table = cell(n,1);
% Q_Table = cell(n,1);
% Re_Table = cell(n,1);

  secWing = repmat(struct( ...
    'b',[], 'cr',[], ...
    'ct',[], 'V_box',[], ...
    'W_box',[], 'Load_w',[]), n,1);

  secFuse = repmat(struct( ...
    'V_max',[], 'W_max',[]), n,1);



apogee = ones(n,1)*17.5;  %for glideDescent.m
wingVars = {'AR_w','Taper_w','Sref_w','Thick_w', "AR_h1"};   % whatever you use
fuseVars = {'Fuse_Mat','Length_f','Dia_f','Amax_f','Abase_f'};
      
switch SensVar
    case wingVars
        for i=1:n
            secWing(i) = computeSecWing(Design_Input(i,:), Airfoil(i, :), Material_Data);
        end

    case fuseVars
        for i = 1:n
            secFuse(i) = computeSecFuse(Design_Input(i,:), Airfoil(i, :), Material_Data);
        end

    otherwise
        error('SensitivityModeling:UnknownVar',...
            'Unknown sensitivity variable "%s"', SensVar);
end
SensitivityData.SensVarRange = SensVarRange;
SensitivityData.ParasiteDragData = Parasite_Drag_Data;
SensitivityData.GlideDescent = GlideData;
% SensitivityData.FFTable = FF_Table; %uncomment these if you need their data
% SensitivityData.QTable = Q_Table;
% SensitivityData.ReTable = Re_Table;
SensitivityData.WingSecondaries = secWing;
SensitivityData.FuseSecondaries = secFuse;


%% Plotting
if Plot_Sensitivity_Data == 1
    switch SensVar
        case wingVars
            % List of wing metrics to plot
            wingMetrics = {'b','cr','ct','V_box','W_box','Load_w'};
            for k = 1:numel(wingMetrics)
                m = wingMetrics{k};
                Y = [secWing.(m)];  % extract data vector
        
                figure('Name',m,'NumberTitle','off');
                plot(SensVarRange, Y, '-o');
                xlabel(SensVar, 'Interpreter','none');
                ylabel(m,      'Interpreter','none');
                title(['Sensitivity of ', m, ' to ', SensVar], 'Interpreter','none');
                grid on;
            end
        case fuseVars
            % List of fuselage metrics to plot
            fuseMetrics = {'V_max','W_max'};
            for k = 1:numel(fuseMetrics)
                m = fuseMetrics{k};
                Y = [secFuse.(m)];  % extract data vector
        
                figure('Name',m,'NumberTitle','off');
                plot(SensVarRange, Y, '-s');
                xlabel(SensVar, 'Interpreter','none');
                ylabel(m,      'Interpreter','none');
                title(['Sensitivity of ', m, ' to ', SensVar], 'Interpreter','none');
                grid on;
            end
    end
end

end

function sec = computeSecFuse(DI, Material_Data)
  %unpacking
  g = 9.81;
  rho = Material_Data.(DI.Fuse_Mat{1});
  L = DI.Length_f;
  d = DI.Dia_f;
  Amax = pi*0.25*(d^2);
  % V = DI.V_o;
  %maximum cylinder volume and weight
  V_max = Amax * L;
  W_max = rho * V_max * g;
  sec = struct( ...
    'V_max', V_max, ...
    'W_max', W_max ...
  );
end

function sec = computeSecWing(DI, AF, Material_Data)
  % constants
  g = 9.81;
  rho = Material_Data.(DI.Wing_Mat{1})(1);  

  % unpack
  S = DI.Sref_w;
  AR = DI.AR_w;
  lambda = DI.Taper_w;
  tc = AF.Thick_w;
  % Sweep = DI.Sweep_w;
  % V = DI.V_o;
  % X_thick = AF.X_thick_w;

  % geometry
  b = sqrt(AR * S);
  cr = 2*S / (b*(1+lambda));
  ct = lambda * cr;

  % volumes
  V_box = 0.5*(cr^2 + ct^2)*tc*b;    % uniform‐thickness box
  W_box = rho * V_box * g;

  % loading
  Load_w = W_box / S;

  sec = struct( ...
    'b', b, ...
    'cr', cr, ...
    'ct', ct, ...
    'V_box', V_box, ...
    'W_box', W_box, ...
    'Load_w',Load_w   ...
  );
end
