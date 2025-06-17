function [SensitivityData] =...
    SensitivityModeling(Design_Input,WingGeo_Data,Airfoil,ATMOS,SensVar, ModelRow, Material_Data,Plot_Sensitivity_Data, WingLiftModel, WingLiftCurve, WingDragCurve, Benchmark, AoA_Count,AirfoilLiftCurve, Component_Data)

%% Sensitivity Data Summary:
% takes in all data and a variable that we will be performing sensitivity
% analysis on, and then shows plots of different variables vs that
% sensitivity variable, where the sensitivity variable varies by +/- 20,
% 40, 60, 80 100% of its original value.
%% outputs

Design_Input = Design_Input(ModelRow,:);
WingGeo_Data = WingGeo_Data(ModelRow,:);
Airfoil = Airfoil(ModelRow,:);

baseVal = double(Design_Input.(SensVar));

%THIS IS THE RANGE YOU CHOSE, please modify this if you want a different
%range (Note: its in percent, ex: -0.8 is -80% of the original value)
pct_sens = (-0.8:0.2:1.0)';
SensVarRange = baseVal*(1+pct_sens);

n = numel(SensVarRange);

Parasite_Drag_Data_Mod = cell(n,1);
FF_Table = cell(n,1);
Q_Table = cell(n,1);
Re_Table = cell(n,1);

  secWing = repmat(struct( ...
    'b',[], 'cr',[], ...
    'ct',[], 'V_box',[], ...
    'W_box',[], 'Load_w',[]), n,1);

  secFuse = repmat(struct( ...
    'V_max',[], 'W_max',[]), n,1);



switch SensVar
    case {'AR_w','Taper_w','Sref_w',}
        Design_Mod = repmat(Design_Input,n,1); %prealocates all of the modified Design_Inputs
        apogee = ones(n,1)*17.5;  %for glideDescent.m
        for i = 1:n
            Design_Mod(i,:).(SensVar)(1) = SensVarRange(i);

            [Parasite_Drag_Data_Mod{i}, FF_Table{i}, Q_Table{i}, Re_Table{i}] = ...
              ParasiteDrag(Design_Mod(i,:), Airfoil, WingGeo_Data, ATMOS, 1, 0);

            [InducedDrag_Data,InducedDrag_Model_Names] = ...
                InducedDrag(Design_Mod(i,:),WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Parasite_Drag_Data_Mod{i},1,Benchmark,0);

            [DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] = ...
                DragPolarV(Parasite_Drag_Data_Mod{i},InducedDrag_Data, InducedDrag_Model_Names,Design_Mod(i,:),AoA_Count,WingLiftCurve,WingDragCurve,AirfoilLiftCurve,Airfoil,Benchmark,1,0);
            
            [LD_mod1,LD_mod2,LD_mod3] = ...
                LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,WingDragCurve,AoA_Count,1,0);

            [Weight_Data,~] = ...
                Weight(Design_Mod(i,:),1,WingGeo_Data,Airfoil,Material_Data,Component_Data,9.81,0);
            
 
            [GlideData_Mod(i,:)] = GlideDescent(LD_mod3, apogee, Design_Mod(i,:), ATMOS, Weight_Data, WingLiftModel, WingLiftCurve,WingDragCurve,1,0); %Must select LD of your best model          

            secWing(i) = computeSecWing(Design_Mod(i,:), Airfoil, Material_Data);
            

        end

    case {'Fuse_Mat', 'Length_f','Dia_f'}
        Design_Mod = repmat(Design_Input,n,1); %prealocates all of the modified Design_Inputs
        apogee = ones(n,1)*17.5; 

        for i = 1:n
            Design_Mod(i,:).(SensVar)(1) = SensVarRange(i);

            [Parasite_Drag_Data_Mod{i}, FF_Table{i}, Q_Table{i}, Re_Table{i}] = ...
              ParasiteDrag(Design_Mod(i,:), Airfoil, WingGeo_Data, ATMOS, 1, 0);

            [InducedDrag_Data,InducedDrag_Model_Names] = ...
                InducedDrag(Design_Mod(i,:),WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Parasite_Drag_Data_Mod{i},1,Benchmark,0);

            [DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] = ...
                DragPolarV(Parasite_Drag_Data_Mod{i},InducedDrag_Data, InducedDrag_Model_Names,Design_Mod(i,:),AoA_Count,WingLiftCurve,WingDragCurve,AirfoilLiftCurve,Airfoil,Benchmark,1,0);
            
            [LD_mod1,LD_mod2,LD_mod3] = ...
                LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,WingDragCurve,AoA_Count,1,0);

            [Weight_Data,~] = ...
                Weight(Design_Mod(i,:),1,WingGeo_Data,Airfoil,Material_Data,Component_Data,9.81,0);
            
 
            [GlideData_Mod(i,:)] = GlideDescent(LD_mod3, apogee, Design_Mod(i,:), ATMOS, Weight_Data, WingLiftModel, WingLiftCurve,WingDragCurve,1,0); %Must select LD of your best model

            secFuse(i) = computeSecFuse(Design_Mod(i,:), Material_Data);
        end
        
    case {'Thick_w'}
        Airfoil_Mod = repmat(Airfoil,n,1); %prealocates all of the modified Design_Inputs
        apogee = ones(n,1)*17.5; 

        for i = 1:n
            Airfoil_Mod(i,:).(SensVar)(1) = SensVarRange(i);

            
            [WingLiftModel,AoA,AoA_Count,AirfoilLiftCurve,WingLiftCurve,WingDragCurve] =...
                                WingLiftDrag(Design_Input,Airfoil_Mod(i,:),1,Benchmark,Plot_Wing_Data); 

            [Parasite_Drag_Data_Mod{i}, FF_Table{i}, Q_Table{i}, Re_Table{i}] = ...
              ParasiteDrag(Design_Input, Airfoil_Mod(i,:), WingGeo_Data, ATMOS, 1, 0);

            [InducedDrag_Data,InducedDrag_Model_Names] = ...
                InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Parasite_Drag_Data_Mod{i},1,Benchmark,0);

            [DragPolar_mod1,DragPolar_mod2,DragPolar_mod3] = ...
                DragPolarV(Parasite_Drag_Data_Mod{i},InducedDrag_Data, InducedDrag_Model_Names,Design_Input,AoA_Count,WingLiftCurve,WingDragCurve,AirfoilLiftCurve,Airfoil_Mod(i,:),Benchmark,1,0);
            
            [LD_mod1,~,~,~] = ...
                LD(Benchmark,DragPolar_mod1,DragPolar_mod2,DragPolar_mod3,WingLiftCurve,WingDragCurve,AoA_Count,1,0);

            [Weight_Data,~] = ...
                Weight(Design_Input,1,WingGeo_Data,Airfoil_Mod(i,:),Material_Data,Component_Data,9.81,0);
            

            [GlideData_Mod(i,:)] = GlideDescent(LD_mod1, apogee,Design_Input, ATMOS, Weight_Data, WingLiftModel, WingLiftCurve,WingDragCurve,1,0); %Must select LD of your best model

            secWing(i) = computeSecWing(Design_Input, Airfoil_Mod(i,:), Material_Data);
        end
      
    otherwise
      error('SensitivityModeling:UnknownVar',...
            'Unknown sensitivity variable "%s"', SensVar);

end
SensitivityData.SensVarRange = SensVarRange;
SensitivityData.ParasiteDragData = Parasite_Drag_Data_Mod;
SensitivityData.GlideDescent = GlideData_Mod;
% SensitivityData.FFTable = FF_Table; %uncomment these if you need their data
% SensitivityData.QTable = Q_Table;
% SensitivityData.ReTable = Re_Table;
SensitivityData.WingSecondaries = secWing;
SensitivityData.FuseSecondaries = secFuse;


%% Plotting
if Plot_Sensitivity_Data == 1
      wingVars = {'AR_w','Taper_w','Sref_w','Thick_w'};   % whatever you use
      fuseVars = {'Fuse_Mat','Length_f','Dia_f','Amax_f','Abase_f'};
      
    if ismember(SensVar, wingVars)
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
        figure(1); clf;    
        n = numel(SensVarRange);
        cmap = colormap([linspace(0,1,n)',linspace(1,0,n)', linspace(1,0,n)']);
        ax = axes;          
        hold(ax, 'on');
        for i = 1:numel(GlideData_Mod.bestGlide) %coppied and pasted code from GlideDescent.m
            plot(ax,...
                [0, GlideData_Mod.bestGlide(i)], [apogee(i), 0],...
                'Color',cmap(i, :), 'DisplayName',sprintf('AR = %.2f best glide', SensVarRange(i)))
        end
        xlabel('Glide distance [m]');
        ylabel('Height [m]');
        title('Best Glide Plots');
        legend();
        grid on
        hold off
    end
    
    elseif ismember(SensVar, fuseVars)
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
        
        figure(1); clf;    
        n = numel(SensVarRange);
        cmap = colormap([linspace(0,1,n)',linspace(1,0,n)', linspace(1,0,n)']);
        ax = axes;          
        hold(ax, 'on');
        for i = 1:numel(GlideData_Mod.bestGlide) %coppied and pasted code from GlideDescent.m
            plot(ax,...
                [0, GlideData_Mod.bestGlide(i)], [apogee(i), 0],...
                'Color',cmap(i, :), 'DisplayName',sprintf('AR = %.2f best glide', SensVarRange(i)))
        end
        xlabel('Glide distance [m]');
        ylabel('Height [m]');
        title('Best Glide Plots');
        legend();
        grid on
        hold off
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
