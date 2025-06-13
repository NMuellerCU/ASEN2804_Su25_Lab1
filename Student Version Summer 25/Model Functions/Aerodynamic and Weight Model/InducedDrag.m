function [InducedDrag_Data, InducedDrag_Model_Names] =...
    InducedDrag(Design_Input,WingLiftModel,WingLiftCurve,WingDragCurve,WingGeo_Data,Parasite_Drag_Data,Count,Benchmark,Plot_Induced_Data)
%% Induced Drag Model Function Summary
% This function evaluates different Oswalds Efficiency Factor models for 
% use in your drag polar model.  It compiles and outputs a variables table
% for three different Oswalds models (mod1, mod2, mod3).  For each Oswalds
% model, the Oswalds Efficiencty Factor (eo_mod1,2,3), and the k1 values
% (k1_mod1,2,3) are outputted in the InducedDrag_Data table.  Note that
% depending on the Oswalds models chosen to evaluate, you may or may not
% need information from teh WingGeo_Data table from the WingGeo function.
% Additionally, this code supports the calculation of the k2 values for
% evaluating non-symmetric airfoil design, but is not required.

%% Outputs:
%
% InducedDrag_Data:
%   Table containing oswalds info and calculated k1 and k2 values for three
%   different models for oswalds (denoted by suffixes _mod1, _mod2, and
%   _mod3)(columns) for each input from the design input spreadsheet (rows)


%% Preallocate variables of interest
eo_mod1 = zeros(Count, 1); % Oswalds for Model #1
eo_mod2 = zeros(Count, 1); % Oswalds for Model #2
eo_mod3 = zeros(Count, 1); % Oswalds for Model #3
k1_mod1 = zeros(Count, 1); % k1 for Model #1
k1_mod2 = zeros(Count, 1); % k1 for Model #2
k1_mod3 = zeros(Count, 1); % k1 for Model #3
k2_mod1 = zeros(Count, 1); % k2 for Model #1
k2_mod2 = zeros(Count, 1); % k2 for Model #2
k2_mod3 = zeros(Count, 1); % k2 for Model #3


% NOTE: k2 values not required if only symmetric airfoils used; however,
% this version of the code includes it as an option

%% Loop through different configurations
calculatek1 = @(eo, AR) 1 / (pi * eo * AR);
calculatek2 = @(k1, min_CL) -2 * k1 * min_CL;
for n = 1:Count 
% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
    %Find CL min Drag value of wing drag polar to estimate k2
    [CD_min,CD_min_index] = min(WingDragCurve{n,:});
    CL_minD = WingLiftCurve{n, CD_min_index};
    AR = Design_Input.AR_w(n);

    %Cavallo Oswalds Model 1 (Baseline required)
    Model1_Name = 'Cavallo'; %Name of first Oswald's Model
    eo_mod1(n) = 1.78 * (1 - 0.045 * AR^0.68) - 0.64;
    k1_mod1(n) = calculatek1(eo_mod1(n), AR);
    k2_mod1(n) = calculatek2(k1_mod1(n), CL_minD);

    %Student Option Oswalds Model 2 
    Model2_Name = 'Kroo'; %Name of second Oswald's Model

    % Find u which is equal to span efficiency factor
    f_taper = 0.0524*Design_Input.Taper_w(n)^4-0.15*Design_Input.Taper_w(n)^3+0.1659*Design_Input.Taper_w(n)^2-0.0706*Design_Input.Taper_w(n)+0.0119;
    u_kroo = 1 / (1 + f_taper * Design_Input.AR_w(n, 1));
    if Design_Input.QuarterSweep_w(n) ~= 0
        u_kroo = u_kroo * cosd(Design_Input.QuarterSweep_w(n));
    end

    % Back calculate wingspan since it is not given
    b_squared = AR * Design_Input.Sref_w(n);
    s_kroo = 1 - 2 * (Design_Input.Dia_f(n)^2 / b_squared);
    Q_kroo = 1 / (u_kroo * s_kroo);
    K = 0.38;
    P = K * Parasite_Drag_Data.CDo(n);

    eo_mod2(n) = 1 / (Q_kroo + P * pi * AR);
    k1_mod2(n) = calculatek1(eo_mod2(n), AR);
    k2_mod2(n) = calculatek2(k1_mod2(n), CL_minD);
   
    %Student Option Oswalds Model 3
    Model3_Name = 'Schaufele'; %Name of third Oswald's Model

    Q_schaufele = 1.03;
    
    eo_mod3(n) = 1 / (Q_schaufele + P * pi * AR);
    k1_mod3(n) = calculatek1(eo_mod3(n), AR);
    k2_mod3(n) = calculatek2(k1_mod3(n), CL_minD);
    
% /////////////////////////////////////////////////////////////////////////

% END OF SECTION TO MODIFY

% ///////////////////////////////////////////////////////////////////////// 





%% Oraganize into table for output

InducedDrag_Data = table(eo_mod1, eo_mod2, eo_mod3, k1_mod1, k1_mod2, k1_mod3, k2_mod1, k2_mod2, k2_mod3);
InducedDrag_Model_Names = {Model1_Name, Model2_Name, Model3_Name};


 %Isolated Induced Drag Coefficients

 CDi_w{n} = (WingLiftCurve{n,:}').^2/ ...
    (pi*WingLiftModel.e(n)*Design_Input.AR_w(1));

 CDi_mod1{n} = (WingLiftCurve{n,:}').^2.*InducedDrag_Data.k1_mod1(n)...

 +InducedDrag_Data.k2_mod1(n).*((WingLiftCurve{n,:}'));

 CDi_mod2{n} = (WingLiftCurve{n,:}').^2.*InducedDrag_Data.k1_mod2(n)...

 +InducedDrag_Data.k2_mod2(n).*((WingLiftCurve{n,:}'));

 CDi_mod3{n} = (WingLiftCurve{n,:}').^2.*InducedDrag_Data.k1_mod3(n)...

 +InducedDrag_Data.k2_mod3(n).*((WingLiftCurve{n,:}'));

 CDi_benchmark{n} = Benchmark.CD-min(Benchmark.CD); % subtract off minCD to get CDi



end

%% Plots for this function (Figure 400 - 499)

if Plot_Induced_Data == 1

 

 % CDi comparison for different Oswalds Models

 for n=1:Count

 figure(399+n)

 hold on

 plot(WingLiftCurve{n,:},CDi_w{n},'--');

 plot(WingLiftCurve{n,:},CDi_mod1{n});

 plot(WingLiftCurve{n,:},CDi_mod2{n});

 plot(WingLiftCurve{n,:},CDi_mod3{n});

 %plot(WingLiftCurve{1,:},CDi_benchmark,'--');


 xlabel('Coefficient of Lift (CL) [ ]');

 ylabel('Induced Drag (CDi) [ ]');

 title(sprintf('Induced Drag (CDi) Modeling Config: %d', n));

 %legend('3D Wing',Model1_Name,Model2_Name,Model3_Name,'Benchmark','Location','southeast');

 legend('3D Wing',Model1_Name,Model2_Name,Model3_Name, 'Location','southeast');

 grid on

 

 hold off

 end

 

 % Reset default color order

 set(0,'DefaultAxesColorOrder','default')

end



end