function [Boost_Initial_Stab,Boost_75_Stab,Boost_50_Stab,Boost_25_Stab,Boost_Empty_Stab,Glide_Stab,STAB_SM_SUMMARY,STAB_Xcg_SUMMARY,STAB_Xnp_SUMMARY,STAB_Vh_SUMMARY,STAB_Vv_SUMMARY,STAB_GLIDE_h1_SUMMARY]...
    = Stability(Design_Input, Count, CG_Data, WingGeo_Data, GlideData, WingLiftModel, Component_Data,Plot_Stability_Data)
%% Weight Model Summary 
% This is a Instructor team- only function for use in calculating static
% margin and other stability terms for designs

% Note that in its current form, this function assumes only h1 and v1

%% Outputs:
%
% <flight phase>Stability:
%   Table containing Static Margin tail coefficients and more for boost or
%   glide

%% Loop over boost and glide
for i = 1:6 %Discrete phases of flight (1 = Boost Full Water, 2 = Boost 75% water, 3 = Boost 50% water, 4 = Boost 25% water, 5 = Boost Empty Water, 6 = Glide Empty & Trim
    % Preallocate for each phase of flight
    Xnp = zeros(Count,1); %Longitudinal distance to neutral point from nose of aircraft [m]
    h_np = zeros(Count,1); %Longitudinal distance to neutral point from nose of aircraft [% MAC]
    Xcg =  zeros(Count,1); %Longitudinal distance to aircraft cg from nose of aircraft [m]
    h_cg = zeros(Count,1); %Longitudinal distance to aircraft cg from nose of aircraft [% MAC]
    Xac_wb = zeros(Count,1); %Longitudinal distance to wing aerodynamic center from nose of aircraft [m]
    h_ac_wb = zeros(Count,1); %Longitudinal distance to wing aerodynamic center from nose of aircraft [% MAC]
    X_ac_h1 = zeros(Count,1); %Longitudinal distance to horz tail#1 aerodynamic center from nose of aircraft [m]
    X_ac_v1 = zeros(Count,1); %Longitudinal distance to vert tail#1 aerodynamic center from nose of aircraft [m]
    X_ac_v2 = zeros(Count,1); %Longitudinal distance to vert tail#2 aerodynamic center from nose of aircraft [m]
    SM = zeros(Count, 1); %Static Margin
    lt_h1 = zeros(Count, 1); %moment arm for horizontal tail 1(aircraft cg to horz tail ac) [m]
    lv_v1 = zeros(Count, 1); %moment arm for vertical tail 1(aircraft cg to ert tail ac) [m]
    lv_v2 = zeros(Count, 1); %moment arm for vertical tail 2 (aircraft cg to ert tail ac) [m]
    Vh = zeros(Count, 1); %Horizontal Tail Volume Coefficient h1
    Vv = zeros(Count, 1); %Vertical Tail Volume Coefficient v1
    Cm_alpha = zeros(Count, 1); %Longitudinal Static Stability Derivative [1/deg]
        i_t = zeros(Count, 1); %tail incidence angle [deg]
    Cm0 = zeros(Count, 1); %Zero lift pitching moment coefficient
    depsilon_dalpha_h1 = zeros(Count, 1); %change in downwash angle w/ change in angle of attack
    eps = zeros(Count, 1); %Downwash Angle [deg]
    at = zeros(Count, 1); %Horizontal stabilizer 3D lift curve slope
    AoA_h1 = zeros(Count, 1); %Angle of attack of horizontal tail 1 [deg]
    AoA_e = zeros(Count,1); %Wing-Aircaft Trim / Equilibrium AoA [deg]
   
    %% Loop through different configurations
    for n = 1:Count
        % Unpack some things
        L = Design_Input.Length_f(n);
        Sref = Design_Input.Sref_w(n);
        b = WingGeo_Data.b_w(n);
        c = WingGeo_Data.MAC_w(n);
        Sref_h1 = Design_Input.Sref_h1(n); %horizontal stab 1 planform
        Sref_v1 = Design_Input.Sref_v1(n); %vertical stab 1 planform
        Sref_v2 = Design_Input.Sref_v2(n); %vertical stab 2 planform
        a = WingLiftModel.a(n);
        X_LE_wing = Component_Data.X_LE_wing(n); %Distance from nose to leading edge of wing root [m]
        X_LE_h1 = Component_Data.X_LE_h1(n); %Distance from nose to leading edge of horz tail #1 root [m]
        X_LE_v1 = Component_Data.X_LE_v1(n); %Distance from nose to leading edge of vert tail #1 root [m]
        X_LE_v2 = Component_Data.X_LE_v2(n); %Distance from nose to leading edge of vert tail #2 root [m]
        
        % Coose options, need to hard code locations
        if i == 1 % Boost intial (full water)
            Xcg(n) = CG_Data.CG_tot(n);
            CL = 0; %Zero lift ballistic ascent assumption
            AoA_e(n)= 0; %Trim AoA set to zero for boost
        elseif i == 2 % Boost 75% water
            Xcg(n)  = CG_Data.CG_tot_75(n);
            CL = 0; %Zero lift ballistic ascent assumption
            AoA_e(n)= 0; %Trim AoA set to zero for boost
        elseif i == 3 % Boost 50% water
            Xcg(n) = CG_Data.CG_tot_50(n);
            CL = 0; %Zero lift ballistic ascent assumption
            AoA_e(n)= 0; %Trim AoA set to zero for boost
        elseif i == 4 % Boost 25% water
            Xcg(n) = CG_Data.CG_tot_25(n);
            CL = 0; %Zero lift ballistic ascent assumption
            AoA_e(n)= 0; %Trim AoA set to zero for boost
        elseif i == 5 % Boost Empty water
            Xcg(n) = CG_Data.CG_empty(n);
            CL = 0; %Zero lift ballistic ascent assumption
            AoA_e(n)= 0; %Trim AoA set to zero for boost
        else % Glide
            Xcg(n) = CG_Data.CG_empty(n);         
            CL = GlideData.CL_LDmax(n); %Coefficient of lift at L/D max
            AoA_e(n) = GlideData.AoA_LDmax(n)-WingLiftModel.AoA_0(n); %Trim Absolute Angle of attack at L/D max (deg)

        end
        
        %% Stabilizer Coefficients
        X_ac_h1(n) = X_LE_h1+(0.25*Design_Input.MAC_h1(n)); % horizontal tail #1 aerodynamic center
        X_ac_v1(n) = X_LE_v1+(0.25*Design_Input.MAC_v1(n)); %  vertical tail #1 aerodynamic center
        X_ac_v2(n) = X_LE_v2+(0.25*Design_Input.MAC_v2(n)); %  vertical tail #2 aerodynamic center
        % Horizontal and Vertical Tail Volume Coefficients for h1 and v1
        lt_h1(n) = abs(X_ac_h1(n) - Xcg(n)); %Moment arm to horz tail 1
        lv_v1(n) = abs(X_ac_v1(n) - Xcg(n)); %Moment arm to vert tail 1
        lv_v2(n) = abs(X_ac_v2(n) - Xcg(n)); %Moment arm to vert tail 2
        Vh(n) = (Sref_h1*lt_h1(n))/(Sref*c); %Horizontal Tail Volume Coefficient
        Vv(n) = (Sref_v1*lv_v1(n)+Sref_v2*lv_v2(n))/(Sref*b); % Combined Vertical tail volume coefficient (all vertical stabilizers)
        
        %%  Spiral Parameter
        % Omitted for now
        % Set desired spiral parameter
        % B = 6;
        % rearrangement of B=((xacv-xcg)*Gamma)/(b*CL);
        % Gamma = B*b*CL/(x_ac_v-xCG); % Dihedral

        %% Neutral point
        Xac_wb(n) = X_LE_wing+WingGeo_Data.x_MAC_w(n);
        h_ac_wb(n) = (Xac_wb(n))/c; % Wing body aero center
                
        %Span Efficiency Model for H1 stabilizer (Foundational source:
        %Horner, S; Fluid Dynamic Drag, Dayton, Ohio, 1951; f_taper
        %polynomical derived from Horner data in Nita, M; Scholz, D;
        %Estimating the Oswals FActor from Basic Aircraft Geometrical
        %Parameters; Hamburg University of Applied Sciences; 2012.
        if Design_Input.Sweep_h1(n) == 0 % If no sweep:
            f_taper = 0.0524*Design_Input.Taper_h1(n)^4-0.15*Design_Input.Taper_h1(n)^3+0.1659*Design_Input.Taper_h1(n)^2-0.0706*Design_Input.Taper_h1(n)+0.0119;
            e_h1(n)=1/(1+f_taper*Design_Input.AR_h1(n));
        else % If there is quarter chord sweep:
            f_taper = 0.0524*Design_Input.Taper_h1(n)^4-0.15*Design_Input.Taper_h1(n)^3+0.1659*Design_Input.Taper_h1(n)^2-0.0706*Design_Input.Taper_h1(n)+0.0119;
            e_h1(n)=1/(1+f_taper*Design_Input.AR_h1(n)).*cosd(Design_Input.Sweep_h1(n));
        end
        at_o = .092; % Approximate 2D flat plate lift curve slope value for horz tail [1/deg].  Assumes all horizontal and vertical stabilizers are flat plates (no airfoil).  Source:Mueller, T.J.; Aerodynamic Measurements at Low Reynolds Numbers for Fixed Wing Micro-Air Vehicles; Hesert Center for Aerospace Research, University of Notre Dame; 1999.
        at(n) = (at_o)/(1+(57.3*at_o)/(pi*e_h1(n)*Design_Input.AR_h1(n))); %3D lift curve slope for h1 stabilizer

        %%Change in Downwash vs AoA Table Lookup based on W.F. Philips,et al, "Estimating the Low-Speed Downwash Angle on an Aft Tail", Journal of Aircraft, Vol. 39, No. 4 (2002), pp. 600-608.
        XT = [0.5 1 2];
        VT = [0.5 0.4 0.35];
        depsilon_dalpha_h1(n) = interp1(XT,VT,lt_h1(n)/(WingGeo_Data.b_w(n)/2),'linear','extrap');
        if Xac_wb(n) < X_ac_h1(n) %Checks if configuration is canard or conventional
            h_np(n) = h_ac_wb(n) + Vh(n)*(at(n)/a)*(1-depsilon_dalpha_h1(n)); %NP in % MAC
            Xnp(n)= h_np(n)*c; %location of NP from nose [m]
        else
            h_np(n) = h_ac_wb(n) -Vh(n)*(at(n)/a); %NP calculation for canard configuration
            Xnp(n) = h_np(n)*c; %location of NP from nose [m]
        end

        %% Static Margin
        h_cg(n) = Xcg(n)/c; %location of CG in % MAC
        SM(n) = h_np(n)-h_cg(n); %Static Margin

        %% CMAlpha
        Cm_alpha(n) = -a*SM(n);

        if Xac_wb(n) < X_ac_h1(n)
            % Wing pitching moment coef and downwash angle
            Cmacwb = 0.00; %= cm at zero lift - assume that of airfoil = zero for flat plate
            %Cmacwb = -0.06; %= cm at zero lift - assume that of airfoil = -0.06 for Clark Y
            eps0 = 0; %Downwash angle on horz tail when wing at zero lift [deg]
            eps(n) = eps0+depsilon_dalpha_h1(n)*AoA_e(n); %Downwash angle [deg]
    
            %set Cmacwb = 0 and rearrange to find incidence angle
            i_t(n) = -(Cmacwb+a*AoA_e(n)*((h_cg(n)-h_ac_wb(n))-at(n)/a*Vh(n)*(1-depsilon_dalpha_h1(n))))/(Vh(n)*at(n))-eps0;
    
            %% Cm0>0
            Cm0(n) = Cmacwb+Vh(n)*at(n)*(i_t(n)+eps0);
    
            %% alpha tail
            AoA_h1(n) = AoA_e(n)-i_t(n)-eps(n);
        else
            % Wing pitching moment coef and downwash angle
            Cmacwb = 0; %= cm at zero lift - assume that of airfoil - zero for flat plate
            eps0 = 0; %Downwash angle on horz tail when wing at zero lift [deg]
            eps(n) = 0; %For canard, set downwash angle to zero [deg]
    
            %set Cmacwb = 0 and rearrange to find canardincidence angle
            i_t(n) = (Cmacwb+a*AoA_e(n)*((h_cg(n)-h_ac_wb(n))+at(n)/a*Vh(n)))/(Vh(n)*at(n)); %Canard incidence angle

            % Cm0>0
            Cm0(n) = Cmacwb-Vh(n)*at(n)*i_t(n); %Zero lift pitch moment for canard
    
            % alpha canard
            AoA_h1(n) = AoA_e(n)-i_t(n); %deg
        end
        
               
    end

    %Write Data to tables
    if i == 1 % Boost intial (full water)
        Boost_Initial_Stab = table(SM,Vh,Vv,AoA_e,AoA_h1,i_t,eps,depsilon_dalpha_h1,at,Cm_alpha,Cm0,Xcg,Xnp,Xac_wb,X_ac_h1,X_ac_v1,h_cg,h_np,h_ac_wb,lt_h1,lv_v1);
    elseif i == 2 % Boost 75% water
        Boost_75_Stab = table(SM,Vh,Vv,AoA_e,AoA_h1,i_t,eps,depsilon_dalpha_h1,at,Cm_alpha,Cm0,Xcg,Xnp,Xac_wb,X_ac_h1,X_ac_v1,h_cg,h_np,h_ac_wb,lt_h1,lv_v1);
    elseif i == 3 % Boost 50% water
        Boost_50_Stab = table(SM,Vh,Vv,AoA_e,AoA_h1,i_t,eps,depsilon_dalpha_h1,at,Cm_alpha,Cm0,Xcg,Xnp,Xac_wb,X_ac_h1,X_ac_v1,h_cg,h_np,h_ac_wb,lt_h1,lv_v1);
    elseif i == 4 % Boost 25% water
        Boost_25_Stab = table(SM,Vh,Vv,AoA_e,AoA_h1,i_t,eps,depsilon_dalpha_h1,at,Cm_alpha,Cm0,Xcg,Xnp,Xac_wb,X_ac_h1,X_ac_v1,h_cg,h_np,h_ac_wb,lt_h1,lv_v1);
    elseif i == 5 % Boost Empty water
        Boost_Empty_Stab = table(SM,Vh,Vv,AoA_e,AoA_h1,i_t,eps,depsilon_dalpha_h1,at,Cm_alpha,Cm0,Xcg,Xnp,Xac_wb,X_ac_h1,X_ac_v1,h_cg,h_np,h_ac_wb,lt_h1,lv_v1);
    else % Glide
        Glide_Stab = table(SM,Vh,Vv,AoA_e,AoA_h1,eps,i_t,depsilon_dalpha_h1,at,Cm_alpha,Cm0,Xcg,Xnp,Xac_wb,X_ac_h1,X_ac_v1,h_cg,h_np,h_ac_wb,lt_h1,lv_v1);
    end
end

STAB_SM_SUMMARY = table(Boost_Initial_Stab.SM,Boost_75_Stab.SM,Boost_50_Stab.SM,Boost_25_Stab.SM,Boost_Empty_Stab.SM,Glide_Stab.SM,'VariableNames',["SM_Full","SM_75","SM_50","SM_25","SM_Empty","SM_Glide"]);
STAB_Xcg_SUMMARY = table(Boost_Initial_Stab.Xcg,Boost_75_Stab.Xcg,Boost_50_Stab.Xcg,Boost_25_Stab.Xcg,Boost_Empty_Stab.Xcg,Glide_Stab.Xcg,'VariableNames',["Xcg_Full","Xcg_75","Xcg_50","Xcg_25","Xcg_Empty","Xcg_Glide"]);
STAB_Xnp_SUMMARY = table(Boost_Initial_Stab.Xnp,Boost_75_Stab.Xnp,Boost_50_Stab.Xnp,Boost_25_Stab.Xnp,Boost_Empty_Stab.Xnp,Glide_Stab.Xnp,'VariableNames',["Xnp_Full","Xnp_75","Xnp_50","Xnp_25","Xnp_Empty","Xnp_Glide"]);
STAB_Vh_SUMMARY = table(Boost_Initial_Stab.Vh,Boost_75_Stab.Vh,Boost_50_Stab.Vh,Boost_25_Stab.Vh,Boost_Empty_Stab.Vh,Glide_Stab.Vh,'VariableNames',["Vh_Full","Vh_75","Vh_50","Vh_25","Vh_Empty","Vh_Glide"]);
STAB_Vv_SUMMARY = table(Boost_Initial_Stab.Vv,Boost_75_Stab.Vv,Boost_50_Stab.Vv,Boost_25_Stab.Vv,Boost_Empty_Stab.Vv,Glide_Stab.Vv,'VariableNames',["Vv_Full","Vv_75","Vv_50","Vv_25","Vv_Empty","Vv_Glide"]);
STAB_GLIDE_h1_SUMMARY = table(Glide_Stab.AoA_e,Glide_Stab.AoA_h1,Glide_Stab.i_t,Glide_Stab.at,'VariableNames',["Trim AOA","h1 AoA","h1 Tail Incidence","h1 Tail Lift Curve Slope"]);

%% Plots for this function (Figure 1100 - 1199)
if Plot_Stability_Data == 1
    
    %Xcg,Xnp,SM shift over flight phases
    figure(1100)
    subplot(1,2,1);
    hold on
     %Standardize colormap for plots to align color to configurations evaluated
        cmap = colormap(hsv(Count)); % Sets color map for the specific number of variables (that last part is important)
        set(0,'DefaultAxesColorOrder',cmap) % overwrites default color order to what we just specified

    for n = 1:Count;
        plot(linspace(1,6,6),STAB_Xcg_SUMMARY{n,:},DisplayName = [Design_Input.Config{n}, ' Xcg'],Marker="o",LineStyle="-",MarkerSize=10);
    end
    for n = 1:Count;
        plot(linspace(1,6,6),STAB_Xnp_SUMMARY{n,:},DisplayName = [Design_Input.Config{n}, ' Xnp'],Marker="+",LineStyle="--",MarkerSize=10);
    end
    title('Longitudinal CG, NP Location Travel');
    ylabel('Distance from Aircraft Nose [m]');
    xticks([1 2 3 4 5 6])
    xticklabels({'Full','75% Water','50% Water','25% Water','Empty Boost','Empty Glide'});
    xlabel('Flight Phase');
    legend();
    grid on
    hold off
    
    subplot(1,2,2);
    hold on
    %Standardize colormap for plots to align color to configurations evaluated
    cmap = colormap(hsv(Count)); % Sets color map for the specific number of variables (that last part is important)
    set(0,'DefaultAxesColorOrder',cmap) % overwrites default color order to what we just specified

    for n = 1:Count;
        plot(linspace(1,6,6),STAB_SM_SUMMARY{n,:},DisplayName = [Design_Input.Config{n}, ' SM'],Marker="x",MarkerSize=10);
    end
    title('Longitudinal Static Margin Shift');
    xticks([1 2 3 4 5 6])
    xticklabels({'Full','75% Water','50% Water','25% Water','Empty Boost','Empty Glide'});
    xlabel('Flight Phase');
    ylabel('Static Margin');
    legend();
    grid on
    hold off
    
    %% Reset default plot color order
    set(0,'DefaultAxesColorOrder','default')

end



end

