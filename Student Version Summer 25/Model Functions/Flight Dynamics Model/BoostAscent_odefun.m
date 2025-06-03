function [dSdt] = BoostAscent_odefun(t,S,consts,thrustVec,Time)
% BOOSTASCENT_odefun
% The goal of this function is to output the derivatives of the current
% state (held in the vector 'S'). These derivatives are calculated using
% fundamental physics, and then are packed together in "dSdt'. 
% 
% ODE45 will pass in time information to this function (since the physics
% are time dependentant), and then will numerically integrate the 
% derivative in the output by something like:
%       S(i+1) = S(i) + dSdt*dt 
% where ODE45 is picking dt.
%
% We need information about the position, velocity, and mass, thus these
% are all of the quantities held in the state vector. Also passed in is a
% constants vector that holds information that the will be needed to
% calculate quanties of interest in this funciton. Finally, there is also
% the trust vector and time vector associated with the thrust vector so
% that we can calculate thrust at any given time using interpolation.

%% Unpack the state vector
Vx = S(1); % inertial velocity in x-direction [m/s]
Vy = S(2); % inertial velocity in y-direction [m/s]
Vz = S(3); % inertial velocity in z-direction [m/s]
x  = S(4); % position in x (inertial) [m]
y  = S(5); % position in y (inertial) [m]
z  = S(6); % position in z (inertial) [m]
m  = S(7); % current total mass [kg]

%% Unpack Constants Vector
% Basic Properties
g       = consts(1); % Accelatation due to gravity [m/s^2]
rho_w   = consts(2); % Density of water [kg/m^3]
rho_a   = consts(3); % Density of air [kg/m^3]
mu_k    = consts(4); % Launch rail coefficient of dynamic friction []
% Vehicle info
A_exit  = consts(5); % Area of bottle outlet [m^2]
C_D     = consts(6); % C_D of the vehicle (assume zero lift) []
S_ref   = consts(7); % Wing reference area [m^2]
m_empty = consts(8); % Weight of the rocket with no water [kg]
% Launch Direction
eliv    = consts(9); % Launch elevation [degrees]
azim    = consts(10); % Launch Azimuth [degrees], measured CW from north when looking down on the map; also known as compass heading

% /////////////////////////////////////////////////////////////////////////
% MODIFY THIS SECTION
% /////////////////////////////////////////////////////////////////////////
%% Magnitude of velocity vector (relative = inertial, assume no wind)
V_mag= sqrt(Vx^2+Vy^2+Vz^2); % Magnitude of velocity [m/s]

%% Set velocity direction & launch rail friction
% Constant if still on the rails
f_rails = 0; % friction due to the rails, preallocate to zero
if sqrt(x^2+y^2+z^2) <= 0.5 % still on rails
    % We need to add 90 degrees to elevation as it is normally measured
    % from the horizon but the spherical to cartiesian coordinate
    % conversion measures from the positive-z (strait down)
    hBod_x = sind(eliv+90)*cosd(azim); % heading of the body (inertial pointing direction of the body)
    hBod_y = sind(eliv+90)*sind(azim);
    hBod_z = cosd(eliv+90);
    f_rails = mu_k*m*g; % making friction non-zero only if on the rails, assuming the force on the rails is the full weight (this is a strange assumption but attempts to account for extra forces during launch)
else %free flight, velocity is into the relative wind
    hVrel = 1/V_mag*[Vx; Vy; Vz]; % relative velocity unit vector (direction of vector)
    hBod_x = hVrel(1); % x component of relative velocity unit vector 
    hBod_y = hVrel(2); % y component of relative velocity unit vector 
    hBod_z = hVrel(3); % z component of relative velocity unit vector 
end

%% Interpolate from the thrust curve
% We want to do this so that we allow ode 45 to choose its own time step size
if t < 0.5
    T = interp1(Time, thrustVec, t);
else
    T = 0;
end

%% Calculate total drag
D = C_D*1/2*rho_a*(V_mag^2)*S_ref;

%% Sum the forces
% Assume that all forces exept gravity act in (or against) the direction
% that the body is pointing. This means that the force in any component
% direction is the sum of the forces multiplied by the component of the
% unit vector in line with the body pointing (hBod). We then assume that
% gravity acts only in -z.

% Sum of the forces
%Since forces T, D, f_rails are always aligned w/ relative velocity vector (body frame),
%components of force in inertial frame can be found my multiplying sum of magnitude of forces by relative velocity unit vector.
Fx = (T-D-f_rails)*hBod_x; %Sum of forces in x inertial frame
Fy = (T-D-f_rails)*hBod_y; %Sum of forces in y inertial frame
Fz = (T-D-f_rails)*hBod_z + m*g; %Sum of forces in z inertial frame

%% Calculate the derivatives we need to pass out
% Change in postition = Acceleration = F/m (Newtons 2nd Law) - This value
% comes from our physics since the sum of the forces is time dependent
dVdt_x = Fx/m;
dVdt_y = Fy/m;
dVdt_z = Fz/m;

% Now we do some error checking:
% If we are not yet producing enough thrust to get a positive component of
% acceleration while very near the base of the the launcher, this means
% that the rocket is trying to slide backwards down the rails. There is
% backing plate that will stop this from happening, so if this backwards
% motion is detected by the following test, set the acceleration to zero
if sqrt(x^2+y^2+z^2) <= 0.05 && dVdt_z > 0% still on rails
    dVdt_x = 0;
    dVdt_y = 0;
    dVdt_z = 0;
end

% Change in postition = Velocity - we already have this! So we can just ask
% MATLAB to multiply our veocity by the timestep
dxdt = Vx;
dydt = Vy;
dzdt = Vz;

% Mass flow rate
if m > m_empty % if water is not  yet exausted as measured by weight
    mDot = -sqrt(rho_w*A_exit*T);
else % ignore any mass change from expulsed air
    mDot = 0;
end
% /////////////////////////////////////////////////////////////////////////
% END OF SECTION TO MODIFY
% /////////////////////////////////////////////////////////////////////////
%% Pass out the derivatives in time for ODE45 to intake
dSdt = [dVdt_x; dVdt_y; dVdt_z; dxdt; dydt; dzdt; mDot];
end