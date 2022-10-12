%% CL, CD vs alpha
% Trevor Burgoyne 4 Oct 2022

% Paths for data loading
ROOT_DIR = "C:/Users/Trevor/Desktop/AEM 4602W/Fluids Lab/Fluids Lab Data/";
FORCE_DIR = ROOT_DIR + "Force Measurements/";
ANGLES = ["-5", "00", "05", "10", "12", "14", "15", "18"];

% Arrays to store CD, CL, and a
CL_arr = zeros(1, length(ANGLES));
CD_arr = zeros(1, length(ANGLES));
a_arr  = zeros(1, length(ANGLES));

% Airfoil properties
c = .254; % m +/- .005m, chord length
b = .670; % m +/- .005m, wing span
S =  b*c; % m^2, approx. wing area

% lb -> N = (lb) * 4.448 N/lb
LB_TO_N = 4.448;

for i = 1:length(ANGLES)
    path = FORCE_DIR + "force_mes_a_" + ANGLES(i) + ".mat";
    data = load(path); % lab data, with P, rho, v, Fn, Fa, a
    
    % q  = .5*rho*v^2, dynamic pressure
    q = .5 * data.rho * data.v^2;
    
    % L  = -Fa*sin(a) + Fn*cos(a), Lift Force
    L = -data.Fa*sind(data.a) + data.Fn*cosd(data.a);
    L = L * LB_TO_N;
    
    % D = Fa*cos(a) + Fn*sin(a), Drag Force
    D = -data.Fa*cosd(data.a) + data.Fn*sind(data.a);
    D = D * LB_TO_N;
    
    % CL = L / q*S, coefficient of lift
    CL = L / (q*S);
    
    % CD = D / (q*S), coefficient of drag
    CD = D / (q*S);
    
    % Store in arrays for graphing
    CL_arr(i) = CL;
    CD_arr(i) = CD;
    a_arr(i)  = data.a;
end

%% CL vs Angle of Attack

plot(a_arr, CL_arr, "-o")
xlabel("Angle of Attack (deg)")
ylabel("CL")

%% CD vs Angle of Attack
plot(a_arr, CD_arr, "-o")
xlabel("Angle of Attack (deg)")
ylabel("CD")
