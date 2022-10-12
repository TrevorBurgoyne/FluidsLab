%% Turbulence Profiles
% Trevor Burgoyne 4 Oct 2022

% Paths for data loading
ROOT_DIR = "C:/Users/Trevor/Desktop/AEM 4602W/Fluids Lab/Fluids Lab Data/";
HOTWIRE_DIR = ROOT_DIR + "Hotwire Measurements/";
ANGLES = ["-5", "00", "05", "20"];
N_DATAPOINTS = [11, 10, 15, 78];
V_AVG = 24.4247; % m/s, mean of all velocities, excluding nan

% Arrays to store data per angle
angle_data_arr = repmat(...
    struct(...
        "len_scale",[],...
        "v_normalized", [],...
        "v_rms", []...
    ), length(ANGLES), 1 ... 
);

% in -> m = (in) * .0254 m/in
IN_TO_M = 0.0254;

% Experiment properties
c =  .254; % m +/- .005m, airfoil chord length
x = .75*c; % m +/- .005m, distance from trailing edge to hotwire
L = 19.5*IN_TO_M; % m +/- 1/32in, distance from hotwire to top of tunnel

% Calibration constants: (E+offset)^2 = A + B*U^n
offset = -8.92; % V, voltage at zero flow
A      = -74.9; % constant from linear regression
B      =  53.7; % constant from linear regression
n      =  0.55; % exponent in King's Law that gave straightest fit


for i = 1:length(ANGLES)
   angle = ANGLES(i);
   angle_data_arr(i) = struct(...
      "len_scale", zeros(1, N_DATAPOINTS(i)),...
      "v_normalized", zeros(1, N_DATAPOINTS(i)),...
      "v_rms", zeros(1,N_DATAPOINTS(i))...
   );
   
   L0 = L - (x)*sind(str2double(angle)); % m, height of LE adjusted for angle of attack
   for j = 1:N_DATAPOINTS(i)
      path = HOTWIRE_DIR + "a_" + ANGLES(i) + "/data_" + j + ".mat";
      data = load(path); % lab data, with P, rho, v, a, y, and V_arr
      
      % len_scale = (L0 - y) / L0
      len_scale = (L0 - (data.y * IN_TO_M)) / L0;
      
      % hotwire velocity
      % from calibration: v_hotwire = (((E + offset)^2 - A)/B)^(1/n)
      v_hotwire_arr = (((data.V_arr + offset).^2 - A)./B).^(1/n);
      
      % v_rms = remove mean from velocities, take the average of their
      % squares, and then take the square root
      v_rms = sqrt( mean( ( v_hotwire_arr - mean(v_hotwire_arr) ).^2 ) );
      
      % v_hotwire / v_freestream
      % NOTE: for some reason, the pitot tube returned a speed of 'nan' for
      % all of our measurements at a = 0. This isn't a huge deal, since the
      % freestream was always set to the same speed, so a good
      % approximation for this case is to use the average of all other
      % velocities we measured
      if(isnan(data.v))
          data.v = V_AVG;
      end
      v_normalized = mean(v_hotwire_arr) / data.v;
      
       
       % Store in global arr
       angle_data_arr(i).len_scale(j) = len_scale;
       angle_data_arr(i).v_normalized(j) = v_normalized;
       angle_data_arr(i).v_rms(j) = v_rms;
   end
end


%% Graph

colors = ["green", "red", "blue", "black"];
shapes = ["-o", "-*", "-x", "-^"];

hold on
for i = 1:length(ANGLES)
    plot(angle_data_arr(i).v_normalized, angle_data_arr(i).len_scale, shapes(i), 'MarkerFaceColor', colors(i))
end
xlabel("v / vinf")
ylabel("y / Lo")
xlim([.3, 1])
legend('a=-5', 'a=0', 'a=5', "a=20", 'AutoUpdate', 'off', 'Location', 'northwest')










