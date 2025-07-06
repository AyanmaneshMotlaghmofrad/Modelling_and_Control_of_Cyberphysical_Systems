clc
clear variables
close all
format compact

%% loading the data

load localization_data.mat

size_D = size(D);
q = size_D(1);
n = size_D(2);

%solution
target_position = 37;
x_tilde = zeros(n,1);
x_tilde(target_position) = 1;

a_support_idx = [1,10,14,16,17];
a_tilde_support = zeros(q,1);
a_tilde_support(a_support_idx) = 1;


% solving the problem
% solver parameters

delta = 1e-10;
delta_my_ISTA = 1e-6;
max_iteration = 1e4;


% to avoid numerical problems component of 
% D have large values compared to I.
%normalizing D

G = normalize([D, eye(q)]);

                       
%ISTA for localization full LASSO +++++++++++++++++++++++++++++++++++++++++
nu = 2.33/(norm(G,2)^2);
lambda_1 = 14;
lambda_2 = 11.4;

% lambda_1 = 10;
% lambda_2 = 10;
%initialization
sensor_position_hat_ISTA = zeros(n + q,max_iteration);


for k = 1:max_iteration
    
    sensor_position_hat_ISTA(:,k+1) = prox_l1_vec(sensor_position_hat_ISTA(:,k) - ... 
                        nu*G'*(G*sensor_position_hat_ISTA(:,k)-y),[nu*lambda_1;nu*lambda_2],n,q);
  
    if(norm(sensor_position_hat_ISTA(:,k+1)-sensor_position_hat_ISTA(:,k))^2<delta)  % condition for maximum iteration number
        break
    end
end
last_col_ISTA = find(any(sensor_position_hat_ISTA, 1), 1, 'last'); %finding the last non-zero column

% cleaning the data %the x_i  =  {0, 1}
senseor_position_hat_ISTA = sensor_position_hat_ISTA(1:n,last_col_ISTA);
attack_position_hat_ISTA = sensor_position_hat_ISTA(n+1:end,last_col_ISTA);

non_zero_x = find(senseor_position_hat_ISTA);
non_zero_a = find(attack_position_hat_ISTA);

cleaning_threshold = 0.4;

% the state of the cell can be {0,1}
% so we have to clean the data
for i = 1:length(non_zero_x)
    if(senseor_position_hat_ISTA(non_zero_x(i))<cleaning_threshold)
        senseor_position_hat_ISTA(non_zero_x(i)) = 0;
    
    else
        senseor_position_hat_ISTA(non_zero_x(i)) = 1;
    end

end
target_position = find(senseor_position_hat_ISTA)
% the check whether position of the attacks are detected correctly.
% [a_support_idx' non_zero_a]

% identifying the attack
a_support_matrix = zeros(q);

for i = 1:length(non_zero_a)
    a_support_matrix(non_zero_a(i),non_zero_a(i)) = 1;
end

a_hat = pinv(a_support_matrix)*(y-D*senseor_position_hat_ISTA);

% just to check the index of the attacks
[find(a_hat) find(attack_position_hat_ISTA)]



% now that the position of the target is clear, we can solve the 
% problem for identifying the attack

%% plotting

figure(1)
bar(1:length(find(a_hat)),a_hat(find(a_hat)),'g')

xticks([1, 2, 3, 4, 5]);  % Specify the x-tick positions
xticklabels({'1', '10', '14', '16','17'});  % Set the labels for the x-ticks
xlabel('attacked sensor')
ylabel('the value of the attack')
title('The value of the attack corrupting the measurement')                    

%% finding the position of the sensors

[~, sensor_position] = max(D, [], 2);

grid_matrix = nan(10,10);

% Sensor cell numbers (each sensor is assigned a cell number)
sensorCells = sensor_position';

% Initialize a 10x10 matrix to count sensors in each cell
gridCounts = zeros(10, 10);

% Count sensors in each cell using linear indexing
for i = 1:length(sensorCells)
    if sensorCells(i) >= 1 && sensorCells(i) <= 100
        gridCounts(sensorCells(i)) = gridCounts(sensorCells(i)) + 1;
    else
        error('Sensor cell number %d is out of range (should be 1-100)', sensorCells(i));
    end
end

% Create a new matrix that defines color mapping for cells:
% 0: empty, 1: one sensor, 2: two or more sensors.
displayGrid = zeros(10,10);
displayGrid(gridCounts == 1) = 1;
displayGrid(gridCounts >= 2) = 2;  % same color for two or more sensors

% Define a custom black and white colormap:
% Row 1: Black for empty cells
% Row 2: White for cells with exactly one sensor
% Row 3: Gray for cells with two or more sensors
cmap = [0 0 0;       % Black
        1 1 1;       % White
        0.5 0.5 0.5];% Gray

% --- Transpose the matrices ---
% When displaying the transposed grid, rows and columns switch.
displayGridT = displayGrid';
gridCountsT = gridCounts';

% Create the figure and plot the transposed grid
figure;
imagesc(displayGridT);
colormap(cmap);
set(gca, 'YDir', 'normal');  % (1,1) will now appear at the lower left corner
axis equal tight;
colorbar off;

% Set ticks and labels
xticks(1:10);
yticks(1:10);
set(gca, 'XTickLabel', 1:10, 'YTickLabel', 1:10);
grid on;
set(gca, 'GridColor', 'k', 'GridAlpha', 1);

% Add text inside each cell for nonempty positions (transposed)
for row = 1:10
    for col = 1:10
        if gridCountsT(row, col) > 0
            text(col, row, num2str(gridCountsT(row, col)), 'HorizontalAlignment', 'center', ...
                 'Color', 'k', 'FontWeight', 'bold', 'FontSize', 12);
        end
    end
end

title('Sensor Distribution on 10x10 Grid');

hold on;
hEmpty = patch(NaN, NaN, cmap(1,:));
hOne   = patch(NaN, NaN, cmap(2,:));
hMulti = patch(NaN, NaN, cmap(3,:));
hold off;

legend([hEmpty, hOne, hMulti], {'Empty Cell', '1 Sensor in the cell', '2 Sensors in the cell'}, ...
       'Location', 'best','FontSize',10);

%% attack-free finger print localization with 1 sensor

for i = 1:20
    D_sensors = D(i,:);
    y_sensors = y(i);
    
    % vecnorm(vector, norm type, dimension)
    [~,position] = min(vecnorm(D_sensors - y_sensors, 2, 1));

    [i position]
   
end

% with two sensors
D_sensors = D(2:3,:);
y_sensors = y(2:3);

% vecnorm(vector, norm type, dimension)
[~,position] = min(vecnorm(D_sensors - y_sensors, 2, 1))
