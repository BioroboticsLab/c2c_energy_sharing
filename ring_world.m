 function [M, S] = ring_world(N)

G.stepwidth       = 4;              % acceleration factor
G.motion_dist_per_step = 0.001;     % in the ring world the max distance is 2  
G.simsteps        = 10000;          %
G.consumption_per_timestep = 0.001; % 
G.share_energy    = 1;
G.plot_world      = 1;
G.iterations      = 1;
G.charge_threshold_full         = 3;        % 
G.charge_from_grid_per_timestep = 0.004;    % 4 time faster charging than spending
G.charge_from_car_per_timestep  = 0.012;    % using the current scheme, half of this is going to be shared at max
G.proportion_taxi               = 0.3;      % how many taxis (never charge on grid)
G.Ntxi                          = round(N * G.proportion_taxi);

G.M = nan(1,N);
G.W = nan(1,N);

S = [];

for iter = 1 : G.iterations
fprintf('\r%d\n', iter);

G = initialize_world(N, G);

% without sharing energy some cars might not reach their first goal
G.energy_left(G.energy_left < 2) = 2; % 

% fprintf('mean energy: %.2f\n', mean(G.energy_left));

    for t = 1 : G.simsteps
        G.t = t;
        % those that arrived, do not need to go (square not circ. region)
        G.lgc_on_grid       = sum( abs(G.positions - G.targets) < 0.01, 2) == 2;
        G.lgc_in_field      = ~G.lgc_on_grid;
        
        % all vehicles that arrived are charged and marked for departure
        % when fully charged 
        G = charge_from_grid(G);
        
        G = determine_who_can_go(G);
        
        G = set_new_targets((G.lgc_on_grid & G.lgc_fully_charged) | (G.lgc_on_grid & G.lgc_taxis), G);      
        
        G.M(t) = sum(G.lgc_can_go) / N;
        fprintf('\b\b\b\b%1.2f',G.M(t));
        
        % move one step forward
        G.positions(G.lgc_can_go,:)   = G.positions(G.lgc_can_go,:) + G.stepwidth * G.motion_dist_per_step * G.directions(G.lgc_can_go,:);

        % step consumes energy
        G.energy_left(G.lgc_can_go) = G.energy_left(G.lgc_can_go) - G.stepwidth * G.consumption_per_timestep;
        G.energy_left(G.energy_left < 0) = 0;

        % determine who is close to whom
        % matrix of booleans
        MP = nearby_pairs(G.positions, 0.05);
        % array with ones denoting the index of interacting agents
        I = find(sum(MP, 2) > 0);

        if G.share_energy
            % for each interaction group
            for i = I'    
                J = [i find(MP(i, :))]; % in J there are all interaction partner indices

                % give proportional to what you have (one value for each car)
                E = G.energy_left(J); 

                % scale the energy amount with how far you go per timestep
                % (stepwidth)
                energy_to_give = G.stepwidth * G.charge_from_car_per_timestep * E / sum(E); 

                % take from everyone
                G.energy_left(J) = G.energy_left(J) - energy_to_give;

                % get mean from everyone
                G.energy_left(J) = G.energy_left(J) + mean(energy_to_give);     
            end
        end
            

        % paint cars and world
        if G.plot_world
            plot_step_results(G, I);
        end   
    end
        
%    fprintf('mean energy after: %.2f\n', mean(energy_left));
    
%     directions  = targets - positions;
%     lengths     = sqrt(sum(directions.^2, 2));
%     M = [M; mean(lengths)];
%     S = [S; std(lengths)];
%     
    
end

function plot_step_results(G, I)
subplot(4, 2, [5:8])
plot( G.positions(:,1), G.positions(:,2), '*')
hold on
plot( G.positions(G.lgc_taxis,1), G.positions(G.lgc_taxis,2), 'om')   % taxis
plot( G.positions(I,1), G.positions(I,2), 'y*')                 % potential sharers
plot( G.positions(~G.lgc_can_go,1),  G.positions(~G.lgc_can_go,2), 'r*')% emtpy batteries    
plot( G.positions( G.lgc_on_grid,1), G.positions( G.lgc_on_grid,2), 'g*')% arrived and charging from grid
hold off
axis equal
axis([-1 1 -1 1])
subplot(4, 2, 1)
title('fraction of active cars')
plot(G.M(1:G.t))

subplot(4, 2, 2)
title('energy distribution')
hist(G.energy_left)

subplot(4, 2, 3)
title('current avg charging time on grid')
plot(G.W(1:G.t), '.')

drawnow

function I = nearby_pairs(P, threshold)
I = squareform(pdist(P) < threshold);

function [G] = initialize_world(N, G)
G.t  = 0;
% start and end locations on the circle (in radians)
angular_positions   = rand(N, 1) * 2*pi; % arrange cars on a circle
angular_targets     = rand(N, 1) * 2*pi;
G.energy_left         = rand(N, 1) * 4;  
% make cartesian coordinates 
G.positions         = [cos(angular_positions) sin(angular_positions)];
G.targets           = [cos(angular_targets) sin(angular_targets)];
% direction vectors, all unit lengths
G.directions        = G.targets - G.positions;
lengths             = sqrt(sum(G.directions.^2, 2));
G.directions        = G.directions./repmat(lengths, 1, 2);
G.lgc_on_grid       = zeros(N, 1) == 1;
G.grid_charging_time = zeros(N, 1);
G.lgc_taxis         = (1:N)' <= G.Ntxi;
G.lgc_fully_charged = zeros(N,1) == 1;

function [G] = charge_from_grid(G)
% add energy to battery of those that arrived (technically they are arrived
% as long as they have not fully charged)
G.energy_left(G.lgc_on_grid)    = G.energy_left(G.lgc_on_grid) + G.stepwidth * G.charge_from_grid_per_timestep;

% count time from arrival 
G.grid_charging_time(G.lgc_on_grid) = G.grid_charging_time(G.lgc_on_grid) + G.stepwidth;

% mark the full ones
G.lgc_fully_charged(G.lgc_on_grid)	= G.energy_left(G.lgc_on_grid) > G.charge_threshold_full; 

% collect statistics
G.W(G.t) = mean(G.grid_charging_time(G.lgc_fully_charged));

% set back to zero for the ones that can go
G.grid_charging_time(G.lgc_fully_charged)    = 0;


function G = determine_who_can_go(G)
G.lgc_can_go = (G.lgc_taxis & G.energy_left>0) | G.lgc_fully_charged | ( G.lgc_in_field & (G.energy_left > 0));

function G = set_new_targets(lgc_idx, G)
angular_targets     = rand(sum(lgc_idx), 1) * 2*pi;
G.targets(lgc_idx,:)    = [cos(angular_targets) sin(angular_targets)];
G.directions(lgc_idx,:) = G.targets(lgc_idx,:) - G.positions(lgc_idx,:);
lengths             = sqrt(sum(G.directions(lgc_idx,:).^2, 2));
G.directions(lgc_idx,:) = G.directions(lgc_idx,:)./repmat(lengths, 1, 2);
    