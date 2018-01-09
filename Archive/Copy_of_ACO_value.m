%This code is for a simple VRP with normal constraints. Inputs are:
%heuristic matrix (nxn), demand array (1xn), time window array (2xn),
%service time array (1xn), max vehicle load, heuristic importance,
%pheromone importance, the number of iterations, and allotted run time.
%(Note that only of the last two is necessary for the code to run, so put
%in a non-negative value for one you'd like to use and zero for the other.
%If you put in non-negative values for both, the algorithm will default to
%run time).

%Inputs: heuristic = matrix of distances/costs between nodes (symmetric)
%demand = array of demands for each node tw = "2xn" array with beginning
%(1st row) and end (2nd row) of time window service = array of service
%times for each node Q = max vehicle load alpha = relative importance of
%heuristic information beta = relative importance of pheromone information
%ro = pheromone evaporation rate (now calculated) min = min pheromone level
%max = max pheromone level

%Potential inputs not currently included q0 = acceptance criteria for city
%selection (Ant System modification) ants = number of ants (to include just
%introduce as input and use as
%       index for outermost for loop

function [full_solution]=ACO_value(d_dist, p, demand, tw, service, Q, beta, run_time, pbest, qo, rho, gp_ratio, num_ants, phi, dist_limit, veh_speed, threat, num_vehicles, theta, current_inventory)

%Start the timers.
tic;
timer_start = cputime;

%p_dist is probabilities of destruction,  just to match the minimization
%structure of the problem. p_dist2 are probabilities of survival.  All are
%in log format so they're additive.
p_dist = -log(1-p);
p_dist2 = -log(p);

%define heuristic2 as distances, heuristic as times, and hx as the number
%of nodes. Heuristic3 is the probabilities of survival on a route.
heuristic3 = p_dist;
heuristic2 = d_dist;
heuristic = d_dist./veh_speed;
hx=size(heuristic,1);

%__________________________________________________________________________
%Test parameters
%pbest
%qo
%rho
%beta
%gp_ratio
%num_ants
%__________________________________________________________________________
%__________________________________________________________________________
%Constants
alpha = 1;
iterations = 15; %50
%cl_size = ceil(hx/5);
cl_size = ceil(hx/2);
%__________________________________________________________________________

%Calculate rho based on number of iterations.  If iteration count isn't
%given, estimate based on allotted run_time.
% if iterations>0
%     sbar = iterations;
% else
%     sbar = run_time*5;
% end

%Stutzle's rule
%rho = ((1-(.05)^(1/(hx-1)))/(((hx-1)/2-1)*(.05)^(1/(hx-1))))^(1/sbar);

%Form the candidate lists cl_size = round(hx/5);
candidate_list = zeros(hx,cl_size);

for i = 1:hx
    candidate_list(1,i)         = i;
end

for i = 2:hx
    [~,temp]                    = sort(heuristic(i,2:hx));
    candidate_list(i,1:cl_size) = temp(2:(cl_size + 1)) + 1;
end

%Initialize global best solution, cost, and array of all costs
sbest                               = [];
costbest                            = sum(sum(heuristic));
antcostcounter                      = 0;
antcosts                            = zeros(1,max((ceil(sqrt(hx-1))*iterations),1));
bestantcostcounter                  = 0;
bestantcosts                        = zeros(1,max(iterations,1));
lscostcounter                       = 0;
lscosts                             = zeros(1,max(iterations,1));
loads_bestofiterationroute          = [];
loads_best                          = [];
start_times_bestofiterationroute    = [];
start_times_best                    = [];
route_cost_bestofiterationroute     = [];
ant_time_index                      = 1;
all_ant_times                       = zeros(50,1);
LS_time_index                       = 1;
all_LS_times                        = zeros(50,1);
route_value_best                    = [];
s_value_best                        = [];
s_value_bestofiteration             = [];
route_value_bestofiteration         = [];
%route_value = [];
route_dist_bestofiterationroute     = [];
all_value_bestofiteration           = [];


%Run a single iteration with a beta of zero, which is equivalent to a
%probabilistic nearest neighbor solution.
[~, cost1, ~] = ant_value(heuristic, heuristic2, heuristic3, candidate_list, cl_size, ones(hx), demand, tw, service, Q, alpha, 0, rho, 0, 1, hx, qo, phi, dist_limit, threat, num_vehicles, theta, current_inventory);

%Use this initial solution to estimate the maximum pheromone amount.
mx = 1/(1-rho)*1/((hx-1)*cost1);

%Using Stutzle's rule, construct min
pdec = pbest^(1/(hx-2));
mn = (2*mx*(1-pdec))/((hx-1)*pdec);

%Form pheromone matrix, setting all pheromones to max value.
pheromone = ones(hx).*mx;

%Define my limit as either time or iterations
if run_time > 0
    limit   = run_time;
    counter = toc;
else
    limit   = iterations;
    counter = 0;
end

iteration_counter = 0;
pheromone_counter = 0;

%Main ACO function
while counter < limit

    %start timer for this iteration of ants
    ant_timer_start   = cputime;
    
    iteration_counter = iteration_counter + 1;
    pheromone_counter = pheromone_counter + 1;
    
    %Set iteration best cost to big M, in this case 10x the probabilistic
    %nearest neighbor solution calculated above.
    bestofiterationcost = cost1*10;
    bestofiterationsolution = [];
    
    %# of ants is equal to the number of customers
    for i = 1:num_ants
        
        %Subrouting constructs an entire solution for 1 ant
        [s, cost, route_cost, pheromone, start_times, loads, route_dist] = ant_value(heuristic, heuristic2, heuristic3, candidate_list, cl_size, pheromone, demand, tw, service, Q, alpha, beta, rho, mn, mx, hx, qo, phi, dist_limit, threat, num_vehicles, theta, current_inventory);
    
        %Calculate the values of the given routes.
        all_value = route_values(s, loads, current_inventory, theta);
        if num_vehicles < size(all_value,1)
            route_value = all_value(1:num_vehicles,:);
        else
            route_value = all_value;
        end

        if size(s,2) > size(loads,2)
            temp = zeros(size(loads,1),(size(s,2)-size(loads,2)));
            loads = [loads,temp];
        end
        
        if size(route_value,1) < num_vehicles
            temp = zeros(num_vehicles - size(route_value,1),size(route_value,2));
            route_value = [route_value;temp];
        end
        
        %Calculate the total value of a solution.
        s_value = sum(sum(route_value));
        
        %Calculate probabilities.
        %Initialize to 0, then initialize the first column to 1. We're assuming
        %the base doesn't get bombed, so probability of survival is 1 prior to
        %departure.
        prob = zeros(size(route_cost_bestofiterationroute));
        %Compute route survival probabilities.
        for i31 = 1:size(s,1)
            for i32 = 2:size(s,2)
                %Check to see if we're at the end of a route.
                if s(i31,i32) == 0
                    break;
                end

                %If not, add probability of survival of reaching that node
                %(i.e., probability of survival associated with previous leg).
                prob(i31,i32) = heuristic3(s(i31,i32-1),s(i31,i32));
            end
        end
        
        %If the first run, then store the solution as the best.
        if i == 1
            %Store as best of iteration
            s_value_bestofiteration = s_value;
            route_value_bestofiteration = route_value;
            bestofiterationcost = cost;
            bestofiterationsolution = s;
            start_times_bestofiterationroute = start_times;
            loads_bestofiterationroute = loads;
            route_cost_bestofiterationroute = route_cost;
            route_dist_bestofiterationroute = route_dist;
            all_value_bestofiteration = all_value;
            prob_bestofiteration = prob;


            
            if counter == 0
                s_value_best = s_value;
                route_value_best = route_value;
                costbest = cost;
                sbest = s;
                start_times_best = start_times;
                loads_best = loads;
                %route_cost_best = route_cost;
                route_dist_best = route_dist;
                %all_value_best = all_value;
                prob_best = prob;
            end
        end
            
       
        %Update global best
        if s_value > s_value_best
            s_value_best = s_value;
            route_value_best = route_value;
            costbest = cost;
            sbest = s;
            start_times_best = start_times;
            loads_best = loads;
%             route_cost_best = route_cost;
            route_dist_best = route_dist;
%             all_value_best = all_value;
            prob_best = prob;
        end

        %Store current cost in array
        antcostcounter = antcostcounter + 1;
        antcosts(antcostcounter) = cost;
        
        if s_value > s_value_bestofiteration
            s_value_bestofiteration = s_value;
            route_value_bestofiteration = route_value;
            bestofiterationcost = cost;
            bestofiterationsolution = s;
            start_times_bestofiterationroute = start_times;
            loads_bestofiterationroute = loads;
            route_cost_bestofiterationroute = route_cost;
            route_dist_bestofiterationroute = route_dist;
            all_value_bestofiteration = all_value;
        end
        
    end
%     
%     bestantcostcounter = bestantcostcounter + 1;
%     bestantcosts(bestantcostcounter) = bestofiterationcost;
%     
%     ant_time = cputime - ant_timer_start;
%     all_ant_times(ant_time_index) = ant_time;
%     ant_time_index = ant_time_index + 1;

%__________________________________________________________________________    
%Local Search based on value function

    feeder_solution = bestofiterationsolution;
    feeder_cost = bestofiterationcost;
    feeder_route_cost = route_cost_bestofiterationroute;
    feeder_start_times = start_times_bestofiterationroute;
    feeder_loads = loads_bestofiterationroute;
    feeder_route_dist = route_dist_bestofiterationroute;
    feeder_s_value = s_value_bestofiteration;
    feeder_route_value = route_value_bestofiteration;
    feeder_all_value = all_value_bestofiteration;
    lowest_num_veh = size(feeder_solution,1);

%     vale = min(num_vehicles, size(feeder_solution,1));
%     if vale > 1
%         change1 = 1;
%         change2 = 1;
%     elseif (size(feeder_solution,1)-num_vehicles) <= 0
%         change1 = 1;
%         change2 = 1;        
%     else
%         change1 = 0;
%         change2 = 0;
%     end
    change3 = 1;
    trigger = 1;

    while trigger == 1
        
        trigger = 0;
% change1=0;%temporarily disable 2-opt*

        
        vale = min(num_vehicles, size(feeder_solution,1));
        if (vale > 1) && (size(feeder_solution,1) > num_vehicles)
            change1 = 1;
            change2 = 1;
        else %(size(feeder_solution,1)-num_vehicles) <= 0
            change1 = 0;
            change2 = 0;        
        end

        while change1 == 1
            
            if (size(feeder_solution,1) - num_vehicles) <= 0
                break;
            end
            
            [feeder_solution, feeder_cost, feeder_route_cost, feeder_start_times, feeder_loads, feeder_route_dist, feeder_s_value, feeder_route_value, feeder_all_value, change1] = two_opt_star_value(feeder_solution, feeder_cost, feeder_route_cost, heuristic, heuristic2, tw, service, Q, feeder_start_times, feeder_loads, dist_limit, feeder_route_dist, num_vehicles, feeder_s_value, feeder_route_value, feeder_all_value, current_inventory, theta);
            if change1 == 1
                trigger = 1;
                change2 = 1;
                change3 = 1;
            end
        end
% change2=0;%temporarily disable Or-opt        
        while change2 == 1

            if (size(feeder_solution,1) - num_vehicles) <= 0
                break;
            end
            
            [feeder_solution, feeder_cost, feeder_route_cost, feeder_start_times, feeder_loads, feeder_route_dist, feeder_s_value, feeder_route_value, feeder_all_value, change2] = Or_opt_v2_value(feeder_solution, feeder_cost, feeder_route_cost, heuristic, heuristic2, tw, service, Q, feeder_start_times, feeder_loads, dist_limit, num_vehicles, feeder_route_dist, feeder_s_value, feeder_route_value, feeder_all_value);
            if change2 == 1
                trigger = 1;
                change1 = 1;
                change3 = 1;
            end
        end

        while change3 == 1

            [feeder_solution, feeder_cost, feeder_route_cost, feeder_start_times, feeder_loads, feeder_route_value, feeder_all_value, change3] = relocate_value(feeder_solution, feeder_cost, feeder_route_cost, heuristic, tw, service, Q, feeder_start_times, feeder_loads, feeder_route_value, feeder_all_value, num_vehicles);
            if change3 == 1
                if num_vehicles > 1
                    trigger = 1;
                    change1 = 1;
                    change2 = 1;
                end
            end
        end
    end

%     feeder_solution2 = feeder_solution;
%     feeder_cost2 = feeder_cost;
%     feeder_start_times2 = feeder_start_times;
%     feeder_loads2 = feeder_loads;
%     feeder_route_dist2 = feeder_route_dist;
    
%     lscostcounter = lscostcounter + 1;
%     lscosts(lscostcounter) = ls_cost;
    
    
    %Update global best
    if (feeder_s_value > s_value_best) || ((feeder_s_value == s_value_best) && (size(feeder_all_value,1) < lowest_num_veh))
        s_value_best = feeder_s_value;
        route_value_best = feeder_route_value;
        costbest = feeder_cost;
        sbest = feeder_solution;
        start_times_best = feeder_start_times;
        loads_best = feeder_loads;
%        route_cost_best = feeder_route_cost;
        route_dist_best = feeder_route_dist;
%        all_value_best = feeder_all_value;
    end
    
%__________________________________________________________________________
% LS_timer_start = cputime;

%Local Search to optimize distances, run only on the appropriate number of
%routes.

%     feeder_solution2 = bestofiterationsolution;
%     feeder_cost2 = bestofiterationcost;
%     feeder_route_cost2 = route_cost_bestofiterationroute;
%     feeder_start_times2 = start_times_bestofiterationroute;
%     feeder_loads2 = loads_bestofiterationroute;
%     feeder_route_dist2 = route_dist_bestofiterationroute;

    %Trim and recalculate as appropriate.
    trim = min(num_vehicles, size(feeder_solution,1));
    feeder_solution = feeder_solution(1:trim,:);
    
    trim = min(num_vehicles, size(feeder_route_cost,1));
    feeder_route_cost = feeder_route_cost(1:trim,:);
    feeder_cost = sum(feeder_route_cost);
    
    trim = min(num_vehicles, size(feeder_start_times,1));
    feeder_start_times = feeder_start_times(1:trim,:);
    
    trim = min(num_vehicles, size(feeder_loads,1));
    feeder_loads = feeder_loads(1:trim,:);

    trim = min(num_vehicles, size(feeder_route_dist,1));
    feeder_route_dist = feeder_route_dist(1:trim,:);
    
    %Initialize to 0, then initialize the first column to 1. We're assuming
    %the base doesn't get bombed, so probability of survival is 1 prior to 
    %departure.
    feeder_prob = zeros(size(feeder_route_cost));
    %Compute route survival probabilities.
    for i31 = 1:size(feeder_solution,1)
        for i32 = 2:size(feeder_solution,2)
            %Check to see if we're at the end of a route.
            if feeder_solution(i31,i32) == 0
                break;
            end
            
            %If not, add probability of survival of reaching that node
            %(i.e., probability of survival associated with previous leg).
            feeder_prob(i31,i32) = heuristic3(feeder_solution(i31,i32-1),feeder_solution(i31,i32));
        end
    end

    change1 = 1;
    change2 = 1;
    change3 = 1;
    trigger = 1;

    while trigger == 1
 
        trigger = 0;
      
        while change1 == 1
            
            [feeder_solution, feeder_cost, feeder_prob, feeder_start_times, feeder_loads, feeder_route_dist, change1] = two_opt_star(feeder_solution, feeder_cost, feeder_prob, heuristic, heuristic2, heuristic3, tw, service, Q, feeder_start_times, feeder_loads, dist_limit, feeder_route_dist);
            if change1 == 1
                trigger = 1;
                change2 = 1;
                change3 = 1;
            end
        end

        while change2 == 1
            
            [feeder_solution, feeder_cost, feeder_prob, feeder_start_times, feeder_loads, feeder_route_dist, change2] = Or_opt_v2(feeder_solution, feeder_cost, feeder_prob, heuristic, heuristic2, heuristic3, tw, service, Q, feeder_start_times, feeder_loads, dist_limit, feeder_route_dist);
            if change2 == 1
                trigger = 1;
                change1 = 1;
                change3 = 1;
            end
        end

        while change3 == 1

            [feeder_solution, feeder_cost, feeder_prob, feeder_start_times, feeder_loads, change3] = relocate(feeder_solution, feeder_cost, feeder_prob, heuristic, heuristic3, tw, service, Q, feeder_start_times, feeder_loads);
            if change3 == 1
                trigger = 1;
                change1 = 1;
                change2 = 1;
            end
        end
    end

    %Recalculate the route_value. Total value is unchanged but customers
    %(and their associated values) were shifted around.
    route_value = route_values(feeder_solution, feeder_loads, current_inventory, theta);
    vale = min(size(route_value,1),num_vehicles);
    route_value = route_value(1:vale,:);
    
    ls_solution = feeder_solution;
    ls_cost = feeder_cost;
    ls_start_times = feeder_start_times;
    ls_loads = feeder_loads;
    ls_dist = feeder_route_dist;
    ls_value = feeder_s_value;
    ls_route_value = route_value;
    ls_prob = feeder_prob;
    
%     lscostcounter = lscostcounter + 1;
%     lscosts(lscostcounter) = ls_cost;
    
%     LS_time = cputime - LS_timer_start;
%     all_LS_times(LS_time_index) = LS_time;
%     LS_time_index = LS_time_index + 1;
%__________________________________________________________________________
    
    %Update global best
%     if ls_value > s_value_best
        costbest = ls_cost;
        route_dist_best = ls_dist;
        sbest = ls_solution;
        start_times_best = ls_start_times;
        loads_best = ls_loads;
        s_value_best = ls_value;
        route_value_best = ls_route_value;
        prob_best = ls_prob;
%     end
    
    %Local pheromone update.  Decrease all by rho.
    pheromone = pheromone.*(1-rho);
    
    %Check against min values.
    for x1 = size(pheromone,1)
        for x2 = size(pheromone,2)
            if pheromone(x1,x2) < mn
                pheromone(x1,x2) = mn;
            end
        end
    end
    
        
    %Global pheromone update using global best
    for k=1:size(sbest,1)-1
        for j=1:size(sbest,2)-1
            
            %Only update if the entry is an actual path in a route
            if (sbest(k,j) ~= 0) && (sbest(k,j+1) ~= 0)
                pheromone(sbest(k,j),sbest(k,j+1)) = pheromone(sbest(k,j),sbest(k,j+1)) + mx/3*gp_ratio;%(1/costbest)*gp_ratio;
                
                %Check pheromones against min/max values.  Only need to
                %check against max because we're just adding positive
                %pheromone at this step.
%                 if pheromone(sbest(k,j),sbest(k,j+1)) < mn
%                     pheromone(sbest(k,j),sbest(k,j+1)) = mn;
%                 elseif pheromone(sbest(k,j),sbest(k,j+1)) > mx
                if pheromone(sbest(k,j),sbest(k,j+1)) > mx
                    pheromone(sbest(k,j),sbest(k,j+1)) = mx;
                end
               
                %Maintain symmetry
                pheromone(sbest(k,j+1),sbest(k,j)) = pheromone(sbest(k,j),sbest(k,j+1));
                
            end
        end
    
    end
    
    %Global pheromone update using iteration best
    for k=1:size(ls_solution,1)-1
        for j=1:size(ls_solution,2)-1
            
            %Only update if the entry is an actual path in a route
            if (ls_solution(k,j) ~= 0) && (ls_solution(k,j+1) ~= 0)
                pheromone(ls_solution(k,j),ls_solution(k,j+1)) = pheromone(ls_solution(k,j),ls_solution(k,j+1)) + mx/3*gp_ratio;%(1/ls_cost)*(1-gp_ratio);
                
                %Check pheromones against min/max values
%                 if pheromone(bestofiterationsolution(k,j),bestofiterationsolution(k,j+1)) < mn
%                     pheromone(bestofiterationsolution(k,j),bestofiterationsolution(k,j+1)) = mn;
%                 elseif pheromone(bestofiterationsolution(k,j),bestofiterationsolution(k,j+1)) > mx
                if pheromone(ls_solution(k,j),ls_solution(k,j+1)) > mx
                    pheromone(ls_solution(k,j),ls_solution(k,j+1)) = mx;
                end
 
                %Maintain symmetry
                pheromone(ls_solution(k,j+1),ls_solution(k,j)) = pheromone(ls_solution(k,j),ls_solution(k,j+1));

            end
        end
    end
    
    
    if run_time > 0
        counter = toc;
    else
        counter=counter+1;
    end
end


%Fix sbest to match Ian's numbering (depot is last, not first).
for j = 1:size(sbest,2)
    for i = 1:size(sbest,1)
        if sbest(i,j) == 1
            sbest(i,j) = hx;
        elseif sbest(i,j) > 1
            sbest(i,j) = sbest(i,j) - 1;
        end
    end
end

%Compute the cumulative probabilities.
cumul_prob = zeros(size(sbest));
for i = 1:size(sbest,1)
    for j = 2:size(sbest,2)
        if sbest(i,j) == 0
            break;
        end
        
        cumul_prob(i,j) = cumul_prob(i,j-1) + p_dist2(sbest(i,j-1),sbest(i,j));
    end
end

%Convert probabilities back from log format to actual probabilities.
prob_best = 1-exp(-prob_best);
prob_best(:,1) = 1;
cumul_prob = exp(-cumul_prob);


%Trim all of the parts of the solution to make sure they're the same size,
%meaning delete any dummy columns.
%Sbest
x = size(sbest,2);
y = sum(sbest,1);
for k = x:-1:1
    if y(k) == 0
        sbest(:,k) = [];
    else
        break;
    end
end

%Start_times
y = size(start_times_best,2);
if y > x
    for k = y:-1:x+1
        start_times_best(:,k) = [];
    end
end

%loads
y = size(loads_best,2);
if y > x
    for k = y:-1:x+1
        loads_best(:,k) = [];
    end
end

%route_value
y = size(route_value_best,2);
if y > x
    for k = y:-1:x+1
        route_value_best(:,k) = [];
    end
end

%prob
y = size(prob_best,2);
if y > x
    for k = y:-1:x+1
        prob_best(:,k) = [];
    end
end

%cumul_prob
y = size(cumul_prob,2);
if y > x
    for k = y:-1:x+1
        cumul_prob(:,k) = [];
    end
end

%Output solution
full_solution{1} = sbest;
full_solution{2} = sum(sum(route_dist));%This is the total distance required to fulfill the solution.
full_solution{3} = start_times_best;%This will give you the total time for each route. The last non-zero value in each row is the time the vehicle arrives back at the depot.
full_solution{4} = loads_best;%This is the load amount for each delivery.
% full_solution{5} = antcosts(1:antcostcounter);%Disregard
% full_solution{6} = bestantcosts(1:bestantcostcounter);%Disregard
% full_solution{7} = lscosts(1:lscostcounter);%Disregard
full_solution{5} = toc;%This is the total time it took to generate all of this.
full_solution{6} = cputime - timer_start;%Same as above, just calculated slightly differently. Should be within a couple seconds of each other.
% full_solution{10} = all_ant_times;%Disregard
% full_solution{11} = all_LS_times;%Disregard
full_solution{7} = route_dist_best;%This is the distance covered by each route.
full_solution{8} = s_value_best;%This is the total value added by executing this solution.
full_solution{9} = route_value_best;%These are the individual values of each delivery.
full_solution{10} = prob_best;%These are the probabilities of survival for each leg of a route.
full_solution{11} = cumul_prob;%These are the cumulative probabilites for each route.

end
