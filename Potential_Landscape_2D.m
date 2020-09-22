%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%-----Equations describing the dynamical system-----%

% dxdt = c1 + a1*(x^n)/(Kdxx^n + (x^n)) + (b1*(Kdyx^n))/(Kdyx^n + (y^n)) - (x*k1)            
% dydt = c1 + a2*(y^n)/(Kdyy^n + (y^n)) + (b2*(Kdxy^n))/(Kdxy^n + (x^n)) - (y*k2)

%NOTE***: All inputs are at the bottom of the code
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all;
tic;  %track run time

[x_lower_lim, y_lower_lim, x_upper_lim, y_upper_lim, xyGridSpacing, numTimeSteps, dt, tol] = program_parameters();
[grid_lines] = graph_parameters();

% Calculate total no. of paths for defined grid spacing
numPaths = 0;
for i = x_lower_lim : xyGridSpacing : x_upper_lim
    for j = y_lower_lim : xyGridSpacing : y_upper_lim
        numPaths = numPaths + 1;
    end
end
 
% Initialize "path" variable matrices
x_path = zeros(numPaths,numTimeSteps);  %x-coord. along path
y_path = zeros(numPaths,numTimeSteps);  %y-coord. along path
pot_path = zeros(numPaths,numTimeSteps);    %pot. along path
 
path_tag = ones(numPaths,1);   %tag denotes the basin of attraction of each path
                               % ** initialized to 1 for all paths **
                               
% Initialize "Path counter" to 1
path_counter = 1;   
 
% Initialize no. of attractors and separatrices (basin boundaries)
num_attractors = 0;  

% Assign array to keep track of attractors and their coordinates; and pot.
attractors_num_X_Y = [];

% Assign array to keep track of no. of paths per attractor
numPaths_att = [];
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Loop over x-y grid
for i = x_lower_lim : xyGridSpacing : x_upper_lim
    for j = y_lower_lim : xyGridSpacing : y_upper_lim
       
        %Initialize coords. (Starting points)
        x0 = i;    
        y0 = j;
 
        
        p0 = 0;   %Initialized to zero for comparison of potential                 
        
        %Initialize "path" variables 
        x_p = x0;
        y_p = y0;
        
        %Initialize accumulators for "potential" along path
        Pot = p0;
        Pot_old = 1.0e7;   %initialize to large number
 
        %Initialize global arrays (time t = 0 counts as "time step #1")
        x_path(path_counter, 1) = x_p;
        y_path(path_counter, 1) = y_p;        
        pot_path(path_counter, 1) = Pot;
        
        %Evaluate potential (Integrate) over trajectory from init cond to  stable steady state
        for n_steps = 2:numTimeSteps
            
            Pot_old = Pot;
    
            %calculation of the time derivatives of the variables
            
            dxdt = calc_first_variable(x_p, y_p);            
            dydt = calc_second_variable(x_p, y_p);
            
                      
            dx = (dxdt)*dt;
            dy = (dydt)*dt;
            
            %update x, y 
            x_p = x_p + dx;
            y_p = y_p + dy;
            
            x_path(path_counter, n_steps) = x_p;
            y_path(path_counter, n_steps) = y_p;
            
            %update "potential" 
            dPot = - (dxdt)*dx - (dydt)*dy;   % signs ensure that "potential" decreases as "velocity" increases
            Pot = Pot_old + dPot; 
            pot_path(path_counter, n_steps) = Pot;
                       
        end    % end integration over path
              
        %check for convergence
        if (abs(Pot - Pot_old) > tol)
            fprintf(1,'Warning: not converged!\n');        
        end
 
        % --- assign path tag to track the basins of attraction of each
        % path
        if (path_counter == 1) %record attractor of first path and its coords 
            num_attractors = num_attractors + 1;
            current_att_num_X_Y = [num_attractors  x_p  y_p  Pot  path_counter]; %create array
            attractors_num_X_Y = [attractors_num_X_Y; current_att_num_X_Y]; % appending vertically
            path_tag(path_counter) = num_attractors;  %initialize path tag 
            numPaths_att = [numPaths_att;  1]; % appending vertically
            
        else    % i.e. if path counter > 1
 
            xp_lastPath = x_path((path_counter - 1), numTimeSteps);
            yp_lastPath = y_path((path_counter - 1), numTimeSteps);
            
            %calculate distance between "end points" of current and previous paths
            endPt_dist_sqr = ((x_p - xp_lastPath)^2 + (y_p - yp_lastPath)^2);    %This is for identifying attractors               
  
            %check if the current path *ended* in a different point compared to previous path
            % (x-y grid spacing used as a "tolerance" for distance)          
            if ( endPt_dist_sqr > (2 * (xyGridSpacing^2)) )  
                
                % --- check if this "different" attractor has been identified before
                new_attr_found = 1;
                
                for k = 1 : num_attractors
                    x_att = attractors_num_X_Y(k,2);
                    y_att = attractors_num_X_Y(k,3);                   
                    if ( (abs(x_p - x_att) < xyGridSpacing) && (abs(y_p - y_att) < xyGridSpacing) )
                        % this attractor has been identified before
                        new_attr_found = 0;
                        path_tag(path_counter) = k;  
                        numPaths_att(k) =  numPaths_att(k) + 1;
                        break;  %exit for-loop
                    end                  
                end
                
                if (new_attr_found == 1)
                    num_attractors = num_attractors + 1;
                    current_att_num_X_Y = [num_attractors  x_p  y_p  Pot  path_counter];%create array
                    attractors_num_X_Y = [attractors_num_X_Y; current_att_num_X_Y];%appended vertically
                    path_tag(path_counter) = num_attractors;  
                    numPaths_att = [numPaths_att; 1]; %Appended array vertically                                              
                    
                end
            
            else  %i.e. current path converged at same pt. as previous path
                
                %%update path tag 
                path_tag(path_counter) = path_tag(path_counter - 1);
                
                %update no. of paths for current attractor
                % (path tag already updated at start of path-counter loop)
                tag = path_tag(path_counter);
                numPaths_att(tag) =  numPaths_att(tag) + 1;
                
            end   % end check for location of path end-pt.  
            
        end
            
        % increment "path counter"
        path_counter = path_counter + 1;                     
        
    end
end
 
fprintf(1,'  Ran path-loop okay!\n'); 
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -- need 1-D "lists" (vectors) to plot all x,y, Pot values along paths --
list_size = numPaths*numTimeSteps; 
 
x_p_list = zeros(list_size,1);
y_p_list = zeros(list_size,1);
pot_p_list = zeros(list_size,1);
 
n_list = 1;
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Normalising potential values to "Align" potential values so all path-potentials end up at same global min.
for n_path = 1:numPaths
    tag = path_tag(n_path);
    del_pot = pot_path(n_path, numTimeSteps) - attractors_num_X_Y(tag, 4);
    
    %noramlise pot. at each time step along path
    for n_steps = 1:numTimeSteps
        pot_old = pot_path(n_path, n_steps);    
        pot_path(n_path, n_steps) = pot_old - del_pot;
        
        %add data point to list
        x_p_list(n_list) = x_path(n_path, n_steps);
        y_p_list(n_list) = y_path(n_path, n_steps);
        pot_p_list(n_list) = pot_path(n_path, n_steps);
        
        n_list = n_list + 1;    % increment n_list
 
    end 
    
end
 
fprintf(1,'  Ran path-alignment okay!\n');  
fprintf(1,'  ************************\n'); 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate surface interpolation grid 
%--- Create X,Y grid to interpolate "potential surface" ---
%No. of grid lines in x- and y- directions
xlin = linspace(min(x_p_list), max(x_p_list), grid_lines);
ylin = linspace(min(y_p_list), max(y_p_list), grid_lines);
[Xgrid,Ygrid] = meshgrid(xlin,ylin);
Zgrid = griddata(x_p_list, y_p_list, pot_p_list, Xgrid, Ygrid);
 
fprintf(1,'  Ran surface grid-interpolation okay!\n');  
fprintf(1,'  ************************************\n'); 

% ------ Plot "potential surface" ------------------

%result 

figure
mesh(Xgrid, Ygrid, Zgrid);
attractors_num_X_Y(:, 5) = [];
attractors_num_X_Y(:, 1) = [];
writematrix(attractors_num_X_Y, 'Landscape_2D_Data.txt');
%--Program ends here. All input related changes have to be made below--%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

toc; %From given time to run the program we can modify our programs parameters

%--Functions used in the program--%


%Function which defines the programs parameters
function [x_lower_lim, y_lower_lim, x_upper_lim, y_upper_lim, grid_spacing, num_time_steps, timestep, tolerance] = program_parameters()

%------Inputs: You can change the following------%
grid_spacing = 2;

x_upper_lim = 20.0;   % upper limit of x, y (zero to ...)
y_upper_lim = 20.0;

x_lower_lim = -20.0;
y_lower_lim = -20.0;


%--No. of time steps for integrating along each path (to ensure uniform
%arrays)--%
num_time_steps = 1400;   %Choose high-enough number for convergence with given dt


%--Time step and tolerance to test for convergence%
timestep = 1.0e-2;  % ** Check convergence for assigned "numTimeSteps" and tol
tolerance = 1.0e-4;   % **

%------All Inputs end-----
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [gridlines] = graph_parameters()
    gridlines = 45;
end



%Functions to describe the dynamical system 
function der_first_var = calc_first_variable(x, y)

%------Inputs: You can change the following------%
%--Parameters of the equation of the first variable of the dynamical system--%%
n = 4;
a1 = 10.0;
a2 = 10.0;
Kdxx = 4;
Kdyx = 4;
Kdyy = 4;
Kdxy = 4;
b1 = 10.0;
b2 = 10.0;
k1 = 1.0;
k2 = 1.0;
c1 = 0;
c2 = 0;

%--LHS is to be changed according to system--%
der_first_var = c1 + a1*(x^n)/(Kdxx^n + (x^n)) + (b1*(Kdyx^n))/(Kdyx^n + (y^n)) - (x*k1);  
%der_first_var = x - y^3;
%------All Inputs end-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function der_second_var = calc_second_variable(x, y)

%------Inputs: You can change the following------%
%--Parameters of the equation of the second variable of the dynamical system--%
n = 4;
a1 = 10.0;
a2 = 10.0;
Kdxx = 4;
Kdyx = 4;
Kdyy = 4;
Kdxy = 4;
b1 = 10.0;
b2 = 10.0;
k1 = 1.0;
k2 = 1.0;
c1 = 0;
c2 = 0;

%--LHS is to be changed according to system--%
der_second_var = c1 + a2*(y^n)/(Kdyy^n + (y^n)) + (b2*(Kdxy^n))/(Kdxy^n + (x^n)) - (y*k2);
%der_second_var = y - x^3;
%------All Inputs end-----%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end