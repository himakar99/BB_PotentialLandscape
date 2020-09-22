%-----Set-up; Inputs: You can change the following-----
xmin = -25;%Lower bound of interval
xmax = 25; %Upper bound of interval
%-----All Inputs end-----

%-----Initial Value; Inputs: You can change the following-----
U0 = 0; %Initialisation of U at xmin
%-----All Inputs end-----

%Solving ODE using ODE Solver
[xSol, USol] = ode45(@(x,U) Func(x,U), [xmin xmax], U0); 

%Generating Data into .txt file
data_landscape = [xSol, USol]; 
plot(xSol, USol)
writematrix(data_landscape, 'Landscape_Data.txt');

%In order to read data
readmatrix('Landscape_Data.txt');