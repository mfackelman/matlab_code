% Clear everything from the symbol table
clear;

% Clear the command window to make it more readable
clc;

% Create and display descriptive information
Assignment = "Solving Multipipe Systems";
Author     = "Molly Fackelman";
Class      = "Fluid Dynamics";
Date       = "11/25/2019";

fprintf( ...
  '%-15s%s\n%-15s%s\n%-15s%s\n%-15s%s\n\n', ...
  "Assignment:", Assignment, "Author:", Author, ...
  "Class:", Class, "Date:", Date ...
);

%%%%%%%%%%%%%
%%  Setup  %%
%%%%%%%%%%%%%

%%%
% Default values that may be changed by the user
%%%


% Use this to automatically accept all default values, which will speed
% up testing when your model becomes more reliable.
Accept_All_Defaults = Get_Bool_YesOrNo("Accept all defaults", 0);

% This is the default verbosity level, which controls how much output
% is displayed on the screen.
%   0 = Minimum
%   1 = A little more...
%   2 = A little more...
% And so on. You control this in the script by deciding which statements
% will be displayed based on the verbosity level in an 'if' statement.
Verbosity_Default = 1;

% Other default values
KGlobe_Count_Default = 100;
KGlobe_Start_Default = 0.01;
KGlobe_End_Default   = 1001.0;
Pipe_Area_Default    = 0.087266;

%%%
% Default values that may _not_ be changed by the user
%%%
%Vdot = 325; % Units: gpm
rho = 60.121;       % Units: lbm/ft3
mu  = 2.04e-4;      % Units: lbm/(ft*s)
d   = 0.3333333333; % Units: ft
ee  = .00015;       % [m] roughness
La  = 10;           % top branch; Units: ft
Lb  = 10;           % bottom branch; Units: ft
Lc  = 100;          % entrance and exit pipe; Units: ft
g   = 32.2;         % Units: ft/s2

% Create K value variables
Kent   = 0.5; % K at entrance
Kex    = 1.0; % K at exit
Kelbow = 0.3; % K for elbow
Ktee   = 0.6; % K for tee
Kheat  = 50;  % K for heat exchanger

% Explicitly declare Va, Vb, and Vc to be symbols
syms Va Vb Vc

if Accept_All_Defaults == 1
    Verbosity = Verbosity_Default;
else
    Verbosity = Get_Positive_Integer( ...
        "Enter a positive integer to set the verbosity level", ...
        Verbosity_Default ...
    );
end

%%%%%%%%%%%%%%
%%  Part 1  %%
%%%%%%%%%%%%%%

fprintf("\n");
fprintf(repmat('-', 1, 12));
fprintf("\n-- PART 1 --\n");
fprintf(repmat('-', 1, 12));
fprintf("\n\n");

% Get input values
if Accept_All_Defaults == 1
    KGlobe_Count = KGlobe_Count_Default;
    KGlobe_Start = KGlobe_Start_Default;
    KGlobe_End   = KGlobe_End_Default;
    Pipe_Area    = Pipe_Area_Default;
else
    [KGlobe_Count, KGlobe_Start, KGlobe_End] = Get_KGlobe_Parameters( ...
        KGlobe_Count_Default, KGlobe_Start_Default, KGlobe_End_Default ...
    );
    Pipe_Area = Get_Float("Pipe Area", Pipe_Area_Default);
end

% Create the list of K values for the globe valve
if Verbosity > 0
    fprintf("\nCreating %d KGlobe values from %0.2f to %0.2f\n", ...
        KGlobe_Count, KGlobe_Start, KGlobe_End);
end
Kglobe = linspace(KGlobe_Start, KGlobe_End, KGlobe_Count);

% Solve and populate a vector of Hp values
if Verbosity > 0
    fprintf("\nLooping through KGlobe values");
    if Verbosity > 1
        fprintf(": ");
    else
        fprintf("\n");
    end
end

% Preallocate Vdot_num and Hp to reduce overhead in the script. This
% will make complex algorithms run much faster especially if the vectors
% are large.
Vdot_num = zeros(1, KGlobe_Count);
Hp = zeros(1, KGlobe_Count);

for i = 1 : KGlobe_Count 
    if Verbosity > 1
        Status(i, KGlobe_Count);
    end
    
    % Calculate Rea, Reb, and Rec
    % <Explain what everything means here>
    Rea = rho * Va * d / mu;
    Reb = rho * Vb * d / mu;
    Rec = rho * Vc * d / mu;

    % Calculate fa, fb, and rc
    % <Explain what everything means here>
    % Using Haaland
    fa = (1/( -1.8*log10(6.9/Rea + ((ee/d)/3.7)^1.11)))^2;
    fb = (1/( -1.8*log10(6.9/Reb + ((ee/d)/3.7)^1.11)))^2;
    fc = (1/( -1.8*log10(6.9/Rec + ((ee/d)/3.7)^1.11)))^2;

    % Calculate Vdot
    % <Explain what it is here>
    Vdot = Vc * Pipe_Area;

    % Calculate hp
    % <Explain what it is here>
    hp = (-6e-5) * Vdot^2 + 0.0138 * Vdot + 31.603; 

    % Calculate EQ1
    % <Explain what it is here>
    eq1 = Vc == Va + Vb;
    
    % Calculate EQ2
    % <Explain what it is here>
    eq2 = ...
        Vb^2 == ...
        (Va^2 * ((fa * La) / d + Kglobe(i) + 2*Kelbow)) ...
        / ((fb * Lb)/d + Kheat + 2*Kelbow);
    
    % Calculate EQ3
    % <Explain what it is here>
    eq3 = ...
        -hp ... 
        + 2 * (Vc^2/(2*g)) * ((fc*Lc)/d + Kent + 2 * Ktee + Kex) ...
        + (Va^2/(2*g)) * ((fa*La)/d ...
        + Kglobe(i)...
        + 2*Kelbow) ...
        == 10;

    % Calculate S
    % <Explain what it is here>
    S = vpasolve([eq1, eq2, eq3], [Va, Vb, Vc], [5,5,10]);

    % Calculate vectors Vdot_num and Hp
    % <Explain what they are here>
    Vdot_num(i) = S.Vc * Pipe_Area * 448.831;
    Hp(i) = ...
        (-6e-5) * Vdot_num(i)^2 ...
        + 0.0138 * Vdot_num(i) ...
        + 31.603;
end

% Print a newline to make sure the output is formatted properly
fprintf("\n");

if Verbosity > 1
    fprintf('Va = %.3f ft/s \n', S.Va);
    fprintf('Vb = %.3f ft/s \n', S.Vb);
    fprintf('Vc = %.3f ft/s \n', S.Vc);
    fprintf('Hp = %.3f ft \n', Hp);
end



%%%%%%%%
% I commented this out because it doesn't seem to work quite right and
% you may want to try to solve the problem again using the new functions
% I created. Using the methods I showed above while programming will
% make your life a lot easier.

%{
% copy and paste equations, run new loop with linspace for the velocities
% at a single K

%% Part 2
clc
Vc = linspace(2,600,10);
Kglobe = 1; %set K for the globe valve
A = 0.087266; %area of pipe

for i = 1:length(Vc)
Rea = rho*Va*d/mu;
Reb = rho*Vb*d/mu;
Rec = rho*Vc(i)*d/mu;

% Using Haaland
fa = (1/( -1.8*log10(6.9/Rea + ((ee/d)/3.7)^1.11)))^2;
fb = (1/( -1.8*log10(6.9/Reb + ((ee/d)/3.7)^1.11)))^2;
fc = (1/( -1.8*log10(6.9/Rec(i) + ((ee/d)/3.7)^1.11)))^2;

eq1 = Vc(i) == Va + Vb;
eq2 = Vb^2 == (Va^2*((fa*La)/d + Kglobe + 2*Kelbow))/((fb*Lb)/d + Kheat + 2*Kelbow);
eq3 = -hp + 2*(Vc(i).^2/(2*g))*((fc(i)*Lc)/d + Kent + 2*Ktee + Kex) + (Va^2/(2*g))*((fa*La)/d + Kglobe + 2*Kelbow) == 10;

S = vpasolve([eq1, eq2, eq3], [Va, Vb, hp], [5,5,10]);

Hp(i) = S.hp;

end

disp(Vc);
fprintf('Hp = %.3f ft \n', Hp);


%% Scratch
clc
Vc = linspace(0,600,10);
disp(Vc);
%}

%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
%% Helper Functions %%
%%%%%%%%%%%%%%%%%%%%%%

% Get the K parameters from the user for the globe valve.
function [Count, Start, End] = Get_KGlobe_Parameters( ...
    Count_Default, Start_Default, End_Default)
    
    % Get the count from the user
    Count = Get_Positive_Integer("Number of KGlobe Values", Count_Default);
    
    % Get the first KGlobe value from the user
    Start = Get_Float("First KGlobe value", Start_Default);
    
    % Get the first KGlobe value from the user
    End = Get_Float("Last KGlobe value", End_Default);
end

% This function asks the user to say yes or no to a prompt and then
% returns a boolean value (0 for false and 1 for true). It will ask
% for more input until the condition is satisfied or until the user
% enters nothing, which is an indication that the default value should
% be used.
function V = Get_Bool_YesOrNo(Prompt, Default)
    Default_String = "N";
    if Default == 1
        Default_String = "Y";
    end
    
    % Print a message showing the prompt with the default value.
    Prompt = sprintf( ...
        "%s [%s]: ", ...
        Prompt, ...
    	Default_String ...
    );

    % Loop forever until the user enters a float or the script
    % is terminated with Ctrl+C or some other method.
    while 1
        V = input(Prompt,'s');

        % If the user entered nothing, use the default.
        if isempty(V)
            V = Default;
            break;
        end

        % Otherwise, the user entered something and it needs to be parsed.
        % If the user entered Y or y, return true.
        if V == "Y" || V == "y"
            V = 1;
            break;
        
        % If the user entered N or n, return false.
        elseif V == "N" || V == "n"
            V = 0;
            break;
        
        % Otherwise, the user entered an invalid choice.
        else
            disp('Please enter Y (yes) or N (no)')
        end
        
    end
end

% This function asks the user for a float. If the input value is
% not a float, it will ask for more input until the condition
% is satisfed. The only exception is if the user enters nothing, which
% is an indication that the default value should be used.
function V = Get_Float(Prompt, Default)
    % Print a message showing the prompt with the default value.
    Prompt = sprintf( ...
        "%s [%d]: ", ...
        Prompt, ...
        Default ...
    );
    
    % Loop forever until the user enters a float or the script
    % is terminated with Ctrl+C or some other method.
    while 1
        V = input(Prompt,'s');

        % If the user entered nothing, use the default.
        if isempty(V)
            V = Default;
            break;
        end
        
        % Convert the input string to a number.
        V = str2double(V);
        
        % Otherwise, the user entered something and it needs to be parsed.
        % If the user entered a non-float, ask for more input.
        if ~isfloat(V)
            disp('Please enter a valid floating point number')

        % Otherwise, the user entered a valid float.
        else
            break;
        end
    end
end

% This function asks the user for a positive integer. If the input value is
% not a positive integer, it will ask for more input until the condition
% is satisfed. The only exception is if the user enters nothing, which
% is an indication that the default value should be used.
function V = Get_Positive_Integer(Prompt, Default)
    % Print a message showing the prompt with the default value.
    Prompt = sprintf( ...
        "%s [%d]: ", ...
        Prompt, ...
        Default ...
    );
    
    % Loop forever until the user enters an integer or the script
    % is terminated with Ctrl+C or some other method.
    while 1
        V = input(Prompt,'s');

        % If the user entered nothing, use the default.
        if isempty(V)
            V = Default;
            break;
        end
        
        % Convert the input string to a number.
        V = str2double(V);
        
        % Otherwise, the user entered something and it needs to be parsed.
        % If the user entered a non-integer, ask for more input.
        if ~isIntegerValue(V) || V < 1
            disp('Please enter a positive integer')

        % Otherwise, the user entered a valid integer.
        else
            break;
        end
    end
end

% This function checks to see if X is an integer. Note that it also works
% with infinity and -infinity whereas many other integer checking
% solutions do not.
function T = isIntegerValue(X)
    T = (mod(X, 1) == 0);
end

% This function will output a status message that will overwrite itself
% to prevent spamming the console window with messages.
function S = Status(Iteration, Total)
    persistent Remove_Count;
    
    if isempty(Remove_Count)
        Remove_Count = 0;
    end
    
    % Clear the last status output value
    fprintf(repmat('\b', 1, Remove_Count));
    
    % Create the new status output value
    Status_String = sprintf( ...
        '%d out of %d (%0.f percent)', ...
        Iteration, Total, Iteration / Total * 100 ...
    ); 
    
    % Store the number of characters to remove next time
    Remove_Count = numel(Status_String);
    
    % Output the new status value
    fprintf(Status_String);
    
    S = 0;
end
