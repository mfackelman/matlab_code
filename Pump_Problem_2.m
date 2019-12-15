%%% Solving multipipe systems %%%

clear;
clc;

syms Va Vb Vc

%Vdot = 325; %[gpm]
rho = 60.121; %[lbm/ft3]
mu = 2.04e-4; %[lbm/(ft*s)]
d = 0.3333333333; %[ft]
ee = .00015; %[m] roughness
La = 10; %[ft] top branch
Lb = 10; %[ft] bottom branch
Lc = 100; %[ft] entrance and exit pipe
g = 32.2; %[ft/s2]
Kent = 0.5; %K at entrance
Kex = 1.0; %K at exit
Kelbow = 0.3; %K for elbow
Ktee = 0.6; %K for tee
Kheat = 50; %K for heat exchanger

%% Part 1
Kglobe = linspace(.01,1001,100); %guess K value for the globe valve (fully open, K=6); FIX to givens
A = 0.087266; %area of pipe

for i = 1:length(Kglobe)
Rea = rho*Va*d/mu;
Reb = rho*Vb*d/mu;
Rec = rho*Vc*d/mu;

% Using Haaland
fa = (1/( -1.8*log10(6.9/Rea + ((ee/d)/3.7)^1.11)))^2;
fb = (1/( -1.8*log10(6.9/Reb + ((ee/d)/3.7)^1.11)))^2;
fc = (1/( -1.8*log10(6.9/Rec + ((ee/d)/3.7)^1.11)))^2;

Vdot = Vc*A;

hp = (-6e-5)*Vdot^2 + 0.0138*Vdot + 31.603; 

eq1 = Vc == Va + Vb;
eq2 = Vb^2 == (Va^2*((fa*La)/d + Kglobe(i) + 2*Kelbow))/((fb*Lb)/d + Kheat + 2*Kelbow);
eq3 = -hp + 2*(Vc^2/(2*g))*((fc*Lc)/d + Kent + 2*Ktee + Kex) + (Va^2/(2*g))*((fa*La)/d + Kglobe(i) + 2*Kelbow) == 10;

S = vpasolve([eq1, eq2, eq3], [Va, Vb, Vc], [5,5,10]);

Vdot_num(i) = S.Vc*A*448.831;
Hp(i) = (-6e-5)*Vdot_num(i)^2 + 0.0138*Vdot_num(i) + 31.603;

plot(Vdot_num,Hp)
end


fprintf('Va = %.3f ft/s \n', S.Va);
fprintf('Vb = %.3f ft/s \n', S.Vb);
fprintf('Vc = %.3f ft/s \n', S.Vc);
fprintf('Hp = %.3f ft \n', Hp);



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
fprintf('Hp = %.3f ft \n', S.Hp);


%% Scratch
clc
Vc = linspace(0,600,10);
disp(Vc);
