%% low-filled composites. (Vf< 13%)
clear
clc

% Import datas
data_file_path = 'xxx';
r = data(:,2)/2; % r (m)
Vf = data(:,1); % Vf
Kf = data(:,4); % Kf (W/m K)
Km = data(:,3); % Km
Ri = data(:,5); % Ri (W m2/K)
K0 = data(:,6);

num_data = size(Vf, 1);
Keff_result = zeros(num_data, 1);

% Define integral variables
    syms x m;
    % Start a parallel pool
if isempty(gcp('nocreate'))
    parpool('local');
end

parfor i = 1:num_data
    
    %% Basic geometry parameters
    l = r(i) * (pi / 3 / Vf(i))^(1/3); % The length of the BCC structural unit side
    a = sqrt(2) * r(i);
    b = sqrt(2/3) * r(i);
    w = 2 * r(i)/sqrt(3);
    
    %% Solve for the thermal resistance of R1 and R3
    % 0<m<r
    ex=(-3*m+sqrt(-3*m^2+4*r(i)^2))/(2*sqrt(2));
    ey=(m+sqrt(-3*m^2+4*r(i)^2))/(2*sqrt(2));
    A1_expr=ey^2+2*int(b/a * sqrt(a^2 - x^2), x, ex, a);
    A1_func = matlabFunction(A1_expr);
    Vf1 = pi * (r(i)^2 - m^2) / 4 ./ A1_func(m);
    Vm1 = 1 - Vf1;
    K1_expr = Vf1 * Kf(i) + Vm1 * Km(i); 
    K1_func=matlabFunction(K1_expr);
    
    R1 =double(integral(@(m) 1 ./ A1_func(m) ./ K1_func(m), 0, r(i)));
    
    R3=R1;

    %% Solve for the thermal resistance of R2
    % r<m<w
    kx=(-3*m-sqrt(-3*m^2+4*r(i)^2))/(2*sqrt(2));
    ky=(m-sqrt(-3*m^2+4*r(i)^2))/(2*sqrt(2));
    A2_expr1=2*int(b/a * sqrt(a^2 - x^2), x, -a, kx)+(ey+ky)*(ex-kx)+2*int(b/a * sqrt(a^2 - x^2), x, ex, a);
    A2_func1 = matlabFunction(A2_expr1);

    % w<m<l-w
    A2_expr2 = pi*a*b;
        
    R2 =double(2*integral(@(m) 1 ./ A2_func1(m) ./ Km(i), r(i), w)+(l-2*w)/(A2_expr2 * Km(i)));

    %% Solve for the thermal resistance of R4
    % 0<m<r
    A4_expr1 = l^2 - A1_expr;
    A4_func1 = matlabFunction(A4_expr1);
    % r<m<w
    A4_expr2 = l^2 - A2_expr1;
    A4_func2 = matlabFunction(A4_expr2);
    % w<m<l-w
    A4_expr3 = l^2 - A2_expr2;
    
    R4 = double(2*integral(@(m) 1 ./ A4_func1(m) ./ Km(i), 0 , r(i))+2*integral(@(m) 1 ./ A4_func2(m) ./Km(i), r(i), w)+(l-2*w)/(A4_expr3 * Km(i)));

    %% Solve for the thermal resistance of RI
    RI = double(Ri(i)/ (pi * r(i)^2 / 2));

    %% Solve for the thermal resistance of Re
    Re = 1 / (1 / (2 * R1 + R2 + 2 * RI) + 1 / R4);

    %% Solve for the thermal resistance of R and Keff
    R = Re / 2;
    curr_Keff = (2 * l) / (2 * l)^2 / R;

    %% Stores the results of the current row in the result matrix
    R1_result(i) = double(R1);
    R2_result(i) = double(R2);
    R3_result(i) = double(R3);
    R4_result(i) = double(R4);
    RI_result(i) = double(RI);
    Keff_result(i) = double(curr_Keff);
end

Keff_result = real(Keff_result);
R1_result = real(R1_result)';
R2_result = real(R2_result)';
R3_result = real(R3_result)';
R4_result = real(R4_result)';
RI_result = real(RI_result)';

% Display the results
disp(R1_result');
disp(R2_result');
disp(RI_result');
disp(Keff_result');

% Calculate the error between each K0 and Keff
errors = abs(K0 - Keff_result) ./ K0 * 100;

% Calculate the average error of the whole
mean_error = mean(errors);

% Displays the error result
disp('Error between each K0 and Keff:');
disp(errors);
disp('Overall average error:');
disp(mean_error);