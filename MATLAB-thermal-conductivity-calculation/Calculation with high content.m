%% high-filled composites. (Vf> 13%)
clear
clc

% Import datas
data_file_path = 'xxxx';
r = data(:,2)/2; 
Vf = data(:,1); 
Kf = data(:,4); 
Km = data(:,3); 
Ri = data(:,5); 
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
    l = r(i) * (pi / 3 / Vf(i))^(1/3); 
    a = sqrt(2) * r(i);
    b = sqrt(2/3) * r(i);
    w = 2 * r(i)/sqrt(3);
    
    %% Solve for the thermal resistance of R1 and R3
    % 0<m<l-w
    ex=(-3*m+sqrt(-3*m^2+4*r(i)^2))/(2*sqrt(2));
    ey=(m+sqrt(-3*m^2+4*r(i)^2))/(2*sqrt(2));
    A1_expr1=ey^2+2*int(b/a * sqrt(a^2 - x^2), x, ex, a);
    A1_func1 = matlabFunction(A1_expr1);
    Vf1_1 = pi * (r(i)^2 - m^2) / 4 ./ A1_func1(m);
    Vm1_1 = 1 - Vf1_1;
    K1_expr1 = Vf1_1 * Kf(i) + Vm1_1 * Km(i); 
    K1_func1=matlabFunction(K1_expr1);
    % l-w<m<l-r
    mx=(3*l-3*m-sqrt(-3*(l-m)^2+4*r(i)^2))/(2*sqrt(2));
    my=(l-m+sqrt(-3*(l-m)^2+4*r(i)^2))/(2*sqrt(2));
    ox=(3*l-3*m+sqrt(-3*(l-m)^2+4*r(i)^2))/(2*sqrt(2));
    oy=(l-m-sqrt(-3*(l-m)^2+4*r(i)^2))/(2*sqrt(2));
    A1_expr2=ey^2+2*int(b/a * sqrt(a^2 - x^2), x, ex, mx)+(oy+my)*(ox-mx)+2*int(b/a * sqrt(a^2 - x^2), x, ox, a);
    A1_func2 = matlabFunction(A1_expr2);
    Vf1_2 = pi * (r(i)^2 - m^2) / 4 ./ A1_func2(m);
    Vm1_2 = 1 - Vf1_2;
    K1_expr2 = Vf1_2 * Kf(i) + Vm1_2 * Km(i); 
    K1_func2=matlabFunction(K1_expr2);

    R1 =double(integral(@(m) 1 ./ A1_func1(m) ./ K1_func1(m), 0, l-w)+integral(@(m) 1 ./ A1_func2(m) ./ K1_func2(m), l-w, l-r(i)));
    
    R3=R1;

    %% Solve for the thermal resistance of R12
    % l-r<m<r
    A2_expr=ey^2+2*int(b/a * sqrt(a^2 - x^2), x, ex, mx)+my^2;
    A2_func = matlabFunction(A2_expr);

    Vf2=pi * (r(i)^2 - m^2 + r(i)^2 - (l-m).^2) / 4 ./ A2_func(m);
    Vm2 = 1 - Vf2;
    K2_expr = Vf2 * Kf(i) + Vm2 * Km(i); 
    K2_func=matlabFunction(K2_expr);

    R2 =double(integral(@(m) 1 ./ A2_func(m) ./ K2_func(m), l-r(i), r(i)));

    %% Solve for the thermal resistance of R4
    % 0<m<l-w
    A4_expr1 = l^2 - A1_expr1;
    A4_func1 = matlabFunction(A4_expr1);
    % l-w<m<l-r
    A4_expr2 = l^2 - A1_expr2;
    A4_func2 = matlabFunction(A4_expr2);
    % l-r<m<r
    A4_expr3 = l^2 - A2_expr;
    A4_func3 = matlabFunction(A4_expr3);
    
    R4 = double(2*integral(@(m) 1 ./ A4_func1(m) ./ Km(i), 0 , l-w)+2*integral(@(m) 1 ./ A4_func2(m) ./Km(i), l-w, l-r(i)))+integral(@(m) 1 ./ A4_func3(m) ./Km(i), l-r(i), r(i));

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

