% BEAM_ANALYSIS                  Rikky Duivenvoorden -- 07/27/2015
%   Calculates the shear force, bending moment, slope and deflection of
%   a slender beam under uniform loading with one or more intermediate
%   supports and no cantilevers, using Euler-Bernoulli Beam theory.
%   The supports are assumed to be level and the beam is constrained to
%   zero deflection at each support. Additionally, the moment is zero at
%   the beam ends (x = 0 and L), and the beam is assumed to have constant
%   bending stiffness EI along the span.
%   elem_length is a vector containing the lengths of individual spans
%   between supports (i.e. sum(elem_length) = total_span), EI is the
%   bending stiffness, and q is the uniform load applied to the beam.
%   Units in inches, lbf.in2, and lbf/in respectively results in shear
%   force in lbf, moment in lbf.in, and deflection in inches.

function [x, w, w_p, M, V, deflections] = beam_analysis(elem_length, EI, q, plot_flag)

% % elem_length = [141.75 276.25];
% elem_length = [141.75 134.625 141.625];
% E = 200000*(25.4^2/9.8/0.454); % Young's Modulus of Elasticity - psi
% I = 48; % Second Moment of Area - in4
% EI = E*I; %lbf.in2
% % EI = 3.28e13/(9.8*0.454*25.4^2); %lbf.in2 for W12x175
% 
% floor_load = 40; %psf
% floor_semispan = 17.5; %ft
% q = floor_load*floor_semispan/12; % uniform load - lbf/in

num_elements = length(elem_length);

total_span = sum(elem_length);

% create matrix of middle elements (see pdf for details)
A = zeros(num_elements*3-1);

for i = 0:num_elements-1
    A(i*3+1,i*3+1) = -elem_length(i+1)^2/6;
    A(i*3+1,i*3+2) = -elem_length(i+1)/2;
    A(i*3+1,i*3+3) = -1;
    A(i*3+2,i*3+1) = -elem_length(i+1)^2/2;
    A(i*3+2,i*3+2) = -elem_length(i+1);
    A(i*3+2,i*3+3) = -1;
    A(i*3+2,i*3+6) = 1;
    A(i*3+3,i*3+1) =  -elem_length(i+1);
    A(i*3+3,i*3+2) = -1;
    A(i*3+3,i*3+5) = 1;
end

% cutout second column and 2nd last row for end elements (zero moment)
A = [A(:,1) A(:,3:end-3)];
A = [A(1:end-2,:); A(end,:)];

% create load vector
B = zeros(num_elements*3-1,1);

for i = 0:num_elements-1
    B(i*3+1) = q/24*elem_length(i+1)^3;
    B(i*3+2) = q/6*elem_length(i+1)^3;
    B(i*3+3) = q/2*elem_length(i+1)^2;
end

B = [B(1:end-2); B(end)];

% solve for integration constants
C = A\B;

% Find shear force, moment, slope and deflection at global spanwise x
% locations
x = linspace(0,total_span,1000);

% initialize deflection, slope, moment and shear force vectors
w = zeros(size(x));
w_p = w;
M = w;
V = w;

% plot shear, moment and deflection

allowable_w = zeros(size(x));

%loop through spanwise locations
for i=1:length(x)
    % keep a sum_len variable that tracks where x(i) should search along
    % the span of the beam
    sum_len = 0;
    
    % search which element (hence which constants of integration) to use
    for j=1:num_elements
        
        sum_len = elem_length(j) + sum_len;
        
        % exit the loop if x(i) is no longer in the jth beam element
        if x(i) < sum_len
            allowable_w(i) = elem_length(j)/360; % target deflection of L/720
            % if j=1, then x(i) is in the first span (c_12 = 0)
            if j==1 
                w(i) = -(x(i)^4*(q/24) + x(i)^3*C(1)/6 + x(i)*C(2))*(1/EI);
                w_p(i) = -(q/6*x(i)^3 + C(1)/2*x(i)^2+C(2))*(1/EI);
                M(i) = -q/2*x(i)^2 - C(1)*x(i);
                V(i) = -q*x(i)-C(1);
                % break the loop for x(i+1)
                break;
            % otherwise, x(i) is in the rest of the beam where c_i2 ~= 0
            else
                local_x = x(i) + elem_length(j) - sum_len;
                w(i) = -(local_x^4*(q/24) + local_x^3*(C((j-1)*3)/6) + local_x^2*C((j-1)*3+1)/2 + local_x*C((j-1)*3+2))*(1/EI);
                w_p(i) = -(q/6*local_x^3 + C((j-1)*3)/2*local_x^2 + C((j-1)*3+1)*local_x + C((j-1)*3+2))*(1/EI);
                M(i) = -q/2*local_x^2 - C((j-1)*3)*local_x - C((j-1)*3+1);
                V(i) = -q*local_x-C((j-1)*3);
                % break the loop for x(i+1)
                break;
            end
        end
        
        
    end
end

% Determine max deflection
deflections = zeros(num_elements,1);

span_location = 0;
prev_location = 1;
sum_len = 0;
for i = 1:num_elements
    local_span = elem_length(i);
    sum_len = elem_length(i) + sum_len;
    
    j = 1;
    while span_location < sum_len
        j = j+1;
        span_location = x(j);
    end
    
    deflections(i) = max(abs(w(prev_location:j)))/local_span*360;
    prev_location = j;
end

% Plot diagrams
if plot_flag == 1
    figure(1)
    plot(x,w)
    hold on
    plot(x, allowable_w, 'r')
    plot(x, -allowable_w, 'r')
    grid on
    % axis equal
    xlabel('')
    ylabel('')
    title('')

%     figure(2)
%     plot(x,w_p)
%     grid on
%     xlabel('')
%     ylabel('')
%     title('')
% 
%     figure(3)
%     plot(x,M)
%     grid on
%     xlabel('')
%     ylabel('')
%     title('')
% 
%     figure(4)
%     plot(x,V)
%     grid on
%     xlabel('')
%     ylabel('')
%     title('')
end

end
    
