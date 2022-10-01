clc; clear all;

fileid = fopen('example.txt');

%// Node data
NON = fscanf(fileid, 'Node = %d', 1);
x = fscanf(fileid,'%f %f %f',[5,NON]);
TDOF = 2*NON;
node = x';

%// Element data
NOE = fscanf(fileid, '\nElement = %d' , 1);
y = fscanf(fileid,'%f %f',[7,NOE]);
element = y';
con_table = element(:,1:3);
elem = zeros(NOE,11);
elem(:,1)=1:NOE;
for i = 1:NOE
    a = con_table(i,2);
    b = con_table(i,3);
    A = find(node(:,1)==a);
    B = find(node(:,1)==b);
    elem(i,2:6) = node(A,:);
    elem(i,7:11) = node(B,:);
end

%// Material and Geometrical Properties
E = y(4,:);
Area = y(5,:);
Coeffthrmlexp = y(6,:);
delta_T = y(7,:);

%// Load Data
F_k = fscanf(fileid, '\n\nKnown Forces = %d' , 1);
load = fscanf(fileid, '%f %f' , [2,F_k]);
load = load';
F = zeros(TDOF,2);
F(:,1) = 1:TDOF;
for i = 1:F_k
    a = load(i,1);
    A = find(load(:,1)==a);
    F(a,2) = load(A,2);
end
load = F(:,2);


%// Boundary Conditions
BC = fscanf(fileid, '\n\n\nBoundary Condition = %d' , 1);
boun_cond = fscanf(fileid, '%f %f' , [2,BC]);
boun_cond = boun_cond';
Q = sym('Q',[TDOF 1]);

for i = 1:TDOF
    if boun_cond(i,2) == 0
        Q(i) = 0;
    end
end
%// Multi-point Constraint
incl_roller = fscanf(fileid, '\n\n\n\nInclined Roller = %d' , 1);
incl_roller_data = fscanf(fileid, '%f %f' , [2,incl_roller]);
incl_roller_node = incl_roller_data(1,:);
incl_roller_theta = incl_roller_data(2,:);
beta_0 = zeros(TDOF,2) ; 
beta_1 = zeros(TDOF,2) ; 
beta_2 = zeros(TDOF,2) ; 
for i = incl_roller_node
    
        beta_1(2*i-1,:) = [-sin(incl_roller_theta * pi/180),1] ; 
        beta_2(2*i,:) = [cos(incl_roller_theta * pi/180),1] ;
   
end

%// Calculate length and the factors l and m for each element:
L = zeros(NOE,1);
l = zeros(NOE,1);
m = zeros(NOE,1);
for i = 1:NOE
    L(i) = sqrt( ( elem(i,9) - elem(i,4) )^2 + ( elem(i,11) - elem(i,6) )^2 ) ; % length, mm
    l(i) = ( elem(i,9) - elem(i,4) ) / L(i) ; % l factor (unitless)
    m(i) = ( elem(i,11) - elem(i,6) ) / L(i) ; % m factor (unitless)
end

%// Temperature load
temp_load = zeros(TDOF,1);
for i = 1:NOE
    T_load = E(i).*Area(i).*Coeffthrmlexp(i).*delta_T(i)*[-l(i);-m(i);l(i);m(i)];
    index = [ elem(i,3);elem(i,5);elem(i,8);elem(i,10) ];
    temp_load(index) = temp_load(index) + T_load;
    
end

%// Construct global stiffness matrix from element local matrices:
K = zeros(TDOF,TDOF);
for i = 1:NOE
    k = E(i) * Area(i) / L(i) .* [l(i)^2 l(i)*m(i) -l(i)^2 -l(i)*m(i);l(i)*m(i) m(i)^2 -l(i)*m(i) -m(i)^2;-l(i)^2 -l(i)*m(i) l(i)^2 l(i)*m(i);-l(i)*m(i) -m(i)^2 l(i)*m(i) m(i)^2 ] ; % element local stiffness matrix (kN/mm)
    index = [ elem(i,3) elem(i,5) elem(i,8) elem(i,10) ] ; 
    K(index,index) = K(index,index) + k ;
    disp(K)
end

C_factor = 10^8;
C = abs(max(K,[],'all')) * C_factor;
q = zeros(TDOF,1); 
modify = false(TDOF,1);

K_p = K; 
F_p = load +temp_load; 

for i = 1:TDOF
    if isreal(Q(i)) % condition for which a displacement is known by the user
        q(i) = Q(i) ; % if condition is met, overwrite q_i with Q_i
        modify(i) = true ; % set the i-th logical entry to true (1) so the program knows to modify
    end
    if modify(i) == true
        K_p(i,i) = K_p(i,i) + C ; % modify i-th diagonal entry K_ii wherever Q_i is known
        F_p(i) = F_p(i) + C * q(i) ; % modify i-th entry F_i wherever Q_i is known
    end
    for j = 1:TDOF
        if beta_1(i,2) == beta_2(j,2) % condition for which i and j are a set of multi-point constraints
            K_p(i,j) = K_p(i,j) + C * beta_1(i,1) * beta_2(j,1) ; % modify K_ij
            K_p(j,i) = K_p(j,i) + C * beta_1(i,1) * beta_2(j,1) ; % modify K_ji
            K_p(i,i) = K_p(i,i) + C * beta_1(i,1)^2 ; % modify K_ii         
            K_p(j,j) = K_p(j,j) + C * beta_2(j,1)^2 ; % modify K_jj
            F_p(i) = F_p(i) + C * beta_0(i) * beta_1(i) ; % modify F_i
            F_p(j) = F_p(j) + C * beta_0(j) * beta_2(j) ; % modify F_j
        end
    end
end

Q = inv(K_p)*F_p;
Q = double(Q);  % convert symbolic array to double precision

%// Solve for unknown reactions
R = zeros(TDOF,1) ; % initialize vector of reaction forces as array of zeros
for i = 1:TDOF
    if modify(i) == true 
      R(i) = R(i) - C * ( Q(i) - q(i) ) ; % modify i-th entry R_i (kN) wherever Q_i was initially known
    end
    for j = 1:TDOF
        if beta_1(i,2) == beta_2(j,2) % modify R_i and R_j for each set of multi-point constraints
            R(i) = R(i) -C * beta_1(i,1) * ( beta_1(i,1) * Q(i) + beta_2(j,1) * Q(j) - beta_0(i) ) ; % (kN)
            R(j) = R(j) -C * beta_2(j,1) * ( beta_1(i,1) * Q(i) + beta_2(j,1) * Q(j) - beta_0(j) ) ; % (kN)
        end
    end
end

%// Determine elemental strains and stresses:
epsilon = zeros(NOE,1);
sigma = zeros(NOE,1);
for i = 1:NOE
    epsilon(i) = 1 / L(i) * [-l(i) -m(i) l(i) m(i)] * [Q(elem(i,3));Q(elem(i,5));Q(elem(i,8));Q(elem(i,10))] ; % element strain (unitless)
    sigma(i) = E(i) * epsilon(i) - E(i)*Coeffthrmlexp(i)*delta_T(i); % element stress (kN/mm^2)
end
disp('Displacements:')
disp(Q)
disp('Element Strain:')
disp(epsilon)
disp('Element Stress:')
disp(sigma)
disp('Reactions:')
disp(R)

