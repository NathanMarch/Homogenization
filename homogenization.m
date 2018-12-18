function Deff = homogenization(D,varargin)
% homogenization computes the effective diffusivity tensor Deff for a
% two-dimensional heterogenous domain using a semi-analytical method.

% homogenization solves the two-dimensional homogenization problem on a
% finite domain consisting of an m by m grid of blocks, where each of the
% m^2 blocks can have a different diffusivity. The code is applicable to
% domains of size (xn-x0) by (xn-x0) in which each of the m^2 blocks is 
% equally sized.

% For domains in which the domain comprises an m by n grid of blocks (where
% m may not equal n) and/or the blocks in the domain are not all equally
% sized, please use the code homogenisation_CD

% homogenization is an implementation of the semi-analytical method 
% proposed by March, Carr and Turner.

%% If you use this code please city the following publication:

% Nathan G. March, Elliot J. Carr , and Ian W. Turner,
% A fast semi-analytical homogenization method for block heterogeneous
% media, Submitted to Journal of Computational Physics, 
% https://arxiv.org/abs/1812.06680


% -------------------------------------------------------------------------
% Details of inputs and outputs
% -------------------------------------------------------------------------

% INPUTS:

% D = matrix of diffusivities, where D(i,j) is the diffusivity in the
% (i,j)th block
% options = options containing the following parameters:
% x0 = the leftmost point of the domain
% xn = the rightmost point of the domain
% Nx = the number of abscissa to use in all integrals
% N = the number of terms to use in the summation

% OUTPUTS:
% Deff = effective diffusivity tensor

% -------------------------------------------------------------------------
% Check inputs:
% -------------------------------------------------------------------------

if nargin<1
    error('Not enough input arguments.')
elseif nargin == 1
    options = struct;
elseif nargin == 2
    options = varargin{1};
else
    error('Too many input arguments.')
end

% Check that D is square
[m,n] = size(D);
if m ~= n
    error('Diffusivity matrix D must be square. Please use HOMOGENISATION_CD for domains in which D is not square.')
end

% Check that all entries of D are positive
if sum(sum(D>0)) ~= m^2
    error('All entries of D must be positive.')
end

% Check that options is a structure
if ~isa(options,'struct')
    error('options must be a structure.')
end

if isfield(options,'x0')
    x0 = options.x0;
else
    x0 = 0; % default value
end

if isfield(options,'xn')
    xn = options.xn;
    if x0>=xn
        error('x0 must be less than xn.')
    end
else
    xn = 1; % default value
end

if isfield(options,'Nx')
    Nx = options.Nx;
    if round(Nx) ~= Nx && Nx < 1
        error('Nx must be an integer greater than or equal to 1.')
    end
else
    Nx = 16; % default value
end

if isfield(options,'N')
    N = options.N;
    if round(N) ~= N || N < 1 || N > 2*Nx-3 || N > 100
        error('N must be an integer greater than or equal to 1 and less than 2Nx-3 and less than 100.')
    elseif N < 2*Nx-3
        warning('Using N < 2Nx-3 may produce inaccurate results')
    end
else
    N = min(2*Nx-3,100); % default value
end



%% Set up
D = D(end:-1:1,:); % Reverse order of diffusivity matrix as y axis is reversed
x = x0+(xn-x0)/m; % Location of interface
xint = linspace(x0+(x-x0)/(2*Nx),x-(x-x0)/(2*Nx),Nx); % Spacing of abscissa
lx = x-x0; % Block width
bx = zeros(m*n*(2*Nx+1),1); % Right hand side vector b used when homogenization problem is solved in x direction
by = zeros(m*n*(2*Nx+1),1); % Right hand side vector b used when homogenization problem is solved in y direction
wx = xint(2)-xint(1); % Weighting of abscissa used in midpoint rule

N_matrix = 14*m^2*Nx^2+8*m^2*Nx-4*Nx+m^2; % Number of entries in linear system

row = zeros(N_matrix,1); % vector of row_start indices
col = zeros(N_matrix,1); % vector of column indices
val = zeros(N_matrix,1); % vector of values

index = 1; % loop through entries in r,c,v


%% Compute terms used in summation
square_vec_1 = wx/(2*lx^2) * (xint - x0).^2;
square_vec_2 = wx/(2*lx^2) * (xint - x(1)).^2;
gamma = zeros(N,1);

cos_prod_sum = zeros(Nx,Nx);
cosh_cos_prod_sum = zeros(Nx,Nx);
cosh_cos_sum_1 = zeros(Nx,Nx);
cosh_cos_sum_2 = zeros(Nx,Nx);
cosh_cos_sum_3 = zeros(Nx,Nx);
cosh_cos_sum_4 = zeros(Nx,Nx);
for k = 1:N
    gamma(k) = k*pi*sinh(k*pi);
    cos_prod_sum = cos_prod_sum + cos((k*pi/lx) * xint)' * cos((k*pi/lx) * xint)/gamma(k);
    cosh_cos_prod_sum = cosh_cos_prod_sum + cosh(k*pi)*cos((k*pi/lx) * xint)' * cos((k*pi/lx) * xint)/gamma(k);
    cosh_cos_sum_1 = cosh_cos_sum_1 + (-1)^k * cosh((k*pi/lx) * xint)' * cos((k*pi/lx) * xint)/gamma(k);
    cosh_cos_sum_2 = cosh_cos_sum_2 + (-1)^k * cosh((k*pi/lx) * (xint-x(1)))' * cos((k*pi/lx) * xint)/gamma(k);
    cosh_cos_sum_3 = cosh_cos_sum_3 + cosh((k*pi/lx) * xint)' * cos((k*pi/lx) * xint)/gamma(k);
    cosh_cos_sum_4 = cosh_cos_sum_4 + cosh((k*pi/lx) * (xint-x(1)))' * cos((k*pi/lx) * xint)/gamma(k);
end
cos_prod_sum = 2*wx*cos_prod_sum;
cosh_cos_prod_sum = 2*wx*cosh_cos_prod_sum;
cosh_cos_sum_1 = 2*wx*cosh_cos_sum_1;
cosh_cos_sum_2 = 2*wx*cosh_cos_sum_2;
cosh_cos_sum_3 = 2*wx*cosh_cos_sum_3;
cosh_cos_sum_4 = 2*wx*cosh_cos_sum_4;





%% Formulate linear system
% Equations are numbered with the first equation being v_{2,1} - v{1,1},
% second equation being v_{2,2} - v_{1,2} etc

% Equation v_{i+1,j}(x^{(p)},y_i) - v_{i,j}(x^{(p)},y_i) = 0
for i = 1:m-1 % Loop through rows
    for j = 1:n % Loop through columns
        z8 = m*n*(2*Nx) + i*n+j; % indicates position of K_{i+1,j)
        z9 = m*n*(2*Nx) + (i-1)*n+j;  % indicates position of K_{i,j)
        for p = 1:Nx
            e = ((i-1)*n+j-1)*Nx+p; % equation number
            
            row(index) = e;
            col(index) = z8;
            val(index) = 1;
            index = index+1;
            
            row(index) = e;
            col(index) = z9;
            val(index) = -1;
            index = index+1;
            
            for r = 1:Nx
                z1 = ((j-1)*m+i-1)*Nx+r; % indicates position in vector of abscissas
                z2 = ((j-1)*m+i)*Nx+r;
                if i == m-1
                    z3 = ((j-1)*m)*Nx+r;
                else
                    z3 = ((j-1)*m+i+1)*Nx+r;
                end
                
                row(index) = e;
                col(index) = z1;       
                val(index) = val(index) + cos_prod_sum(p,r)/D(i,j);         
                index = index+1;
                
                row(index) = e;
                col(index) = z2;   
                val(index) = val(index)  - (wx/2) * (1/D(i,j) + 1/D(i+1,j)) - cosh_cos_prod_sum(p,r) * (1/D(i,j) + 1/D(i+1,j));      
                index = index+1;
                
                row(index) = e;
                col(index) = z3;  
                val(index) = val(index) + cos_prod_sum(p,r)/D(i+1,j);
                index = index+1;
            end
            
            for r = 1:Nx  % indicates position in vector of abscissas
                z4 = m*n*Nx+((i-1)*n+j-1)*Nx+r;
                z6 = m*n*Nx+(i*n+j-1)*Nx+r;
                if j == n
                    z5 = m*n*Nx+(i-1)*n*Nx+r;
                    z7 = m*n*Nx+(i*n)*Nx+r;
                else
                    z5 = m*n*Nx+((i-1)*n+j)*Nx+r;
                    z7 = m*n*Nx+(i*n+j)*Nx+r;
                end
                
                row(index) = e;
                col(index) = z4;
                val(index) = val(index) + (square_vec_2(p)/D(i,j))  + (cosh_cos_sum_2(p,r)/D(i,j));               
                index = index+1;
                
                row(index) = e;
                col(index) = z5;    
                val(index) = val(index) - square_vec_1(p)/D(i,j) - cosh_cos_sum_1(p,r)/D(i,j);
                index = index+1;
                
                row(index) = e;
                col(index) = z6;
                val(index) = val(index) - square_vec_2(p)/D(i+1,j) - cosh_cos_sum_4(p,r)/D(i+1,j);
                index = index+1;
                
                row(index) = e;
                col(index) = z7;
                val(index) = val(index) + square_vec_1(p)/D(i+1,j) + cosh_cos_sum_3(p,r)/D(i+1,j);
                index = index+1;
            end
        end
    end
end

% Equation v_{i,j+1}(x_j,y_^{(p)}) - v_{i,j}(x_j,y^{(p)}) = 0
for i = 1:m % Loop through rows
    for j = 1:n-1 % Loop through columns
        z8 = m*n*(2*Nx) + (i-1)*n+j+1; % indicates position of K_{i,j+1)
        z9 = m*n*(2*Nx) + (i-1)*n+j;  % indicates position of K_{i,j)
        for p = 1:Nx
            e = (m-1)*n*Nx + ((j-1)*m+i-1)*Nx+p;  % equation number
            row(index) = e;
            col(index) = z8;
            val(index) = 1;
            index = index+1;
            
            row(index) = e;
            col(index) = z9;
            val(index) = -1;
            index = index+1;
            for r = 1:Nx
                z1 = ((j-1)*m+i-1)*Nx+r; % indicates position in vector of abscissas
                z3 = (j*m+i-1)*Nx+r;
                if i == m
                    z2 = (j-1)*m*Nx+r;
                    z4 = j*m*Nx+r;
                else
                    z2 = ((j-1)*m+i)*Nx+r;
                    z4 = (j*m+i)*Nx+r;
                end
                
                row(index) = e;
                col(index) = z1;
                val(index) = val(index) + square_vec_2(p)/D(i,j) + cosh_cos_sum_2(p,r)/D(i,j);
                index = index+1;
                
                row(index) = e;
                col(index) = z2;
                val(index) = val(index) - square_vec_1(p)/D(i,j) - cosh_cos_sum_1(p,r)/D(i,j);
                index = index+1;
                
                row(index) = e;
                col(index) = z3;
                val(index) = val(index) - square_vec_2(p)/D(i,j+1) - cosh_cos_sum_4(p,r)/D(i,j+1);
                index = index+1;
                
                row(index) = e;
                col(index) = z4;
                val(index) = val(index) + square_vec_1(p)/D(i,j+1)+ cosh_cos_sum_3(p,r)/D(i,j+1);
                index = index+1;
                
            end
            
            
            
            for r = 1:Nx  % indicates position in vector of abscissas
                z5 = m*n*Nx+((i-1)*n+j-1)*Nx+r;
                z6 = m*n*Nx+((i-1)*n+j)*Nx+r;
                if j == n-1
                    z7 = m*n*Nx+(i-1)*n*Nx+r;
                else
                    z7 = m*n*Nx+((i-1)*n+j+1)*Nx+r;
                end
                
                row(index) = e;
                col(index) = z5;
                val(index) = val(index) + cos_prod_sum(p,r)/D(i,j);
                index = index+1;
                
                row(index) = e;
                col(index) = z6;
                val(index) = val(index)  - (wx/2) * (1/D(i,j) + 1/D(i,j+1)) - cosh_cos_prod_sum(p,r) * (1/D(i,j) + 1/D(i,j+1));
                index = index+1;
                
                row(index) = e;
                col(index) = z7;
                val(index) = val(index) + cos_prod_sum(p,r)/D(i,j+1);
                index = index+1;      
            end
        end
    end
end

% Equation v_{m,j}(x^{(p)},y_m) - v_{1,j}(x^{(p)},y_0) = y_m - y_0
for j = 1:n % Loop through columns
    z8 = m*n*(2*Nx) + (m-1)*n+j; % indicates position of K_{m,j)
    z9 = m*n*(2*Nx) + j;  % indicates position of K_{1,j)
    for p = 1:Nx
        e = (m-1)*n*Nx+(n-1)*m*Nx+(j-1)*Nx+p; % equation number
        
        row(index) = e;
        col(index) = z8;
        val(index) = 1;
        index = index+1;
        
        row(index) = e;
        col(index) = z9;
        val(index) = -1;
        index = index+1;
        
        by(e) = xn-x0; % right hand side vector
        for r = 1:Nx
            z1 = ((j-1)*m)*Nx+r; % indicates position in vector of abscissas
            z2 = ((j-1)*m+1)*Nx+r;
            z3 = (j*m-1)*Nx+r;
            
            row(index) = e;
            col(index) = z1;
            val(index) = val(index)  + (wx/2) * (1/D(1,j) + 1/D(m,j))  + cosh_cos_prod_sum(p,r) * (1/D(1,j) + 1/D(m,j));
            index = index+1;
            
            row(index) = e;
            col(index) = z2;
            val(index) = val(index) - cos_prod_sum(p,r)/D(1,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z3; 
            val(index) = val(index) - cos_prod_sum(p,r)/D(m,j);        
            index = index+1;  
        end
        
        for r = 1:Nx  % indicates position in vector of abscissas
            z4 = m*n*Nx+(j-1)*Nx+r;
            z6 = m*n*Nx+((m-1)*n+j-1)*Nx+r;
            if j == n
                z5 = m*n*Nx+r;
                z7 = m*n*Nx+(m-1)*n*Nx+r;
            else
                z5 = m*n*Nx+j*Nx+r;
                z7 = m*n*Nx+((m-1)*n+j)*Nx+r;
            end
            row(index) = e;
            col(index) = z4;
            val(index) = val(index) + square_vec_2(p)/D(1,j) + cosh_cos_sum_4(p,r)/D(1,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z5;
            val(index) = val(index) - square_vec_1(p)/D(1,j) - cosh_cos_sum_3(p,r)/D(1,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z6;
            val(index) = val(index) - square_vec_2(p)/D(m,j) - cosh_cos_sum_2(p,r)/D(m,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z7;
            val(index) = val(index) + square_vec_1(p)/D(m,j) + cosh_cos_sum_1(p,r)/D(m,j);
            index = index+1; 
        end
    end
end

% Equation v_{i,n}(x_n,y_^{(p)}) - v_{i,1}(x_0,y^{(p)}) = x_n-x_0
for i = 1:m % Loop through rows
    z8 = m*n*(2*Nx) + i*n; % indicates position of K_{i,n)
    z9 = m*n*(2*Nx) + (i-1)*n+1;  % indicates position of K_{i,1)
    for p = 1:Nx
        e = m*n*Nx+(n-1)*m*Nx+(i-1)*Nx+p; % equation number
        row(index) = e;
        col(index) = z8;
        val(index) = 1;
        index = index+1;
        
        row(index) = e;
        col(index) = z9;
        val(index) = -1;
        index = index+1;
        bx(e) = xn - x0; % right hand side vector
        for r = 1:Nx
            z1 = (i-1)*Nx+r; % indicates position in vector of abscissas
            z3 = ((n-1)*m+i-1)*Nx+r;
            if i == m
                z2 = r;
                z4 = ((n-1)*m)*Nx+r;
            else
                z2 = i*Nx+r;
                z4 =((n-1)*m+i)*Nx+r;
            end
            
            row(index) = e;
            col(index) = z1;
            val(index) = val(index) + square_vec_2(p)/D(i,1) + cosh_cos_sum_4(p,r)/D(i,1);
            index = index+1;
            
            row(index) = e;
            col(index) = z2;
            val(index) = val(index) - square_vec_1(p)/D(i,1) - cosh_cos_sum_3(p,r)/D(i,1);
            index = index+1;
            
            row(index) = e;
            col(index) = z3; 
            val(index) = val(index) - square_vec_2(p)/D(i,n) - cosh_cos_sum_2(p,r)/D(i,n);
            index = index+1;
            
            row(index) = e;
            col(index) = z4;
            val(index) = val(index) + square_vec_1(p)/D(i,n) + cosh_cos_sum_1(p,r)/D(i,n);
            index = index+1;
        end
        
        for r = 1:Nx  % indicates position in vector of abscissas
            z5 = m*n*Nx+((i-1)*n)*Nx+r;
            z6 = m*n*Nx+((i-1)*n+1)*Nx+r;
            z7 = m*n*Nx+(i*n-1)*Nx+r;
            
            row(index) = e;
            col(index) = z5;
            val(index) = val(index)  + (wx/2) * (1/D(i,1) + 1/D(i,n))  + cosh_cos_prod_sum(p,r) * (1/D(i,1) + 1/D(i,n));
            index = index+1;
            
            row(index) = e;
            col(index) = z6;
            val(index) = val(index) - cos_prod_sum(p,r)/D(i,1);
            index = index+1;
            
            row(index) = e;
            col(index) = z7;
            val(index) = val(index) - cos_prod_sum(p,r)/D(i,n);
            index = index+1;
        end
    end
end

%% Formulate rows of linear system based on zero net flux condition
for i = 1:m-1
    for j = 1:n-1
        e = m*n*(2*Nx)+(i-1)*(n-1)+j; % equation number
        for r = 1:Nx
            z1 = ((j-1)*m+i-1)*Nx+r; % indicates position in vector of abscissas
            z2 = ((j-1)*m+i)*Nx+r; % indicates position in vector of abscissas
            
            row(index) = e;
            col(index) = z1;
            val(index) = wx/D(i,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z2;
            val(index) = -wx/D(i,j);
            index = index+1;
        end
        
        for r = 1:Nx
            z3 = m*n*Nx+((i-1)*n+j-1)*Nx+r; % indicates position in vector of abscissas
            z4 = m*n*Nx+((i-1)*n+j)*Nx+r; % indicates position in vector of abscissas
            
            row(index) = e;
            col(index) = z3;
            val(index) = wx/D(i,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z4;
            val(index) = -wx/D(i,j);
            index = index+1;
        end
    end
end

for j = 1:n-1
    e = m*n*(2*Nx)+(m-1)*(n-1)+j; % equation number
    for r = 1:Nx
        z1 = (j*m-1)*Nx+r; % indicates position in vector of abscissas
        z2 = ((j-1)*m)*Nx+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z1;
        val(index) = wx/D(m,j);
        index = index+1;
        
        row(index) = e;
        col(index) = z2;
        val(index) = -wx/D(m,j);
        index = index+1;
    end
    
    for r = 1:Nx
        z3 = m*n*Nx+((m-1)*n+j-1)*Nx+r; % indicates position in vector of abscissas
        z4 = m*n*Nx+((m-1)*n+j)*Nx+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z3;
        val(index)  = wx/D(m,j);      
        index = index+1;
        
        row(index) = e;
        col(index) = z4;
        val(index) = -wx/D(m,j);
        index = index+1;
    end
end

for i = 1:m-1
    e = m*n*(2*Nx)+m*(n-1)+i; % equation number
    for r = 1:Nx
        z1 = ((n-1)*m+i-1)*Nx+r; % indicates position in vector of abscissas
        z2 = ((n-1)*m+i)*Nx+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z1;
        val(index) = wx/D(i,n);
        index = index+1;
        
        row(index) = e;
        col(index) = z2;
        val(index) = -wx/D(i,n);
        index = index+1;
        
    end
    
    for r = 1:Nx
        z3 = m*n*Nx+(i*n-1)*Nx+r; % indicates position in vector of abscissas
        z4 = m*n*Nx+((i-1)*n)*Nx+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z3;
        val(index) = wx/D(i,n);  
        index = index+1;
        
        row(index) = e;
        col(index) = z4;
        val(index) = -wx/D(i,n);
        index = index+1;
    end
end

e = m*n*(2*Nx+1); % equation number

row(N_matrix-m*n+1:N_matrix) = e*ones(m*n,1);
col(N_matrix-m*n+1:N_matrix) = (m*n*(2*Nx)+1:m*n*(2*Nx+1))';
val(N_matrix-m*n+1:N_matrix) = ones(m*n,1);
bx(e) = (xn+x0)/2;
by(e) = (xn+x0)/2;

%% Solve linear system and assign flux values to vectors

A = sparse(row,col,val,m*n*(2*Nx+1),m*n*(2*Nx+1));
qg = A\[bx by];
qgx = qg(:,1);
qgy = qg(:,2);

Deff = zeros(2,2); % effective diffusvity matrix
for direction = 0:1
    
    if direction == 0
        
        q = qgx(1:m*n*Nx);
        g = qgx(m*n*Nx+1:m*n*(2*Nx));
        
    else
        q = qgy(1:m*n*Nx);
        g = qgy(m*n*Nx+1:m*n*(2*Nx));
    end
    q = reshape(q,[Nx,m*n]);
    g = reshape(g,[Nx,m*n]);
    
    %% Build solution
    
    % Compute polynomial coefficients
    a0 = zeros(m,n);
    b0 = zeros(m,n);
    c0 = zeros(m,n);
    d0 = zeros(m,n);
    for i = 1:m-1
        for j = 1:n-1
            if direction == 0
                a0(i,j) = (2*wx/lx) * sum(g(:,(i-1)*n+j))/D(i,j);
                b0(i,j) = (2*wx/lx) * sum(g(:,(i-1)*n+j+1))/D(i,j);
            end
            c0(i,j) = (2*wx/lx) * sum(q(:,(j-1)*m+i))/D(i,j);
            d0(i,j) = (2*wx/lx) * sum(q(:,(j-1)*m+i+1))/D(i,j);
        end
    end
    
    for j = 1:n-1
        if direction == 0
            a0(m,j) = (2*wx/lx)  * sum(g(:,(m-1)*n+j))/D(m,j);
            b0(m,j) = (2*wx/lx)  * sum(g(:,(m-1)*n+j+1))/D(m,j);
        end
        c0(m,j) = (2*wx/lx)*sum(q(:,j*m))/D(m,j);
        d0(m,j) = (2*wx/lx)*sum(q(:,(j-1)*m+1))/D(m,j);
    end
    
    
    for i = 1:m-1
        if direction == 0
            a0(i,n) = (2*wx/lx)*sum(g(:,i*n))/D(i,n);
            b0(i,n) = (2*wx/lx)*sum(g(:,(i-1)*n+1))/D(i,n);
        end
        c0(i,n) = (2*wx/lx)*sum(q(:,(n-1)*m+i))/D(i,n);
        d0(i,n) = (2*wx/lx)*sum(q(:,(n-1)*m+i+1))/D(i,n);
    end
    if direction == 0
        a0(m,n) = (2*wx/lx)*sum(g(:,m*n))/D(m,n);
        b0(m,n) = (2*wx/lx)*sum(g(:,(m-1)*n+1))/D(m,n);
    end
    c0(m,n) = (2*wx/lx)*sum(q(:,m*n))/D(m,n);
    d0(m,n) = (2*wx/lx)*sum(q(:,(n-1)*m+1))/D(m,n);
    
    % Compute trigonometric function coefficients
    
    ak = zeros(m,n,k);
    bk = zeros(m,n,k);
    ck = zeros(m,n,k);
    dk = zeros(m,n,k);
    
    for i = 1:m-1
        for j = 1:n-1
            for k = 1:2:N
                ak(i,j,k) = (2*wx/(lx*D(i,j))) * sum(g(:,(i-1)*n+j) .* cos(k*pi*xint'/lx));
                bk(i,j,k) = (2*wx/(lx*D(i,j))) * sum(g(:,(i-1)*n+j+1) .* cos(k*pi*xint'/lx));
                if direction  == 0
                    ck(i,j,k) = (2*wx/(lx*D(i,j))) * sum(q(:,(j-1)*m+i) .* cos(k*pi*xint'/lx));
                    dk(i,j,k) = (2*wx/(lx*D(i,j))) * sum(q(:,(j-1)*m+i+1) .* cos(k*pi*xint'/lx));
                end
            end
        end
    end
    
    
    
    for j = 1:n-1
        for k = 1:2:N
            ak(m,j,k) = (2*wx/(lx*D(m,j))) * sum(g(:,(m-1)*n+j) .* cos(k*pi*xint'/lx));
            bk(m,j,k) = (2*wx/(lx*D(m,j))) * sum(g(:,(m-1)*n+j+1) .* cos(k*pi*xint'/lx));
            if direction  == 0 
                ck(m,j,k) = (2*wx/(lx*D(m,j))) * sum(q(:,j*m) .* cos(k*pi*xint'/lx));
                dk(m,j,k) = (2*wx/(lx*D(m,j))) * sum(q(:,(j-1)*m+1) .* cos(k*pi*xint'/lx));
            end
        end
    end
    
    
    for i = 1:m-1
        for k = 1:2:N
            ak(i,n,k) = (2*wx/(lx*D(i,n))) * sum(g(:,i*n) .* cos(k*pi*xint'/lx));
            bk(i,n,k) = (2*wx/(lx*D(i,n))) * sum(g(:,(i-1)*n+1) .* cos(k*pi*xint'/lx));
            if direction  == 0
                ck(i,n,k) = (2*wx/(lx*D(i,n))) * sum(q(:,(n-1)*m+i) .* cos(k*pi*xint'/lx));
                dk(i,n,k) = (2*wx/(lx*D(i,n))) * sum(q(:,(n-1)*m+i+1) .* cos(k*pi*xint'/lx));
            end
        end
    end
    
    for k = 1:2:N
        ak(m,n,k) = (2*wx/(lx*D(m,n))) * sum(g(:,m*n) .* cos(k*pi*xint'/lx));
        bk(m,n,k) = (2*wx/(lx*D(m,n))) * sum(g(:,(m-1)*n+1) .* cos(k*pi*xint'/lx));
        if direction == 0
            ck(m,n,k) = (2*wx/(lx*D(m,n))) * sum(q(:,m*n) .* cos(k*pi*xint'/lx));
            dk(m,n,k) = (2*wx/(lx*D(m,n))) * sum(q(:,(n-1)*m+1).* cos(k*pi*xint'/lx));
        end
    end
    
    
    %% Compute effective diffusivity
    for i = 1:m
        for j = 1:n
            if direction  == 0
                Deff(1,1) = Deff(1,1) + D(i,j)* lx^2 * (a0(i,j)+b0(i,j))/4;
                
                for k = 1:2:N
                    Deff(1,1) = Deff(1,1) + D(i,j)*lx^2 * 2 * (ck(i,j,k)-dk(i,j,k))/((k*pi)^2); 
                end
                Deff(2,1) = Deff(2,1) + D(i,j)* lx^2 * (c0(i,j)+d0(i,j))/4;
                
                for k = 1:2:N
                    Deff(2,1) = Deff(2,1) + D(i,j)* lx^2 * 2 * (ak(i,j,k)-bk(i,j,k))/((k*pi)^2); 
                end
                
            else
                Deff(2,2) = Deff(2,2) + D(i,j)* lx^2 * (c0(i,j)+d0(i,j))/4;
   
                for k = 1:2:N
                    Deff(2,2) = Deff(2,2) + D(i,j)* lx^2 * 2 * (ak(i,j,k)-bk(i,j,k))/((k*pi)^2);
                end
            end
        end
    end
end

Deff(1,2) = Deff(2,1);
Deff = Deff/((xn-x0)*(xn-x0));