function Deff = homogenization_CD(D,x0,xn,x,y0,ym,y,varargin)
% homogenization_CD computes the effective diffusivity tensor Deff for a
% two-dimensional heterogenous domain using a semi-analytical method.

% homogenization_CD solves the two-dimensional homogenization problem on a
% finite domain consisting of an m by mn grid of blocks, where each of the
% mn blocks can have a different diffusivity. The code is applicable to
% domains of size (xn-x0) by (ym-y0) in which each of the mn blocks can be
% of different sizes.

% For domains in which the domain comprises an m by m grid of equally sized
% blocks, please use the code homogenization as it contains optimisations
% applicable to those domains that will decrease the computational time.

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
% x0 = the leftmost point of the domain
% xn = the rightmost point of the domain
% x = vector of coordinates of vertical interfaces
% y0 = the topmost point of the domain
% ym = the bottommost point of the domain
% y = vector of coordinates of horizontal interfaces


% options = options containing the following parameters:
% Nx = the number of abscissa to use in all integrals of vertical interface
% functions
% Ny = the number of abscissa to use in all integrals of horizontal
% interface functions
% N = the number of terms to use in the summation

% OUTPUTS:
% Deff = effective diffusivity tensor

% -------------------------------------------------------------------------
% Check inputs:
% -------------------------------------------------------------------------
if nargin<7
    error('Not enough input arguments.')
elseif nargin == 7
    options = struct;
elseif nargin == 8
    options = varargin{1};
else
    error('Too many input arguments.')
end

[m,n] = size(D);

% Check that all entries of D are positive
if sum(sum(D>0)) ~= m*n
    error('All entries of D must be positive.')
end

% Check that options is a structure
if ~isa(options,'struct')
    error('options must be a structure.')
end

if sum(x == sort(x)) ~= n-1
    error('Coordinates of vertical interfaces must be listed in ascending order.')
end

if sum(y == sort(y)) ~= m-1
    error('Coordinates of horizontal interfaces must be listed in ascending order.')
end

if sum(x0 < x)  ~= n-1
    error('x0 must be less than coordinates of all vertical interfaces.')
end

if sum(y0 < y)  ~= m-1
    error('y0 must be less than coordinates of all vertical interfaces.')
end

if sum(xn > x)  ~= n-1
    error('xn must be greater than coordinates of all vertical interfaces.')
end

if sum(ym > y)  ~= m-1
    error('y0 must be less than coordinates of all vertical interfaces.')
end

if isfield(options,'Nx')
    Nx = options.Nx;
    if round(Nx) ~= Nx && Nx < 1
        error('Nx must be an integer greater than or equal to 1.')
    end
else
    Nx = 16; % default value
end

if isfield(options,'Ny')
    Ny = options.Ny;
    if round(Ny) ~= Ny && Ny < 1
        error('Ny must be an integer greater than or equal to 1.')
    end
else
    Ny = 16; % default value
end

max_N = min(2*Nx-3,2*Ny-3);

if isfield(options,'N')
    N = options.N;
    if round(N) ~= N || N < 1 || N > max_N || N > 100
        error('N must be an integer greater than or equal to 1 and less than min(2*Nx-3,2*Ny-3) and less than 100.')
    elseif N < max_N
        warning('Using N < min(2*Nx-3,2*Ny-3) may produce inaccurate results')
    end
else
    N = min(max_N,100); % default value
end

%% Set up
D = D(end:-1:1,:); % Reverse order of diffusivity matrix as y axis is reversed
xint = zeros(Nx,n); % abscissas in x direction
yint = zeros(Ny,m); % abscissas in y direction
yint(:,1) = linspace(y0+(y(1)-y0)/(2*Ny),y(1)-(y(1)-y0)/(2*Ny),Ny);
yint(:,m) = linspace(y(m-1)+(ym-y(m-1))/(2*Ny),ym-(ym-y(m-1))/(2*Ny),Ny);
for i = 2:m-1
    yint(:,i) = linspace(y(i-1)+(y(i)-y(i-1))/(2*Ny),y(i)-(y(i)-y(i-1))/(2*Ny),Ny);
end
xint(:,1) = linspace(x0+(x(1)-x0)/(2*Nx),x(1)-(x(1)-x0)/(2*Nx),Nx);
xint(:,n) = linspace(x(n-1)+(xn-x(n-1))/(2*Nx),xn-(xn-x(n-1))/(2*Nx),Nx);
for j = 2:n-1
    xint(:,j) = linspace(x(j-1)+(x(j)-x(j-1))/(2*Nx),x(j)-(x(j)-x(j-1))/(2*Nx),Nx);
end
lx = zeros(n,1); % layer widths in x direction
ly = zeros(m,1); % layer widths in y direction
lx(1) = x(1)-x0;
ly(1) = y(1)-y0;
for j = 2:n-1
    lx(j) = x(j) - x(j-1);
end
for i = 2:m-1
    ly(i) = y(i) - y(i-1);
end
lx(n) = xn-x(n-1);
ly(m) = ym - y(m-1);
eigsx = zeros(n,N); % Eigenvalues in x direction
eigsy = zeros(m,N); % Eigenvalues in y direction
for k = 1:N
    for i = 1:m
        eigsy(i,k) = k*pi/ly(i);
    end
    for j = 1:n
        eigsx(j,k) = k*pi/lx(j);
    end
end
bx = zeros(m*n*(Nx+Ny+1),1); % Right hand side vector b
by = zeros(m*n*(Nx+Ny+1),1); % Right hand side vector b
hx = xint(2,:)-xint(1,:); % Distance betwen abscissas in x direction
hy = yint(2,:)-yint(1,:); % Distance betwen abscissas in x direction
midpointx = ones(Nx,1);
wx = midpointx*hx; % Weighting of abscissa used in x direction in midpoint rule
midpointy = ones(Ny,1);
wy = midpointy*hy; % Weighting of abscissa used in y direction in midpoint rule

N_matrix = m*n*Nx*(3*Nx+4*Ny+2) + m*n*Ny*(3*Ny+4*Nx+2) + (m*n-1)*(2*(Nx+Ny))+m*n;

row = zeros(N_matrix,1); % vector of row_start indices
col = zeros(N_matrix,1); % vector of column indices
val = zeros(N_matrix,1); % vector of values

index = 1; % loop through entries in r,c,v


%% Compute constants gamma_{i,j,k} and mu_{i,j,k}
gamma = zeros(m,n,N);
mu = zeros(m,n,N);
for k = 1:N
    for i = 1:m
        for j = 1:n
            gamma(i,j,k) = k*pi*sinh(k*pi*lx(j)/ly(i));
            mu(i,j,k) = k*pi*sinh(k*pi*ly(i)/lx(j));
        end
    end
end

%% Formulate linear system
% Equations are numbered with the first equation being v_{2,1} - v{1,1},
% second equation being v_{2,2} - v_{1,2} etc

% Equation v_{i+1,j}(x^{(p)},y_i) - v_{i,j}(x^{(p)},y_i) = 0
for i = 1:m-1 % Loop through rows
    if i == 1
        bb = y0; % bottom boundary
    else
        bb = y(i-1);
    end
    for j = 1:n % Loop through columns
        if j == 1 % Need to account for x0
            lb = x0; % left boundary
        else
            lb = x(j-1);
        end
        if j == n
            rb = xn; % right boundary
        else
            rb = x(j);
        end
        z8 = m*n*(Nx+Ny) + i*n+j; % indicates position of K_{i+1,j)
        z9 = m*n*(Nx+Ny) + (i-1)*n+j;  % indicates position of K_{i,j)
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
                for k = 1:N
                    val(index) = val(index) + 2*wx(r,j)/(D(i,j)*mu(i,j,k))...
                        * cos(eigsx(j,k)*(xint(p,j)-lb))*cos(eigsx(j,k)*(xint(r,j)-lb));
                end             
                index = index+1;
                
                row(index) = e;
                col(index) = z2;
                val(index) = val(index)  - (wx(r,j)/(2*lx(j))) * (ly(i)/D(i,j) + ly(i+1)/D(i+1,j)); 
                for k = 1:N
                    val(index) = val(index)  - 2*wx(r,j) * cos(eigsx(j,k)*(xint(p,j)-lb))*cos(eigsx(j,k)*(xint(r,j)-lb)) ...
                        * (cosh(eigsx(j,k)*ly(i))/(D(i,j)*mu(i,j,k)) + cosh(eigsx(j,k)*ly(i+1))/(D(i+1,j)*mu(i+1,j,k)));
                end   
                index = index+1;
                
                row(index) = e;
                col(index) = z3;              
                for k = 1:N
                    val(index) = val(index) + 2*wx(r,j)/(D(i+1,j)*mu(i+1,j,k))...
                        * cos(eigsx(j,k)*(xint(p,j)-lb))*cos(eigsx(j,k)*(xint(r,j)-lb));
                end
                index = index+1;  
            end
            
            for r = 1:Ny  % indicates position in vector of abscissas
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
                val(index) = val(index) + wy(r,i)*(xint(p,j)-rb)^2/(2*D(i,j)*lx(j)*ly(i));
                for k = 1:N
                    val(index) = val(index) + 2*wy(r,i)*(-1)^k/(D(i,j)*gamma(i,j,k))... 
                        * cosh(eigsy(i,k)*(xint(p,j)-rb))*cos(eigsy(i,k)*(yint(r,i)-bb));
                end              
                index = index+1;
                
                row(index) = e;
                col(index) = z5;                
                val(index) = val(index) - wy(r,i)*(xint(p,j)-lb)^2/(2*D(i,j)*lx(j)*ly(i)); 
                for k = 1:N
                    val(index) = val(index) - 2*wy(r,i)*(-1)^k/(D(i,j)*gamma(i,j,k))... 
                        * cosh(eigsy(i,k)*(xint(p,j)-lb))*cos(eigsy(i,k)*(yint(r,i)-bb));
                end              
                index = index+1;
                
                row(index) = e;
                col(index) = z6;               
                val(index) = val(index) - wy(r,i+1)*(xint(p,j)-rb)^2/(2*D(i+1,j)*lx(j)*ly(i+1)); 
                for k = 1:N
                    val(index) = val(index) - 2*wy(r,i+1)/(D(i+1,j)*gamma(i+1,j,k))... 
                        * cosh(eigsy(i+1,k)*(xint(p,j)-rb))*cos(eigsy(i+1,k)*(yint(r,i+1)-y(i)));
                end              
                index = index+1;
                
                row(index) = e;
                col(index) = z7;             
                val(index) = val(index) + wy(r,i+1)*(xint(p,j)-lb)^2/(2*D(i+1,j)*lx(j)*ly(i+1)); 
                for k = 1:N
                    val(index) = val(index) + 2*wy(r,i+1)/(D(i+1,j)*gamma(i+1,j,k))...
                        * cosh(eigsy(i+1,k)*(xint(p,j)-lb))*cos(eigsy(i+1,k)*(yint(r,i+1)-y(i)));
                end
                index = index+1;
            end
        end
    end
end

% Equation v_{i,j+1}(x_j,y_^{(p)}) - v_{i,j}(x_j,y^{(p)}) = 0
for i = 1:m % Loop through rows
    if i == 1
        bb = y0; % bottom boundary
    else
        bb = y(i-1);
    end
    if i == m
        ub = ym; % upper boundary
    else
        ub = y(i);
    end
    for j = 1:n-1 % Loop through columns
        if j == 1 % Need to account for x0
            lb = x0; % left boundary
        else
            lb = x(j-1);
        end
        z8 = m*n*(Nx+Ny) + (i-1)*n+j+1; % indicates position of K_{i,j+1)
        z9 = m*n*(Nx+Ny) + (i-1)*n+j;  % indicates position of K_{i,j)
        for p = 1:Ny
            e = (m-1)*n*Nx + ((j-1)*m+i-1)*Ny+p;  % equation number
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
                val(index) = val(index) + wx(r,j)*(yint(p,i)-ub)^2/(2*D(i,j)*lx(j)*ly(i));              
                for k = 1:N
                    val(index) = val(index) + 2*wx(r,j)*(-1)^k/(D(i,j)*mu(i,j,k))... 
                        * cosh(eigsx(j,k)*(yint(p,i)-ub))*cos(eigsx(j,k)*(xint(r,j)-lb));    
                end
                
                index = index+1;                
                row(index) = e;
                col(index) = z2;               
                val(index) = val(index) - wx(r,j)*(yint(p,i)-bb)^2/(2*D(i,j)*lx(j)*ly(i));
                for k = 1:N
                    val(index) = val(index) - 2*wx(r,j)*(-1)^k/(D(i,j)*mu(i,j,k))... 
                        * cosh(eigsx(j,k)*(yint(p,i)-bb))*cos(eigsx(j,k)*(xint(r,j)-lb));
                end
                index = index+1; 
                
                row(index) = e;
                col(index) = z3;              
                val(index) = val(index) - wx(r,j)*(yint(p,i)-ub)^2/(2*D(i,j+1)*lx(j+1)*ly(i));
                for k = 1:N
                    val(index) = val(index) - 2*wx(r,j)/(D(i,j+1)*mu(i,j+1,k))... 
                        * cosh(eigsx(j+1,k)*(yint(p,i)-ub))*cos(eigsx(j+1,k)*(xint(r,j+1)-x(j)));
                end              
                index = index+1;
                
                row(index) = e;
                col(index) = z4; 
                val(index) = val(index) + wx(r,j)*(yint(p,i)-bb)^2/(2*D(i,j+1)*lx(j+1)*ly(i));
                for k = 1:N
                    val(index) = val(index) + 2*wx(r,j)/(D(i,j+1)*mu(i,j+1,k))... 
                        * cosh(eigsx(j+1,k)*(yint(p,i)-bb))*cos(eigsx(j+1,k)*(xint(r,j+1)-x(j)));
                end
                index = index+1;
                
            end
            
            
            
            for r = 1:Ny  % indicates position in vector of abscissas
                z5 = m*n*Nx+((i-1)*n+j-1)*Nx+r;
                z6 = m*n*Nx+((i-1)*n+j)*Nx+r;
                if j == n-1
                    z7 = m*n*Nx+(i-1)*n*Nx+r;
                else
                    z7 = m*n*Nx+((i-1)*n+j+1)*Nx+r;
                end
                
                row(index) = e;
                col(index) = z5;
                for k = 1:N
                    val(index) = val(index) + 2*wy(r,i)/(D(i,j)*gamma(i,j,k))... 
                        * cos(eigsy(i,k)*(yint(p,i)-bb))*cos(eigsy(i,k)*(yint(r,i)-bb));
                end        
                index = index+1;
                
                row(index) = e;
                col(index) = z6;
                val(index) = val(index) - (wy(r,i)/(2*ly(i))) * (lx(j)/D(i,j) + lx(j+1)/D(i,j+1));             
                for k = 1:N
                    val(index) = val(index)  - 2*wy(r,i) * cos(eigsy(i,k)*(yint(p,i)-bb))*cos(eigsy(i,k)*(yint(r,i)-bb)) ...
                        * (cosh(eigsy(i,k)*lx(j))/(D(i,j)*gamma(i,j,k)) + cosh(eigsy(i,k)*lx(j+1))/(D(i,j+1)*gamma(i,j+1,k)));
                end               
                index = index+1;
                
                row(index) = e;
                col(index) = z7;          
                for k = 1:N
                    val(index) = val(index) + 2*wy(r,i)/(D(i,j+1)*gamma(i,j+1,k))... 
                        * cos(eigsy(i,k)*(yint(p,i)-bb))*cos(eigsy(i,k)*(yint(r,i)-bb));
                end              
                index = index+1;            
            end
        end
    end
end

% Equation v_{m,j}(x^{(p)},y_m) - v_{1,j}(x^{(p)},y_0) = y_m - y_0


for j = 1:n % Loop through columns
    if j == 1 % Need to account for x0
        lb = x0; % left boundary
    else
        lb = x(j-1);
    end
    if j == n
        rb = xn; % right boundary
    else
        rb = x(j);
    end
    z8 = m*n*(Nx+Ny) + (m-1)*n+j; % indicates position of K_{m,j)
    z9 = m*n*(Nx+Ny) + j;  % indicates position of K_{1,j)
    for p = 1:Nx
        e = (m-1)*n*Nx+(n-1)*m*Ny+(j-1)*Nx+p; % equation number
        
        row(index) = e;
        col(index) = z8;
        val(index) = 1;
        index = index+1;
        
        row(index) = e;
        col(index) = z9;
        val(index) = -1;
        index = index+1;
        
        by(e) = ym-y0; % right hand side vector
        for r = 1:Nx
            z1 = ((j-1)*m)*Nx+r; % indicates position in vector of abscissas
            z2 = ((j-1)*m+1)*Nx+r;
            z3 = (j*m-1)*Nx+r;
            
            row(index) = e;
            col(index) = z1;  
            val(index) = val(index) + (wx(r,j)/(2*lx(j))) * (ly(1)/D(1,j) + ly(m)/D(m,j)); 
            for k = 1:N   
                val(index) = val(index)  + 2*wx(r,j) * cos(eigsx(j,k)*(xint(p,j)-lb))*cos(eigsx(j,k)*(xint(r,j)-lb)) ...
                    * (cosh(eigsx(j,k)*ly(1))/(D(1,j)*mu(1,j,k)) + cosh(eigsx(j,k)*ly(m))/(D(m,j)*mu(m,j,k))); 
            end           
            index = index+1;
            
            row(index) = e;
            col(index) = z2;
            for k = 1:N
                val(index) = val(index) - 2*wx(r,j)/(D(1,j)*mu(1,j,k))... 
                    * cos(eigsx(j,k)*(xint(p,j)-lb))*cos(eigsx(j,k)*(xint(r,j)-lb));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z3; 
            for k = 1:N 
                val(index) = val(index) - 2*wx(r,j)/(D(m,j)*mu(m,j,k))... 
                    * cos(eigsx(j,k)*(xint(p,j)-lb))*cos(eigsx(j,k)*(xint(r,j)-lb));
            end 
            index = index+1;
        end
        
        for r = 1:Ny  % indicates position in vector of abscissas
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
            val(index) = val(index) + wy(r,1)*(xint(p,j)-rb)^2/(2*D(1,j)*lx(j)*ly(1)); 
            for k = 1:N  
                val(index) = val(index) + 2*wy(r,1)/(D(1,j)*gamma(1,j,k))... 
                    * cosh(eigsy(1,k)*(xint(p,j)-rb))*cos(eigsy(1,k)*(yint(r,1)-y0));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z5;
            val(index) = val(index) - wy(r,1)*(xint(p,j)-lb)^2/(2*D(1,j)*lx(j)*ly(1)); 
            for k = 1:N

                
                val(index) = val(index) - 2*wy(r,1)/(D(1,j)*gamma(1,j,k))... 
                    * cosh(eigsy(1,k)*(xint(p,j)-lb))*cos(eigsy(1,k)*(yint(r,1)-y0));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z6;
            val(index) = val(index) - wy(r,m)*(xint(p,j)-rb)^2/(2*D(m,j)*lx(j)*ly(m));
            for k = 1:N   
                val(index) = val(index) - 2*wy(r,m)*(-1)^k/(D(m,j)*gamma(m,j,k))... 
                    * cosh(eigsy(m,k)*(xint(p,j)-rb))*cos(eigsy(m,k)*(yint(r,m)-y(m-1)));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z7;
            val(index) = val(index) + wy(r,m)*(xint(p,j)-lb)^2/(2*D(m,j)*lx(j)*ly(m));
            for k = 1:N
                val(index) = val(index) + 2*wy(r,m)*(-1)^k/(D(m,j)*gamma(m,j,k))... 
                    * cosh(eigsy(m,k)*(xint(p,j)-lb))*cos(eigsy(m,k)*(yint(r,m)-y(m-1)));
            end
            index = index+1; 
        end
    end
end

% Equation v_{i,n}(x_n,y_^{(p)}) - v_{i,1}(x_0,y^{(p)}) = x_n-x_0
for i = 1:m % Loop through rows
    if i == 1
        bb = y0; % bottom boundary
    else
        bb = y(i-1);
    end
    if i == m
        ub = ym; % upper boundary
    else
        ub = y(i);
    end
    z8 = m*n*(Nx+Ny) + i*n; % indicates position of K_{i,n)
    z9 = m*n*(Nx+Ny) + (i-1)*n+1;  % indicates position of K_{i,1)
    for p = 1:Ny
        e = m*n*Nx+(n-1)*m*Ny+(i-1)*Ny+p; % equation number
        
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
            val(index) = val(index) + wx(r,1)*(yint(p,i)-ub)^2/(2*D(i,1)*lx(1)*ly(i));
            for k = 1:N
                val(index) = val(index) + 2*wx(r,1)/(D(i,1)*mu(i,1,k))... 
                    * cosh(eigsx(1,k)*(yint(p,i)-ub))*cos(eigsx(1,k)*(xint(r,1)-x0));  
            end
            
            index = index+1;
            
            row(index) = e;
            col(index) = z2;
            val(index) = val(index) - wx(r,1)*(yint(p,i)-bb)^2/(2*D(i,1)*lx(1)*ly(i));
            for k = 1:N
                val(index) = val(index) - 2*wx(r,1)/(D(i,1)*mu(i,1,k))... 
                    * cosh(eigsx(1,k)*(yint(p,i)-bb))*cos(eigsx(1,k)*(xint(r,1)-x0));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z3;
            val(index) = val(index) - wx(r,n)*(yint(p,i)-ub)^2/(2*D(i,n)*lx(n)*ly(i)); 
            for k = 1:N
                val(index) = val(index) - 2*wx(r,n)*(-1)^k/(D(i,n)*mu(i,n,k))... 
                    * cosh(eigsx(n,k)*(yint(p,i)-ub))*cos(eigsx(n,k)*(xint(r,n)-x(n-1)));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z4;
            val(index) = val(index) + wx(r,n)*(yint(p,i)-bb)^2/(2*D(i,n)*lx(n)*ly(i)); 
            for k = 1:N
                val(index) = val(index) + 2*wx(r,n)*(-1)^k/(D(i,n)*mu(i,n,k))... 
                    * cosh(eigsx(n,k)*(yint(p,i)-bb))*cos(eigsx(n,k)*(xint(r,n)-x(n-1)));
            end
            index = index+1;
        end
        
        for r = 1:Ny  % indicates position in vector of abscissas
            z5 = m*n*Nx+((i-1)*n)*Nx+r;
            z6 = m*n*Nx+((i-1)*n+1)*Nx+r;
            z7 = m*n*Nx+(i*n-1)*Nx+r;
            
            row(index) = e;
            col(index) = z5;
            
            val(index) = val(index) + (wy(r,i)/(2*ly(i))) * (lx(n)/D(i,n) + lx(1)/D(i,1));
            for k = 1:N
                    val(index) = val(index)  + 2*wy(r,i) * cos(eigsy(i,k)*(yint(p,i)-bb))*cos(eigsy(i,k)*(yint(r,i)-bb)) ...
                        * (cosh(eigsy(i,k)*lx(n))/(D(i,n)*gamma(i,n,k)) + cosh(eigsy(i,k)*lx(1))/(D(i,1)*gamma(i,1,k)));
            end
            index = index+1;
            
            row(index) = e;
            col(index) = z6;
            for k = 1:N
                val(index) = val(index) - 2*wy(r,i)/(D(i,1)*gamma(i,1,k))...
                    * cos(eigsy(i,k)*(yint(p,i)-bb))*cos(eigsy(i,k)*(yint(r,i)-bb));
            end        
            index = index+1;
            
            row(index) = e;
            col(index) = z7;
            for k = 1:N
                val(index) = val(index) - 2*wy(r,i)/(D(i,n)*gamma(i,n,k))...
                    * cos(eigsy(i,k)*(yint(p,i)-bb))*cos(eigsy(i,k)*(yint(r,i)-bb));
            end
            index = index+1;
        end
    end
end

%% Formulate rows of linear system based on zero net flux condition
for i = 1:m-1
    for j = 1:n-1
        e = m*n*(Nx+Ny)+(i-1)*(n-1)+j; % equation number
        for r = 1:Nx
            z1 = ((j-1)*m+i-1)*Nx+r; % indicates position in vector of abscissas
            z2 = ((j-1)*m+i)*Nx+r; % indicates position in vector of abscissas
            
            row(index) = e;
            col(index) = z1;
            val(index) = wx(r,j)/D(i,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z2;
            val(index) = -wx(r,j)/D(i,j);
            index = index+1;
        end
        
        for r = 1:Ny
            z3 = m*n*Nx+((i-1)*n+j-1)*Ny+r; % indicates position in vector of abscissas
            z4 = m*n*Nx+((i-1)*n+j)*Ny+r; % indicates position in vector of abscissas
            
            row(index) = e;
            col(index) = z3;
            val(index) = wy(r,i)/D(i,j);
            index = index+1;
            
            row(index) = e;
            col(index) = z4;
            val(index) = -wy(r,i)/D(i,j);
            index = index+1;
        end
    end
end

for j = 1:n-1
    e = m*n*(Nx+Ny)+(m-1)*(n-1)+j; % equation number
    for r = 1:Nx
        z1 = (j*m-1)*Nx+r; % indicates position in vector of abscissas
        z2 = ((j-1)*m)*Nx+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z1;
        val(index) = wx(r,j)/D(m,j);
        index = index+1;
        
        row(index) = e;
        col(index) = z2;
        val(index) = -wx(r,j)/D(m,j);
        index = index+1;
    end
    
    for r = 1:Ny
        z3 = m*n*Nx+((m-1)*n+j-1)*Ny+r; % indicates position in vector of abscissas
        z4 = m*n*Nx+((m-1)*n+j)*Ny+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z3;
        val(index)  = wy(r,m)/D(m,j);
        index = index+1;
        
        row(index) = e;
        col(index) = z4;
        val(index) = -wy(r,m)/D(m,j);
        index = index+1;
    end
end

for i = 1:m-1
    e = m*n*(Nx+Ny)+m*(n-1)+i; % equation number
    for r = 1:Nx
        z1 = ((n-1)*m+i-1)*Nx+r; % indicates position in vector of abscissas
        z2 = ((n-1)*m+i)*Nx+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z1;
        val(index) = wx(r,n)/D(i,n);
        index = index+1;
        
        row(index) = e;
        col(index) = z2;
        val(index) = -wx(r,n)/D(i,n);
        index = index+1;
        
    end
    
    for r = 1:Ny
        z3 = m*n*Nx+(i*n-1)*Ny+r; % indicates position in vector of abscissas
        z4 = m*n*Nx+((i-1)*n)*Ny+r; % indicates position in vector of abscissas
        
        row(index) = e;
        col(index) = z3;
        val(index) = wy(r,i)/D(i,n);
        index = index+1;
        
        row(index) = e;
        col(index) = z4;
        val(index) = -wy(r,i)/D(i,n);
        index = index+1;
    end
end

e = m*n*(Nx+Ny+1); % equation number

row(N_matrix-m*n+1:N_matrix) = e*ones(m*n,1);
col(N_matrix-m*n+1:N_matrix) = (m*n*(Nx+Ny)+1:m*n*(Nx+Ny+1))';
val(N_matrix-m*n+1:N_matrix) = ones(m*n,1);
bx(e) = (xn+x0)/2;
by(e) = (ym+y0)/2;


%% Solve linear system and assign flux values to vectors

A = sparse(row,col,val,m*n*(Nx+Ny+1),m*n*(Nx+Ny+1));

qg = A\[bx by];
qgx = qg(:,1);
qgy = qg(:,2);

Deff = zeros(2); % effective diffusvity matrix

for direction = 0:1
    
    if direction == 0
        
        q = qgx(1:m*n*Nx);
        g = qgx(m*n*Nx+1:m*n*(Nx+Ny));
        
    else
        q = qgy(1:m*n*Nx);
        g = qgy(m*n*Nx+1:m*n*(Nx+Ny));
    end
    
    
    q = reshape(q,[Nx,m*n]);
    g = reshape(g,[Ny,m*n]);
    
    %% Build solution
    
    % Compute polynomial coefficients
    a0 = zeros(m,n);
    b0 = zeros(m,n);
    c0 = zeros(m,n);
    d0 = zeros(m,n);
    for i = 1:m-1
        for j = 1:n-1
            if direction == 0
                a0(i,j) = (2/ly(i)) * wy(:,i)'*g(:,(i-1)*n+j)/D(i,j);
                b0(i,j) = (2/ly(i)) * wy(:,i)'*g(:,(i-1)*n+j+1)/D(i,j);
            end
            c0(i,j) = (2/lx(j)) * wx(:,j)'*q(:,(j-1)*m+i)/D(i,j);
            d0(i,j) = (2/lx(j)) * wx(:,j)'*q(:,(j-1)*m+i+1)/D(i,j);
        end
    end
    
    for j = 1:n-1
        if direction == 0
            
            a0(m,j) = (2/ly(m)) * wy(:,m)'*g(:,(m-1)*n+j)/D(m,j);
            b0(m,j) = (2/ly(m)) * wy(:,m)'*g(:,(m-1)*n+j+1)/D(m,j);
        end
        c0(m,j) = (2/lx(j)) * wx(:,j)'*q(:,j*m)/D(m,j);
        d0(m,j) = (2/lx(j)) * wx(:,j)'*q(:,(j-1)*m+1)/D(m,j);
    end
    
    
    for i = 1:m-1
        if direction == 0
            
            a0(i,n) = (2/ly(i)) * wy(:,i)'*g(:,i*n)/D(i,n);
            b0(i,n) = (2/ly(i)) * wy(:,i)'*g(:,(i-1)*n+1)/D(i,n);
        end
        c0(i,n) = (2/lx(n)) * wx(:,n)'*q(:,(n-1)*m+i)/D(i,n);
        d0(i,n) = (2/lx(n)) * wx(:,n)'*q(:,(n-1)*m+i+1)/D(i,n);
    end
    if direction == 0
        
        a0(m,n) = (2/ly(m)) * wy(:,m)'*g(:,m*n)/D(m,n);
        b0(m,n) = (2/ly(m)) * wy(:,m)'*g(:,(m-1)*n+1)/D(m,n);
    end
    c0(m,n) = (2/lx(n)) * wx(:,n)'*q(:,m*n)/D(m,n);
    d0(m,n) = (2/lx(n)) * wx(:,n)'*q(:,(n-1)*m+1)/D(m,n);
    
    % Compute trigonometric function coefficients
    
    ak = zeros(m,n,k);
    bk = zeros(m,n,k);
    ck = zeros(m,n,k);
    dk = zeros(m,n,k);
    
    for i = 1:m-1
        if i == 1
            bb = y0;
        else
            bb = y(i-1);
        end
        for j = 1:n-1
            if j == 1
                lb = x0;
            else
                lb = x(j-1);
            end
            for k = 1:2:N
                ak(i,j,k) = (2/ly(i)) * wy(:,i)'*((g(:,(i-1)*n+j)/D(i,j)) .* cos(eigsy(i,k)*(yint(:,i)-bb)));
                bk(i,j,k) = (2/ly(i)) * wy(:,i)'*((g(:,(i-1)*n+j+1)/D(i,j)) .* cos(eigsy(i,k)*(yint(:,i)-bb)));
                if direction  == 0
                    
                    ck(i,j,k) = (2/lx(j)) * wx(:,j)'*((q(:,(j-1)*m+i)/D(i,j)) .* cos(eigsx(j,k)*(xint(:,j)-lb)));
                    dk(i,j,k) = (2/lx(j)) * wx(:,j)'*((q(:,(j-1)*m+i+1)/D(i,j)) .* cos(eigsx(j,k)*(xint(:,j)-lb)));
                end
            end
        end
    end
    
    for j = 1:n-1
        if j == 1
            lb = x0;
        else
            lb = x(j-1);
        end
        for k = 1:2:N
            ak(m,j,k) = (2/ly(m)) * wy(:,m)'*((g(:,(m-1)*n+j)/D(m,j)) .* cos(eigsy(m,k)*(yint(:,m)-y(m-1))));
            bk(m,j,k) = (2/ly(m)) * wy(:,m)'*((g(:,(m-1)*n+j+1)/D(m,j)) .* cos(eigsy(m,k)*(yint(:,m)-y(m-1))));
            if direction  == 0
                ck(m,j,k) = (2/lx(j)) * wx(:,j)'*((q(:,j*m)/D(m,j)) .* cos(eigsx(j,k)*(xint(:,j)-lb)));
                dk(m,j,k) = (2/lx(j)) * wx(:,j)'*((q(:,(j-1)*m+1)/D(m,j)) .* cos(eigsx(j,k)*(xint(:,j)-lb)));
            end
        end
    end
    
    for i = 1:m-1
        if i == 1
            bb = y0;
        else
            bb = y(i-1);
        end
        for k = 1:2:N
            ak(i,n,k) = (2/ly(i)) * wy(:,i)'*((g(:,i*n)/D(i,n)) .* cos(eigsy(i,k)*(yint(:,i)-bb)));
            bk(i,n,k) = (2/ly(i)) * wy(:,i)'*((g(:,(i-1)*n+1)/D(i,n)) .* cos(eigsy(i,k)*(yint(:,i)-bb)));
            if direction  == 0
                
                ck(i,n,k) = (2/lx(n)) * wx(:,n)'*((q(:,(n-1)*m+i)/D(i,n)) .* cos(eigsx(n,k)*(xint(:,n)-x(n-1))));
                dk(i,n,k) = (2/lx(n)) * wx(:,n)'*((q(:,(n-1)*m+i+1)/D(i,n)) .* cos(eigsx(n,k)*(xint(:,n)-x(n-1))));
            end
        end
    end
    
    for k = 1:2:N
        ak(m,n,k) = (2/ly(m)) * wy(:,m)'*((g(:,m*n)/D(m,n)) .* cos(eigsy(m,k)*(yint(:,m)-y(m-1))));
        bk(m,n,k) = (2/ly(m)) * wy(:,m)'*((g(:,(m-1)*n+1)/D(m,n)) .* cos(eigsy(m,k)*(yint(:,m)-y(m-1))));
        if direction == 0
            ck(m,n,k) = (2/lx(n)) * wx(:,n)'*((q(:,m*n)/D(m,n)) .* cos(eigsx(n,k)*(xint(:,n)-x(n-1))));
            dk(m,n,k) = (2/lx(n)) * wx(:,n)'*((q(:,(n-1)*m+1)/D(m,n)) .* cos(eigsx(n,k)*(xint(:,n)-x(n-1))));
        end
    end
    
    
    %% Compute effective diffusivity
    for i = 1:m
        for j = 1:n
            if direction  == 0
                Deff(1,1) = Deff(1,1) + D(i,j)*(lx(j)*ly(i) * (a0(i,j)+b0(i,j)))/4;
                for k = 1:2:N
                    Deff(1,1) = Deff(1,1) + D(i,j)*lx(j)^2 * 2 * (ck(i,j,k)-dk(i,j,k))/((k*pi)^2);
                end
                Deff(2,1) = Deff(2,1) + D(i,j)*(lx(j)*ly(i) * (c0(i,j)+d0(i,j)))/4;
                for k = 1:2:N
                    Deff(2,1) = Deff(2,1) + D(i,j)*ly(i)^2 * 2 * (ak(i,j,k)-bk(i,j,k))/((k*pi)^2);
                end
                
            else
                Deff(2,2) = Deff(2,2) + D(i,j)*(lx(j)*ly(i) * (c0(i,j)+d0(i,j)))/4;
                for k = 1:2:N
                    Deff(2,2) = Deff(2,2) + D(i,j)*ly(i)^2 * 2 * (ak(i,j,k)-bk(i,j,k))/((k*pi)^2);
                end
            end
        end
    end
end
Deff(1,2) = Deff(2,1);
Deff = Deff/((xn-x0)*(ym-y0));
