# Fractals-from-the-Lorenz-Equations

Determine the fractal that arises from using Newton's method to compute the fixed-point solutions of the Lorenz equations.  Use the parameter values   and  Initial values are taken on a grid in the -  plane with always .  For assessment purposes, the computational grid and the graphics code will be given in the Learner Template. To pass the assessment, every pixel in your figure needs to be colored correctly.
(Hint: Some grid points may require as many as 33 Newton iterations to converge while others may require as few as three.  Unfortunately, if you uniformly use 33 Newton iterations at every grid point, the MATLAB Grader may time out.  You can accelerate your code by using a while loop instead of a for loop.) 

r=28; sigma=10; beta=8/3; 
x1=0; y1=0; z1=0;
x2=sqrt(beta*(r-1)); y2=sqrt(beta*(r-1)); z2=r-1;
x3=-sqrt(beta*(r-1)); y3=-sqrt(beta*(r-1)); z3=r-1;
nx=500; nz=500;
xmin=-40; xmax=40; zmin=-40; zmax=40;
x_grid=linspace(xmin,xmax,nx); z_grid=linspace(zmin,zmax,nz);
[X,Z]=meshgrid(x_grid,z_grid);
% Set initial y = 3*sqrt(2) for all points
Y = 3 * sqrt(2) * ones(size(X));

% Now apply Newton's method
tol = 1e-6;  % Tolerância para convergência
max_iter = 40;  % Máximo de iterações para segurança

for i = 1:nx
    for j = 1:nz
        x = X(j,i);
        y = Y(j,i);
        z = Z(j,i);
        
        err = 1;
        iter = 0;
        
        while err > tol && iter < max_iter
            % Funções (f1, f2, f3)
            f1 = sigma*(y - x);
            f2 = x*(r - z) - y;
            f3 = x*y - beta*z;
            
            % Jacobiana (J)
            J = [-sigma, sigma, 0;
                  r - z, -1, -x;
                  y, x, -beta];
            
            F = [f1; f2; f3];
            
            % Atualizar (Delta = -J\F)
            delta = -J\F;
            x = x + delta(1);
            y = y + delta(2);
            z = z + delta(3);
            
            % Atualizar erro
            err = norm(delta);
            iter = iter + 1;
        end
        
        % Atualizar as matrizes após convergência
        X(j,i) = x;
    end
end
eps=1.e-03;
X1 = abs(X-x1) < eps; X2 = abs(X-x2) < eps; X3 = abs(X-x3) < eps;
X4 = ~(X1+X2+X3);
figure; 
map = [1 0 0; 0 1 0; 0 0 1; 0 0 0]; colormap(map); %[red;green;blue;black]
X=(X1+2*X2+3*X3+4*X4); 
image([xmin xmax], [zmin zmax], X); set(gca,'YDir','normal');
xlabel('$x$', 'Interpreter', 'latex', 'FontSize',14);
ylabel('$z$', 'Interpreter', 'latex', 'FontSize',14);
title('Fractal from the Lorenz Equations', 'Interpreter', 'latex','FontSize', 16)  
