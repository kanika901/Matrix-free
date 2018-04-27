% Evaluates
%  -Laplacian u - lambda*exp(u) = 0,  0 < x,y < 1,
% with boundary condition u = 0  for  x = 0, x = 1, y = 0, y = 1.
% lambda is nonlinearity parameter, 0 <= lambda <= 6.81
% http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex5s.c.html
function F = Bratu(x, lambda, mx, my)
    hx = 1/(mx-1);
    hy = 1/(my-1);
    hxdhy = hx / hy;
    hydhx = hy / hx;
    sc = hx*hy*lambda;
    F = zeros(mx*my, 1);
    for i = 1:mx
        for j = 1:my
            row = i + (j-1) * mx;
            if i == 1 || j == 1 || i == mx || j == my
                F(row) = x(row);
				continue
            end
            u = x(row);
            uxx    = (2*u - x(row-1)  - x(row+1))*hydhx;
            uyy    = (2*u - x(row-mx) - x(row+mx))*hxdhy;
            F(row) = uxx + uyy - sc*exp(u);
        end
    end
end
