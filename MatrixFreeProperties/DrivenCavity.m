% Nonlinear driven cavity problem
% Based off
% http://www.mcs.anl.gov/petsc/petsc-current/src/snes/examples/tutorials/ex19.c.html
% lidvelocity - dimensionless velocity of lid
% grashof     - dimensionless temperature gradent
% prandtl     - dimensionless thermal/momentum diffusity ratio
% Returns:
% A 3-tensor F where
% F(i,j,k) is the (x,y) position of the grid for the kth measurement, where
% F(:,:,1) => u
% F(:,:,1) => v
% F(:,:,1) => omega
% F(:,:,1) => temp
function F = DrivenCavity(xu, xv, xomega, xtemp, mx, my, lidvelocity, grashof, prandtl)
	Fu = zeros(mx,my);
	Fv = zeros(mx,my);
	Fomega = zeros(mx,my);
	Ftemp = zeros(mx,my);
	dhx   = mx-1;
	dhy   = my-1;
	hx    = 1/dhx;
	hy    = 1/dhy;
	hxdhy = hx*dhy;
	hydhx = hy*dhx;
	for j = 1:my
		for i = 1:mx
			%row = i + (j-1) * mx;
			if j == 1                             % Bottom edge
				Fu(i,j)     = xu(i,j);
				Fv(i,j)     = xv(i,j);
				Fomega(i,j) = xomega(i,j) + (xu(i,j+1) - xu(i,j))*dhy;
				Ftemp(i,j)  = xtemp(i,j) - xtemp(i,j+1);
			elseif j == mx                        % Top edge
				Fu(i,j)     = xu(i,j) - lidvelocity;
				Fv(i,j)     = xv(i,j);
				Fomega(i,j) = xomega(i,j) + (xu(i,j) - xu(i,j-1))*dhy;
				Ftemp(i,j)  = xtemp(i,j) - xtemp(i,j-1);
            end

			if i == 1                             % Left edge
				Fu(i,j) = xu(i,j);
				Fv(i,j) = xv(i,j);
				Fomega(i,j) = xomega(i,j) - (xv(i+1,j) - xv(i,j))*dhx;
				Ftemp(i,j) = xtemp(i,j);
			elseif i == mx                        % Right edge
				Fu(i,j) = xu(i,j);
				Fv(i,j) = xv(i,j);
				Fomega(i,j) = xomega(i,j) - (xv(i,j) - xv(i-1,j))*dhx;
				Ftemp(i,j) = xtemp(i,j) - (grashof>0);
			end

			if i > 1 && i < mx && j > 1 && j < mx % Interior
				% convective coefficients for upwinding
				vx  = xu(i,j);
				avx = abs(vx);
				vxp = .5*(vx+avx);
				vxm = .5*(vx-avx);
				vy  = xv(i,j);
				avy = abs(vy);
				vyp = .5*(vy+avy);
				vym = .5*(vy-avy);

				% U velocity
				u         = xu(i,j);
				uxx       = (2.0*u - xu(i-1,j) - xu(i+1,j))*hydhx;
				uyy       = (2.0*u - xu(i,j-1) - xu(i,j+1))*hxdhy;
				Fu(i,j)   = uxx + uyy - .5*(xomega(i,j+1) - xomega(i,j-1))*hx;

				% V velocity
				u         = xv(i,j);
				uxx       = (2.0*u - xv(i-1,j) - xv(i+1,j))*hydhx;
				uyy       = (2.0*u - xv(i,j-1) - xv(i,j+1))*hxdhy;
				Fv(i,j)   = uxx + uyy + .5*(xomega(i+1,j)-xomega(i-1,j))*hy;

				% Omega
				u             = xomega(i,j);
				uxx           = (2.0*u - xomega(i-1,j) - xomega(i+1,j))*hydhx;
				uyy           = (2.0*u - xomega(i,j-1) - xomega(j+1,i))*hxdhy;
				Fomega(i,j)   = uxx + uyy + (vxp*(u - xomega(i-1,j)) + vxm*(xomega(i+1,j) - u))*hy + ...
				                (vyp*(u - xomega(i,j-1)) + vym*(xomega(i,j+1) - u))*hx - ...
				                .5*grashof*(xtemp(i+1,j) - xtemp(i-1,j))*hy;

				% Temperature
				u            = xtemp(i,j);
				uxx          = (2.0*u - xtemp(i-1,j) - xtemp(i+1,j))*hydhx;
				uyy          = (2.0*u - xtemp(i,j-1) - xtemp(i,j+1))*hxdhy;
				Ftemp(i,j)   =  uxx + uyy  + prandtl*((vxp*(u - xtemp(i-1,j)) + vxm*(xtemp(i+1,j) - u))*hy + ...
				                                      (vyp*(u - xtemp(i,j-1)) + vym*(xtemp(i,j+1) - u))*hx);
			end
		end
	end
	F = [reshape(Fu,[],1), reshape(Fv,[],1), reshape(Fomega,[],1), reshape(Ftemp,[],1)];
end
