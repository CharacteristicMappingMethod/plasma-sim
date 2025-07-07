function [f] = Advect(f,Efield,grid,dt)

f = Adv_x(f,grid,dt/2);
f = Adv_vx(f,grid,Efield,dt);
f = Adv_x(f,grid,dt/2);

end

function f = Adv_x(f,grid,dt)
X = grid.X;
V = grid.V;

%X_new = mod(X - V*dt,grid.Lx-grid.dx);
X_new = X - V*dt;
for j = 1:size(X_new,1)
    %f(j,:) = interp1(X(j,:),f(j,:),X_new(j,:),grid.method);
    f(j,:) = lagrange1d_local_interp_periodic(X_new(j,:),X(j,:),f(j,:),3);
end

end

function f = Adv_vx(f,grid,Efield,dt)
V = grid.V;
Lv = grid.Lv;
%V_new = mod(V + Efield*dt+Lv,2*grid.Lv-grid.dv)-Lv;
V_new = V + Efield*dt;
for j = 1:size(V,2)
    %f(:,j) = interp1(V(:,j),f(:,j),V_new(:,j),grid.method);
    f(:,j) = lagrange1d_local_interp_periodic(V_new(:,j)+Lv,V(:,j)+Lv,f(:,j),3);
end

end