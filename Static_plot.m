function U_vector = Static_plot(x, T, le, U, Fext, type)

U = real(U);

U_vector = Convert2vector(x, T, U);

nodes1 = Fext(1:8,1);
nodes2 = Fext(9:end,1);

disp_wing1 = U_vector(nodes1,:);
disp_wing2 = U_vector(nodes2,:);

msg1 = [type ' - Vertical displacement along the span'];
msg2 = [type ' - Horizontal displacement along the span'];

figure;
subplot(2,1,1);
plot(x(nodes2,1),disp_wing2(:,3)*1000); hold on;
plot(x(nodes1,1),disp_wing1(:,3)*1000); 
xlabel('Span [m]','Interpreter','latex','Fontsize',12);
ylabel('Displacement in z-direction [mm]','Interpreter','latex','Fontsize',12);
title(msg1,'Interpreter','latex','Fontsize',14);
legend ({'Leading edge','Trailing edge'},'Interpreter','latex','Fontsize',12)
ylim([ min(disp_wing2(:,3)*1000)-1  max(disp_wing2(:,3)*1000)+1])

subplot(2,1,2);
plot(x(nodes2,1),disp_wing2(:,2)*1000); hold on;
plot(x(nodes1,1),disp_wing1(:,2)*1000); 
xlabel('Span [m]','Interpreter','latex','Fontsize',12);
ylabel('Displacement in y-direction [mm]','Interpreter','latex','Fontsize',12);
title(msg2,'Interpreter','latex','Fontsize',14);
legend ({'Leading edge','Trailing edge'},'Interpreter','latex','Fontsize',12,'location','Southeast')


nelements = size(T,1);

u = U;

for i=1:nelements
        unit(1,i) = u(T(i,1)*(3-2));
        unit(2,i) = u(T(i,1)*(3-1));
        unit(3,i) = u(T(i,1)*(3));
        unit(7,i) = u(T(i,2)*(3-2));
        unit(8,i) = u(T(i,2)*(3-1));
        unit(9,i) = u(T(i,2)*(3));
end

unit = zeros(12,nelements);
n = zeros(1,nelements);
qy = zeros(1,nelements);
qz = qy;
t = qz;
my = zeros(2,nelements);
mz = zeros(2,nelements);

plotWing(x,T,le,U,unit,n,qy,qz,t,my,mz)
