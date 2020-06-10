function Aircraft_plot(x,T,le,U);

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