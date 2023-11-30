function [T,tburn]=Thrust_ASAS_13(t)

Point1=[0,11120];
Point2=[6,19000];
Point3=[14,17792];
Point4=[14.3,0];
tburn=14.4;

if t<Point2(1)
    a=(Point2(2)-Point1(2))/(Point2(1)-Point1(1));
    b=Point1(2)-a*Point1(1);
    T=a*t+b;
elseif t<Point3(1)
    a=(Point3(2)-Point2(2))/(Point3(1)-Point2(1));
    b=Point2(2)-a*Point2(1);
    T=a*t+b;
elseif t<Point4(1)
    a=(Point4(2)-Point3(2))/(Point4(1)-Point3(1));
    b=Point3(2)-a*Point3(1);
    T=a*t+b;
else
    T=0;
end

end