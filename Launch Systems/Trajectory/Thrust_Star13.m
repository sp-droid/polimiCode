function [T,tburn]=Thrust_Star13(t)

Point1=[0,7784];
Point2=[4,6227];
Point3=[14,9600];
Point4=[15,0];
tburn=15;

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