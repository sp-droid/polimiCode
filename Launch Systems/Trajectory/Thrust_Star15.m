function [T,tburn]=Thrust_Star15(t)

Point1=[0,6672];
Point2=[15,12900];
Point3=[20,4448];
Point4=[30,3559];
Point5=[35,0];
tburn=35;

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
elseif t<Point5(1)
    a=(Point5(2)-Point4(2))/(Point5(1)-Point4(1));
    b=Point4(2)-a*Point4(1);
    T=a*t+b;
else
    T=0;
end

end