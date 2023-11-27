% Trust over time
function [T]=Thrust(t)

Point1=[0,70000];
Point2=[5,85000];
Point3=[20,80000];
Point4=[35,50000];
Point5=[42,40000];
Point6=[52,30000];
Point7=[55,0];

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
elseif t<Point6(1)
    a=(Point6(2)-Point5(2))/(Point6(1)-Point5(1));
    b=Point5(2)-a*Point5(1);
    T=a*t+b;
elseif t<Point7(1)
    a=(Point7(2)-Point6(2))/(Point7(1)-Point6(1));
    b=Point6(2)-a*Point6(1);
    T=a*t+b;
else
    T=0;
end

end