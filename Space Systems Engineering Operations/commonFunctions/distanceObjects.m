function distance = distanceObjects(object1, object2)
    distance = sqrt((object1(:,1)-object2(:,1)).^2+(object1(:,2)-object2(:,2)).^2+(object1(:,3)-object2(:,3)).^2);
end