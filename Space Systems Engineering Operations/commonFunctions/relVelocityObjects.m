function relVelocity = relVelocityObjects(object1, object2)
    relVelocity = sqrt((object1(:,4)-object2(:,4)).^2+(object1(:,5)-object2(:,5)).^2+(object1(:,6)-object2(:,6)).^2);
end