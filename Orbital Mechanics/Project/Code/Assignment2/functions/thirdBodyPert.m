function a = thirdBodyPert(r, rThirdBody, muThirdBody)

a = muThirdBody*((rThirdBody-r)/norm(rThirdBody-r)^3-rThirdBody/norm(rThirdBody)^3);
end