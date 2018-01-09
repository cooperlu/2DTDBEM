function [ DISTC ] = distc(XP,YP,X1,Y1,X2,Y2,SR)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
PROY=abs((X2-X1)*(X2-XP)+(Y2-Y1)*(Y2-YP))/SR;
      DIST=sqrt((X2-XP)^2+(Y2-YP)^2);
if SR<PROY
     DISTC= sqrt(DIST^2-SR*(2*PROY-SR));
     return
else
    DISTC = sqrt(DIST^2-PROY^2);
    return
end
end

