function [AgeFeatureVector] = CreateAgeFeature(AgeVector)
% Input: AgeVector - Vector with Ages
% Output: AgeFeatureVector - A Vector in the same size of AgeVector with
% ranks depend on age:        Age < 40: 1, 
%                       40 <= Age < 45: 2,
%                       45 <= Age < 50: 3,
%                       50 <= Age < 55: 4,
%                       55 <= Age < 60: 5,
%                       60 <= Age < 65: 6,
%                       65 <= Age < 70: 7,
%                       70 <= Age     : 8,

% Initialize Output Vector
AgeFeatureVector = zeros(size(AgeVector));

% Create Ranking
AgeFeatureVector(AgeVector < 40) = 1;
AgeFeatureVector((AgeVector >= 40) & (AgeVector < 45)) = 2;
AgeFeatureVector((AgeVector >= 45) & (AgeVector < 50))  = 3;
AgeFeatureVector((AgeVector >= 50) & (AgeVector < 55))  = 4;
AgeFeatureVector((AgeVector >= 55) & (AgeVector < 60))  = 5;
AgeFeatureVector((AgeVector >= 60) & (AgeVector < 65))  = 6;
AgeFeatureVector((AgeVector >= 65) & (AgeVector < 70))  = 7;
AgeFeatureVector(AgeVector >= 70) = 8;
end

