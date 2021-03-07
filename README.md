# FastSineCos
These classes calculate fast mathematical Sine/Cos for limited accuracy. If you need more speed for calculating Sine and reduced accuracy for std::sin() is ok then this class can be 80%-320% faster than std::sin() (according to my current tests - which might be faulty, so if you want to be sure do your own tests :)).

FastSin is a class to calculate fast mathematical Sine for a given angle in radians.
It uses MiniMax polynomial approximation and the degree of the polynomial approximation can be chosen. Smaller degree gives faster results.

Currently two degrees (7 and 9) can be used, but it is easy to add more degrees.

Maximum error for Degree 7: 9.39101e-07
Maximum error for Degree 9: 5.31399e-09

According to my testings FastSin seems to be 80%-340% faster than std::sin(). 

NOTE: FastSin is only fast if you call it so that your consequent angles
are close (about 2Pi) to each others. So for example calling with angles: 1.521, 1.540, 1.600, 1.425. If you pass random angles consequently it should still be
faster than std::sin() but not much. 
So FastSin is good for calculating rotation angles because when when rotating normally consequent
angles are close to each others.

Usage example 1:
```C++
FastSin fastSin1, fastSin2;
auto sin1 = fastSin1(0.268);
auto sin2 = fastSin2(55.689);
```

Usage example 2:
Creating a degree 9 polynomial approximation (more accurate than degree 7) with double type:
```C++
FastSin<double, 9> fastSin3;
auto sin3 = fastSin1(2.2351);
```

Usage example 3:
Creating a float type approximation:
```C++
FastSin<float> fastSin4;
auto sin4 = fastSin1(1.85111);
```
  
This is based on the MinMax values found from:
https://github.com/publik-void/sin-cos-approximations
