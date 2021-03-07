///////////////////////////////////////////////////////////////////////////////
//
// MIT License
//
// Copyright (c) 2021 Juha Kettunen
// Contact: cpekkak ( at ) gmail.com
//
// This algorithm is based on the article:
// "Fast MiniMax Polynomial Approximations of Sine and Cosine"
// https://gist.github.com/publik-void/067f7f2fef32dbe5c27d6e215f824c91
// From that website you can also find more degrees for polynomial approximation.
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so.
//
// I have tested this a lot and I am pretty confident it works but please note
// that it is not yet fully tested so I can not promise it works 100%.
// Especially for extreme values (like huge values, or very small values near zero)
// it is not fully tested.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
///////////////////////////////////////////////////////////////////////////////
//
// Version info
// 07/03/21: Juha Kettunen
// First version. class FastSin added.
//

#ifndef __FAST_SIN__
#define __FAST_SIN__

// FastSin: A class to calculate mathematical sin for a given angle in radians.
// T: The type of the calculations/return value (double/float)
// Degree: the degree of the polynomial approximation used when approximation Sin.
// Can be 7 or 9 (9 is more accurate).
// Maximum error for Degree 7: 9.39101e-07
// Maximum error for Degree 9: 5.31399e-09
// According to my testings FastSin seems to be 80%-340% faster than std::sin(). 
//   NOTE: FastSin is only fast if you call it so that your consequent angles
// are close (about 2*Pi) to each others. So for example calling with angles: 1.521, 1.540, 1.600, 1.425.
// If you pass random angles it should still be faster than std::sin() but not much. 
// So FastSin is good for calculating rotation angles because when rotating normally consequent
// angles are close each others.
//
// Usage example 1:
// FastSin fastSin1, fastSin2;
// auto sin1 = fastSin1(0.268);
// auto sin2 = fastSin2(55.689);
//
// Usage example 2:
// Creating a degree 9 polynomial approximation (more accurate than degree 7) with double type:
// FastSin<double, 9> fastSin3;
// auto sin3 = fastSin1(2.2351);
//
// Usage example 3:
// Creating a float type approximation:
// FastSin<float> fastSin4;
// auto sin4 = fastSin1(1.85111);
//
template<typename T = double, int Degree = 7>
class FastSin
{
public:
    // angle: in radians
    // returns: Mathematical Sine for the angle @angle using template 
    // argument Degree level of polynomial approximation.
    T operator()(T angle);

private:
    // constants used for speedy calculation of the (next) approximation
    inline const static double FAST_SIN_PI{ 3.141592653589793 };
    inline const static double PI_DIV_2{ FAST_SIN_PI / 2.0 };
    inline const static double PI_MULT_3_DIV_2{ FAST_SIN_PI * 3.0 / 2.0 };
    inline const static double PI_MULT_2{ 2.0 * FAST_SIN_PI };
    inline const static double PI_MULT_4{ 4.0 * FAST_SIN_PI };
    // Variables to store information about the previous Sine calculation. These
    // can then be used to calculate fast the next Sine value.
    bool m_hasValidPreviousAngle{ false };
    T m_previousAngle;
    int m_previousFullCyckles;
    double m_previousFullCycklesAngle;
};

template<typename T, int Degree>
T FastSin<T, Degree>::operator()(const T angle)
{
    double angleShort;
    // If previous angle is "near" (near is about 2*Pi) use it as an 
    // advantage to calculate the new angle - it is faster to calculate
    // knowing the information about the last angle values.
    if (m_hasValidPreviousAngle)
    {
        const double diff = angle - m_previousAngle;
        angleShort = angle - m_previousFullCycklesAngle;
        if (diff > 0.0)
        {
            if (angleShort > PI_MULT_2)
            {
                if (angleShort <= PI_MULT_4)
                {
                    ++m_previousFullCyckles;
                    m_previousFullCycklesAngle = m_previousFullCyckles * PI_MULT_2;
                    angleShort = angle - m_previousFullCycklesAngle;
                }
                else
                    m_hasValidPreviousAngle = false;
            }
        }
        else
        {
            if (angleShort < 0.0)
            {
                if (angleShort >= -PI_MULT_2)
                {
                    --m_previousFullCyckles;
                    m_previousFullCycklesAngle = m_previousFullCyckles * PI_MULT_2;
                    angleShort = angle - m_previousFullCycklesAngle;
                }
                else
                    m_hasValidPreviousAngle = false;
            }
        }
    }
    // If we do not have previous angle (to calculate fast), just use a formula
    // which works for all angles but is slower.
    if (!m_hasValidPreviousAngle)
    {
        const double div = angle / PI_MULT_2; // quite slow
        m_previousFullCyckles = div;
        m_previousFullCycklesAngle = m_previousFullCyckles * PI_MULT_2;
        angleShort = (div - static_cast<int>(div)) * PI_MULT_2; // quite slow
    }
    // The polynomial approximation only knows the values from the first quarter section (0 - Pi/2) of the radians unit
    // circle (0 - 2*Pi), so if the angle is on the other 3 quarter sections of the unit circle (Pi/2 - 2*Pi) we need
    // to find the corresponding value (or its negation value) on the first section. Note: If we know all the values 
    // from the first quarter or the unit circle, then we can get the value also for other sections 
    // (Pi/2 - 2*Pi).
    double sign = 1.0;
    if (angleShort > PI_DIV_2 && angleShort <= FAST_SIN_PI)
        angleShort = FAST_SIN_PI - angleShort;
    else if (angleShort > FAST_SIN_PI && angleShort <= PI_MULT_3_DIV_2)
    {
        angleShort = angleShort - FAST_SIN_PI;
        sign = -1.0;
    }
    else if (angleShort > PI_MULT_3_DIV_2 && angleShort <= PI_MULT_2)
    {
        angleShort = PI_MULT_2 - angleShort;
        sign = -1.0;
    }

    const double x1 = angleShort;
    const double x2 = angleShort * angleShort;

    m_previousAngle = angle;
    m_hasValidPreviousAngle = true;

    // Use the polynomial degree according to the template argument @Degree.
    //
    // Below, tests done using a for loop:
    // (*) for (long long i{ -3610000000 }; i < 3610000000; ++i)
    //      angle =  i / 10000000.0;
    //
    // degree 7: x1*(0.999999060898976336474926982596043563 + x2*(-0.166655540927576933646197607200949732 + x2*(0.00831189980138987918776159520367912155 - 0.000184881402886071911033139680005197992*x2)))
    // degree 9: x1*(0.999999994686007336752316120259640318 + x2*(-0.166666566840071513590695269999128453 + x2*(0.00833302513896936729848481553136180314 + x2*(-0.000198074187274269708745741141088641071 + 2.60190306765146018582500885337773154e-6*x2))))
    if constexpr (Degree == 7)
    {
        // degree 7 - Maximum error (*): 9.39101e-07
        return sign * x1 * (static_cast<T>(0.999999060898976) + x2 * (static_cast<T>(-0.166655540927576) +
            x2 * (static_cast<T>(0.00831189980138987) - static_cast<T>(0.000184881402886071 * x2))));
    }
    else if constexpr (Degree == 9)
    {
        // degree 9 - Maximum error (*): 5.31399e-09
        return sign * x1 * (static_cast<T>(0.999999994686007) + x2 * (static_cast<T>(-0.166666566840071) +
            x2 * (static_cast<T>(0.00833302513896936) + x2 * (static_cast<T>(-0.000198074187274269) +
                static_cast<T>(2.601903067651460e-6) * x2))));
    }
}

#endif // __FAST_SIN__
