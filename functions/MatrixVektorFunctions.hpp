/// MATVEK C++ Library
/// Copyright (c) 2017 Eetu Asikainen
/// 
/// This software is provided 'as-is', without any express or implied
/// warranty. In no event will the authors be held liable for any damages
/// arising from the use of this software.
/// 
/// Permission is granted to anyone to use this software for any purpose,
/// including commercial applications, and to alter it and redistribute it
/// freely, subject to the following restrictions:
/// 
/// 1. The origin of this software must not be misrepresented; you must not
///    claim that you wrote the original software. If you use this software
///    in a product, an acknowledgement in the product documentation would be
///    appreciated but is not required.
/// 2. Altered source versions must be plainly marked as such, and must not be
///    misrepresented as being the original software.
/// 3. This notice may not be removed or altered from any source distribution.

#ifndef MV_MVF_HPP
#define MV_MVF_HPP

//Cross product between Mat3 and Vek3
template <typename Type>
static MV_API Mat<3u, 3u, Type> Cross(const Mat<3u, 3u, Type>& m, const Vektor<3u, Type>& v)
{
    return Multiply(m, Tilde(v));
}
//////////////////////////////////////////////////////////


//Cross product between Vek3 and Mat3
template <typename Type>
static MV_API Mat<3u, 3u, Type> Cross(const Vektor<3u, Type>& v, const Mat<3u, 3u, Type>& m)
{
    return Multiply(Tilde(v), m);
}
//////////////////////////////////////////////////////////


//Returns composite axis around which this rotation matrix is happening
template <typename Type>
static MV_API Vektor<3u, Type> GetCompositeRotationAxis(const Mat<3u, 3u, Type>& m)
{
    return Vektor<3u, Type>(
        m.at(1, 2) - m.at(2, 1),
        m.at(2, 0) - m.at(0, 2),
        m.at(0, 1) - m.at(1, 0)
    );
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<4u, 4u, Type> LookAt(const Vektor<3u, Type>& eye, const Vektor<3u, Type>& target, const Vektor<3u, Type>& up)
{
    const Vektor<3u, Type> Z(GetUnitVektor(eye - target));
    const Vektor<3u, Type> X(GetUnitVektor(Cross(up, Z)));
    const Vektor<3u, Type> Y(Cross(Z, X));
    
    return mv::Mat<4u, 4u, Type>(
        X[0], Y[0], Z[0], static_cast<Type>(0),
        X[1], Y[1], Z[1], static_cast<Type>(0),
        X[2], Y[2], Z[2], static_cast<Type>(0),
        -Dot(X, eye), -Dot(Y, eye), -Dot(Z, eye), static_cast<Type>(1)
    );
}
//////////////////////////////////////////////////////////


#if MV_EXPERIMENTAL
//Combined rotation matrix
//Roll, Pitch, Roll (ZXZ)
//Longitude, Latitude, Track direction
template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeCompositeOrbitEphemerisRotationM(const Vektor<3u, Type>& rads)
{
    using math::Cos;
    using math::Sin;
    
    //mu, epsilon, ro (ZXZ)
    const Type cm = Cos(rads[0]);
    const Type sm = Sin(rads[0]);
    const Type ce = Cos(rads[1]);
    const Type se = Sin(rads[1]);
    const Type cr = Cos(rads[2]);
    const Type sr = Sin(rads[2]);
    
    return Mat<3u, 3u, Type>(
        cm*ce, sm*ce, se,
        -cm*se*sr - sm*cr, -sm*se*sr + cm*cr, ce*sr,
        -cm*se*cr + sm*sr, -sm*se*cr - cm*sr, ce*cr
    );
}
//////////////////////////////////////////////////////////
#endif


#if MV_EXPERIMENTAL
//Combined rotation matrix
//Roll, Pitch, Roll (ZXZ)
//to Ascending node, Orbit inclination, Argument of Latitude
template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeCompositeOrbitRotationM(const Vektor<3u, Type>& rads)
{
    using math::Cos;
    using math::Sin;
    
    //omega, iota, nu (ZXZ)
    const Type co = Cos(rads[0]);
    const Type so = Sin(rads[0]);
    const Type ci = Cos(rads[1]);
    const Type si = Sin(rads[1]);
    const Type cn = Cos(rads[2]);
    const Type sn = Sin(rads[2]);
    
    return Mat<3u, 3u, Type>(
        co*cn - so*ci*sn, so*cn + co*ci*sn, si*sn,
        -co*sn - so*ci*cn, -so*sn + co*ci*cn, si*cn,
        so*si, -co*si, ci
    );
}
//////////////////////////////////////////////////////////
#endif


#if MV_EXPERIMENTAL
//Combined rotation matrix
//Pitch, yaw, roll (XYZ)
template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeCompositeRotationM3(const Vektor<3u, Type>& rads)
{
    using math::Cos;
    using math::Sin;
    
    //pitch, yaw, roll
    const Type cp = Cos(rads[0]);
    const Type sp = Sin(rads[0]);
    const Type cy = Cos(rads[1]);
    const Type sy = Sin(rads[1]);
    const Type cr = Cos(rads[2]);
    const Type sr = Sin(rads[2]);
    
    return Mat<3u, 3u, Type>(
        cr*cy, sp*cy, -sy,
        cr*sy*sp - sr*cp, sr*sy*sp + cr*cp, cy*sp,
        cr*sy*cp + sr*sp, sr*sy*cp - cr*sp, cy*cp
    );
}
//////////////////////////////////////////////////////////
#endif

#if MV_EXPERIMENTAL
//Combined rotation matrix
//Pitch, yaw, roll (XYZ)
template <typename Type>
static MV_API Mat<4u, 4u, Type> MakeCompositeRotationM4(const Vektor<3u, Type>& rads)
{
    using math::Cos;
    using math::Sin;
    const Type z = static_cast<Type>(0);
    const Type o = static_cast<Type>(1);
    return Multiply(
        Multiply(
            Mat<4u, 4u, Type>(
                o, z, z, z,
                z, Cos(rads[0]), -Sin(rads[0]), z,
                z, Sin(rads[0]), Cos(rads[0]), z,
                z, z, z, o
            ),
            Mat<4u, 4u, Type>(
                Cos(rads[1]), z, Sin(rads[1]), z,
                z, o, z, z,
                -Sin(rads[1]), z, Cos(rads[1]), z,
                z, z, z, o
            )
        ),
        Mat<4u, 4u, Type>(
            Cos(rads[2]), -Sin(rads[2]), z, z,
            Sin(rads[2]), Cos(rads[2]), z, z,
            z, z, o, z,
            z, z, z, o
        )
    );
}
//////////////////////////////////////////////////////////
#endif


//Vektor v rotated theta radians
//Requires standard reference frame (ijk-aligned)
template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeRotationM(const Vektor<3u, Type>& v, const Type theta)
{
    const Type o = static_cast<Type>(1);
    const Type c = math::Cos(theta);
    const Type s = math::Sin(theta);
    
    return Mat<3u, 3u, Type>(
        v[0] * v[0] + (v[1] * v[1] + v[2] * v[2])*c, v[0] * v[1] * (o - c) - v[2] * s, v[0] * v[2] * (o - c) + v[1] * s,
        v[0] * v[1] * (o - c) + v[2] * s, v[1] * v[1] + (v[2] * v[2] + v[0] * v[0])*c, v[1] * v[2] * (o - c) - v[0] * s,
        v[2] * v[0] * (o - c) - v[1] * s, v[1] * v[2] * (o - c) + v[0] * s, v[2] * v[2] + (v[0] * v[0] + v[1] * v[1])*c
    );
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeScaling3(const Vektor<2u, Type>& s)
{
    const Type z = static_cast<Type>(0);
    return Mat<3u, 3u, Type>(
        s[0], z, z,
        z, s[1], z,
        z, z, static_cast<Type>(1)
    );
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<4u, 4u, Type> MakeScaling4(const Vektor<3u, Type>& s)
{
    const Type z = static_cast<Type>(0);
    return Mat<4u, 4u, Type>(
        s[0], z, z, z,
        z, s[1], z, z,
        z, z, s[2], z,
        z, z, z, static_cast<Type>(1)
    );
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<4u, 4u, Type> MakeScalingViaPoint(const Vektor<3u, Type>& s, const Vektor<3u, Type>& p)
{
    const Type z = static_cast<Type>(0);
    return Mat<4u, 4u, Type>(
        s[0], z, z, p[0] * (1 - s[0]),
        z, s[1], z, p[1] * (1 - s[1]),
        z, z, s[2], p[2] * (1 - s[2]),
        z, z, z, static_cast<Type>(1)
    );
}
//////////////////////////////////////////////////////////


//Returns translation matrix3
template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeTranslation3(const Vektor<2, Type>& v)
{
    const Type z = static_cast<Type>(0);
    const Type o = static_cast<Type>(1);
    return Mat<3u, 3u, Type>(
        o, z, v[0],
        z, o, v[1],
        z, z, o
    );
}
//////////////////////////////////////////////////////////


//Returns translation matrix4
template <typename Type>
static MV_API Mat<4u, 4u, Type> MakeTranslation4(const Vektor<3, Type>& v)
{
    const Type z = static_cast<Type>(0);
    const Type o = static_cast<Type>(1);
    return Mat<4u, 4u, Type>(
        o, z, z, v[0],
        z, o, z, v[1],
        z, z, o, v[2],
        z, z, z, o
    );
}
//////////////////////////////////////////////////////////


//Multiply matrix and vektor
template <UInt8 Rows, UInt16 Size, typename Type>
static MV_API Mat<Rows, 1u, Type> Multiply(const Mat<Rows, (UInt8)Size, Type>& m, const Vektor<Size, Type>& v)
{
    return Multiply(m, ToMatrix(v));
}
//////////////////////////////////////////////////////////


//Multiply vektor and matrix
template <UInt8 Cols, UInt16 Size, typename Type>
static MV_API Mat<1u, Cols, Type> Multiply(const Vektor<Size, Type>& v, const Mat<(UInt8)Size, Cols, Type>& m)
{
    return Multiply(ToMatH(v), m);
}
//////////////////////////////////////////////////////////


//v = axis
template <typename Type>
static MV_API Mat<4u, 4u, Type> RotateAroundAxis(const Vektor<3u, Type>& v, const Type rad)
{
    using math::Cos;
    using math::Sin;
    
    const Type K = static_cast<Type>(1) - Cos(rad);
    return Mat<4u, 4u, Type>(
        (v[0] * v[0]) * K + Cos(rad), v[0] * v[1] * K - v[2] * Sin(rad), v[0] * v[2] * K + v[1] * Sin(rad),
        v[0] * v[1] * K + v[2] * Sin(rad), (v[1] * v[1]) * K + Cos(rad), v[1] * v[2] * K - v[0] * Sin(rad),
        v[0] * v[2] * K - v[1] * Sin(rad), v[1] * v[2] * K + v[0] * Sin(rad), (v[2] * v[2]) * K + Cos(rad)
    );
}


//Vektor to skew-symmetric matrix (for Mat-Vek Cross)
template <typename Type>
static MV_API Mat<3u, 3u, Type> Tilde(const Vektor<3u, Type>& v)
{
    const Type z = static_cast<Type>(0);
    return Mat<3u, 3u, Type>(
        z, -v[2], v[1],
        v[2], z, -v[0],
        -v[1], v[0], z
    );
}
//////////////////////////////////////////////////////////


//Construct mat34 from mat3 and an optional vek3
template <typename Type>
static MV_API Mat<3u, 4u, Type> ToMat34(const Mat<3u, 3u, Type>& m, const Vektor<3u, Type>& v = Vektor<3u, Type>())
{
    return Mat<3u, 4u, Type>(
        m[0], m[1], m[2], v[0],
        m[3], m[4], m[5], v[1],
        m[6], m[7], m[8], v[2]
    );
}
//////////////////////////////////////////////////////////


//Construct mat4 from mat3 and an optional vek3
template <typename Type>
static MV_API Mat<4u, 4u, Type> ToMat4(const Mat<3u, 3u, Type>& m, const Vektor<3u, Type>& v = Vektor<3u, Type>())
{
    const Type z = static_cast<Type>(0);
    return Mat<4u, 4u, Type>(
        m[0], m[1], m[2], v[0],
        m[3], m[4], m[5], v[1],
        m[6], m[7], m[8], v[2],
        z, z, z, static_cast<Type>(1)
    );
}
//////////////////////////////////////////////////////////


//Returns horizontal matrix
template <UInt8 Size, typename Type>
static Mat<1u, Size, Type> ToMatH(const Vektor<Size, Type>& v)
{
    return Mat<1u, Size, Type>(v.data());
}
//////////////////////////////////////////////////////////


//Returns vertical matrix
template <UInt8 Size, typename Type>
static MV_API Mat<Size, 1u, Type> ToMatV(const Vektor<Size, Type>& v)
{
    return Mat<Size, 1u, Type>(v.data());
}
//////////////////////////////////////////////////////////


template <UInt8 MS, typename Type, UInt8 TS = MS>
static MV_API Vektor<TS, Type> ToVek(const Mat<1u, MS, Type>& m)
{
    return ToVek<TS, Type>(Vektor<MS, Type>(m.data()));
}
template <UInt8 MS, typename Type, UInt8 TS = MS>
static MV_API Vektor<TS, Type> ToVek(const Mat<MS, 1u, Type>& m)
{
    return ToVek<TS, Type>(Vektor<MS, Type>(m.data()));
}
//////////////////////////////////////////////////////////


#if !MV_CONSTEXPR
namespace detail
{
    template <bool, UInt8 Rows, UInt8 Cols, typename Type>
    struct ToVektorHelper final
    {
        // Row by row
        MV_DISABLE_CLASS(ToVektorHelper);
        static Vektor<Rows * Cols, Type> ToVektor(const Mat<Rows, Cols, Type>& m)
        {
            std::array<Type, Rows * Cols> arr;
            for (UInt8 i = 0u; i < Rows; ++i)
            {
                for (UInt8 j = 0u; j < Cols; ++i)
                {
                    arr[(i * Cols) + j] = m.at(i, j);
                }
            }
            return Vektor<Rows * Cols, Type>(arr);
        }
    };
    template <UInt8 Rows, UInt8 Cols, typename Type>
    struct ToVektorHelper<false, Rows, Cols, Type> final
    {
        // Column by column
        MV_DISABLE_CLASS(ToVektorHelper);
        static Vektor<Rows * Cols, Type> ToVektor(const Mat<Rows, Cols, Type>& m)
        {
            std::array<Type, Rows * Cols> arr;
            UInt16 index = 0u;
            for (UInt8 j = 0u; j < Cols; ++j)
            {
                for (UInt8 i = 0u; i < Rows; ++i)
                {
                    arr[index] = m.at(i, j);
                    ++index;
                }
            }
            return Vektor<Rows * Cols, Type>(arr);
        }
    };
} // detail
#endif // !MV_CONSTEXPR ToVektorHelper

template <bool FlattenByRow, UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Vektor<Rows * Cols, Type> ToVektor(const Mat<Rows, Cols, Type>& m)
{
    #if MV_CONSTEXPR
    std::array<Type, Rows * Cols> arr;
    if MV_API (FlattenByRow) // Row by row
    {
        for (UInt8 i = 0u; i < Rows; ++i)
        {
            for (UInt8 j = 0u; j < Cols; ++i)
            {
                arr[(i * Cols) + j] = m.at(i, j);
            }
        }
    }
    else // Column by column
    {
        UInt16 index = 0u;
        for (UInt8 j = 0u; j < Cols; ++j)
        {
            for (UInt8 i = 0u; i < Rows; ++i)
            {
                arr[index] = m.at(i, j);
                ++index;
            }
        }
    }
    return Vektor<Rows * Cols, Type>(arr);
    #else
    return detail::ToVektorHelper<FlattenByRow, Rows, Cols, Type>::ToVektor(m);
    #endif
}


#if MV_EXPERIMENTAL
template <typename Type>
static MV_API Vektor<4u, Type> TransformDirection(const Mat<4u, 4u, Type>& m, const Vektor<4u, Type>& v)
{
    return ToVek(Multiply(m, v));
}
//////////////////////////////////////////////////////////
#endif


#if 0
//Search for pattern
//Returns pattern starting coordinates if found, out of bound coordinates otherwise
//Exactness: 0-ignored, 1-exact requirement, 2+ required somewhere, use EV
template <
UInt8 TarRows, UInt8 TarCols,
UInt8 PatRows, UInt8 PatCols,
typename Type, UInt8 EVSize = 0u>
static Mat<1u, 2u, UInt8> SearchForPattern(
    const Mat<TarRows, TarCols, Type>& target,
    const Mat<PatRows, PatCols, Type>& pattern,
    const Mat<PatRows, PatCols, UInt8>* exactness = nullptr,
    const Vektor<EVSize, UInt8>* EVCounts = nullptr,
    const Vektor<EVSize, Type>* EVValues = nullptr
)
{
    static_assert(PatCols > 0u, "Too small pattern to search for");
    static_assert(PatRows > 0u, "Too small pattern to search for");
    static_assert(PatRows < TarRows, "Pattern is larger than searchable area");
    static_assert(PatCols < TarCols, "Pattern is larger than searchable area");
    
    const UInt8 lastMatch = (PatRows * PatCols) - 1u;
    bool lastWasFalse = true;
    
    Mat<PatRows, PatCols, bool> match;
    match.Fill(false);
    
    Vektor<EVSize, UInt8> evmatchcounter;
    evmatchcounter.Fill(0u);
    
    
    Mat<PatRows, PatCols, UInt8> exact;
    if (exactness)
    exact = *exactness;
    else
    exact.Fill(1u);
    
    
    for (UInt8 tarY = 0u; tarY < TarRows; ++tarY)
    {
        for (UInt8 tarX = 0u; tarX < TarCols; ++tarX)
        {
            for (UInt8 patY = 0u; patY < PatRows; ++patY)
            {
                for (UInt8 patX = 0u; patX < PatCols; ++patX)
                {
                    if (tarY + patY >= TarRows || tarX + patX >= TarCols)
                    {
                        break;
                    }
                    else if ((target.at(tarY + patY, tarX + patX) == pattern.at(patY, patX) && exact.at(patY, patX) == 1u) || (exact.at(patY, patX) == 0u))
                    {
                        match.at(patY, patX) = true;
                        lastWasFalse = false;
                    }
                    else
                    {
                        for (UInt8 ev = 0u; ev < EVSize; ++ev)
                        {
                            if (target.at(tarY + patY, tarX + patX) == (*EVValues)[ev])
                            {
                                ++evmatchcounter[ev];
                            }
                        }
                        match.at(patY, patX) = true;
                        lastWasFalse = false;
                        if (lastWasFalse)
                        {
                            break;
                        }
                    }
                }
                if (lastWasFalse)
                {
                    break;
                }
                else
                {
                    lastWasFalse = true;
                }
            }
            if (match[lastMatch] == true)
            {
                if (std::all_of(match.data().begin() + 1u, match.data().end(), [](bool b) {return b == true; }))
                {
                    bool allTrue = true;
                    for (UInt8 m = 0u; m < EVSize; ++m)
                    {
                        if ((*EVCounts)[m] > evmatchcounter[m])
                        {
                            evmatchcounter.Fill(0u);
                            allTrue = false;
                            break;
                        }
                    }
                    if (allTrue)
                    {
                        return Mat<1u, 2u, UInt8>(tarY, tarX);
                    }
                }
            }
            else
            {
                match.Fill(false);
                lastWasFalse = true;
            }
        }
    }
    
    //Pattern not found
    return Mat<1u, 2u, UInt8>(TarRows, TarCols);
}
#endif // if 0 Pattern


#endif // MV_MVF_HPP