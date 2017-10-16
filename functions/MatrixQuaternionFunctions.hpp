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

#ifndef MV_MQF_HPP
#define MV_MQF_HPP

#if MV_EXPERIMENTAL
//Experimental
//Returns derivative of a quaternion operator
template <typename Type>
static MV_API Mat<3u, 3u, Type> DerivativeOfQuaternionOperator(const Quat<Type>& q)
{
    const Type t = static_cast<Type>(2);
    const Vektor<3u, Type> w(AngularRateVektor(q._v));
    const Mat<3u, 3u, Type> T(ToMat3(q));
    return Transpose(Mat<3u,3u,Type>(
        Cross(Vektor<3u, Type>(T.at(0, 0)*t, T.at(1, 0)*t, T.at(2, 0)*t), w),
        Cross(Vektor<3u, Type>(T.at(0, 1)*t, T.at(1, 1)*t, T.at(2, 1)*t), w),
        Cross(Vektor<3u, Type>(T.at(0, 2)*t, T.at(1, 2)*t, T.at(2, 2)*t), w)
    ));
}
//////////////////////////////////////////////////////////
#endif


//Converts quaternion to equivalent matrix
template <typename Type>
static MV_API Mat<3u, 3u, Type> ToMat3(const Quat<Type>& q)
{
    const Type o = static_cast<Type>(1);
    const Type t = static_cast<Type>(2);
    const Quat<Type> u(GetUnitQuat(q));
    return Mat<3u, 3u, Type>(
        t*(u.S*u.S + u.X*u.X) - o, t*(u.X*u.Y - u.S*u.Z), t*(u.X*u.Z + u.S*u.Y),
        t*(u.X*u.Y + u.S*u.Z), t*(u.S*u.S + u.Y*u.Y) - o, t*(u.Y*u.Z - u.S*u.X),
        t*(u.X*u.Z - u.S*u.Y), t*(u.Y*u.Z + u.S*u.X), t*(u.S*u.S + u.Z*u.Z) - o
    );
}
//////////////////////////////////////////////////////////


//Returns quaternion
template <typename Type>
static MV_API Quat<Type> ToQuat(const Mat<3u, 3u, Type>& m)
{
    using math::Sqrt;
    
    const MV_TYPE sq = static_cast<MV_TYPE>(m[0] + m.at(1, 1) + m.at(2, 2));
    if MV_API (sq > 0)
    {
        const MV_TYPE s = static_cast<MV_TYPE>(0.5f) / Sqrt(sq + static_cast<MV_TYPE>(1.f));
        return Quat<Type>(
            static_cast<Type>(0.25f / s),
            static_cast<Type>((m.at(2, 1) - m.at(1, 2)) * s),
            static_cast<Type>((m.at(0, 2) - m.at(2, 0)) * s),
            static_cast<Type>((m.at(1, 0) - m.at(0, 1)) * s)
        );
    }
    else
    {
        if MV_API (m[0] > m.at(1, 1) && m[0] > m.at(2, 2))
        {
            const MV_TYPE s = 2.0f * Sqrt(1.0f + m[0] - m.at(1, 1) - m.at(2, 2));
            return Quat<Type>(
                static_cast<Type>((m.at(2, 1) - m.at(1, 2)) / s),
                static_cast<Type>(0.25f * s),
                static_cast<Type>((m.at(0, 1) + m.at(1, 0)) / s),
                static_cast<Type>((m.at(0, 2) + m.at(2, 0)) / s)
            );
        }
        else if MV_API (m.at(1, 1) > m.at(2, 2))
        {
            const MV_TYPE s = 2.0f * Sqrt(1.0f + m.at(1, 1) - m[0] - m.at(2, 2));
            return Quat<Type>(
                static_cast<Type>((m.at(0, 2) - m.at(2, 0)) / s),
                static_cast<Type>((m.at(0, 1) + m.at(1, 0)) / s),
                static_cast<Type>(0.25f * s),
                static_cast<Type>((m.at(1, 2) + m.at(2, 1)) / s)
            );
        }
        else
        {
            const MV_TYPE s = 2.0f * Sqrt(1.0f + m.at(2, 2) - m[0] - m.at(1, 1));
            return Quat<Type>(
                static_cast<Type>((m.at(1, 0) - m.at(0, 1)) / s),
                static_cast<Type>((m.at(0, 2) + m.at(2, 0)) / s),
                static_cast<Type>((m.at(1, 2) + m.at(2, 1)) / s),
                static_cast<Type>(0.25f * s)
            );
        }
    }
}
//////////////////////////////////////////////////////////


#endif // MV_MQF_HPP
