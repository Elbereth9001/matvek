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

#ifndef MV_MF_HPP
#define MV_MF_HPP



namespace detail
{
    //Determinants
    template <UInt8, UInt8, typename Type>
    MV_API static Type Det(const Mat<2u, 2u, Type>& m)
    {
        return (m.at(0, 0) * m.at(1, 1)) - (m.at(0, 1) * m.at(1, 0));
    }
    template <UInt8, UInt8, typename Type>
    static MV_API Type Det(const Mat<3u, 3u, Type>& m)
    {
        return (
            m.at(0, 0) * Det<2u, 2u, Type>(GetMinor(m, 0, 0)) -
            m.at(0, 1) * Det<2u, 2u, Type>(GetMinor(m, 0, 1)) +
            m.at(0, 2) * Det<2u, 2u, Type>(GetMinor(m, 0, 2))
        );
    }
    template <UInt8 Rows, UInt8 Cols, typename Type>
    static MV_API Type Det(const Mat<Rows, Cols, Type>& m)
    {
        Type det = static_cast<Type>(0);
        for (UInt8 j = 0u; j < Cols; ++j)
        {
            det += (
                math::Pow<Type>(static_cast<Type>(-1), j) *
                (m.at(0, j) * Det<Rows - 1u, Cols - 1u, Type>(GetMinor(m, 0, j)))
            );
        }
        return det;
    }
    //////////////////////////////////////////////////////////
    
    
    //Cofactors
    template <UInt8, UInt8, typename Type>
    static MV_API Mat<2u, 2u, Type> Cof(const Mat<2u, 2u, Type>& m)
    {
        return Mat<2u, 2u, Type>(m[3], -m[2], -m[1], m[0]);
    }
    template <UInt8 Rows, UInt8 Cols, typename Type>
    static MV_API Mat<Rows, Cols, Type> Cof(const Mat<Rows, Cols, Type>& m)
    {
        Mat<Rows, Cols, Type> mat;
        for (UInt8 i = 0u; i < Rows; ++i)
        {
            for (UInt8 j = 0u; j < Cols; ++j)
            {
                mat.at(i, j) = (
                    math::Pow<Type>(static_cast<Type>(-1), i + j) *
                    Det<Rows - 1u, Cols - 1u, Type>(GetMinor(m, i, j))
                );
            }
        }
        return mat;
    }
    //////////////////////////////////////////////////////////
    
} //Detail


//Calculate adjunct matrix
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Mat<Cols, Rows, Type> Adjunct(const Mat<Rows, Cols, Type>& m)
{
    return Transpose(Cofactor(m));
}
//////////////////////////////////////////////////////////


//Calculate cofactor matrix
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Mat<Rows, Cols, Type> Cofactor(const Mat<Rows, Cols, Type>& m)
{
    return detail::Cof<Rows, Cols, Type>(m);
}
//////////////////////////////////////////////////////////


//Calculate determinant
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Type Determinant(const Mat<Rows, Cols, Type>& m)
{
    return detail::Det<Rows, Cols, Type>(m);
}
//////////////////////////////////////////////////////////


//Returns total rotation around the composite axis
template <typename Type>
static MV_API Type GetCompositeRotationAngle(const Mat<3u, 3u, Type>& m)
{
    return math::ACos((Trace(m) - static_cast<Type>(1)) / static_cast<Type>(2));
}
//////////////////////////////////////////////////////////


//Returns minor
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Mat<Rows - 1u, Cols - 1u, Type> GetMinor(const Mat<Rows, Cols, Type>& parent, const UInt8 exRow, const UInt8 exCol)
{
    MV_ASSERT(Rows > 1u, "Too small matrix to retrieve minor");
    MV_ASSERT(Cols > 1u, "Too small matrix to retrieve minor");
    
    Mat<Rows - 1u, Cols - 1u, Type> m;
    
    for (UInt8 i = 0u, mi = 0u; i < Rows; ++i, ++mi)
    {
        if (i == exRow)
        {
            --mi;
            continue;
        }
        
        for (UInt8 j = 0u, mj = 0u; j < Cols; ++j, ++mj)
        {
            if (j == exCol)
            {
                --mj;
                continue;
            }
            m.at(mi, mj) = parent.at(i, j);
        }
    }
    return m;
}
//////////////////////////////////////////////////////////


//Inverse matrix
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Mat<Rows, Cols, Type> Inverse(const Mat<Rows, Cols, Type>& m)
{
    Type det = Determinant(m);
    if MV_API (math::Epsilon(det))
    {
        MV_ASSERT(false, "Matrix division by zero");
        det = static_cast<Type>(detail::MV_BOUNDARY);
    }
    
    Mat<Cols, Rows, Type> mat(Adjunct(m));
    
    for (auto& itr : mat.data())
    {
        itr /= det;
    }
    return mat;
}
//////////////////////////////////////////////////////////


//Make identity matrix
template <UInt8 Size, typename Type = MV_TYPE>
static MV_API Mat<Size, Size, Type> MakeIdentity()
{
    Mat<Size, Size, Type> m;
    for (UInt8 i = 0u; i < Size; ++i)
    m.at(i, i) = static_cast<Type>(1);
    return m;
}
//////////////////////////////////////////////////////////


//Make diagonal matrix
template <UInt8 Rows, UInt8 Columns = Rows, typename Type = MV_TYPE, typename ... Args>
static MV_API Mat<Rows, Columns, Type> MakeDiagonal(Args&& ... args)
{
    Mat<Rows, Columns, Type> m;
    UInt8 coord = 0u;
    for (const auto& arg : { args... })
    {
        m.at(coord, coord) = arg;
        ++coord;
    }
    return m;
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<4u, 4u, Type> MakeOrtographic(
    const Type left,
    const Type right,
    const Type top,
    const Type bottom,
    const Type near,
    const Type far
)
{
    const Type z = static_cast<Type>(0);
    const Type t = static_cast<Type>(2);
    
    return Mat<4u, 4u, Type>(
        t / (right - left), z, z, z,
        z, t / (top - bottom), z, z,
        z, z, -t / (far - near), z,
        -(right + left) / (right - left), -(top + bottom) / (top - bottom), -(far + near) / (far - near), z
    );
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<4u, 4u, Type> MakePerspective(const Type fovRad, const Type aspectRatio, const Type near, const Type far)
{
    const Type z = static_cast<Type>(0);
    const Type o = static_cast<Type>(1);
    const Type halfTanFovY = math::Tan(fovRad / static_cast<Type>(2));
    
    return Mat<4u, 4u, Type>(
        o / (aspectRatio * halfTanFovY), z, z, z,
        z, o / halfTanFovY, z, z,
        z, z, -(far + near) / (far - near), -o,
        z, z, -(static_cast<Type>(2) * far * near) / (far - near), z
    );
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<3u, 3u, Type> MakeRotationSingle3(const Type rad, const AXIS axis)
{
    using math::Cos;
    using math::Sin;
    
    const Type z = static_cast<Type>(0);
    const Type o = static_cast<Type>(1);
    
    switch (axis)
    {
        case mv::AXIS::X:
        return Mat<3u, 3u, Type>(
            o, z, z,
            z, Cos(rad), -Sin(rad),
            z, Sin(rad), Cos(rad)
        );
        case mv::AXIS::Y:
        return Mat<3u, 3u, Type>(
            Cos(rad), z, Sin(rad),
            z, o, z,
            -Sin(rad), z, Cos(rad)
        );
        case mv::AXIS::Z:
        return Mat<3u, 3u, Type>(
            Cos(rad), -Sin(rad), z,
            Sin(rad), Cos(rad), z,
            z, z, o
        );
        default:
        MV_ASSERT(false, "Invalid rotation axis");
    }
    return Mat<3u, 3u, Type>();
}
//////////////////////////////////////////////////////////


template <typename Type>
static MV_API Mat<4u, 4u, Type> MakeRotationSingle4(const Type rad, const AXIS axis)
{
    using math::Cos;
    using math::Sin;
    
    const Type z = static_cast<Type>(0);
    const Type o = static_cast<Type>(1);
    
    switch (axis)
    {
        case mv::AXIS::X:
        return Mat<4u, 4u, Type>(
            o, z, z, z,
            z, Cos(rad), -Sin(rad), z,
            z, Sin(rad), Cos(rad), z,
            z, z, z, o
        );
        case mv::AXIS::Y:
        return Mat<4u, 4u, Type>(
            Cos(rad), z, Sin(rad), z,
            z, o, z, z,
            -Sin(rad), z, Cos(rad), z,
            z, z, z, o
        );
        case mv::AXIS::Z:
        return Mat<4u, 4u, Type>(
            Cos(rad), -Sin(rad), z, z,
            Sin(rad), Cos(rad), z, z,
            z, z, o, z,
            z, z, z, o
        );
        default:
        MV_ASSERT(false, "Invalid rotation axis");
    }
    return Mat<4u, 4u, Type>();
}
//////////////////////////////////////////////////////////


//Multiply matrices
template <typename Type, UInt8 RowsA, UInt8 ColsARowsB, UInt8 ColsB>
static MV_API Mat<RowsA, ColsB, Type> Multiply(
    const Mat<RowsA, ColsARowsB, Type>& matA,
    const Mat<ColsARowsB, ColsB, Type>& matB
)
{
    Mat<RowsA, ColsB, Type> m;
    for (UInt8 iRow = 0u; iRow < RowsA; ++iRow)
    {
        for (UInt8 iCol = 0u; iCol < ColsB; ++iCol)
        {
            for (UInt8 restARow = 0u; restARow < ColsARowsB; ++restARow)
            {
                m.at(iRow, iCol) += matA.at(iRow, restARow) * matB.at(restARow, iCol);
            }
        }
    }
    return m;
}
//////////////////////////////////////////////////////////


template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Mat<Rows, Cols, Type> Multiply(Mat<Rows, Cols, Type> copy, const Type scalar)
{
    copy.scale(scalar);
    return copy;
}
//////////////////////////////////////////////////////////


#if MV_DEBUG
//Print matrix
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API void Print(const Mat<Rows, Cols, Type>& m)
{
    UInt8 i = 0u;
    for(const auto& itr : m.data())
    {
        if(++i == Cols)
        {
            _msg(itr, '\n');
            i = 0u;
        }
        else
        {
            _msg(itr, ',');
        }
    }
}
//////////////////////////////////////////////////////////
#endif


template <UInt8 NewRows, UInt8 NewCols, typename Type, UInt8 OldRows, UInt8 OldCols>
static MV_API Mat<NewRows, NewCols, Type> ResizeMatrix(const Mat<OldRows, OldCols, Type>& m)
{
    
    static_assert(NewRows == NewCols, "Only square matrices supported are currently supported for resizing");
    static_assert(OldRows == OldCols, "Only square matrices supported are currently supported for resizing");
    
    static_assert(NewRows > 0u, "Cannot resize matrix to zero");
    static_assert(NewCols > 0u, "Cannot resize matrix to zero");
    static_assert(OldRows > 0u, "Cannot resize a matrix with 0 elements");
    static_assert(OldCols > 0u, "Cannot resize a matrix with 0 elements");
    if MV_API (NewRows == OldRows && NewCols == OldCols)
    {
        return m;
    }
    else if MV_API (NewRows == NewCols && OldRows == OldCols && NewRows < OldRows)
    {
        // square smaller
        std::array<Type, NewRows * NewCols> arr;
        for (UInt8 i = 0u; i < NewRows; ++i)
        {
            for (UInt8 j = 0u; j < NewCols; ++j)
            {
                arr[(i * NewCols) + j] = m.at(i, j);
            }
        }
        return Mat<NewRows, NewCols, Type>(arr);
    }
    else if MV_API (NewRows == NewCols && OldRows == OldCols)
    {
        // Square larger
        std::array<Type, NewRows * NewCols> arr;
        for (UInt8 i = 0u; i < OldRows; ++i)
        {
            for (UInt8 j = 0u; j < OldCols; ++j)
            {
                arr[(i * NewCols) + j] = m.at(i, j);
            }
        }
        for (UInt16 i = OldRows * OldCols; i < NewRows * NewCols; ++i)
        {
            arr[i] = static_cast<Type>(0);
        }
        return Mat<NewRows, NewCols, Type>(arr);
    }
    // TODO: non-square matrices
    return Mat<NewRows, NewCols, Type>();
}


//Scale matrix
template <UInt8 Rows, UInt8 Cols, typename Type>
static MV_API Mat<Rows, Cols, Type> Scale(Mat<Rows, Cols, Type> copy, const Type scalar)
{
    copy.scale(scalar);
    return copy;
}
//////////////////////////////////////////////////////////


//Cuts Mat3 down to Mat2
template <typename Type>
static MV_API Mat<2u, 2u, Type> ToMat2(const Mat<3u, 3u, Type>& m)
{
    return Mat<2u, 2u, Type>(
        m.at(0, 0), m.at(0, 1),
        m.at(1, 0), m.at(1, 1)
    );
}
//////////////////////////////////////////////////////////


//Cut mat4 down to mat3
template <typename Type>
static MV_API Mat<3u, 3u, Type> ToMat3(const Mat<4u, 4u, Type>& m)
{
    return Mat<3u, 3u, Type>(
        m[0], m[1], m[2],
        m[4], m[5], m[6],
        m[8], m[9], m[10]
    );
}
//////////////////////////////////////////////////////////


//Get trace of a matrix (sum of diagonal elements)
template <UInt8 RowsCols, typename Type>
static MV_API Type Trace(const Mat<RowsCols, RowsCols, Type>& m)
{
    Type t = static_cast<Type>(0);
    for (UInt8 i = 0u; i < RowsCols; ++i)
    {
        t += m.at(i, i);
    }
    return t;
}
//////////////////////////////////////////////////////////


//Transpose a matrix
template <UInt8 NewRows, UInt8 NewCols, typename Type>
static MV_API Mat<NewRows, NewCols, Type> Transpose(const Mat<NewCols, NewRows, Type>& m)
{
    Mat<NewRows, NewCols, Type> n;
    
    UInt8 oldRow = 0u;
    UInt8 oldCol = 0u;
    for (UInt16 i = 0u; i < NewRows * NewCols; ++i)
    {
        n.at(oldCol, oldRow) = m.at(oldRow, oldCol);
        
        if MV_API (oldCol == (NewRows - 1u))
        {
            oldCol = 0u;
            ++oldRow;
        }
        else
        {
            ++oldCol;
        }
    }
    return n;
}
//////////////////////////////////////////////////////////


#endif // MV_MF_HPP
