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

#ifndef MV_MATRIX_HPP
#define MV_MATRIX_HPP

#include "../functions/Utilities.hpp"

namespace mv
{
    
    template <UInt8 Rows, UInt8 Columns = Rows, typename Type = MV_TYPE>
    class Mat
    {
        
        //Matrix data
        union
        {
            std::array<Type, Rows * Columns> __data;
            std::array<std::array<Type, Columns>, Rows> _data;
        };
        
        
    public:
        
        //Ctors
        
        //Default ctor (zero matrix)
        MV_API Mat() : __data()
        {
            for (auto& itr : __data)
            {
                itr = static_cast<Type>(0);
            }
        }
        //Array ctor
        MV_API Mat(const std::array<Type, Rows * Columns>& arr) : __data(arr)
        {
            
        }
        //Single argument ctor
        template <typename T,
        typename std::enable_if<detail::AcceptedType<Type, T>::value>::type* = nullptr
        >
        MV_API Mat(T arg) : __data({ std::forward<T>(arg) })
        {
            
        }
        //Variadic ctor
        template <typename T, typename ... Args,
        typename std::enable_if<detail::AcceptedType<Type, Args...>::value>::type* = nullptr,
        typename std::enable_if<detail::AcceptedType<Type, T>::value>::type* = nullptr
        >
        MV_API Mat(T arg, Args&& ... args) : __data({ std::forward<T>(arg), std::forward<Args>(args)... })
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Dtor
        ~Mat() = default;
        //////////////////////////////////////////////////////////
        
        
        //Move ctor
        MV_API Mat(Mat&& other) : __data(std::move(other.__data))
        {
            //__data = std::move(other.__data);
        }
        //////////////////////////////////////////////////////////
        
        
        //Copy-ctor
        MV_API Mat(const Mat& other) : __data(other.__data)
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Assignment
        MV_API Mat& operator=(Mat&& other)
        {
            __data = other.__data;
            // swap(*this, other);
            return *this;
        }
        //////////////////////////////////////////////////////////

        MV_API Mat& operator=(const Mat& other)
        {
            __data = other.__data;
            return *this;
        }
        
        
        //Swap
        MV_API friend void swap(Mat& first, Mat& second)
        {
            using std::swap;
            swap(first.__data, second.__data);
        }
        //////////////////////////////////////////////////////////
        
        
        //Type conversion
        template <typename T, typename std::enable_if<detail::MV_IsAllowed<T, Type>::value>::type* = nullptr>
        MV_API operator Mat<Rows, Columns, T>()
        {
            std::array<T, Rows * Columns> arr;
            for (UInt16 i = 0u; i < Rows * Columns; ++i)
            {
                arr[i] = static_cast<T>(__data[i]);
            }
            return Mat<Rows, Columns, T>(arr);
        }
        //////////////////////////////////////////////////////////
        
        
        //Return number of elements in matrix
        MV_API UInt32 size() const
        {
            return Rows * Columns;
        }
        //////////////////////////////////////////////////////////
        
        
        //Fill entire matrix with one data value
        MV_API void fill(Type data)
        {
            for (auto& itr : __data)
            {
                itr = data;
            }
        }
        //////////////////////////////////////////////////////////
        
        
        //Get const ref to data
        MV_API const std::array<Type, Rows * Columns>& data() const
        {
            return __data;
        }
        //Get ref to data
        MV_API std::array<Type, Rows * Columns>& data()
        {
            return __data;
        }
        //////////////////////////////////////////////////////////
        
        
        //Get const ref to data
        MV_API const std::array<Type, Columns>& dataRow(const UInt8 row) const
        {
            MV_ASSERT(row < Rows, "Index out of range");
            return _data[row];
        }
        //Get ref to data
        MV_API std::array<Type, Columns>& dataRow(const UInt8 row)
        {
            MV_ASSERT(row < Rows, "Index out of range");
            return _data[row];
        }
        //////////////////////////////////////////////////////////
        
        //Return data from raw index
        MV_API Type operator[](const UInt16 index) const
        {
            MV_ASSERT(index < (Rows * Columns), "Raw index out of range");
            return __data[index];
        }
        MV_API Type& operator[](const UInt16 index)
        {
            MV_ASSERT(index < (Rows * Columns), "Raw index out of range");
            return __data[index];
        }
        MV_API Type at(const UInt16 index) const
        {
            MV_ASSERT(index < (Rows * Columns), "Raw index out of range");
            return __data[index];
        }
        MV_API Type& at(const UInt16 index)
        {
            MV_ASSERT(index < (Rows * Columns), "Raw index out of range");
            return __data[index];
        }
        //////////////////////////////////////////////////////////
        
        
        //Return data from index
        MV_API Type operator[](const Index index) const
        {
            MV_ASSERT(index.Row < Rows, "Row index out of range");
            MV_ASSERT(index.Col < Columns, "Column index out of range");
            return __data[index.Row][index.Col];
        }
        MV_API Type& operator[](const Index index)
        {
            MV_ASSERT(index.Row < Rows, "Row index out of range");
            MV_ASSERT(index.Col < Columns, "Column index out of range");
            return __data[index.Row][index.Col];
        }
        MV_API Type at(const UInt8 row, const UInt8 column) const
        {
            MV_ASSERT(row < Rows, "Row index out of range");
            MV_ASSERT(column < Columns, "Column index out of range");
            return _data[row][column];
        }
        MV_API Type& at(const UInt8 row, const UInt8 column)
        {
            MV_ASSERT(row < Rows, "Row index out of range");
            MV_ASSERT(column < Columns, "Column index out of range");
            return _data[row][column];
        }
        //////////////////////////////////////////////////////////
        
        
        //Increment
        MV_API Mat& operator+=(const Mat& rhs)
        {
            for (UInt16 i = 0u; i < (Rows * Columns); ++i)
            {
                __data[i] += rhs.__data[i];
            }
            return *this;
        }
        MV_API friend Mat operator+(Mat lhs, const Mat& rhs)
        {
            lhs += rhs;
            return lhs;
        }
        //////////////////////////////////////////////////////////
        
        
        //Substract
        MV_API Mat& operator-=(const Mat& rhs)
        {
            for (UInt16 i = 0u; i < (Rows * Columns); ++i)
            {
                __data[i] -= rhs.__data[i];
            }
            return *this;
        }
        MV_API friend Mat operator-(Mat lhs, const Mat& rhs)
        {
            lhs -= rhs;
            return lhs;
        }
        //////////////////////////////////////////////////////////
        
    }; //Mat
    
    #include "../functions/MatrixFunctions.hpp"
    
    #ifdef MV_VEKTOR_HPP
    #include "../functions/MatrixVektorFunctions.hpp"
    #endif
    #ifdef MV_QUATERNION_HPP
    #include "../functions/MatrixQuaternionFunctions.hpp"
    #endif
    
} // mv

#endif // MV_MATRIX_HPP
