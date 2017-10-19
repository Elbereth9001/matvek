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

#ifndef MV_VEKTOR_HPP
#define MV_VEKTOR_HPP

#include "../utilities/Utilities.hpp"

#include <limits>

#define MV_VEKTORDATA_CTORS(className, paramName, Type, Size) \
/* Default ctor (zero vektor) */ \
MV_API className() : paramName() { for (UInt16 i = 0u; i < Size; ++i) this->paramName[i] = static_cast<Type>(0); } \
/* Array ctor */ \
MV_API className(const std::array<Type, Size>& data) : paramName(data) {} \
template <typename T> /* Single argument ctor */ \
MV_API className(T arg) : paramName{{ std::forward<T>(arg) }} {} \
template <typename T, typename ... Args> /* Variadic ctor */ \
MV_API className(T arg, Args&& ... args) : paramName{{ std::forward<T>(arg), std::forward<Args>(args)... }} {}

namespace mv
{
    namespace detail
    {
        template <typename SizeType, SizeType Size, typename Type>
        class VektorImpl
        {
            static_assert(std::is_same<SizeType, UInt16>::value, "Tried to create non-UInt16 Vektor");
        };
        
        template <UInt16 Size, typename Type>
        struct VektorData
        {
            std::array<Type, Size> _data;
            MV_VEKTORDATA_CTORS(VektorData, _data, Type, Size);
        };
        template <typename Type>
        struct VektorData<2u, Type>
        {
            union
            {
                struct
                {
                    std::array<Type, 2u> _data;
                };
                struct
                {
                    Type x;
                    Type y;
                };
            };
            MV_VEKTORDATA_CTORS(VektorData, _data, Type, 2u);
        };
        template <typename Type>
        struct VektorData<3u, Type>
        {
            union
            {
                struct
                {
                    std::array<Type, 3u> _data;
                };
                struct
                {
                    Type x;
                    Type y;
                    Type z;
                };
            };
            MV_VEKTORDATA_CTORS(VektorData, _data, Type, 3u);
        };
        template <typename Type>
        struct VektorData<4u, Type>
        {
            union
            {
                struct
                {
                    std::array<Type, 4u> _data;
                };
                struct
                {
                    Type x;
                    Type y;
                    Type z;
                    Type w;
                };
            };
            MV_VEKTORDATA_CTORS(VektorData, _data, Type, 4u);
        };
        
        
        template <UInt16 Size, typename Type>
        struct VektorImpl<UInt16, Size, Type> : public detail::VektorData<Size, Type>
        {
            //std::array<Type, Size> _data;
            
            using VData = detail::VektorData<Size, Type>;
            
        public:
            
            // Default ctor (zero vektor)
            MV_API VektorImpl() : VData() {}
            // Array ctor
            MV_API VektorImpl(const std::array<Type, Size>& data) : VData(data) {}
            // Single argument ctor
            template <typename T>
            MV_API VektorImpl(T arg) : VData(std::forward<T>(arg)) {}
            // Variadic ctor
            template <typename T, typename ... Args>
            MV_API VektorImpl(T arg, Args&& ... args) : VData(std::forward<T>(arg), std::forward<Args>(args)... )
            {
                static_assert(sizeof...(Args) == (Size - 1u), "Invalid vektor constructor argument count");
            }
            
            
            //Type conversion
            template <typename T>
            MV_API operator VektorImpl<UInt16, Size, T>()
            {
                std::array<T, Size> arr;
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    arr[i] = static_cast<T>(this->_data[i]);
                }
                return VektorImpl<UInt16, Size, T>(arr);
            }
            //////////////////////////////////////////////////////////
            
            
            // Calls abs to all elements of the vektor
            MV_API void absolute()
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] = math::Abs(this->_data[i]);
                }
            }
            //////////////////////////////////////////////////////////
            
            
            // Fill entire vektor with the same value
            MV_API void fill(const Type value)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] = value;
                }
            }
            //////////////////////////////////////////////////////////
            
            
            // Return vektor data
            MV_API const std::array<Type, Size>& data() const
            {
                return this->_data;
            }
            MV_API std::array<Type, Size>& data()
            {
                return this->_data;
            }
            //////////////////////////////////////////////////////////
            
            
            // Checks if all elements are between min and max
            MV_API bool isWithin(const VektorImpl& min, const VektorImpl& max) const
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    if (min[i] < this->_data[i] || max[i] < this->_data[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            //////////////////////////////////////////////////////////
            
            
            // Checks if all elements are within offset of target
            MV_API bool isWithin(const VektorImpl& target, const Type& offset)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    if ((target[i] - offset) < this->_data[i] || (target[i] + offset) < this->_data[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            
            
            MV_API Type length() const
            {
                return math::Sqrt(lengthSquared());
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API Type lengthSquared() const
            {
                Type l = static_cast<Type>(0);
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    l += this->_data[i] * this->_data[i];
                }
                return l;
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void normalize()
            {
                Type l = length();
                if(math::Epsilon(l))
                {
                    MV_ASSERT(false, "VektorImpl normalization by zero");
                    l = static_cast<Type>(detail::MV_BOUNDARY);
                }
                l = static_cast<Type>(1) / l;
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] *= l;
                }
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void invert()
            {
                static_assert(std::numeric_limits<Type>::is_signed, "Tried to call VektorImpl reversal with unsigned type");
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] = -this->_data[i];
                }
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void scale(const VektorImpl& v)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] *= v[i];
                }
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void scale(const Type s)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] *= s;
                }
            }
            //////////////////////////////////////////////////////////
            
            
            //Returns the number of elements
            MV_API UInt16 size() const
            {
                return Size;
            }
            //////////////////////////////////////////////////////////
            
            
            //Access data in index
            MV_API const Type& at(const UInt16 index) const
            {
                MV_ASSERT(index < Size && index >= 0u, "VektorImpl index out of range");
                return this->_data[index];
            }
            MV_API Type& at(const UInt16 index)
            {
                MV_ASSERT(index < Size && index >= 0u, "VektorImpl index out of range");
                return this->_data[index];
            }
            //////////////////////////////////////////////////////////
            
            
            //Access data in index
            MV_API const Type& operator[](const UInt16 index) const
            {
                return this->_data[index];
            }
            MV_API Type& operator[](const UInt16 index)
            {
                return this->_data[index];
            }
            //////////////////////////////////////////////////////////
            
            template <UInt16 Index = 0u>
            MV_API Type get() const
            {
                static_assert(Index < Size, "Invalid Index for VektorImpl.get<>()");
                return std::get<Index>(this->_data);
            }
            
            
            //Increment
            MV_API VektorImpl& operator+=(const VektorImpl& rhs)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] += rhs._data[i];
                }
                return *this;
            }
            MV_API friend VektorImpl operator+(VektorImpl lhs, const VektorImpl& rhs)
            {
                lhs += rhs;
                return lhs;
            }
            //////////////////////////////////////////////////////////
            
            
            //Substract
            MV_API VektorImpl& operator-=(const VektorImpl& rhs)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    this->_data[i] -= rhs._data[i];
                }
                return *this;
            }
            MV_API friend VektorImpl operator-(VektorImpl lhs, const VektorImpl& rhs)
            {
                lhs -= rhs;
                return lhs;
            }
            //////////////////////////////////////////////////////////
            
            
            //Returns true if elements are smaller is lhs than rhs (elementwise)
            MV_API friend bool operator<(const VektorImpl& lhs, const VektorImpl& rhs)
            {
                for (UInt16 i = 0u; i < Size; ++i)
                {
                    if (lhs._data[i] >= rhs._data[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            //////////////////////////////////////////////////////////
            
            
            //Returns true if elements are larger in lhs than rhs (elementwise)
            MV_API friend bool operator>(const VektorImpl& lhs, const VektorImpl& rhs)
            {
                return rhs < lhs;
            }
            //////////////////////////////////////////////////////////
            
            
            //Returns true if elements are smaller or equal in lhs than rhs (elementwise)
            MV_API friend bool operator<=(const VektorImpl& lhs, const VektorImpl& rhs)
            {
                return !(lhs > rhs);
            }
            //////////////////////////////////////////////////////////
            
            
            //Returns true if elements are larger or equal in lhs than rhs (elementwise)
            MV_API friend bool operator>=(const VektorImpl& lhs, const VektorImpl& rhs)
            {
                return !(lhs < rhs);
            }
            //////////////////////////////////////////////////////////
            
            
            //True if all elements are equal (elementwise)
            MV_API friend bool operator==(const VektorImpl& lhs, const VektorImpl& rhs)
            {
                for(UInt16 i = 0u; i < Size; ++i)
                {
                    if (lhs._data[i] != rhs._data[i])
                    {
                        return false;
                    }
                }
                return true;
            }
            //////////////////////////////////////////////////////////
            
            
            //True if not all elements are equal (elementwise)
            MV_API friend bool operator!=(const VektorImpl& lhs, const VektorImpl& rhs)
            {
                return !(lhs, rhs);
            }
            //////////////////////////////////////////////////////////
            
            
        }; //VektorImpl<UInt16, Size, Type> : VektorImplData
        
    } // detail
    
    #include "../functions/VektorFunctions.hpp"
    
    #ifdef MV_MATRIX_HPP
    #include "../functions/MatrixVektorFunctions.hpp"
    #endif
    #ifdef MV_QUATERNION_HPP
    #include "../functions/QuaternionVektorFunctions.hpp"
    #endif
    
} // mv

#undef MV_VEKTOR_CTORS
#undef MV_VEKTORDATA_CTORS

#endif // MV_VEKTOR_HPP
