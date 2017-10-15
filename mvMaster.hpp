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

/// MATVEK version 1.1 'constexpr'

#pragma once

#pragma warning(push)
#pragma warning(disable: 4201)

#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <limits>
#include <type_traits>

#if defined(MV_DEBUG)
#include <iostream>
#include <string>
#include <sstream>
#endif

#if MV_CONSTEXPR
#define MV_API constexpr
#else
#define MV_API
#endif

#define MV_DISABLE_CLASS(name) \
name() = delete; \
~name() = delete; \
name(name&&) = delete; \
name(const name&) = delete; \
name& operator=(name&&) = delete; \
name& operator=(const name&) = delete;

namespace mv
{
    #if defined(MV_DEBUG)
    
    namespace detail
    {
        template <typename M>
        static std::string _msgImpl(const M& m)
        {
            std::ostringstream ss;
            ss << m;
            return ss.str();
        }
        
        template <typename M, typename ... Ms>
        static std::string _msgImpl(const M& m, const Ms& ... ms)
        {
            return _msgImpl(m) + _msgImpl(ms...);
        }
    } //detail
    
    template <typename ... Ms>
    static void _msg(const Ms& ... ms)
    {
        std::cout << detail::_msgImpl(ms...);
    }
    
    namespace detail
    {
        static void assertion(const bool expression, const char* message)
        {
            if (!expression)
            {
                _msg(message);
                #if defined (_DEBUG)
                __debugbreak();
                #endif
                std::exit(EXIT_FAILURE);
            }
        }
    } //detail
    
    #define MV_ASSERT(expression, message) ::mv::detail::assertion(expression, message)
    #else
    #define MV_ASSERT(expression, message)
    
    template <typename ... Ms>
    static void _msg(const Ms& ...) {}
    #endif
    
    
    using Int8 = signed char;
    using UInt8 = unsigned char;
    
    using Int16 = signed short int;
    using UInt16 = unsigned short int;
    
    using Int32 = signed int;
    using UInt32 = unsigned int;
    
    namespace detail
    {
        
        //Default type for matrices and vectors
        #if defined (MV_USE_LONG_DOUBLE)
        #define MV_TYPE long double
        #elif defined (MV_USE_DOUBLE)
        #define MV_TYPE double
        #else
        #define MV_TYPE float
        #endif
        
        #if defined (MV_BOUNDARY_STRICT)
        static const MV_TYPE MV_BOUNDARY = static_cast<MV_TYPE>(0.0001);
        #elif defined (MV_BOUNDARY_LOOSE)
        static const MV_TYPE MV_BOUNDARY = static_cast<MV_TYPE>(0.01);
        #else
        static const MV_TYPE MV_BOUNDARY = static_cast<MV_TYPE>(0.001);
        #endif
        
        //How strict constructors should be
        #if defined (MV_MAT_STRICT)
        template <typename T1, typename T2>
        using MV_IsAllowed = std::is_same<T1, T2>;
        #else
        template <typename T1, typename T2>
        using MV_IsAllowed = std::is_convertible<T1, T2>;
        #endif
        
        
        //Base
        template <typename, typename...>
        struct AcceptedType{};
        //Last arg
        template <typename Type, typename T>
        struct AcceptedType<Type, T>
        { enum : bool { value = MV_IsAllowed<T, Type>::value }; };
        //Check constructor arguments recursively
        template <typename Type, typename T, typename ... Rest>
        struct AcceptedType<Type, T, Rest...>
        {
            enum : bool { value = MV_IsAllowed<T, Type>::value };
            static_assert(value && AcceptedType<Type, Rest...>::value, "One or more generic matrix constructor parameter is incompatible with desired type.");
        };
        
    } //detail
    
    
    //For accessing matrix via []-operator with two arguments
    struct Index
    {
        UInt8 Row;
        UInt8 Col;
        Index(const UInt8 row, const UInt8 col) : Row(row), Col(col) {}
    };
    
    enum class AXIS : Int8
    {
        Invalid = -1,
        X = 0, Pitch = 0,
        Y = 1, Yaw = 1,
        Z = 2, Roll = 2,
    };
    
    
    //Math
    namespace math
    {
        MV_API const double MV_PI = 3.14159265358979;
        
        template <typename Type>
        inline static Type Abs(const Type val)
        {
            return static_cast<Type>(std::abs(val));
        }
        template<typename Type>
        inline static Type ACos(const Type rad)
        {
            return static_cast<Type>(std::acos(rad));
        }
        template<typename Type>
        inline static Type ASin(const Type rad)
        {
            return static_cast<Type>(std::asin(rad));
        }
        template<typename Type>
        inline static Type ATan(const Type rad)
        {
            return static_cast<Type>(std::atan(rad));
        }
        template<typename Type>
        inline static Type ATan2(const Type y, const Type x)
        {
            return static_cast<Type>(std::atan2(y, x));
        }
        template<typename Type>
        inline static Type Cos(const Type rad)
        {
            return static_cast<Type>(std::cos(rad));
        }
        template <typename Type = MV_TYPE>
        static Type Pi()
        {
            return static_cast<Type>(MV_PI);
        }
        template<typename Type, typename T>
        inline static Type Pow(const Type val, const T exp)
        {
            return static_cast<Type>(std::pow(val, exp));
        }
        template<typename Type>
        inline static Type Sin(const Type rad)
        {
            return static_cast<Type>(std::sin(rad));
        }
        template<typename Type>
        inline static Type Sqrt(const Type val)
        {
            return static_cast<Type>(std::sqrt(val));
        }
        template<typename Type>
        inline static Type Tan(const Type rad)
        {
            return static_cast<Type>(std::tan(rad));
        }
        
        //Returns clamped value
        template <typename Type>
        static Type Clamp(const Type val, const Type min, const Type max)
        {
            return val < min ? min : max < val ? max : val;
        }
        //Returns true if val is close to zero
        template <typename Type>
        static bool Epsilon(const Type val)
        {
            return
            -detail::MV_BOUNDARY < static_cast<MV_TYPE>(val) &&
            static_cast<MV_TYPE>(val) < detail::MV_BOUNDARY;
        }
        //Returns smallest value
        template <typename Type>
        static Type GetLargest(const Type a, const Type b, const Type c)
        {
            return (b < a) ? ((c < a) ? a : c) : ((c < b) ? b : c);
        }
        //Returns index of smallest element
        template <typename Type>
        static UInt8 GetLargest(const Type(&arr)[3u])
        {
            const Type max = GetLargest(arr[0], arr[1], arr[2]);
            return arr[0] == max ? 0u : arr[1] == max ? 1u : 2u;
        }
        //Returns smallest value
        template <typename Type>
        static Type GetSmallest(const Type a, const Type b, const Type c)
        {
            return (a < b) ? ((c < a) ? c : a) : ((c < b) ? c : b);
        }
        //Returns index of the smallest element
        template <typename Type>
        static UInt8 GetSmallest(const Type(&arr)[3u])
        {
            const Type max = GetSmallest(arr[0], arr[1], arr[2]);
            return arr[0] == max ? 0u : arr[1] == max ? 1u : 2u;
        }
        //Runge-Kutta4
        template <typename Type, typename Function>
        static Type RungeKutta4(const Type x, const Type y, const Type dx, Function f)
        {
            const Type t = static_cast<Type>(2);
            const Type k0 = f(x, y);
            const Type k1 = f(x + (dx / t), y + (k0 / t) * dx);
            const Type k2 = f(x + (dx / t), y + (k1 / t) * dx);
            const Type k3 = f(x + dx, y + k2 * dx);
            return y + ((dx * (k0 + t*k1 + t*k2 + k3)) / static_cast<Type>(6));
        }
        //////////////////////////////////////////////////////////
        
    } //math
    
    #define MV_VEKTOR_CTORS_IMPL(className, paramName, Type, Size) \
    /* Array ctor */ \
    MV_API className(const std::array<Type, Size>& data) : paramName(data) {} \
    template <typename T> /* Single argument ctor */ \
    MV_API className(T arg) : paramName({ std::forward<T>(arg) }) {} \
    template <typename T, typename ... Args> /* Variadic ctor */ \
    MV_API className(T arg, Args&& ... args) : paramName({ std::forward<T>(arg), std::forward<Args>(args)... }) {}
    
    #define MV_VEKTORDATA_CTORS(className, paramName, Type, Size) \
    /* Default ctor (zero vektor) */ \
    MV_API className() : paramName() { for (UInt8 i = 0u; i < Size; ++i) this->paramName[i] = static_cast<Type>(0); } \
    MV_VEKTOR_CTORS_IMPL(className, paramName, Type, Size)
    
    
    #define MV_VEKTOR_CTORS(className, paramName, Type, Size) \
    /* Default ctor (zero vektor) */ \
    MV_API className() : paramName() {} \
    MV_VEKTOR_CTORS_IMPL(className, paramName, Type, Size)
    
    namespace detail
    {
        template <typename SizeType, SizeType Size, typename Type>
        class VektorImpl
        {
            static_assert(std::is_same<SizeType, UInt8>::value, "Tried to create non-UInt8 Vektor");
        };
        
        template <UInt8 Size, typename Type>
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
        
        
        //template <UInt8 Size, typename Type>
        //using Vektor = detail::VektorImpl<UInt8, Size, Type>;
        
        template <UInt8 Size, typename Type>
        struct VektorImpl<UInt8, Size, Type> : public detail::VektorData<Size, Type>
        {
            //std::array<Type, Size> _data;
            
            using VData = detail::VektorData<Size, Type>;
            
        public:
            
            MV_VEKTOR_CTORS(VektorImpl, VData, Type, Size);
            
            /*
            //Default ctor (zero vektor)
            MV_API VektorImpl() : VektorData()
            {
                for (UInt8 i = 0u; i < Size; ++i)
                this->_data[i] = static_cast<Type>(0);
            }
            
            //Array ctor
            MV_API VektorImpl(const std::array<Type, Size>& data) :
            _data(data)
            {
                
            }
            
            //Single argument ctor
            template <typename T>
            MV_API VektorImpl(T arg) :
            _data({ std::forward<T>(arg) })
            {
                
            }
            //Variadic ctor
            template <typename T, typename ... Args>
            MV_API VektorImpl(T arg, Args&& ... args) :
            _data({ std::forward<T>(arg), std::forward<Args>(args)... })
            {
                
            }
            //////////////////////////////////////////////////////////
            
            */
            
            //Type conversion
            template <typename T>
            MV_API operator VektorImpl<UInt8, Size, T>()
            {
                std::array<T, Size> arr;
                for (UInt8 i = 0u; i < Size; ++i)
                {
                    arr[i] = static_cast<T>(this->_data[i]);
                }
                return VektorImpl<UInt8, Size, T>(arr);
            }
            //////////////////////////////////////////////////////////
            
            
            // Calls abs to all elements of the vektor
            MV_API void absolute()
            {
                for(auto& itr : this->_data)
                {
                    itr = math::Abs(itr);
                }
            }
            //////////////////////////////////////////////////////////
            
            
            // Fill entire vektor with the same value
            MV_API void fill(const Type value)
            {
                for (auto& itr : this->_data)
                {
                    itr = value;
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
                for (UInt8 i = 0u; i < Size; ++i)
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
                for (UInt8 i = 0u; i < Size; ++i)
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
                return static_cast<Type>(math::Sqrt(lengthSquared()));
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API Type lengthSquared() const
            {
                Type l = static_cast<Type>(0);
                for(const auto& itr : this->_data)
                {
                    l += itr * itr;
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
                for(auto& itr : this->_data)
                {
                    itr /= l;
                }
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void reverse()
            {
                static_assert(std::numeric_limits<Type>::is_signed, "Tried to call VektorImpl reversal with unsigned type");
                for (UInt8 i = 0u; i < Size; ++i)
                {
                    this->_data[i] = -this->_data[i];
                }
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void scale(const VektorImpl& v)
            {
                for (UInt8 i = 0u; i < Size; ++i)
                {
                    this->_data[i] *= v[i];
                }
            }
            //////////////////////////////////////////////////////////
            
            
            MV_API void scale(const Type t)
            {
                for (auto& itr : this->_data)
                {
                    itr *= t;
                }
            }
            //////////////////////////////////////////////////////////
            
            
            //Returns the number of elements
            MV_API UInt32 count() const
            {
                return Size;
            }
            //////////////////////////////////////////////////////////
            
            //Returns the number of elements
            MV_API UInt32 size() const
            {
                return Size;
            }
            //////////////////////////////////////////////////////////
            
            //Access data in index
            MV_API Type at(const UInt8 index) const
            {
                MV_ASSERT(index < Size && index >= 0u, "VektorImpl index out of range");
                return this->_data[index];
            }
            MV_API Type& at(const UInt8 index)
            {
                MV_ASSERT(index < Size && index >= 0u, "VektorImpl index out of range");
                return this->_data[index];
            }
            //////////////////////////////////////////////////////////
            
            
            //Access data in index
            MV_API const Type& operator[](const UInt8 index) const
            {
                return this->_data[index];
            }
            MV_API Type& operator[](const UInt8 index)
            {
                return this->_data[index];
            }
            //////////////////////////////////////////////////////////
            
            template <UInt8 Index = 0u>
            MV_API Type get() const
            {
                static_assert(Index < Size, "Invalid Index for VektorImpl.get<>()");
                return std::get<Index>(this->_data);
            }
            
            
            //Increment
            MV_API VektorImpl& operator+=(const VektorImpl& rhs)
            {
                for (UInt8 i = 0u; i < Size; ++i)
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
                for (UInt8 i = 0u; i < Size; ++i)
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
                for(UInt8 i = 0u; i < Size; ++i)
                {
                    if(lhs._data[i] >= rhs._data[i])
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
                for(UInt8 i = 0u; i < Size; ++i)
                {
                    if(lhs._data[i] != rhs._data[i])
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
            
            
        }; //VektorImpl<UInt8, Size, Type> : VektorImplData
        
    } // detail
    
    template <UInt8 Size, typename Type = MV_TYPE>
    using Vektor = detail::VektorImpl<UInt8, Size, Type>;
    
    
    //Returns copy of v with each element gone through abs
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Absolute(const Vektor<Size, Type>& v)
    {
        std::array<Type, Size> arr(v.data());
        for(auto& itr : arr)
        itr = math::Abs(itr);
        return Vektor<Size, Type>(arr);
    }
    //////////////////////////////////////////////////////////
    
    
    //Rotates target rad radians around axis
    template <typename Type>
    static MV_API Vektor<3u, Type> AxisAngleRotation(const Vektor<3u, Type>& target, const Vektor<3u, Type>& axis, const Type rad)
    {
        const Vektor<3u, Type> n(GetUnitVektor(axis));
        return Vektor<3u, Type>(
            Multiply(target, math::Cos(rad)) +
            Multiply(Multiply(n, Dot(target, n)), (static_cast<Type>(1) - math::Cos(rad))) +
            Multiply(Cross(n, target), math::Sin(rad))
        );
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns a vektor which has the components of two vektors multiplied
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> ComponentProduct(const Vektor<Size, Type>& v1, const Vektor<Size, Type>& v2)
    {
        std::array<Type, Size> arr;
        for(UInt8 i = 0u; i < Size; ++i)
        arr[i] = v1[i] * v2[i];
        return Vektor<Size, Type>(arr);
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns a copy of the vektor
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Copy(const Vektor<Size, Type>& v)
    {
        return Vektor<Size, Type>(v);
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns cross product of two vektors (only applicable on vektors with a size of 3)
    template <typename Type>
    static MV_API Vektor<3u, Type> Cross(const Vektor<3u, Type>& v1, const Vektor<3u, Type>& v2)
    {
        return Vektor<3u, Type>(
            ((v1[1] * v2[2]) - (v1[2] * v2[1])),
            -((v1[0] * v2[2]) - (v1[2] * v2[0])),
            ((v1[0] * v2[1]) - (v1[1] * v2[0]))
        );
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns dot product of two vektors
    template <UInt8 Size, typename Type>
    static MV_API Type Dot(const Vektor<Size, Type>& v1, const Vektor<Size, Type>& v2)
    {
        Type res = static_cast<Type>(0);
        for (UInt8 i = 0u; i < Size; ++i)
        res += v1[i] * v2[i];
        return res;
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns angle between two vektors in radians
    template <UInt8 Size, typename Type>
    static MV_API Type GetAngleBetween(const Vektor<Size, Type>& v1, const Vektor<Size, Type>& v2)
    {
        return static_cast<Type>(math::ACos(Dot(v1, v2) / (Length(v1) * Length(v2))));
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns radians of angle
    template <UInt8 Size, typename Type, typename std::enable_if<Size == 3u>::type* = nullptr>
    static MV_API Type GetDirectionCosine(const Vektor<Size, Type>& v, AXIS axis)
    {
        switch (axis)
        {
            case mv::AXIS::X: return static_cast<Type>(math::ACos(v[0] / Length(v)));
            case mv::AXIS::Y: return static_cast<Type>(math::ACos(v[1] / Length(v)));
            case mv::AXIS::Z: return static_cast<Type>(math::ACos(v[2] / Length(v)));
            default: MV_ASSERT(false, "Invalid rotation angle");
        }
        return static_cast<Type>(0);
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns unitvektor of v as new vektor
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> GetUnitVektor(const Vektor<Size, Type>& v)
    {
        Type l = v.length();
        if(math::Epsilon(l))
        {
            MV_ASSERT(false, "Vektor normalization by zero");
            l = static_cast<Type>(detail::MV_BOUNDARY);
        }
        
        Vektor<Size, Type> u;
        for (UInt8 i = 0u; i < Size; ++i)
        {
            u[i] = v[i] / l;
        }
        return u;
    }
    //////////////////////////////////////////////////////////
    
    
    //Linear interpolation between v1 and v2 in a timeframe t 0...1
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> InterpolateLERP(const Vektor<Size, Type>& v1, const Vektor<Size, Type>& v2, const Type t)
    {
        if (t <= static_cast<Type>(0))
        return v1;
        if (t >= static_cast<Type>(1))
        return v2;
        
        return Vektor<Size, Type>(v1 + Multiply((v2 - v1), t));
    }
    //////////////////////////////////////////////////////////
    
    
    //Spherical linear interpolation between v1 and v2 in a timeframe t 0...1
    //Uses linear interpolation if theta is very small
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> InterpolateSLERP(const Vektor<Size, Type>& v1, const Vektor<Size, Type>& v2, const Type t)
    {
        using math::Sin;
        
        if (t <= static_cast<Type>(0))
        return v1;
        if (t >= static_cast<Type>(1))
        return v2;
        
        const Type theta = GetAngleBetween(v1, v2);
        
        //Use linear interpolation if theta is very small
        return
        static_cast<MV_TYPE>(theta) < static_cast<MV_TYPE>(0.001f) &&
        static_cast<MV_TYPE>(theta) > static_cast<MV_TYPE>(-0.001f) ?
        InterpolateLERP(v1, v2, t) :
        Vektor<Size, Type>(
            Multiply(v1, (Sin((static_cast<Type>(1) - t) * theta)) / (Sin(theta))) +
            Multiply(v2, (Sin(t * theta)) / (Sin(theta)))
        );
        
    }
    //////////////////////////////////////////////////////////
    
    
    // Checks if all elements in target are within min and max
    template <UInt8 Size, typename Type>
    static MV_API bool IsWithin(const Vektor<Size, Type>& target, const Vektor<Size, Type>& min, const Vektor<Size, Type>& max)
    {
        return target.isWithin(min, max);
    }
    //////////////////////////////////////////////////////////
    
    
    // Checks if all elements in target are +- offset of minmax
    template <UInt8 Size, typename Type>
    static MV_API bool IsWithin(const Vektor<Size, Type>& target, const Vektor<Size, Type>& minmax, const Type& offset)
    {
        return target.isWithin(minmax, offset);
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns length of a vektor
    template <UInt8 Size, typename Type>
    static MV_API Type Length(const Vektor<Size, Type>& v)
    {
        return v.length();
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns squared length of a vektor
    template <UInt8 Size, typename Type>
    static MV_API Type LengthSquared(const Vektor<Size, Type>& v)
    {
        return v.lengthSquared();
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns length between two vektors
    template <UInt8 Size, typename Type>
    static MV_API Type Length(const Vektor<Size, Type>& v1, const Vektor<Size, Type>& v2)
    {
        return (v1 - v2).length();
    }
    //////////////////////////////////////////////////////////
    
    
    //Modifies A and B so that they are perpendicular, and return C, which is perpendicular to both A and B. Direction of A will not change.
    template <typename Type>
    static MV_API Vektor<3u, Type> MakeOrthoNormal(Vektor<3u, Type>& a, Vektor<3u, Type>& b)
    {
        a.normalize();
        Vektor<3u, Type> c(Cross(a, b));
        if(math::Epsilon(c.lengthSquared()))
        {
            MV_ASSERT(false, "Tried to orthonormalize vektors that are parallel");
            return c;
        }
        c.normalize();
        b = Cross(c, a);
        return c;
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Multiply(const Vektor<Size, Type>& v, const Type s)
    {
        Vektor<Size, Type> v2(v);
        for (auto& itr : v2.data())
        itr *= s;
        return v2;
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Normalize(const Vektor<Size, Type>& v)
    {
        Vektor<Size, Type> copy(v);
        copy.normalize();
        return copy;
    }
    //////////////////////////////////////////////////////////
    
    
    // Project a onto b
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Project(const Vektor<Size, Type>& a, const Vektor<Size, Type>& b)
    {
        return Multiply(b, Dot(a, b) / b.lengthSquared());
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Reverse(Vektor<Size, Type> copy)
    {
        static_assert(std::numeric_limits<Type>::is_signed, "Tried to call reverse on Vektor of unsigned type");
        copy.reverse();
        return copy;
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Scale(Vektor<Size, Type> v, const Type scalar)
    {
        v.scale(scalar);
        return v;
    }
    //////////////////////////////////////////////////////////
    
    
    //Multiplies every element in v1 with the corresponding element in v2
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> Scale(Vektor<Size, Type> copy, const Vektor<Size, Type>& scalars)
    {
        for (UInt8 i = 0u; i < Size; ++i)
        copy[i] *= scalars[i];
        return copy;
    }
    //////////////////////////////////////////////////////////
    
    
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
        MV_API Mat(Mat&& other) : __data(std::move(other._data))
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Copy-ctor
        MV_API Mat(const Mat& other) : __data(other._data)
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Assignment
        MV_API Mat& operator=(Mat other)
        {
            swap(*this, other);
            return *this;
        }
        //////////////////////////////////////////////////////////
        
        
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
            arr[i] = static_cast<T>(__data[i]);
            return Mat<Rows, Columns, T>(arr);
        }
        //////////////////////////////////////////////////////////
        
        
        //Return number of elements - 1 in matrix
        MV_API UInt32 size() const
        {
            return Rows * Columns - 1u;
        }
        //////////////////////////////////////////////////////////
        
        
        //Fill entire matrix with one data value
        MV_API void fill(Type data)
        {
            for (auto& itr : __data)
            itr = data;
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
            MV_ASSERT(index < (Rows * Cols), "Raw index out of range");
            return __data[index];
        }
        MV_API Type& operator[](const UInt16 index)
        {
            MV_ASSERT(index < (Rows * Cols), "Raw index out of range");
            return __data[index];
        }
        MV_API Type at(const UInt16 index) const
        {
            MV_ASSERT(index < (Rows * Cols), "Raw index out of range");
            return __data[index];
        }
        MV_API Type& at(const UInt16 index)
        {
            MV_ASSERT(index < (Rows * Cols), "Raw index out of range");
            return __data[index];
        }
        //////////////////////////////////////////////////////////
        
        
        //Return data from index
        MV_API Type operator[](const Index index) const
        {
            MV_ASSERT(index.Row < Rows, "Row index out of range");
            MV_ASSERT(index.Col < Cols, "Column index out of range");
            return __data[index.Row][index.Col];
        }
        MV_API Type& operator[](const Index index)
        {
            MV_ASSERT(index.Row < Rows, "Row index out of range");
            MV_ASSERT(index.Col < Cols, "Column index out of range");
            return __data[index.Row][index.Col];
        }
        MV_API Type at(const UInt8 row, const UInt8 column) const
        {
            MV_ASSERT(row < Rows, "Row index out of range");
            MV_ASSERT(column < Cols, "Column index out of range");
            return _data[row][column];
        }
        MV_API Type& at(const UInt8 row, const UInt8 column)
        {
            MV_ASSERT(row < Rows, "Row index out of range");
            MV_ASSERT(column < Cols, "Column index out of range");
            return _data[row][column];
        }
        //////////////////////////////////////////////////////////
        
        
        //Increment
        MV_API Mat& operator+=(const Mat& rhs)
        {
            for (UInt16 i = 0u; i < (Rows * Cols); ++i)
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
            for (UInt16 i = 0u; i < (Rows * Cols); ++i)
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
    
    
    namespace detail
    {
        //Determinants
        template <UInt8, UInt8, typename Type>
        MV_API static Type Det(const Mat<2u, 2u, Type>& m)
        {
            return m.at(0, 0) * m.at(1, 1) - m.at(0, 1) * m.at(1, 0);
        }
        template <UInt8, UInt8, typename Type>
        static MV_API Type Det(const Mat<3u, 3u, Type>& m)
        {
            return
            m.at(0, 0) * Det<2u, 2u, Type>(GetMinor(m, 0, 0)) -
            m.at(0, 1) * Det<2u, 2u, Type>(GetMinor(m, 0, 1)) +
            m.at(0, 2) * Det<2u, 2u, Type>(GetMinor(m, 0, 2));
        }
        template <UInt8 Rows, UInt8 Cols, typename Type>
        static MV_API Type Det(const Mat<Rows, Cols, Type>& m)
        {
            Type det = static_cast<Type>(0);
            for (UInt8 j = 0u; j < Cols; ++j)
            {
                det += math::Pow<Type>(static_cast<Type>(-1), j) *
                (m.at(0, j) * Det<Rows - 1u, Cols - 1u, Type>(GetMinor(m, 0, j)));
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
            for (UInt8 j = 0u; j < Cols; ++j)
            mat.at(i, j) = math::Pow<Type>(static_cast<Type>(-1), i + j) *
            Det<Rows - 1u, Cols - 1u, Type>(GetMinor(m, i, j));
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
    
    
    //Multiply matrices
    template <typename Type, UInt8 RowsA, UInt8 ColsARowsB, UInt8 ColsB>
    static MV_API Mat<RowsA, ColsB, Type> Multiply(
        const Mat<RowsA, ColsARowsB, Type>& matA,
        const Mat<ColsARowsB, ColsB, Type>& matB
    )
    {
        Mat<RowsA, ColsB, Type> m;
        for (UInt8 iRow = 0u; iRow < RowsA; ++iRow)
        for (UInt8 iCol = 0u; iCol < ColsB; ++iCol)
        for (UInt8 restARow = 0u; restARow < ColsARowsB; ++restARow)
        m.at(iRow, iCol) += matA.at(iRow, restARow) * matB.at(restARow, iCol);
        return m;
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Rows, UInt8 Cols, typename Type>
    static MV_API Mat<Rows, Cols, Type> Multiply(const Mat<Rows, Cols, Type>& m, const Type scalar)
    {
        Mat<Rows, Cols, Type> mNew(m);
        for (auto& itr : mNew.data())
        itr *= scalar;
        return mNew;
    }
    //////////////////////////////////////////////////////////
    
    
    //Scale matrix
    template <UInt8 Rows, UInt8 Cols, typename Type>
    static MV_API Mat<Rows, Cols, Type>& Scale(Mat<Rows, Cols, Type>& m, const Type scalar)
    {
        for (auto& itr : m.data())
        itr *= scalar;
        return m;
    }
    //////////////////////////////////////////////////////////
    
    
    //Get trace of a matrix (sum of diagonal elements)
    template <UInt8 RowsCols, typename Type>
    static MV_API Type Trace(const Mat<RowsCols, RowsCols, Type>& m)
    {
        Type t = static_cast<Type>(0);
        for (UInt8 i = 0u; i < RowsCols; ++i)
        t += m.at(i, i);
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
    
    
    template <typename Type = MV_TYPE>
    class Quat
    {
        
    public:
        
        union
        {
            struct
            {
                Vektor<3u, Type> _v;
                Type _s;
            };
            
            struct
            {
                Type X;
                Type Y;
                Type Z;
                
                Type S;
            };
        };
        
        
        MV_API Quat() :
        S(static_cast<Type>(1)),
        X(static_cast<Type>(0)),
        Y(static_cast<Type>(0)),
        Z(static_cast<Type>(0))
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Ctor
        MV_API Quat(Type s, Type x, Type y, Type z) :
        S(s),
        X(x),
        Y(y),
        Z(z)
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Ctor
        MV_API Quat(Type s, const Vektor<3u, Type>& v) :
        _s(s),
        _v(v)
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        //Ctor (pure quaternion)
        MV_API Quat(const Vektor<3u, Type>& v) :
        _s(static_cast<Type>(0)),
        _v(v)
        {
            
        }
        //////////////////////////////////////////////////////////
        
        
        template <typename T>
        MV_API operator Quat<T>()
        {
            return Quat<T>(
                static_cast<T>(S),
                static_cast<T>(X),
                static_cast<T>(Y),
                static_cast<T>(Z)
            );
        }
        //////////////////////////////////////////////////////////
        
        
        MV_API Type length() const
        {
            return math::Sqrt(S*S + X*X + Y*Y + Z*Z);
        }
        //////////////////////////////////////////////////////////
        
        
        MV_API void normalize()
        {
            Type l = length();
            
            if(math::Epsilon(l))
            {
                MV_ASSERT(false, "Quaternion division by zero");
                l = static_cast<Type>(detail::MV_BOUNDARY);
            }
            S /= l;
            X /= l;
            Y /= l;
            Z /= l;
        }
        //////////////////////////////////////////////////////////
        
        
        //Increment
        MV_API Quat& operator+=(const Quat& rhs)
        {
            S += rhs.S;
            X += rhs.X;
            Y += rhs.Y;
            Z += rhs.Z;
            return *this;
        }
        MV_API friend Quat operator+(Quat lhs, const Quat& rhs)
        {
            lhs += rhs;
            return lhs;
        }
        //////////////////////////////////////////////////////////
        
        
        //Substract
        MV_API Quat& operator-=(const Quat& rhs)
        {
            S -= rhs.S;
            X -= rhs.X;
            Y -= rhs.Y;
            Z -= rhs.Z;
            return *this;
        }
        MV_API friend Quat operator-(Quat lhs, const Quat& rhs)
        {
            lhs -= rhs;
            return lhs;
        }
        //////////////////////////////////////////////////////////
        
        
        //Multiply
        MV_API Quat& operator*=(const Quat& rhs)
        {
            const Type s = S * rhs.S - Dot(_v, rhs._v);
            _v = ScaleNew(rhs._v, _s) + ScaleNew(_v, rhs._s) + Cross(_v, rhs._v);
            _s = s;
            return *this;
        }
        MV_API friend Quat operator*(Quat lhs, const Quat& rhs)
        {
            lhs *= rhs;
            return lhs;
        }
        //////////////////////////////////////////////////////////
        
    }; //quat
    
    
    //Get complex conjugate (keeps S, inverts V)
    template <typename Type>
    static MV_API Quat<Type> ComplexConjugate(const Quat<Type>& q)
    {
        return Quat<Type>(q.S, -q.X, -q.Y, -q.Z);
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns total rotation around the composite axis
    template <typename Type>
    static MV_API Type GetCompositeRotationAngle(const Quat<Type>& q)
    {
        return static_cast<Type>(2)*math::ACos(q.S);
    }
    //////////////////////////////////////////////////////////
    
    
    //Divides q by its length, no effect if already normalized
    template <typename Type>
    static MV_API Quat<Type> GetUnitQuat(const Quat<Type>& q)
    {
        Type l = Length(q);
        if MV_API (math::Epsilon(l))
        {
            MV_ASSERT(false, "Quaternion division by zero");
            l = static_cast<Type>(detail::MV_BOUNDARY);
        }
        return Quat<Type>(q.S / l, q.X / l, q.Y / l, q.Z / l);
    }
    //////////////////////////////////////////////////////////
    
    
    //Returns SLERP-interpolated position between positions q1 and q2, from a timeframe of t (0...1)
    template <typename Type>
    static MV_API Quat<Type> InterpolateSLERP(const Quat<Type>& q1, const Quat<Type>& q2, const Type t)
    {
        using math::Sin;
        
        if (t <= static_cast<Type>(0))
        return q1;
        if (t >= static_cast<Type>(1))
        return q2;
        
        const Type theta = math::ACos(
            (q1.S*q2.S + q1.X*q2.X + q1.Y*q2.Y + q1.Z*q2.Z) /
            (q1.length() * q2.length())
        );
        
        return 
        Multiply(q1, (Sin((static_cast<Type>(1) - t) * theta)) / (Sin(theta))) +
        Multiply(q2, (Sin(t * theta)) / (Sin(theta)));
    }
    //////////////////////////////////////////////////////////
    
    
    //Invert a quaternion (get complex conjugate (keeps S, inverts V))
    template <typename Type>
    static MV_API Quat<Type> Inverse(const Quat<Type>& q)
    {
        return Multiply(ComplexConjugate(q), (static_cast<Type>(1) / (Length(q) * Length(q))));
    }
    //////////////////////////////////////////////////////////
    
    
    //Norm
    template <typename Type>
    static MV_API Type Length(const Quat<Type>& q)
    {
        return math::Sqrt(q.S*q.S + q.X*q.X + q.Y*q.Y + q.Z*q.Z);
    }
    //////////////////////////////////////////////////////////
    
    
    template <typename Type>
    static MV_API Quat<Type> Multiply(const Quat<Type>& q1, const Quat<Type>& q2)
    {
        return Quat<Type>(
            (q1.S * q2.S) - Dot(q1._v, q2._v),
            Multiply(q2._v, q1._s) + Multiply(q1._v, q2._s) + Cross(q1._v, q2._v)
        );
    }
    //////////////////////////////////////////////////////////
    
    
    template <typename Type>
    static MV_API Quat<Type> Multiply(const Quat<Type>& q, const Type scalar)
    {
        return Quat<Type>(
            q.S * scalar,
            q.X * scalar,
            q.Y * scalar,
            q.Z * scalar
        );
    }
    //////////////////////////////////////////////////////////
    
    
    template <typename Type>
    static MV_API Quat<Type> Scale(Quat<Type> copy, const Type scalar)
    {
        copy.scale(scalar);
        return copy;
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Size, typename Type>
    static MV_API Type AngularRateAboutVektor(const Vektor<Size, Type>& v)
    {
        MV_ASSERT(false, "Math under construction");
        return static_cast<Type>(1);
    }
    //////////////////////////////////////////////////////////
    
    
    template <UInt8 Size, typename Type>
    static MV_API Vektor<Size, Type> AngularRateVektor(const Vektor<Size, Type>& v)
    {
        return Multiply(v, AngularRateAboutVektor(v));
    }
    //////////////////////////////////////////////////////////
    
    
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
    
    
    #if MV_EXPERIMENTAL
    //Experimental
    //Returns derivative of quaternion
    template <typename Type>
    static MV_API Quat<Type> DerivativeOfQuaternion(const Quat<Type>& q)
    {
        return Multiply(q, AngularRateVektor(q._v));
    }
    //////////////////////////////////////////////////////////
    #endif
    
    #if MV_EXPERIMENTAL
    //Experimental
    //Returns derivative of quaternions conjugate (don't pass conjugate in)
    template <typename Type>
    static MV_API Quat<Type> DerivativeOfQuaternionConjugate(const Quat<Type>& q)
    {
        return Multiply(ComplexConjugate(q), -AngularRateVektor(q._v));
    }
    //////////////////////////////////////////////////////////
    #endif
    
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
    
    #if MV_EXPERIMENTAL
    //Experimental
    //Decomposes result into return value
    //Un-normalizes result (does this have any applications)
    template <typename Type>
    static MV_API Quat<Type> FactorsDistinctPrincipalAxis(const Quat<Type>& result)
    {
        const Quat<Type> p(GetUnitQuat(result));
        const Type A = p.S*p.X + p.Y*p.Z;
        const Type B = -p.S*p.S + p.Y*p.Y;
        const Type D = p.X*p.X - (p.Z*p.Z);
        const Type theta = math::ATan(-(static_cast<Type>(2)*A) / (B + D));
        const Quat<Type> c(
            math::Cos(theta / static_cast<Type>(2)),
            math::Sin(theta / static_cast<Type>(2)),
            static_cast<Type>(0),
            static_cast<Type>(0));
            return Quat<Type>(
                p.S*c.S + p.X*c.X,
                p.X*c.S - p.S*c.X,
                p.Y*c.S - p.Z*c.X,
                p.Z*c.S + p.Y*c.X
            );
        }
        //////////////////////////////////////////////////////////
        #endif
        
        #if MV_EXPERIMENTAL
        //Experimental
        template <typename Type>
        static MV_API Quat<Type> FactorsRepeatedPrincipalAxis(const Quat<Type>& result)
        {
            const Quat<Type> p(GetUnitQuat(result));
            const Type A = p.X*(p.S + p.Z);
            const Type B = -p.S*p.X - p.S*p.Y;
            const Type D = p.Z*(p.X + p.Y);
            const Type E = -p.Y*p.Z - p.Y*p.S;
            
            const Type o = static_cast<Type>(1);
            const Type t = static_cast<Type>(2);
            
            const Type part1 = (B + D) / (t*E);
            const Type part2 = math::Sqrt(o - (A*t*t*E) / (math::Pow(B + D, t)));
            
            const Type res1 = part1 * (-o + part2);
            const Type res2 = part1 * (-o - part2);
            
            
            //a = resx * b
            //b = a / resx
            //sqrt(a^2 + (a/resx)^2) = 1 = sqrt(b^2 + (resx*b)^2)
            //-> 8 solutions for res1
            //-> 8 solutions for res2
            //Correct res?
            
            return result;
        }
        //////////////////////////////////////////////////////////
        #endif
        
        #if MV_EXPERIMENTAL
        //Experimental
        //Decomposes result into q1 and q2 which make up result (in standard reference frame)
        //Use empty pre-initialized quaternions for q1 and q2 for all data in them will be overwritten
        template <typename Type>
        static MV_API void FactorsThreePrincipalAxis(const Quat<Type>& result, Quat<Type>& q1, Quat<Type>& q2)
        {
            const Type at = math::ATan(result.Y / result.S);
            
            q1.S = math::Cos(at);
            q1.X = static_cast<Type>(0);
            q1.Y = math::Sin(at);
            q1.Z = static_cast<Type>(0);
            
            q2.S = q1.S*result.S + q1.Y*result.Y;
            q2.X = q1.S*result.X - q1.Y*result.Z;
            q2.Y = static_cast<Type>(0);
            q2.Z = q1.Y*result.X + q1.S*result.Z;
        }
        //////////////////////////////////////////////////////////
        #endif
        
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
        
        
        //Returns composite axis around which this rotation quaternion is happening
        template <typename Type>
        static MV_API Vektor<3u, Type> GetCompositeRotationAxis(const Quat<Type>& q)
        {
            return Vektor<3u, Type>(q._v);
        }
        //////////////////////////////////////////////////////////
        
        
        //Returns rotations from a quaternion
        //Pitch Yaw, Roll (XYZ)
        template <typename Type>
        static MV_API Vektor<3u, Type> GetRotations(const Quat<Type>& q)
        {
            const Type o = static_cast<Type>(1);
            const Type t = static_cast<Type>(2);
            return Vektor<3u, Type>(
                math::ATan2(t*(q.S*q.X + q.Y*q.Z), o - t*(q.X*q.X + q.Y*q.Y)),
                math::ASin(t*(q.S*q.Y - q.Z*q.X)),
                math::ATan2(t*(q.S*q.Z + q.X*q.Y), o - t*(q.Y*q.Y + q.Z*q.Z))
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
        //Combined rotation quaternion
        //Roll, Pitch, Roll (ZXZ)
        //Longitude, Latitude, Track direction
        template <typename Type>
        static MV_API Quat<Type> MakeCompositeOrbitEphemerisRotationQ(const Vektor<3u, Type>& rads)
        {
            using math::Cos;
            using math::Sin;
            
            const Type t = static_cast<Type>(2);
            //mu, epsilon, ro
            const Type cm = Cos(rads[0] / t);
            const Type sm = Sin(rads[0] / t);
            const Type ce = Cos(rads[1] / t);
            const Type se = Sin(rads[1] / t);
            const Type cr = Cos(rads[2] / t);
            const Type sr = Sin(rads[2] / t);
            
            return Quat<Type>(
                cm*ce*cr - sm*se*sr,
                cm*ce*sr + sm*se*cr,
                sm*ce*sr - cm*se*cr,
                cm*se*sr + sm*ce*cr
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
        //Combined rotation quaternion
        //Roll, Pitch, Roll (ZXZ)
        //to Ascending node, Orbit inclination, Argument of Latitude
        template <typename Type>
        static MV_API Quat<Type> MakeCompositeOrbitRotationQ(const Vektor<3u, Type>& rads)
        {
            using math::Cos;
            using math::Sin;
            
            //omega, beta, gamma
            const Type co = Cos(rads[0] * static_cast<Type>(0.5f));
            const Type so = Sin(rads[0] * static_cast<Type>(0.5f));
            const Type cb = Cos(rads[1] * static_cast<Type>(0.5f));
            const Type sb = Sin(rads[1] * static_cast<Type>(0.5f));
            const Type cg = Cos(rads[2] * static_cast<Type>(0.5f));
            const Type sg = Sin(rads[2] * static_cast<Type>(0.5f));
            
            return Quat<Type>(
                co*cb*cg - so*cb*sg,
                co*sb*cg + so*sb*sg,
                -co*sb*sg + so*sb*cg,
                co*cb*sg + so*cb*cg
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
            return Multiply(Multiply(
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
                )),
                Mat<4u, 4u, Type>(
                    Cos(rads[2]), -Sin(rads[2]), z, z,
                    Sin(rads[2]), Cos(rads[2]), z, z,
                    z, z, o, z,
                    z, z, z, o
                ));
            }
            //////////////////////////////////////////////////////////
            #endif
            
            //Combined rotation quaternion
            //Pitch, Yaw, Roll (XYZ)
            template <typename Type>
            static MV_API Quat<Type> MakeCompositeRotationQ(const Vektor<3u, Type>& rads)
            {
                using math::Cos;
                using math::Sin;
                
                //pitch, yaw, roll
                const Type cp = Cos(rads[0] * static_cast<Type>(0.5f));
                const Type sp = Sin(rads[0] * static_cast<Type>(0.5f));
                const Type cy = Cos(rads[1] * static_cast<Type>(0.5f));
                const Type sy = Sin(rads[1] * static_cast<Type>(0.5f));
                const Type cr = Cos(rads[2] * static_cast<Type>(0.5f));
                const Type sr = Sin(rads[2] * static_cast<Type>(0.5f));
                
                return Quat<Type>(
                    cr*cy*cp + sr*sy*sp,
                    cr*cy*sp - sr*sy*cp,
                    cr*sy*cp + sr*cy*sp,
                    sr*cy*cp - cr*sy*sp
                );
            }
            //////////////////////////////////////////////////////////
            
            
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
            
            
            //Returns rotation quaternion rad radians around axis v
            template <typename Type>
            static MV_API Quat<Type> MakeRotationQ(const Vektor<3u, Type>& v, const Type rad)
            {
                return Quat<Type>(
                    math::Cos(rad / static_cast<Type>(2)),
                    mv::Multiply(mv::GetUnitVektor(v), math::Sin(rad / static_cast<Type>(2)))
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
            
            
            //Rotates v around a rotation of q
            template <typename Type>
            static MV_API Vektor<3u, Type> Multiply(const Quat<Type>& q, const Vektor<3u, Type>& v)
            {
                return Multiply(Multiply(q, Quat<Type>(v)), Inverse(q))._v;
            }
            //////////////////////////////////////////////////////////
            
            
            //Multiply matrix and vektor
            template <UInt8 Rows, UInt8 Size, typename Type>
            static MV_API Mat<Rows, 1u, Type> Multiply(const Mat<Rows, Size, Type>& m, const Vektor<Size, Type>& v)
            {
                return Multiply(m, ToMatV(v));
            }
            //////////////////////////////////////////////////////////
            
            
            //Multiply vektor and matrix
            template <UInt8 Cols, UInt8 Size, typename Type>
            static MV_API Mat<1u, Cols, Type> Multiply(const Vektor<Size, Type>& v, const Mat<Size, Cols, Type>& m)
            {
                return Multiply(ToMatH(v), m);
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
            
            
            #if MV_DEBUG
            //Print quaternion
            template <typename Type>
            static MV_API void Print(const Quat<Type>& q)
            {
                _msg(q.S, ',', q.X, ',', q.Y, ',', q.Z, '\n');
            }
            //////////////////////////////////////////////////////////
            #endif
            
            
            #if MV_DEBUG
            //Print vektor
            template <UInt8 Size, typename Type>
            static MV_API void Print(const Vektor<Size, Type>& v)
            {
                for(const auto& itr : v.data())
                _msg(itr, ',');
                _msg('\n');
            }
            //////////////////////////////////////////////////////////
            #endif
            
            #if MV_EXPERIMENTAL
            //Experimental
            template <typename Type>
            static MV_API Quat<Type> QuaternionPerturbation(Quat<Type>& q)
            {
                if (q.S == static_cast<Type>(0))
                return q;
                
                MV_ASSERT(false, "Math under construction");
                return q;
            }
            //////////////////////////////////////////////////////////
            #endif
            
            
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
            template <typename Type>
            static MV_API Vektor<3u, Type> RotateAroundAxis(const Quat<Type>& rot, const Vektor<3u, Type>& point)
            {
                return Multiply(rot, point);
            }
            template <typename Type>
            static MV_API Vektor<3u, Type> RotateAroundAxis(const Vektor<3u, Type>& axis, const Vektor<3u, Type>& point, const Type rad)
            {
                return Multiply(
                    Quat<Type>(math::Cos(rad / static_cast<Type>(2)),
                    Multiply(GetUnitVektor(axis), math::Sin(rad / static_cast<Type>(2)))),
                    point);
                }
                //////////////////////////////////////////////////////////
                
                
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
                
                
                template <typename Type>
                static MV_API Type ToDeg(const Type rad)
                {
                    return static_cast<Type>(static_cast<MV_TYPE>(rad)*static_cast<MV_TYPE>(180.0) / static_cast<MV_TYPE>(math::MV_PI));
                }
                template <typename Type>
                static MV_API Type ToRad(const Type deg)
                {
                    return static_cast<Type>(static_cast<MV_TYPE>(deg)*static_cast<MV_TYPE>(math::MV_PI) / static_cast<MV_TYPE>(180.0));
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
                
                
                //Returns vektor of target size (TS)
                template <UInt8 TS, typename Type, UInt8 VS>
                static MV_API Vektor<TS, Type> ToVek(const Vektor<VS, Type>& v, const Vektor<TS-VS, Type>* fill = nullptr)
                {
                    std::array<Type, TS> arr;
                    if MV_API (TS <= VS)
                    for(UInt8 i = 0u; i < TS; ++i)
                    arr[i] = v[i];
                    else
                    {
                        for (UInt8 i = 0u; i < VS; ++i)
                        arr[i] = v[i];
                        if (fill)
                        for(UInt8 i = 0u; i < (TS-VS); ++i)
                        arr[i + VS] = (*fill)[i];
                        else
                        for(UInt8 i = VS; i < TS; ++i)
                        arr[i] = static_cast<Type>(0);
                    }
                    return Vektor<TS, Type>(arr);
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
                
                
                template <UInt8 TS, typename Type, UInt8 VS>
                static MV_API Vektor<TS, Type> ToVek(const Vektor<VS, Type>& v, const Type arg)
                {
                    const Vektor<1u, Type> vt(arg);
                    return ToVek<TS>(v, &vt);
                }
                //////////////////////////////////////////////////////////
                
                
                template <UInt8 TS, typename Type, UInt8 VS, typename ... Args>
                static MV_API Vektor<TS, Type> ToVek(const Vektor<VS, Type>& v, Args...args)
                {
                    static_assert(sizeof...(Args) == (TS-VS), "ToVek: Invalid variadic pack size! Should be TS-VS");
                    //const std::array<Type, sizeof...(Args)> arr = { std::forward<Args>(args)... };
                    //const Vektor<sizeof...(Args)> vt(arr);
                    const Vektor<sizeof...(Args), Type> vt(std::forward<Args>(args)...);
                    return ToVek<TS>(v, &vt);
                }
                //////////////////////////////////////////////////////////
                
                
                template <typename Type>
                static MV_API Vektor<4u, Type> TransformDirection(const Mat<4u, 4u, Type>& m, const Vektor<4u, Type>& v)
                {
                    return ToVek(Multiply(m, v));
                }
                //////////////////////////////////////////////////////////
                
                
                
                //    //Translate matrix 4
                //    template <typename Type>
                //    static Mat<4u, 4u, Type> Translate(const Mat<4u, 4u, Type>& m, const Vektor<3u, Type>& v)
                //    {
                    //        return Multiply(m, MakeTranslation(v));
                    //    }
                    //    //////////////////////////////////////////////////////////
                    
                    
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
                    
                    
                    
                    template <UInt8 NewSize, typename Type, UInt8 OldSize>
                    static MV_API Vektor<NewSize, Type> ResizeVektor(const Vektor<OldSize, Type>& v)
                    {
                        if MV_API (NewSize == OldSize)
                        {
                            return v;
                        }
                        else if MV_API (NewSize < OldSize)
                        {
                            std::array<Type, NewSize> arr;
                            for (UInt8 i = 0u; i < NewSize; ++i)
                            {
                                arr[i] = v[i];
                            }
                            return Vektor<NewSize, Type>(arr);
                        }
                        else
                        {
                            std::array<Type, NewSize> arr;
                            for (UInt8 i = 0u; i < OldSize; ++i)
                            {
                                arr[i] = v[i];
                            }
                            for (UInt8 i = OldSize; i < NewSize; ++i)
                            {
                                arr[i] = static_cast<Type>(0);
                            }
                            return Vektor<NewSize, Type>(arr);
                        }
                    }
                    
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
                    
                }
                
                #pragma warning(pop) // 4201
                
                #undef MV_API
                #undef MV_ASSERT
                #undef MV_DISABLE_CLASS
                #undef MV_TYPE