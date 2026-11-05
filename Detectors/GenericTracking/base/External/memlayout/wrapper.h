// Copyright 2019-2020 CERN and copyright holders of ALICE O2.
// See https://alice-o2.web.cern.ch/copyright for details of the copyright holders.
// All rights not expressly granted are reserved.
//
// This software is distributed under the terms of the GNU General Public
// License v3 (GPL Version 3), copied verbatim in the file "COPYING".
//
// In applying this license CERN does not waive the privileges and immunities
// granted to it by virtue of its status as an Intergovernmental Organization
// or submit itself to any jurisdiction.

#ifndef MEMLAYOUT_H
#define MEMLAYOUT_H

// Upstream forward-declares std::iterator_traits here.  libc++ declares it
// differently, so the two are ambiguous in any translation unit that also
// includes a standard container.  Including the real header instead.
#include <iterator>
#include <type_traits>

#ifndef MEMLAYOUT_HOST_DEVICE
#ifdef GPUhdi
#define MEMLAYOUT_HOST_DEVICE GPUhdi()
#else
#define MEMLAYOUT_HOST_DEVICE
#endif
#endif

// Local addition: the element access -- operator[] and what it calls -- is
// always inlined. It builds a proxy of one reference per column; inlined,
// the references a caller never reads fold away and `view[i].x` is one load,
// but a compiler that keeps it out of line (GCC at -O2 does, for a record of
// many columns) computes and stores every one of them on each call. nvcc and
// HIP's clang take the attribute as they take __forceinline__.
//
// It only reaches as far as the proxy does inline: a proxy returned by, or
// passed to, a function the compiler keeps out of line is built whole again.
// Code that needs a few columns of a row is better served reading them from
// the view, `view.x[i]`, which builds no proxy at all.
#ifndef MEMLAYOUT_INLINE
#if defined(__GNUC__) || defined(__clang__)
#define MEMLAYOUT_INLINE __attribute__((always_inline))
#else
#define MEMLAYOUT_INLINE
#endif
#endif

namespace memlayout
{

using size_t = decltype(sizeof(0));
using ptrdiff_t = decltype(static_cast<int*>(nullptr) - static_cast<int*>(nullptr));

template <class T>
using value = T;
template <class T>
using reference = T&;
template <class T>
using const_reference = const T&;
template <class T>
using pointer = T*;
template <class T>
using const_pointer = const T*;

// Local addition: a default value for a field of a struct templated on its
// storage.  In the record form, Container = value, the field takes `v`; in
// every other form it is value-initialised -- a null pointer, an empty
// container -- so one declaration serves as record and as columns:
//
//   C<float> x = memlayout::default_value<C>(0.F);
//
// The reference forms have no default; that is safe because a default member
// initializer of a class template is instantiated only where it is used, and
// the wrapper always initialises every member.
template <template <class> class Container, class T>
MEMLAYOUT_HOST_DEVICE constexpr Container<T> default_value(T v)
{
  if constexpr (std::is_same_v<Container<T>, T>) {
    return v;
  } else {
    return Container<T>{};
  }
}

struct Dummy {
};

// The concepts are stated over a type, Container<Dummy>, and the SoA
// specialisations are constrained through those. nvcc (13.1, 13.2) mishandles
// concepts over a template template parameter used in a partial
// specialisation's constraint: it picks the mutable random-access
// specialisation for const_pointer -- or its front end, cudafe++, aborts on a
// larger translation unit. The type forms mean the same and keep the
// subsumption the specialisations are ordered by: the const one's constraint
// contains the random-access one.
template <class It>
concept random_access_like_type = requires(It it, size_t n) { it + n; it - n; it - it; it[n]; ++it; --it; *it; };

template <class It>
concept writable_iterator_type = requires(It it) { *it = Dummy{}; };

template <class It>
concept const_random_access_like_type = random_access_like_type<It> && !writable_iterator_type<It>;

template <template <class> class Container>
concept random_access_like = random_access_like_type<Container<Dummy>>;

template <template <class> class Iterator>
concept non_const_iterator_like = writable_iterator_type<Iterator<Dummy>>;

template <class ReturnType>
struct RandomAccessAt {
  memlayout::size_t i;
  template <class... Args>
  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(Args&... args) const
  {
    return {args[i]...};
  }
  template <class... Args>
  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(const Args&... args) const
  {
    return {args[i]...};
  }
};

template <class ReturnType>
struct GetPointer {
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(Args&... args) const
  {
    return {&args...};
  }
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(const Args&... args) const
  {
    return {&args...};
  }
};

template <class ReturnType>
struct AggregateConstructor {
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(Args&... args) const
  {
    return {args...};
  }
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(const Args&... args) const
  {
    return {args...};
  }
};

struct FirstMember {
  template <class T, class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr T& operator()(T& t, Args&... args) const
  {
    return t;
  }
  template <class T, class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr const T& operator()(const T& t, const Args&... args) const
  {
    return t;
  }
};

template <class ReturnType>
struct PreIncrement {
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(Args&... args) const
  {
    return {++args...};
  }
};

template <class ReturnType>
struct PreDecrement {
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(Args&... args) const
  {
    return {--args...};
  }
};

template <class ReturnType>
struct Advance {
  ptrdiff_t i;
  template <class... Args>
  MEMLAYOUT_HOST_DEVICE constexpr ReturnType operator()(const Args&... args) const
  {
    return {(args + i)...};
  }
};

struct CopyAssignment {
  template <class Left, class Right>
  MEMLAYOUT_HOST_DEVICE constexpr Left& operator()(Left& left, const Right& right) const
  {
    return left = right;
  }
};

template <class ReturnType>
struct Begin {
  template <class... Args>
  constexpr ReturnType operator()(Args&... args) const
  {
    return {args.begin()...};
  }
  template <class... Args>
  constexpr ReturnType operator()(const Args&... args) const
  {
    return {args.cbegin()...};
  }
};

template <class ReturnType>
struct End {
  template <class... Args>
  constexpr ReturnType operator()(Args&... args) const
  {
    return {args.end()...};
  }
  template <class... Args>
  constexpr ReturnType operator()(const Args&... args) const
  {
    return {args.cend()...};
  }
};

enum class Layout { aos = 0,
                    soa = 1 };

template <
  template <template <class> class> class Struct,
  template <class> class Container,
  Layout L = Layout::soa>
struct Wrapper;

template <
  template <template <class> class> class Struct,
  template <class> class Container>
struct Wrapper<Struct, Container, Layout::aos> : public Container<Struct<value>> {
  static constexpr Layout layout_type = Layout::aos;
  using Base = Container<Struct<value>>;
  using Base::Base;

  constexpr Wrapper() = default;
  template <template <class> class other_Container>
  constexpr Wrapper(Wrapper<Struct, other_Container, Layout::aos>& other) : Base(other)
  {
  }

  MEMLAYOUT_INLINE constexpr Wrapper<Struct, reference> operator[](size_t i) { return Base::operator[](i); }
  MEMLAYOUT_INLINE constexpr Wrapper<Struct, const_reference> operator[](size_t i) const { return Base::operator[](i); }

  constexpr Wrapper<Struct, reference> operator*() { return operator[](0); }
  constexpr Wrapper<Struct, const_reference> operator*(ptrdiff_t) const { return operator[](0); }
};

template <template <template <class> class> class Struct>
struct Wrapper<Struct, pointer, Layout::aos> {
  static constexpr Layout layout_type = Layout::aos;

  pointer<Struct<value>> data;
  using Data = pointer<Struct<value>>;

  operator Data&() { return data; }
  operator const Data&() const { return data; }

  constexpr bool operator==(const Wrapper& other) const { return data == other.data; }
  constexpr bool operator!=(const Wrapper& other) const { return data != other.data; }
  constexpr bool operator<(const Wrapper& other) const { return data < other.data; }

  constexpr Wrapper operator+(ptrdiff_t i) const { return {data + i}; }
  constexpr Wrapper operator-(ptrdiff_t i) const { return {data - i}; }
  constexpr ptrdiff_t operator-(const Wrapper& other) const { return {data - other.data}; }

  constexpr Wrapper operator++() { return {++data}; }
  constexpr Wrapper operator+=(ptrdiff_t i) { return {data += i}; }
  constexpr Wrapper operator--() { return {--data}; }
  constexpr Wrapper operator-=(ptrdiff_t i) { return {data -= i}; }

  MEMLAYOUT_INLINE constexpr Wrapper<Struct, reference> operator[](size_t i) { return data[i]; }
  MEMLAYOUT_INLINE constexpr Wrapper<Struct, const_reference> operator[](size_t i) const { return data[i]; }

  constexpr Wrapper<Struct, reference> operator*() { return *data; }
  constexpr Wrapper<Struct, const_reference> operator*(ptrdiff_t) const { return *data; }

  constexpr Wrapper<Struct, reference> operator->() { return operator[](0); }
  constexpr Wrapper<Struct, const_reference> operator->() const { return operator[](0); }
};

template <template <template <class> class> class Struct>
struct Wrapper<Struct, const_pointer, Layout::aos> {
  static constexpr Layout layout_type = Layout::aos;

  const_pointer<Struct<value>> data;
  using Data = const_pointer<Struct<value>>;

  operator const Data&() const { return data; }

  constexpr bool operator==(const Wrapper& other) const { return data == other.data; }
  constexpr bool operator!=(const Wrapper& other) const { return data != other.data; }
  constexpr bool operator<(const Wrapper& other) const { return data < other.data; }

  constexpr Wrapper operator+(ptrdiff_t i) const { return {data + i}; }
  constexpr Wrapper operator-(ptrdiff_t i) const { return {data - i}; }
  constexpr ptrdiff_t operator-(const Wrapper& other) const { return {data - other.data}; }

  constexpr Wrapper operator++() { return {++data}; }
  constexpr Wrapper operator+=(ptrdiff_t i) { return {data += i}; }
  constexpr Wrapper operator--() { return {--data}; }
  constexpr Wrapper operator-=(ptrdiff_t i) { return {data -= i}; }

  MEMLAYOUT_INLINE constexpr Wrapper<Struct, const_reference> operator[](size_t i) const { return data[i]; }
  constexpr Wrapper<Struct, const_reference> operator*() const { return *data; }
  constexpr Wrapper<Struct, const_reference> operator->() const { return operator[](0); }
};

template <
  template <template <class> class> class Struct,
  template <class> class Container>
struct Wrapper<Struct, Container, Layout::soa> : public Struct<Container> {
  static constexpr Layout layout_type = Layout::soa;
  using Base = Struct<Container>;
  using Base::Base;
  template <class T>
  using iterator = typename Container<T>::iterator;
  template <class T>
  using const_iterator = typename Container<T>::const_iterator;

  constexpr Wrapper() = default;
  constexpr Wrapper(Base b) : Base{static_cast<Base&&>(b)} {}
  template <template <class> class other_Container>
  constexpr Wrapper(Struct<other_Container>& other) : Base{other.apply(AggregateConstructor<Base>{})}
  {
  }
  template <template <class> class other_Container>
  constexpr Wrapper(const Struct<other_Container>& other) : Base{other.apply(AggregateConstructor<Base>{})}
  {
  }

  MEMLAYOUT_INLINE constexpr Wrapper<Struct, reference> operator[](size_t i) { return Base::apply(RandomAccessAt<Struct<reference>>{i}); }
  MEMLAYOUT_INLINE constexpr Wrapper<Struct, const_reference> operator[](size_t i) const { return Base::apply(RandomAccessAt<Struct<const_reference>>{i}); }

  constexpr Wrapper<Struct, reference> operator*() { return operator[](0); }
  constexpr Wrapper<Struct, const_reference> operator*(ptrdiff_t) const { return operator[](0); }

  constexpr Wrapper<Struct, iterator> begin() { return Base::apply(Begin<Struct<iterator>>{}); }
  constexpr Wrapper<Struct, iterator> end() { return Base::apply(End<Struct<iterator>>{}); }
  constexpr Wrapper<Struct, const_iterator> cbegin() const { return Base::apply(Begin<Struct<const_iterator>>{}); }
  constexpr Wrapper<Struct, const_iterator> cend() const { return Base::apply(End<Struct<const_iterator>>{}); }
};

template <
  template <template <class> class> class Struct,
  template <class> class Container>
  requires random_access_like_type<Container<Dummy>>
struct Wrapper<Struct, Container, Layout::soa> : public Struct<Container> {
  using iterator_category = std::iterator_traits<Container<Dummy>>::iterator_category;
  using difference_type = ptrdiff_t;
  using value_type = Wrapper<Struct, memlayout::value>;
  using pointer = void;
  using reference = Wrapper<Struct, memlayout::reference>;

  static constexpr Layout layout_type = Layout::soa;
  using Base = Struct<Container>;

  constexpr Wrapper() = default;
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper(Base b) : Base{static_cast<Base&&>(b)} {}

  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, memlayout::reference> operator[](size_t i) { return Base::apply(RandomAccessAt<Struct<memlayout::reference>>{i}); }
  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr const Wrapper<Struct, memlayout::const_reference> operator[](size_t i) const { return Base::apply(RandomAccessAt<Struct<memlayout::const_reference>>{i}); }

  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, memlayout::reference> operator*() { return operator[](0); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, memlayout::const_reference> operator*() const { return operator[](0); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, memlayout::reference> operator->() { return operator[](0); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, memlayout::const_reference> operator->() const { return operator[](0); }

  MEMLAYOUT_HOST_DEVICE constexpr bool operator==(const Wrapper& other) const { return Base::apply(FirstMember{}) == other.apply(FirstMember{}); }
  MEMLAYOUT_HOST_DEVICE constexpr bool operator!=(const Wrapper& other) const { return !this->operator==(other); }
  MEMLAYOUT_HOST_DEVICE constexpr bool operator<(const Wrapper& other) const { return Base::apply(FirstMember{}) < other.apply(FirstMember{}); }

  MEMLAYOUT_HOST_DEVICE constexpr Wrapper operator+(ptrdiff_t i) const { return Base::apply(Advance<Base>{i}); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper operator-(ptrdiff_t i) const { return operator+(-i); }
  MEMLAYOUT_HOST_DEVICE constexpr ptrdiff_t operator-(const Wrapper& other) const { return Base::apply(FirstMember{}) - other.apply(FirstMember{}); }

  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator++()
  {
    Base::apply(PreIncrement<Base>{});
    return *this;
  }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator+=(ptrdiff_t i) { return *this = *this + i; }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator--()
  {
    Base::apply(PreDecrement<Base>{});
    return *this;
  }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator-=(ptrdiff_t i) { return *this = *this - i; }
};

template <
  template <template <class> class> class Struct,
  template <class> class Container>
  requires const_random_access_like_type<Container<Dummy>>
struct Wrapper<Struct, Container, Layout::soa> : public Struct<Container> {
  using iterator_category = std::iterator_traits<Container<Dummy>>::iterator_category;
  using difference_type = ptrdiff_t;
  using value_type = Wrapper<Struct, memlayout::value>;
  using pointer = void;
  using reference = Wrapper<Struct, memlayout::const_reference>;

  static constexpr Layout layout_type = Layout::soa;
  using Base = Struct<Container>;

  constexpr Wrapper() = default;
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper(Base b) : Base{static_cast<Base&&>(b)} {}

  template <template <class> class OtherContainer>
    requires random_access_like_type<OtherContainer<Dummy>>
  constexpr Wrapper(const Struct<OtherContainer>& other) : Base(other.apply(AggregateConstructor<Base>{}))
  {
  }

  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, const_reference> operator[](size_t i) const { return Base::apply(RandomAccessAt<Struct<const_reference>>{i}); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, const_reference> operator*() const { return operator[](0); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, const_reference> operator->() const { return operator[](0); }

  MEMLAYOUT_HOST_DEVICE constexpr bool operator==(const Wrapper& other) const { return Base::apply(FirstMember{}) == other.apply(FirstMember{}); }
  MEMLAYOUT_HOST_DEVICE constexpr bool operator!=(const Wrapper& other) const { return !this->operator==(other); }
  MEMLAYOUT_HOST_DEVICE constexpr bool operator<(const Wrapper& other) const { return Base::apply(FirstMember{}) < other.apply(FirstMember{}); }

  MEMLAYOUT_HOST_DEVICE constexpr Wrapper operator+(ptrdiff_t i) const { return Base::apply(Advance<Base>{i}); }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper operator-(ptrdiff_t i) const { return operator+(-i); }
  MEMLAYOUT_HOST_DEVICE constexpr ptrdiff_t operator-(const Wrapper& other) const { return Base::apply(FirstMember{}) - other.apply(FirstMember{}); }

  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator++()
  {
    Base::apply(PreIncrement<Base>{});
    return *this;
  }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator+=(ptrdiff_t i) { return *this = *this + i; }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator--()
  {
    Base::apply(PreDecrement<Base>{});
    return *this;
  }
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper& operator-=(ptrdiff_t i) { return *this = *this - i; }
};

template <template <template <class> class> class Struct>
struct Wrapper<Struct, value, Layout::soa> : public Struct<value> {
  using Base = Struct<value>;

  constexpr Wrapper() = default;
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper(Base b) : Base{static_cast<Base&&>(b)} {}
  constexpr Wrapper(const Struct<reference>& other) : Base(other.apply(AggregateConstructor<Base>{})) {}
  constexpr Wrapper(const Struct<const_reference>& other) : Base(other.apply(AggregateConstructor<Base>{})) {}
};

template <template <template <class> class> class Struct>
struct Wrapper<Struct, reference, Layout::soa> : public Struct<reference> {
  using Base = Struct<reference>;

  constexpr Wrapper() = delete;
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper(Base b) : Base{static_cast<Base&&>(b)} {}
  constexpr Wrapper(Struct<value>& other) : Base(other.apply(AggregateConstructor<Base>{})) {}

  constexpr Wrapper(const Wrapper& other) = default;

  constexpr Wrapper& operator=(const Wrapper<Struct, value>& other)
  {
    Base::apply(other, CopyAssignment{});
    return *this;
  }
  constexpr Wrapper& operator=(const Wrapper& other)
  {
    Base::apply(other, CopyAssignment{});
    return *this;
  }
  constexpr Wrapper& operator=(const Wrapper<Struct, const_reference>& other)
  {
    Base::apply(other, CopyAssignment{});
    return *this;
  }

  constexpr Wrapper(Wrapper&& other) = default;

  constexpr Wrapper& operator=(Wrapper&& other) { return operator=(other); }

  constexpr Wrapper<Struct, pointer> operator&() { return Base::apply(GetPointer<Struct<pointer>>{}); }
  // constexpr Wrapper<Struct, const_pointer> operator&() const { return Base::apply(GetPointer<Struct<const_pointer>>{}); }
  constexpr pointer<Wrapper<Struct, reference>> operator->() { return this; }
};

template <template <template <class> class> class Struct>
struct Wrapper<Struct, const_reference, Layout::soa> : public Struct<const_reference> {
  using Base = Struct<const_reference>;

  constexpr Wrapper() = delete;
  MEMLAYOUT_HOST_DEVICE constexpr Wrapper(Base b) : Base{static_cast<Base&&>(b)} {}
  constexpr Wrapper(const Struct<value>& other) : Base(other.apply(AggregateConstructor<Base>{})) {}
  constexpr Wrapper(const Struct<reference>& other) : Base(other.apply(AggregateConstructor<Base>{})) {}

  MEMLAYOUT_HOST_DEVICE constexpr Wrapper<Struct, const_pointer> operator&() const { return Base::apply(GetPointer<Struct<const_pointer>>{}); }
  MEMLAYOUT_HOST_DEVICE constexpr const_pointer<Wrapper<Struct, const_reference>> operator->() const { return this; }
};

} // namespace memlayout

#define MEMLAYOUT_PARENS ()

#define MEMLAYOUT_EXPAND(...) MEMLAYOUT_EXPAND4(MEMLAYOUT_EXPAND4(MEMLAYOUT_EXPAND4(__VA_ARGS__)))
#define MEMLAYOUT_EXPAND4(...) MEMLAYOUT_EXPAND3(MEMLAYOUT_EXPAND3(MEMLAYOUT_EXPAND3(__VA_ARGS__)))
#define MEMLAYOUT_EXPAND3(...) MEMLAYOUT_EXPAND2(MEMLAYOUT_EXPAND2(MEMLAYOUT_EXPAND2(__VA_ARGS__)))
#define MEMLAYOUT_EXPAND2(...) MEMLAYOUT_EXPAND1(MEMLAYOUT_EXPAND1(MEMLAYOUT_EXPAND1(__VA_ARGS__)))
#define MEMLAYOUT_EXPAND1(...) __VA_ARGS__

#define MEMLAYOUT_FOR_EACH(macro, ...) \
  __VA_OPT__(MEMLAYOUT_EXPAND(MEMLAYOUT_FOR_EACH_HELPER(macro, __VA_ARGS__)))
#define MEMLAYOUT_FOR_EACH_HELPER(macro, a1, ...) \
  macro(a1)                                       \
    __VA_OPT__(, MEMLAYOUT_FOR_EACH_AGAIN MEMLAYOUT_PARENS(macro, __VA_ARGS__))
#define MEMLAYOUT_FOR_EACH_AGAIN() MEMLAYOUT_FOR_EACH_HELPER

// per-member transforms: bare name for the unary case, `f(m, other.m)` for binary
#define MEMLAYOUT_UNARY_ONE(m) m
#define MEMLAYOUT_BINARY_ONE(m) f(m, other.m)

// One macro per struct, member list written exactly once. Generates a column
// enum in the same order as the apply arguments. The final enumerator gives
// the number of columns; member names must not use ColumnCount.
// Also generates the four `apply(other, f)` overloads and the two unary
// `apply(f)` overloads.
#define MEMLAYOUT_APPLY(STRUCT_NAME, ...)                                                    \
  enum class Column : memlayout::size_t { __VA_ARGS__,                                       \
                                          ColumnCount };                                     \
  template <class Function>                                                                  \
  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr auto apply(Function&& f)                  \
  {                                                                                          \
    return f(MEMLAYOUT_FOR_EACH(MEMLAYOUT_UNARY_ONE, __VA_ARGS__));                          \
  }                                                                                          \
  template <class Function>                                                                  \
  MEMLAYOUT_INLINE MEMLAYOUT_HOST_DEVICE constexpr auto apply(Function&& f) const            \
  {                                                                                          \
    return f(MEMLAYOUT_FOR_EACH(MEMLAYOUT_UNARY_ONE, __VA_ARGS__));                          \
  }                                                                                          \
  template <template <class> class other_Container, class Function>                          \
  constexpr STRUCT_NAME apply(STRUCT_NAME<other_Container>& other, Function&& f)             \
  {                                                                                          \
    return {MEMLAYOUT_FOR_EACH(MEMLAYOUT_BINARY_ONE, __VA_ARGS__)};                          \
  }                                                                                          \
  template <template <class> class other_Container, class Function>                          \
  constexpr STRUCT_NAME apply(STRUCT_NAME<other_Container>& other, Function&& f) const       \
  {                                                                                          \
    return {MEMLAYOUT_FOR_EACH(MEMLAYOUT_BINARY_ONE, __VA_ARGS__)};                          \
  }                                                                                          \
  template <template <class> class other_Container, class Function>                          \
  constexpr STRUCT_NAME apply(const STRUCT_NAME<other_Container>& other, Function&& f)       \
  {                                                                                          \
    return {MEMLAYOUT_FOR_EACH(MEMLAYOUT_BINARY_ONE, __VA_ARGS__)};                          \
  }                                                                                          \
  template <template <class> class other_Container, class Function>                          \
  constexpr STRUCT_NAME apply(const STRUCT_NAME<other_Container>& other, Function&& f) const \
  {                                                                                          \
    return {MEMLAYOUT_FOR_EACH(MEMLAYOUT_BINARY_ONE, __VA_ARGS__)};                          \
  }

#endif // MEMLAYOUT_H
