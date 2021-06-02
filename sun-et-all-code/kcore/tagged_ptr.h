#ifndef TAGGED_PTR_HPP
#define TAGGED_PTR_HPP

#include <cstddef>
#include <limits>

template <typename T>
class tagged_ptr
{
public:
  using tag_t = std::size_t;

  tagged_ptr() noexcept: ptr(nullptr), tag(0) {}
  explicit tagged_ptr(T* p, tag_t t = 0) : ptr(p), tag(t) {}

  tag_t next_tag() const
  {
    return (tag + 1) & (std::numeric_limits<tag_t>::max)();
  }

  void set_tag(tag_t t)
  {
    tag = t;
  }

  T* get_ptr() const
  {
    return ptr;
  }

  void set_ptr(T* p)
  {
    ptr = p;
  }

  T& operator*() const
  {
    return *ptr;
  }

  T* operator->() const
  {
    return ptr;
  }

  operator bool(void) const
  {
    return ptr != nullptr;
  }

  bool operator==(volatile tagged_ptr<T>& other)
  {
    return (ptr == other.ptr) && (tag == other.tag);
  }

  bool operator!= (volatile tagged_ptr<T>& other)
  {
    return !(operator==(other));
  }

private:
  T* ptr;
  tag_t tag;
};

#endif /* TAGGED_PTR_HPP */