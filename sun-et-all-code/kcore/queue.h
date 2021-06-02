#ifndef LOCKFREE_QUEUE_H
#define LOCKFREE_QUEUE_H

#include <atomic>   /* std::atomic */

#include "tagged_ptr.h"

/*
 * lqueue provides a multi-reader/multi-writer queue with lock-free push and pop
 * operations. Construction and destruction must be synchronized. Based off of
 * the lock-free algorithm from "Simple, Fast, and Practical Non-Blocking and
 * Blocking Concurrent Queue Algorithms" by MM Michael and ML Scott (1996).
 */
template<typename T>
class lqueue {
public:

  /*
   * Constructs an empty queue.
   */
  lqueue() {
    node *dummy = new node();
    head = tagged_ptr<node>(dummy, 0);
    tail = tagged_ptr<node>(dummy, 0);
  }

  /*
   * Pushes object t into the queue. Returns true if the push succeeds.
   */
  bool push(const T &t) {
    node *n = new node(t);
    for (;;) {
      auto old_tail = tail.load(std::memory_order_acquire);
      auto old_next = old_tail->next.load(std::memory_order_acquire);
      if (old_tail==tail.load(std::memory_order_acquire)) {
        if (!old_next) {
          tagged_ptr<node> new_next(n, old_next.next_tag());
          if (old_tail->next.compare_exchange_weak(old_next, new_next)) {
            tagged_ptr<node> new_tail(n, old_tail.next_tag());
            tail.compare_exchange_strong(old_tail, new_tail);
            return true;
          }
        } else {
          tagged_ptr<node> new_tail(old_next.get_ptr(), old_tail.next_tag());
          tail.compare_exchange_strong(old_tail, new_tail);
        }
      }
    }
  }

  /*
   * Pops an object from the queue and stores it in ret. Returns true if the
   * pop succeeds.
   */
  bool pop(T &ret) {
    for (;;) {
      auto old_head = head.load(std::memory_order_acquire);
      auto old_tail = tail.load(std::memory_order_acquire);
      auto old_next = old_head->next.load(std::memory_order_acquire);
      if (old_head==head.load(std::memory_order_acquire)) {
        if (old_head.get_ptr()==old_tail.get_ptr()) {
          if (!old_next) {
            return false;
          }
          tagged_ptr<node> new_tail(old_next.get_ptr(), old_tail.next_tag());
          tail.compare_exchange_strong(old_tail, new_tail);
        } else {
          ret = old_next->data;
          tagged_ptr<node> new_head(old_next.get_ptr(), old_head.next_tag());
          if (head.compare_exchange_strong(old_head, new_head)) {
            //delete old_head.get_ptr();
            return true;
          }
        }
      }
    }
  }

private:
  class node {
  public:
    node() noexcept {
      next = tagged_ptr<node>(nullptr);
    }
    node(const T &v) : data(v) {
      auto old_next = next.load(std::memory_order_relaxed);
      tagged_ptr<node> new_next(nullptr, old_next.next_tag());
      next.store(new_next, std::memory_order_release);
    }

    std::atomic<tagged_ptr<node>> next;
    T data;
  };

  std::atomic<tagged_ptr<node>> head, tail;
};

#endif /* LOCKFREE_QUEUE_H */