#ifndef ZISA_TREE_HPP
#define ZISA_TREE_HPP

#include <algorithm>
#include <memory>

#include <zisa/config.hpp>

namespace zisa {

template <class T>
struct integer_tree_traits {
  static constexpr T sentinel = std::numeric_limits<T>::max();
};

template <class Int>
struct tree_traits
    : public integer_tree_traits<
          typename std::enable_if<std::is_integral<Int>::value, Int>::type> {};

template <class T, int WIDTH>
class Tree {
private:
  class Layer;

  class Node {
  public:
    T value = tree_traits<T>::sentinel;

  public:
    Node() = default;
    explicit Node(const T &value)
        : value(value), next_layer(std::make_unique<Layer>()){};

    template <class F>
    void insert(const T &value, const F &f) {
      if (is_leaf()) {
        if (is_empty()) {
          this->value = value;
        } else {
          next_layer = std::make_unique<Layer>();
          next_layer->insert(value, f);
        }
      } else {
        next_layer->insert(value, f);
      }
    }

    template <class F>
    const T &locate(const F &f) const {
      LOG_ERR_IF(is_empty(), "Can't locate anything in an empty node.");
      return is_leaf() ? this->value : next_layer->locate(f);
    }

    bool is_leaf() const { return next_layer == nullptr; }
    bool is_empty() const { return value == tree_traits<T>::sentinel; }

  private:
    std::unique_ptr<Layer> next_layer;
  };

  class Layer {
  public:
    template <class F>
    void insert(const T &value, const F &distance) {
      auto &node = find_best_node</* AllowEmpty = */ true>(
          [&value, &distance](const T &query) {
            return distance(value, query);
          });
      node.insert(value, distance);
    }

    template <class F>
    const T &locate(const F &f) const {
      const auto &node = find_best_node</* AllowEmpty = */ false>(f);
      return node.locate(f);
    }

    bool empty() const {
      return std::all_of(std::begin(nodes),
                         std::end(nodes),
                         [](const Node &node) { return node.is_empty(); });
    }

  protected:
    const Node &first_leaf_node() const {
      for (const auto &node : nodes) {
        if (node.is_leaf()) {
          return node;
        }
      }

      LOG_ERR("Failed to find an empty node.");
    }

    Node &first_leaf_node() {
      return const_cast<Node &>(std::as_const(*this).first_leaf_node());
    }

    template <bool AllowEmpty, class F>
    const Node &find_best_node(const F &f) const {
      auto index = 0;
      auto best_distance
          = std::numeric_limits<decltype(f(nodes[index].value))>::max();

      for (int i = 0; i < WIDTH; ++i) {
        if (nodes[i].is_empty()) {
          if constexpr (AllowEmpty) {
            return nodes[i];
          } else {
            break;
          }
        }

        auto current_distance = f(nodes[i].value);
        if (current_distance < best_distance) {
          index = i;
          best_distance = current_distance;
        }
      }

      return nodes[index];
    }

    template <bool AllowEmpty, class F>
    Node &find_best_node(const F &f) {
      return const_cast<Node &>(
          std::as_const(*this).template find_best_node<AllowEmpty>(f));
    }

    bool is_leaf_layer() const {
      for (const auto &node : nodes) {
        if (node.is_leaf()) {
          return true;
        }
      }
      return false;
    }

  private:
    Node nodes[WIDTH];
  };

public:
  template <class F>
  void insert(const T &value, const F &distance) {
    top_layer.insert(value, distance);
  }

  template <class F>
  const T &locate(const F &f) const {
    return top_layer.locate(f);
  }

  bool empty() const { return top_layer.empty(); }

private:
  Layer top_layer;
};
}

#endif // ZISA_TREE_HPP
