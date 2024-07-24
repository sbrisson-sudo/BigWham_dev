#include <string>
#include <stdexcept>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

#include <vector>

#include "il/Array.h"

#include "io/bigwham_io_gen.h"
#include "io/bigwham_io_helper.h"

using namespace bigwham;
template <typename T>
using jlarray = jlcxx::ArrayRef<T>;

/* -------------------------------------------------------------------------- */

template <typename T>
inline jlarray<T> as_jlarray(il::Array<T> &seq)
{
    return jlcxx::make_julia_array(seq.Data(), seq.size());
}
/* -------------------------------------------------------------------------- */

template <typename T>
inline il::ArrayView<T> as_array_view(const jlarray<T> &c)
{
    il::ArrayView<T> d{&c[0], c.size()};
    return d;
}
/* -------------------------------------------------------------------------- */

std::string greet()
{
    return "hello, world";
}

void test_array_set(jlcxx::ArrayRef<double> a, const int64_t i, const double v)
{
    a[i] = v;
}

struct World
{
    World(const std::string &message = "default hello") : msg(message) {}
    void set(const std::string &msg) { this->msg = msg; }
    std::string greet() { return msg; }
    std::string msg;
    il::Array<double> a;
    il::Array<double> matvec(const il::ArrayView<double> &b)
    {
        this->a.Resize(b.size());
        for (il::int_t i = 0; i < b.size(); i++)
        {
            this->a[i] = b[i] * 2;
        }
        return this->a;
    }
    ~World() { std::cout << "Destroying World with message " << msg << std::endl; }
};

JLCXX_MODULE define_julia_module(jlcxx::Module &mod)
{
    mod.method("greet", &greet, "documentation for greet");
    mod.method("test_array_set", &test_array_set, "documentation for test_array_set");
    mod.method("mat_vec", [](const jlarray<double> b)
               {
                   auto c = as_array_view<double>(b);
                   auto d = il::Array<double>(c.size(), 0.);
                   for (il::int_t i = 0; i < c.size(); i++)
                   {
                       d[i] = c[i];
                       std::cout << b[i] << c[i] << d[i] << std::endl;
                   }
                   auto e = jlcxx::make_julia_array(d.Data(), d.size());
                   for (il::int_t i = 0; i < d.size(); i++)
                   {
                       std::cout << e[i] << std::endl;
                   }
                   return e;
                   // return as_jlarray<double>(d);
               },
               "documentation for as_jlarray");
    mod.add_type<World>("World")
        .constructor<const std::string &>()
        .method("set", &World::set)
        .method("greet", &World::greet)
        .method("matvec", [](World &w, const jlarray<double> b)
                {
      auto c = as_array_view<double>(b);
      w.matvec(c);
      for (il::int_t i = 0; i < w.a.size(); i++) {
        std::cout << b[i] << c[i] << w.a[i] << std::endl;
      }
      return as_jlarray<double>(w.a); })
        .method("get_a", [](World &w)
                { return as_jlarray<double>(w.a); });
    mod.add_type<BigWhamIOGen>("BEMatrix")
        .constructor([](const jlarray<double> x, const jlarray<int> y, const std::string &z, const jlarray<double> w)
                     {
                         auto x_ = std::vector<double>(x.data(), x.data() + x.size());
                         auto y_ = std::vector<int>(y.data(), y.data() + y.size());
                         auto w_ = std::vector<double>(w.data(), w.data() + w.size());
                         return new BigWhamIOGen(x_, y_, z, w_);
                     })
        .method("build_hierarchical_matrix", &BigWhamIOGen::BuildHierarchicalMatrix)
        .method("get_collocation_points", &BigWhamIOGen::GetCollocationPoints)
        .method("matvec", [](BigWhamIOGen &w, const jlarray<double> x){
      auto c = as_array_view<double>(x);
      w.MatVec(c);
      return as_jlarray<double>(w.m_yout_); });
    ;
}