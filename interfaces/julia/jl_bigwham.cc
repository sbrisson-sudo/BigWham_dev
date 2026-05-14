#include <string>
#include <stdexcept>

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/stl.hpp"

#include <vector>
#include <climits>

#include "il/Array.h"
#include "io/bigwham_io.h"
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
    il::ArrayView<T> d{&c[0], static_cast<int>(c.size())};
    return d;
}
/* -------------------------------------------------------------------------- */

template <typename T>
inline il::ArrayEdit<T> as_array_edit(jlarray<T> &c)
{
    il::ArrayEdit<T> d{&c[0], static_cast<int>(c.size())};
    return d;
}
/* -------------------------------------------------------------------------- */

JLCXX_MODULE define_julia_module(jlcxx::Module &mod)
{
    mod.add_type<BigWhamIO>("BigWhamIO")
        .constructor([](const jlarray<double> x, const jlarray<int> y, const std::string &z, const jlarray<double> w, const int num_threads)
                     {
                         auto x_ = std::vector<double>(x.data(), x.data() + x.size());
                         auto y_ = std::vector<int>(y.data(), y.data() + y.size());
                         auto w_ = std::vector<double>(w.data(), w.data() + w.size());
                         return new BigWhamIO(x_, y_, z, w_, num_threads, num_threads, true, false,false,-1); })
        .method("build_hierarchical_matrix", &BigWhamIO::BuildHierarchicalMatrix)
        .method("get_collocation_points", &BigWhamIO::GetCollocationPoints)
        .method("get_element_normals", &BigWhamIO::GetElementNormals)
        .method("get_rotation_matrix", &BigWhamIO::GetRotationMatrix)
        .method("get_omp_threads", &BigWhamIO::GetOmpThreads)
        .method("matvec!", [](BigWhamIO &w, const jlarray<double> xin, jlarray<double> xout)
                {
                    auto cin = as_array_view<double>(xin);
                    auto cout  = as_array_edit<double>(xout);
                    w.MatVecVoid(cin, cout);
                })
        .method("matvec", [](BigWhamIO &w, const jlarray<double> x)
                {
                    auto c = as_array_view<double>(x);
                    w.MatVecVoid(c);
                    return as_jlarray<double>(w.m_yout_); 
                });
    ;
}
