#include <cxxwrap/julia.h>
#include "bigwham_io_gen.h"

JULIA_CPP_MODULE_BEGIN(registry)
    using namespace JLCxx;

    Module& mod = registry.create_module("BigWhamIO");

    mod.add_type<BigWhamIOGen>("BigWhamIOGen")
        .constructor<const std::vector<double>&, const std::vector<int>&, const std::string&, const std::vector<double>&>()
        .method("some_method", &BigWhamIOGen::some_method)  // Replace with actual methods
        // Add more methods as needed
        ;
JULIA_CPP_MODULE_END
