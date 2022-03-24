
#include <LTemplate.h>

// Dummy class for "free" functions
struct Manager {

    // Release a HMatExpr
    void releaseHMatExpr(mint id) {
        int err = mma::libData->releaseManagedLibraryExpression("HMatExpr", id);
        if (err)
            throw mma::LibraryError("Managed library expression does not exist.");
    }
};
