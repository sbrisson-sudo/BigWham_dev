#include <vector>
#include <memory>
#include <iostream>

#include <il/Array2D.h>

#include "hmat/hmatrix/contiguous_array2d_vector.hpp"

int main(){

    // We create a vector of Array2D
    std::vector<std::unique_ptr<il::Array2D<double>>> full_rank_blocks_;

    // We populate it 
    const int N_BLOCKS = 5;
    full_rank_blocks_.reserve(N_BLOCKS);
    for (int n(0); n<N_BLOCKS; n++){
        auto new_block = std::make_unique<il::Array2D<double>>(n+1,n+1); 
        il::Array2DEdit<double> M = new_block->Edit();
        for (int i(0); i<n+1; i++){
            for (int j(0); j<n+1; j++){
                M(i,j) = static_cast<int>(n);
            }
        }
        full_rank_blocks_.push_back(std::move(new_block));
    }

    // Try to access an element
    std::cout << "Before copy : " << std::endl;
    for (int n(0); n<N_BLOCKS; n++){
        std::cout << "Block # " << n << " = [ ";
        auto block_view = full_rank_blocks_[n]->view();
        for (int i(0); i<n+1; i++){
            for (int j(0); j<n+1; j++){
                std::cout << block_view(i,j) << ", ";
            }
        }
        std::cout << "]" << std::endl;
        std::cout << "Stored at (rel to first) : " << full_rank_blocks_[n]->data() -  full_rank_blocks_[0]->data()<< std::endl;
    }

    // Try the contiguous operator 
    ContiguousArray2DVector<double> full_blocks_cont(full_rank_blocks_);
    auto& full_blocks_cont_vec = full_blocks_cont.blocks();

    // Delete previous 

    std::cout << "Copying vector done" << std::endl;

    std::cout << "After copy : " << std::endl;
    for (int n(0); n<N_BLOCKS; n++){
        std::cout << "Block # " << n << " = [ ";
        auto block_view = full_blocks_cont_vec[n]->view();
        for (int i(0); i<n+1; i++){
            for (int j(0); j<n+1; j++){
                std::cout << block_view(i,j) << ", ";
            }
        }
        std::cout << "]" << std::endl;
        std::cout << "Stored at (rel to first) : " << full_blocks_cont_vec[n]->data() - full_blocks_cont_vec[0]->data() << std::endl;

    }




    return 0;

}