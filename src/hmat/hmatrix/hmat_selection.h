#ifndef HMAT_SELECTION
#define HMAT_SELECTION

#include "hmat.h"

namespace bigwham {

template <typename T>
class HmatSelection : public Hmat<T> {
private:
    il::Array<int> row_indices_;
    il::Array<int> col_indices_;
    il::Array<int> row_indices_perm_;
    il::Array<int> col_indices_perm_;
    const Hmat<T>& base_hmat_;  // Reference to the original matrix

    il::StaticArray<il::int_t, 2> base_size_;
    
    // References to the block containers from the base matrix
    const std::vector<std::unique_ptr<bigwham::LowRank<T>>>& low_rank_blocks_ref_;
    const std::vector<std::shared_ptr<il::Array2D<T>>>& full_rank_blocks_ref_;

    // List of blocks to compute on 
    std::vector<int> lr_blocks_selected_;
    std::vector<int> fr_blocks_selected_;

    // Utilities
    void validateIndices() const;
    void blockSelection();

public:

    // Constrcutor from an Hmat object 
    HmatSelection(const Hmat<T>& base_hmat, 
                  const il::Array<int>& row_selection, 
                  const il::Array<int>& col_selection)
        : row_indices_(row_selection), col_indices_(col_selection), base_hmat_(base_hmat),
          low_rank_blocks_ref_(base_hmat.get_low_rank_blocks()),
          full_rank_blocks_ref_(base_hmat.get_full_rank_blocks()) {
        
        // Copy/share the hierarchical representation from base
        this->hr_ = base_hmat.get_hr();
        this->dof_dimension_ = base_hmat.get_dof_dimension();

        // Store "full" size
        this->base_size_[0] = base_hmat.size(0);
        this->base_size_[1] = base_hmat.size(1);

        // Validate selection indices
        validateIndices();
        
        // Update size to reflect the selection
        this->size_[0] = static_cast<il::int_t>(row_indices_.size()) * this->dof_dimension_;
        this->size_[1] = static_cast<il::int_t>(col_indices_.size()) * this->dof_dimension_;

        // Copy other attributes
        this->isBuilt_ = base_hmat.get_isBuilt();
        this->isBuilt_LR_ = base_hmat.get_isBuilt_LR();
        this->isBuilt_FR_ = base_hmat.get_isBuilt_FR();
        this->n_openMP_threads_ = base_hmat.get_n_openMP_threads();
        this->frb_chunk_size_ = base_hmat.get_frb_chunk_size();
        this->lrb_chunk_size_ = base_hmat.get_lrb_chunk_size();
        this->verbose_ = base_hmat.get_verbose();
        this->fixed_rank_ = base_hmat.get_fixed_rank();

        // Finally we call the block selection method
        blockSelection();
    }

    ~HmatSelection() = default;

    // Overridding the matvec function
    il::Array<T> matvec(il::ArrayView<T> x) override;
    il::Array<T> matvecOriginal(il::ArrayView<T> x) override;
    il::Array<T> matvec_full(il::ArrayView<T> x);

    // Getting the full blocks for preconditionners 
    void fullBlocksOriginal(il::io_t, il::Array<T> & val_list,il::Array<int> & pos_list) override;
};

}

#endif