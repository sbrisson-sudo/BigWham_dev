#include "hmat_selection.h"

// #define DEBUG

namespace bigwham {

// Check that asked spans are within range
// Important : passed as range on elements
template <typename T>
void HmatSelection<T>::validateIndices() const {

    const int dim_dof = this->dof_dimension_;

    #ifdef DEBUG 
    std::cout << "row_indices_ = [";
    for (int i=0; i<row_indices_.size(); i++) std::cout << row_indices_[i] << ", ";
    std::cout << "]" <<  std::endl; 
    std::cout << "col_indices_ = [";
    for (int i=0; i<col_indices_.size(); i++) std::cout << i <<":"<< col_indices_[i] << ", ";
    std::cout << "]" <<  std::endl; 

    std::cout << "col_indices_.size() = " << col_indices_.size() << std::endl;
    std::cout << "row_indices_.size() = " << row_indices_.size() << std::endl;
    std::cout << "this->base_size_[0] = " << this->base_size_[0] << std::endl;
    std::cout << "this->base_size_[1] = " << this->base_size_[1] << std::endl;
    std::cout << "dim_dof = " << dim_dof << std::endl;
    #endif

    for (int i=0; i<row_indices_.size(); i++){
        if (row_indices_[i] < 0 || row_indices_[i] >= this->base_size_[0]/dim_dof) {
            throw std::out_of_range("Row index out of bounds");
        }
    }
    
    for (int i=0; i<col_indices_.size(); i++){
        if (col_indices_[i] < 0 || col_indices_[i] >= this->base_size_[1]/dim_dof) {
            throw std::out_of_range("Column index out of bounds");
        }
    }
}

// Perform block selection
template <typename T>
void HmatSelection<T>::blockSelection() {

    const int dim_dof = this->dof_dimension_;

    // 1. We apply permutation to the input spans

    il::Array<il::int_t> inv_permutation{this->hr_->permutation_1_.size()};
    for (size_t i = 0; i < this->hr_->permutation_1_.size(); ++i) {
        inv_permutation[this->hr_->permutation_1_[i]] = i;
    }

    row_indices_perm_.Resize(row_indices_.size());
    col_indices_perm_.Resize(col_indices_.size());

    auto row_indices_perm_edit = row_indices_perm_.Edit();
    auto col_indices_perm_edit = col_indices_perm_.Edit();

    for (int i=0; i<row_indices_.size(); i++) {
        row_indices_perm_edit[i] = inv_permutation[row_indices_[i]];
    }
    for (int i=0; i<col_indices_.size(); i++) {
        col_indices_perm_edit[i] = inv_permutation[col_indices_[i]];
    }

    // 2. We define boolean array (1 true, 0 false)
    il::Array<int> x_test{static_cast<il::int_t>(this->base_size_[0]/this->dof_dimension_), 0};
    il::Array<int> y_test{static_cast<il::int_t>(this->base_size_[0]/this->dof_dimension_), 0};

    auto x_test_edit = x_test.Edit();
    auto y_test_edit = y_test.Edit();

    for (int i=0; i<row_indices_.size(); i++) {
        y_test_edit[row_indices_perm_[i]] = 1;
    }
    for (int i=0; i<col_indices_.size(); i++) {
        x_test_edit[col_indices_perm_[i]] = 1;
    }

    // 3. We loop on the blocks, if both x and y spans have one non zero element > we append the block

    auto x_test_view = x_test.view();
    auto y_test_view = y_test.view();
    int x_test_sum, y_test_sum;

    for (il::int_t i = 0; i < this->hr_->pattern_.n_FRB; i++) {
        auto i0 = this->hr_->pattern_.FRB_pattern(1, i);
        auto j0 = this->hr_->pattern_.FRB_pattern(2, i);
        auto iend = this->hr_->pattern_.FRB_pattern(3, i);
        auto jend = this->hr_->pattern_.FRB_pattern(4, i);

        x_test_sum = 0;
        for (int k=j0; k<jend; k++){
            x_test_sum += x_test[k];
        }
        y_test_sum = 0;
        for (int k=i0; k<iend; k++){
            y_test_sum += y_test[k];
        }

        if ((x_test_sum > 0) && (y_test_sum > 0)){
            fr_blocks_selected_.push_back(i);
        }
    }


    for (il::int_t i = 0; i < this->hr_->pattern_.n_LRB; i++) {
        auto i0 = this->hr_->pattern_.LRB_pattern(1, i);
        auto j0 = this->hr_->pattern_.LRB_pattern(2, i);
        auto iend = this->hr_->pattern_.LRB_pattern(3, i);
        auto jend = this->hr_->pattern_.LRB_pattern(4, i);

        x_test_sum = 0;
        for (int k=j0; k<jend; k++){
            x_test_sum += x_test[k];
        }
        y_test_sum = 0;
        for (int k=i0; k<iend; k++){
            y_test_sum += y_test[k];
        }

        if ((x_test_sum > 0) && (y_test_sum > 0)){
            lr_blocks_selected_.push_back(i);
        }
    }

    #ifdef DEBUG 
    std::cout << "fr_blocks_selected_ = [";
    for (auto i : fr_blocks_selected_) std::cout << i << ", ";
    std::cout << "]" <<  std::endl; 

    std::cout << "lr_blocks_selected_ = [";
    for (auto i : lr_blocks_selected_) std::cout << i << ", ";
    std::cout << "]" <<  std::endl; 
    #endif
}

template <typename T> 
il::Array<T> HmatSelection<T>::matvec_full(il::ArrayView<T> x_full) {

    int dof_dimension = this->dof_dimension_;

    il::Array<T> y_full(this->base_size_[0], 0.0, il::align_t(), 64);

    #pragma omp parallel num_threads(this->n_openMP_threads_)
    {

        il::Array<T> yprivate(this->base_size_[0], 0.0, il::align_t(), 64);

        #pragma omp for schedule(guided) nowait

            for (auto i : fr_blocks_selected_) {
                auto i0 = this->hr_->pattern_.FRB_pattern(1, i);
                auto j0 = this->hr_->pattern_.FRB_pattern(2, i);
                auto iend = this->hr_->pattern_.FRB_pattern(3, i);
                auto jend = this->hr_->pattern_.FRB_pattern(4, i);

                auto a = (*full_rank_blocks_ref_[i]).view();
                auto xs = x_full.view(il::Range{j0 * dof_dimension, jend * dof_dimension});
                auto ys = yprivate.Edit(il::Range{i0 * dof_dimension, iend * dof_dimension});

                il::blas(1.0, a, xs, 1.0, il::io, ys);
            }

        #pragma omp for schedule(guided) nowait
            for (auto ii : lr_blocks_selected_) {
                auto i0 = this->hr_->pattern_.LRB_pattern(1, ii);
                auto j0 = this->hr_->pattern_.LRB_pattern(2, ii);
                auto iend = this->hr_->pattern_.LRB_pattern(3, ii);
                auto jend = this->hr_->pattern_.LRB_pattern(4, ii);

                auto a = low_rank_blocks_ref_[ii]->A.view();
                auto b = low_rank_blocks_ref_[ii]->B.view();

                auto xs = x_full.view(il::Range{j0 * dof_dimension, jend * dof_dimension});
                auto ys =yprivate.Edit(il::Range{i0 * dof_dimension, iend * dof_dimension});
                auto r = a.size(1);
                il::Array<double> tmp{r, 0.0};

                il::blas(1.0, b, il::Dot::None, xs, 0.0, il::io, tmp.Edit());
                il::blas(1.0, a, tmp.view(), 1.0, il::io, ys);
            }

        #pragma omp critical
        il::blas(1., yprivate.view(), il::io_t{}, y_full.Edit());

    }

    return y_full;
}

template <typename T> 
il::Array<T> HmatSelection<T>::matvec(il::ArrayView<T> x) {
    IL_EXPECT_FAST(this->isBuilt_);
    IL_EXPECT_FAST(x.size() == this->size_[1]);

    int dof_dimension = this->dof_dimension_;

    // 1. We copy the values we need in a full x vector
    // Note : x is already permuted here
    il::Array<T> x_full(this->base_size_[0], 0.0, il::align_t(), 64);
    auto x_full_edit = x_full.Edit(); 
    auto x_view = x.view(); 
    for (int i=0; i<col_indices_.size(); i++) {
        for (int j=0; j<dof_dimension; j++){
            x_full_edit[dof_dimension*col_indices_[i]+j] = x_view[dof_dimension*i+j];
        } 
    }

    il::Array<T> y_full = matvec_full(x_full.view());

    // Now we extract only the values we need
    il::Array<T> y(this->size_[0], 0.0, il::align_t(), 64);
    auto y_edit = y.Edit(); 
    auto y_full_view = y_full.view(); 
    for (int i=0; i<row_indices_.size(); i++) {
        for (int j=0; j<dof_dimension; j++){
            y_edit[dof_dimension*i+j] = y_full_view[dof_dimension*row_indices_[i]+j]; 
        } 
    }

    return y;
}

// Same but with permutation :
template <typename T> 
il::Array<T> HmatSelection<T>::matvecOriginal(il::ArrayView<T> x) {

    IL_EXPECT_FAST(this->isBuilt_);
    IL_EXPECT_FAST(x.size() == this->size_[1]);

    int dof_dimension = this->dof_dimension_;

    // 1. We copy the values we need in a full x vector
    // Note : here we also permute x
    il::Array<T> x_full(this->base_size_[0], 0.0, il::align_t(), 64);
    auto x_full_edit = x_full.Edit(); 
    auto x_view = x.view(); 
    for (int i=0; i<col_indices_perm_.size(); i++) {
        for (int j=0; j<dof_dimension; j++){
            x_full_edit[dof_dimension*col_indices_perm_[i]+j] = x_view[dof_dimension*i+j];
        } 
    }

    il::Array<T> y_full = matvec_full(x_full.view());

    // Now we extract only the values we need
    il::Array<T> y(this->size_[0], 0.0, il::align_t(), 64);
    auto y_edit = y.Edit(); 
    auto y_full_view = y_full.view(); 
    for (int i=0; i<row_indices_perm_.size(); i++) {
        for (int j=0; j<dof_dimension; j++){
            y_edit[dof_dimension*i+j] = y_full_view[dof_dimension*row_indices_perm_[i]+j]; 
        } 
    }

    return y;
}


/*
@brief This metohd returns the FR blocks data that falls wihin the selection
*/
template <typename T> 
void HmatSelection<T>::fullBlocksOriginal(il::io_t, il::Array<T> & val_list,il::Array<int> & pos_list){

    const int dim_dof = this->dof_dimension_;

    // 1. We compute the number of entries
    // 1.a. We define boolean array (1 true, 0 false)
    il::Array<int> x_test{static_cast<il::int_t>(this->base_size_[0]/this->dof_dimension_), 0};
    il::Array<int> y_test{static_cast<il::int_t>(this->base_size_[0]/this->dof_dimension_), 0};

    auto x_test_edit = x_test.Edit();
    auto y_test_edit = y_test.Edit();

    for (int i=0; i<row_indices_.size(); i++) {
        y_test_edit[row_indices_perm_[i]] = 1;
    }
    for (int i=0; i<col_indices_.size(); i++) {
        x_test_edit[col_indices_perm_[i]] = 1;
    }

    // 1.b. We loop on the blocks, if both x and y spans have one non zero element > we append the block

    auto x_test_view = x_test.view();
    auto y_test_view = y_test.view();
    int nb_entry_col, nb_entry_row;

    int nbfentry = 0;

    for (auto i : fr_blocks_selected_) {
        auto i0 = this->hr_->pattern_.FRB_pattern(1, i);
        auto j0 = this->hr_->pattern_.FRB_pattern(2, i);
        auto iend = this->hr_->pattern_.FRB_pattern(3, i);
        auto jend = this->hr_->pattern_.FRB_pattern(4, i);

        nb_entry_col = 0;
        for (int k=j0; k<jend; k++){
            nb_entry_col += x_test[k];
        }
        nb_entry_row = 0;
        for (int k=i0; k<iend; k++){
            nb_entry_row += y_test[k];
        }

        nbfentry += nb_entry_col*nb_entry_row * dim_dof*dim_dof;
    }

    #ifdef DEBUG 
    std::cout << "nbfentry = " << nbfentry << std::endl;
    #endif

    // prepare outputs
    pos_list.Resize(nbfentry * 2);
    val_list.Resize(nbfentry);

    // inverse the (partial) permutations
    il::Array<int> permut_col_inv{dim_dof * this->hr_->permutation_0_.size(), -1};
    il::Array<int> permut_row_inv{dim_dof * this->hr_->permutation_0_.size(), -1};
    for (int i=0; i<col_indices_perm_.size(); i++) {
        for (int j=0; j<this->dof_dimension_; j++){
            // std::cout << "col_indices_perm_[i]*dim_dof+j = " << col_indices_perm_[i]*dim_dof+j << std::endl;
            permut_col_inv[col_indices_perm_[i]*dim_dof+j] = i*dim_dof+j;
        }
    }
    for (int i=0; i<row_indices_perm_.size(); i++) {
        for (int j=0; j<this->dof_dimension_; j++){
            // std::cout << "row_indices_perm_[i]*dim_dof+j = " << row_indices_perm_[i]*dim_dof+j << std::endl;
            permut_row_inv[row_indices_perm_[i]*dim_dof+j] = i*dim_dof+j;
        }
    }

    // copy the data 
    int counter = 0;

    for (auto i : fr_blocks_selected_) {
        auto i0 = this->hr_->pattern_.FRB_pattern(1, i);
        auto j0 = this->hr_->pattern_.FRB_pattern(2, i);
        auto iend = this->hr_->pattern_.FRB_pattern(3, i);
        auto jend = this->hr_->pattern_.FRB_pattern(4, i);

        auto a = (*full_rank_blocks_ref_[i]).view();

        for (int ii=i0; ii<iend; ii++){
            for (int jj=j0; jj<jend; jj++){
                if (x_test[jj]+y_test[ii] == 2){

                    // We take these values
                    for (int iii=0; iii<dim_dof; iii++){
                        for (int jjj=0; jjj<dim_dof; jjj++){

                            // std::cout << "counter = " << counter << std::endl;
                        
                            // We take this value
                            val_list[counter] = a((ii-i0)*dim_dof + iii, (jj-j0)*dim_dof + jjj);

                            // We set its position (in the selected, original dof ordering)
                            // std::cout << "ii*dim_dof + iii = " << ii*dim_dof + iii << std::endl;
                            // std::cout << "jj*dim_dof + jjj = " << jj*dim_dof + jjj << std::endl;
                            pos_list[2*counter] = permut_row_inv[ii*dim_dof + iii];
                            pos_list[2*counter+1] = permut_col_inv[jj*dim_dof + jjj];

                            // We increment the counter 
                            counter ++;
                        }
                    }
                }
            }
        }
    }

}

template class HmatSelection<double>;

}