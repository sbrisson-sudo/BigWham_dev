;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

(
 ;; (python-mode . ((pyvenv-workon . "myscipy")))
 (nil . ((compile-command . (concat "source bigwham_venv/bin/activate;. /opt/intel/mkl/bin/mklvars.sh intel64; cd build; cmake -DCMAKE_CXX_COMPILER=g++ -DBUILD_GOOGLE_TESTS=1 -DBUILD_PYTHON_BINDINGS=0 -DUSE_INTEL=0 -DIL_OPENMP=1 -DIL_OPENBLAS=0 -DIL_MKL=1  ..;make -j2; ./BigWhamElastUnitTest"))))

;; (nil
;;   (format-all-formatters
;;    ("C++" clang-format)
;;    ("C" clang-format)
;;    ("Python" (yapf "--style" "google"))))
 )
