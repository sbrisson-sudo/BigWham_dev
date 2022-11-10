;;; Directory Local Variables
;;; For more information see (info "(emacs) Directory Variables")

(
 ;; (python-mode . ((pyvenv-workon . "myscipy")))
 (nil . ((compile-command . (concat "cd build; source ~/pythevirons/myscipy/bin/activate; cmake -DLGFM_ROOT=" "\"~/epfl/projects/lgfm/build\"" " -DCMAKE_BUILD_TYPE=Release -DHLIBPRO_ROOT=" "\"~/epfl/projects/hlibpro-2.8.1\"" " -DCXX11_ABI=true ..; make -j2"))))
;; (nil
;;   (format-all-formatters
;;    ("C++" clang-format)
;;    ("C" clang-format)
;;    ("Python" (yapf "--style" "google"))))
 )
