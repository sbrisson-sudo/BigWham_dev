(* ::Package:: *)

(* Specify build settings such as compiler and linker flags, libraries to be linked, etc. here *)



(*- SPECIFYING MKL AND TBB ROOT here -- customize as needed *)
Switch[$SystemID, "Linux-x86-64", {
	$MKLROOT = "@MKLROOT@";
  $TBBROOT = "/opt/intel/compilers_and_libraries_2018.2.199/linux/tbb";
  $TBBLIBPATH="/opt/intel/compilers_and_libraries_2018.2.199/linux/tbb/lib/intel64_lin/gcc4.7/libtbb.so.2";
  }
  , "MacOSX-x86-64", {
  $MKLROOT = "@MKLROOT@";
  $TBBROOT = "/opt/intel/compilers_and_libraries_2020/mac/tbb";
  $TBBLIBPATH= "/opt/intel/compilers_and_libraries_2020/mac/tbb/lib/libtbb.dylib";
      }];


(* PATH of the BigWham root folder *)
(* edit as see fit *)
$BigWhamDirectory =  "@CMAKE_BINARY_DIR@";



(*
Switch[$OperatingSystem,
  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = {
    "CompileOptions" -> {(* "-std=c++14", "-mmacosx-version-min=10.9" *)}

    (*
    , "IncludeDirectories" -> {}
    , "LibraryDirectories" -> {}
    *)
  },

  "Unix", (* Compilation settings for Linux *)
  $buildSettings = {
    "CompileOptions" -> {(* "-std=c++14" *)}

    (*
    , "IncludeDirectories" -> {}
    , "LibraryDirectories" -> {}
    *)
  },

  "Windows", (* Compilation settings for Windows *)
  $buildSettings = { 
    "CompileOptions" -> {"/EHsc", "/wd4244", "/DNOMINMAX"}

    (*
    , "IncludeDirectories" -> {}
    , "LibraryDirectories" -> {}
    *)
  }
]
*)
