(* ::Package:: *)

(* Specify build settings such as compiler and linker flags, libraries to be linked, etc. here *)

(* PATH of the BigWham root folder *)
(* edit as see fit *)
$BigWhamDirectory =  "@PROJECT_BINARY_DIR@";

(*
Switch[$OperatingSystem,
  "MacOSX", (* Compilation settings for OS X *)
  $buildSettings = {
    "CompileOptions" -> {(* "-std=c++17", "-mmacosx-version-min=10.9" *)}

    (*
    , "IncludeDirectories" -> {}
    , "LibraryDirectories" -> {}
    *)
  },

  "Unix", (* Compilation settings for Linux *)
  $buildSettings = {
    "CompileOptions" -> {(* "-std=c++17" *)}

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
