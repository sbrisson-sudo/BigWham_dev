//
// This file is part of BigWham.
//
// Created by Brice Lecampion on 22.12.19.
// Copyright (c) EPFL (Ecole Polytechnique Fédérale de Lausanne) , Switzerland,
// Geo-Energy Laboratory, 2016-2020.  All rights reserved. See the LICENSE.TXT
// file for more details.
//

#include <iostream>
#include <string>

#include <src/BigWham.h>

/////////////////////////////////////////////////////////////////////////////////
int main() {

  std::cout << "++++++++++++++++++++\n";

  for (;;) {
    std::cout << "Please enter one of: new, load, save, or quit:\n";
    std::string option;
    std::cin >> option;
    switch (hash_djb2a(option)) {
    case "new"_sh:
      std::cout << "You entered \'new\'\n";
      break;
    case "load"_sh:
      std::cout << "You entered \'load\'\n";
      break;
    case "save"_sh:
      std::cout << "You entered \'save\'\n";
      break;
    case "quit"_sh:
      std::cout << "You entered \'quit\'\n";
      return 0;
    default:
      std::cout << "Command not recognized!\n";
      break;
    }
  }

  std::cout << "\n End of BigWham - exe "   << "\n";
}

/////////////////////////////////////////////////////////////////////////////////
