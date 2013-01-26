#include "moleculesystem.h"

#include <iostream>

using namespace std;

int main()
{
    MoleculeSystem system;
    system.load("systems/argon/states/initial.xyz");

    return 0;
}

