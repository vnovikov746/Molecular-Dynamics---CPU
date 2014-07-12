/*
 * Molecular Dynamics Project
 * Author: Vladimir Novikov
 * This file contains the ........
 */

#ifndef LISTS_HANDLER_H_
#define LISTS_HANDLER_H_

#include "Structures.h"
#include "Configurations.h"
#include "Constants.h"
#include <string>
#include <fstream>

int readConfig(string,configurations*);
int readInput(string, allLists*, configurations*);
int listsAlloc(configurations*, allLists*);
int listsDestroy(configurations*, allLists*);


#endif //LISTS_HANDLER_H_