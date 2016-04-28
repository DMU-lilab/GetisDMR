#ifndef FINDFILECURRENT_H
#define FINDFILECURRENT_H


#include <iostream>
#include <string>
#include <fstream>

#define BOOST_FILESYSTEM_NO_DEPRECATED
#include "boost/filesystem.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/progress.hpp"
#include <iostream>

#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

string FindFileHome(string exp);

#endif
