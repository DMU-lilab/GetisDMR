#ifndef BetaBinomialReg_H
#define BetaBinomialReg_H

#include <vector>

class Regression;
bool BetaBinomialReg(Regression &r, std::vector<double> initial_parameters = std::vector<double> ());
#endif 
