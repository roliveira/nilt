/*
 * Verification - standard test functions for numerical inverse Laplace transforms.
 *
 * Evaluates all three inversion algorithms (Stehfest, Talbot, De Hoog) against
 * 10 Laplace transform pairs with known analytical solutions. The first five
 * are from Stehfest (1970) and the remaining five from Abate & Whitt (2006).
 *
 * Output: one CSV per (function, method) combination, containing columns
 *   t, fta (analytical), ftn (numerical), err (relative error).
 *
 * References:
 *   Stehfest, H. (1970). Commun. ACM 13(1), 47-49.
 *   Abate, J. & Whitt, W. (2006). INFORMS J. Comput. 18(4), 408-421.
 */
#include <cstdlib>
#include <iostream>

#include "nilt.hpp"
#include "testfunctions.hpp"
#include "utils.hpp"


int main()
{
    create_results_table("verification");
    return EXIT_SUCCESS;
}
