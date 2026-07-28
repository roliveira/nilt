#ifndef CONFTEST_HEADER
#define CONFTEST_HEADER

#include <iostream>
#include <fstream>

#include <nlohmann/json.hpp>


static nlohmann::json load_tol_file() {
    std::ifstream f("tolerances.json");
    return f.is_open() ? nlohmann::json::parse(f) : nlohmann::json::object();
}

static const nlohmann::json TOL = load_tol_file();


#endif // CONFTEST_HEADER
