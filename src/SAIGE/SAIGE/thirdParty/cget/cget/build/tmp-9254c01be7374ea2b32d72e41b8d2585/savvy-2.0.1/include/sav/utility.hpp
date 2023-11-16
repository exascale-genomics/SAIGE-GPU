/*
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/.
 */

#ifndef SAVVY_SAV_UTILITY_HPP
#define SAVVY_SAV_UTILITY_HPP

#include "savvy/region.hpp"

#include <string>
#include <vector>
#include <unordered_set>

savvy::genomic_region string_to_region(const std::string& s);
std::string join_vector_to_string(const std::vector<std::string>& vec, std::string delim);
std::vector<std::string> split_string_to_vector(const char* in, char delim);
std::unordered_set<std::string> split_string_to_set(const char* in, char delim);
std::unordered_set<std::string> split_file_to_set(const char* in);
std::vector<std::string> split_file_to_vector(const char* in, std::size_t size_hint = 0);

#endif //SAVVY_SAV_UTILITY_HPP
