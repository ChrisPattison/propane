#pragma once
#include <string>
#include <vector>

namespace propane::util
{
/** Basic implemention of assert used for validating inputs.
*/
void Check(bool cond, const char* err);

enum StringSplitOptions {
    kStringSplitOptions_None,
    kStringSplitOptions_RemoveEmptyEntries
};

/** Splits string along seperator.
 * Does not return seperator.
 * Optionally removes empty entries.
 */
std::vector<std::string> Split(std::string target, char seperator, StringSplitOptions options = kStringSplitOptions_None);
/** Split with the option of multiple seperators.
*/
std::vector<std::string> Split(std::string target, std::vector<char> seperator, StringSplitOptions options = kStringSplitOptions_None);
}