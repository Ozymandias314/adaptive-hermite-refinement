# Usage: ./format.sh <[0|1]>
# if set to 1, the script only checks formatting and fails on error
# if set to 0, the script formats the code in-place
CHECK=${1:-0} # in-place by default
if [ $CHECK -eq 0 ]; then
  FLAGS="-i"
else
  FLAGS="--dry-run --Werror"
fi
find cpp/bench cpp/test/ cpp/lib -name "*.cpp" -o -name "*.h" -o -name "*.hpp" | xargs clang-format-15 $FLAGS