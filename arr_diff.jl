using PrettyTables

# Read two arrays from standard input, compute their difference and print it to standard output.
# Input: Unless passed as cli args, first 2 lines contain n and m, the sizes of the 2D array
# Next n lines contain n integers each, representing the elements of the first array
# Next n lines contain n integers each, representing the elements of the second array
# Output: n lines containing n integers each, representing the elements of the difference of the two arrays

function arr_diff(n, m, arr1, arr2)
    diff = Array{Float64}(undef, n, m)
    reldiff = Array{Float64}(undef, n, m)
    for i in 1:n
        for j in 1:m
            diff[i, j] = abs(arr1[i, j] - arr2[i, j])
            max_abs = max(abs(arr1[i, j]), abs(arr2[i, j]))
            # avoid division by zero
            reldiff[i, j] = max_abs != 0 ? diff[i, j] / max_abs : 0
        end
    end
    return diff, reldiff
end

# Read input
if length(ARGS) == 2
    n = parse(Int64, ARGS[1])
    m = parse(Int64, ARGS[2])
elseif length(ARGS) == 1
    n = parse(Int64, ARGS[1])
    m = parse(Int64, ARGS[1])
else
    n = parse(Int64, readline()) # rows
    m = parse(Int64, readline()) # columns
end

# Trim, replace ", " with "," and split by " "
read_array_line() = split(replace(strip(readline()), Pair(", ", ",")), " ")

lines_1 = [read_array_line() for _ in 1:n]
lines_2 = [read_array_line() for _ in 1:n]
@assert all(length(line) == m for line in lines_1)
@assert all(length(line) == m for line in lines_2)

# Convert to 2D array and transpose (lists are column vectors)
lines_1 = permutedims(hcat(lines_1...))
lines_2 = permutedims(hcat(lines_2...))
println(size(lines_1))
println(size(lines_2))

function parse_complex(str)
    # Remove parentheses
    str = strip(str, ['(', ')'])
    # Split by comma
    parts = split(str, ",")
    # Convert to Float64 and create Complex number
    return Complex{Float64}(parse(Float64, parts[1]), parse(Float64, parts[2]))
end

# Numbers of form (a,b) are converted to complex numbers a + bi
# Otherwise, they are converted to real numbers

function parse_array(lines)
    if lines[1, 1][1] == '('
        return [parse_complex(lines[i, j]) for i in 1:n, j in 1:m]
    else
        return [parse(Float64, lines[i, j]) for i in 1:n, j in 1:m]
    end
end

arr1 = parse_array(lines_1)
arr2 = parse_array(lines_2)

println("Array 1:")
pretty_table(arr1)
println("Array 2:")
pretty_table(arr2)

diff, reldiff = arr_diff(n, m, arr1, arr2)
println("Difference:")
pretty_table(diff)
println("Relative difference:")
pretty_table(reldiff)

argmax_diff = argmax(diff)
print("Max diff at indices $(argmax_diff): $(diff[argmax_diff]) ")
println("(values=($(arr1[argmax_diff]), $(arr2[argmax_diff]))")
argmax_reldiff = argmax(reldiff)
print("Max relative diff at indices $(argmax_reldiff): $(reldiff[argmax_reldiff]) ")
println("(values=($(arr1[argmax_reldiff]), $(arr2[argmax_reldiff]))")
