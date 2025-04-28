# Switch on test reporting
if !debug
    redirect_stdout(oldstd)
end;

test_keys = collect(keys(test_results));

# Define current and test structs
current = results;
test = test_results;

# Test names and values
try
    @testset "Results tests" begin
        for i in test_keys

            # Grab the key, replace NaN's
            i_v1 = replace(current[i], NaN => 0.0)
            i_v2 = replace(test[i], NaN => 0.0)

            # Run tests
            result = wrapper_test(i_v1, i_v2)

            # Warn if test failed
            if !(result isa Test.Pass)
                @warn "Test failed for field $i"
            end
        end
    end
catch e
    global all_passed = false
end;

# Switch off all other reporting
if !debug
    redirect_stdout(devnull)
end;
