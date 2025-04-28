# Switch on test reporting
if !debug
    redirect_stdout(oldstd)
end;

# Define current and test structs
current = ss_full;
test = test_ss_full;

# Test names and values
try
    @testset "Steady state tests" begin

        # Find shared and different names
        names_both = intersection(test, current)

        # Iterate over the common fields
        for i in names_both
            i_v1 = getproperty(test, i)
            i_v2 = getproperty(current, i)

            # If the field is a struct, go deeper, else run test
            if is_named_struct(i_v1) || is_named_struct(i_v2)

                # Find shared and different names
                i_names_both = intersection(i_v1, i_v2)

                # Iterate over the common fields
                for j in i_names_both
                    j_v1 = getproperty(i_v1, j)
                    j_v2 = getproperty(i_v2, j)

                    # If the field is a struct, go deeper, else run test
                    if is_named_struct(j_v1) || is_named_struct(j_v2)

                        # Find shared and different names
                        j_names_both = intersection(j_v1, j_v2)

                        # Iterate over the common fields
                        for k in j_names_both
                            k_v1 = getproperty(j_v1, k)
                            k_v2 = getproperty(j_v2, k)

                            if is_named_struct(k_v1) || is_named_struct(k_v2)
                                @error "You need to go deeper in the testing!"
                                exit(1) # Exit with a nonzero code to make GitHub Actions fail
                            end

                            # Run test
                            result = wrapper_test(k_v1, k_v2)

                            # Warn if test failed
                            if !(result isa Test.Pass)
                                @warn "Test failed for field $i at $j at $k"
                            end
                        end
                    else
                        # Run test
                        result = wrapper_test(j_v1, j_v2)

                        # Warn if test failed
                        if !(result isa Test.Pass)
                            @warn "Test failed for field $i at $j"
                        end
                    end
                end

            else
                # Run test
                result = wrapper_test(i_v1, i_v2)

                # Warn if test failed
                if !(result isa Test.Pass)
                    @warn "Test failed for field $i"
                end
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
