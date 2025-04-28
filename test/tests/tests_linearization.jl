# Switch on test reporting
if !debug
    redirect_stdout(oldstd)
end;

# Tolerance across systems
if Sys.isapple()
    thistol = 1e-4
elseif Sys.islinux()
    thistol = 1e-3
elseif Sys.iswindows()
    thistol = 1e-4
end;

# Define current and test structs
current = lr_full;
test = test_lr_full;

# Test names and values
try
    @testset "Linearization tests" begin

        # Find shared and different names
        names_both = intersection(test, current)

        # Iterate over the common fields
        for i in names_both
            i_v1 = getproperty(test, i)
            i_v2 = getproperty(current, i)

            if is_named_struct(i_v1) || is_named_struct(i_v2)
                @error "You need to go deeper in the testing!"
                exit(1) # Exit with a nonzero code to make GitHub Actions fail
            end

            # Run test
            result = wrapper_test(i_v1, i_v2; tol = thistol)

            # Warn if test failed
            if !(result isa Test.Pass)
                @warn "Test failed for field $i at tolerance $thistol"
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
