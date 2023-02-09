for f in readdir(pwd())
    if endswith(f,"test.jl")
        @info "\033[94mreading $f\033[0m"
        try
            local t = @elapsed include(f)
            @info "\033[92m$f  $t(s)\033[0m"
        catch
            @error "Problem with "*f
            continue 
        end
    end
end

