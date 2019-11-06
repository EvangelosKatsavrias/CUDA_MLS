function hide_Approximation(obj, approximationIndices)

for approximationIndex = approximationIndices
    delete(obj.handles_approximations(approximationIndex)); obj.handles_approximations(approximationIndex)=[];
end

end