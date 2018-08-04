try
echo on
mex pinpolyhedron.cpp tmain.cpp kodtree.cpp
echo off
catch exception
    if(~isempty(exception))
        fprintf(['\n Error during compilation, be sure :\n'...
            'i)  You have C compiler installed (prefered compiler are MSVC/Intel/GCC)\n'...
            'ii) You did "mex -setup" in matlab prompt before running mexpinpolyhedron.cpp'...
            ' tmain.cpp kodtree.cpp\n']);
    end    
end