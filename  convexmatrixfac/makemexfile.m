if ispc
    
    mex -g -O qpmex.cpp m_dw.c dwqp.c solveqp.c  ...
        -lmwlapack  -largeArrayDims ...
        D:\lapack\pthreadVC2.lib -ID:\lapack
else
    !gcc -O -c -fpic  m_dw.c
    !gcc -O -c -fpic  dwqp.c
    !gcc -O -c -fpic  solveqp.c
    
    mex -O qpmex.cpp dwqp.o solveqp.o m_dw.o -largeArrayDims -lmwlapack  -lpthread
    !rm *.o    
end

% mex -O m_dwmex.cpp m_dw.c dwqp.c solveqp.c  -lmwlapack  -largeArrayDims  D:\lapack\pthreadVC2.lib -ID:\lapack
