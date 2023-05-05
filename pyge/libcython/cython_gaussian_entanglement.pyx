'''Cython implementation of the G' computation'''
cimport cython

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef void cross(
    double[:] x, 
    double[:] a, 
    double[:] b):
    '''Computes the cross product'''
    x[0] = a[1]*b[2] - a[2]*b[1]
    x[1] = a[2]*b[0] - a[0]*b[2]
    x[2] = a[0]*b[1] - a[1]*b[0]

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef double dot(
    double[:] a, 
    double[:] b):
    '''Compute dot product'''
    return a[0]*b[0] + a[1]*b[1] + a[2]*b[2]

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)    # Deactivate Python division
def cython_gaussian_entanglement(
    Py_ssize_t i1, 
    Py_ssize_t i2, 
    Py_ssize_t j1, 
    Py_ssize_t j2, 
    double[:, :] R, 
    double[:, :] deltaR):
    '''Returns G' for loop starting from residue i1 to i2 and 
    thread starting from residue j1 to j2'''

    cdef Py_ssize_t i, j, k
    cdef double ge = 0
    cdef double PI = 3.141592653589793

    cdef double dR_C[3] # C array implementation
    cdef double[:] dR = dR_C # memoryview for efficient memory access
    cdef double ddR_C[3]
    cdef double[:] ddR = ddR_C

    for i in range(i1, i2):
        for j in range(j1, j2):
            
            for k in range(3):
                dR[k] = R[i][k] - R[j][k]

            cross(ddR, deltaR[i], deltaR[j])
            ge += dot(dR, ddR)/(dot(dR, dR)**(1.5))

    return ge/(4*PI)