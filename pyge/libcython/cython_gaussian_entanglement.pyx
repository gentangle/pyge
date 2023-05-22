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
cdef double cython_gaussian_entanglement(
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

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)    # Deactivate Python division
def cython_loops_iteration(
    long int[:, :] loops,
    double[:, :] midpos,
    double[:, :] bonds,
    Py_ssize_t thr_min_len,
    Py_ssize_t len_chain,
    long int[:, :] loop_res,
    long int[:, :] thr_n_res,
    double[:] ge_n_res,
    long int[:, :] thr_c_res,
    double[:] ge_c_res
    ):

    cdef double GE_max_n, GE_max_c, GE_loop
    cdef Py_ssize_t i1, i2, j1, j2, j1_max_n, j2_max_n, j1_max_c, j2_max_c, increment, idx

    cdef Py_ssize_t len_loops = loops.shape[0]

    for idx in range(len_loops):
        # loop extrema
        i1 = loops[idx][0]
        i2 = loops[idx][1]

        GE_max_n = 0
        j1_max_n = 0
        j2_max_n = 0
        GE_max_c = 0
        j1_max_c = 0
        j2_max_c = 0
        # Search before the loop: N terminus
        for j1 in range(i1 - thr_min_len):
            for increment in range(i1 - thr_min_len - j1):
                j2 = j1 + thr_min_len + increment
                GE_loop = cython_gaussian_entanglement(i1, i2, j1, j2, midpos, bonds)
                if abs(GE_loop) > abs(GE_max_n):
                    GE_max_n = GE_loop
                    j1_max_n = j1
                    j2_max_n = j2
        # search after the loop
        for j1 in range(i2 + 1, len_chain - thr_min_len):
            for increment in range(len_chain - thr_min_len - j1):
                j2 = j1 + thr_min_len + increment
                GE_loop = cython_gaussian_entanglement(i1, i2, j1, j2, midpos, bonds)
                if abs(GE_loop) > abs(GE_max_c):
                    GE_max_c = GE_loop
                    j1_max_c = j1
                    j2_max_c = j2
        # loop's GE is the maximum in modulus
        loop_res[idx][0] = i1
        loop_res[idx][1] = i2
        thr_n_res[idx][0] = j1_max_n
        thr_n_res[idx][1] = j2_max_n
        ge_n_res[idx] = GE_max_n
        thr_c_res[idx][0] = j1_max_c
        thr_c_res[idx][1] = j2_max_c
        ge_c_res[idx] = GE_max_c
        
