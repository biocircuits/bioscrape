from cython.operator cimport dereference as deref, preincrement as inc

cdef extern from "<vector>" namespace "std" nogil:
    cdef cppclass vector[T]:
        cppclass iterator:
            T operator*()
            iterator operator++()
            bint operator==(iterator)
            bint operator!=(iterator)
        vector()
        void push_back(T&)
        T& operator[](int)
        T& at(int)
        iterator begin()
        iterator end()
        void clear()
        size_t size()


