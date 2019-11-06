#include <vector>
#include <iostream>

template <class T>
void wrapArrayInVector( T *sourceArray, size_t arraySize, std::vector<T, std::allocator<T> > &targetVector ) {
  typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
    (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)((void *) &targetVector);
  vectorPtr->_M_start = sourceArray;
  vectorPtr->_M_finish = vectorPtr->_M_end_of_storage = vectorPtr->_M_start + arraySize;
}

template <class T>
void releaseVectorWrapper( std::vector<T, std::allocator<T> > &targetVector ) {
  typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *vectorPtr =
        (typename std::_Vector_base<T, std::allocator<T> >::_Vector_impl *)((void *) &targetVector);
  vectorPtr->_M_start = vectorPtr->_M_finish = vectorPtr->_M_end_of_storage = NULL;
}

/*
int main() {

  int tests[6] = { 1, 2, 3, 6, 5, 4 };
  std::vector<int> targetVector;
  wrapArrayInVector( tests, 6, targetVector);

  std::cout << std::hex << &tests[0] << ": " << std::dec
            << tests[1] << " " << tests[3] << " " << tests[5] << std::endl;

  std::cout << std::hex << &targetVector[0] << ": " << std::dec
            << targetVector[1] << " " << targetVector[3] << " " << targetVector[5] << std::endl;

  releaseVectorWrapper( targetVector );
}

*/
