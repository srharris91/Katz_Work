template <typename T>
class Array1D
{
 public:
  // constructor
  Array1D(){
    d1 = 0;
    A  = 0;
  }

  Array1D(unsigned l1){
    d1 = l1;
    A  = 0;
    if (d1 > 0) A = new T[d1];
  }
 
  // destructor
  ~Array1D(){
    d1 = 0;
    if (A) delete[] A;
    A  = 0;
  }

  // allocate memory
  void allocate(unsigned l1){
    d1 = l1;
    if (!A && d1 > 0) A = new T[d1];
  }

  // deallocate memory
  void deallocate(){
    d1 = 0;
    if (A) delete[] A;
    A = 0;
  }
  
  // indexing (parenthesis operator), two of them (for const correctness)
  const T& operator () (unsigned i) const {return A[i];}

  T& operator () (unsigned i) {return A[i];}
  
  // get dims
  unsigned getD1() const {return d1;}
  
  // private data members
 private:

  unsigned d1;
  T*       A;


  // to prevent unwanted copying:
  Array1D(const Array1D<T>&);
  Array1D& operator = (const Array1D<T>&);
};
