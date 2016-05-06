template <typename T>
class Array2D
{
 public:
  // constructor
  Array2D(){
    d1 = 0;
    d2 = 0;
    A  = 0;
  }

  Array2D(unsigned l1,
	  unsigned l2){
    d1 = l1;
    d2 = l2;
    A  = 0;
    if (d1 > 0 &&
	d2 > 0) A = new T[d1*d2];
  }
 
  // destructor
  ~Array2D(){
    d1 = 0;
    d2 = 0;
    if (A) delete[] A;
    A  = 0;
  }

  // allocate memory
  void allocate(unsigned l1,
		unsigned l2){
    d1 = l1;
    d2 = l2;
    if (!A &&
	d1 > 0 &&
	d2 > 0) A = new T[d1*d2];
  }

  // deallocate memory
  void deallocate(){
    d1 = 0;
    d2 = 0;
    if (A) delete[] A;
    A = 0;
  }
  
  // indexing (parenthesis operator), two of them (for const correctness)
  const T& operator () (unsigned i,
			unsigned j) const{
    return A[d1*j+i];}

  T& operator () (unsigned i,
		  unsigned j){
    return A[d1*j+i];}
  
  // get dims
  unsigned getD1() const {return d1;}
  unsigned getD2() const {return d2;}

  
  // private data members
 private:

  unsigned d1;
  unsigned d2;
  T*       A;


  // to prevent unwanted copying:
  Array2D(const Array2D<T>&);
  Array2D& operator = (const Array2D<T>&);
};
