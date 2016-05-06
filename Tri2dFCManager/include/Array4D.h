template <typename T>
class Array4D
{
 public:
  // constructor
  Array4D(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    d4 = 0;
    A  = 0;
  }

  Array4D(unsigned l1,
	  unsigned l2,
	  unsigned l3,
	  unsigned l4){
    d1 = l1;
    d2 = l2;
    d3 = l3;
    d4 = l4;
    A  = 0;
    if (d1 > 0 &&
	d2 > 0 &&
	d3 > 0 &&
	d4 > 0) A = new T[d1*d2*d3*d4];
  }
 
  // destructor
  ~Array4D(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    d4 = 0;
    if (A) delete[] A;
    A  = 0;
  }

  // allocate memory
  void allocate(unsigned l1,
		unsigned l2,
		unsigned l3,
		unsigned l4){
    d1 = l1;
    d2 = l2;
    d3 = l3;
    d4 = l4;
    if (!A &&
	d1 > 0 &&
	d2 > 0 &&
	d3 > 0 &&
	d4 > 0) A = new T[d1*d2*d3*d4];
  }

  // deallocate memory
  void deallocate(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    d4 = 0;
    if (A) delete[] A;
    A = 0;
  }
  
  // indexing (parenthesis operator), two of them (for const correctness)
  const T& operator () (unsigned i,
			unsigned j,
			unsigned k,
			unsigned l) const{
    return A[d4*(d3*(d2*i+j)+k)+l];}

  T& operator () (unsigned i,
		  unsigned j,
		  unsigned k,
		  unsigned l){
    return A[d4*(d3*(d2*i+j)+k)+l];}
  
  // get dims
  unsigned getD1() const {return d1;}
  unsigned getD2() const {return d2;}
  unsigned getD3() const {return d3;}
  unsigned getD4() const {return d4;}


  // private data members
 private:

  unsigned d1;
  unsigned d2;
  unsigned d3;
  unsigned d4;
  T*       A;


  // to prevent unwanted copying:
  Array4D(const Array4D<T>&);
  Array4D& operator = (const Array4D<T>&);
};
