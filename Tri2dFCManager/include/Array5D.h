template <typename T>
class Array5D
{
 public:
  // constructor
  Array5D(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    d4 = 0;
    d5 = 0;
    A  = 0;
  }

  Array5D(unsigned l1,
	  unsigned l2,
	  unsigned l3,
	  unsigned l4,
	  unsigned l5){
    d1 = l1;
    d2 = l2;
    d3 = l3;
    d4 = l4;
    d5 = l5;
    A  = 0;
    if (d1 > 0 &&
	d2 > 0 &&
	d3 > 0 &&
	d4 > 0 &&
	d5 > 0) A = new T[d1*d2*d3*d4*d5];
  }
 
  // destructor
  ~Array5D(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    d4 = 0;
    d5 = 0;
    if (A) delete[] A;
    A  = 0;
  }

  // allocate memory
  void allocate(unsigned l1,
		unsigned l2,
		unsigned l3,
		unsigned l4,
		unsigned l5){
    d1 = l1;
    d2 = l2;
    d3 = l3;
    d4 = l4;
    d5 = l5;
    if (!A &&
	d1 > 0 &&
	d2 > 0 &&
	d3 > 0 &&
	d4 > 0 &&
	d5 > 0) A = new T[d1*d2*d3*d4*d5];
  }

  // deallocate memory
  void deallocate(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    d4 = 0;
    d5 = 0;
    if (A) delete[] A;
    A = 0;
  }
  
  // indexing (parenthesis operator), two of them (for const correctness)
  const T& operator () (unsigned i,
			unsigned j,
			unsigned k,
			unsigned l,
			unsigned m) const{
    return A[d5*(d4*(d3*(d2*i+j)+k)+l)+m];}

  T& operator () (unsigned i,
		  unsigned j,
		  unsigned k,
		  unsigned l,
		  unsigned m){
    return A[d5*(d4*(d3*(d2*i+j)+k)+l)+m];}
  
  // get dims
  unsigned getD1() const {return d1;}
  unsigned getD2() const {return d2;}
  unsigned getD3() const {return d3;}
  unsigned getD4() const {return d4;}
  unsigned getD5() const {return d5;}

  
  // private data members
 private:

  unsigned d1;
  unsigned d2;
  unsigned d3;
  unsigned d4;
  unsigned d5;
  T*       A;


  // to prevent unwanted copying:
  Array5D(const Array5D<T>&);
  Array5D& operator = (const Array5D<T>&);
};
