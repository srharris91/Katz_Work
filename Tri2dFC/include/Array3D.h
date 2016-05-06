template <typename T>
class Array3D
{
 public:
  // constructor
  Array3D(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    A  = 0;
  }

  Array3D(unsigned l1,
	  unsigned l2,
	  unsigned l3){
    d1 = l1;
    d2 = l2;
    d3 = l3;
    A  = 0;
    if (d1 > 0 &&
	d2 > 0 &&
	d3 > 0) A = new T[d1*d2*d3];
  }

  // destructor
  ~Array3D(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    if (A) delete[] A;
    A  = 0;
  }

  // allocate memory
  void allocate(unsigned l1,
		unsigned l2,
		unsigned l3){
    d1 = l1;
    d2 = l2;
    d3 = l3;
    if (!A &&
	d1 > 0 &&
	d2 > 0 &&
	d3 > 0) A = new T[d1*d2*d3];
  }

  // deallocate memory
  void deallocate(){
    d1 = 0;
    d2 = 0;
    d3 = 0;
    if (A) delete[] A;
    A = 0;
  }
  
  void set(T n){
      for(int i=0;i<d1*d2*d3;i++){
          A[i] = n;}
  }

  // indexing (parenthesis operator), two of them (for const correctness)
  const T& operator () (unsigned i,
			unsigned j,
			unsigned k) const{
    return A[d3*(d2*i+j)+k];}

  T& operator () (unsigned i,
		  unsigned j,
		  unsigned k){
    return A[d3*(d2*i+j)+k];}
  
  // get dims
  unsigned getD1() const {return d1;}
  unsigned getD2() const {return d2;}
  unsigned getD3() const {return d3;}

  
  // private data members
 private:

  unsigned d1;
  unsigned d2;
  unsigned d3;
  T*       A;


  // to prevent unwanted copying:
  Array3D(const Array3D<T>&);
  Array3D& operator = (const Array3D<T>&);
};
