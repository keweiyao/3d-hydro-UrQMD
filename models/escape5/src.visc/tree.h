class TTree ;

class MyTree{
 TTree *tree ;
 Float_t *X, *Y, *Z, *T, *Px, *Py, *Pz, *E ;
 Int_t *Id, *MId ;
 Short_t *Chrg, *Bar, *Strg ;
 Int_t nfill ;
public:
 MyTree(char *name) ;
 void fill(int iev) ;
} ;
